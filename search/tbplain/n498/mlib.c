// mlib.c : the M part of QH(498): sign lists of the m-part (build_block, gen_list), exact four-list join (join4), partial-agreement counts,
// the Gaussian proxy of the number of completions.  #included by band.c; also built as a shared object for test_mlib.py:
//   gcc -O2 -march=native -shared -fPIC -o mlib.so mlib.c -lm
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
typedef uint64_t u64; typedef int16_t i16; typedef int8_t i8; typedef uint32_t u32;
#define P 41
#define NC 13
#define LANES 16
static int orb_of[P], NO=0, OREP[9], NEGO[9];
static int NN[9][9][9];
static int WQ[NC][27][27];
static u64 WH[LANES];
static int inited=0;
static u64 sm(u64*x){ u64 z=(*x+=0x9E3779B97F4A7C15ULL); z=(z^(z>>30))*0xBF58476D1CE4E5B9ULL; z=(z^(z>>27))*0x94D049BB133111EBULL; return z^(z>>31);}
static const int MC_B[NC]={0,0,0,0,1,1,1,1,1,1,1,1,1};
static const int MC_O[NC]={1,2,3,6,0,1,2,3,4,5,6,7,8};
void minit(void){
  if(inited) return; inited=1;
  for(int z=0;z<P;z++) orb_of[z]=-1;
  for(int z=0;z<P;z++) if(orb_of[z]<0){ int w=z; do{ orb_of[w]=NO; w=(w*10)%P; }while(w!=z); NO++; }
  for(int o=0;o<NO;o++){ for(int z=0;z<P;z++) if(orb_of[z]==o){ OREP[o]=z; break; } }
  for(int o=0;o<NO;o++) NEGO[o]=orb_of[(P-OREP[o])%P];
  memset(NN,0,sizeof NN);
  for(int op=0;op<NO;op++){ int g=OREP[op]; for(int z=0;z<P;z++) NN[op][orb_of[z]][orb_of[(z+g)%P]]++; }
  memset(WQ,0,sizeof WQ);
  for(int c=0;c<NC;c++) for(int b=0;b<3;b++) for(int o=0;o<9;o++) for(int o2=0;o2<9;o2++) WQ[c][9*b+o][9*((b+MC_B[c])%3)+o2]+=NN[MC_O[c]][o][o2];
  u64 s=777; for(int j=0;j<LANES;j++) WH[j]=sm(&s)|1;
}
typedef struct { int n; int nv[32][2]; int np[32]; i16 diag[LANES]; i8 K[32][32][LANES]; } Block;   // K[i][j][c] symmetric, zero on diagonal
// typ 0: C/D type (free positions where p==0), typ 1: symmetric (A/B)
int build_block(int typ, const int *p /*27 ints: p[9b+o]*/, Block *B){
  memset(B,0,sizeof *B);
  int n=0; int pos[32][2]; int np_[32];
  if(typ==0){ for(int b=0;b<3;b++) for(int o=0;o<9;o++) if(p[9*b+o]==0){ pos[n][0]=9*b+o; np_[n]=1; n++; } }
  else {
    for(int o=0;o<9;o++){ int q=NEGO[o]; if(o<=q && p[o]==0){ pos[n][0]=o; if(q!=o){ pos[n][1]=q; np_[n]=2; } else np_[n]=1; n++; } }
    for(int o=0;o<9;o++) if(p[9+o]==0){ pos[n][0]=9+o; pos[n][1]=18+NEGO[o]; np_[n]=2; n++; }
  }
  if(n>30) return -1;
  B->n=n;
  for(int c=0;c<NC;c++){
    for(int i=0;i<n;i++) for(int j=0;j<n;j++){
      int q=0; for(int a=0;a<np_[i];a++) for(int b=0;b<np_[j];b++) q+=WQ[c][pos[i][a]][pos[j][b]];
      if(i==j) B->diag[c]+=q; else { B->K[i][j][c]+=q; B->K[j][i][c]+=q; }   // symmetrised: K[i][j] = Q[i][j]+Q[j][i]
    }
  }
  for(int i=0;i<n;i++){ B->nv[i][0]=pos[i][0]; B->nv[i][1]=(np_[i]>1)?pos[i][1]:-1; B->np[i]=np_[i]; }
  return n;
}
// all 2^n sign vectors in Gray order: list index m <-> pattern gray(m) = m ^ (m>>1) (bit i set <=> s_i = -1)
void gen_list(const Block *B, i16 *out){
  int n=B->n; int s[32]; for(int i=0;i<n;i++) s[i]=1;
  i16 v[LANES]; memset(v,0,sizeof v);
  for(int c=0;c<NC;c++){ int x=B->diag[c]; for(int i=0;i<n;i++) for(int j=i+1;j<n;j++) x+=B->K[i][j][c]; v[c]=(i16)x; }
  u64 N=(u64)1<<n;
  memcpy(out,v,sizeof v);
  for(u64 m=1;m<N;m++){
    int i=__builtin_ctzll(m);
    i16 acc[LANES]; memset(acc,0,sizeof acc);
    for(int j=0;j<n;j++){ if(j==i) continue; int sj=s[j]; const i8*k=B->K[i][j]; for(int c=0;c<LANES;c++) acc[c]+=sj*k[c]; }
    int si=s[i];
    for(int c=0;c<LANES;c++) v[c]-=(i16)(2*si*acc[c]);
    s[i]=-si;
    memcpy(out+m*LANES,v,sizeof v);
  }
}
static inline u64 hvec(const i16*v){ u64 h=0; for(int c=0;c<NC;c++) h+=WH[c]*(u64)(int64_t)v[c]; return h; }
// 4-list exact join: V0[a]+V1[b]+V2[c]+V3[d] == T (NC columns) ; pairing (0,2) vs (1,3).  returns number of solutions (stored up to maxout as 4 ints), *nmatch = key matches
long join4(const i16*V0,long n0,const i16*V1,long n1,const i16*V2,long n2,const i16*V3,long n3,const i16*T,int*out,long maxout,long*nmatch){
  u64*h0=malloc(n0*8),*h1=malloc(n1*8),*h2=malloc(n2*8),*h3=malloc(n3*8);
  for(long i=0;i<n0;i++) h0[i]=hvec(V0+i*LANES); for(long i=0;i<n1;i++) h1[i]=hvec(V1+i*LANES);
  for(long i=0;i<n2;i++) h2[i]=hvec(V2+i*LANES); for(long i=0;i<n3;i++) h3[i]=hvec(V3+i*LANES);
  u64 hT=hvec(T);
  u64 npairs=(u64)n0*n2; int tb=4; while(((u64)1<<tb)<2*npairs+16) tb++;
  size_t cap=(size_t)1<<tb;
  u64*tk=malloc(cap*8); u64*tv=malloc(cap*8); if(!tk||!tv) return -1;
  for(size_t i=0;i<cap;i++) tv[i]=~0ULL;
  for(long a=0;a<n0;a++){ u64 ha=h0[a]; for(long c=0;c<n2;c++){ u64 key=ha+h2[c]; size_t p=(size_t)((key*0x9E3779B97F4A7C15ULL)>>(64-tb)); while(tv[p]!=~0ULL) p=(p+1)&(cap-1); tk[p]=key; tv[p]=(u64)a*(u64)n2+(u64)c; } }
  long ns=0,nm=0;
  for(long b=0;b<n1;b++){ u64 hb=hT-h1[b]; for(long d=0;d<n3;d++){ u64 key=hb-h3[d]; size_t p=(size_t)((key*0x9E3779B97F4A7C15ULL)>>(64-tb));
      while(tv[p]!=~0ULL){ if(tk[p]==key){ nm++; u64 id=tv[p]; long a=(long)(id/(u64)n2), c=(long)(id%(u64)n2); int ok=1;
          for(int j=0;j<NC;j++){ int x=V0[a*LANES+j]+V1[b*LANES+j]+V2[c*LANES+j]+V3[d*LANES+j]; if(x!=T[j]){ok=0;break;} }
          if(ok){ if(ns<maxout){ out[4*ns]=(int)a; out[4*ns+1]=(int)b; out[4*ns+2]=(int)c; out[4*ns+3]=(int)d; } ns++; } }
        p=(p+1)&(cap-1); } } }
  *nmatch=nm; free(h0);free(h1);free(h2);free(h3);free(tk);free(tv); return ns;
}
// number of quadruples agreeing on the first j columns of an arbitrary column order (cols[0..j-1]); counting only (hash multiset)
long count_partial(const i16*V0,long n0,const i16*V1,long n1,const i16*V2,long n2,const i16*V3,long n3,const int*cols,int j,const i16*T){
  u64 W[LANES]; for(int t=0;t<j;t++) W[t]=WH[cols[t]];
  u64*h[4]; const i16*V[4]={V0,V1,V2,V3}; long nn[4]={n0,n1,n2,n3};
  for(int k=0;k<4;k++){ h[k]=malloc(nn[k]*8); for(long i=0;i<nn[k];i++){ u64 x=0; for(int t=0;t<j;t++) x+=W[t]*(u64)(int64_t)V[k][i*LANES+cols[t]]; h[k][i]=x; } }
  u64 hT=0; for(int t=0;t<j;t++) hT+=W[t]*(u64)(int64_t)T[cols[t]];
  u64 npairs=(u64)n0*n2; int tb=4; while(((u64)1<<tb)<2*npairs+16) tb++;
  size_t cap=(size_t)1<<tb; u64*tk=malloc(cap*8); u32*tc=calloc(cap,4);
  for(long a=0;a<n0;a++){ u64 ha=h[0][a]; for(long c=0;c<n2;c++){ u64 key=ha+h[2][c]; size_t p=(size_t)((key*0x9E3779B97F4A7C15ULL)>>(64-tb)); while(tc[p]&&tk[p]!=key) p=(p+1)&(cap-1); tk[p]=key; tc[p]++; } }
  long tot=0;
  for(long b=0;b<n1;b++){ u64 hb=hT-h[1][b]; for(long d=0;d<n3;d++){ u64 key=hb-h[3][d]; size_t p=(size_t)((key*0x9E3779B97F4A7C15ULL)>>(64-tb)); while(tc[p]){ if(tk[p]==key){ tot+=tc[p]; break; } p=(p+1)&(cap-1); } } }
  for(int k=0;k<4;k++) free(h[k]); free(tk); free(tc); return tot;
}

// Gaussian proxy (lattice cell 1) of the number of M-completions: log10 E_t = log10( 2^(sum n) * density at 0 of the 13-column sum )
#include <math.h>
double proxy4(const Block*B0,const Block*B1,const Block*B2,const Block*B3,double*zz,double*logdet_out){
  const Block*Bs[4]={B0,B1,B2,B3};
  double d[NC]; double cov[NC][NC]; int ntot=0;
  memset(d,0,sizeof d); memset(cov,0,sizeof cov);
  for(int k=0;k<4;k++){ const Block*B=Bs[k]; ntot+=B->n; for(int c=0;c<NC;c++) d[c]+=B->diag[c];
    for(int i=0;i<B->n;i++) for(int j=i+1;j<B->n;j++){ const i8*kk=B->K[i][j]; for(int c=0;c<NC;c++){ double x=kk[c]; if(x==0) continue; for(int c2=c;c2<NC;c2++) cov[c][c2]+=x*kk[c2]; } } }
  for(int c=0;c<NC;c++) for(int c2=0;c2<c;c2++) cov[c][c2]=cov[c2][c];
  // Cholesky
  double L[NC][NC]; memset(L,0,sizeof L);
  for(int i=0;i<NC;i++){ for(int j=0;j<=i;j++){ double s=cov[i][j]; for(int k=0;k<j;k++) s-=L[i][k]*L[j][k];
      if(i==j){ if(s<=1e-9) return -1e30; L[i][i]=sqrt(s); } else L[i][j]=s/L[j][j]; } }
  double y[NC]; double z2=0, ld=0;
  for(int i=0;i<NC;i++){ double s=-d[i]; for(int k=0;k<i;k++) s-=L[i][k]*y[k]; y[i]=s/L[i][i]; z2+=y[i]*y[i]; ld+=2*log(L[i][i]); }
  if(zz) *zz=z2; if(logdet_out) *logdet_out=ld;
  double ln=ntot*log(2.0)-0.5*(NC*log(2*M_PI)+ld)-0.5*z2;
  return ln/log(10.0);
}

// Monte-Carlo count: sum over nsamp random (a,b) of #{(c,d): V2[c]+V3[d] = T - V0[a] - V1[b]} on the columns cols[0..j-1]
long sample_count(const i16*V0,long n0,const i16*V1,long n1,const i16*V2,long n2,const i16*V3,long n3,const int*cols,int j,const i16*T,long nsamp,u64 seed,long*nev){
  u64 W[LANES]; for(int t=0;t<j;t++) W[t]=WH[cols[t]];
  #define HH(V,i) ({ u64 x_=0; for(int t_=0;t_<j;t_++) x_+=W[t_]*(u64)(int64_t)(V)[(long)(i)*LANES+cols[t_]]; x_; })
  u64*h2=malloc(n2*8); for(long i=0;i<n2;i++) h2[i]=HH(V2,i);
  u64 hT=0; for(int t=0;t<j;t++) hT+=W[t]*(u64)(int64_t)T[cols[t]];
  int tb=4; while(((u64)1<<tb)<2*(u64)n3+16) tb++; size_t cap=(size_t)1<<tb;
  u64*tk=malloc(cap*8); u32*tc=calloc(cap,4);
  for(long d=0;d<n3;d++){ u64 key=HH(V3,d); size_t p=(size_t)((key*0x9E3779B97F4A7C15ULL)>>(64-tb)); while(tc[p]&&tk[p]!=key) p=(p+1)&(cap-1); tk[p]=key; tc[p]++; }
  u64 rs=seed; long tot=0, ev=0;
  for(long s=0;s<nsamp;s++){ long a=sm(&rs)%n0, b=sm(&rs)%n1; u64 need=hT-HH(V0,a)-HH(V1,b); long here=0;
    for(long c=0;c<n2;c++){ u64 key=need-h2[c]; size_t p=(size_t)((key*0x9E3779B97F4A7C15ULL)>>(64-tb)); while(tc[p]){ if(tk[p]==key){ here+=tc[p]; break; } p=(p+1)&(cap-1); } }
    tot+=here; if(here) ev++; }
  *nev=ev; free(h2); free(tk); free(tc); return tot;
}
