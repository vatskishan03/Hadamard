// tower.c : resolve ONE (R,w)-tuple of the QH(358) search completely (G = Z2^4 x Z11; see README.md and engine.py).
//   input : a (R,w)-tuple = the complete member lists of A_t, B_t (each member: its 35 PAF class values, its line-1 and line-2 labels a1, a2 on
//           z = 0..5), and for C_t, D_t the row sums r (6 classes) and the column weights w (11 values); the targets T (35 columns), tau1, tau2.
//   level 1: line-1 labels a1(z) of C_t, D_t (label spaces from the column weights), classes of A_t,B_t members by a1; exact join on the 6 line-1
//            PAF columns -> children
//   level 2: inside a child, line-2 labels a2(z); classes of the child's A',B' members by a2; exact join on the 6 line-2 PAF columns -> grandchildren
//   leaf   : fibres of A'',B'' (member lists), C'',D'' (DFS enumeration of the letter words with row sums r, exact PAF class vectors by the table K),
//            exact 4-list join on all 35 columns (pair sums (A,C) hashed, (B,D) looked up, every key match re-checked on the 35 columns);
//            full matches are printed as "HIT iA iB letters(C) letters(D)".
//   Nothing in the tuple is sampled: a hit is a solution, a miss in a resolved tuple proves there is none.
// build: gcc -O2 -march=native -o tower tower.c
// usage: tower SPEC.bin OUT.txt
// SPEC: int32 m(35), nA, nB; per A member: int16 paf[m], int8 a1[6], int8 a2[6]; same for B; int32 T[m], tau1[6], tau2[6]; int16 K[64*64*6]; int32 colmap[36];
//       int8 LET[64*6]; int32 row1[6], row2[6]; for C,D: int32 r[6], w[11].
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <time.h>
typedef uint64_t u64; typedef int16_t i16; typedef uint32_t u32; typedef int32_t i32; typedef int8_t i8; typedef uint8_t u8;
#define NZ 11
static int M;
static i16 *K; static i32 colmap[36]; static i8 LET[64*6]; static int wL[64], a1L[64], a2L[64]; static i32 row1[6], row2[6];
static i32 T[64], tau1[6], tau2[6];
static u64 HW[64]; static u64 HW6[6];
static u64 sm(u64*x){ u64 z=(*x+=0x9E3779B97F4A7C15ULL); z=(z^(z>>30))*0xBF58476D1CE4E5B9ULL; z=(z^(z>>27))*0x94D049BB133111EBULL; return z^(z>>31);}
typedef struct { int n; i16 *paf; i8 *a1; i8 *a2; } Mem;
static Mem MA, MB;
static FILE *fo;
// ---------------- labels
typedef struct { int n; int nopt[NZ]; int optv[NZ][4]; i8 (*lab)[NZ]; i32 (*paf)[6]; } LabSp;
static void paf6_of(const i8*a11, i32*out){ for(int s=0;s<6;s++){ int t=0; for(int z=0;z<NZ;z++) t+=a11[z]*a11[(z+s)%NZ]; out[s]=t; } }
// allowed[z] = list (letters), val = value table; wfix: weight list
static void labsp_build(LabSp*ls,int (*allowed)[64],int *nallowed,const int*val){
  ls->n=1;
  for(int z=0;z<NZ;z++){ int k=0; for(int i=0;i<nallowed[z];i++){ int v=val[allowed[z][i]]; int f=0; for(int j=0;j<k;j++) if(ls->optv[z][j]==v){f=1;break;} if(!f){ if(k>=4){ fprintf(stderr,"too many options\n"); exit(9);} ls->optv[z][k++]=v; } }
    ls->nopt[z]=k; ls->n*=k; }
  ls->lab=malloc((size_t)ls->n*NZ); ls->paf=malloc((size_t)ls->n*6*4);
  for(int idx=0;idx<ls->n;idx++){ int x=idx; i8 a[NZ]; for(int z=0;z<NZ;z++){ a[z]=(i8)ls->optv[z][x%ls->nopt[z]]; x/=ls->nopt[z]; } memcpy(ls->lab[idx],a,NZ); paf6_of(a,ls->paf[idx]); }
}
static void labsp_free(LabSp*ls){ free(ls->lab); free(ls->paf); }
// ---------------- label join
typedef struct { int ia,ib,ic,id; } Quad;
static u64 h6(const i32*v){ u64 h=0; for(int j=0;j<6;j++) h+=HW6[j]*(u64)(int64_t)(v[j]+4096); return h; }
static int label_join(int nA,i32 (*PA)[6],int nB,i32 (*PB)[6],int nC,i32 (*PC)[6],int nD,i32 (*PD)[6],const i32*tau,Quad**out){
  int tb=1; while((1<<tb)<2*nD+16 && tb<30) tb++; int cap=1<<tb; int *head=malloc(sizeof(int)*(size_t)cap); int *nxt=malloc(sizeof(int)*((size_t)(nD>0?nD:0)+1)); for(int i=0;i<cap;i++) head[i]=-1;
  for(int d=0;d<nD;d++){ u64 h=h6(PD[d]); int p=(int)((h*0x9E3779B97F4A7C15ULL)>>(64-tb)); nxt[d]=head[p]; head[p]=d; }
  int nout=0,ocap=64; Quad*o=malloc(sizeof(Quad)*ocap);
  for(int a=0;a<nA;a++) for(int b=0;b<nB;b++){ i32 need[6]; for(int j=0;j<6;j++) need[j]=tau[j]-PA[a][j]-PB[b][j];
    for(int c=0;c<nC;c++){ i32 v[6]; for(int j=0;j<6;j++) v[j]=need[j]-PC[c][j]; u64 h=h6(v); int p=(int)((h*0x9E3779B97F4A7C15ULL)>>(64-tb));
      for(int d=head[p];d>=0;d=nxt[d]){ if(memcmp(v,PD[d],24)==0){ if(nout==ocap){ ocap*=2; o=realloc(o,sizeof(Quad)*ocap);} o[nout].ia=a;o[nout].ib=b;o[nout].ic=c;o[nout].id=d; nout++; } } } }
  free(head); free(nxt); *out=o; return nout;
}
// ---------------- member classes by label (a1 or a2 of 6 values)
static const i8 *g_lab; 
static int cmp6(const void*x,const void*y){ int i=*(const int*)x,j=*(const int*)y; int c=memcmp(g_lab+(size_t)i*6,g_lab+(size_t)j*6,6); return c?c:(i-j); }
typedef struct { int n; int *idx; int ncl; int *start; } Classes;      // idx sorted by label; class c = idx[start[c]..start[c+1])
static void classes_build(Classes*cl,const int*members,int n,const i8*lab6 /* indexed by member id: lab6[member*6] */){
  cl->n=n; cl->idx=malloc(sizeof(int)*(n+1)); memcpy(cl->idx,members,sizeof(int)*n); g_lab=lab6; qsort(cl->idx,n,sizeof(int),cmp6);
  cl->start=malloc(sizeof(int)*(n+2)); int nc=0; for(int i=0;i<n;i++){ if(i==0||memcmp(lab6+(size_t)cl->idx[i]*6,lab6+(size_t)cl->idx[i-1]*6,6)!=0) cl->start[nc++]=i; } cl->start[nc]=n; cl->ncl=nc; }
static void classes_free(Classes*cl){ free(cl->idx); free(cl->start); }
static void class_paf6(const Classes*cl,const i8*lab6,i32 (*out)[6]){ static const int sym[NZ]={0,1,2,3,4,5,5,4,3,2,1}; for(int c=0;c<cl->ncl;c++){ const i8*l=lab6+(size_t)cl->idx[cl->start[c]]*6; i8 a[NZ]; for(int z=0;z<NZ;z++) a[z]=l[sym[z]]; paf6_of(a,out[c]); } }
// ---------------- DFS fibre enumeration
typedef struct { i16 *V; u8 *Lt; long n, cap; } Fib;
static void fib_reset(Fib*f){ f->n=0; }
static void fib_emit(Fib*f,const i32*acc,const u8*let){
  if(f->n==f->cap){ f->cap=f->cap?f->cap*2:256; f->V=realloc(f->V,(size_t)f->cap*M*2); f->Lt=realloc(f->Lt,(size_t)f->cap*NZ); if(!f->V||!f->Lt){ fprintf(stderr,"oom\n"); exit(7);} }
  i16*v=f->V+(size_t)f->n*M; for(int c=0;c<6;c++) for(int s=0;s<6;s++){ int j=colmap[c*6+s]; if(j>=0) v[j]=(i16)acc[c*6+s]; }
  memcpy(f->Lt+(size_t)f->n*NZ,let,NZ); f->n++; }
static int (*g_allowed)[64]; static int *g_nallowed;
static void dfs(Fib*f,int z,const i32*rem,const i32*acc,u8*let){
  if(z==NZ){ fib_emit(f,acc,let); return; }
  for(int o=0;o<g_nallowed[z];o++){ int L=g_allowed[z][o]; i32 nr[6]; int ok=1; int lim=NZ-1-z;
    for(int c=0;c<6;c++){ nr[c]=rem[c]-LET[L*6+c]; if(nr[c]>lim||nr[c]<-lim){ ok=0; break; } }
    if(!ok) continue;
    i32 na[36]; memcpy(na,acc,sizeof na);
    for(int c=0;c<6;c++) na[c*6+0]+=K[((size_t)L*64+L)*6+c];
    for(int zp=0;zp<z;zp++){ int Lp=let[zp]; int d=z-zp; int s=d<NZ-d?d:NZ-d; for(int c=0;c<6;c++) na[c*6+s]+=K[((size_t)L*64+Lp)*6+c]; }
    let[z]=(u8)L; dfs(f,z+1,nr,na,let); }
}
// ---------------- deep join of four small lists (vectors int16[M]) ; calls emit on match
static i16 *gVA,*gVB,*gVC,*gVD; 
static long deep_join(int nA,const int*ia_list,const i16*VA,int nB,const int*ib_list,const i16*VB,int nC,const i16*VC,int nD,const i16*VD,int (*cb)(int,int,int,int)){
  // VA indexed through ia_list into global member array; VC,VD are fibre arrays
  double c1=(double)nA*nC+(double)nB*nD, c2=(double)nA*nD+(double)nB*nC;
  int sw=(c2<c1);
  // sets X=(A, C or D) stored ; Y=(B, D or C) probed
  int nX=nA, nY=nB, nXc=sw?nD:nC, nYc=sw?nC:nD; const i16*VXc=sw?VD:VC; const i16*VYc=sw?VC:VD;
  double pairs=(double)nX*nXc; int passbits=0; while(pairs/(double)((u64)1<<passbits) > 6e6 && passbits<12) passbits++;
  int tb=1; while(((u64)1<<tb) < 2.4*(pairs/(double)((u64)1<<passbits))+64) tb++;
  size_t cap=(size_t)1<<tb; static u64 *tk=0,*tv=0; static size_t tcap=0; if(tcap<cap){ free(tk); free(tv); tk=malloc(cap*8); tv=malloc(cap*8); tcap=cap; }
  u64 hT=0; for(int j=0;j<M;j++) hT+=HW[j]*(u64)(int64_t)T[j];
  long nsol=0;
  // hashes
  u64 *hA=malloc(sizeof(u64)*(nA+1)),*hB=malloc(sizeof(u64)*(nB+1)),*hXc=malloc(sizeof(u64)*(nXc+1)),*hYc=malloc(sizeof(u64)*(nYc+1));
  for(int i=0;i<nA;i++){ u64 x=0; const i16*v=VA+(size_t)ia_list[i]*M; for(int j=0;j<M;j++) x+=HW[j]*(u64)(int64_t)v[j]; hA[i]=x; }
  for(int i=0;i<nB;i++){ u64 x=0; const i16*v=VB+(size_t)ib_list[i]*M; for(int j=0;j<M;j++) x+=HW[j]*(u64)(int64_t)v[j]; hB[i]=x; }
  for(int i=0;i<nXc;i++){ u64 x=0; const i16*v=VXc+(size_t)i*M; for(int j=0;j<M;j++) x+=HW[j]*(u64)(int64_t)v[j]; hXc[i]=x; }
  for(int i=0;i<nYc;i++){ u64 x=0; const i16*v=VYc+(size_t)i*M; for(int j=0;j<M;j++) x+=HW[j]*(u64)(int64_t)v[j]; hYc[i]=x; }
  for(u64 pass=0;pass<((u64)1<<passbits);pass++){
    for(size_t i=0;i<cap;i++) tv[i]=~0ULL;
    for(int x=0;x<nX;x++) for(int c=0;c<nXc;c++){ u64 key=hA[x]+hXc[c]; if(passbits && (key>>(64-passbits))!=pass) continue; size_t p=(size_t)((key*0x9E3779B97F4A7C15ULL)>>(64-tb)); while(tv[p]!=~0ULL) p=(p+1)&(cap-1); tk[p]=key; tv[p]=((u64)x<<32)|(u64)c; }
    for(int y=0;y<nY;y++) for(int c=0;c<nYc;c++){ u64 key=hT-hB[y]-hYc[c]; if(passbits && (key>>(64-passbits))!=pass) continue; size_t p=(size_t)((key*0x9E3779B97F4A7C15ULL)>>(64-tb));
      while(tv[p]!=~0ULL){ if(tk[p]==key){ u64 id=tv[p]; int x=(int)(id>>32), cx=(int)(id&0xffffffffULL);
          const i16*a=VA+(size_t)ia_list[x]*M,*b=VB+(size_t)ib_list[y]*M,*cc=VXc+(size_t)cx*M,*dd=VYc+(size_t)c*M; int ok=1; for(int j=0;j<M;j++) if(a[j]+b[j]+cc[j]+dd[j]!=T[j]){ ok=0; break; }
          if(ok){ int ci=sw?c:cx, di=sw?cx:c; nsol+=cb(ia_list[x],ib_list[y],ci,di); } }
        p=(p+1)&(cap-1); } } }
  free(hA);free(hB);free(hXc);free(hYc); return nsol;
}
static Fib FC, FD; static long g_hits=0;
static int cb_hit(int a,int b,int c,int d){ fprintf(fo,"HIT %d %d",a,b); for(int z=0;z<NZ;z++) fprintf(fo," %d",FC.Lt[(size_t)c*NZ+z]); for(int z=0;z<NZ;z++) fprintf(fo," %d",FD.Lt[(size_t)d*NZ+z]); fprintf(fo,"\n"); g_hits++; return 1; }
int main(int argc,char**argv){
  FILE*f=fopen(argv[1],"rb"); fo=fopen(argv[2],"w"); i32 hd[3]; if(fread(hd,4,3,f)!=3) return 2; M=hd[0]; MA.n=hd[1]; MB.n=hd[2];
  Mem*mm[2]={&MA,&MB}; for(int q=0;q<2;q++){ Mem*m=mm[q]; m->paf=malloc((size_t)m->n*M*2+2); m->a1=malloc((size_t)m->n*6+8); m->a2=malloc((size_t)m->n*6+8);
    for(int i=0;i<m->n;i++){ if(fread(m->paf+(size_t)i*M,2,M,f)!=(size_t)M) return 3; if(fread(m->a1+(size_t)i*6,1,6,f)!=6) return 3; if(fread(m->a2+(size_t)i*6,1,6,f)!=6) return 3; } }
  if(fread(T,4,M,f)!=(size_t)M) return 4; if(fread(tau1,4,6,f)!=6) return 4; if(fread(tau2,4,6,f)!=6) return 4;
  K=malloc((size_t)64*64*6*2); if(fread(K,2,(size_t)64*64*6,f)!=(size_t)64*64*6) return 5;
  if(fread(colmap,4,36,f)!=36) return 6; if(fread(LET,1,64*6,f)!=64*6) return 6; if(fread(row1,4,6,f)!=6) return 6; if(fread(row2,4,6,f)!=6) return 6;
  i32 rC[6],rD[6],wC[NZ],wD[NZ]; if(fread(rC,4,6,f)!=6||fread(wC,4,NZ,f)!=NZ||fread(rD,4,6,f)!=6||fread(wD,4,NZ,f)!=NZ) return 7; fclose(f);
  u64 sd=12345; for(int j=0;j<64;j++) HW[j]=sm(&sd)|1; for(int j=0;j<6;j++) HW6[j]=sm(&sd)|1;
  static const i8 W6[6]={1,3,3,3,3,3};
  for(int L=0;L<64;L++){ int w=0,a1=0,a2=0; for(int c=0;c<6;c++){ w+=W6[c]*LET[L*6+c]; a1+=row1[c]*LET[L*6+c]; a2+=row2[c]*LET[L*6+c]; } wL[L]=w; a1L[L]=a1; a2L[L]=a2; }
  clock_t t0=clock(); double cover=0; long nchild=0,ngrand=0,nleaf=0,membC=0; double pairops=0;
  // ---- level 1
  int *allm=malloc(sizeof(int)*(MA.n+MB.n+2)); for(int i=0;i<MA.n;i++) allm[i]=i; Classes clA1,clB1; classes_build(&clA1,allm,MA.n,MA.a1); for(int i=0;i<MB.n;i++) allm[i]=i; classes_build(&clB1,allm,MB.n,MB.a1);
  i32 (*PA1)[6]=malloc(sizeof(i32[6])*(clA1.ncl+1)),(*PB1)[6]=malloc(sizeof(i32[6])*(clB1.ncl+1)); class_paf6(&clA1,MA.a1,PA1); class_paf6(&clB1,MB.a1,PB1);
  int allowedC[NZ][64],nalC[NZ],allowedD[NZ][64],nalD[NZ];
  for(int z=0;z<NZ;z++){ nalC[z]=0; nalD[z]=0; for(int L=0;L<64;L++){ if(wL[L]==wC[z]) allowedC[z][nalC[z]++]=L; if(wL[L]==wD[z]) allowedD[z][nalD[z]++]=L; } }
  LabSp lsC1,lsD1; labsp_build(&lsC1,allowedC,nalC,a1L); labsp_build(&lsD1,allowedD,nalD,a1L);
  Quad*ch; int nch=label_join(clA1.ncl,PA1,clB1.ncl,PB1,lsC1.n,lsC1.paf,lsD1.n,lsD1.paf,tau1,&ch);
  fprintf(fo,"# level1: A classes %d B classes %d C labels %d D labels %d -> %d children\n",clA1.ncl,clB1.ncl,lsC1.n,lsD1.n,nch);
  for(int ci=0;ci<nch;ci++){
    Quad q=ch[ci]; nchild++;
    int nA1=clA1.start[q.ia+1]-clA1.start[q.ia], nB1=clB1.start[q.ib+1]-clB1.start[q.ib]; int *mA1=clA1.idx+clA1.start[q.ia], *mB1=clB1.idx+clB1.start[q.ib];
    // level 2 label spaces: letters with (w, a1 = label value)
    int alC2[NZ][64],nlC2[NZ],alD2[NZ][64],nlD2[NZ];
    for(int z=0;z<NZ;z++){ nlC2[z]=0; nlD2[z]=0; for(int i=0;i<nalC[z];i++){ int L=allowedC[z][i]; if(a1L[L]==lsC1.lab[q.ic][z]) alC2[z][nlC2[z]++]=L; } for(int i=0;i<nalD[z];i++){ int L=allowedD[z][i]; if(a1L[L]==lsD1.lab[q.id][z]) alD2[z][nlD2[z]++]=L; } }
    LabSp lsC2,lsD2; labsp_build(&lsC2,alC2,nlC2,a2L); labsp_build(&lsD2,alD2,nlD2,a2L);
    Classes clA2,clB2; classes_build(&clA2,mA1,nA1,MA.a2); classes_build(&clB2,mB1,nB1,MB.a2);
    i32 (*PA2)[6]=malloc(sizeof(i32[6])*(clA2.ncl+1)),(*PB2)[6]=malloc(sizeof(i32[6])*(clB2.ncl+1)); class_paf6(&clA2,MA.a2,PA2); class_paf6(&clB2,MB.a2,PB2);
    Quad*gr; int ngr=label_join(clA2.ncl,PA2,clB2.ncl,PB2,lsC2.n,lsC2.paf,lsD2.n,lsD2.paf,tau2,&gr);
    for(int gi=0;gi<ngr;gi++){
      Quad g=gr[gi]; ngrand++;
      int alC3[NZ][64],nlC3[NZ],alD3[NZ][64],nlD3[NZ];
      for(int z=0;z<NZ;z++){ nlC3[z]=0; nlD3[z]=0; for(int i=0;i<nlC2[z];i++){ int L=alC2[z][i]; if(a2L[L]==lsC2.lab[g.ic][z]) alC3[z][nlC3[z]++]=L; } for(int i=0;i<nlD2[z];i++){ int L=alD2[z][i]; if(a2L[L]==lsD2.lab[g.id][z]) alD3[z][nlD3[z]++]=L; } }
      i32 acc[36]; memset(acc,0,sizeof acc); u8 let[16];
      fib_reset(&FC); g_allowed=alC3; g_nallowed=nlC3; dfs(&FC,0,rC,acc,let);
      if(FC.n==0) continue;
      fib_reset(&FD); g_allowed=alD3; g_nallowed=nlD3; dfs(&FD,0,rD,acc,let);
      if(FD.n==0) continue;
      int nA2=clA2.start[g.ia+1]-clA2.start[g.ia], nB2=clB2.start[g.ib+1]-clB2.start[g.ib]; const int*mA2=clA2.idx+clA2.start[g.ia],*mB2=clB2.idx+clB2.start[g.ib];
      nleaf++; membC+=FC.n+FD.n; cover+=(double)nA2*nB2*FC.n*FD.n; pairops+=(double)nA2*FC.n+(double)nB2*FD.n;
      deep_join(nA2,mA2,MA.paf,nB2,mB2,MB.paf,(int)FC.n,FC.V,(int)FD.n,FD.V,cb_hit);
    }
    free(gr); free(PA2); free(PB2); classes_free(&clA2); classes_free(&clB2); labsp_free(&lsC2); labsp_free(&lsD2);
  }
  fprintf(fo,"# STATS children %ld grandchildren(consistent) %ld leaves(non-empty) %ld enumerated_members %ld cover %.6g pairops %.6g hits %ld time %.2fs\n",nchild,ngrand,nleaf,membC,cover,pairops,g_hits,(double)(clock()-t0)/CLOCKS_PER_SEC);
  fclose(fo); return 0;
}
