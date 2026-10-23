// qspin.c: quaternary SPIN skew GS over Z_v: C = B o mu, D = B o mu^2 (mu^3 in H), A twisted-skew; H-orbit level.  Exhaustive and exact:
// the list of all skew blocks A (2 * 4^((n-1)/2), hashed by the PAF vector) against every B (4^(n-1), normalised B(0) = 1), with the
// row-sum filter |R_A|^2 + 3 |R_B|^2 = 4v and a first-class filter.  The plain reference version (n <= 13); qspin_fast.c is the
// same search for up to n = 21 orbits (equal counts: v = 7 <1>, mu = 2: COUNT 1820; v = 13, mu = 3: COUNT 110448).
// input (spin_drv.spin_shape): v n ; sz ; neg ; cperm[n] (class permutation induced by mu, classes 1..n-1, cperm[0]=0) ; operm[n] (orbit permutation) ; N tensor (n-1 classes)
// usage: qspin file maxprint
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
typedef unsigned long long u64;
static int v,n,d; static int sz[64],neg_[64],cperm[64],operm[64]; static int *Nt; static u64 P[128];
static inline int m4(int x){return x&3;}
static void pafc(const int*e,int c,int*re,int*im){
  int M[4]={0,0,0,0}; const int*N=Nt+(size_t)c*n*n;
  for(int a=0;a<n;a++){ const int*row=N+a*n; int ea=e[a]; for(int b=0;b<n;b++){ int x=row[b]; if(x) M[(ea-e[b])&3]+=x; } }
  *re=M[0]-M[2]; *im=M[1]-M[3];
}
typedef struct{ u64 key; u64 ex; int rre,rim; } El;
int main(int argc,char**argv){
  FILE*f=fopen(argv[1],"r"); int maxprint=argc>2?atoi(argv[2]):3;
  if(fscanf(f,"%d %d",&v,&n)!=2) return 1; d=n-1;
  for(int i=0;i<n;i++) if(fscanf(f,"%d",&sz[i])!=1) return 1;
  for(int i=0;i<n;i++) if(fscanf(f,"%d",&neg_[i])!=1) return 1;
  for(int i=0;i<n;i++) if(fscanf(f,"%d",&cperm[i])!=1) return 1;
  for(int i=0;i<n;i++) if(fscanf(f,"%d",&operm[i])!=1) return 1;
  Nt=malloc(sizeof(int)*(size_t)n*n*n);
  for(int c=1;c<n;c++) for(int a=0;a<n;a++) for(int b=0;b<n;b++) if(fscanf(f,"%d",&Nt[(size_t)c*n*n+a*n+b])!=1) return 1;
  fclose(f);
  u64 s=88172645463325252ULL; for(int i=0;i<128;i++){ s^=s<<13;s^=s>>7;s^=s<<17; P[i]=s|1ULL; }
  // A list
  int pr[64],np=0; for(int a=1;a<n;a++) if(a<neg_[a]) pr[np++]=a;
  long nA=2; for(int i=0;i<np;i++) nA*=4;
  El*A=malloc(sizeof(El)*nA);
  static char okR[1<<24]; memset(okR,0,sizeof(okR)); // |R_A|^2 values <= 4v
  long total=4L*v;
  // first-class bitmap of A's PAF at class 1 : we use value (re,im) offset by v
  int W=2*v+1; char*bm=calloc((size_t)W*W,1);
  int vec[256];
  for(long idx=0;idx<nA;idx++){
    long t=idx; int o0=t%2; t/=2; int e[64]; e[0]=o0?3:0;
    for(int i=0;i<np;i++){ int ej=t%4; t/=4; e[pr[i]]=ej; e[neg_[pr[i]]]=m4(1-ej); }
    u64 k=0; for(int c=1;c<n;c++){ int re,im; pafc(e,c,&re,&im); k+=(u64)(long long)re*P[2*(c-1)]+(u64)(long long)im*P[2*(c-1)+1]; if(c==1){ bm[(size_t)(-im+v)*W+(-re+v)]=1; } }
    int rr=0,ri=0; u64 ex=0; for(int a=0;a<n;a++){ int tt=m4(e[a]); ex|=((u64)tt)<<(2*a); switch(tt){case 0:rr+=sz[a];break;case 1:ri+=sz[a];break;case 2:rr-=sz[a];break;case 3:ri-=sz[a];break;} }
    A[idx].key=k;A[idx].ex=ex;A[idx].rre=rr;A[idx].rim=ri;
    long nn=(long)rr*rr+(long)ri*ri; if(nn<=total) okR[nn]=1;
  }
  long tsz=1; while(tsz<2*nA) tsz<<=1; unsigned *tab=calloc(tsz,sizeof(unsigned));
  for(long i=0;i<nA;i++){ u64 h=(A[i].key*0x9E3779B97F4A7C15ULL)>>20; h&=tsz-1; while(tab[h]) h=(h+1)&(tsz-1); tab[h]=(unsigned)i+1; }
  fprintf(stderr,"n=%d |A|=%ld\n",n,nA);
  // enumerate B with e0 = 0
  int e[64]; for(int i=0;i<n;i++) e[i]=0;
  long long count=0; int printed=0; long long nB=0,passDC=0,pass1=0;
  long long total_B=1; for(int i=1;i<n;i++) total_B*=4;
  // class reps: we need P(c) for c in {1, cperm[1], cperm[cperm[1]]}
  int c1=1,c2=cperm[1],c3=cperm[c2];
  for(long long idx=0;idx<total_B;idx++){
    long long t=idx; e[0]=0; for(int i=1;i<n;i++){ e[i]=t&3; t>>=2; }
    int rr=0,ri=0; for(int a=0;a<n;a++){ switch(e[a]){case 0:rr+=sz[a];break;case 1:ri+=sz[a];break;case 2:rr-=sz[a];break;case 3:ri-=sz[a];break;} }
    long nn=3L*((long)rr*rr+(long)ri*ri); if(nn>total) continue; if(!okR[total-nn]) continue;
    passDC++;
    int re1,im1,re2,im2,re3,im3; pafc(e,c1,&re1,&im1); pafc(e,c2,&re2,&im2); pafc(e,c3,&re3,&im3);
    int sre=re1+re2+re3, sim=im1+im2+im3; // S(c1)
    if(sre<-v||sre>v||sim<-v||sim>v) continue;
    if(!bm[(size_t)(sim+v)*W+(sre+v)]) continue;   // need PAF_A(c1) = -S(c1): bitmap stored at (-im,-re) offset => index (-imA+v)*W+(-reA+v) = (S im + v)*W+(S re+v)
    pass1++;
    // full
    int pv[64][2]; for(int c=1;c<n;c++){ pafc(e,c,&pv[c][0],&pv[c][1]); }
    u64 k=0; for(int c=1;c<n;c++){ int c2_=cperm[c],c3_=cperm[c2_]; int sr=pv[c][0]+pv[c2_][0]+pv[c3_][0], si=pv[c][1]+pv[c2_][1]+pv[c3_][1]; k+=(u64)(long long)sr*P[2*(c-1)]+(u64)(long long)si*P[2*(c-1)+1]; }
    u64 kk=-k; u64 h=(kk*0x9E3779B97F4A7C15ULL)>>20; h&=tsz-1;
    while(tab[h]){ if(A[tab[h]-1].key==kk){
        // verify exactly: A + B + B o mu + B o mu^2
        long ai=tab[h]-1; int ea[64]; for(int i=0;i<n;i++) ea[i]=(A[ai].ex>>(2*i))&3;
        int ec[64],ed[64]; for(int i=0;i<n;i++){ ec[i]=e[operm[i]]; ed[i]=e[operm[operm[i]]]; }
        int ok=1; for(int c=1;c<n&&ok;c++){ int r0,i0,r1,i1,r2,i2,r3,i3; pafc(ea,c,&r0,&i0); pafc(e,c,&r1,&i1); pafc(ec,c,&r2,&i2); pafc(ed,c,&r3,&i3); if(r0+r1+r2+r3||i0+i1+i2+i3) ok=0; }
        if(ok){ count++; if(printed<maxprint){ printed++; printf("SOL A="); for(int i=0;i<n;i++) printf("%d",ea[i]); printf(" B="); for(int i=0;i<n;i++) printf("%d",e[i]); printf("\n"); fflush(stdout);} }
      }
      h=(h+1)&(tsz-1); }
  }
  printf("COUNT %lld  (B normalised e0=0; passDC %lld pass1 %lld of %lld)\n",count,passDC,pass1,total_B);
  return 0;
}
