// tri_join.c: partitioned four-list join for the triangular F1 family v = 1 + k(k+1)/2 (search/NOTES.md section 88; section 92).
// v prime = 1 + k(k+1)/2, H <= Z_v* of order k, e = (v-1)/k cosets C_j = g^j H; quaternary blocks constant on cosets.
// Pattern: A twisted-skew (u(-x) = i conj u(x), u(0) in {1,-i}) with coset sum S_A = m(1+i), m = +-1 aligned (big, |R_A|^2 = 4v-3),
// B, C, D with u(0) = 1 and coset sum 0 (small).  Condition: sum of the four complex PAFs = 0 at every s != 0.
// The join is partitioned by a linear residue pi(PAF) mod M: left pairs (A,B) with pi = r, right pairs (C,D) with pi = -r.
// usage: tri_join v k g M r0 r1 maxsol      (g a primitive root mod v; residues r0..r1-1)
// output: SOL lines with the coset exponents (u(0), c_0..c_{e-1}) of A B C D; STAT lines.
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
typedef unsigned long long u64; typedef long long i64;
static int v,k,e,n,np_;            // n = e+1 orbits (0 = {0}, j+1 = C_j); np_ = e/2 complex PAF classes
static int *orb;                   // orb[x] in 0..e
static int *Ncnt;                  // Ncnt[c][a][b] for the np_ class representatives c
static int crep[64];
static u64 KR[128], QR[128];
typedef struct { u64 key; u64 pi; u64 ex; } El;   // ex: 2 bits per orbit (n <= 32)
static int paf(const int *x, int *out){ // x[0..n-1] exponents
  for(int t=0;t<np_;t++){ int M[4]={0,0,0,0}; const int *N=Ncnt+(size_t)t*n*n;
    for(int a=0;a<n;a++) for(int b=0;b<n;b++){ int c=N[a*n+b]; if(c) M[(x[a]-x[b])&3]+=c; }
    out[2*t]=M[0]-M[2]; out[2*t+1]=M[1]-M[3]; }
  return 0;
}
static void fill(El *E, const int *x, u64 Mmask){
  int p[128]; paf(x,p); u64 kk=0,q=0; for(int i=0;i<2*np_;i++){ kk+=(u64)(i64)p[i]*KR[i]; q+=(u64)(i64)p[i]*QR[i]; }
  E->key=kk; E->pi=q&Mmask; u64 ex=0; for(int a=0;a<n;a++) ex|=((u64)(x[a]&3))<<(2*a); E->ex=ex;
}
static El *A,*L; static long nA=0,nL=0;
static void genL(int pos,int *x,int c0,int c1,int c2,int c3,int h, u64 Mm){ // counts so far; need c0=c2,c1=c3 at the end
  if(pos==n){ if(c0==c2&&c1==c3){ fill(&L[nL++],x,Mm);} return; }
  int rem=n-pos;
  for(int t=0;t<4;t++){ int d0=c0+(t==0),d1=c1+(t==1),d2=c2+(t==2),d3=c3+(t==3);
    if(abs(d0-d2)+abs(d1-d3)>rem-1) continue; x[pos]=t; genL(pos+1,x,d0,d1,d2,d3,h,Mm); }
}
int main(int argc,char**argv){
  v=atoi(argv[1]); k=atoi(argv[2]); int g=atoi(argv[3]); long M=atol(argv[4]); long r0=atol(argv[5]), r1=atol(argv[6]); long maxsol=atol(argv[7]);
  e=(v-1)/k; n=e+1; np_=e/2; u64 Mm=(u64)M-1;
  orb=malloc(sizeof(int)*v); orb[0]=0;
  long hgen=1; for(int i=0;i<e;i++) hgen=hgen*g%v; // g^e generates H
  long gj=1; for(int j=0;j<e;j++){ long y=gj; for(int i=0;i<k;i++){ orb[y]=j+1; y=y*hgen%v; } gj=gj*g%v; }
  for(int x=1;x<v;x++) if(orb[x]<1){ fprintf(stderr,"orbit error\n"); return 1; }
  // class representatives c = g^t, t = 0..np_-1 (classes t and t+e/2 are conjugate)
  Ncnt=calloc((size_t)np_*n*n,sizeof(int));
  gj=1; for(int t=0;t<np_;t++){ crep[t]=(int)gj; for(int x=0;x<v;x++) Ncnt[(size_t)t*n*n+orb[x]*n+orb[(x+gj)%v]]++; gj=gj*g%v; }
  u64 s=88172645463325252ULL; for(int i=0;i<128;i++){ s^=s<<13;s^=s>>7;s^=s<<17; KR[i]=s|1ULL; s^=s<<13;s^=s>>7;s^=s<<17; QR[i]=s; }
  // A list: u0 = 1 with m = +1, u0 = -i (exp 3) with m = -1; pairs (j, j+e/2): c in {0,1} -> +(1+i), {2,3} -> -(1+i); partner 1-c
  int h=e/2; long cap=2; for(int i=0;i<h;i++) cap*=4; A=malloc(sizeof(El)*cap);
  int x[64];
  for(int o0=0;o0<2;o0++){ int m=o0?-1:1; int npl=(h+m)/2;
    long tot=1; for(int i=0;i<h;i++) tot*=4;
    for(long idx=0;idx<tot;idx++){ long t=idx; int pl=0; for(int i=0;i<h;i++){ int c=t%4; t/=4; x[1+i]=c; x[1+i+h]=(1-c)&3; if(c<2) pl++; }
      if(pl!=npl) continue; x[0]=o0?3:0; fill(&A[nA++],x,Mm); } }
  // L list: u0 = 1, coset sum 0
  long capL=1; { // C(e,h)^2 upper bound
    double c=1; for(int i=1;i<=h;i++) c=c*(e-h+i)/i; capL=(long)(c*c)+16; }
  L=malloc(sizeof(El)*capL); x[0]=0; genL(1,x,0,0,0,0,h,Mm);
  fprintf(stderr,"v=%d k=%d e=%d |A|=%ld |L|=%ld M=%ld\n",v,k,e,nA,nL,M);
  // bucket L by pi
  long *cnt=calloc(M+1,sizeof(long)); for(long i=0;i<nL;i++) cnt[L[i].pi+1]++; for(long r=0;r<M;r++) cnt[r+1]+=cnt[r];
  long *ord=malloc(sizeof(long)*nL); long *pos=malloc(sizeof(long)*M); memcpy(pos,cnt,sizeof(long)*M);
  for(long i=0;i<nL;i++) ord[pos[L[i].pi]++]=i;
  long found=0; long long probes=0;
  for(long r=r0;r<r1 && found<maxsol;r++){
    // left table: (a,b), pi(a)+pi(b) = r
    long nl=0; for(long i=0;i<nA;i++){ u64 rb=((u64)r-A[i].pi)&Mm; nl+=cnt[rb+1]-cnt[rb]; }
    long tsz=1; while(tsz<2*nl+2) tsz<<=1;
    u64 *tk=calloc(tsz,8); unsigned *ta=calloc(tsz,4), *tb=calloc(tsz,4);
    for(long i=0;i<nA;i++){ u64 rb=((u64)r-A[i].pi)&Mm; for(long q=cnt[rb];q<cnt[rb+1];q++){ long j=ord[q]; u64 kk=A[i].key+L[j].key;
        u64 hh=(kk*0x9E3779B97F4A7C15ULL)>>17; hh&=tsz-1; while(ta[hh]) hh=(hh+1)&(tsz-1); tk[hh]=kk; ta[hh]=(unsigned)i+1; tb[hh]=(unsigned)j+1; } }
    // right: (c,d) with pi(c)+pi(d) = -r
    for(long ic=0;ic<nL && found<maxsol;ic++){ u64 rd=((u64)(-r)-L[ic].pi)&Mm;
      for(long q=cnt[rd];q<cnt[rd+1];q++){ long id=ord[q]; u64 kk=-(L[ic].key+L[id].key); probes++;
        u64 hh=(kk*0x9E3779B97F4A7C15ULL)>>17; hh&=tsz-1;
        while(ta[hh]){ if(tk[hh]==kk){ // verify exactly
            int xs[4][64]; El *E4[4]={&A[ta[hh]-1],&L[tb[hh]-1],&L[ic],&L[id]}; int ok=1; int sum[128]; memset(sum,0,sizeof(sum));
            for(int z=0;z<4;z++){ for(int a=0;a<n;a++) xs[z][a]=(E4[z]->ex>>(2*a))&3; int p[128]; paf(xs[z],p); for(int i=0;i<2*np_;i++) sum[i]+=p[i]; }
            for(int i=0;i<2*np_;i++) if(sum[i]) ok=0;
            if(ok){ found++; printf("SOL r=%ld",r); for(int z=0;z<4;z++){ printf(" "); for(int a=0;a<n;a++) printf("%d",xs[z][a]); } printf("\n"); fflush(stdout); if(found>=maxsol) break; }
          } hh=(hh+1)&(tsz-1); } }
    }
    printf("STAT r=%ld left=%ld probes=%lld found=%ld\n",r,nl,probes,found); fflush(stdout);
    free(tk); free(ta); free(tb);
  }
  return 0;
}
