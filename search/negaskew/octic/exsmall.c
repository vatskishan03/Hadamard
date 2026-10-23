// exsmall.c: EXHAUSTIVE count of the octic skew Goethals-Seidel quadruples of negacyclic length n = 4v <= 24 (order 16v, v = 1, 3, 5), no multiplier group:
// A palindromic (a_0 = 1), B, C, D free +-1 sequences, counted through the canonical representatives of their orbits under the monomials +-x^s
// (meet in the middle on the keys of the negaperiodic autocorrelations of lags 1 .. n/2-1; every hit is a solution).  Prints
//   COUNT                     ordered (A, B, C, D) with B, C, D canonical representatives
//   TOTAL_ORDERED_FULL        the number of ordered quadruples of all +-1 sequences, weighting each representative by its orbit size
//   DOUBLING_TEST tested/doubled   among every dstride-th solution those made of two pairs (X, m theta(X)), theta: x -> -x, m a monomial
//                                  (with adjoints), i.e. the images of order-8v quadruples with a centre-palindromic block (NOTES.md section 89)
// usage: exsmall n [maxsols [dstride]]     n = 4, 12, 20; maxsols = number of SOL lines printed (default 2000); dstride = 1: test every solution
//   (default), k: every k-th, 0: no test.   Expected: n = 4  TOTAL_ORDERED_FULL 6144;  n = 12  2057674752  (COUNT 148848);
//   n = 20  75652141056000  (COUNT 1182064704; a few minutes, about 650 MB).
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
typedef uint64_t u64; typedef uint32_t u32; typedef uint16_t u16;
static int n,h; static u64 R[32];
static inline u32 rot1(u32 a){ // multiply by x: (x a)_k = a_{k-1}, (x a)_0 = -a_{n-1}
  u32 top=(a>>(n-1))&1; u32 r=((a<<1)&((1u<<n)-1)); r|=(top^1); return r; }
static u32 canon(u32 a){ u32 best=a, c=a, mask=(1u<<n)-1; for(int i=0;i<n;i++){ c=rot1(c); if(c<best) best=c; u32 d=c^mask; if(d<best) best=d; } if((a^mask)<best) best=a^mask; return best; }

static u32 theta_b(u32 a){ u32 m=0; for(int k=1;k<n;k+=2) m|=1u<<k; return a^m; }   // x -> -x : flips odd positions
static u32 adj_b(u32 a){ u32 r=a&1u; for(int k=1;k<n;k++){ u32 b=(a>>(n-k))&1; r|=(b^1u)<<k; } return r; }
static int paired(u32 X,u32 Y){ u32 t=theta_b(X); return canon(t)==canon(Y) || canon(adj_b(t))==canon(Y); }
static int is_doubled(u32 A,u32 B,u32 C,u32 D){ u32 ta=theta_b(A); u32 c1=canon(ta), c2=canon(adj_b(ta)); u32 xs[3]={B,C,D};
  for(int i=0;i<3;i++){ if(canon(xs[i])==c1||canon(xs[i])==c2){ u32 Y=xs[(i+1)%3],Z=xs[(i+2)%3]; if(paired(Y,Z)) return 1; } } return 0; }
static void paf(u32 a,int*p){ // p[s-1] for s=1..h ; entries +-1: bit=1 -> -1
  int x[32]; for(int i=0;i<n;i++) x[i]=((a>>i)&1)?-1:1;
  for(int s=1;s<=h;s++){ int t=0; for(int i=0;i<n;i++){ if(i+s<n) t+=x[i]*x[i+s]; else t-=x[i]*x[i+s-n]; } p[s-1]=t; } }
static u64 keyof(u32 a){ int p[32]; paf(a,p); u64 k=0; for(int s=0;s<h;s++) k+=(u64)(int64_t)p[s]*R[s]; return k; }
typedef struct { u64 key; u32 ib; u16 ia; u16 used; } Ent;
int main(int argc,char**argv){
  n=atoi(argv[1]); long maxs=argc>2?atol(argv[2]):2000; h=n/2-1;
  u64 s=88172645463325252ULL; for(int i=0;i<32;i++){ s^=s<<13;s^=s>>7;s^=s<<17; R[i]=s|1ULL; }
  u32 N=1u<<n; u32*cl=malloc(sizeof(u32)*N); long nc=0;
  for(u32 a=0;a<N;a++) if(canon(a)==a) cl[nc++]=a;
  u64*ck=malloc(sizeof(u64)*nc); int*osz=malloc(sizeof(int)*nc); for(long i=0;i<nc;i++){ u32 a=cl[i],c=a,mask=(1u<<n)-1; int stab=0; for(int t=0;t<n;t++){ if(c==a) stab++; if((c^mask)==a) stab++; c=rot1(c);} osz[i]=(2*n)/stab; } for(long i=0;i<nc;i++) ck[i]=keyof(cl[i]);
  // palindromic list: bits s=1..n/2 free, a_{n-s}=a_s, a_0 = 0 (bit 0 => +1)
  int hn=n/2; long na=1L<<hn; u32*al=malloc(sizeof(u32)*na); u64*ak=malloc(sizeof(u64)*na);
  for(long t=0;t<na;t++){ u32 a=0; for(int q=1;q<=hn;q++){ int b=(t>>(q-1))&1; if(b){ a|=1u<<q; a|=1u<<(n-q);} } al[t]=a; ak[t]=keyof(a); }
  fprintf(stderr,"n=%d canonical=%ld palindromic=%ld\n",n,nc,na);
  long m=na*nc; long tsz=1; while(tsz<2*m) tsz<<=1;
  Ent*tab=calloc(tsz,sizeof(Ent)); if(!tab){fprintf(stderr,"alloc\n");return 1;}
  for(long ia=0;ia<na;ia++) for(long ib=0;ib<nc;ib++){ u64 k=ak[ia]+ck[ib]; u64 hh=(k*0x9E3779B97F4A7C15ULL)>>20; hh&=tsz-1; while(tab[hh].used) hh=(hh+1)&(tsz-1); tab[hh].key=k;tab[hh].ia=(u16)ia;tab[hh].ib=(u32)ib;tab[hh].used=1; }
  long long count=0; long long tested=0,dbl=0; long dstride=argc>3?atol(argv[3]):1; unsigned __int128 tot=0; long printed=0; u64 rs=12345;
  for(long ic=0;ic<nc;ic++) for(long id=0;id<nc;id++){
    u64 k=-(ck[ic]+ck[id]); u64 hh=(k*0x9E3779B97F4A7C15ULL)>>20; hh&=tsz-1;
    while(tab[hh].used){ if(tab[hh].key==k){ count++; if(dstride>0 && (count%dstride)==0){ tested++; if(is_doubled(al[tab[hh].ia],cl[tab[hh].ib],cl[ic],cl[id])) dbl++; } tot+=(unsigned __int128)osz[tab[hh].ib]*osz[ic]*osz[id];
        // reservoir-free sampling: print the first maxs and then every so often
        if(printed<maxs){ printf("SOL %u %u %u %u\n",al[tab[hh].ia],cl[tab[hh].ib],cl[ic],cl[id]); printed++; }
        else { rs^=rs<<13;rs^=rs>>7;rs^=rs<<17; if((rs%(u64)(count))<(u64)maxs){ /* replace not possible in stream; skip */ } }
      } hh=(hh+1)&(tsz-1); }
  }
  printf("DOUBLING_TEST tested %lld doubled %lld\n",tested,dbl); printf("TOTAL_ORDERED_FULL %llu\n",(unsigned long long)tot); printf("COUNT %lld  (ordered (A,B,C,D), A palindromic, B,C,D canonical monomial-orbit reps; nc=%ld)\n",count,nc);
  return 0; }
