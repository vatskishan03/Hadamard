// tri2.c : cache-resident partitioned four-list join for the triangular family v = 1 + k(k+1)/2 prime, k = 3 mod 4 (search/NOTES.md section 88 and section 92; README.md of this directory).
// H <= Z_v* of order k, e = (v-1)/k cosets C_j = g^j H (orbit 0 = {0}, orbit j+1 = C_j); the quaternary blocks are constant on the cosets, values i^x; the four complex PAFs must sum to zero at every shift != 0.
// Compared with tri_join.c (modulus a power of two, 64-bit random linear key): the partition is by a residue modulo a PRIME M, the left table of one partition is sized to the pairs of that partition and stays in the cache
// (time per probe: README.md), every hit is verified in integer arithmetic, and red=1 (below) divides the work by the symmetry group of PROOFS.md P5.
// mode 3 (k = 3 mod 8): A = twisted-skew BIG block (u0 in {1,-i}, S_A = m(1+i), m = +-1 aligned), B,C,D free small (u0 = 1, S = 0).      left pairs (A,B), right pairs (C,D)
// mode 7 (k = 7 mod 8): A = twisted-skew SMALL block (S_A = 0), B = free BIG (u0 = 1, S = 1+i), C,D free small.  left pairs (A,B), right pairs (C,D)
// (the row-sum pattern 4v-3, 1, 1, 1 is PROVED in PROOFS.md, P2; it assigns the big block.)
// red=1 (modes 3 and 7): A over representatives of the Z_e shift orbits (weight = orbit size), C,D over tilde-class representatives (tilde: c_j -> conj c_{j+e/2}),
//   unordered {C,D}; B only S = 1+i (tilde gives S = 1-i).  Total normalised count = (position of the big block: 3) * 2 * sum of weights  (in mode 7);
//   in mode 3 and in mode 7 with red=0 every solution (ordered A,B,C,D as in qskew's normalised count) is counted once with weight 1
//   (mode 7, red=0: counts B big with S = 1+i only; total = 3 * 2 * count).
// The join is partitioned by a linear residue pi(PAF) mod M (M prime): left pairs with pi = r, right pairs with pi = T - r.   Every solution lies in exactly one partition,
//   and a solution without a stabiliser comes with 31 others in the same partition (PROOFS.md, P5: tilde of B, C, D, swap of C and D, A -> -A), so independent hits are clusters of 32 (16 in mode 7, where the tilde of the big block B leaves the list).
// Limits: the digits of a block are packed 2 bits each into 32 bits, so n = e + 1 <= 16 (e <= 15: v = 7, 29, 67, 191, 277, 379); the program refuses larger e (a guard added when the code entered the repository, no computation changed).
// usage: tri2 v k g mode red M r0 step cnt idx0 maxsol seed      (g a primitive root mod v)
//   partitions i = idx0..cnt-1 with r_i = (r0 + i*step) mod M   (resumable by idx0; with step = 1, r0 = 0, cnt = M all partitions are processed)
//   output: SOL lines "SOL i=I r=R w=W A B C D" (digit strings u0, c_0..c_{e-1} in coset order; w = the weight of a reduced hit), one STAT line per partition, a TOTAL line.
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
typedef unsigned long long u64; typedef long long i64; typedef unsigned int u32;
static int v,k,e,n,np_,h,mode,red;
static int *orb; static int *Ncnt;
static u64 KR[128]; static i64 QR[128]; static i64 M;
typedef struct { u64 key; u32 pi; u32 ex; } El;
static int paf(const int *x,int *out){
  for(int t=0;t<np_;t++){ int Mm[4]={0,0,0,0}; const int *N=Ncnt+(size_t)t*n*n;
    for(int a=0;a<n;a++){ const int *Na=N+a*n; int xa=x[a]; for(int b=0;b<n;b++){ int c=Na[b]; if(c) Mm[(xa-x[b])&3]+=c; } }
    out[2*t]=Mm[0]-Mm[2]; out[2*t+1]=Mm[1]-Mm[3]; }
  return 0;
}
static u64 TK; static i64 TPI;
static void fill(El *E,const int *x){
  int p[128]; paf(x,p); u64 kk=0; i64 q=0;
  for(int t=0;t<np_;t++){
    if(!(p[2*t]&1)||(p[2*t+1]&1)){ fprintf(stderr,"parity assumption violated\n"); exit(2); }
    i64 a=(p[2*t]-1)/2, b=p[2*t+1]/2;
    kk+=(u64)a*KR[2*t]+(u64)b*KR[2*t+1]; q+=a*QR[2*t]+b*QR[2*t+1]; }
  q%=M; if(q<0) q+=M;
  E->key=kk; E->pi=(u32)q; u64 ex=0; for(int a=0;a<n;a++) ex|=((u64)(x[a]&3))<<(2*a); E->ex=(u32)ex;
}
static void unpack(u32 ex,int *x){ for(int a=0;a<n;a++) x[a]=(ex>>(2*a))&3; }
static u32 tilde_ex(u32 ex){ int x[64],y[64]; unpack(ex,x); y[0]=0; for(int j=0;j<e;j++) y[1+j]=(-x[1+((j+h)%e)])&3; u32 r=0; for(int a=0;a<n;a++) r|=((u32)y[a]&3)<<(2*a); return r; }
static u32 shift_ex(u32 ex,int s){ int x[64],y[64]; unpack(ex,x); y[0]=x[0]; for(int j=0;j<e;j++) y[1+j]=x[1+((j+s)%e)]; u32 r=0; for(int a=0;a<n;a++) r|=((u32)y[a]&3)<<(2*a); return r; }
// ---------- list generators
static El *Lst; static long nL=0, capL;
static void genSmall(int pos,int *x,int c0,int c1,int c2,int c3){ // u0 = 1, coset sum 0
  if(pos==n){ if(c0==c2&&c1==c3) fill(&Lst[nL++],x); return; }
  int rem=n-pos;
  for(int t=0;t<4;t++){ int d0=c0+(t==0),d1=c1+(t==1),d2=c2+(t==2),d3=c3+(t==3);
    if(abs(d0-d2)+abs(d1-d3)>rem-1) continue; x[pos]=t; genSmall(pos+1,x,d0,d1,d2,d3); }
}
static El *Big; static long nBig=0;
static void genBig(int pos,int *x,int c0,int c1,int c2,int c3){ // u0 = 1, coset sum = 1+i
  if(pos==n){ if(c0-c2==1&&c1-c3==1) fill(&Big[nBig++],x); return; }
  int rem=n-pos;
  for(int t=0;t<4;t++){ int d0=c0+(t==0),d1=c1+(t==1),d2=c2+(t==2),d3=c3+(t==3);
    if(abs(d0-d2-1)+abs(d1-d3-1)>rem-1) continue; x[pos]=t; genBig(pos+1,x,d0,d1,d2,d3); }
}
static El *Alst; static long nA=0;
static double binom(int a,int b){ double c=1; for(int i=1;i<=b;i++) c=c*(a-b+i)/i; return c; }
static void genA(int bigflag){
  // pairs (j, j+h), j<h: sign + (c in {0,1}) or - (c in {2,3}); partner (1-c)&3.  big: m=+1 with u0=1, m=-1 with u0=-i (exp 3);  small: sum of signs 0, u0 either.
  long cap=2; for(int i=0;i<h;i++) cap*=4; Alst=malloc(sizeof(El)*cap);
  int x[64];
  for(int o0=0;o0<2;o0++){
    int need = bigflag ? (o0?-1:1) : 0;   // (#plus - #minus)
    for(long mask=0;mask<(1L<<h);mask++){ int pl=__builtin_popcountl(mask); if(2*pl-h!=need) continue;
      for(long chs=0;chs<(1L<<h);chs++){
        for(int j=0;j<h;j++){ int plus=(mask>>j)&1; int c=((chs>>j)&1)+(plus?0:2); x[1+j]=c; x[1+j+h]=(1-c)&3; }
        x[0]=o0?3:0; fill(&Alst[nA++],x);
      } } }
}
// ---------- sorting by pi
static int cmpPi(const void*a,const void*b){ const El*x=a,*y=b; if(x->pi!=y->pi) return x->pi<y->pi?-1:1; return x->ex<y->ex?-1:(x->ex>y->ex); }
typedef struct { u64 *key; u32 *pi; u32 *ex; long *off; long cnt; } Sorted;
static void makeSorted(Sorted *S,El *E,long cnt){
  qsort(E,cnt,sizeof(El),cmpPi);
  S->cnt=cnt; S->key=malloc(8*cnt); S->pi=malloc(4*cnt); S->ex=malloc(4*cnt); S->off=calloc(M+2,sizeof(long));
  for(long i=0;i<cnt;i++){ S->key[i]=E[i].key; S->pi[i]=E[i].pi; S->ex[i]=E[i].ex; S->off[E[i].pi+1]++; }
  for(i64 r=0;r<M;r++) S->off[r+1]+=S->off[r];
}
int main(int argc,char**argv){
  if(argc<13){ fprintf(stderr,"usage: tri2 v k g mode red M r0 step cnt idx0 maxsol seed\n"); return 1; }
  v=atoi(argv[1]); k=atoi(argv[2]); int g=atoi(argv[3]); mode=atoi(argv[4]); red=atoi(argv[5]); M=atol(argv[6]);
  long r0=atol(argv[7]), step=atol(argv[8]), cnt=atol(argv[9]), idx0=atol(argv[10]), maxsol=atol(argv[11]); u64 seed=strtoull(argv[12],0,10);
  e=(v-1)/k; n=e+1; np_=e/2; h=e/2;
  if(n>16){ fprintf(stderr,"e = %d: the digits of a block no longer fit into 32 bits (n = e+1 <= 16)\n",e); return 1; }
  if(mode==3 && (k%8)!=3){ fprintf(stderr,"mode 3 needs k=3 mod 8\n"); return 1; }
  if(mode==7 && (k%8)!=7){ fprintf(stderr,"mode 7 needs k=7 mod 8\n"); return 1; }
  orb=malloc(sizeof(int)*v); orb[0]=0;
  long hgen=1; for(int i=0;i<e;i++) hgen=hgen*g%v;
  long gj=1; for(int j=0;j<e;j++){ long y=gj; for(int i=0;i<k;i++){ orb[y]=j+1; y=y*hgen%v; } gj=gj*g%v; }
  for(int x=1;x<v;x++) if(orb[x]<1){ fprintf(stderr,"orbit error\n"); return 1; }
  Ncnt=calloc((size_t)np_*n*n,sizeof(int));
  gj=1; for(int t=0;t<np_;t++){ for(int x=0;x<v;x++) Ncnt[(size_t)t*n*n+orb[x]*n+orb[(x+gj)%v]]++; gj=gj*g%v; }
  u64 s=seed*0x9E3779B97F4A7C15ULL+88172645463325252ULL; for(int i=0;i<128;i++){ s^=s<<13;s^=s>>7;s^=s<<17; KR[i]=s|1ULL; s^=s<<13;s^=s>>7;s^=s<<17; QR[i]=(i64)(s%(u64)M); }
  TK=0; TPI=0; for(int t=0;t<np_;t++){ TK+=(u64)(-2LL)*KR[2*t]; TPI+=-2*QR[2*t]; } TPI%=M; if(TPI<0) TPI+=M;
  int x[64];
  // lists
  long expL=1; { double c=binom(e,h); expL=(long)(c*c)+16; }
  Lst=malloc(sizeof(El)*expL); x[0]=0; genSmall(1,x,0,0,0,0);
  fprintf(stderr,"v=%d k=%d e=%d mode=%d red=%d M=%ld |L|=%ld\n",v,k,e,mode,red,(long)M,nL);
  // reduced small list
  if(red){ long m2=0; for(long i=0;i<nL;i++){ u32 t=tilde_ex(Lst[i].ex); if(Lst[i].ex<=t) Lst[m2++]=Lst[i]; } nL=m2; fprintf(stderr,"|L'| (tilde reps) = %ld\n",nL); }
  // right list R = small
  static Sorted R,XB; makeSorted(&R,Lst,nL);
  // A list
  genA(mode==3); fprintf(stderr,"|A|=%ld\n",nA);
  long *wA=NULL;
  if(red){ // shift-orbit representatives
    long m2=0; wA=malloc(sizeof(long)*nA);
    for(long i=0;i<nA;i++){ u32 ex=Alst[i].ex; u32 best=ex; int sz=0; u32 seen[64]; for(int sft=0;sft<e;sft++){ u32 y=shift_ex(ex,sft); int dup=0; for(int q=0;q<sz;q++) if(seen[q]==y){dup=1;break;} if(!dup) seen[sz++]=y; if(y<best) best=y; }
      if(best==ex){ wA[m2]=sz; Alst[m2++]=Alst[i]; } }
    nA=m2; fprintf(stderr,"A shift-orbit reps = %ld\n",nA);
  }
  // XB: second left list
  Sorted *XBp;
  if(mode==3){ XBp=&R; }
  else { Big=malloc(sizeof(El)*((size_t)(binom(e,h)*binom(e,h+1)*1.05)+1024)); x[0]=0; genBig(1,x,0,0,0,0); fprintf(stderr,"|Big(S=1+i)|=%ld\n",nBig); makeSorted(&XB,Big,nBig); XBp=&XB; }
  // pi histogram of R
  { long mx=0; for(i64 r=0;r<M;r++){ long c=R.off[r+1]-R.off[r]; if(c>mx) mx=c; } fprintf(stderr,"R buckets: mean %.2f max %ld\n",(double)nL/M,mx); }
  // A sorted by pi isn't needed (loop over all A)
  size_t maxleft=0; // allocate tables lazily
  u64 *tkey=NULL; u32 *tia=NULL,*tib=NULL; size_t tcap=0; u64 *bm=NULL; size_t bmcap=0;
  long long tot_probes=0, tot_hits=0; double wsum=0;
  for(long idx=idx0; idx<cnt; idx++){
    i64 r=( (i64)r0 + (i64)idx*step ) % M; i64 rho=(TPI - r)%M; if(rho<0) rho+=M;
    // count left
    long nl=0; for(long ia=0;ia<nA;ia++){ i64 b=(r-(i64)Alst[ia].pi)%M; if(b<0) b+=M; nl+=XBp->off[b+1]-XBp->off[b]; }
    size_t tsz=1; while(tsz<(size_t)(2*nl+16)) tsz<<=1;
    if(tsz>tcap){ free(tkey);free(tia);free(tib); tcap=tsz; tkey=malloc(8*tsz); tia=malloc(4*tsz); tib=malloc(4*tsz); }
    memset(tia,0,4*tsz);
    int bbits=1; while(((size_t)1<<bbits)<(size_t)(16*nl+64)) bbits++; size_t bmsz=((size_t)1<<bbits)/64; if(bmsz>bmcap){ free(bm); bmcap=bmsz; bm=malloc(8*bmsz); } memset(bm,0,8*bmsz);
    u64 tmask=tsz-1;
    for(long ia=0;ia<nA;ia++){ i64 b=(r-(i64)Alst[ia].pi)%M; if(b<0) b+=M;
      for(long q=XBp->off[b];q<XBp->off[b+1];q++){ u64 kk=Alst[ia].key+XBp->key[q]; u64 hh=(kk*0x9E3779B97F4A7C15ULL)>>20; hh&=tmask; while(tia[hh]) hh=(hh+1)&tmask; tkey[hh]=kk; tia[hh]=(u32)ia+1; tib[hh]=(u32)q+1;
        u64 hb=(kk*0xC2B2AE3D27D4EB4FULL)>>(64-bbits); bm[hb>>6]|=1ULL<<(hb&63); } }
    long long probes=0; long hits=0; double wpart=0;
    for(long pc=0;pc<R.cnt;pc++){
      i64 pd=rho-(i64)R.pi[pc]; if(pd<0) pd+=M;
      long q0=R.off[pd], q1=R.off[pd+1]; if(red && q0<pc) q0=pc;   // unordered pairs: position of d >= position of c
      u64 kc=TK-R.key[pc];
      for(long q=q0;q<q1;q++){
        u64 kk=kc-R.key[q]; probes++;
        u64 hb=(kk*0xC2B2AE3D27D4EB4FULL)>>(64-bbits);
        if(!((bm[hb>>6]>>(hb&63))&1)) continue;
        u64 hh=(kk*0x9E3779B97F4A7C15ULL)>>20; hh&=tmask;
        while(tia[hh]){ if(tkey[hh]==kk){ // exact verification
            int xs[4][64]; long ia=tia[hh]-1, ib=tib[hh]-1; u32 exs[4]={Alst[ia].ex,XBp->ex[ib],R.ex[pc],R.ex[q]}; int sum[128]; memset(sum,0,sizeof(sum)); int ok=1;
            for(int z=0;z<4;z++){ unpack(exs[z],xs[z]); int p[128]; paf(xs[z],p); for(int i=0;i<2*np_;i++) sum[i]+=p[i]; }
            for(int i=0;i<2*np_;i++) if(sum[i]) ok=0;
            if(ok){ double w=1; if(red){ w=(double)wA[ia]; int sc=(R.ex[pc]==tilde_ex(R.ex[pc]))?1:2, sd=(R.ex[q]==tilde_ex(R.ex[q]))?1:2; w*=sc*sd; if(pc!=q) w*=2; if(mode==3){ u32 eb=XBp->ex[ib]; w*=(eb==tilde_ex(eb))?1:2; } }
              hits++; wpart+=w;
              if(hits+tot_hits<=maxsol){ printf("SOL i=%ld r=%ld w=%.0f",idx,(long)r,w); for(int z=0;z<4;z++){ printf(" "); for(int a=0;a<n;a++) printf("%d",xs[z][a]); } printf("\n"); } }
          } hh=(hh+1)&tmask; }
      }
    }
    tot_probes+=probes; tot_hits+=hits; wsum+=wpart;
    printf("STAT i=%ld r=%ld left=%ld probes=%lld hits=%ld wsum=%.0f\n",idx,(long)r,nl,probes,hits,wpart); fflush(stdout);
  }
  printf("TOTAL parts=%ld probes=%lld hits=%lld wsum=%.0f\n",cnt-idx0,tot_probes,tot_hits,wsum);
  return 0;
}
