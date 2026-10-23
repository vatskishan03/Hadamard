// qskew.c: orbit-level exhaustive join for GS quadruples over Z_v (v odd) invariant under a multiplier group H.
// mode 0: REAL skew GS (values +-1, A(0)=1, A(-g)=-A(g)), PAF sum = 0 off 0  (ordered count = 8 * printed count)
// mode 1: QUATERNARY skew GS (values i^e, A twisted-skew u(-g)=i conj(u(g)), u(0) in {1,-i}); ordered count = 64 * printed count
// input file: v n ; sz[n] ; neg[n] ; then for c=1..n-1: n*n matrix N_c(a,b)=#{g in O_a : g+s_c in O_b}
// usage: qskew file mode maxprint [firstonly]
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
typedef unsigned long long u64;
static int v,n,d; static int sz[64],neg_[64]; static int *Nt; // Nt[c*n*n+a*n+b]
static u64 P[128];
typedef struct { u64 key; int rre,rim; u64 ex; } El;
static El *L,*A; static long nL,nA;
static int mode;
static inline int mod4(int x){return x&3;}
// compute PAF vector (re,im per class pair) and key
static void paf(const int *e, int *out /* size 2*d */){
  for(int c=1;c<n;c++){
    int M[4]={0,0,0,0};
    const int *N=Nt+(size_t)c*n*n;
    for(int a=0;a<n;a++) for(int b=0;b<n;b++){ int x=N[a*n+b]; if(x) M[mod4(e[a]-e[b])]+=x; }
    out[2*(c-1)]=M[0]-M[2]; out[2*(c-1)+1]=M[1]-M[3];
  }
}
static u64 keyof(const int *vec){ u64 k=0; for(int i=0;i<2*(n-1);i++) k+=(u64)(long long)vec[i]*P[i]; return k; }
static void fill(El *E,const int *e){
  int vec[256]; paf(e,vec); E->key=keyof(vec);
  int rr=0,ri=0; u64 ex=0;
  for(int a=0;a<n;a++){ int t=mod4(e[a]); ex|=((u64)t)<<(2*a); switch(t){case 0:rr+=sz[a];break;case 1:ri+=sz[a];break;case 2:rr-=sz[a];break;case 3:ri-=sz[a];break;} }
  E->rre=rr;E->rim=ri;E->ex=ex;
}
static int cmpR(const void*x,const void*y){const El*a=x,*b=y; if(a->rre!=b->rre) return a->rre<b->rre?-1:1; if(a->rim!=b->rim) return a->rim<b->rim?-1:1; return 0;}
typedef struct {int rre,rim; long lo,hi;} Bk;
static Bk *bL,*bA; static long nbL,nbA;
static Bk* buckets(El*E,long m,long*nb){ Bk*b=malloc(sizeof(Bk)*(m+1)); long k=0; for(long i=0;i<m;){ long j=i; while(j<m&&E[j].rre==E[i].rre&&E[j].rim==E[i].rim) j++; b[k].rre=E[i].rre;b[k].rim=E[i].rim;b[k].lo=i;b[k].hi=j;k++; i=j;} *nb=k; return b;}
// hash table for left pairs
typedef struct { u64 key; unsigned int ia, ib; } Ent;
int main(int argc,char**argv){
  FILE*f=fopen(argv[1],"r"); mode=atoi(argv[2]); int maxprint=argc>3?atoi(argv[3]):5; int firstonly=argc>4?atoi(argv[4]):0; long stride=argc>5?atol(argv[5]):1;
  if(fscanf(f,"%d %d",&v,&n)!=2) return 1; d=n-1;
  for(int i=0;i<n;i++) if(fscanf(f,"%d",&sz[i])!=1) return 1;
  for(int i=0;i<n;i++) if(fscanf(f,"%d",&neg_[i])!=1) return 1;
  Nt=malloc(sizeof(int)*(size_t)n*n*n);
  for(int c=1;c<n;c++) for(int a=0;a<n;a++) for(int b=0;b<n;b++) if(fscanf(f,"%d",&Nt[(size_t)c*n*n+a*n+b])!=1) return 1;
  fclose(f);
  u64 s=88172645463325252ULL; for(int i=0;i<128;i++){ s^=s<<13;s^=s>>7;s^=s<<17; P[i]=s|1ULL; }
  // list L for B,C,D: w_0 = 1 (e0=0); values: mode0 e in {0,2}, mode1 e in 0..3
  int base=mode?4:2, step=mode?1:2;
  nL=1; for(int i=1;i<n;i++) nL*=base;
  L=malloc(sizeof(El)*nL);
  int e[64];
  for(long idx=0;idx<nL;idx++){ long t=idx; e[0]=0; for(int i=1;i<n;i++){ e[i]=(t%base)*step; t/=base; } fill(&L[idx],e); }
  // skew list A
  int pr[64],np=0; for(int a=1;a<n;a++) if(a<neg_[a]) pr[np++]=a;
  int a0opts=2; nA=a0opts; for(int i=0;i<np;i++) nA*=base;
  A=malloc(sizeof(El)*nA);
  for(long idx=0;idx<nA;idx++){ long t=idx; int o0=t%2; t/=2; e[0]= mode? (o0?3:0) : 0; if(!mode && o0) { /* real: only o0=0 valid */ }
    for(int i=0;i<np;i++){ int ej=(t%base)*step; t/=base; e[pr[i]]=ej; e[neg_[pr[i]]]= mode? mod4(1-ej) : mod4(ej+2); }
    e[0]= mode ? (o0?3:0) : 0;
    fill(&A[idx],e); if(!mode && o0) A[idx].key=0, A[idx].rre=1<<30; }
  if(!mode){ long k=0; for(long i=0;i<nA;i++) if(A[i].rre!=(1<<30)) A[k++]=A[i]; nA=k; }
  fprintf(stderr,"n=%d d=%d |L|=%ld |A|=%ld\n",n,d,nL,nA);
  qsort(L,nL,sizeof(El),cmpR); qsort(A,nA,sizeof(El),cmpR);
  bL=buckets(L,nL,&nbL); bA=buckets(A,nA,&nbA);
  long total=4L*v; long long count=0; int printed=0; long incomplete=0;
  // loop over (RA,RB) build left table, then (RC,RD)
  for(long ia=0;ia<nbA;ia++) for(long ib=0;ib<nbL;ib++){
    long na_=bA[ia].hi-bA[ia].lo, nb_=bL[ib].hi-bL[ib].lo; long nl=na_*nb_;
    long nA2=(long)bA[ia].rre*bA[ia].rre+(long)bA[ia].rim*bA[ia].rim, nB2=(long)bL[ib].rre*bL[ib].rre+(long)bL[ib].rim*bL[ib].rim;
    if(nA2+nB2>total) continue;
    { static long fc=0; fc++; if(stride>1 && (fc%stride)!=0) continue; }
    // chunk the A bucket so that the left table has at most 4e7 entries
    long chunk=na_; if(nl>40000000L){ chunk=40000000L/nb_; if(chunk<1) chunk=1; }
    for(long alo=bA[ia].lo; alo<bA[ia].hi; alo+=chunk){
    long ahi=alo+chunk; if(ahi>bA[ia].hi) ahi=bA[ia].hi;
    long nlc=(ahi-alo)*nb_;
    long tsz=1; while(tsz<2*nlc) tsz<<=1;
    int built=0; Ent*tab=NULL;
    for(long ic=0;ic<nbL;ic++){
      long nC2=(long)bL[ic].rre*bL[ic].rre+(long)bL[ic].rim*bL[ic].rim; if(nA2+nB2+nC2>total) continue;
      for(long id=0;id<nbL;id++){
        long nD2=(long)bL[id].rre*bL[id].rre+(long)bL[id].rim*bL[id].rim; if(nA2+nB2+nC2+nD2!=total) continue;
        if(!built){ tab=calloc(tsz,sizeof(Ent)); if(!tab){ fprintf(stderr,"calloc fail\n"); incomplete++; break; } for(long x=alo;x<ahi;x++) for(long y=bL[ib].lo;y<bL[ib].hi;y++){ u64 k=A[x].key+L[y].key; u64 h=(k*0x9E3779B97F4A7C15ULL)>>20; h&=tsz-1; while(tab[h].ia|tab[h].ib?1:0){ h=(h+1)&(tsz-1);} tab[h].key=k;tab[h].ia=(unsigned)x+1;tab[h].ib=(unsigned)y+1; } built=1; }
        for(long x=bL[ic].lo;x<bL[ic].hi;x++) for(long y=bL[id].lo;y<bL[id].hi;y++){
          u64 k=-(L[x].key+L[y].key); u64 h=(k*0x9E3779B97F4A7C15ULL)>>20; h&=tsz-1;
          while(tab[h].ia|tab[h].ib){ if(tab[h].key==k){ if(printed>=maxprint){ count++; h=(h+1)&(tsz-1); continue; } // verify exactly
              int e1[64],e2[64],e3[64],e4[64]; for(int i=0;i<n;i++){ e1[i]=(A[tab[h].ia-1].ex>>(2*i))&3; e2[i]=(L[tab[h].ib-1].ex>>(2*i))&3; e3[i]=(L[x].ex>>(2*i))&3; e4[i]=(L[y].ex>>(2*i))&3; }
              int v1[256],v2[256],v3[256],v4[256]; paf(e1,v1);paf(e2,v2);paf(e3,v3);paf(e4,v4); int ok=1; for(int i=0;i<2*d;i++) if(v1[i]+v2[i]+v3[i]+v4[i]){ok=0;break;}
              if(ok){ count++; if(printed<maxprint){ printed++; printf("SOL"); for(int t=0;t<4;t++){ int*ee=t==0?e1:t==1?e2:t==2?e3:e4; printf(" "); for(int i=0;i<n;i++) printf("%d",ee[i]); } printf("\n"); fflush(stdout);} if(firstonly&&printed>=maxprint){ printf("COUNT_AT_STOP %lld\n",count); return 0;} }
            }
            h=(h+1)&(tsz-1); }
        }
      }
    }
    if(built) free(tab);
    }
  }
  printf("INCOMPLETE_TUPLES %ld\n",incomplete); printf("COUNT %lld (normalized: B,C,D have w0=1; multiply by %d for ordered; stride %ld)\n",count, mode?64:8, stride);
  return 0;
}
