// oskew.c: orbit-level exhaustive join for OCTIC skew Goethals-Seidel quadruples over Z_v (v odd), alphabet A = {c0+c1 w+c2 w^2+c3 w^3, c_i=+-1},
// blocks constant on orbits of a multiplier group H (-1 not in H).  Ring Z[w][Z_v] = Z[x]/(x^{4v}+1);  sum_X U_X U_X^* = 16 v delta.
// A skew: U(0) = (1,c,d,c), U(-g) = -conj U(g).  B,C,D normalised by w^k: U(0) in {(1,1,1,1),(1,-1,1,1)}.
// code e in 0..15: bit j set <=> c_j = -1.
// input: v n ; sz[n] ; neg[n] ; then for c=1..n-1 the n*n matrix N_c(a,b) = #{g in O_a : g + s_c in O_b}
// usage: oskew file maxprint [firstonly]
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
typedef uint64_t u64; typedef int64_t i64;
static int v,n,ncls; static int sz[32],neg_[32]; static int *Nt; static int cls[32];
static int cb[16][4];                  // coefficients
static int Tp[16][16][4];              // c_p * conj(c_q)
static u64 Rw[32][4]; static u64 *W;   // W[((a*n+b)*16+x)*16+y]
#define WI(a,b,x,y) W[((((size_t)(a)*n+(b))*16+(x))*16+(y))]
static void zmul(const int*a,const int*b,int*o){ int r[4]={0,0,0,0}; for(int i=0;i<4;i++)for(int j=0;j<4;j++){ int k=i+j; if(k>=4) r[k-4]-=a[i]*b[j]; else r[k]+=a[i]*b[j]; } for(int i=0;i<4;i++) o[i]=r[i]; }
static void zconj(const int*a,int*o){ o[0]=a[0]; o[1]=-a[3]; o[2]=-a[2]; o[3]=-a[1]; }
typedef struct { u64 key; u64 pk; } El;       // pk = dcid<<40 | ex  (ex: 4 bits per orbit, n<=10)
static El *L,*A; static long nL,nA;
// DC classes
typedef struct { int R[4]; int P,Q; } Dc; static Dc *dc; static int ndc=0; static u64 *dctab; static int *dcidx; static long dcsz=1<<22;
static int dc_id(const int*R){ u64 p=0; for(int j=0;j<4;j++) p|=((u64)(R[j]+4096))<<(13*j); u64 h=(p*0x9E3779B97F4A7C15ULL)>>40; h&=dcsz-1;
  while(dctab[h]){ if(dctab[h]-1==p) return dcidx[h]; h=(h+1)&(dcsz-1);} 
  dctab[h]=p+1; dcidx[h]=ndc; for(int j=0;j<4;j++) dc[ndc].R[j]=R[j]; dc[ndc].P=R[0]*R[0]+R[1]*R[1]+R[2]*R[2]+R[3]*R[3]; dc[ndc].Q=R[0]*R[1]+R[1]*R[2]+R[2]*R[3]-R[0]*R[3]; return ndc++; }
static void paf_vec(const int*e,int*out /*ncls*4*/){ for(int ci=0;ci<ncls;ci++){ int c=cls[ci]; int M[4]={0,0,0,0}; const int*N=Nt+(size_t)c*n*n;
    for(int a=0;a<n;a++)for(int b=0;b<n;b++){ int x=N[a*n+b]; if(x) for(int j=0;j<4;j++) M[j]+=x*Tp[e[b]][e[a]][j]; } for(int j=0;j<4;j++) out[ci*4+j]=M[j]; } }
static u64 keyof(const int*e){ u64 k=0; for(int a=0;a<n;a++)for(int b=0;b<n;b++) k+=WI(a,b,e[a],e[b]); return k; }
static int gen_mode; static long gpos; static int ee[32];
static void leaf(El*E){ u64 ex=0; int R[4]={0,0,0,0}; for(int a=0;a<n;a++){ ex|=((u64)ee[a])<<(4*a); for(int j=0;j<4;j++) R[j]+=sz[a]*cb[ee[a]][j]; }
  E->key=keyof(ee); E->pk=((u64)dc_id(R)<<40)|ex; }
static void genL(int a){ if(a==n){ leaf(&L[gpos++]); return; } if(a==0){ int opts[2]={0,2}; for(int t=0;t<2;t++){ ee[0]=opts[t]; genL(1);} return; } for(int e=0;e<16;e++){ ee[a]=e; genL(a+1);} }
static int pr[32],np_;
static int negcode(int e){ int b0=e&1,b1=(e>>1)&1,b2=(e>>2)&1,b3=(e>>3)&1; return (b0^1)|(b3<<1)|(b2<<2)|(b1<<3); }
static void genA(int i){ if(i==np_){ leaf(&A[gpos++]); return; } for(int e=0;e<16;e++){ ee[pr[i]]=e; ee[neg_[pr[i]]]=negcode(e); genA(i+1);} }
static int cmpPk(const void*x,const void*y){ u64 a=((const El*)x)->pk,b=((const El*)y)->pk; return a<b?-1:a>b; }
typedef struct { int id; long lo,hi; } Bk;
static Bk*gb;
static int cmpPQ(const void*x,const void*y){ int a=*(const int*)x,b=*(const int*)y; Dc*da=&dc[gb[a].id],*db=&dc[gb[b].id]; if(da->P!=db->P) return da->P<db->P?-1:1; if(da->Q!=db->Q) return da->Q<db->Q?-1:1; return a<b?-1:a>b; }
static Bk* buckets(El*E,long m,int*nb){ Bk*b=malloc(sizeof(Bk)*(m+1)); int k=0; for(long i=0;i<m;){ long j=i; u64 id=E[i].pk>>40; while(j<m&&(E[j].pk>>40)==id) j++; b[k].id=(int)id;b[k].lo=i;b[k].hi=j;k++; i=j;} *nb=k; return b; }
typedef struct { u64 key; unsigned int ia, ib; } Ent;
int main(int argc,char**argv){
  FILE*f=fopen(argv[1],"r"); int maxprint=argc>2?atoi(argv[2]):3; int firstonly=argc>3?atoi(argv[3]):0;
  if(fscanf(f,"%d %d",&v,&n)!=2) return 1; for(int i=0;i<n;i++) if(fscanf(f,"%d",&sz[i])!=1) return 1; for(int i=0;i<n;i++) if(fscanf(f,"%d",&neg_[i])!=1) return 1;
  Nt=malloc(sizeof(int)*(size_t)n*n*n); for(int c=1;c<n;c++)for(int a=0;a<n;a++)for(int b=0;b<n;b++) if(fscanf(f,"%d",&Nt[(size_t)c*n*n+a*n+b])!=1) return 1; fclose(f);
  for(int a=1;a<n;a++) if(neg_[a]==a){ printf("DEAD_SELFNEG\nCOUNT 0\n"); return 0; }
  for(int e=0;e<16;e++)for(int j=0;j<4;j++) cb[e][j]=((e>>j)&1)?-1:1;
  for(int p=0;p<16;p++)for(int q=0;q<16;q++){ int cq[4]; zconj(cb[q],cq); zmul(cb[p],cq,Tp[p][q]); }
  ncls=0; for(int c=1;c<n;c++) if(c<neg_[c]) cls[ncls++]=c;
  u64 s=88172645463325252ULL; for(int c=0;c<32;c++)for(int j=0;j<4;j++){ s^=s<<13;s^=s>>7;s^=s<<17; Rw[c][j]=s|1ULL; }
  W=calloc((size_t)n*n*256,sizeof(u64));
  // PAF_c = sum_{a,b} N_c(a,b) U_b conj(U_a) = sum N_c(a,b) Tp[e_b][e_a]
  for(int a=0;a<n;a++)for(int b=0;b<n;b++)for(int x=0;x<16;x++)for(int y=0;y<16;y++){ u64 acc=0; for(int ci=0;ci<ncls;ci++){ int c=cls[ci]; int N=Nt[(size_t)c*n*n+a*n+b]; if(!N) continue; for(int j=0;j<4;j++) acc+=(u64)(i64)(N*Tp[y][x][j])*Rw[ci][j]; } WI(a,b,x,y)=acc; }
  dc=malloc(sizeof(Dc)*(1<<21)); dctab=calloc(dcsz,sizeof(u64)); dcidx=malloc(sizeof(int)*dcsz);
  nL=2; for(int i=1;i<n;i++) nL*=16; L=malloc(sizeof(El)*nL); gpos=0; genL(0);
  np_=0; for(int a=1;a<n;a++) if(a<neg_[a]) pr[np_++]=a; nA=4; for(int i=0;i<np_;i++) nA*=16; A=malloc(sizeof(El)*nA); gpos=0;
  for(int u=0;u<4;u++){ int c=u&1,d=(u>>1)&1; ee[0]=(c?(2|8):0)|(d?4:0); // U(0)=(1,c,d,c): bits: c1=bit1,c2=bit2,c3=bit3
    long save=gpos; (void)save; genA(0); }
  if(gpos!=nA){ fprintf(stderr,"nA mismatch %ld %ld\n",gpos,nA); return 1; }
  fprintf(stderr,"n=%d v=%d |L|=%ld |A|=%ld ndc=%d ncls=%d\n",n,v,nL,nA,ndc,ncls);
  qsort(L,nL,sizeof(El),cmpPk); qsort(A,nA,sizeof(El),cmpPk);
  int nbL,nbA; Bk*bL=buckets(L,nL,&nbL); Bk*bA=buckets(A,nA,&nbA);
  long T=16L*v; long long count=0, countw=0; int printed=0; long incomplete=0; long long lookups=0, tuples=0;
  // PQ map for L classes: sort L-class indices by (P,Q)
  int*ord=malloc(sizeof(int)*nbL); for(int i=0;i<nbL;i++) ord[i]=i;
  gb=bL;      // sort the class indices by (P,Q), P + sqrt2 Q = |R|^2 of the class
  qsort(ord,nbL,sizeof(int),cmpPQ);
  int*posPQ=ord;
  for(int ia=0;ia<nbA;ia++){ Dc*da=&dc[bA[ia].id]; for(int ib=0;ib<nbL;ib++){ Dc*db=&dc[bL[ib].id]; if(da->P+db->P>T) continue;
    // find any valid (cc>=ib,cd>=cc) quickly: loop cc
    long na_=bA[ia].hi-bA[ia].lo, nb_=bL[ib].hi-bL[ib].lo; long nl=na_*nb_; long chunk=na_; if(nl>30000000L){ chunk=30000000L/nb_; if(chunk<1) chunk=1; }
    for(long alo=bA[ia].lo; alo<bA[ia].hi; alo+=chunk){ long ahi=alo+chunk; if(ahi>bA[ia].hi) ahi=bA[ia].hi; long nlc=(ahi-alo)*nb_; long tsz=1; while(tsz<2*nlc) tsz<<=1; Ent*tab=NULL; int built=0;
      for(int ic=ib;ic<nbL;ic++){ Dc*dcc=&dc[bL[ic].id]; long Pc=T-da->P-db->P-dcc->P; long Qc=-(da->Q+db->Q+dcc->Q); if(Pc<0) continue;
        // binary search in posPQ for (Pc,Qc)
        long lo=0,hi=nbL; while(lo<hi){ long mid=(lo+hi)/2; Dc*dm=&dc[bL[posPQ[mid]].id]; if(dm->P<Pc||(dm->P==Pc&&dm->Q<Qc)) lo=mid+1; else hi=mid; }
        for(long t=lo;t<nbL;t++){ int id=posPQ[t]; Dc*dd=&dc[bL[id].id]; if(dd->P!=Pc||dd->Q!=Qc) break; if(id<ic) continue;
          int mult = (ib==ic&&ic==id)?1: (ib==ic||ic==id)?3:6; tuples++;
          if(!built){ tab=calloc(tsz,sizeof(Ent)); if(!tab){ incomplete++; break; } for(long x=alo;x<ahi;x++)for(long y=bL[ib].lo;y<bL[ib].hi;y++){ u64 k=A[x].key+L[y].key; u64 h=(k*0x9E3779B97F4A7C15ULL)>>20; h&=tsz-1; while(tab[h].ia|tab[h].ib) h=(h+1)&(tsz-1); tab[h].key=k; tab[h].ia=(unsigned)x+1; tab[h].ib=(unsigned)y+1; } built=1; }
          for(long x=bL[ic].lo;x<bL[ic].hi;x++)for(long y=bL[id].lo;y<bL[id].hi;y++){ u64 k=-(L[x].key+L[y].key); u64 h=(k*0x9E3779B97F4A7C15ULL)>>20; h&=tsz-1; lookups++;
            while(tab[h].ia|tab[h].ib){ if(tab[h].key==k){ // verify exactly
                int e1[32],e2[32],e3[32],e4[32]; for(int i=0;i<n;i++){ e1[i]=(A[tab[h].ia-1].pk>>(4*i))&15; e2[i]=(L[tab[h].ib-1].pk>>(4*i))&15; e3[i]=(L[x].pk>>(4*i))&15; e4[i]=(L[y].pk>>(4*i))&15; }
                int v1[512],v2[512],v3[512],v4[512]; paf_vec(e1,v1);paf_vec(e2,v2);paf_vec(e3,v3);paf_vec(e4,v4); int ok=1; for(int i=0;i<4*ncls;i++) if(v1[i]+v2[i]+v3[i]+v4[i]){ok=0;break;}
                if(ok){ count++; countw+=mult; if(printed<maxprint){ printed++; printf("SOL"); for(int t2=0;t2<4;t2++){ int*q=t2==0?e1:t2==1?e2:t2==2?e3:e4; printf(" "); for(int i=0;i<n;i++) printf("%x",q[i]); } printf("\n"); fflush(stdout); if(firstonly&&printed>=maxprint){ printf("COUNT_AT_STOP %lld\n",count); return 0; } } } }
              h=(h+1)&(tsz-1); } }
        }
      }
      if(built) free(tab);
    } } }
  printf("INCOMPLETE_TUPLES %ld\nDC_TUPLES %lld LOOKUPS %lld\n",incomplete,tuples,lookups);
  printf("COUNT %lld  (class-ordered B<=C<=D; ORDERED total of (A,B,C,D) with B,C,D normalised by w^k = %lld; times 8^3 for all ordered)\n",count,countw);
  return 0; }
