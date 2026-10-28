// band.c : QH(498) search, sparse-band exhaustive canonical run.  Structure: G = Z2 x Z3 x Z41, H = <(1,1,10)> of order 5, Q = Z3, T = all
// characters, N = 249, t = 3 (see README.md).  Four +-1 sequences A, B (symmetric), C, D on G with sum PAF = -12 on {b = 0} \ {0}, 0 elsewhere.
//   The Z2 coordinate splits every sequence into p = (x0 + x1)/2 and m = (x0 - x1)/2 on Z123 = Z3 x Z41: p has values in {-1,0,1} with a
//   balanced-ternary slice structure ("P part", 14 PAF columns, listed by pcv()), m is a sign vector on the complement of the support of p
//   ("M part", 13 columns, target 0).  Level 1: P-tuples (A, B, C, D objects with sum of the P class vectors = target), found by a
//   bucketed exact join of the complete lists A (canonical under the multipliers of Z41), B, and the band lists C, D (objects with exactly
//   n_C resp. n_D free signs, CBAND / DBAND); keys are the P columns 6, 7, 8 (buckets s0 = key 0 of A + C, then s1, s2).  Level 2: for every
//   P-tuple with at least NMIN free signs in total, the sign lists of the four M parts (2^n each, Gray order), pruned by the single-list
//   power-spectrum bound PSD_k(chi) <= 249 - (sum of the other lists' minima), then the exact four-list join on the 13 M columns.
//   Symmetry: A canonical under the multiplier group (mod 4 of the primitive root 6 of Z41), B canonical under Stab(A), C and D reflected
//   (z -> -z) to the smaller id, C <-> D ordered; every orbit of tuples is represented (check_sym.py), and a hit stands for 16 sign copies.
//   Every solution is printed as "SOL ida idb idc idd ga gb gc gd" (object ids and Gray indices); verify_sol.py rebuilds and checks it on G.
// build: gcc -O2 -march=native -o band band.c -lm        (needs AVX-512F: immintrin.h, _mm512_*)
// usage: [CBAND=lo:hi] [DBAND=lo:hi] [MAXB=pairs] band NMIN OUTPREFIX [S0LO S0HI]
//   NMIN       gate: P-tuples with fewer than NMIN free signs in total are not completed (41 in all runs)
//   OUTPREFIX  writes OUTPREFIX.sol (solutions), .prog (one line per finished s0 bucket; a restart skips them), .skip (tuples with a block of
//              more than NMAXVAR = 22 free signs, not treated), and prints progress to stderr
//   S0LO S0HI  only the buckets s0 in [S0LO, S0HI] (default all); buckets are taken in the order -2, -3, -1, -4, 0, ...
//   CBAND, DBAND  free-sign counts of C and D (e.g. 16:16 = exactly 16 each; DBAND defaults to CBAND); MAXB  pairs per pass (memory)
// switches for the checks: PLANTM=1 plants a random M completion per P-tuple and counts whether the join recovers it (printed as PLANT ok/fail);
//   DUMPT=1 prints every P-tuple (PTUPLE lines); NOCANA, NOCANB, NOREFL, NODEDUP switch the four symmetry reductions off.
#include "mlib.c"
#include <time.h>
#include <immintrin.h>
static float ST[13][16] __attribute__((aligned(64)));
static double tm_gen=0,tm_psd=0,tm_prune=0,tm_join=0,tm_build=0;
static inline double nowsec(void){ struct timespec ts; clock_gettime(CLOCK_MONOTONIC,&ts); return ts.tv_sec+1e-9*ts.tv_nsec; }
#define NV 1107
#define NU 2214
#define NPC 14
static i8 Uv[NU][9]; static int RF[NU]; static int OSIZE[9]; static int ZCNT[NU];
static i8 (*G1)[NU][9]; static i8 G0[NU][9];
static const int CLS0[5]={0,1,2,3,6};
static const int PT[14]={243,-6,-6,-6,-6,0,0,0,0,0,0,0,0,0};
static u64 PW[14];
static float S14[14][LANES];
static void pinit(void){
  for(int o=0;o<9;o++) OSIZE[o]=0; for(int z=0;z<P;z++) OSIZE[orb_of[z]]++;
  int cnt=0; int v[8];
  for(int code=0;code<6561;code++){ int c=code,s=0; for(int i=0;i<8;i++){ v[i]=c%3-1; c/=3; s+=v[i]; } if(s) continue;
    for(int i=0;i<8;i++){ Uv[cnt][1+i]=v[i]; Uv[NV+cnt][1+i]=v[i]; } Uv[cnt][0]=0; Uv[NV+cnt][0]=1; cnt++; }
  if(cnt!=NV){ fprintf(stderr,"bad NV\n"); exit(1); }
  for(int i=0;i<NU;i++){ i8 r[9]; for(int o=0;o<9;o++) r[o]=Uv[i][NEGO[o]]; int f=-1; for(int j=0;j<NU;j++) if(!memcmp(Uv[j],r,9)){ f=j; break; } RF[i]=f; }
  G1=malloc((size_t)NU*NU*9);
  for(int a=0;a<NU;a++) for(int b=0;b<NU;b++) for(int op=0;op<9;op++){ int s=0; for(int o=0;o<9;o++) if(Uv[a][o]) for(int o2=0;o2<9;o2++) s+=Uv[a][o]*Uv[b][o2]*NN[op][o][o2]; G1[a][b][op]=(i8)s; }
  for(int a=0;a<NU;a++) memcpy(G0[a],G1[a][a],9);
  for(int i=0;i<NU;i++){ int z=0; for(int o=0;o<9;o++) if(Uv[i][o]==0) z++; ZCNT[i]=z; }
  u64 s=12345; for(int j=0;j<14;j++) PW[j]=sm(&s)|1;
  // PSD matrix for the 13 M columns and 14 character classes of Z123
  int colb[13],colo[13];
  for(int c=0;c<13;c++){ int b=MC_B[c],o=MC_O[c]; int b2=(3-b)%3, o2=NEGO[o]; if(b2<b||(b2==b&&o2<o)){ b=b2; o=o2; } colb[c]=b; colo[c]=o; }
  int seen[3][9]; memset(seen,0,sizeof seen); int nr=0;
  memset(S14,0,sizeof S14);
  for(int beta=0;beta<3;beta++) for(int j=0;j<P;j++){
    int o=orb_of[j], o2=orb_of[(P-j)%P]; int kb=beta, ko=o; int kb2=(3-beta)%3, ko2=o2; if(kb2<kb||(kb2==kb&&ko2<ko)){ kb=kb2; ko=ko2; }
    if(seen[kb][ko]) continue; seen[kb][ko]=1;
    for(int b=0;b<3;b++) for(int z=0;z<P;z++){
      int oo=orb_of[z]; int cb=b, co=oo; int cb2=(3-b)%3, co2=NEGO[oo]; if(cb2<cb||(cb2==cb&&co2<co)){ cb=cb2; co=co2; }
      if(cb==0&&co==0) continue;
      int col=-1; for(int c=0;c<13;c++) if(colb[c]==cb&&colo[c]==co) col=c; if(col<0){ fprintf(stderr,"bad class\n"); exit(1); }
      double ang=2*M_PI*((double)(beta*b)/3.0+(double)((long)j*z%P)/(double)P);
      S14[nr][col]+=(float)cos(ang); }
    nr++; }
  if(nr!=14){ fprintf(stderr,"bad nr %d\n",nr); exit(1); }
}
static void stinit(void){ for(int c=0;c<13;c++) for(int r=0;r<16;r++) ST[c][r]=(r<14)?S14[r][c]:0.0f; }
static inline void pcv(int u0,int u1,int u2,i8*out){
  for(int k=0;k<5;k++){ int o=CLS0[k]; out[k]=G0[u0][o]+G0[u1][o]+G0[u2][o]; }
  for(int o=0;o<9;o++) out[5+o]=G1[u0][u1][o]+G1[u1][u2][o]+G1[u2][u0][o];
  out[14]=out[15]=0;
}
static inline u64 hv14(const i8*v){ u64 h=0; for(int j=0;j<14;j++) h+=PW[j]*(u64)(int64_t)v[j]; return h; }
// ---------- lists
static int KC[3]={6,7,8}; static int kmin[3], kmax[3], kr[3]; static long NK;
static inline u32 keyof(const i8*v){ u32 k=0; for(int t=0;t<3;t++) k=k*kr[t]+(u32)(v[KC[t]]-kmin[t]); return k; }
typedef struct { long n; u32*id; u64*h; long*gs; } List;
static int sym[2][19];
static inline u32 rand_cid0(u64*rs){ u64 a=sm(rs)%NV; u64 b=sm(rs)%NV; u64 c=sm(rs)%NV; return (u32)((a*NV+b)*NV+c); }
static inline int zeros_of(u32 id){ int c2=id%NV, c1=(id/NV)%NV, c0=id/(NV*NV); return ZCNT[NV+c0]+ZCNT[c1]+ZCNT[c2]; }
static inline u32 rand_cid_f(u64*rs,int nmin){ for(;;){ u32 id=rand_cid0(rs); if(nmin<=0||zeros_of(id)>=nmin) return id; } }
static void obj_c(u32 id,i8*v){ int c2=id%NV, c1=(id/NV)%NV, c0=id/(NV*NV); pcv(NV+c0,c1,c2,v); }
static void obj_s(int cen,u32 id,i8*v){ int s=id/NV, y=id%NV; pcv(sym[cen][s],y,RF[y],v); }
// ---------- tuple processing
static double st_E=0, st_Eby[64]; static double st_E2[27*27]; static long st_N2[27*27]; static long hist4[22*22*22*22]; static long st_tuples=0, st_gate=0, st_sols=0, st_pairs=0, st_skipbig=0, st_prune_empty=0; static long hist_n[64];
static FILE*fsol; static int NMIN; static int NMAXVAR=22; static FILE*fskip; static long MAXB=10000000L;
static void slices_of(int typ,u32 id,int*p){ // fill p[27] (p[9b+o])
  int u0,u1,u2;
  if(typ==0||typ==1){ int cen=typ; int s=id/NV,y=id%NV; u0=sym[cen][s]; u1=y; u2=RF[y]; }   // typ 0 = B-type (center 0), 1 = A-type (center 1)
  else { int c2=id%NV, c1=(id/NV)%NV, c0=id/(NV*NV); u0=NV+c0; u1=c1; u2=c2; }
  for(int o=0;o<9;o++){ p[o]=Uv[u0][o]; p[9+o]=Uv[u1][o]; p[18+o]=Uv[u2][o]; }
}
static int plant_mode=0; static long plant_ok=0, plant_fail=0;
typedef struct { u64 key; u32 x,y; } PEnt;
#define NPART 512

static int dedupe=0, norefl=0;
static int SIG[4][NU]; static int symidx[2][NU];
static void sinit(void){
  int g=6; { int x=1,ord=0; do{ x=(x*g)%P; ord++; }while(x!=1); if(ord!=40){ fprintf(stderr,"g=6 not primitive\n"); exit(1); } }
  for(int k=0;k<4;k++){ int gk=1; for(int t=0;t<k;t++) gk=(gk*g)%P; int sg[9]; for(int o=0;o<9;o++) sg[o]=orb_of[(gk*OREP[o])%P];
    for(int u=0;u<NU;u++){ i8 v2[9]; for(int o=0;o<9;o++) v2[sg[o]]=Uv[u][o]; int f=-1; for(int j=0;j<NU;j++) if(!memcmp(Uv[j],v2,9)){ f=j; break; } if(f<0){ fprintf(stderr,"sigma lookup fail\n"); exit(1); } SIG[k][u]=f; } }
  for(int c=0;c<2;c++){ for(int u=0;u<NU;u++) symidx[c][u]=-1; for(int s=0;s<19;s++) symidx[c][sym[c][s]]=s; }
}
static int canonA(u32 id){ int s=id/NV, y=id%NV; int u0=sym[1][s]; for(int k=1;k<4;k++){ int u0p=SIG[k][u0], yp=SIG[k][y]; int sp=symidx[1][u0p]; if(sp<0||yp>=NV){ fprintf(stderr,"canonA map error\n"); exit(1); } u32 idp=(u32)(sp*NV+yp); if(idp<id) return 0; } return 1; }
static unsigned char stabA[19*NV];
static int sigsym(int k,int cen,u32 id){ int s=id/NV, y=id%NV; int u0=sym[cen][s]; int u0p=SIG[k][u0], yp=SIG[k][y]; int sp=symidx[cen][u0p]; return sp*NV+yp; }
static void init_stab(void){ for(u32 i=0;i<19*NV;i++){ unsigned char m=0; for(int k=1;k<4;k++) if((u32)sigsym(k,1,i)==i) m|=(1<<k); stabA[i]=m; } }
static int canonB(u32 aid,u32 bid){ unsigned char m=stabA[aid]; if(!m) return 1; for(int k=1;k<4;k++) if(m&(1<<k)){ if((u32)sigsym(k,0,bid)<bid) return 0; } return 1; }
static inline u32 reflid(u32 id){ int c2=id%NV, c1=(id/NV)%NV, c0=id/(NV*NV); int u0=NV+c0; int r0=RF[u0]-NV, r1=RF[c2], r2=RF[c1]; return (u32)(((long)r0*NV+r1)*NV+r2); }
static void build_band(int lo,int hi,u32**out,long*nout){
  long cap=1<<20, n=0; u32*ids=malloc(cap*4);
  for(int c0=0;c0<NV;c0++){ int z0=ZCNT[NV+c0];
    for(int c1=0;c1<NV;c1++){ int z01=z0+ZCNT[c1]; if(z01+9<lo) continue;
      for(int c2=0;c2<NV;c2++){ int z=z01+ZCNT[c2]; if(z<lo||z>hi) continue; u32 id=(u32)(((long)c0*NV+c1)*NV+c2); if(!norefl && reflid(id)<id) continue; if(n>=cap){ cap*=2; ids=realloc(ids,cap*4); } ids[n++]=id; } } }
  *out=ids; *nout=n;
}
static Block BL[4];
static i16 *LV[4]; static long LN[4]; static float *LP[4];
static i16 *SV[4]; static long SN[4]; static u32 *SI[4];          // survivors: vectors and original Gray indices
static long capv[4]={0,0,0,0};
static int dumpt=0;
static void process_tuple(u32 ida,u32 idb,u32 idc,u32 idd){
  st_tuples++;
  if(dumpt) fprintf(fsol,"PTUPLE %u %u %u %u\n",ida,idb,idc,idd);
  double tb0=nowsec(); int p[27]; int ntot=0; int nn[4]; int wm[4];
  u32 ids[4]={ida,idb,idc,idd};
  // block order: 0 = A (sym, center 1), 1 = B (sym, center 0), 2 = C, 3 = D
  for(int k=0;k<4;k++){
    int typ=(k==0)?1:(k==1)?0:2;
    slices_of(typ,ids[k],p);
    int wp=0; for(int b=0;b<3;b++) for(int o=0;o<9;o++) if(p[9*b+o]) wp+=OSIZE[o];
    wm[k]=123-wp;
    int bt=(k<2)?1:0;
    int n=build_block(bt,p,&BL[k]); if(n<0){ st_skipbig++; return; }
    nn[k]=n; ntot+=n;
  }
  hist_n[ntot<64?ntot:63]++;
  if(ntot<NMIN) return;
  for(int k=0;k<4;k++) if(nn[k]>NMAXVAR){ st_skipbig++; if(fskip){ fprintf(fskip,"SKIP %u %u %u %u n=%d,%d,%d,%d\n",ida,idb,idc,idd,nn[0],nn[1],nn[2],nn[3]); fflush(fskip); } return; }
  tm_build+=nowsec()-tb0; st_gate++;
  { double zz,ld; double pv=proxy4(&BL[0],&BL[1],&BL[2],&BL[3],&zz,&ld); if(pv>-1e29){ double Et=pow(10.0,pv)*8192.0; st_E+=Et; st_Eby[ntot<64?ntot:63]+=Et; st_E2[nn[2]*27+nn[3]]+=Et; st_N2[nn[2]*27+nn[3]]++; } hist4[((nn[0]*22+nn[1])*22+nn[2])*22+nn[3]]++; }
  for(int k=0;k<4;k++){
    long N=1L<<nn[k]; if(N>capv[k]){ capv[k]=N; LV[k]=realloc(LV[k],N*LANES*2); { if(LP[k]) free(LP[k]); LP[k]=aligned_alloc(64,(size_t)N*64); } SV[k]=realloc(SV[k],N*LANES*2); SI[k]=realloc(SI[k],N*4); }
    double tg0=nowsec(); gen_list(&BL[k],LV[k]); LN[k]=N; tm_gen+=nowsec()-tg0; double tp0=nowsec();
    // PSD of each member for the 14 character classes
    { __m512 w0=_mm512_set1_ps((float)wm[k]); for(long m=0;m<N;m++){ const i16*v=LV[k]+m*LANES; __m512 acc=w0; for(int c=0;c<13;c++) acc=_mm512_fmadd_ps(_mm512_load_ps(ST[c]),_mm512_set1_ps((float)v[c]),acc); _mm512_store_ps(LP[k]+m*16,acc); } }
    tm_psd+=nowsec()-tp0;
  }
  double tpr0=nowsec(); float PSI[14]; for(int r=0;r<14;r++) PSI[r]=249.0f;
  i16 T[LANES]; memset(T,0,sizeof T); long pm[4]={0,0,0,0};
  if(plant_mode){ int wsum=0; for(int k=0;k<4;k++) wsum+=wm[k]; for(int k=0;k<4;k++){ pm[k]=rand()%LN[k]; for(int c=0;c<LANES;c++) T[c]+=LV[k][pm[k]*LANES+c]; }
    for(int r=0;r<14;r++){ float s2=(float)wsum; for(int c=0;c<13;c++) s2+=S14[r][c]*T[c]; PSI[r]=s2; } }
  // iterative pruning: PSD_k(chi) <= 249 - sum of the other lists' minima, also >= 0 trivially
  char *alive[4]; for(int k=0;k<4;k++){ alive[k]=malloc(LN[k]); memset(alive[k],1,LN[k]); }
  float mn[4][16]; long cntalive[4];
  for(int k=0;k<4;k++) cntalive[k]=LN[k];
  for(int it=0;it<6;it++){
    for(int k=0;k<4;k++){ __m512 mv=_mm512_set1_ps(1e30f); for(long m=0;m<LN[k];m++) if(alive[k][m]) mv=_mm512_min_ps(mv,_mm512_load_ps(LP[k]+m*16)); _mm512_storeu_ps(mn[k],mv); for(int r=14;r<16;r++) mn[k][r]=0; }
    int changed=0;
    for(int k=0;k<4;k++){ float bound[16]; for(int r=0;r<14;r++){ float s=PSI[r]; for(int j=0;j<4;j++) if(j!=k) s-=mn[j][r]; bound[r]=s+0.01f; } bound[14]=bound[15]=1e30f; __m512 bv=_mm512_loadu_ps(bound);
      long c=0; for(long m=0;m<LN[k];m++){ if(!alive[k][m]) continue; __mmask16 bad=_mm512_cmp_ps_mask(_mm512_load_ps(LP[k]+m*16),bv,_CMP_GT_OQ); if(bad){ alive[k][m]=0; changed=1; } else c++; }
      cntalive[k]=c; }
    if(!changed) break;
  }
  tm_prune+=nowsec()-tpr0; int empty=0; for(int k=0;k<4;k++) if(cntalive[k]==0) empty=1;
  if(empty){ st_prune_empty++; for(int k=0;k<4;k++) free(alive[k]); return; }
  for(int k=0;k<4;k++){ long c=0; for(long m=0;m<LN[k];m++) if(alive[k][m]){ memcpy(SV[k]+c*LANES,LV[k]+m*LANES,LANES*2); SI[k][c]=(u32)m; c++; } SN[k]=c; free(alive[k]); }
  // choose the pairing with the fewest pair operations: (x,z)|(y,w)
  static const int PR[3][4]={{0,1,2,3},{0,2,1,3},{0,1,3,2}};   // join4(V[a],V[b],V[c],V[d]) pairs (a,c) with (b,d)
  int best=0; double bc=1e300;
  for(int t=0;t<3;t++){ double c1=(double)SN[PR[t][0]]*SN[PR[t][2]], c2=(double)SN[PR[t][1]]*SN[PR[t][3]]; double c=c1>c2?c1:c2; if(c<bc){ bc=c; best=t; } }
  st_pairs+=(long)(bc*2);
  int*out=malloc(4*sizeof(int)*1000); long nm=0;
  const int*q=PR[best];
  double tj0=nowsec(); long ns=join4(SV[q[0]],SN[q[0]],SV[q[1]],SN[q[1]],SV[q[2]],SN[q[2]],SV[q[3]],SN[q[3]],T,out,1000,&nm);
  tm_join+=nowsec()-tj0;
  if(plant_mode){ int found=0; for(long s=0;s<ns&&s<1000;s++){ u32 gi[4]; gi[q[0]]=SI[q[0]][out[4*s]]; gi[q[1]]=SI[q[1]][out[4*s+1]]; gi[q[2]]=SI[q[2]][out[4*s+2]]; gi[q[3]]=SI[q[3]][out[4*s+3]]; if(gi[0]==pm[0]&&gi[1]==pm[1]&&gi[2]==pm[2]&&gi[3]==pm[3]) found=1; }
    if(found) plant_ok++; else plant_fail++; free(out); return; }
  if(ns>0){
    for(long s=0;s<ns&&s<1000;s++){ u32 gi[4]; gi[q[0]]=SI[q[0]][out[4*s]]; gi[q[1]]=SI[q[1]][out[4*s+1]]; gi[q[2]]=SI[q[2]][out[4*s+2]]; gi[q[3]]=SI[q[3]][out[4*s+3]];
      fprintf(fsol,"SOL %u %u %u %u %u %u %u %u\n",ida,idb,idc,idd,gi[0],gi[1],gi[2],gi[3]); }
    fflush(fsol); st_sols+=ns;
  }
  free(out);
}
int main(int argc,char**argv){
  if(argc<3){ fprintf(stderr,"usage: band NMIN OUTPREFIX [S0LO S0HI]   (see the header of band.c)\n"); return 1; }
  NMIN=atoi(argv[1]); const char*pref=argv[2];
  int lo0=argc>4?atoi(argv[3]):-100000, hi0=argc>4?atoi(argv[4]):100000;
  if(getenv("MAXB")) MAXB=atol(getenv("MAXB"));
  minit(); pinit(); stinit(); plant_mode = getenv("PLANTM")!=0; dumpt = getenv("DUMPT")!=0; srand(12345);
  char fn[512]; snprintf(fn,sizeof fn,"%s.sol",pref); fsol=fopen(fn,"a");
  char fprog[512]; snprintf(fprog,sizeof fprog,"%s.prog",pref); FILE*fp=fopen(fprog,"a"); { char fs[512]; snprintf(fs,sizeof fs,"%s.skip",pref); fskip=fopen(fs,"a"); }
  int ns2[2]={0,0};
  for(int i=0;i<NU;i++){ int ok=1; for(int o=0;o<9;o++) if(Uv[i][o]!=Uv[i][NEGO[o]]) ok=0; if(ok) sym[(int)Uv[i][0]][ns2[(int)Uv[i][0]]++]=i; }
  sinit(); init_stab(); int nocanb=getenv("NOCANB")!=0;
  int cl=18,ch=99,dl=18,dh=99; { const char*e=getenv("CBAND"); if(e) sscanf(e,"%d:%d",&cl,&ch); e=getenv("DBAND"); if(e) sscanf(e,"%d:%d",&dl,&dh); else { dl=cl; dh=ch; } }
  int nocana=getenv("NOCANA")!=0; norefl=getenv("NOREFL")!=0; int nodedup=getenv("NODEDUP")!=0;
  u32*cids,*dids; long ncid,ndid; build_band(cl,ch,&cids,&ncid); int same=(cl==dl&&ch==dh); if(same){ dids=cids; ndid=ncid; } else build_band(dl,dh,&dids,&ndid); dedupe=same&&!nodedup;
  fprintf(stderr,"bands C %d:%d (%ld canonical objects)  D %d:%d (%ld)  same=%d\n",cl,ch,ncid,dl,dh,ndid,same);
  // A (center 1, canonical under the multiplier group) and B (center 0): lists
  List A,B; long nAall=19L*NV; A.id=malloc(nAall*4); B.id=malloc(nAall*4); A.h=malloc(nAall*8); B.h=malloc(nAall*8);
  static i8 (*AV)[16]; static i8 (*BV)[16]; AV=malloc(nAall*16); BV=malloc(nAall*16);
  long na=0; for(long i=0;i<nAall;i++){ if(!nocana && !canonA((u32)i)) continue; A.id[na]=(u32)i; obj_s(1,(u32)i,AV[na]); A.h[na]=hv14(AV[na]); na++; }
  A.n=na; B.n=nAall; for(long i=0;i<nAall;i++){ B.id[i]=(u32)i; obj_s(0,(u32)i,BV[i]); B.h[i]=hv14(BV[i]); }
  fprintf(stderr,"A canonical objects %ld of %ld, B %ld\n",A.n,nAall,B.n);
  for(int t=0;t<3;t++){ kmin[t]=127; kmax[t]=-128; }
  for(long i=0;i<A.n;i++){ for(int t=0;t<3;t++){ int x=AV[i][KC[t]]; if(x<kmin[t]) kmin[t]=x; if(x>kmax[t]) kmax[t]=x; } }
  for(long i=0;i<B.n;i++){ for(int t=0;t<3;t++){ int x=BV[i][KC[t]]; if(x<kmin[t]) kmin[t]=x; if(x>kmax[t]) kmax[t]=x; } }
  for(int side=0;side<2;side++){ u32*ids=side?dids:cids; long n=side?ndid:ncid; for(long i=0;i<n;i++){ i8 v[16]; obj_c(ids[i],v); for(int t=0;t<3;t++){ int x=v[KC[t]]; if(x<kmin[t]) kmin[t]=x; if(x>kmax[t]) kmax[t]=x; } } }
  for(int t=0;t<3;t++) kr[t]=kmax[t]-kmin[t]+1;
  NK=(long)kr[0]*kr[1]*kr[2]; fprintf(stderr,"key ranges %d %d %d NK=%ld\n",kr[0],kr[1],kr[2],NK);
  List C={ncid,0,0,0}, D={ndid,0,0,0};
  List*AB[2]={&A,&B}; i8 (*ABV[2])[16]={AV,BV};
  for(int t=0;t<2;t++){ List*L=AB[t]; long n=L->n; u32*key=malloc(n*4); L->gs=calloc(NK+1,sizeof(long)); for(long i=0;i<n;i++){ key[i]=keyof(ABV[t][i]); L->gs[key[i]+1]++; } for(long k=0;k<NK;k++) L->gs[k+1]+=L->gs[k];
    long*pos=malloc(NK*sizeof(long)); memcpy(pos,L->gs,NK*sizeof(long)); u32*nid=malloc(n*4); u64*nh=malloc(n*8); for(long i=0;i<n;i++){ long q=pos[key[i]]++; nid[q]=L->id[i]; nh[q]=L->h[i]; } free(pos); free(key); free(L->id); free(L->h); L->id=nid; L->h=nh; }
  for(int side=0;side<2;side++){ if(side==1&&same){ D=C; break; } List*L=side?&D:&C; u32*ids=side?dids:cids; long n=L->n; if(n<=0){ fprintf(stderr,"empty band\n"); return 4; } u32*key=malloc(n*4); L->gs=calloc(NK+1,sizeof(long));
    for(long i=0;i<n;i++){ i8 v[16]; obj_c(ids[i],v); key[i]=keyof(v); L->gs[key[i]+1]++; }
    for(long k=0;k<NK;k++) L->gs[k+1]+=L->gs[k];
    long*pos=malloc(NK*sizeof(long)); memcpy(pos,L->gs,NK*sizeof(long)); L->id=malloc(n*4); L->h=malloc(n*8);
    for(long i=0;i<n;i++){ i8 v[16]; obj_c(ids[i],v); long q=pos[key[i]]++; L->id[q]=ids[i]; L->h[q]=hv14(v); }
    free(pos); free(key); }
  free(cids); if(!same) free(dids);
  fprintf(stderr,"lists built: nC=%ld nD=%ld\n",C.n,D.n);
  // nonempty A,B groups
  long*gA=malloc(NK*sizeof(long)),*gB=malloc(NK*sizeof(long)); long nga=0,ngb=0;
  for(long k=0;k<NK;k++){ if(A.gs[k+1]>A.gs[k]) gA[nga++]=k; if(B.gs[k+1]>B.gs[k]) gB[ngb++]=k; }
  int *ga0=malloc(nga*4),*ga1=malloc(nga*4),*ga2=malloc(nga*4),*gb0=malloc(ngb*4),*gb1=malloc(ngb*4),*gb2=malloc(ngb*4);
  for(long i=0;i<nga;i++){ ga2[i]=(int)(gA[i]%kr[2])+kmin[2]; ga1[i]=(int)((gA[i]/kr[2])%kr[1])+kmin[1]; ga0[i]=(int)(gA[i]/((long)kr[1]*kr[2]))+kmin[0]; }
  for(long i=0;i<ngb;i++){ gb2[i]=(int)(gB[i]%kr[2])+kmin[2]; gb1[i]=(int)((gB[i]/kr[2])%kr[1])+kmin[1]; gb0[i]=(int)(gB[i]/((long)kr[1]*kr[2]))+kmin[0]; }
  int Tk[3]={PT[KC[0]],PT[KC[1]],PT[KC[2]]};
  u64 hT=0; for(int j=0;j<14;j++) hT+=PW[j]*(u64)(int64_t)PT[j];
  size_t tcap=1<<12; u64*tk=malloc(tcap*8); u64*tv=malloc(tcap*8);
  PEnt*pbuf_ac=0,*pbuf_bd=0; long pbuf_ac_cap=0,pbuf_bd_cap=0; static long pcnt_ac[NPART], pcnt_bd[NPART];
  int s0min=kmin[0]*2, s0max=kmax[0]*2, s1min=kmin[1]*2, s1max=kmax[1]*2, s2min=kmin[2]*2, s2max=kmax[2]*2;
  long npairs_total=0; time_t t0=time(0);
  static char doneS0[4096]; { FILE*pf=fopen(fprog,"r"); if(pf){ char line[4096]; while(fgets(line,sizeof line,pf)){ int v; if(sscanf(line,"s0=%d ",&v)==1 && v-s0min>=0 && v-s0min<4096) doneS0[v-s0min]=1; } fclose(pf); } }
  int nS0=s0max-s0min+1; int*s0ord=malloc(nS0*sizeof(int)); { int cnt0=0; int cen=-2; for(int d=0; cnt0<nS0; d++){ int a=cen-d, b2=cen+d; if(a>=s0min&&a<=s0max) s0ord[cnt0++]=a; if(d>0&&b2>=s0min&&b2<=s0max) s0ord[cnt0++]=b2; if(d>400) break; } nS0=cnt0; }
  for(int is0=0;is0<nS0;is0++){ int s0=s0ord[is0]; if(s0<lo0||s0>hi0||doneS0[s0-s0min]) continue;
   for(int s1=s1min;s1<=s1max;s1++) for(int s2=s2min;s2<=s2max;s2++){
    long cnt=0;
    for(long ig=0;ig<nga;ig++){ int c0=s0-ga0[ig]-kmin[0], c1=s1-ga1[ig]-kmin[1], c2=s2-ga2[ig]-kmin[2]; if(c0<0||c0>=kr[0]||c1<0||c1>=kr[1]||c2<0||c2>=kr[2]) continue;
      long kc=((long)c0*kr[1]+c1)*kr[2]+c2; cnt+=(A.gs[gA[ig]+1]-A.gs[gA[ig]])*(C.gs[kc+1]-C.gs[kc]); }
    if(!cnt) continue;
    int t0k=Tk[0]-s0, t1k=Tk[1]-s1, t2k=Tk[2]-s2; long cnt2=0;
    for(long ig=0;ig<ngb;ig++){ int d0=t0k-gb0[ig]-kmin[0], d1=t1k-gb1[ig]-kmin[1], d2=t2k-gb2[ig]-kmin[2]; if(d0<0||d0>=kr[0]||d1<0||d1>=kr[1]||d2<0||d2>=kr[2]) continue;
      long kd=((long)d0*kr[1]+d1)*kr[2]+d2; cnt2+=(B.gs[gB[ig]+1]-B.gs[gB[ig]])*(D.gs[kd+1]-D.gs[kd]); }
    if(!cnt2) continue;
    npairs_total+=cnt+cnt2;
    long npass=(cnt+MAXB-1)/MAXB; if(npass<1) npass=1;
    for(long pass=0;pass<npass;pass++){
      long exp_ac=cnt/npass+1, exp_bd=cnt2/npass+1;
      int retry=0; double slack=1.35;
      for(;;){
        long capR_ac=(long)((double)exp_ac/NPART*slack)+256, capR_bd=(long)((double)exp_bd/NPART*slack)+256;
        if(pbuf_ac_cap<capR_ac*NPART){ free(pbuf_ac); pbuf_ac_cap=capR_ac*NPART; pbuf_ac=malloc(pbuf_ac_cap*sizeof(PEnt)); }
        if(pbuf_bd_cap<capR_bd*NPART){ free(pbuf_bd); pbuf_bd_cap=capR_bd*NPART; pbuf_bd=malloc(pbuf_bd_cap*sizeof(PEnt)); }
        memset(pcnt_ac,0,sizeof pcnt_ac); memset(pcnt_bd,0,sizeof pcnt_bd);
        int overflow=0;
        for(long ig=0;ig<nga&&!overflow;ig++){ int c0=s0-ga0[ig]-kmin[0], c1=s1-ga1[ig]-kmin[1], c2=s2-ga2[ig]-kmin[2]; if(c0<0||c0>=kr[0]||c1<0||c1>=kr[1]||c2<0||c2>=kr[2]) continue;
          long kc=((long)c0*kr[1]+c1)*kr[2]+c2; long ka=gA[ig];
          for(long x=A.gs[ka];x<A.gs[ka+1]&&!overflow;x++){ u64 ha=A.h[x];
            for(long y=C.gs[kc];y<C.gs[kc+1];y++){ u64 key=ha+C.h[y]; if(npass>1 && (long)((key>>20)%(u64)npass)!=pass) continue; int pr=(int)((key>>44)&(NPART-1));
              long c=pcnt_ac[pr]; if(c>=capR_ac){ overflow=1; break; } PEnt*e=&pbuf_ac[(long)pr*capR_ac+c]; e->key=key; e->x=(u32)x; e->y=(u32)y; pcnt_ac[pr]=c+1; } } }
        for(long ig=0;ig<ngb&&!overflow;ig++){ int d0=t0k-gb0[ig]-kmin[0], d1=t1k-gb1[ig]-kmin[1], d2=t2k-gb2[ig]-kmin[2]; if(d0<0||d0>=kr[0]||d1<0||d1>=kr[1]||d2<0||d2>=kr[2]) continue;
          long kd=((long)d0*kr[1]+d1)*kr[2]+d2; long kb=gB[ig];
          for(long x=B.gs[kb];x<B.gs[kb+1]&&!overflow;x++){ u64 hb=hT-B.h[x];
            for(long y=D.gs[kd];y<D.gs[kd+1];y++){ u64 key=hb-D.h[y]; if(npass>1 && (long)((key>>20)%(u64)npass)!=pass) continue; int pr=(int)((key>>44)&(NPART-1));
              long c=pcnt_bd[pr]; if(c>=capR_bd){ overflow=1; break; } PEnt*e=&pbuf_bd[(long)pr*capR_bd+c]; e->key=key; e->x=(u32)x; e->y=(u32)y; pcnt_bd[pr]=c+1; } } }
        if(overflow){ slack*=1.6; if(++retry>6){ fprintf(stderr,"partition overflow\n"); return 3; } continue; }
        for(int pr=0;pr<NPART;pr++){ long nac=pcnt_ac[pr], nbd=pcnt_bd[pr]; if(!nac||!nbd) continue;
          size_t nb2=4; while(((size_t)1<<nb2)<(size_t)nac*2+8) nb2++; size_t cp2=(size_t)1<<nb2;
          if(cp2>tcap){ free(tk); free(tv); tcap=cp2; tk=malloc(tcap*8); tv=malloc(tcap*8); }
          for(size_t i=0;i<cp2;i++) tv[i]=~0ULL;
          PEnt*pa=&pbuf_ac[(long)pr*capR_ac]; for(long i=0;i<nac;i++){ u64 key=pa[i].key; size_t p=(size_t)((key*0x9E3779B97F4A7C15ULL)>>(64-nb2)); while(tv[p]!=~0ULL) p=(p+1)&(cp2-1); tk[p]=key; tv[p]=((u64)pa[i].x<<32)|(u64)pa[i].y; }
          PEnt*pb=&pbuf_bd[(long)pr*capR_bd];
          for(long i=0;i<nbd;i++){ u64 key=pb[i].key; size_t p=(size_t)((key*0x9E3779B97F4A7C15ULL)>>(64-nb2));
            while(tv[p]!=~0ULL){ if(tk[p]==key){ long xa=(long)(tv[p]>>32), yc=(long)(tv[p]&0xffffffffu); long xb=pb[i].x, yd=pb[i].y;
                i8 va[16],vb[16],vc[16],vd[16]; obj_s(1,A.id[xa],va); obj_s(0,B.id[xb],vb); obj_c(C.id[yc],vc); obj_c(D.id[yd],vd); int ok=1; for(int j=0;j<14;j++) if(va[j]+vb[j]+vc[j]+vd[j]!=PT[j]) ok=0;
                if(ok && !(dedupe && C.id[yc]>D.id[yd]) && (nocanb || canonB(A.id[xa],B.id[xb]))) process_tuple(A.id[xa],B.id[xb],C.id[yc],D.id[yd]); }
              p=(p+1)&(cp2-1); } } }
        break;
      }
    } // pass
   }
   fprintf(fp,"s0=%d tuples=%ld gate=%ld sols=%ld pairs=%ld mpairs=%ld skipbig=%ld pruneempty=%ld t=%ld Esum=%.6g\n",s0,st_tuples,st_gate,st_sols,npairs_total,st_pairs,st_skipbig,st_prune_empty,(long)(time(0)-t0),st_E); fflush(fp);
   if(s0%4==0){ fprintf(stderr,"s0=%d tuples=%ld gate=%ld sols=%ld t=%lds\n",s0,st_tuples,st_gate,st_sols,(long)(time(0)-t0)); }
  }
  fprintf(fp,"TIMERS build=%.1f gen=%.1f psd=%.1f prune=%.1f join=%.1f\n",tm_build,tm_gen,tm_psd,tm_prune,tm_join);
  fprintf(fp,"PLANT ok=%ld fail=%ld\n",plant_ok,plant_fail); fprintf(fp,"E2(nC,nD:Esum,count):"); for(int a=0;a<27;a++) for(int b=0;b<27;b++) if(st_N2[a*27+b]) fprintf(fp," %d,%d:%.5g,%ld",a,b,st_E2[a*27+b],st_N2[a*27+b]); fprintf(fp,"\n"); fprintf(fp,"EBYN:"); for(int i=0;i<64;i++) if(st_Eby[i]>0) fprintf(fp," %d:%.5g",i,st_Eby[i]); fprintf(fp,"\n"); fprintf(fp,"DONE tuples=%ld gate=%ld sols=%ld hist:",st_tuples,st_gate,st_sols); for(int i=0;i<64;i++) if(hist_n[i]) fprintf(fp," %d:%ld",i,hist_n[i]); fprintf(fp,"\n"); fclose(fp); fclose(fsol);
  return 0;
}
