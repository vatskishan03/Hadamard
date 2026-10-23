// qsa.c: simulated annealing for the quaternary skew GS orbit-level system (4 blocks).  usage: qsa shapefile seed seconds T0 T1 itersPerRun
// shapefile: as written by qdrv.shape (e.g. shapes/t_191.txt after `python3 -c "from qdrv import *; shape(191,[180],'shapes/t_191.txt')"`);
// each run anneals from a random start, geometric cooling T0 -> T1 over itersPerRun moves; prints a SOL line (four orbit-exponent
// strings) if the sum of the squared PAF values reaches 0.  It does not converge even at n = 11 (v = 191): NOTES.md section 88.
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
typedef unsigned long long u64;
static int v,n,d; static int sz[64],neg_[64]; static int *Nt;
static u64 rs;
static inline u64 rnd(){ rs^=rs<<13; rs^=rs>>7; rs^=rs<<17; return rs; }
static inline double urand(){ return (rnd()>>11)*(1.0/9007199254740992.0); }
static int e[4][64]; static int pr[64],np;
static int pv[4][64][2]; // PAF per block per class (re,im)
static int tot[64][2];
static void paf_block(const int*ee,int out[][2]){
  for(int c=1;c<n;c++){ int M[4]={0,0,0,0}; const int*N=Nt+(size_t)c*n*n; for(int a=0;a<n;a++){const int*row=N+a*n; int ea=ee[a]; for(int b=0;b<n;b++){int x=row[b]; if(x) M[(ea-ee[b])&3]+=x;}} out[c][0]=M[0]-M[2]; out[c][1]=M[1]-M[3]; }
}
static long fval(){ long f=0; for(int c=1;c<n;c++) f+=(long)tot[c][0]*tot[c][0]+(long)tot[c][1]*tot[c][1]; return f; }
int main(int argc,char**argv){
  FILE*f=fopen(argv[1],"r"); rs=88172645463325252ULL^(u64)atoll(argv[2])*0x9E3779B97F4A7C15ULL; if(!rs) rs=1; for(int i=0;i<20;i++) rnd();
  double secs=atof(argv[3]),T0=atof(argv[4]),T1=atof(argv[5]); long iters=atol(argv[6]);
  if(fscanf(f,"%d %d",&v,&n)!=2) return 1; d=n-1;
  for(int i=0;i<n;i++) if(fscanf(f,"%d",&sz[i])!=1) return 1;
  for(int i=0;i<n;i++) if(fscanf(f,"%d",&neg_[i])!=1) return 1;
  Nt=malloc(sizeof(int)*(size_t)n*n*n);
  for(int c=1;c<n;c++) for(int a=0;a<n;a++) for(int b=0;b<n;b++) if(fscanf(f,"%d",&Nt[(size_t)c*n*n+a*n+b])!=1) return 1;
  fclose(f);
  for(int a=1;a<n;a++) if(a<neg_[a]) pr[np++]=a;
  time_t t0=time(NULL); long runs=0; long bestf=1L<<60;
  while(difftime(time(NULL),t0)<secs){
    runs++;
    for(int X=0;X<4;X++){ e[X][0]=0; for(int a=1;a<n;a++) e[X][a]=rnd()&3; }
    e[0][0]=(rnd()&1)?3:0; for(int i=0;i<np;i++) e[0][neg_[pr[i]]]=(1-e[0][pr[i]])&3;
    for(int X=0;X<4;X++) paf_block(e[X],pv[X]);
    for(int c=1;c<n;c++){ tot[c][0]=tot[c][1]=0; for(int X=0;X<4;X++){ tot[c][0]+=pv[X][c][0]; tot[c][1]+=pv[X][c][1]; } }
    long cur=fval();
    for(long it=0;it<iters;it++){
      double T=T0*pow(T1/T0,(double)it/iters);
      int X=rnd()&3; int newpv[64][2]; int old[64]; int a,na; 
      if(X==0){
        int r=rnd()%(np+1); memcpy(old,e[0],sizeof(int)*n);
        if(r==np){ e[0][0]=e[0][0]==0?3:0; } else { int p=pr[r]; int nv=(e[0][p]+1+rnd()%3)&3; e[0][p]=nv; e[0][neg_[p]]=(1-nv)&3; }
        paf_block(e[0],newpv);
      } else {
        a=1+rnd()%(n-1); na=(e[X][a]+1+rnd()%3)&3; memcpy(old,e[X],sizeof(int)*n); e[X][a]=na; 
        paf_block(e[X],newpv); // simple full recompute
      }
      long nf=0; int nt[64][2];
      for(int c=1;c<n;c++){ nt[c][0]=tot[c][0]-pv[X][c][0]+newpv[c][0]; nt[c][1]=tot[c][1]-pv[X][c][1]+newpv[c][1]; nf+=(long)nt[c][0]*nt[c][0]+(long)nt[c][1]*nt[c][1]; }
      long df=nf-cur;
      if(df<=0 || urand()<exp(-df/T)){ cur=nf; memcpy(pv[X],newpv,sizeof(int)*128); for(int c=1;c<n;c++){tot[c][0]=nt[c][0];tot[c][1]=nt[c][1];} if(cur==0) break; }
      else { memcpy(e[X],old,sizeof(int)*n); }
    }
    if(cur<bestf){ bestf=cur; fprintf(stderr,"run %ld best f %ld\n",runs,bestf); }
    if(cur==0){ printf("SOL"); for(int X=0;X<4;X++){ printf(" "); for(int a=0;a<n;a++) printf("%d",e[X][a]); } printf("\n"); fflush(stdout); break; }
  }
  fprintf(stderr,"runs %ld bestf %ld\n",runs,bestf);
  return 0;
}
