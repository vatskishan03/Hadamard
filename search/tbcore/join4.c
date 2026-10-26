// exact 4-list join: V_A[a]+V_B[b]+V_C[c]+V_D[d] == T, by hashing pair sums (A,C) and (B,D), in key chunks.
// input file: int32 nA nB nC nD m; then int16 VA[nA*m], VB, VC, VD; then int32 T[m]
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <time.h>
typedef uint64_t u64; typedef int16_t i16; typedef int32_t i32; typedef uint32_t u32;
static u64 splitmix(u64 *x){ u64 z=(*x+=0x9E3779B97F4A7C15ULL); z=(z^(z>>30))*0xBF58476D1CE4E5B9ULL; z=(z^(z>>27))*0x94D049BB133111EBULL; return z^(z>>31);}
int main(int argc,char**argv){
  if(argc<3){fprintf(stderr,"usage: join4 file nch_log2 [seed]\n");return 1;}
  FILE*f=fopen(argv[1],"rb"); int lch=atoi(argv[2]); u64 seed=argc>3?atoll(argv[3]):12345;
  i32 hdr[5]; if(fread(hdr,4,5,f)!=5) return 2;
  i32 n[4]={hdr[0],hdr[1],hdr[2],hdr[3]}, m=hdr[4];
  i16 *V[4]; for(int k=0;k<4;k++){V[k]=malloc((size_t)n[k]*m*2); if(fread(V[k],2,(size_t)n[k]*m,f)!=(size_t)n[k]*m) return 3;}
  i32 *T=malloc(m*4); if(fread(T,4,m,f)!=(size_t)m) return 4; fclose(f);
  u64 rnd[64]; u64 s=seed; for(int j=0;j<m;j++) rnd[j]=splitmix(&s)|1;
  u64 *h[4]; for(int k=0;k<4;k++){ h[k]=malloc((size_t)n[k]*8); for(i32 i=0;i<n[k];i++){ u64 x=0; for(int j=0;j<m;j++) x+=rnd[j]*(u64)(int64_t)V[k][(size_t)i*m+j]; h[k][i]=x; } }
  u64 hT=0; for(int j=0;j<m;j++) hT+=rnd[j]*(u64)(int64_t)T[j];
  // order: left pairs (A,C), right pairs (B,D). choose left = smaller of nA*nC, nB*nD? keep as is
  u64 nL=(u64)n[0]*n[2], nR=(u64)n[1]*n[3];
  int nch=1<<lch; long nsol=0, nmatch=0;
  u64 per=(nL/nch)+1; int tb=1; while(((u64)1<<tb) < 2*per+16) tb++;
  size_t cap=(size_t)1<<tb;
  u64 *tk=malloc(cap*8); u32 *tv=malloc(cap*4);
  if(!tk||!tv){fprintf(stderr,"alloc fail cap=%zu\n",cap);return 5;}
  clock_t t0=clock();
  for(int ch=0;ch<nch;ch++){
    memset(tv,0xff,cap*4);
    u64 ins=0;
    for(i32 a=0;a<n[0];a++){ u64 ha=h[0][a];
      for(i32 c=0;c<n[2];c++){ u64 key=ha+h[2][c];
        if(lch && (int)(key>>(64-lch))!=ch) continue;
        size_t p=(size_t)((key*0x9E3779B97F4A7C15ULL)>>(64-tb));
        while(tv[p]!=0xffffffffu) p=(p+1)&(cap-1);
        tk[p]=key; tv[p]=(u32)((u64)a*n[2]+c); ins++; } }
    for(i32 b=0;b<n[1];b++){ u64 hb=hT-h[1][b];
      for(i32 d=0;d<n[3];d++){ u64 key=hb-h[3][d];
        if(lch && (int)(key>>(64-lch))!=ch) continue;
        size_t p=(size_t)((key*0x9E3779B97F4A7C15ULL)>>(64-tb));
        while(tv[p]!=0xffffffffu){ if(tk[p]==key){ nmatch++; u64 id=tv[p]; i32 a=id/n[2], c=id%n[2];
            int ok=1; for(int j=0;j<m;j++){ int x=V[0][(size_t)a*m+j]+V[1][(size_t)b*m+j]+V[2][(size_t)c*m+j]+V[3][(size_t)d*m+j]; if(x!=T[j]){ok=0;break;} }
            if(ok){ nsol++; printf("SOL %d %d %d %d\n",a,b,c,d); fflush(stdout);} }
          p=(p+1)&(cap-1);} } }
    if(nch>1 && ch%1==0) fprintf(stderr,"chunk %d/%d ins %llu matches %ld sols %ld t=%.1fs\n",ch+1,nch,(unsigned long long)ins,nmatch,nsol,(double)(clock()-t0)/CLOCKS_PER_SEC);
  }
  fprintf(stderr,"pairs L=%llu R=%llu keymatches=%ld solutions=%ld time=%.1fs\n",(unsigned long long)nL,(unsigned long long)nR,nmatch,nsol,(double)(clock()-t0)/CLOCKS_PER_SEC);
  printf("DONE %ld %ld\n",nmatch,nsol);
  return 0;
}
