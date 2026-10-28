// pjoin.c : bucketed exact four-list join with prefix-agreement histogram (stage 1 of the QH(358) tower, small_analogue.py, test_pjoin.py).
// Lists V0..V3 (int16, n_i x m, columns already in the chosen order), target T[m].  A quadruple (a,b,c,d) is KEPT when the sums agree with T on the first K0 columns;
// for every kept quadruple the length l of the agreeing prefix (K0 <= l <= m) is histogrammed; quadruples with l >= LSTORE are printed.
// Key: h(v) = sum_{j<K0} w_j v_j (mod 2^64), w_j random odd.  (A,C) pairs with h_a+h_c in bucket beta are put into an open-addressing table, (B,D) pairs with T-h_b-h_d in the same bucket probe it.
// every key match is re-verified on the K0 columns.  Buckets = top NBITS bits of the key; part/nparts splits the bucket range between processes.
// build: gcc -O2 -march=native -o pjoin pjoin.c
// usage: pjoin LISTS.bin K0 NBITS SEED PART NPARTS LSTORE OUT.txt   (LISTS.bin as written by pio.write_lists)
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <time.h>
typedef uint64_t u64; typedef int16_t i16; typedef uint32_t u32; typedef int32_t i32;
static u64 sm(u64*x){ u64 z=(*x+=0x9E3779B97F4A7C15ULL); z=(z^(z>>30))*0xBF58476D1CE4E5B9ULL; z=(z^(z>>27))*0x94D049BB133111EBULL; return z^(z>>31);}
typedef struct { u64 h; u32 i; } HI;
static int cmp(const void*a,const void*b){ u64 x=((const HI*)a)->h, y=((const HI*)b)->h; return x<y?-1:(x>y?1:0); }
static int m, K0; static i16*V[4]; static i32 n[4]; static i32*T;
static HI*S[4];
// first index with h >= s
static i32 lb(const HI*a,i32 na,u64 s){ i32 lo=0,hi=na; while(lo<hi){ i32 mid=(lo+hi)>>1; if(a[mid].h<s) lo=mid+1; else hi=mid; } return lo; }
// iterate over indices with h in the circular interval [s, s+W) ; calls via macro-like loop: returns up to two ranges
static void ranges(const HI*a,i32 na,u64 s,u64 W,i32*r0,i32*r1,i32*r2,i32*r3){
  *r2=0; *r3=0;
  if(W==0){ *r0=0; *r1=na; return; }
  u64 e=s+W;
  if(e>s){ *r0=lb(a,na,s); *r1=lb(a,na,e); }
  else { *r0=lb(a,na,s); *r1=na; *r2=0; *r3=lb(a,na,e); }
}
int main(int argc,char**argv){
  if(argc<9){ fprintf(stderr,"usage\n"); return 1; }
  FILE*f=fopen(argv[1],"rb"); K0=atoi(argv[2]); int nbits=atoi(argv[3]); u64 seed=strtoull(argv[4],0,10); int part=atoi(argv[5]), nparts=atoi(argv[6]); int lstore=atoi(argv[7]); FILE*fo=fopen(argv[8],"w");
  i32 hdr[5]; if(fread(hdr,4,5,f)!=5) return 2; m=hdr[0]; for(int k=0;k<4;k++) n[k]=hdr[1+k];
  for(int k=0;k<4;k++){ V[k]=malloc((size_t)n[k]*m*2); if(fread(V[k],2,(size_t)n[k]*m,f)!=(size_t)n[k]*m) return 3; }
  T=malloc(4*m); if(fread(T,4,m,f)!=(size_t)m) return 4; fclose(f);
  u64 w[64]; u64 sd=seed; for(int j=0;j<m;j++) w[j]=sm(&sd)|1;
  for(int k=0;k<4;k++){ S[k]=malloc(sizeof(HI)*n[k]); for(i32 i=0;i<n[k];i++){ u64 x=0; for(int j=0;j<K0;j++) x+=w[j]*(u64)(int64_t)V[k][(size_t)i*m+j]; S[k][i].h=x; S[k][i].i=i; } qsort(S[k],n[k],sizeof(HI),cmp); }
  u64 hT=0; for(int j=0;j<K0;j++) hT+=w[j]*(u64)(int64_t)T[j];
  u64 NB=(u64)1<<nbits; u64 W=(nbits==0)?0:((u64)1<<(64-nbits));   // W==0 means a single bucket = whole circle
  double pairsAC=(double)n[0]*n[2]; u64 per=(u64)(pairsAC/(double)NB*1.15)+1024; int tb=1; while(((u64)1<<tb)<2*per) tb++;
  size_t cap=(size_t)1<<tb; u64*tk=malloc(cap*8), *tv=malloc(cap*8); if(!tk||!tv){ fprintf(stderr,"alloc\n"); return 5; }
  long long hist[80]; memset(hist,0,sizeof hist); long long nmatchkey=0,nfalse=0,npairs_ac=0,npairs_bd=0,overflow=0;
  u64 b0=(NB*(u64)part)/nparts, b1=(NB*(u64)(part+1))/nparts;
  time_t t0=time(0);
  for(u64 beta=b0;beta<b1;beta++){
    u64 low=beta*W;   // bucket [low, low+W)
    for(size_t i=0;i<cap;i++) tv[i]=~0ULL;
    size_t filled=0;
    for(i32 a=0;a<n[0];a++){
      u64 s=low-S[0][a].h; i32 r0,r1,r2,r3; ranges(S[2],n[2],s,W,&r0,&r1,&r2,&r3);
      for(int pass=0;pass<2;pass++){
        i32 lo=pass?r2:r0, hi=pass?r3:r1;
        for(i32 c=lo;c<hi;c++){
          u64 key=S[0][a].h+S[2][c].h; size_t p=(size_t)((key*0x9E3779B97F4A7C15ULL)>>(64-tb));
          while(tv[p]!=~0ULL) p=(p+1)&(cap-1);
          tk[p]=key; tv[p]=((u64)S[0][a].i<<32)|(u64)S[2][c].i; filled++; npairs_ac++;
          if(filled*2>cap){ overflow++; goto next_a; } } }
      next_a:;
    }
    // (B,D) pairs with key' = hT - hB - hD in the bucket: hD in [hT - hB - (low+W) + 1, hT - hB - low]
    for(i32 b=0;b<n[1];b++){
      u64 s=hT-S[1][b].h-(low+W)+1; if(W==0) s=0; i32 r0,r1,r2,r3; ranges(S[3],n[3],s,W,&r0,&r1,&r2,&r3);
      for(int pass=0;pass<2;pass++){
        i32 lo=pass?r2:r0, hi=pass?r3:r1;
        for(i32 d=lo;d<hi;d++){
          u64 key=hT-S[1][b].h-S[3][d].h; npairs_bd++;
          size_t p=(size_t)((key*0x9E3779B97F4A7C15ULL)>>(64-tb));
          while(tv[p]!=~0ULL){
            if(tk[p]==key){ nmatchkey++; u64 id=tv[p]; i32 a=(i32)(id>>32), c=(i32)(id&0xffffffffULL), bb=S[1][b].i, dd=S[3][d].i;
              int l=0; const i16*pa=V[0]+(size_t)a*m,*pb=V[1]+(size_t)bb*m,*pc=V[2]+(size_t)c*m,*pd=V[3]+(size_t)dd*m;
              while(l<m && pa[l]+pb[l]+pc[l]+pd[l]==T[l]) l++;
              if(l<K0) nfalse++; else { hist[l]++; if(l>=lstore) fprintf(fo,"Q %d %d %d %d %d\n",l,a,bb,c,dd); } }
            p=(p+1)&(cap-1); } } }
    }
    if(((beta-b0)&63)==63 || beta+1==b1){ fprintf(fo,"# bucket %llu/%llu done, %lds\n",(unsigned long long)(beta+1-b0),(unsigned long long)(b1-b0),(long)(time(0)-t0)); fflush(fo); }
  }
  fprintf(fo,"# STATS nAC_pairs %lld nBD_pairs %lld keymatches %lld false %lld overflow %lld\n",npairs_ac,npairs_bd,nmatchkey,nfalse,overflow);
  fprintf(fo,"HIST"); for(int l=K0;l<=m;l++) fprintf(fo," %lld",hist[l]); fprintf(fo,"\n"); fclose(fo); return 0;
}
