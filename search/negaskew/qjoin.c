// qjoin.c: exact 4-list join for quaternary skew GS over an orbit shape with EXPLICIT lists (one binary file of packed exponent
// vectors per block; 2 bits per orbit, orbit i at bits 2i..2i+1, n <= 32).  Blocks: list 0 = A (twisted skew, given explicitly),
// lists 1..3 = B, C, D.  Condition: sum of the 4 complex PAF vectors over classes 1..n-1 = 0 (exact, verified on every hit).
// Pairing: table from lists (p0,p1), probes from lists (p2,p3), chosen by the driver as argv.  DC (row sum) buckets: sum |R|^2 = 4v.
// usage: qjoin shape.txt l0.bin l1.bin l2.bin l3.bin p0 p1 p2 p3 [maxprint] [maxtab]
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
typedef unsigned long long u64;
static int v,n,d; static int sz[64],neg_[64]; static int *Nt; static u64 P[128];
typedef struct { u64 key; int rre,rim; u64 ex; } El;
static inline int m4(int x){return x&3;}
static void paf(const int*e,int*out){
  for(int c=1;c<n;c++){ int M[4]={0,0,0,0}; const int*N=Nt+(size_t)c*n*n;
    for(int a=0;a<n;a++){ const int*row=N+a*n; int ea=e[a]; for(int b=0;b<n;b++){ int x=row[b]; if(x) M[(ea-e[b])&3]+=x; } }
    out[2*(c-1)]=M[0]-M[2]; out[2*(c-1)+1]=M[1]-M[3]; }
}
static void fill(El*E,u64 ex){
  int e[64]; for(int i=0;i<n;i++) e[i]=(ex>>(2*i))&3;
  int vec[256]; paf(e,vec); u64 k=0; for(int i=0;i<2*d;i++) k+=(u64)(long long)vec[i]*P[i];
  int rr=0,ri=0; for(int a=0;a<n;a++){ switch(e[a]){case 0:rr+=sz[a];break;case 1:ri+=sz[a];break;case 2:rr-=sz[a];break;case 3:ri-=sz[a];break;} }
  E->key=k;E->rre=rr;E->rim=ri;E->ex=ex;
}
static int cmpR(const void*x,const void*y){const El*a=x,*b=y; if(a->rre!=b->rre) return a->rre<b->rre?-1:1; if(a->rim!=b->rim) return a->rim<b->rim?-1:1; return 0;}
typedef struct {int rre,rim; long lo,hi;} Bk;
static Bk* buckets(El*E,long m,long*nb){ Bk*b=malloc(sizeof(Bk)*(m+1)); long k=0; for(long i=0;i<m;){ long j=i; while(j<m&&E[j].rre==E[i].rre&&E[j].rim==E[i].rim) j++; b[k].rre=E[i].rre;b[k].rim=E[i].rim;b[k].lo=i;b[k].hi=j;k++; i=j;} *nb=k; return b;}
typedef struct { u64 key; unsigned int ia, ib; } Ent;
static El* load(const char*fn,long*m){
  FILE*f=fopen(fn,"rb"); if(!f){fprintf(stderr,"cannot open %s\n",fn); exit(1);} fseek(f,0,SEEK_END); long sz_=ftell(f); fseek(f,0,SEEK_SET);
  long cnt=sz_/8; u64*raw=malloc(sz_); if(fread(raw,8,cnt,f)!=(size_t)cnt) exit(1); fclose(f);
  El*E=malloc(sizeof(El)*cnt); for(long i=0;i<cnt;i++) fill(&E[i],raw[i]); free(raw); *m=cnt; return E;
}
int main(int argc,char**argv){
  FILE*f=fopen(argv[1],"r"); if(fscanf(f,"%d %d",&v,&n)!=2) return 1; d=n-1;
  for(int i=0;i<n;i++) if(fscanf(f,"%d",&sz[i])!=1) return 1;
  for(int i=0;i<n;i++) if(fscanf(f,"%d",&neg_[i])!=1) return 1;
  Nt=malloc(sizeof(int)*(size_t)n*n*n);
  for(int c=1;c<n;c++) for(int a=0;a<n;a++) for(int b=0;b<n;b++) if(fscanf(f,"%d",&Nt[(size_t)c*n*n+a*n+b])!=1) return 1;
  fclose(f);
  u64 s=88172645463325252ULL; for(int i=0;i<128;i++){ s^=s<<13;s^=s>>7;s^=s<<17; P[i]=s|1ULL; }
  El*Ls[4]; long Nl[4]; for(int i=0;i<4;i++){ Ls[i]=load(argv[2+i],&Nl[i]); qsort(Ls[i],Nl[i],sizeof(El),cmpR); }
  int p0=atoi(argv[6]),p1=atoi(argv[7]),p2=atoi(argv[8]),p3=atoi(argv[9]); int maxprint=argc>10?atoi(argv[10]):5; long maxtab=argc>11?atol(argv[11]):40000000L;
  Bk*Bs[4]; long Nb[4]; for(int i=0;i<4;i++) Bs[i]=buckets(Ls[i],Nl[i],&Nb[i]);
  fprintf(stderr,"n=%d sizes %ld %ld %ld %ld ; buckets %ld %ld %ld %ld ; pairing (%d,%d | %d,%d)\n",n,Nl[0],Nl[1],Nl[2],Nl[3],Nb[0],Nb[1],Nb[2],Nb[3],p0,p1,p2,p3);
  long total=4L*v; long long count=0; int printed=0; long long probes=0, tabs=0;
  El*X=Ls[p0],*Y=Ls[p1],*Z=Ls[p2],*W=Ls[p3]; Bk*bX=Bs[p0],*bY=Bs[p1],*bZ=Bs[p2],*bW=Bs[p3]; long nX=Nb[p0],nY=Nb[p1],nZ=Nb[p2],nW=Nb[p3];
  for(long ix=0;ix<nX;ix++) for(long iy=0;iy<nY;iy++){
    long nx_=bX[ix].hi-bX[ix].lo, ny_=bY[iy].hi-bY[iy].lo;
    long n1=(long)bX[ix].rre*bX[ix].rre+(long)bX[ix].rim*bX[ix].rim, n2=(long)bY[iy].rre*bY[iy].rre+(long)bY[iy].rim*bY[iy].rim;
    if(n1+n2>total) continue;
    // any (z,w) bucket pair valid?
    long chunk=nx_; if(nx_*ny_>maxtab){ chunk=maxtab/ny_; if(chunk<1) chunk=1; }
    int built=0;
    for(long xlo=bX[ix].lo; xlo<bX[ix].hi; xlo+=chunk){
      long xhi=xlo+chunk; if(xhi>bX[ix].hi) xhi=bX[ix].hi; long nlc=(xhi-xlo)*ny_;
      long tsz=1; while(tsz<2*nlc) tsz<<=1; Ent*tab=NULL; built=0;
      for(long iz=0;iz<nZ;iz++){
        long n3=(long)bZ[iz].rre*bZ[iz].rre+(long)bZ[iz].rim*bZ[iz].rim; if(n1+n2+n3>total) continue;
        for(long iw=0;iw<nW;iw++){
          long n4=(long)bW[iw].rre*bW[iw].rre+(long)bW[iw].rim*bW[iw].rim; if(n1+n2+n3+n4!=total) continue;
          if(!built){ tab=calloc(tsz,sizeof(Ent)); if(!tab){ fprintf(stderr,"calloc fail\n"); return 2; } tabs++;
            for(long x=xlo;x<xhi;x++) for(long y=bY[iy].lo;y<bY[iy].hi;y++){ u64 k=X[x].key+Y[y].key; u64 h=((k*0x9E3779B97F4A7C15ULL)>>20)&(tsz-1); while(tab[h].ia|tab[h].ib) h=(h+1)&(tsz-1); tab[h].key=k;tab[h].ia=(unsigned)x+1;tab[h].ib=(unsigned)y+1; }
            built=1; }
          for(long z=bZ[iz].lo;z<bZ[iz].hi;z++) for(long w=bW[iw].lo;w<bW[iw].hi;w++){
            probes++;
            u64 k=-(Z[z].key+W[w].key); u64 h=((k*0x9E3779B97F4A7C15ULL)>>20)&(tsz-1);
            while(tab[h].ia|tab[h].ib){ if(tab[h].key==k){
                int ee[4][64]; u64 exs[4]; exs[p0]=X[tab[h].ia-1].ex; exs[p1]=Y[tab[h].ib-1].ex; exs[p2]=Z[z].ex; exs[p3]=W[w].ex;
                int vv[4][256]; int ok=1; for(int t=0;t<4;t++){ for(int i=0;i<n;i++) ee[t][i]=(exs[t]>>(2*i))&3; paf(ee[t],vv[t]); }
                for(int i=0;i<2*d;i++) if(vv[0][i]+vv[1][i]+vv[2][i]+vv[3][i]){ok=0;break;}
                if(ok){ count++; if(printed<maxprint){ printed++; printf("SOL"); for(int t=0;t<4;t++){ printf(" "); for(int i=0;i<n;i++) printf("%d",ee[t][i]); } printf("\n"); fflush(stdout);} }
              }
              h=(h+1)&(tsz-1); }
          }
        }
      }
      if(built) free(tab);
    }
  }
  printf("COUNT %lld  (tables %lld probes %lld)\n",count,tabs,probes);
  return 0;
}
