// qspin_fast.c: quaternary SPIN skew GS over Z_v: C = B o mu, D = B o mu^2 (mu^3 in H), A twisted-skew; H-orbit level; up to n = 21 orbits.
// Same search as qspin.c (exhaustive, exact) with three speed-ups: the row-sum condition |R_A|^2 + 3 |R_B|^2 = 4v is met in the middle
// over the low and high digits of B (no 4^(n-1) enumeration), a two-class presence filter against the A list, and a bilinear fast
// path for the class sums of B (about 10 table lookups per B that passes the row sums).  Counts and statistics equal those of qspin.c.
// input (spin_drv.spin_shape): v n ; sz ; neg ; cperm[n] (class permutation induced by mu, classes 1..n-1, cperm[0]=0) ; operm[n] (orbit permutation) ; N tensor (n-1 classes)
// usage: qspin_fast file maxprint [part nparts]      (part/nparts: split the low-digit buckets over independent runs; the COUNT lines add up)
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <time.h>
static double now(){ struct timespec ts; clock_gettime(CLOCK_MONOTONIC,&ts); return ts.tv_sec+ts.tv_nsec*1e-9; }
static double T0;
typedef unsigned long long u64;
static int v,n,d; static int sz[64],neg_[64],cperm[64],operm[64]; static int *Nt; static u64 P[128];
static inline int m4(int x){return x&3;}
static void pafc(const int*e,int c,int*re,int*im){
  int M[4]={0,0,0,0}; const int*N=Nt+(size_t)c*n*n;
  for(int a=0;a<n;a++){ const int*row=N+a*n; int ea=e[a]; for(int b=0;b<n;b++){ int x=row[b]; if(x) M[(ea-e[b])&3]+=x; } }
  *re=M[0]-M[2]; *im=M[1]-M[3];
}
typedef struct{ u64 key; u64 ex; int rre,rim; } El;
typedef struct { int rr,ri; u64 ex; } Hf;
static int cmpH(const void*x,const void*y){ const Hf*a=x,*b=y; if(a->rr!=b->rr) return a->rr<b->rr?-1:1; if(a->ri!=b->ri) return a->ri<b->ri?-1:1; return 0; }
int main(int argc,char**argv){ T0=now();
  FILE*f=fopen(argv[1],"r"); int maxprint=argc>2?atoi(argv[2]):3;
  if(fscanf(f,"%d %d",&v,&n)!=2) return 1; d=n-1;
  for(int i=0;i<n;i++) if(fscanf(f,"%d",&sz[i])!=1) return 1;
  for(int i=0;i<n;i++) if(fscanf(f,"%d",&neg_[i])!=1) return 1;
  for(int i=0;i<n;i++) if(fscanf(f,"%d",&cperm[i])!=1) return 1;
  for(int i=0;i<n;i++) if(fscanf(f,"%d",&operm[i])!=1) return 1;
  Nt=malloc(sizeof(int)*(size_t)n*n*n);
  for(int c=1;c<n;c++) for(int a=0;a<n;a++) for(int b=0;b<n;b++) if(fscanf(f,"%d",&Nt[(size_t)c*n*n+a*n+b])!=1) return 1;
  fclose(f);
  u64 s=88172645463325252ULL; for(int i=0;i<128;i++){ s^=s<<13;s^=s>>7;s^=s<<17; P[i]=s|1ULL; }
  // A list
  int pr[64],np=0; for(int a=1;a<n;a++) if(a<neg_[a]) pr[np++]=a;
  long nA=2; for(int i=0;i<np;i++) nA*=4;
  El*A=malloc(sizeof(El)*nA);
  static char okR[1<<24]; memset(okR,0,sizeof(okR)); // |R_A|^2 values <= 4v
  long total=4L*v;
  // first-class bitmap of A's PAF at class 1 : we use value (re,im) offset by v
  int W=2*v+1; char*bm=calloc((size_t)W*W,1);
  int vec[256];
  for(long idx=0;idx<nA;idx++){
    long t=idx; int o0=t%2; t/=2; int e[64]; e[0]=o0?3:0;
    for(int i=0;i<np;i++){ int ej=t%4; t/=4; e[pr[i]]=ej; e[neg_[pr[i]]]=m4(1-ej); }
    u64 k=0; for(int c=1;c<n;c++){ int re,im; pafc(e,c,&re,&im); k+=(u64)(long long)re*P[2*(c-1)]+(u64)(long long)im*P[2*(c-1)+1]; if(c==1){ bm[(size_t)(-im+v)*W+(-re+v)]=1; } }
    int rr=0,ri=0; u64 ex=0; for(int a=0;a<n;a++){ int tt=m4(e[a]); ex|=((u64)tt)<<(2*a); switch(tt){case 0:rr+=sz[a];break;case 1:ri+=sz[a];break;case 2:rr-=sz[a];break;case 3:ri-=sz[a];break;} }
    A[idx].key=k;A[idx].ex=ex;A[idx].rre=rr;A[idx].rim=ri;
    long nn=(long)rr*rr+(long)ri*ri; if(nn<=total) okR[nn]=1;
  }
  // second class (different pi-orbit) for a 2-class presence filter
  int cb=0; { int c1_=1,c2_=cperm[1],c3_=cperm[c2_]; for(int c=1;c<n;c++) if(c!=c1_&&c!=c2_&&c!=c3_){ cb=c; break; } }
  int cb2=cperm[cb], cb3=cperm[cb2];
  long ps=1; while(ps<4*nA) ps<<=1; u64*pset=calloc(ps,sizeof(u64));
  for(long idx=0;idx<nA;idx++){
    int ea[64]; for(int i=0;i<n;i++) ea[i]=(A[idx].ex>>(2*i))&3;
    int r1,i1,r2,i2; pafc(ea,1,&r1,&i1); pafc(ea,cb,&r2,&i2);
    u64 k2=((u64)(unsigned)(-r1)&0xFFFF)|(((u64)(unsigned)(-i1)&0xFFFF)<<16)|(((u64)(unsigned)(-r2)&0xFFFF)<<32)|(((u64)(unsigned)(-i2)&0xFFFF)<<48);
    k2|=1ULL<<63; u64 hq=(k2*0x9E3779B97F4A7C15ULL)>>20; hq&=ps-1; while(pset[hq]&&pset[hq]!=k2) hq=(hq+1)&(ps-1); pset[hq]=k2;
  }
  long tsz=1; while(tsz<2*nA) tsz<<=1; unsigned *tab=calloc(tsz,sizeof(unsigned));
  for(long i=0;i<nA;i++){ u64 h=(A[i].key*0x9E3779B97F4A7C15ULL)>>20; h&=tsz-1; while(tab[h]) h=(h+1)&(tsz-1); tab[h]=(unsigned)i+1; }
  fprintf(stderr,"n=%d |A|=%ld  t=%.1f\n",n,nA,now()-T0);
  // meet-in-the-middle over the DC condition: B = (lo digits 1..h, hi digits h+1..n-1); only (lo,hi) bucket pairs whose R sum passes DC are processed
  int part=argc>3?atoi(argv[3]):0, nparts=argc>4?atoi(argv[4]):1;
  int h=(n-1)/2, g=n-1-h; // lo has h digits (orbits 1..h), hi has g digits (orbits h+1..n-1)
  long nlo=1, nhi=1; for(int i=0;i<h;i++) nlo*=4; for(int i=0;i<g;i++) nhi*=4;
  Hf*LO=malloc(sizeof(Hf)*nlo), *HI=malloc(sizeof(Hf)*nhi);
  static const int UX[4]={1,0,-1,0}, UY[4]={0,1,0,-1};
  for(long idx=0;idx<nlo;idx++){ long t=idx; int rr=0,ri=0; u64 ex=0; for(int i=1;i<=h;i++){ int dgt=t&3; t>>=2; ex|=((u64)dgt)<<(2*i); rr+=sz[i]*UX[dgt]; ri+=sz[i]*UY[dgt]; } LO[idx].rr=rr;LO[idx].ri=ri;LO[idx].ex=ex; }
  for(long idx=0;idx<nhi;idx++){ long t=idx; int rr=0,ri=0; u64 ex=0; for(int i=h+1;i<=n-1;i++){ int dgt=t&3; t>>=2; ex|=((u64)dgt)<<(2*i); rr+=sz[i]*UX[dgt]; ri+=sz[i]*UY[dgt]; } HI[idx].rr=rr;HI[idx].ri=ri;HI[idx].ex=ex; }
  qsort(LO,nlo,sizeof(Hf),cmpH); qsort(HI,nhi,sizeof(Hf),cmpH);
  // buckets
  long *loS=malloc(sizeof(long)*(nlo+1)), nlb=0; for(long i=0;i<nlo;){ long j=i; while(j<nlo&&LO[j].rr==LO[i].rr&&LO[j].ri==LO[i].ri) j++; loS[nlb++]=i; i=j; } loS[nlb]=nlo;
  long *hiS=malloc(sizeof(long)*(nhi+1)), nhb=0; for(long i=0;i<nhi;){ long j=i; while(j<nhi&&HI[j].rr==HI[i].rr&&HI[j].ri==HI[i].ri) j++; hiS[nhb++]=i; i=j; } hiS[nhb]=nhi;
  int e[64]; for(int i=0;i<n;i++) e[i]=0;
  long long count=0; int printed=0; long long passDC=0,pass1=0;
  int c1=1,c2=cperm[1],c3=cperm[c2];
  int Wd=2*v+1; char*okRB=calloc((size_t)Wd*Wd,1);
  for(int x=-v;x<=v;x++) for(int y=-v;y<=v;y++){ long nn=3L*((long)x*x+(long)y*y); if(nn<=total && okR[total-nn]) okRB[(size_t)(x+v)*Wd+(y+v)]=1; }
  fprintf(stderr,"lo buckets %ld hi buckets %ld  t=%.1f\n",nlb,nhb,now()-T0);
  fprintf(stderr,"setup done t=%.1f\n",now()-T0);
  // ---- FAST PATH: bilinear decomposition of the class-sum PAF over (lo, hi) digit halves ----
  int cls3[2][3]={{c1,c2,c3},{cb,cb2,cb3}};
  int *HHre[2],*HHim[2]; for(int gi=0;gi<2;gi++){ HHre[gi]=malloc(sizeof(int)*nhi); HHim[gi]=malloc(sizeof(int)*nhi); } unsigned char *HIdig=malloc((size_t)nhi*g);
  for(long ih=0;ih<nhi;ih++){
    u64 exh=HI[ih].ex; int eh[64]; for(int b=0;b<g;b++){ int dg=(exh>>(2*(h+1+b)))&3; HIdig[(size_t)ih*g+b]=dg; eh[b]=dg; }
    for(int gi=0;gi<2;gi++){ int sr=0,si=0;
      for(int ci=0;ci<3;ci++){ const int*N=Nt+(size_t)cls3[gi][ci]*n*n; int M[4]={0,0,0,0};
        for(int a=0;a<g;a++) for(int b=0;b<g;b++){ int x=N[(h+1+a)*n+(h+1+b)]; if(x) M[(eh[a]-eh[b])&3]+=x; }
        sr+=M[0]-M[2]; si+=M[1]-M[3]; }
      HHre[gi][ih]=sr; HHim[gi][ih]=si; }
  }
  int nc1=n-1; int *HHall=malloc(sizeof(int)*(size_t)nhi*nc1*2);
  for(long ih=0;ih<nhi;ih++){ u64 exh=HI[ih].ex; int eh[64]; for(int b=0;b<g;b++) eh[b]=(exh>>(2*(h+1+b)))&3;
    for(int c=1;c<n;c++){ const int*N=Nt+(size_t)c*n*n; int M[4]={0,0,0,0};
      for(int a=0;a<g;a++) for(int b=0;b<g;b++){ int x=N[(h+1+a)*n+(h+1+b)]; if(x) M[(eh[a]-eh[b])&3]+=x; }
      HHall[((size_t)ih*nc1+(c-1))*2]=M[0]-M[2]; HHall[((size_t)ih*nc1+(c-1))*2+1]=M[1]-M[3]; } }
  static const int RX[4]={1,0,-1,0}, RY[4]={0,1,0,-1};
  long *vbh=malloc(sizeof(long)*(nhb+1));
  for(long bl=0;bl<nlb;bl++){
    if(bl%nparts!=part) continue;
    int lrr=LO[loS[bl]].rr, lri=LO[loS[bl]].ri;
    long nv=0;
    for(long bh=0;bh<nhb;bh++){
      int rr=lrr+HI[hiS[bh]].rr+sz[0], ri=lri+HI[hiS[bh]].ri;
      if(rr<-v||rr>v||ri<-v||ri>v) continue;
      if(!okRB[(size_t)(rr+v)*Wd+(ri+v)]) continue;
      vbh[nv++]=bh;
    }
    if(nv==0) continue;
    for(long il=loS[bl];il<loS[bl+1];il++){
      int el[64]; el[0]=0; u64 exl=LO[il].ex; for(int i=1;i<=h;i++) el[i]=(exl>>(2*i))&3;
      int base_re[2]={0,0}, base_im[2]={0,0};
      for(int gi=0;gi<2;gi++) for(int ci=0;ci<3;ci++){ const int*N=Nt+(size_t)cls3[gi][ci]*n*n; int M[4]={0,0,0,0};
        for(int a=0;a<=h;a++) for(int b=0;b<=h;b++){ int x=N[a*n+b]; if(x) M[(el[a]-el[b])&3]+=x; }
        base_re[gi]+=M[0]-M[2]; base_im[gi]+=M[1]-M[3]; }
      int T[2][32][4][2];
      int tallBuilt=0; static int Tall[64][32][4][2]; static int baseall[64][2];
      for(int gi=0;gi<2;gi++) for(int bb=0;bb<g;bb++){
        int ob=h+1+bb; int Ure=0,Uim=0,Vre=0,Vim=0;
        for(int ci=0;ci<3;ci++){ const int*N=Nt+(size_t)cls3[gi][ci]*n*n;
          for(int a=0;a<=h;a++){ int n1=N[a*n+ob], n2=N[ob*n+a];
            if(n1){ Ure+=n1*RX[el[a]]; Uim+=n1*RY[el[a]]; }
            if(n2){ int k2=(4-el[a])&3; Vre+=n2*RX[k2]; Vim+=n2*RY[k2]; } } }
        for(int d=0;d<4;d++){
          int kd=(4-d)&3; int ar=Ure*RX[kd]-Uim*RY[kd], ai=Ure*RY[kd]+Uim*RX[kd];
          int br=Vre*RX[d]-Vim*RY[d], bi=Vre*RY[d]+Vim*RX[d];
          T[gi][bb][d][0]=ar+br; T[gi][bb][d][1]=ai+bi; }
      }
      for(long iv=0;iv<nv;iv++){ long bh=vbh[iv];
        passDC+=hiS[bh+1]-hiS[bh];
        for(long ih=hiS[bh];ih<hiS[bh+1];ih++){
          int sre=base_re[0]+HHre[0][ih], sim=base_im[0]+HHim[0][ih]; const unsigned char*dg=HIdig+(size_t)ih*g;
          for(int bb=0;bb<g;bb++){ sre+=T[0][bb][dg[bb]][0]; sim+=T[0][bb][dg[bb]][1]; }
          if(sre<-v||sre>v||sim<-v||sim>v) continue;
          if(!bm[(size_t)(sim+v)*W+(sre+v)]) continue;
          int s2r=base_re[1]+HHre[1][ih], s2i=base_im[1]+HHim[1][ih];
          for(int bb=0;bb<g;bb++){ s2r+=T[1][bb][dg[bb]][0]; s2i+=T[1][bb][dg[bb]][1]; }
          { u64 k2=((u64)(unsigned)(sre)&0xFFFF)|(((u64)(unsigned)(sim)&0xFFFF)<<16)|(((u64)(unsigned)(s2r)&0xFFFF)<<32)|(((u64)(unsigned)(s2i)&0xFFFF)<<48);
            k2|=1ULL<<63; u64 hq=(k2*0x9E3779B97F4A7C15ULL)>>20; hq&=ps-1; int found=0; while(pset[hq]){ if(pset[hq]==k2){found=1;break;} hq=(hq+1)&(ps-1); }
            if(!found) continue; }
          u64 ex=LO[il].ex|HI[ih].ex; for(int i=0;i<n;i++) e[i]=(ex>>(2*i))&3;
          pass1++;
          if(!tallBuilt){ tallBuilt=1;
            for(int c=1;c<n;c++){ const int*N=Nt+(size_t)c*n*n; int M[4]={0,0,0,0};
              for(int a=0;a<=h;a++) for(int b=0;b<=h;b++){ int x=N[a*n+b]; if(x) M[(el[a]-el[b])&3]+=x; }
              baseall[c][0]=M[0]-M[2]; baseall[c][1]=M[1]-M[3];
              for(int bb=0;bb<g;bb++){ int ob=h+1+bb; int Ure=0,Uim=0,Vre=0,Vim=0;
                for(int a=0;a<=h;a++){ int n1=N[a*n+ob], n2=N[ob*n+a];
                  if(n1){ Ure+=n1*RX[el[a]]; Uim+=n1*RY[el[a]]; }
                  if(n2){ int k2=(4-el[a])&3; Vre+=n2*RX[k2]; Vim+=n2*RY[k2]; } }
                for(int d=0;d<4;d++){ int kd=(4-d)&3; int ar=Ure*RX[kd]-Uim*RY[kd], ai=Ure*RY[kd]+Uim*RX[kd]; int br=Vre*RX[d]-Vim*RY[d], bi=Vre*RY[d]+Vim*RX[d];
                  Tall[c][bb][d][0]=ar+br; Tall[c][bb][d][1]=ai+bi; } } } }
          int pv[64][2];
          for(int c=1;c<n;c++){ int pr=baseall[c][0]+HHall[((size_t)ih*nc1+(c-1))*2], pi=baseall[c][1]+HHall[((size_t)ih*nc1+(c-1))*2+1];
            for(int bb=0;bb<g;bb++){ pr+=Tall[c][bb][dg[bb]][0]; pi+=Tall[c][bb][dg[bb]][1]; }
            pv[c][0]=pr; pv[c][1]=pi; }
          u64 k=0; for(int c=1;c<n;c++){ int c2_=cperm[c],c3_=cperm[c2_]; int sr=pv[c][0]+pv[c2_][0]+pv[c3_][0], si=pv[c][1]+pv[c2_][1]+pv[c3_][1]; k+=(u64)(long long)sr*P[2*(c-1)]+(u64)(long long)si*P[2*(c-1)+1]; }
          u64 kk=-k; u64 hh=(kk*0x9E3779B97F4A7C15ULL)>>20; hh&=tsz-1;
          while(tab[hh]){ if(A[tab[hh]-1].key==kk){
              long ai=tab[hh]-1; int ea[64]; for(int i=0;i<n;i++) ea[i]=(A[ai].ex>>(2*i))&3;
              int ec[64],ed[64]; for(int i=0;i<n;i++){ ec[i]=e[operm[i]]; ed[i]=e[operm[operm[i]]]; }
              int ok=1; for(int c=1;c<n&&ok;c++){ int r0,i0,r1,i1,r2,i2,r3,i3; pafc(ea,c,&r0,&i0); pafc(e,c,&r1,&i1); pafc(ec,c,&r2,&i2); pafc(ed,c,&r3,&i3); if(r0+r1+r2+r3||i0+i1+i2+i3) ok=0; }
              if(ok){ count++; if(printed<maxprint){ printed++; printf("SOL A="); for(int i=0;i<n;i++) printf("%d",ea[i]); printf(" B="); for(int i=0;i<n;i++) printf("%d",e[i]); printf("\n"); fflush(stdout);} }
            }
            hh=(hh+1)&(tsz-1); }
        }
      }
    }
    if(bl%200==part) fprintf(stderr,"lo bucket %ld/%ld passDC %lld pass1 %lld count %lld\n",bl,nlb,passDC,pass1,count);
  }
  fprintf(stderr,"loop done t=%.1f\n",now()-T0);
  printf("COUNT %lld  (part %d/%d passDC %lld pass1 %lld)\n",count,part,nparts,passDC,pass1);
  return 0;
}
