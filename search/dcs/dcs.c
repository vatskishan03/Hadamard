/* Palindromic-type DCS search.  a real +-1; b over {+-1,+-i,+-s}, type pattern symmetric:
   i eps-palindromic on Omega_i, u arbitrary on Omega_u (symmetric), b_j = u_j (real) if u_j u_{n-1-j} = eps else u_j s.
   Then the DCS conditions reduce to N_a + N_u + N_i = 2n delta (aperiodic).  Meet in the middle on N_a.
   usage: dcs n eps [maxprint]  */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
typedef unsigned long long u64;
static int n, eps, h;
static u64 *hk; static uint32_t *ha; static u64 hmask;
static u64 hashv(const int *T){ u64 x=1469598103934665603ULL; for(int k=1;k<n;k++){ x^=(u64)(T[k]+64); x*=1099511628211ULL; x^=x>>29;} return x|1; }
static int issq[4096];
int main(int argc,char**argv){
  n=atoi(argv[1]); eps=atoi(argv[2]); long maxprint=argc>3?atol(argv[3]):20; h=(n-1)/2;
  for(int x=0;x*x<4096;x++) issq[x*x]=1;
  /* phase 1 */
  u64 cap=1; while(cap < (2ULL<<(n-1))) cap<<=1; hmask=cap-1;
  hk=calloc(cap,8); ha=calloc(cap,4); if(!hk||!ha){fprintf(stderr,"oom\n");return 1;}
  long ins=0; int a[40],T[40];
  for(uint32_t bits=0; bits < (1u<<(n-1)); bits++){
    a[0]=1; int sa=1; for(int j=1;j<n;j++){ a[j]=((bits>>(j-1))&1)?-1:1; sa+=a[j]; }
    if(sa*sa>2*n) continue;
    for(int k=1;k<n;k++){ int s=0; for(int j=0;j+k<n;j++) s+=a[j]*a[j+k]; T[k]=s; }
    u64 key=hashv(T), p=key&hmask;
    while(hk[p] && hk[p]!=key) p=(p+1)&hmask;
    if(!hk[p]){ hk[p]=key; ha[p]=bits; ins++; }
  }
  fprintf(stderr,"n=%d eps=%d table %ld distinct\n",n,eps,ins);
  long nsol=0; double nlook=0, nstep=0;
  int u[40], iv[40], Nu[40], Ni[40], pos[40];
  for(uint32_t P=0; P < (1u<<h); P++){
    for(int mid=0; mid<=(eps==1); mid++){
      int npi=__builtin_popcount(P), nis=npi+mid;
      int nu=0; for(int j=0;j<n;j++){ int inI = (j<h)? ((P>>j)&1) : (j>h? ((P>>(n-1-j))&1) : mid); if(!inI) pos[nu++]=j; }
      if(nis==0 && nu==0) continue;
      int ifirst = (nis>0);           /* fix first i sign */
      for(uint32_t is=0; is < (nis? (1u<<(nis-1)) : 1u); is++){
        memset(iv,0,sizeof iv); int c=0, si=0;
        for(int j=0;j<h;j++) if((P>>j)&1){ int sg = (c==0)?1:(((is>>(c-1))&1)?-1:1); iv[j]=sg; iv[n-1-j]=eps*sg; c++; }
        if(mid){ int sg=(c==0)?1:(((is>>(c-1))&1)?-1:1); iv[h]=sg; c++; }
        for(int j=0;j<n;j++) si+=iv[j];
        for(int k=1;k<n;k++){ int s=0; for(int j=0;j+k<n;j++) s+=iv[j]*iv[j+k]; Ni[k]=s; }
        /* u over pos[0..nu-1], first fixed +1, Gray code over the remaining nu-1 */
        memset(u,0,sizeof u); int su=0; for(int t=0;t<nu;t++){ u[pos[t]]=1; su++; }
        for(int k=1;k<n;k++){ int s=0; for(int j=0;j+k<n;j++) s+=u[j]*u[j+k]; Nu[k]=s; }
        u64 G = nu? (1ULL<<(nu-1)) : 1;
        for(u64 g=0; g<G; g++){
          if(g){ int t=__builtin_ctzll(g)+1; int j=pos[t]; int uj=u[j];
            for(int k=1;k<n;k++){ int d=0; if(j+k<n) d+=u[j+k]; if(j-k>=0) d+=u[j-k]; Nu[k]-=2*uj*d; }
            u[j]=-uj; su-=2*uj; }
          nstep++;
          int r=2*n-su*su-si*si; if(r<0||!issq[r]) continue;
          if(((int)( (r>0)?1:0 ))==0 && (n&1)) continue; /* sa odd for odd n */
          for(int k=1;k<n;k++) T[k]=-(Nu[k]+Ni[k]);
          nlook++;
          u64 key=hashv(T), p=key&hmask;
          while(hk[p] && hk[p]!=key) p=(p+1)&hmask;
          if(hk[p]==key){
            uint32_t bits=ha[p]; a[0]=1; for(int j=1;j<n;j++) a[j]=((bits>>(j-1))&1)?-1:1;
            int ok=1; for(int k=1;k<n&&ok;k++){ int s=0; for(int j=0;j+k<n;j++) s+=a[j]*a[j+k]; if(s!=T[k]) ok=0; }
            if(!ok) continue;
            nsol++;
            if(nsol<=maxprint){
              printf("a=");
              for(int j=0;j<n;j++) printf("%s%d",j?",":"",a[j]);
              printf(" b=");
              for(int j=0;j<n;j++){ const char*sg; if(iv[j]){ printf("%s%si",j?",":"",iv[j]>0?"":"-"); }
                else { printf("%s%s%s",j?",":"",u[j]>0?"":"-", (u[j]*u[n-1-j]==eps)?"1":"s"); } }
              printf("\n"); fflush(stdout);
            }
          }
        }
      }
    }
  }
  fprintf(stderr,"n=%d eps=%d steps %.3g lookups %.3g solutions %ld\n",n,eps,nstep,nlook,nsol);
  printf("# n=%d eps=%d solutions %ld (up to sign of u and i; a mod negation and reversal)\n",n,eps,nsol);
  return 0;
}
