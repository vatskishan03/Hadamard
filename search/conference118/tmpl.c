// For a cubic template graph on n vertices (edge list on the command line), enumerate D in {+-1}^n with sum 0, edge
// signs s_e and phases c_e in Z_3 (tree edges c = 0), and test W = D + 2 Phi, Phi_uv = s_e w^(c_e), W^2 = 13 I exactly
// over Z[w]. Prints the number of solutions and the solvable classes (S, D).
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef struct {int x,y;} Z; // x + y w, w^2 = -1 - w
static Z zmul(Z a,Z b){Z r={a.x*b.x-a.y*b.y, a.x*b.y+a.y*b.x-a.y*b.y};return r;}
static Z zconj(Z a){Z r={a.x-a.y,-a.y};return r;}
static Z WP[3]={{1,0},{0,1},{-1,-1}};
int n,m; int eu[16],ev[16]; int tree[16];
int main(int argc,char**argv){
  n=atoi(argv[1]); m=3*n/2; for(int e=0;e<m;e++){eu[e]=atoi(argv[2+2*e]);ev[e]=atoi(argv[3+2*e]);}
  // spanning tree by union-find
  int par[16]; for(int i=0;i<n;i++)par[i]=i;
  #define FIND(x) ({int _x=(x); while(par[_x]!=_x)_x=par[_x]; _x;})
  int nt=0; for(int e=0;e<m;e++){int a=FIND(eu[e]),b=FIND(ev[e]); if(a!=b){par[a]=b;tree[e]=1;nt++;}else tree[e]=0;}
  int fe[16],nf=0; for(int e=0;e<m;e++) if(!tree[e]) fe[nf++]=e;
  long long total=0; long long sols=0;
  // adjacency lookup
  int adj[16][16]; memset(adj,-1,sizeof(adj)); for(int e=0;e<m;e++){adj[eu[e]][ev[e]]=e;adj[ev[e]][eu[e]]=e;}
  int nsp=1<<m; int nph=1; for(int i=0;i<nf;i++) nph*=3;
  // signs pattern histogram of solvable (D,s): record sign-vector + D mask
  static char seen[1<<12][1<<8]; memset(seen,0,sizeof(seen));
  int nsolvable=0;
  for(int sm=0;sm<nsp;sm++){
    for(int dm=0;dm<(1<<n);dm++){
      if(__builtin_popcount(dm)!=n/2) continue;
      int D[16]; for(int i=0;i<n;i++) D[i]=(dm>>i&1)?1:-1;
      for(int ph=0;ph<nph;ph++){
        int c[16]; for(int e=0;e<m;e++)c[e]=0; int t=ph; for(int i=0;i<nf;i++){c[fe[i]]=t%3;t/=3;}
        // build Phi_ij (i->j) = s_e w^{c_e} for i=eu, conj for reverse
        Z Phi[16][16]; memset(Phi,0,sizeof(Phi));
        for(int e=0;e<m;e++){int s=(sm>>e&1)?-1:1; Z v={WP[c[e]].x*s,WP[c[e]].y*s}; Phi[eu[e]][ev[e]]=v; Phi[ev[e]][eu[e]]=zconj(v);}
        int ok=1;
        for(int i=0;i<n&&ok;i++)for(int j=i+1;j<n&&ok;j++){
          // (W^2)_ij = 2(Di+Dj)Phi_ij + 4 sum_k Phi_ik Phi_kj
          Z s={0,0};
          for(int k=0;k<n;k++){ if(k==i||k==j) continue; Z pr=zmul(Phi[i][k],Phi[k][j]); s.x+=pr.x;s.y+=pr.y; }
          Z q={2*(D[i]+D[j])*Phi[i][j].x+4*s.x, 2*(D[i]+D[j])*Phi[i][j].y+4*s.y};
          if(q.x||q.y) ok=0;
        }
        total++;
        if(ok){sols++; if(!seen[sm][dm]){seen[sm][dm]=1;nsolvable++;}}
      }
    }
  }
  printf("n=%d tested=%lld solutions=%lld distinct (signs,D) classes solvable=%d\n",n,total,sols,nsolvable);
  // print the solvable classes
  int cnt=0;
  for(int sm=0;sm<nsp;sm++)for(int dm=0;dm<(1<<n);dm++) if(seen[sm][dm]){ if(cnt<400){ printf("S="); for(int e=0;e<m;e++)printf("%c",(sm>>e&1)?'-':'+'); printf(" D="); for(int i=0;i<n;i++)printf("%c",(dm>>i&1)?'+':'-'); printf("\n");} cnt++; }
  return 0;
}
