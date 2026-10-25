// Exact cover of the 78 pairs of points of Z_13 by component instances from insts.bin (inst.cpp): each pair covered
// once, each (point, D) used at most twice. Prints the first cover. Argument: bitmask of templates (1 = K4, 2 = prism,
// 4 = K33, 8 = Q3).
#include <bits/stdc++.h>
using namespace std;
typedef unsigned __int128 u128; typedef unsigned long long u64;
struct Inst{u128 mask;u64 usage;};
vector<Inst> rows; vector<array<uint8_t,26>> us; vector<int> rowtype;
vector<vector<int>> bypair(78);
long long nodes=0; int sols=0; int maxdepth=0;
u128 covered=0; int capr[26];
vector<int> chosen;
int NP=78;
static inline bool valid(int r){
  if(rows[r].mask&covered) return false;
  for(int k=0;k<26;k++) if(us[r][k]>capr[k]) return false;
  return true;
}
bool dfs(int depth){
  nodes++;
  if(depth>maxdepth){maxdepth=depth;}
  if((nodes&((1<<22)-1))==0) {fprintf(stderr,"nodes=%lld depth=%d maxdepth=%d\n",nodes,depth,maxdepth);}
  if(covered==(((u128)1<<78)-1)){ sols++; printf("SOLUTION with %zu components\n",chosen.size()); for(int r:chosen) printf(" type%d row%d\n",rowtype[r],r); fflush(stdout); return true; }
  // choose pair with min valid rows
  int bestp=-1; int bestc=INT_MAX; 
  for(int p=0;p<NP;p++){
    if((covered>>p)&1) continue;
    int c=0;
    for(int r:bypair[p]) { if(valid(r)){ c++; if(c>=bestc) break; } }
    if(c<bestc){bestc=c;bestp=p; if(c==0) return false;}
  }
  for(int r:bypair[bestp]){
    if(!valid(r)) continue;
    covered|=rows[r].mask; for(int k=0;k<26;k++) capr[k]-=us[r][k]; chosen.push_back(r);
    bool f=dfs(depth+1);
    chosen.pop_back(); covered&=~rows[r].mask; for(int k=0;k<26;k++) capr[k]+=us[r][k];
    if(f) return true;
  }
  return false;
}
int main(int argc,char**argv){
  int use=argc>1?atoi(argv[1]):15; // bitmask of types to include: 1=K4,2=prism,4=K33,8=Q3
  FILE*f=fopen("insts.bin","rb");
  for(int t=0;t<4;t++){ uint32_t n; if(fread(&n,4,1,f)!=1) return 1; vector<Inst> v(n); if(fread(v.data(),sizeof(Inst),n,f)!=n) return 1; if(!((use>>t)&1)) continue; for(auto&x:v){rows.push_back(x);rowtype.push_back(t);} }
  fclose(f);
  us.resize(rows.size());
  for(size_t r=0;r<rows.size();r++){ u64 u=rows[r].usage; for(int k=0;k<26;k++){us[r][k]=u%3;u/=3;} for(int p=0;p<NP;p++) if((rows[r].mask>>p)&1) bypair[p].push_back(r); }
  fprintf(stderr,"rows %zu\n",rows.size());
  for(int k=0;k<26;k++)capr[k]=2;
  bool r=dfs(0);
  printf("done found=%d nodes=%lld maxdepth=%d\n",(int)r,nodes,maxdepth);
}
