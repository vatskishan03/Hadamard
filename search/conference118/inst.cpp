// Component instances: for each template class (S, D) of templates.txt, every placement of the template's vertices on
// points of Z_13 whose edge signs are the Paley signs chi(p_v - p_u), with distinct pairs and each (point, D) used at
// most twice; unique (pair mask, usage) written to insts.bin.
#include <bits/stdc++.h>
using namespace std;
typedef unsigned __int128 u128; typedef unsigned long long u64;
int P=13; int chi[64];
int pidx[13][13];
struct Tmpl{string name;int n,m;vector<pair<int,int>>E;vector<pair<string,string>>cls;};
vector<Tmpl> T;
struct Inst{u128 mask;u64 usage;};
struct KeyHash{size_t operator()(const pair<u128,u64>&k)const{return (size_t)(k.first^(k.first>>64)^(k.second*0x9E3779B97F4A7C15ULL));}};
unordered_set<pair<u128,u64>,KeyHash> seenset[4];
vector<Inst> insts[4];
int cur_n,cur_m; vector<pair<int,int>> curE; string curS,curD;
int pt[8]; int cnt[13][2];
u128 usedmask;
vector<vector<int>> nbrE; // edges incident
int ti;
void dfs(int v){
  if(v==cur_n){
    u64 us=0; for(int i=P-1;i>=0;i--)for(int d=1;d>=0;d--) us=us*3+cnt[i][d];
    pair<u128,u64> key(usedmask,us);
    if(seenset[ti].insert(key).second) insts[ti].push_back({usedmask,us});
    return;
  }
  int d=(curD[v]=='+')?0:1;
  for(int p=0;p<P;p++){
    if(cnt[p][d]>=2) continue;
    // check edges to earlier vertices
    u128 add=0; bool ok=true;
    for(int e: nbrE[v]){
      int a=curE[e].first,b=curE[e].second; int u=(a==v)?b:a; if(u>v) continue;
      if(pt[u]==p){ok=false;break;}
      int s=(curS[e]=='+')?1:-1;
      if(chi[((p-pt[u])%P+P)%P]!=s){ok=false;break;}
      int i=min(p,pt[u]),j=max(p,pt[u]); u128 bit=(u128)1<<pidx[i][j];
      if((usedmask|add)&bit){ok=false;break;}
      add|=bit;
    }
    if(!ok) continue;
    pt[v]=p; cnt[p][d]++; usedmask|=add;
    dfs(v+1);
    usedmask&=~add; cnt[p][d]--;
  }
}
int main(){
  for(int x=0;x<P;x++)chi[x]=0; for(int x=1;x<P;x++)chi[(x*x)%P]=1; for(int x=1;x<P;x++) if(!chi[x])chi[x]=-1;
  int k=0; for(int i=0;i<P;i++)for(int j=i+1;j<P;j++){pidx[i][j]=pidx[j][i]=k++;}
  ifstream f("templates.txt"); string tag;
  while(f>>tag){ Tmpl t; int nc; f>>t.name>>t.n>>t.m>>nc; t.E.resize(t.m); for(auto&e:t.E)f>>e.first>>e.second; for(int i=0;i<nc;i++){string S,D;f>>S>>D;t.cls.push_back({S,D});} T.push_back(t);}
  for(ti=0;ti<(int)T.size();ti++){
    Tmpl&t=T[ti]; cur_n=t.n;cur_m=t.m;curE=t.E; nbrE.assign(cur_n,{}); for(int e=0;e<cur_m;e++){nbrE[curE[e].first].push_back(e);nbrE[curE[e].second].push_back(e);}
    for(auto&c:t.cls){ curS=c.first;curD=c.second; memset(cnt,0,sizeof(cnt)); usedmask=0; dfs(0);} 
    fprintf(stderr,"%s: classes %zu instances(unique mask,usage) %zu\n",t.name.c_str(),t.cls.size(),insts[ti].size());
  }
  FILE*o=fopen("insts.bin","wb"); for(int ti2=0;ti2<(int)T.size();ti2++){ uint32_t n=insts[ti2].size(); fwrite(&n,4,1,o); fwrite(insts[ti2].data(),sizeof(Inst),n,o);} fclose(o);
}
