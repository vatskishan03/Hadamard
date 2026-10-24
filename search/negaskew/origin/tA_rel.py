import sys
sys.path.insert(0,'.')
from common import *
from collections import Counter
m=int(sys.argv[1]); q=int(sys.argv[2]); hK=eval(sys.argv[3]); fn=sys.argv[4]
S=Shape(m,hK); orb=S.orb; nu=S.nu
g=2 if m!=7 else 3
# class index via discrete log base g' (primitive root) mod number of classes e=(nu-1)
from sympy import primitive_root
gr=primitive_root(m); e=nu-1
def ind(k):
    x=1
    for j in range(m-1):
        if x==k%m: return j
        x=x*gr%m
cidx={}
for k in range(1,m): cidx[orb[k]]=ind(k)%e
inv={c:o for o,c in cidx.items()}
sols=load_solutions(fn)
cnt=Counter(); tot=0
for s in sols:
    a=[s[0][0][inv[c]] for c in range(e)]; v=[s[0][1][inv[c]] for c in range(e)]
    found=[]
    for conj in (0,1):
        for sh in range(e):
            for r in range(4):
                if all(((v[c]-r-((-1 if conj else 1)*a[(c+sh)%e]))%4)==0 for c in range(e)): found.append((conj,sh))
    cnt[tuple(sorted(set(found)))]+=1
    tot+=1
print(tot)
for k,v in cnt.most_common(): print(v,k)
