"""Defect-identity feasibility for Paley-type shapes of the quaternary skew problem, G = Z_K x F_q, H = {1} x QR(q) x H_K.
EXACT necessary condition: sum over blocks X and k in K of [ (q-1)|e|^2 + 2 Re(x0 conj e) ] = 4|K|, weights per defect K-position:
 full (x1=x2=z): q+1 (x0=z), q-3 (x0=-z), q-1 (x0=+-iz);  mixed (x2=+-i x1): (q+1)/2 or (q-3)/2;  none: 0.
 A (twisted skew) has a forced mixed position at k=0 and its other defects come in (k,-k) pairs of equal weight.
usage: python3 idfeas.py K q [h ...]    one shape: K, q, generators of H_K in Z_K*   (e.g. python3 idfeas.py 5 59)
       python3 idfeas.py                every shape of runs/paley_shapes.json (written by paley_shapes.py)
Prints IDENTITY-FEASIBLE or infeasible (EXACTLY empty) for each."""
import sys, json
from math import gcd
def korbits(K, HK):
    seen=[-1]*K; orbs=[]
    for k in range(K):
        if seen[k]<0:
            o=sorted(set((k*h)%K for h in HK)); 
            for x in o: seen[x]=len(orbs)
            orbs.append(o)
    return orbs, seen
def feas(K,q,HK,verbose=False):
    orbs,idx=korbits(K,HK); kap=len(orbs)
    sz=[len(o) for o in orbs]
    neg=[idx[(-o[0])%K] for o in orbs]
    Wfull=[q+1,q-3,q-1]; Wmix=[(q+1)//2,(q-3)//2]
    total=4*K
    # items: (block, orbit) slots; A: orbit 0 forced mixed; A other orbits in pairs (i,neg i) i<neg i  (weight counted twice)
    # DP over multiset of (contribution values): reachable set of totals
    reach={0}
    # A@0
    reach={w for w in Wmix}
    def add_item(reach,vals):
        out=set()
        for r in reach:
            out.add(r)
            for c in vals:
                if r+c<=total: out.add(r+c)
        return out
    # A pairs
    for i in range(1,kap):
        if i<neg[i]:
            vals=[2*sz[i]*w for w in Wfull+Wmix]
            reach=add_item(reach,vals)
    # B,C,D: all orbits
    for blk in range(3):
        for i in range(kap):
            vals=[sz[i]*w for w in Wfull+Wmix]
            reach=add_item(reach,vals)
    return total in reach, kap, sz
if __name__=="__main__":
    from enum_shapes import gen_sub
    if len(sys.argv)>=3:
        K=int(sys.argv[1]); q=int(sys.argv[2]); hk=[int(x) for x in sys.argv[3:]]
        HK=sorted(gen_sub(K,hk)) if hk else [1]
        ok,kk,sz=feas(K,q,HK)
        print('K',K,'q',q,'|HK|',len(HK),'kappa',kk,'sizes',sorted(set(sz)),'IDENTITY-FEASIBLE' if ok else 'infeasible')
        sys.exit(0)
    r=json.load(open('runs/paley_shapes.json'))
    for v,m,p,K,n,o,g,kap in r:
        if K==1: continue
        H=gen_sub(v,g)
        HK=sorted(set(h%K for h in H))
        ok,kk,sz=feas(K,p,HK)
        print(v,m,'K',K,'q',p,'|HK|',len(HK),'kappa',kk,'sizes',sorted(set(sz)),'IDENTITY-FEASIBLE' if ok else 'infeasible', flush=True)
