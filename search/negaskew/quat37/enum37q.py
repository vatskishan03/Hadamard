"""enumerate quaternary plug-in shapes (K=Z_m cyclic, q prime = 3 mod 4, gcd=1, H_K odd-order subgroup of Z_m*) for a list of v.
usage: python3 enum37q.py OUT.json VLO VHI NUMAX [census|targets|all]   (census: the census-undecided v in (VLO, VHI] with 4v not in skew/; targets: Table 5 entries m >= 4; origin-defect combinations are listed per shape as `oof`)"""
import json, os, sys, itertools
HERE = os.path.dirname(os.path.abspath(__file__)); REPO = os.path.join(HERE, '..', '..', '..')
from math import gcd
from sympy import isprime
OUT=sys.argv[1]; VLO=int(sys.argv[2]); VHI=int(sys.argv[3]); NUMAX=int(sys.argv[4]); MODE=sys.argv[5] if len(sys.argv)>5 else 'census'
def odd_subgroups(m):
    U=[u for u in range(1,m) if gcd(u,m)==1]
    def order(u):
        x,k=u,1
        while x!=1: x=x*u%m; k+=1
        return k
    def closure(gens):
        H={1}; fr=[1]
        while fr:
            x=fr.pop()
            for g in gens:
                y=x*g%m
                if y not in H: H.add(y); fr.append(y)
        return frozenset(H)
    subs={frozenset([1]):[1]}
    for u in U:
        if order(u)%2==1 and u!=1:
            subs.setdefault(closure([u]),[u])
    changed=True
    while changed:
        changed=False
        items=list(subs.items())
        for (H1,g1),(H2,g2) in itertools.combinations(items,2):
            if H1<=H2 or H2<=H1: continue
            g=g1+g2; H=closure(g)
            if H not in subs: subs[H]=g; changed=True
    return subs
def orbits(m,H):
    seen=[False]*m; orbs=[]
    for k in range(m):
        if not seen[k]:
            o=sorted({k*h%m for h in H})
            for x in o: seen[x]=True
            orbs.append(o)
    return orbs
def weight_sets(q,ells,negpairs):
    cs=sorted({q+1,q-3,q-1,(q+1)//2,(q-3)//2})
    reach={0}
    for l in ells: reach={r+x for r in reach for x in [0]+[l*c for c in cs if c>0]}
    WB=reach
    reach={(q+1)//2,(q-3)//2}
    for l in negpairs: reach={r+x for r in reach for x in [0]+[2*l*c for c in cs if c>0]}
    return reach,WB
have=set(int(x) for x in open(os.path.join(HERE, 'data', 'r36_have.txt')).read().split())
for d in os.listdir(os.path.join(REPO, 'skew')):
    if d.isdigit(): have.add(int(d))
if MODE=='census':
    cen=[int(l.split()[0]) for l in open(os.path.join(REPO, 'search', 'data', 'skew_census_undecided.txt')) if l.strip() and not l.startswith('#')]
    vs=[v for v in cen if 4*v not in have and VLO<v<=VHI]
else:
    T5={int(k):v for k,v in json.load(open(os.path.join(HERE, '..', 't5.json'))).items()}
    vs=sorted(v for v,m in T5.items() if m>=4 and 4*v not in have and 8*v not in have and VLO<v<=VHI)
print(len(vs),'values',file=sys.stderr)
out=[]; cache={}
for v in vs:
    for q in range(3,v):
        if v%q or not isprime(q) or q%4!=3: continue
        m=v//q
        if m%2==0 or m<3 or gcd(m,q)!=1: continue
        if q>8*m+3: continue          # origin half defect (q-3)/2 <= 4m
        if m not in cache: cache[m]=odd_subgroups(m)
        for H,gens in cache[m].items():
            orbs=orbits(m,H); nu=len(orbs)
            if nu<3 or nu>NUMAX: continue
            ells=[len(o) for o in orbs]
            idx={}
            for i,o in enumerate(orbs):
                for x in o: idx[x]=i
            pairs=[]; seen=set()
            for i in range(1,nu):
                if i in seen: continue
                j=idx[(-orbs[i][0])%m]; seen|={i,j}; pairs.append(len(orbs[i]))
            WA,WB=weight_sets(q,ells,pairs)
            tot=4*m
            WBs=WB; pair2={a+b for a in WB for b in WB}
            ident=any((tot-a-s) in WB for a in WA for s in pair2)
            # origin-defect family
            A0s=[(q+1)//2,(q-3)//2]; Os=sorted({0,(q+1)//2,(q-3)//2,q+1,q-3,q-1})
            oof=[(a,b,c,d) for a in A0s for b in Os for c in Os if c>=b for d in Os if d>=c and a+b+c+d==tot]
            out.append(dict(v=v,K=m,q=q,gens=sorted(gens),Hsize=len(H),nu=nu,ells=ells,ident=ident,oof=oof))
json.dump(out,open(OUT,'w'))
from collections import Counter
print(len(out),'shapes; identity-feasible',sum(s['ident'] for s in out),'; origin-defect feasible',sum(1 for s in out if s['oof']),Counter(s['nu'] for s in out if s['ident']))
