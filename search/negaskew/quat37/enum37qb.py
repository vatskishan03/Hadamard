"""non-cyclic shapes: K = product of cyclic groups (Kmods) with diagonal unit H_K (odd order), q = prime power = 3 mod 4 (p^e e odd, p=3 mod 4),
v = |K| q; for Table-5 targets (m>=4, v<=999) and census-undecided v (<=2499) not in repository.  Writes shapes37b.json"""
import json, os, sys, itertools
HERE = os.path.dirname(os.path.abspath(__file__)); REPO = os.path.join(HERE, '..', '..', '..')
from math import gcd
from sympy import isprime, factorint
sys.argv=sys.argv[:1]
sys.path.insert(0, HERE)
from qshape37q import kshape
have=set(int(x) for x in open(os.path.join(HERE, 'data', 'r36_have.txt')).read().split())
for d in os.listdir(os.path.join(REPO, 'skew')):
    if d.isdigit(): have.add(int(d))
T5={int(k):v for k,v in json.load(open(os.path.join(HERE, '..', 't5.json'))).items()}
cen=[int(l.split()[0]) for l in open(os.path.join(REPO, 'search', 'data', 'skew_census_undecided.txt')) if l.strip() and not l.startswith('#')]
targets={v for v,m in T5.items() if m>=4 and 4*v not in have and 8*v not in have and v<=999}
census={v for v in cen if 4*v not in have and v>999}
vs=sorted(targets|census)
def ppq(q):
    f=factorint(q)
    if len(f)!=1: return None
    p,e=next(iter(f.items()))
    return (p,e) if p%4==3 and e%2==1 else None
def Kcands(m):
    # factorisations of K into cyclic groups: Z_m cyclic handled elsewhere (only when q is a prime power with e>=3); non-cyclic: Z_a x Z_b with a|b, both odd >=3
    out=[]
    for a in range(3,int(m**0.5)+2,2):
        if m%a==0:
            b=m//a
            if b%a==0 and b>=a: out.append([a,b])
            # also non-cyclic a x b with gcd(a,b)>1 not dividing is cyclic-iso if coprime
    return out
rows=[]
for v in vs:
    for q in range(3,v):
        if v%q or not ppq(q): continue
        p,e=ppq(q)
        m=v//q
        if m%2==0 or m<3 or gcd(m,q)!=1: continue
        if q>8*m+3: continue
        Ks=[]
        if e>=3: Ks.append([m])
        Ks+=Kcands(m)
        for Km in Ks:
            # diagonal unit tuples of odd order
            L=1
            for a in Km: L=L*a//gcd(L,a)
            units=[u for u in itertools.product(*[[x for x in range(1,a) if gcd(x,a)==1] for a in Km])]
            def order(u):
                x=u; k=1
                while any(xx%a!=1 for xx,a in zip(x,Km)): x=tuple(xx*uu%a for xx,uu,a in zip(x,u,Km)); k+=1
                return k
            odd=[u for u in units if order(u)%2==1]
            subs={}
            def clos(gens):
                H={tuple(1 for _ in Km)}; fr=list(H)
                while fr:
                    x=fr.pop()
                    for g in gens:
                        y=tuple(xx*gg%a for xx,gg,a in zip(x,g,Km))
                        if y not in H: H.add(y); fr.append(y)
                return frozenset(H)
            for u in odd: subs.setdefault(clos([u]),[u])
            items=list(subs.items())
            for (H1,g1),(H2,g2) in itertools.combinations(items,2):
                H=clos(g1+g2)
                if H not in subs: subs[H]=g1+g2
            subs.setdefault(clos([]),[])
            for H,g in subs.items():
                try: sh=kshape(Km,[tuple(x) for x in g])
                except ValueError: continue
                if sh['nu']>9: continue
                rows.append(dict(v=v,Km=Km,q=q,p=p,e=e,gens=[list(x) for x in g],Hsize=len(H),nu=sh['nu'],ells=sh['ell'],m5=T5.get(v)))
json.dump(rows,open('shapes37b.json','w'))
print(len(rows))
for r in sorted(rows,key=lambda r:(r['nu'],r['v'])): print(r['nu'],r['v'],r['Km'],r['q'],'e',r['e'],r['gens'],'|H|',r['Hsize'],r['ells'],'T5',r['m5'])
