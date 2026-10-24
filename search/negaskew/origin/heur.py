"""HEURISTIC expected number of origin-defect solutions per shape (calibrated on 299, 589, 247):
 E = [ L * 16^{#E} * 4^7 / ((5.2 m) * (12 m)^2) ]^p   (z=0: L = Lock states per pair; z=1 analogously), only for shapes whose row-sum test passes.
 usage: python3 heur.py MMAX NUMAX"""
import sys, itertools, math
sys.path.insert(0,'.')
from oo_shapes import shapes, combos
from lockcount import states
from rowfeas import feasible
from sympy import isprime
Ae=[(1+1j)/2,-(1+1j)/2]
Me=[(1+1j)/2,(1-1j)/2,(-1+1j)/2,(-1-1j)/2]
def eps_choices(cost,q):
    if cost==0: return None   # E
    if cost in ((q+1)//2,(q-3)//2): return Me
    s=cost-(q-1)   # q+1 -> 2 : eps=+1 ; q-1 -> 0 : eps=+-i ; q-3 -> -2 : eps=-1
    return {2:[1],0:[1j,-1j],-2:[-1]}[s]
from functools import lru_cache
@lru_cache(maxsize=None)
def Lcount_sig(sig, z):
    chs=[]; kinds=[]
    for i,c in enumerate(sig):
        if i==0: chs.append(Ae); kinds.append('A')
        elif c is None: continue
        else: chs.append(tuple(c)); kinds.append('F')
    best=0
    for eps in itertools.product(*chs): best=max(best,len(states(list(eps),kinds,z)))
    return best
def Lcount(cost, q, z):
    sig=tuple((None if (i>0 and c==0) else (tuple(eps_choices(c,q)) if i>0 else None)) for i,c in enumerate(cost))
    return Lcount_sig(sig, z)
def Lcount_old(cost, q, z):
    chs=[]; kinds=[]
    for i,c in enumerate(cost):
        if i==0:
            chs.append(Ae); kinds.append('A')   # A's e = +-(1+i)/2 ; which sign is consistent with cost: both x0 options exist
        else:
            ec=eps_choices(c,q)
            if ec is None: continue
            chs.append(ec); kinds.append('F')
    best=0
    # for A restrict e sign by cost via x0 in {1,-i}: cost (q+1)/2 or (q-3)/2 both reachable by some (x0,e); take all
    for eps in itertools.product(*chs):
        best=max(best,len(states(list(eps),kinds,z)))
    return best
def calibrated(L,nE,m,p):
    return (L*16**nE*4**7/((5.2*m)*(12*m)**2))**p
if __name__=='__main__':
    MMAX=int(sys.argv[1]); NUMAX=int(sys.argv[2])
    sh=shapes(7,MMAX,NUMAX,3); sh.sort(key=lambda s:(s['v']))
    rows=[]
    for s in sh:
        m,q,h,nu=s['m'],s['q'],s['hsize'],s['nu']
        if not isprime(m): continue
        p=(nu-1)//2
        best=0; bc=None
        for cb in s['combos']:
            nE=sum(1 for c in cb[1:] if c==0)
            for z in (0,1):
                L=Lcount(cb,q,z)
                if L==0: continue
                if not all(feasible(m,q,h,p,cb,z)): continue
                E=calibrated(L,nE,m,p)
                if E>best: best=E; bc=(cb,z,L,nE)
        rows.append((s['v'],m,q,h,nu,best,bc))
    for r in rows:
        print('v %6d  m %3d q %4d h %2d nu %2d  E_pred %.3g  %s'%r)
