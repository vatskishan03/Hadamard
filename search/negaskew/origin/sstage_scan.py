import sys, os, json
sys.path.insert(0,'.')
from sstage import *
from oo_shapes import shapes
import numpy as np
def parse(cb,q):
    """cb=(A,B,C,D) costs (B<=C<=D sorted); return list of (epsA,xA,epsB,epsD) variants for the (M,F,F,E) pattern, else None"""
    A=cb[0]; fr=list(cb[1:])
    if fr.count(0)!=1: return None
    f=[c for c in fr if c!=0]
    if any(c not in (q+1,q-1,q-3) for c in f): return None
    if A==(q+1)//2: eA=(1+1j)/2
    elif A==(q-3)//2: eA=-(1+1j)/2
    else: return None
    def e_of(c): return {q+1:[1],q-1:[1j,-1j],q-3:[-1]}[c]
    return [(eA,1,eb,ed) for eb in e_of(f[0]) for ed in e_of(f[1])]
if __name__=='__main__':
    NUMAX=int(sys.argv[1]); MMAX=int(sys.argv[2])
    sh=shapes(13,MMAX,NUMAX,5); sh.sort(key=lambda s:(s['nu'],s['v']))
    for s in sh:
        m,q,h,nu=s['m'],s['q'],s['hsize'],s['nu']
        from sympy import isprime
        if not isprime(m): continue
        tot=0; nv=0
        for cb in s['combos']:
            vs=parse(cb,q)
            if not vs: continue
            for (eA,xA,eb,ed) in vs:
                r=run(m,json.loads(json.dumps(s['gens'])),eA,xA,eb,ed,verbose=False)
                if r: tot+=r[2]; nv+=1
        print('nu',nu,'m',m,'q',q,'h',h,'(M,F,F,E) combos',sum(1 for cb in s['combos'] if parse(cb,q)),'variants',nv,'(Lock+S)-solutions total',tot,flush=True)
