"""is there a skew block A with PAF = 3 off the origin (a Legendre-type, flat A)?  Exhaustive over the 4^pairs skew blocks of the shape: python3 leg.py v "[H generators]"
(v = 1 + k(k+1)/2 prime, H the subgroup of order k of Z_v*)"""
import numpy as np, itertools
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))      # ../qdrv.py
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from qdrv import shape
v=int(sys.argv[1]); gens=eval(sys.argv[2])
d=shape(v,gens); n=d['n']; orb=d['orb']; neg=d['neg']
I=np.array([1,1j,-1,-1j])
chi=np.array([0]+[1 if pow(x,(v-1)//2,v)==1 else -1 for x in range(1,v)]); L=chi.copy(); L[0]=1
P=lambda u: np.array([np.vdot(np.roll(u,-s),u) for s in range(v)])  # sum conj(u(x+s)) u(x)
print('v',v,'Legendre PAF off 0:',set(np.round(P(L.astype(complex))[1:].real,6)))
pairs=[a for a in range(1,n) if a<neg[a]]
best=None; hits=0
for u0 in (0,3):
    for ch in itertools.product(range(4),repeat=len(pairs)):
        ex=[None]*n; ex[0]=u0
        for a,c in zip(pairs,ch): ex[a]=c; ex[neg[a]]=(1-c)%4
        u=I[np.array([ex[orb[x]] for x in range(v)])]
        p=P(u)[1:]
        dev=np.abs(p-3).max()
        if best is None or dev<best[0]: best=(dev,ex)
        if dev<1e-9: hits+=1
print('A with PAF=3 off 0:',hits,'best max deviation',round(best[0],3))
