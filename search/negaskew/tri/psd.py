"""power spectrum of every block of a solution file on the cosets: python3 psd.py v "[H generators]" FILE MAXSHOW
(v = 1 + k(k+1)/2 prime, H the subgroup of order k of Z_v*; FILE: SOL lines of ../qskew or tri_join; the blocks are quaternary exponent strings, one digit per H-orbit)"""
import collections, numpy as np
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))      # ../qdrv.py
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from qdrv import shape
from sympy import primitive_root
v=int(sys.argv[1]); gens=eval(sys.argv[2]); fn=sys.argv[3]; maxshow=int(sys.argv[4])
d=shape(v,gens); n=d['n']; orb=d['orb']; e=n-1
g=primitive_root(v); cosrep=[pow(g,j,v) for j in range(e)]
I=np.array([1,1j,-1,-1j])
ns=collections.Counter(); shown=0
for l in open(fn):
    if not l.startswith('SOL'): continue
    E=[[int(c) for c in p] for p in l.split()[1:]]
    rows=[]
    for ex in E:
        u=I[np.array([ex[orb[x]] for x in range(v)])]
        F=np.fft.fft(u)   # F[t]=sum u(x) w^{-tx}
        P=np.abs(F)**2
        vals=np.round([P[t] for t in cosrep],4)
        rows.append(vals)
        ns[len(set(np.round(vals,2)))]+=1
    if shown<maxshow:
        shown+=1; print('SOL',l.split()[1:])
        for X,r in enumerate(rows): print('   ',X,'P(0)=%.1f'%(abs(sum(I[ex[orb[x]]] for x in range(v)))**2 if False else 0), r)
print('distinct-value counts per block', sorted(ns.items()))
