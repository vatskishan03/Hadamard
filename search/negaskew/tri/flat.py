"""power spectra of the blocks of a solution file restricted to coset representatives: number of distinct values per block: python3 flat.py v "[H generators]" FILE
(v = 1 + k(k+1)/2 prime, H the subgroup of order k of Z_v*; FILE: SOL lines of ../qskew or tri_join; the blocks are quaternary exponent strings, one digit per H-orbit)"""
import collections, numpy as np
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))      # ../qdrv.py
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from qdrv import shape
from sympy import primitive_root
v=int(sys.argv[1]); gens=eval(sys.argv[2]); fn=sys.argv[3]
d=shape(v,gens); n=d['n']; orb=d['orb']; e=n-1
g=primitive_root(v); cosrep=[pow(g,j,v) for j in range(e)]
I=np.array([1,1j,-1,-1j])
pat=collections.Counter(); ex1={}
for l in open(fn):
    if not l.startswith('SOL'): continue
    E=[[int(c) for c in p] for p in l.split()[1:]]
    key=[]
    for X,ex in enumerate(E):
        u=I[np.array([ex[orb[x]] for x in range(v)])]
        P=np.abs(np.fft.fft(u))**2
        vals=sorted(set(np.round([P[t] for t in cosrep],3)))
        key.append((len(vals),tuple(vals) if len(vals)<=2 else None))
    k=tuple(key); pat[k]+=1; ex1.setdefault(k,l.split()[1:])
for k,c in pat.most_common(12): print(c,k,ex1[k])
