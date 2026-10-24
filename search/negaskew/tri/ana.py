"""row-sum patterns (|R_X|^2 per block) of the solutions in a qskew SOL file: python3 ana.py v "[H generators]" FILE
(v = 1 + k(k+1)/2 prime, H the subgroup of order k of Z_v*; FILE: SOL lines of ../qskew or tri_join; the blocks are quaternary exponent strings, one digit per H-orbit)"""
import collections, numpy as np
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))      # ../qdrv.py
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from qdrv import shape
v=int(sys.argv[1]); gens=eval(sys.argv[2]); fn=sys.argv[3]
d=shape(v,gens); n=d['n']; sz=d['sz']
I=[1,1j,-1,-1j]
pat=collections.Counter(); sols=[]
for l in open(fn):
    if not l.startswith('SOL'): continue
    E=[[int(c) for c in p] for p in l.split()[1:]]
    R=[sum(I[e[o]]*sz[o] for o in range(n)) for e in E]
    big=tuple(int(abs(r)**2) for r in R)
    pat[big]+=1; sols.append((E,R))
print('v',v,'k',len(d['H']),'n',n,'sols',len(sols))
for k_,c in pat.most_common(20): print(k_,c)
