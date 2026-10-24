"""relations between the small blocks of a solution file (equal up to unit and the tilde B(g) -> conj B(-g), self-dual blocks): python3 rel.py v "[H generators]" FILE
(v = 1 + k(k+1)/2 prime, H the subgroup of order k of Z_v*; FILE: SOL lines of ../qskew or tri_join; the blocks are quaternary exponent strings, one digit per H-orbit)"""
import collections, numpy as np
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))      # ../qdrv.py
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from qdrv import shape
v=int(sys.argv[1]); gens=eval(sys.argv[2]); fn=sys.argv[3]
d=shape(v,gens); n=d['n']; neg=d['neg']; sz=d['sz']
I=[1,1j,-1,-1j]
def canon(ex):
    # block-wise symmetries: multiply by i^a, and tilde: B(g) -> conj B(-g)
    out=[]
    for t in (0,1):
        b=list(ex) if t==0 else [(-ex[neg[o]])%4 for o in range(n)]
        for a in range(4): out.append(tuple((x+a)%4 for x in b))
    return min(out)
def selfdual(ex):
    td=[(-ex[neg[o]])%4 for o in range(n)]
    return any(tuple((x+a)%4 for x in td)==tuple(ex) for a in range(4))
cnt=collections.Counter(); tot=0
for l in open(fn):
    if not l.startswith('SOL'): continue
    E=[[int(c) for c in p] for p in l.split()[1:]]; tot+=1
    R=[sum(I[e[o]]*sz[o] for o in range(n)) for e in E]
    small=[X for X in range(4) if abs(abs(R[X])**2-1)<1e-9]
    cs=[canon(E[X]) for X in small]
    eqp=sum(1 for a in range(len(cs)) for b in range(a+1,len(cs)) if cs[a]==cs[b])
    sd=sum(selfdual(E[X]) for X in small)
    cnt[('small',tuple(small),'equalpairs',eqp,'selfdual_small',sd)]+=1
print('v',v,'total',tot)
for k_,c in sorted(cnt.items(), key=lambda x:-x[1])[:15]: print(c,k_)
