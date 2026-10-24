"""elementary row-sum proof of the q = 3, nu = 11 EMPTY verdicts of data/RESULTS.txt (and of the other shapes whose level 1 is empty by row sums alone): in the real class every H-orbit of G pairs with its negative, so A has row sum 1,
a free block has row sum the signed sum of its orbit sizes (l, l(q-1)/2, l(q-1)/2), and the four squares must add to 4v.  usage: python3 dc_check.py   (no engine involved; prints the v whose verdict this explains)"""
import re, json, sys, os, subprocess, itertools
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, '..'))
from qshape import kshape
rows=[]
for l in open(os.path.join(HERE, 'data', 'RESULTS.txt')):
    m=re.match(r'(EMPTY|TIMEOUT)\s+v=(\d+)\s+(\S+)\s+q=(\d+)\s+sec=\S+\s+l1=(\S+)',l)
    if m: rows.append(m.groups())
def reach(sizes_with_mult):
    # set of achievable signed sums of the given multiset of orbit sizes
    S={0}
    for s in sizes_with_mult:
        S={x+s for x in S}|{x-s for x in S}
    return S
out=[]
for status,v,name,q,l1 in rows:
    mm=re.match(r'[cn]\d+_K(\d+)_q(\d+)_h([\d_]+)$',name)
    if not mm: continue
    m=int(mm.group(1)); q=int(mm.group(2)); gens=[int(x) for x in mm.group(3).split('_') if x]
    v=int(v)
    if m*q!=v or q not in (3,7,11,19,23,31,43,47,59,71): continue
    try:
        sh=kshape([m],gens)
    except Exception as e:
        continue
    ell=sh['ell']
    sizes=[]
    for o,e in enumerate(ell):
        sizes += [e, e*(q-1)//2, e*(q-1)//2] if True else []
    # origin orbit: (0,0) size 1, (0,f in sq) size (q-1)/2 ...: K-orbit 0 has ell=1 so same formula
    S=reach(sizes)
    sq={x for x in S}
    S2=sorted({x*x for x in S})
    Sset=set(S2)
    target=4*v-1
    pairs={a+b for a in S2 for b in S2}
    ok=any((target-a) in pairs for a in S2)
    out.append((status,v,name,l1,ok))
for o in out: print(o)
print('DC-infeasible (explained by row sums alone):',[o[1] for o in out if not o[4]])
