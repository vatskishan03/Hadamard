"""Kinds of blocks of the solutions: positions of the big block, number of real cosets (digits 0 or 2) per kind of block, and whether the small free blocks are all real.
usage: python3 stats1.py v k FILE        (FILE: SOL lines of ../qskew or of tri2)"""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from trilib import *
v=int(sys.argv[1]); k=int(sys.argv[2]); fn=sys.argv[3]
d,sols=load_sols(fn,v,k); e=d['e']
big=lambda b: abs(sum(I4[c] for c in b[1]))>1e-9
cntR=collections.Counter(); allreal=collections.Counter(); bigpos=collections.Counter()
typ=collections.Counter()
for s in sols:
    bp=[X for X in range(4) if big(s[X])]; bigpos[tuple(bp)]+=1
    for X in range(4):
        b=s[X]; nr=sum(1 for c in b[1] if c%2==0)
        cntR[(('big' if X in bp else ('skew' if X==0 else 'small')),nr)]+=1
    sm=[X for X in range(4) if X not in bp and X!=0]
    allreal[tuple(sorted(sum(1 for c in s[X][1] if c%2==1)==0 for X in sm))]+=1
print('v',v,'sols',len(sols),'big positions',dict(bigpos))
print('(kind, #real cosets) counts',sorted(cntR.items()))
print('free small blocks all-real? (sorted tuple over the 2 or 3 small free blocks)',dict(allreal))
