"""Number of flat small blocks (|F(m)|^2 = v + 1 at every coset) per solution and the number of distinct values of |F(m)|^2 per kind of block (search/NOTES.md section 92: a flat free block only at 29 and 67).
usage: python3 stats2.py v k FILE        (FILE: SOL lines of ../qskew or of tri2)"""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from trilib import *
v=int(sys.argv[1]); k=int(sys.argv[2]); fn=sys.argv[3]
d,sols=load_sols(fn,v,k); e=d['e']
big=lambda b: abs(sum(I4[c] for c in b[1]))>1e-9
nflat=collections.Counter(); flatkinds=collections.Counter(); distinct=collections.Counter()
for s in sols:
    sp=[spectrum(d,b) for b in s]
    fl=[]
    for X in range(4):
        R,F=sp[X]; P=np.abs(F)**2
        if not big(s[X]) and np.abs(P-(v+1)).max()<1e-6: fl.append(X)
    nflat[len(fl)]+=1
    for X in range(4):
        R,F=sp[X]; P=np.abs(F)**2
        distinct[(('big' if big(s[X]) else ('skew' if X==0 else 'small')),len(set(np.round(P,4))))]+=1
print('v',v,'num flat small blocks per solution',dict(nflat))
print('(kind, #distinct |F|^2 values over cosets):',sorted(distinct.items()))
