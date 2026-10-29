"""Floating-point (FFT) check of the SOL lines of a tri2 run at v = 379, k = 27 (or the v, k given; blocks in coset order): the four PAF vectors sum to 4v at 0 and to 0 elsewhere (tolerance 1e-6), and the row-sum pattern |R_X|^2 of each solution (search/NOTES.md section 92: (1513, 1, 1, 1)).  The integer verification is done inside tri2 itself.
usage: python3 check379.py FILE [NMAX [v k]]        (NMAX = number of solutions printed, default 5; v k default 379 27, e.g. 277 23 for the reduced 277 log)"""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from trilib import *
v,k=(int(sys.argv[3]),int(sys.argv[4])) if len(sys.argv)>4 else (379,27)      # optional v k (default 379 27)
fn=sys.argv[1]; nmax=int(sys.argv[2]) if len(sys.argv)>2 else 5
d=setup(v,k); e=d['e']; print('g',d['g'],'e',e)
okc=0; n=0; pats=collections.Counter()
for l in open(fn):
    if not l.startswith('SOL'): continue
    p=l.split(); blocks=[]
    for s in p[4:8]:
        x=[int(c) for c in s]; blocks.append((x[0],tuple(x[1:])))
    P=sum(paf(d,b) for b in blocks)
    ok=(np.abs(P[1:]).max()<1e-6 and abs(P[0]-4*v)<1e-6)
    R=[spectrum(d,b)[0] for b in blocks]; pats[tuple(int(round(abs(r)**2)) for r in R)]+=1
    n+=1; okc+=ok
    if n<=nmax: print('SOL',p[1:4],'verified (FFT, floating) =',ok, 'row sums |R|^2',[int(round(abs(r)**2)) for r in R])
print('checked',n,'all ok',okc==n,'; row-sum patterns',dict(pats))
