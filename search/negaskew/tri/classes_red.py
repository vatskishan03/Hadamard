"""Equivalence classes of the reduced hits of a tri2 run with red=1 under the group Z_e x tilde^3 x S_3 x (A -> -A) of PROOFS.md P5 (canonical form: the minimum over the e shifts and the sign change of A, the free blocks sorted by min(block, tilde)); prints the number of classes and the count implied if every orbit has the regular size 96 e.
usage: python3 classes_red.py v k FILE        (FILE: a tri2 run with red=1, e.g. data/c277_red.out)"""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from trilib import *
v,k=int(sys.argv[1]),int(sys.argv[2]); fn=sys.argv[3]
d=setup(v,k); e=d['e']; h=e//2
def shift(b,a): return (b[0],tuple(b[1][(j+a)%e] for j in range(e)))
def tilde(b): return (b[0], tuple((-b[1][(j+h)%e])%4 for j in range(e)))
def negA(b): return ((3 if b[0]==0 else 0), tuple((c+2)%4 for c in b[1]))
def canon(s):
    best=None
    for a in range(e):
        for na in (0,1):
            A=shift(s[0],a); A=negA(A) if na else A
            fr=sorted(min(b,tilde(b)) for b in [shift(x,a) for x in s[1:]])
            key=(A,tuple(fr))
            if best is None or key<best: best=key
    return best
cl=collections.Counter(); cw=collections.Counter(); n=0; wsum=0
for l in open(fn):
    if not l.startswith('SOL'): continue
    p=l.split(); n+=1; wsum+=float(p[3][2:])
    bl=[]
    for s in p[4:8]:
        x=[int(c) for c in s]; bl.append((x[0],tuple(x[1:])))
    cl[canon(tuple(bl))]+=1; cw[canon(tuple(bl))]+=float(p[3][2:])
print('reduced hits',n,'weighted sum',wsum,'distinct classes under the full group (Z_e x tilde^3 x S_3 x negA):',len(cl),'; reduced hits per class',sorted(collections.Counter(cl.values()).items()))
print('implied normalised count if all orbits regular (order 2*e*8*6 =',2*e*48,'):',len(cl)*2*e*48,' vs 6*wsum =',6*wsum)
print('class weight sums (orbit size / 6 in mode 7, / 1 in mode 3; a regular class has 16 * 12 = 192 at e = 12):',sorted(collections.Counter(cw.values()).items()),' total',sum(cw.values()))   # new in the port: output only
