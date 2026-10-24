"""is 0 in the coset of the lattice spanned by the differences of the PAF vectors of sampled big (A) and small (B, C, D) blocks, i.e. is the PAF equation Sum_X PAF_X = 0 solvable over Z modulo the lattice?  For v, k = (67, 11), (191, 19), (379, 27), (631, 35), (947, 43): python3 coset.py   (needs sympy; ~1 min)"""
import numpy as np, sys
import tri_est as T
from sympy import Matrix
def member(basis, x):
    # basis rows in echelon form (pivot columns increasing); reduce x
    x=list(map(int,x))
    for b in basis:
        p=next(i for i,t in enumerate(b) if t)
        if x[p]%b[p]: return False
        q=x[p]//b[p]; x=[xi-q*bi for xi,bi in zip(x,b)]
    return not any(x)
for v,k in [(67,11),(191,19),(379,27),(631,35),(947,43)]:
    d=T.setup(v,k)
    PA=np.array([T.pafvec(d,x) for x in T.sample_big(d,400)]); PS=np.array([T.pafvec(d,x) for x in T.sample_small(d,400)])
    diffs=[x-PA[0] for x in PA[1:]]+[x-PS[0] for x in PS[1:]]
    B=T.int_basis(diffs)
    base=PA[0]+3*PS[0]
    print(v,'rank',len(B),'0 in coset:',member(B,-base))
