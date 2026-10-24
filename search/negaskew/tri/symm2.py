"""The symmetry group Z_e x tilde^3 x S_3 x (A -> -A) of order 96 e (PROOFS.md P5): closure of a solution set under the sign change of A, PAF(A) = PAF(-A) on 200 random solutions, and the number of classes with the distribution of the class sizes (search/NOTES.md section 92: 24, 13, 105 classes at v = 29, 67, 191).
usage: python3 symm2.py v k FILE        (FILE: SOL lines of ../qskew or of tri2)"""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from trilib import *
v=int(sys.argv[1]); k=int(sys.argv[2]); fn=sys.argv[3]
d,sols=load_sols(fn,v,k); e=d['e']; h=e//2
S=set(sols)
def shift(b,a): return (b[0],tuple(b[1][(j+a)%e] for j in range(e)))
def tilde(b): return (b[0], tuple((-b[1][(j+h)%e])%4 for j in range(e)))
def negA(b):  # c -> -c ; u0: 0 <-> 3
    return ((3 if b[0]==0 else 0), tuple((c+2)%4 for c in b[1]))
print('v',v,'#sols',len(sols),'closed under negA:',all((negA(s[0]),)+s[1:] in S for s in sols))
# PAF equality of A and negA(A)
import random
ok=all(np.allclose(paf(d,s[0]),paf(d,negA(s[0]))) for s in random.sample(sols,min(200,len(sols))))
print('P_A = P_negA(A) on 200 random solutions:',ok)
def canon(s):
    best=None
    for a in range(e):
        for na in (0,1):
            A=shift(s[0],a); A=negA(A) if na else A
            fr=sorted(min(b,tilde(b)) for b in [shift(x,a) for x in s[1:]])
            key=(A,tuple(fr))
            if best is None or key<best: best=key
    return best
cl=collections.Counter(canon(s) for s in sols)
print('equivalence classes under Z_e x tilde^3 x S_3 x negA:',len(cl),' class size distribution',sorted(collections.Counter(cl.values()).items()))
