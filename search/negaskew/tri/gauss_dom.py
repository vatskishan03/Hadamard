"""Sparsity of the Z_e DFT (c-hat, the Gauss-sum form of PROOFS.md P6) of the coset sequences of the solutions: support sizes per kind of block (A, big, small) against random small blocks; spectral twins (two free blocks with equal |c-hat|); whether the sum of two PAFs is real (HEURISTIC statistics of EXACT solution sets; search/NOTES.md section 92: no sparsity).
usage: python3 gauss_dom.py v k FILE [NSAMP]        (FILE: SOL lines of ../qskew or of tri2; NSAMP = number of solutions sampled, default 2000)"""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from trilib import *
import random
v=int(sys.argv[1]); k=int(sys.argv[2]); fn=sys.argv[3]; nsamp=int(sys.argv[4]) if len(sys.argv)>4 else 2000
d,sols=load_sols(fn,v,k); e=d['e']; h=e//2
om=np.exp(2j*np.pi/e)
def chat(b):
    c=I4[np.array(b[1])]; return np.array([c@om**(-l*np.arange(e)) for l in range(e)])/e
def supp(b): return int((np.abs(chat(b))>1e-9).sum())
random.seed(1)
pool=random.sample(sols,min(nsamp,len(sols)))
sup=collections.defaultdict(collections.Counter)
for s in pool:
    for X in range(4):
        kind='A' if X==0 else ('big' if abs(sum(I4[c] for c in s[X][1]))>1e-9 else 'small')
        sup[kind][supp(s[X])]+=1
print('v',v,'support size of chat (Z_e DFT of coset sequence) in solutions:')
for kd in sup: print('  ',kd,sorted(sup[kd].items()))
# random baselines
rng=np.random.default_rng(3)
def rnd_small():
    a=rng.integers(0,h+1) if False else None
    while True:
        c=rng.integers(0,4,size=e)
        if abs(I4[c].sum())<1e-9: return (0,tuple(c))
base=collections.Counter(supp(rnd_small()) for _ in range(2000))
print('   random small free blocks',sorted(base.items()))
# spectral twins: |chat| equal between two free blocks?
tw=0; tw2=0
for s in pool:
    M=[np.abs(chat(s[X])) for X in range(1,4)]
    for i in range(3):
        for j in range(i+1,3):
            if np.allclose(M[i],M[j],atol=1e-9): tw+=1
            if np.allclose(np.sort(M[i]),np.sort(M[j]),atol=1e-9): tw2+=1
print('pairs of free blocks with identical |chat| vectors:',tw,' identical as multisets:',tw2,' of',3*len(pool),'pairs')
# imaginary part of P_A+P_B, P_B+P_C etc
def Pc(b):
    P=paf(d,b); return P
im=collections.Counter()
for s in pool[:500]:
    Ps=[Pc(b) for b in s]
    for (i,j) in [(0,1),(0,2),(0,3),(1,2),(1,3),(2,3)]:
        z=Ps[i]+Ps[j]
        im[(i,j,'real' if np.abs(z.imag).max()<1e-6 else 'cplx')]+=1
print('P_X+P_Y real?',sorted(im.items()))
