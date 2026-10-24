"""Closure of a solution set under the shifts, the tilde of each free block and the swaps of blocks; classes under (shift, tilde, permutation of the free blocks) and the pools of distinct blocks used.  (symm2.py adds the sign change of A and is the one quoted in search/NOTES.md section 92.)
usage: python3 symm.py v k FILE        (FILE: SOL lines of ../qskew or of tri2)"""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from trilib import *
v=int(sys.argv[1]); k=int(sys.argv[2]); fn=sys.argv[3]
d,sols=load_sols(fn,v,k); e=d['e']; h=e//2
S=set(sols)
def shift(b,a): return (b[0],tuple(b[1][(j+a)%e] for j in range(e)))
def tilde(b): return (b[0] if b[0]==0 else (-b[0])%4, tuple((-b[1][(j+h)%e])%4 for j in range(e)))
big=lambda b: abs(sum(I4[c] for c in b[1]))>1e-9
# closure tests
def closed(f):
    return all(f(s) in S for s in sols)
print('v',v,'#sols',len(sols))
print('closed under shifts:',closed(lambda s: tuple(shift(b,1) for b in s)))
for X in (1,2,3):
    print('closed under tilde of block',X,closed(lambda s,X=X: tuple(tilde(b) if i==X else b for i,b in enumerate(s))))
print('closed under swap B,C:',closed(lambda s:(s[0],s[2],s[1],s[3])), 'swap C,D', closed(lambda s:(s[0],s[1],s[3],s[2])))
# orbit count by canonical form
def canon(s):
    best=None
    for a in range(e):
        sa=[shift(b,a) for b in s]
        # tilde choices: canonical minimal over tilde of each free block, then sort free blocks
        fr=[min(b,tilde(b)) for b in sa[1:]]
        key=(sa[0],tuple(sorted(fr)))
        if best is None or key<best: best=key
    return best
cl=collections.Counter(canon(s) for s in sols)
print('equivalence classes (shift, tilde, permutation):',len(cl),' orbit-size distribution:',sorted(collections.Counter(cl.values()).items()))
# pools
poolA=set(s[0] for s in sols); poolF=set(b for s in sols for b in s[1:])
print('distinct A used:',len(poolA),' distinct free blocks used (up to nothing):',len(poolF))
poolF2=set(min(b,tilde(b)) for b in poolF)
print('free blocks up to tilde:',len(poolF2))
# canonical up to shift+tilde of single free blocks
def canon_block(b):
    return min(min(shift(b,a),shift(tilde(b),a)) for a in range(e))
print('free blocks up to shift+tilde:',len(set(canon_block(b) for b in poolF)),' A up to shift:',len(set(min(shift(a,t) for t in range(e)) for a in poolA)))
