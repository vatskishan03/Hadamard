"""Free blocks of a solution related by a shift or by the tilde (and by the complex conjugation): the number of related pairs per solution (search/NOTES.md section 92: none beyond e = 6).
usage: python3 spinrel.py v k FILE        (FILE: SOL lines of ../qskew or of tri2)"""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from trilib import *
v=int(sys.argv[1]); k=int(sys.argv[2]); fn=sys.argv[3]
d,sols=load_sols(fn,v,k); e=d['e']; h=e//2
def shift(b,a): return (b[0],tuple(b[1][(j+a)%e] for j in range(e)))
def tilde(b): return (b[0] if b[0]==0 else (-b[0])%4, tuple((-b[1][(j+h)%e])%4 for j in range(e)))
def conjb(b): return ((-b[0])%4, tuple((-x)%4 for x in b[1]))
def orbit(b): 
    return set([shift(x,a) for x in (b,tilde(b)) for a in range(e)])
def orbit2(b):   # also conj and conj-tilde (PAF changes to its conjugate)
    return set([shift(x,a) for x in (b,tilde(b),conjb(b),tilde(conjb(b))) for a in range(e)])
cnt=collections.Counter(); cnt2=collections.Counter(); samesh=0
for s in sols:
    fr=s[1:]
    O=[orbit(b) for b in fr]
    key=tuple(sorted([(fr[i]in O[j]) for i in range(3) for j in range(i+1,3)]))
    cnt[sum((fr[i] in O[j]) for i in range(3) for j in range(i+1,3))]+=1
    O2=[orbit2(b) for b in fr]
    cnt2[sum((fr[i] in O2[j]) for i in range(3) for j in range(i+1,3))]+=1
print('v',v,'#sols',len(sols),'free blocks related by shift/tilde: #related pairs ->',dict(cnt),'; also allowing conj ->',dict(cnt2))
