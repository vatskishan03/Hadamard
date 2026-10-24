"""Falsify the 'uniform cyclotomic' conjecture: take the class-value assignments (order-4 cyclotomic classes, p=2) of the 13,056 solutions at (13,23)
 (or 1,536 at (13,19)) and evaluate the q-free K-level equations (S),(T),(Q) at other primes m' = 5 mod 8 with the SAME class values.
 python3 unif_test.py FILE m' ...  (m=13, H_K=<3>)"""
import sys
sys.path.insert(0,'.')
from common import *
fn=sys.argv[1]; mps=[int(x) for x in sys.argv[2:]]
S13=Shape(13,[3]); orb13=S13.orb
def primroot(p):
    for g in range(2,p):
        if len({pow(g,k,p) for k in range(p-1)})==p-1: return g
def ind_table(p):
    g=primroot(p); tb={}
    x=1
    for k in range(p-1): tb[x]=k; x=x*g%p
    return tb
t13=ind_table(13)
cls13={}  # orbit index -> class mod 4
for k in range(1,13): cls13[orb13[k]]=t13[k]%4
print('orbit->class at 13',cls13)
sols=load_solutions(fn)
print(len(sols),'solutions')
def blocks_at(sol,mp,tb,rel):
    res=[]
    for b in sol:
        # class value arrays indexed by class j (0..3): find orbit with class j
        o_of={cl:o for o,cl in cls13.items()}
        def f(arr,k):
            if k==0: return U[arr[0]]
            return U[arr[o_of[(rel*tb[k])%4]]]
        x0=np.array([f(b[0],k) for k in range(mp)]); x1=np.array([f(b[1],k) for k in range(mp)]); x2=np.array([f(b[2],k) for k in range(mp)])
        e=(x1+x2)/2; s=x0-e; t=(x1-x2)/2
        res.append((s,t))
    return res
def pafc(x): 
    F=np.fft.fft(x); return np.fft.ifft(np.abs(F)**2)  # P(k)=sum conj(x(j)) x(j+k): conj convention irrelevant for zero test
def ccc(a,b):
    return np.fft.ifft(np.conj(np.fft.fft(a))*np.fft.fft(b))
for mp in mps:
    tb=ind_table(mp)
    best=[1e9,1e9,1e9]; nzero=0
    for rel in (1,3):
        for sol in sols:
            bl=blocks_at(sol,mp,tb,rel)
            Ps=sum(pafc(s) for s,t in bl); Pt=sum(pafc(t) for s,t in bl)
            Q=sum(ccc(t,s)-ccc(s,t) for s,t in bl)
            rs=np.max(np.abs(Ps[1:])); rt=np.max(np.abs(Pt[1:])); rq=np.max(np.abs(Q))
            best[0]=min(best[0],rs); best[1]=min(best[1],rt); best[2]=min(best[2],rq)
            if rs<1e-9 and rt<1e-9 and rq<1e-9: nzero+=1
    print('m',mp,'min residual max|sum P_s|,|sum P_t|,|Q| over tested solutions:',[round(x,3) for x in best],'exact solutions:',nzero,flush=True)
