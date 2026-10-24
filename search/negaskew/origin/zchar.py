import sys
sys.path.insert(0,'.')
from common import *
from collections import Counter
m=int(sys.argv[1]); hK=eval(sys.argv[2]); fn=sys.argv[3]
S=Shape(m,hK); orb=S.orb
sols=load_solutions(fn)
def cc(a,b): return np.array([np.sum(a*np.conj(np.roll(b,-k))) for k in range(m)])
legendre=lambda k: 0 if k%m==0 else (1 if pow(k,(m-1)//2,m)==1 else -1)
# orbit classes
cls=[None]*S.nu
for k in range(1,m): cls[orb[k]]=legendre(k)
print('orbit quad class',cls)
cnt=Counter(); aff=Counter(); tot=Counter()
for sol in sols:
    Z1=np.zeros(m,complex)
    xs=[]
    for b in sol:
        x0=np.array([U[b[0][orb[k]]] for k in range(m)]); x1=np.array([U[b[1][orb[k]]] for k in range(m)]); x2=np.array([U[b[2][orb[k]]] for k in range(m)])
        e=(x1+x2)/2; s=x0-e
        Z1+=cc(e,s)+cc(s,e); xs.append((b,e,s))
    z=-round(Z1[1].real)
    for X,(b,e,s) in enumerate(xs):
        # x0 pattern on non-origin orbits
        d={}
        for o in range(S.nu):
            if cls[o] is None: continue
            d.setdefault(cls[o],set()).add(b[0][o])
        isaff=all(len(v)==1 for v in d.values())
        # also t pattern (x1 values) 
        d1={}
        for o in range(S.nu):
            if cls[o] is None: continue
            d1.setdefault(cls[o],set()).add((b[1][o],b[2][o]))
        isaff1=all(len(v)==1 for v in d1.values())
        tot[(z,X)]+=1; aff[(z,X)]+= isaff
        aff[(z,X,'x1x2')]+= isaff1
for z in (0,1):
    print('z',z,[ (X,aff[(z,X)],aff[(z,X,'x1x2')],tot[(z,X)]) for X in range(4)])
