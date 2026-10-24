import sys
sys.path.insert(0,'.')
from common import *
from collections import Counter
m=int(sys.argv[1]); q=int(sys.argv[2]); hK=eval(sys.argv[3]); fn=sys.argv[4]
S=Shape(m,hK); orb=S.orb
sols=load_solutions(fn)
def typ(b):
    d=(b[2][0]-b[1][0])%4
    e=(U[b[1][0]]+U[b[2][0]])/2; x=U[b[0][0]]
    c=round((q-1)*abs(e)**2+2*(x*e.conjugate()).real)
    return {2:'E',0:'F',1:'M',3:'M'}[d]+str(c)
cnt=Counter()
step=max(1,len(sols)//4000)
for sol in sols[::step]:
    for X,b in enumerate(sol):
        x0=np.array([U[b[0][orb[k]]] for k in range(m)]); x1=np.array([U[b[1][orb[k]]] for k in range(m)]); x2=np.array([U[b[2][orb[k]]] for k in range(m)])
        e=(x1+x2)/2; s=x0-e; t=(x1-x2)/2
        Rs=s.sum(); Rt=t.sum()
        cnt[('ABCD'[X] if X==0 else 'free',typ(b),round(abs(Rs)**2,2),round(abs(Rt)**2,2))]+=1
for k,v in sorted(cnt.items()): print(k,v)
