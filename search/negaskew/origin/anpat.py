"""pattern / decoupling analysis of a solution file: python3 anpat.py m q "[hK]" solfile"""
import sys
sys.path.insert(0,'.')
from common import *
from collections import Counter
m=int(sys.argv[1]); q=int(sys.argv[2]); hK=eval(sys.argv[3]); fn=sys.argv[4]
S=Shape(m,hK); orb=S.orb; nu=S.nu
sols=load_solutions(fn)
def fullfun(b):
    x0=np.array([U[b[0][orb[k]]] for k in range(m)]); x1=np.array([U[b[1][orb[k]]] for k in range(m)]); x2=np.array([U[b[2][orb[k]]] for k in range(m)])
    e=(x1+x2)/2; s=x0-e; t=(x1-x2)/2
    return s,t,e
def pafc(x): return np.array([np.sum(x*np.conj(np.roll(x,-k))) for k in range(m)])
def cc(a,b): return np.array([np.sum(a*np.conj(np.roll(b,-k))) for k in range(m)])
def typ(b):
    d=(b[2][0]-b[1][0])%4
    e=(U[b[1][0]]+U[b[2][0]])/2; x=U[b[0][0]]
    c=round((q-1)*abs(e)**2+2*(x*e.conjugate()).real)
    return {2:'E',0:'F',1:'M',3:'M'}[d]+str(c)
res=Counter()
step=max(1,len(sols)//3000)
for sol in sols[::step]:
    Ps=np.zeros(m,complex); Pt=np.zeros(m,complex); Z1=np.zeros(m,complex)
    offdef=[]
    for X,b in enumerate(sol):
        for o in range(1,nu):
            if (b[2][o]-b[1][o])%4!=2: offdef.append((X,o,(b[2][o]-b[1][o])%4))
        s,t,e=fullfun(b)
        Ps+=pafc(s); Pt+=pafc(t); Z1+=cc(e,s)+cc(s,e)
    dec=np.max(np.abs(Ps[1:]))<1e-9
    res[(typ(sol[0]),tuple(sorted(typ(b) for b in sol[1:])),bool(dec),len(offdef),round(Ps[0].real,2),round(Pt[0].real,2))]+=1
for k,v in sorted(res.items()): print(k,v)
