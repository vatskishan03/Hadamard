import sys
sys.path.insert(0,'.')
from common import *
from collections import Counter
m=int(sys.argv[1]); q=int(sys.argv[2]); hK=eval(sys.argv[3]); fn=sys.argv[4]
S=Shape(m,hK); orb=S.orb
sols=load_solutions(fn)
def pafc(x): return np.array([np.sum(x*np.conj(np.roll(x,-k))) for k in range(m)])
def cc(a,b): return np.array([np.sum(a*np.conj(np.roll(b,-k))) for k in range(m)])
cnt=Counter(); cs=Counter()
step=max(1,len(sols)//2000)
for sol in sols[::step]:
    Z1=np.zeros(m,complex); Ps=np.zeros(m,complex)
    for b in sol:
        x0=np.array([U[b[0][orb[k]]] for k in range(m)]); x1=np.array([U[b[1][orb[k]]] for k in range(m)]); x2=np.array([U[b[2][orb[k]]] for k in range(m)])
        e=(x1+x2)/2; s=x0-e
        Z1+=cc(e,s)+cc(s,e); Ps+=pafc(s)
    # per class lags
    vals=tuple(sorted(set((round(Z1[k].real,3),round(Z1[k].imag,3)) for k in range(1,m))))
    cnt[vals]+=1
for k,v in sorted(cnt.items(), key=lambda kv:-kv[1])[:10]: print(v,k)
