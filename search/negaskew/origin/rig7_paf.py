import sys
sys.path.insert(0,'.')
from common import *
from collections import Counter
m=int(sys.argv[1]); fn=sys.argv[2]
S=Shape(m,[1]); orb=np.array(S.orb)
sols=load_solutions(fn)
def pafc(x): return np.array([np.sum(x*np.conj(np.roll(x,-k))) for k in range(m)])
cntpat=Counter()
for sol in sols[:int(sys.argv[3])]:
    row=[]
    for X,b in enumerate(sol):
        x0=np.array(U)[np.array(b[0])[orb]]; x1=np.array(U)[np.array(b[1])[orb]]; x2=np.array(U)[np.array(b[2])[orb]]
        e=(x1+x2)/2; s=x0-e; t=(x1-x2)/2
        ps=pafc(s)[1:4]; pt=pafc(t)[1:4]
        row.append((tuple(np.round(ps,2)),tuple(np.round(pt,2))))
    cntpat[tuple(r[0] for r in row)]+=1
for k,v in cntpat.most_common(6): print(v,k)
