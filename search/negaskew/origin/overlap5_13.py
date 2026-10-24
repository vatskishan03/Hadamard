import sys
sys.path.insert(0,'.')
from common import *
def primroot(p):
    for g in range(2,p):
        if len({pow(g,k,p) for k in range(p-1)})==p-1: return g
def classmap(m,hK):
    S=Shape(m,hK); g=primroot(m); ind={}; x=1
    for k in range(m-1): ind[x]=k; x=x*g%m
    cls={}
    for k in range(1,m): cls[S.orb[k]]=ind[k]%4
    return S,cls
def classvals(sol,cls):
    out=[]
    for b in sol:
        arr=[None]*4
        for o,cl in cls.items():
            arr[cl]=(b[0][o],b[1][o])    # (x0,x1) ; x2=-x1 off origin
        out.append(tuple(arr))
    return out
def canon(blocks):
    best=None
    for r in range(4):
        for sg in (1,-1):
            nb=[]
            for arr in blocks:
                nb.append(tuple(arr[(sg*j+r)%4] for j in range(4)))
            # unit multiplication of a block: multiply x0 and x1 by i^a (global unit of block) is a symmetry of free blocks; normalise x0 at class 0 -> 0
            nb2=[]
            for arr in nb:
                a=arr[0][0]
                nb2.append(tuple(((x0-a)%4,(x1-a)%4) for (x0,x1) in arr))
            key=tuple(sorted(nb2))
            if best is None or key<best: best=key
    return best
sets={}
for nm,(m,hK,fn) in {'5_11':(5,[1],'runs/5_11_full.out'),'13_27':(13,[3],'runs/13_27_all.out'),'13_23':(13,[3],'runs/299_all.out'),'13_19':(13,[3],'runs/13_19_all.out')}.items():
    S,cls=classmap(m,hK)
    sols=load_solutions(fn)
    K=set()
    for s in sols: K.add(canon(classvals(s,cls)))
    sets[nm]=K
    print(nm,len(sols),'solutions ->',len(K),'canonical class-value quadruples')
import itertools
for a,b in itertools.combinations(sets,2):
    print(a,b,'common',len(sets[a]&sets[b]))
