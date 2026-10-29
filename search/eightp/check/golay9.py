# Complex (quaternary) Golay pairs of length n = 2..9 by exhaustive search over 4^n sequences (independent code, search/NOTES.md section 92): none of length 9, which disposes of the 'closure of the
# complex Golay lengths under products' route of ../census/recursive8p.py
import numpy as np, itertools, sys, collections
def golay(n):
    T=np.array(list(itertools.product(range(4),repeat=n)),dtype=np.int8)
    X=np.array([1,0,-1,0])[T]; Y=np.array([0,1,0,-1])[T]
    cols=[]
    for s in range(1,n):
        # sum_j a_{j+s} conj(a_j)
        re=(X[:,s:]*X[:,:n-s]+Y[:,s:]*Y[:,:n-s]).sum(1); im=(Y[:,s:]*X[:,:n-s]-X[:,s:]*Y[:,:n-s]).sum(1)
        cols+= [re,im]
    V=np.stack(cols,1).astype(np.int64)
    w=np.random.default_rng(3).integers(1,2**61,size=V.shape[1],dtype=np.int64)
    h=V@w
    d=collections.defaultdict(list)
    for i,x in enumerate(h): d[int(x)].append(i)
    cnt=0
    for i,x in enumerate(h):
        for j in d.get(int(-x),[]):
            if np.array_equal(V[i],-V[j]): cnt+=1
    return cnt
for n in (2,3,4,5,6,7,8,9):
    print(n, golay(n))
