# The odd-|G| Spence construction rebuilt independently of ../oddspence/oddspence_test.py (independent code of the judge of search/NOTES.md section 92). For cyclic Z_v (v odd): four +-1 sequences with
# sum PAF = -4 off 0; H = [[-S, S x 1^T],[S x 1, M]], M = Goethals-Seidel array with circulant blocks; checks H H^T = 4(v+1) I for v = 3, 5, 7, 9, 11.
import numpy as np, itertools, sys, collections
def circ(a):
    v=len(a); return np.array([[a[(j-i)%v] for j in range(v)] for i in range(v)],dtype=np.int64)
def find(v):
    X=np.array(list(itertools.product((1,-1),repeat=v)),dtype=np.int64)
    # PAF all shifts
    P=np.stack([(np.roll(X,-s,axis=1)*X).sum(1) for s in range(1,v)],1)
    U,inv=np.unique(P,axis=0,return_inverse=True)
    first={}
    for i,q in enumerate(inv.ravel()): first.setdefault(q,i)
    rng=np.random.default_rng(2); w=rng.integers(1,2**60,size=v-1,dtype=np.int64)
    h=U@w; tgt=-4*w.sum()
    hp=(h[:,None]+h[None,:])
    d={}
    L=len(U)
    for i in range(L):
        for j in range(i,L):
            d.setdefault(int(hp[i,j]),[]).append((i,j))
    for key,lst in d.items():
        need=int(tgt-key)
        if need in d:
            for (i,j) in lst:
                for (k,l) in d[need]:
                    if np.array_equal(U[i]+U[j]+U[k]+U[l],np.full(v-1,-4)):
                        return [X[first[q]] for q in (i,j,k,l)]
    return None
def build(seqs):
    v=len(seqs[0]); A,B,C,D=[circ(s) for s in seqs]
    a,b,c,d=[int(s.sum()) for s in seqs]
    assert a*a+b*b+c*c+d*d==4 and all(abs(x)==1 for x in (a,b,c,d)), (a,b,c,d)
    R=np.eye(v,dtype=np.int64)[::-1]
    M=np.block([[A,B@R,C@R,D@R],[-B@R,A,D.T@R,-C.T@R],[-C@R,-D.T@R,A,B.T@R],[-D@R,C.T@R,-B.T@R,A]])
    S=np.array([[a,b,c,d],[-b,a,d,-c],[-c,-d,a,b],[-d,c,-b,a]],dtype=np.int64)
    ones=np.ones((v,1),dtype=np.int64)
    top=np.hstack([-S,np.kron(S,ones.T)])
    bot=np.hstack([np.kron(S,ones),M])
    return np.vstack([top,bot])
for v in (3,5,7,9,11):
    seqs=find(v)
    if seqs is None: print('v=%d: no quadruple'%v); continue
    H=build(seqs); n=H.shape[0]
    print('v=%d order %d: entries +-1 %s, H H^T = %d I: %s, rowsums %s'%(v,n,bool((np.abs(H)==1).all()),n,bool(np.array_equal(H@H.T,n*np.eye(n,dtype=np.int64))),[int(s.sum()) for s in seqs]))
