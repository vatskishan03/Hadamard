"""Validate the Spence construction for ODD |G| = v: four +-1 sequences on Z_v with sum PAF = -4 off 0 and row sums (+-1,+-1,+-1,+-1)
give a Hadamard matrix of order 4(v+1) = [[K, Q x 1^T],[P x 1, M]] with P = Q = S (the 4x4 Hadamard matrix of row sums), K = -S."""
import numpy as np, itertools, sys
def circ_paf(x):
    v=len(x); return np.array([int(x@np.roll(x,-s)) for s in range(v)])
def GS(A,B,C,D,v):
    R=np.zeros((v,v),dtype=np.int64); R[np.arange(v),(-np.arange(v))%v]=1  # group inversion on Z_v
    sub=lambda x: np.array([[x[(h-g)%v] for h in range(v)] for g in range(v)])   # X[g,h]=x(h-g)
    A_,B_,C_,D_=[sub(x) for x in (A,B,C,D)]
    M=np.block([[A_,B_@R,C_@R,D_@R],[-B_@R,A_,D_.T@R,-C_.T@R],[-C_@R,-D_.T@R,A_,B_.T@R],[-D_@R,C_.T@R,-B_.T@R,A_]])
    return M
def build(seqs,v):
    a,b,c,d=[int(s.sum()) for s in seqs]
    S=np.array([[a,b,c,d],[-b,a,d,-c],[-c,-d,a,b],[-d,c,-b,a]])
    assert (S@S.T==4*np.eye(4)).all() and np.abs(S).max()==1
    M=GS(*seqs,v)
    one=np.ones((1,v),dtype=np.int64)
    P=S.copy(); Q=S.copy(); K=-(Q@S.T@P)//4
    H=np.vstack([np.hstack([K,np.kron(Q,one)]),np.hstack([np.kron(P,one.T),M])])
    return H
def find(v, tries=2000000, seed=1):
    rng=np.random.default_rng(seed)
    # exhaustive-ish: enumerate all +-1 sequences with row sum -1 or +1 (up to negation, take -1), group by PAF
    k=(v-1)//2
    lst=[]
    for comb in itertools.combinations(range(v),k):
        x=-np.ones(v,dtype=np.int64); x[list(comb)]=1   # k pluses -> sum = 2k - v = -1
        lst.append(x)
    P=np.array([circ_paf(x) for x in lst])
    target=np.full(v,-4); target[0]=4*v
    # meet in the middle pairs
    from collections import defaultdict
    h=defaultdict(list)
    for i in range(len(lst)):
        for j in range(i,len(lst)):
            h[tuple(P[i]+P[j])].append((i,j))
    for key,v1 in h.items():
        need=tuple(target-np.array(key))
        if need in h:
            i,j=v1[0]; k_,l=h[need][0]
            return [lst[i],lst[j],lst[k_],lst[l]]
    return None
for v in [3,5,7,9,11,13]:
    sol=find(v)
    if sol is None: print(v,'no solution'); continue
    H=build(sol,v); n=4*(v+1)
    ok=(H@H.T==n*np.eye(n,dtype=np.int64)).all() and np.abs(H).max()==1
    print('v',v,'order',n,'HHT=nI',ok, 'rowsums',[int(s.sum()) for s in sol])
