"""Lists for the n=955 structure G = Z2^3 x Z29, H = <(tau, 16)>, T = nontrivial chars of Z2^3 (pattern 2)."""
import itertools, numpy as np, math, time, sys
from math import comb
from zq import *
p=29; u=16; g=2
tau3 = mat([[0,0,1],[1,0,1],[0,1,0]])
Tn = [k for k in itertools.product(range(2),repeat=3) if any(k)]
H = block_diag(tau3, mat([[u]]))
st = ZS([2,2,2,p], [H], phi=lambda c: c[:, :3] % 2, T=Tn, qmods=[2,2,2], name='955S1')
O = [[pow(g,t+4*i,p) for i in range(7)] for t in range(4)]
ordz = {}
for t in range(4):
    for i,z in enumerate(O[t]): ordz[z]=(t,i)
uinv = pow(u,-1,p)
# e_a = tau^a e0
e0 = np.array([0,0,1]); E=[]
v=e0.copy()
for a in range(7):
    E.append(v.copy()); v=(tau3@v)%2
assert len({tuple(x) for x in E})==7
eidx=[int(x[0]*4+x[1]*2+x[2]) for x in E]
# permutation: for each a and z: source index z' = u^{-a} z
uia=[pow(uinv,a,p) for a in range(7)]
def build(f0,f1):
    """f0,f1: arrays (K,29) of +-1/0 ; returns X (K,232) int8 with X[e,z]"""
    K=f0.shape[0]
    X=np.zeros((K,8,p),dtype=np.int8)
    X[:,0,:]=f0
    zs=np.arange(p)
    for a in range(7):
        src=(zs*uia[a])%p
        X[:,eidx[a],:]=f1[:,src]
    return X.reshape(K,8*p)
Hperm=st.perms[0]
def check_inv(X):
    return (X[:,Hperm]==X).all()
# ---- A-type enumeration (all 2^17)
def enum_A():
    out=[]
    bits7=np.array(list(itertools.product((1,-1),repeat=7)),dtype=np.int8)   # 128 x 7
    N=2*4*128*128
    f0=np.zeros((N,p),dtype=np.int8); f1=np.zeros((N,p),dtype=np.int8)
    k=0
    idx=np.arange(N)
    a_s=np.where(idx%2==0,1,-1)
    al0=np.where((idx//2)%2==0,1,-1); al1=np.where((idx//4)%2==0,1,-1)
    b0=bits7[(idx//8)%128]; b1=bits7[(idx//1024)%128]
    # f0
    for t in range(4):
        al = al0 if t%2==0 else al1
        for z in O[t]: f0[:,z]=al
    f0[:,0]=0
    f1[:,0]=a_s
    for i,z in enumerate(O[0]): f1[:,z]=b0[:,i]
    for i,z in enumerate(O[1]): f1[:,z]=b1[:,i]
    for i in range(7):
        f1[:,O[2][(i+3)%7]]=b0[:,i]
        f1[:,O[3][(i+3)%7]]=b1[:,i]
    return f0,f1
def zq_compress(X):
    K=X.shape[0]; R=X.reshape(K,8,p).astype(np.int64)
    ze=R.sum(1)             # sum over e for each z -> (K,29)
    qe=R.sum(2)             # sum over z for each e -> (K,8)
    z5=np.stack([ze[:,0]]+[ze[:,O[t][0]] for t in range(4)],1)
    return z5, qe[:,0], qe[:,eidx[0]]
if __name__=="__main__":
    t0=time.time()
    f0,f1=enum_A()
    XA=build(f0,f1)
    print("A enumerated",XA.shape,"inv",check_inv(XA[:2000]),time.time()-t0)
    # symmetric check
    neg=st.neg
    print("A symmetric:", (XA[:2000][:,neg]==XA[:2000]).all())
    z5,q0,q1=zq_compress(XA)
    print("Q compress sets:",set(zip(q0.tolist(),q1.tolist())))
