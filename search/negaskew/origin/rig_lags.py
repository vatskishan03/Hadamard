"""where does the rigid (Lock+S) count collapse?  count (C,D) [C,D unit seqs, D full] with P_C+P_D = T on the first j lags, j=1..half; m, q, nA, seed"""
import sys, itertools, time
import numpy as np
m=int(sys.argv[1]); q=int(sys.argv[2]); ns=int(sys.argv[3]); seed=int(sys.argv[4])
U=np.array([1,1j,-1,-1j]); half=(m-1)//2
def paf_fft(S):
    F=np.fft.fft(S,axis=1); P=np.fft.ifft(np.abs(F)**2,axis=1); return np.conj(P[:,1:half+1])
idx=np.array(list(itertools.product(range(4),repeat=m)),dtype=np.int8)
seqs=U[idx]; del idx
P=paf_fft(seqs); del seqs
P2r=np.rint(2*P.real).astype(np.int64); P2i=np.rint(2*P.imag).astype(np.int64); del P
cost=(q-3)//2 if 2*((q-3)//2)==4*m else (q+1)//2
if cost==(q-3)//2: eA=-(1+1j)/2; eB=(-1+1j)/2
else: eA=(1+1j)/2; eB=(1+1j)/2
SA=1-eA; SB=1-eB; lam=-2*np.conj(eA)*eB
rng=np.random.default_rng(seed)
states=rng.choice(4**half,size=ns,replace=False)
base=4*m+1; off=2*m
def enc(R,I,j):
    key=np.zeros(R.shape[0],dtype=np.int64)
    for l in range(j):
        key=key*base+(R[:,l]+off)
        key=key*base+(I[:,l]+off)
    return key
tabs={}
for j in range(1,half+1):
    uk,cnt=np.unique(enc(P2r,P2i,j),return_counts=True); tabs[j]=(uk,cnt)
res=np.zeros((ns,half),dtype=np.int64)
for ia,a in enumerate(states):
    ex=[(a//4**(half-1-j))%4 for j in range(half)]
    u=U[ex]
    sA=np.zeros(m,complex); sA[0]=SA; sA[1:half+1]=u
    for k in range(1,half+1): sA[m-k]=1j*np.conj(sA[k])
    sB=sA.copy(); sB[1:]*=lam; sB[0]=SB
    PA=paf_fft(sA[None,:])[0]+paf_fft(sB[None,:])[0]
    TR=np.rint(-2*PA.real).astype(np.int64); TI=np.rint(-2*PA.imag).astype(np.int64)
    for j in range(1,half+1):
        nr=(TR[None,:j]-P2r[:,:j]); ni=(TI[None,:j]-P2i[:,:j])
        ok=np.all(np.abs(nr)<=off,axis=1)&np.all(np.abs(ni)<=off,axis=1)
        if not ok.any(): continue
        kneed=enc(nr[ok],ni[ok],j)
        uk,cnt=tabs[j]
        pos=np.searchsorted(uk,kneed); pos[pos>=len(uk)]=0
        hit=(uk[pos]==kneed)
        res[ia,j-1]=int(np.sum(cnt[pos][hit]))
print('m',m,'q',q,'per-A pair counts matching the first j lags (j=1..%d):'%half)
for ia in range(ns): print(list(res[ia]))
print('mean',res.mean(axis=0))
