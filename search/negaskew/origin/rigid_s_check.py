"""independent re-implementation (FFT PAF, python-tuple dict lookup, FULL D, no canonicalisation, no linear keys) of the rigid (Lock+S) count for sampled A-states.
 python3 rigid_s_check.py m q nsample seed"""
import sys, itertools, time
import numpy as np
m=int(sys.argv[1]); q=int(sys.argv[2]); ns=int(sys.argv[3]); seed=int(sys.argv[4])
U=np.array([1,1j,-1,-1j]); half=(m-1)//2; n=m-1
def paf_fft(S):    # S (N,m) -> P(k)=sum_j s(j) conj s(j+k), k=1..half   via FFT:  P = IFFT( conj(F) * F ) with the right index convention
    F=np.fft.fft(S,axis=1)
    P=np.fft.ifft(np.abs(F)**2,axis=1)       # P[k] = sum_j conj(s(j)) s(j+k)   (= conj of ours)
    return np.conj(P[:,1:half+1])
idx=np.array(list(itertools.product(range(4),repeat=n)),dtype=np.int8)
seqs=np.ones((idx.shape[0],m),complex); seqs[:,1:]=U[idx]
P=paf_fft(seqs)
P2=np.rint(2*P.real).astype(np.int64), np.rint(2*P.imag).astype(np.int64)
keys=[tuple(r) for r in np.concatenate([P2[0],P2[1]],axis=1)]
from collections import defaultdict
tab=defaultdict(int)
for kk in keys: tab[kk]+=1
rng=np.random.default_rng(seed)
# one origin variant: A: x0=1, eps=(-1-i)/2 if cost (q-3)/2 else (1+i)/2 ; B: x0=1 eps_B=(-1+i)/2 / (1+i)/2
cost=(q-3)//2 if 2*((q-3)//2)==4*m else (q+1)//2
if cost==(q-3)//2: eA=-(1+1j)/2; eB=(-1+1j)/2
else: eA=(1+1j)/2; eB=(1+1j)/2
SA=1-eA; SB=1-eB; lam=-2*np.conj(eA)*eB
tot=0; t0=time.time()
states=range(4**half) if ns<=0 else rng.choice(4**half,size=ns,replace=False)
for a in states:
    ex=[(a//4**(half-1-j))%4 for j in range(half)]
    u=U[ex]
    sA=np.zeros(m,complex); sA[0]=SA; sA[1:half+1]=u
    for k in range(1,half+1): sA[m-k]=1j*np.conj(sA[k])
    sB=sA.copy(); sB[1:]*=lam; sB[0]=SB
    PA=paf_fft(sA[None,:])[0]+paf_fft(sB[None,:])[0]
    T2=np.rint(-2*PA).astype(np.int64)    # need 2P_C + 2P_D = T2 ; with D over all seqs
    cnt=0
    # vectorised over D: need = T2 - 2P_D
    needR=T2.real.astype(np.int64)[None,:]-P2[0]; needI=T2.imag.astype(np.int64)[None,:]-P2[1] if False else None
    TR=np.rint((-2*PA).real).astype(np.int64); TI=np.rint((-2*PA).imag).astype(np.int64)
    needR=TR[None,:]-P2[0]; needI=TI[None,:]-P2[1]
    need=np.concatenate([needR,needI],axis=1)
    # only consider D whose need is within the range of C PAF values (|.|<=2m)
    ok=np.all(np.abs(need)<=2*m,axis=1)
    for d in np.nonzero(ok)[0]:
        cnt+=tab.get(tuple(need[d]),0)
    tot+=cnt
print('m',m,'q',q,'sampled A states',len(states),'full-D (Lock+S) count',tot,'time',round(time.time()-t0,1))
