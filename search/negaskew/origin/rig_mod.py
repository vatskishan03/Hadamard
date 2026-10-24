"""2-adic / modular test of the rigid (Lock+S) equation: count (A,C,D) with 2P_C+2P_D = T2 (mod M) vs exact.
 python3 rig_mod.py m M nsampleA seed   (q chosen = 4m-1 variant 1: cost (q+1)/2 ... we use whichever variant rigid_s_check uses with q=4m+3 (cost (q-3)/2 = 2m) )"""
import sys, itertools, time
import numpy as np
m=int(sys.argv[1]); M=int(sys.argv[2]); ns=int(sys.argv[3]); seed=int(sys.argv[4]); q=int(sys.argv[5]) if len(sys.argv)>5 else 4*m+3
U=np.array([1,1j,-1,-1j]); half=(m-1)//2; n=m-1
def paf_fft(S):
    F=np.fft.fft(S,axis=1); P=np.fft.ifft(np.abs(F)**2,axis=1); return np.conj(P[:,1:half+1])
idx=np.array(list(itertools.product(range(4),repeat=m)),dtype=np.int8)
seqs=U[idx]
P=paf_fft(seqs)
P2r=np.rint(2*P.real).astype(np.int64); P2i=np.rint(2*P.imag).astype(np.int64)
cost=(q-3)//2 if 2*((q-3)//2)==4*m else (q+1)//2
assert cost*2==4*m
if cost==(q-3)//2: eA=-(1+1j)/2; eB=(-1+1j)/2
else: eA=(1+1j)/2; eB=(1+1j)/2
SA=1-eA; SB=1-eB; lam=-2*np.conj(eA)*eB
def enc(R,I,Mod):
    key=np.zeros(R.shape[0],dtype=np.int64)
    for j in range(R.shape[1]):
        key=key*Mod+(R[:,j]%Mod)
    for j in range(I.shape[1]):
        key=key*Mod+(I[:,j]%Mod)
    return key
rng=np.random.default_rng(seed)
states=list(range(4**half)) if ns<=0 else rng.choice(4**half,size=min(ns,4**half),replace=False)
if M==0:
    keysC=None
keysC=enc(P2r,P2i,M if M else 1<<20) if M else None
tot=0; nz=0
if M:
    uk,cnt=np.unique(keysC,return_counts=True)
t0=time.time()
for a in states:
    ex=[(a//4**(half-1-j))%4 for j in range(half)]
    u=U[ex]
    sA=np.zeros(m,complex); sA[0]=SA; sA[1:half+1]=u
    for k in range(1,half+1): sA[m-k]=1j*np.conj(sA[k])
    sB=sA.copy(); sB[1:]*=lam; sB[0]=SB
    PA=paf_fft(sA[None,:])[0]+paf_fft(sB[None,:])[0]
    TR=np.rint(-2*PA.real).astype(np.int64); TI=np.rint(-2*PA.imag).astype(np.int64)
    if M:
        # need key(c)+key(d)=T mod M digitwise: enumerate unique keys c, compute needed d key
        # decode digits is expensive; instead recompute from P2 arrays on unique reps
        pass
    # vectorised: need = T - P2(C) for every C, then lookup count in table
    nr=(TR[None,:]-P2r); ni=(TI[None,:]-P2i)
    kneed=enc(nr,ni,M)
    pos=np.searchsorted(uk,kneed); pos[pos>=len(uk)]=0
    hit=(uk[pos]==kneed)
    c=int(np.sum(cnt[pos][hit]))
    tot+=c; nz+= (c>0)
print('m',m,'M',M,'q',q,'A states',len(states),'A with >0 solutions mod M',nz,'pairs count total',tot,'time',round(time.time()-t0,1),flush=True)
