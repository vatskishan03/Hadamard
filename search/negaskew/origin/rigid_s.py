"""Exact (Lock+S) count for the rigid pattern (A,M|E,E), K = Z_m, H_K trivial (all orbits singletons), z = 0 forced when q > 4(m-2) + 2 sum|S_X|.
 s_A = [S_A, u_A],  u_A(-k) = i conj u_A(k);  B = lambda A on K\{0}, lambda = -2 conj(eps_A) eps_B;  C, D free unit sequences (E type, s(0) unit).
 count (A, C, D) with  P_{s_A} + P_{s_B} + P_{s_C} + P_{s_D} = 0  at all lags k = 1..(m-1)/2.   D taken modulo translation and global unit (the equations are invariant).
 usage: python3 rigid_s.py m q  [maxD]"""
import sys, itertools, time
import numpy as np
U=np.array([1,1j,-1,-1j])
m=int(sys.argv[1]); q=int(sys.argv[2])
half=(m-1)//2
def paf_all(seqs):   # seqs (N, m) complex ; returns (N, half) complex  P(k)=sum_j s(j) conj s(j+k)
    out=np.zeros((seqs.shape[0],half),complex)
    for k in range(1,half+1):
        out[:,k-1]=(seqs*np.conj(np.roll(seqs,-k,axis=1))).sum(axis=1)
    return out
def lin(P2):  # P2 integer (N,half) complex scaled by 2 -> int64 keys, base 64
    re=np.rint(P2.real).astype(np.int64); im=np.rint(P2.imag).astype(np.int64)
    key=np.zeros(P2.shape[0],dtype=np.int64)
    for i in range(half):
        key+=re[:,i]*(64**(2*i))+im[:,i]*(64**(2*i+1))
    return key
t0=time.time()
# E-type sequences: s(0)=1, u(1..m-1) free : 4^(m-1)
n=m-1
idx=np.array(list(itertools.product(range(4),repeat=n)),dtype=np.int8)
seqs=np.ones((idx.shape[0],m),complex); seqs[:,1:]=U[idx]
PC=paf_all(seqs)          # (N, half) complex integers
keyC=lin(2*PC)
order=np.argsort(keyC); keyC_sorted=keyC[order]
# canonical D reps modulo translation (+ unit renormalisation): translation by a: s'(j)=s(j-a); renormalise s'(0)=1 -> multiply by conj(s(-a))
def canon_mask(seqs):
    N=seqs.shape[0]
    # integer code of a normalised sequence: base-4 exponents of entries 1..m-1
    def code(S):
        e=np.rint(np.angle(S[:,1:])/(np.pi/2)).astype(np.int64)%4
        c=np.zeros(N,dtype=np.int64)
        for j in range(n): c=c*4+e[:,j]
        return c
    c0=code(seqs); best=c0.copy()
    for a in range(1,m):
        S=np.roll(seqs,a,axis=1)
        S=S*np.conj(S[:,[0]])
        best=np.minimum(best,code(S))
    return c0==best
mask=canon_mask(seqs)
Dseq=seqs[mask]; PD=PC[mask]; keyD=lin(2*PD)
print('E sequences',seqs.shape[0],'canonical D reps',Dseq.shape[0],'time',round(time.time()-t0,1),flush=True)
# A origin variants: cost (q-3)/2 or (q+1)/2 with 2 blocks A,B s.t. cA+cB = 4m
res=[]
def origin_variants(cost,skew):
    out=[]
    for x0 in ([1,-1j] if skew else [1]):
        for x1 in U:
            for x2 in U:
                if skew and abs(x2-1j*np.conj(x1))>1e-9: continue
                e=(x1+x2)/2
                c=(q-1)*abs(e)**2+2*(x0*np.conj(e)).real
                if abs(c-cost)<1e-9: out.append((x0,e))
    return out
cases=[]
for cA in ((q+1)//2,(q-3)//2):
    for cB in ((q+1)//2,(q-3)//2):
        if cA+cB==4*m: cases.append((cA,cB))
print('cost cases',cases,flush=True)
total=0
for cA,cB in cases:
    vA=origin_variants(cA,True); vB=origin_variants(cB,False)
    # distinct (S,eps) pairs only
    seen=set(); pairsA=[]
    for x0,e in vA:
        k=(round(x0.real),round(x0.imag),round(e.real*2),round(e.imag*2))
        if k not in seen: seen.add(k); pairsA.append((x0,e))
    seen=set(); pairsB=[]
    for x0,e in vB:
        k=(round(x0.real),round(x0.imag),round(e.real*2),round(e.imag*2))
        if k not in seen: seen.add(k); pairsB.append((x0,e))
    # A skew sequences: u(1..half) free, u(m-k)= i conj u(k)
    aidx=np.array(list(itertools.product(range(4),repeat=half)),dtype=np.int8)
    for (x0A,eA) in pairsA:
        for (x0B,eB) in pairsB:
            SA=x0A-eA; SB=x0B-eB
            lam=-2*np.conj(eA)*eB
            sA=np.zeros((aidx.shape[0],m),complex); sA[:,0]=SA
            sA[:,1:half+1]=U[aidx]
            for k in range(1,half+1): sA[:,m-k]=1j*np.conj(sA[:,k])
            sB=sA.copy(); sB[:,1:]*=lam; sB[:,0]=SB
            PA=paf_all(sA)+paf_all(sB)           # (nA, half)
            T2=np.rint(-2*PA)                    # need 2P_C + 2P_D = T2
            cnt=0; ex=None
            keyT=lin(T2)
            for a in range(aidx.shape[0]):
                need=keyT[a]-keyD                 # linear keys of 2P_C
                lo=np.searchsorted(keyC_sorted,need,side='left'); hi=np.searchsorted(keyC_sorted,need,side='right')
                c=hi-lo
                nz=np.nonzero(c)[0]
                if len(nz):
                    # exact verification (guards against base-64 aliasing)
                    for d in nz:
                        for ci in order[lo[d]:hi[d]]:
                            tot=PA[a]+PC[ci]+PD[d]
                            if np.max(np.abs(tot))<1e-9: cnt+=1; ex=ex or (a,d,ci)
            print('costs',cA,cB,'A(x0,eps)',(x0A,eA),'B eps',eB,'(Lock+S) solutions with D canonical:',cnt,'time',round(time.time()-t0,1),flush=True)
            total+=cnt
print('TOTAL (D mod translation)',total)
