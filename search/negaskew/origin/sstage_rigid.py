"""(Lock+S) count for the rigid pattern (A,M | E,E): A skew M-type, B M-type free (x0=1), C,D E-type.  K=Z_m, H_K.
 python3 sstage_rigid.py m "[hK]" q  -> counts for z=0 (sum P_s = 0) and z=1 (sum P_s = q) over all origin variants of the identity (A cost, B cost with A+B = 4m)"""
import sys, os, itertools, json
sys.path.insert(0,'.'); sys.path.insert(0,'..')
import numpy as np
from lockcount import states
from qshape import kshape
U=[1,1j,-1,-1j]
def counts(m,hK,q,epsA,xA,epsB,z,maxC=4**11):
    sh=kshape([m],hK); nu=sh['nu']; neg=sh['neg']; N=np.array(sh['N'])
    pairs=[o for o in range(1,nu) if o<neg[o]]; p=len(pairs)
    LS=states([epsA,epsB],['A','F'],z); nL=len(LS)
    if nL==0: return 0
    idx=np.array(list(itertools.product(range(nL),repeat=p))); nc=len(idx)
    uA=np.zeros((nc,nu),complex); uB=np.zeros((nc,nu),complex)
    for t,o in enumerate(pairs):
        po=neg[o]
        a0=np.array([LS[i][0][0] for i in range(nL)]); a1=np.array([LS[i][0][1] for i in range(nL)])
        b0=np.array([LS[i][1][0] for i in range(nL)]); b1=np.array([LS[i][1][1] for i in range(nL)])
        uA[:,o]=a0[idx[:,t]]; uA[:,po]=a1[idx[:,t]]; uB[:,o]=b0[idx[:,t]]; uB[:,po]=b1[idx[:,t]]
    SA=xA-epsA; SB=1-epsB; SC=1.0
    def pafs(u,S):
        out=np.zeros((u.shape[0],nu-1),complex)
        for j in range(1,nu):
            M=N[j][1:,1:].astype(float)
            out[:,j-1]=np.einsum('na,ab,nb->n',u[:,1:],M,np.conj(u[:,1:]))+S*np.conj(u[:,j])+np.conj(S)*u[:,neg[j]]
        return out
    PA=pafs(uA,SA)+pafs(uB,SB)
    allC=np.array(list(itertools.product(range(4),repeat=nu-1)))
    uC=np.zeros((len(allC),nu),complex); uC[:,1:]=np.array(U)[allC]
    PC=pafs(uC,SC)
    # need P_C + P_D = q z - PA   (k != 0)
    key=lambda P:(np.rint(2*P.real).astype(np.int64),np.rint(2*P.imag).astype(np.int64))
    Rc,Ic=key(PC)
    # hash C by tuple
    from collections import defaultdict
    d=defaultdict(int)
    for i in range(len(allC)): d[(tuple(Rc[i]),tuple(Ic[i]))]+=1
    tot=0
    Rd,Id=Rc,Ic
    for a in range(nc):
        tr=np.rint(2*(q*z-PA[a]).real).astype(np.int64); ti=np.rint(2*(q*z-PA[a]).imag).astype(np.int64)
        # for each D: need C = target - D
        needR=tr[None,:]-Rd; needI=ti[None,:]-Id
        for i in range(len(allC)):
            tot+=d.get((tuple(needR[i]),tuple(needI[i])),0)
    return tot
if __name__=='__main__':
    m=int(sys.argv[1]); hK=json.loads(sys.argv[2]); q=int(sys.argv[3])
    Ae=[(1+1j)/2,-(1+1j)/2]
    # A cost (q+1)/2 or (q-3)/2 ; B (M-type free x0=1): cost (q+1)/2: e=(1+-i)/2 ; (q-3)/2: e=(-1+-i)/2 ; identity cA+cB=4m
    cA_opts={ (q+1)//2:(0.5+0.5j), (q-3)//2:-(0.5+0.5j)}
    cB_opts={ (q+1)//2:[(1+1j)/2,(1-1j)/2], (q-3)//2:[(-1+1j)/2,(-1-1j)/2]}
    for cA,eA in cA_opts.items():
        for cB,eBs in cB_opts.items():
            if cA+cB!=4*m: continue
            for eB in eBs:
                for z in (0,1):
                    n=counts(m,hK,q,eA,1,eB,z)
                    print('m',m,'q',q,'costs',cA,cB,'epsB',eB,'z',z,'(Lock+S) solutions',n,flush=True)
