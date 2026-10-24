"""S-stage solutions + t row-sum feasibility given R^s.  python3 sstage3.py m "[hK]" q epsA xA epsB epsD"""
import sys, os, itertools, json
sys.path.insert(0,'.')
from sstage import *
from rowfeas import gaussian_vals
import numpy as np
from collections import defaultdict, Counter
def sstage_solutions(m,hK,epsA,xA,epsB,epsD):
    sh=kshape([m],hK); nu=sh['nu']; neg=sh['neg']; N=np.array(sh['N']); ell=sh['ell']
    pairs=[o for o in range(1,nu) if o<neg[o]]; p=len(pairs)
    LS=local_states(epsA,epsB,epsD); nL=len(LS)
    idx=np.array(list(itertools.product(range(nL),repeat=p))); ncomb=len(idx)
    def build(which):
        u=np.zeros((ncomb,nu),complex)
        for t,o in enumerate(pairs):
            po=neg[o]
            a0=np.array([LS[i][which][0] for i in range(nL)]); a1=np.array([LS[i][which][1] for i in range(nL)])
            u[:,o]=a0[idx[:,t]]; u[:,po]=a1[idx[:,t]]
        return u
    uA,uB,uD=build(0),build(1),build(2)
    allC=np.array(list(itertools.product(range(4),repeat=nu-1)))
    Uarr=np.array(U); uC=np.zeros((len(allC),nu),complex); uC[:,1:]=Uarr[allC]
    SA=xA-epsA; SB=1-epsB; SD=1-epsD; SC=1.0
    def pafs(u,S):
        out=np.zeros((u.shape[0],nu-1),complex)
        for j in range(1,nu):
            M=N[j][1:,1:].astype(float)
            out[:,j-1]=np.einsum('na,ab,nb->n',u[:,1:],M,np.conj(u[:,1:]))+S*np.conj(u[:,j])+np.conj(S)*u[:,neg[j]]
        return out
    PA,PB,PD,PC=pafs(uA,SA),pafs(uB,SB),pafs(uD,SD),pafs(uC,SC)
    tot=PA+PB+PD
    def key(P): return (np.rint(2*P.real).astype(np.int64), np.rint(2*P.imag).astype(np.int64))
    Rt,It=key(-tot); Rc,Ic=key(PC)
    d=defaultdict(list)
    for i in range(len(allC)): d[(tuple(Rc[i]),tuple(Ic[i]))].append(i)
    sols=[]
    for i in range(ncomb):
        for j in d.get((tuple(Rt[i]),tuple(It[i])),[]): sols.append((i,j))
    rs=lambda u,S: S+sum(ell[o]*u[o] for o in range(1,nu))
    out=[]
    for (i,j) in sols:
        out.append((rs(uA[i],SA),rs(uB[i],SB),rs(uC[j],SC),rs(uD[i],SD)))
    return out,p,ell
def t_feasible(Rs,taus,h,p,ct,eps_dummy=None):
    """exists R^t_X in tau_X + h w (w in (1+i)Z[i], |w|<=2p; A: (1-i)J) with sum|R^t|^2=ct and sum Im(R^s conj R^t)=0 ; taus: list of possible tau per block"""
    ws=gaussian_vals(0,2*p); Js=[J for J in range(-p,p+1) if (J-p)%2==0]
    sets=[]
    for X in range(4):
        s=set()
        for tau in taus[X]:
            Rts=[tau+h*(1-1j)*J for J in Js] if X==0 else [tau+h*w for w in ws]
            for rt in Rts: s.add((int(round(4*abs(rt)**2)), int(round(4*(Rs[X]*rt.conjugate()).imag))))
        sets.append(s)
    L={(a[0]+b[0],a[1]+b[1]) for a in sets[0] for b in sets[1]}
    R={(a[0]+b[0],a[1]+b[1]) for a in sets[2] for b in sets[3]}
    tgt=(int(round(4*ct)),0)
    return any((tgt[0]-x[0],tgt[1]-x[1]) in R for x in L)
if __name__=='__main__':
    m=int(sys.argv[1]); hK=json.loads(sys.argv[2]); q=int(sys.argv[3]); c=lambda s: complex(s.replace(' ',''))
    epsA,xA,epsB,epsD=c(sys.argv[4]),c(sys.argv[5]),c(sys.argv[6]),c(sys.argv[7])
    sols,p,ell=sstage_solutions(m,hK,epsA,xA,epsB,epsD)
    h=ell[1]
    n2=abs(epsA)**2+abs(epsB)**2+abs(epsD)**2
    ct=4*m-n2
    # tau options: A: +-i eps ; F: 0 ; E: any unit
    taus=[[1j*epsA,-1j*epsA],[0],[1,1j,-1,-1j],[0]]
    # B and D are F; C is E
    taus=[[1j*epsA,-1j*epsA],[0.0],[1,1j,-1,-1j],[0.0]]
    ok=0; seen=Counter()
    cache={}
    for Rs in sols:
        k=tuple((round(r.real,3),round(r.imag,3)) for r in Rs)
        if k not in cache: cache[k]=t_feasible(Rs,taus,h,p,ct)
        ok+=cache[k]
    print('m',m,'h',h,'p',p,'(S)-solutions',len(sols),' row-sum(t,R3)-feasible',ok)
