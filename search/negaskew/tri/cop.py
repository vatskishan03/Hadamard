"""Gaussian-copula / Gamma-marginal correction of the hybrid estimate of hyb.py (k = 3 mod 8, search/NOTES.md section 92).  HEURISTIC model, validated against exact per-A counts.
usage: python3 cop.py v k [FILE]        (FILE: the SOL lines of ../qskew or of tri2 of the exact set, for the per-A comparison; needs scipy)"""
import sys, math, collections, numpy as np
from scipy import stats
import os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from tricore import *
from trilib import load_sols
def run(v,k,fn=None,Ns=400000,verbose=True,seed=5):
    F=Fam(v,k,seed); e=F.e; nL=comb(e,F.h)**2
    # --- hybrid (Gaussian) per-A estimate, as in hyb.py
    A=F.allA(); UA=np.array([F.blk2U(b) for b in A]); PA=F.paf(UA)[:,F.order]
    M=200000
    PS=np.concatenate([F.paf(F.s_small(20000)) for _ in range(M//20000)])[:,F.order]
    mu3=3*PS.mean(0); C3=3*np.cov(PS.T)
    B=EA.int_basis([(PS[i]-PS[0]) for i in range(1,3000)]); Bm=Matrix(B); covS=math.sqrt(float((Bm*Bm.T).det()))
    base3=3*PS[0]
    def in_coset(y):
        x=list(map(int,y-base3))
        for b in B:
            p=next(i for i,t in enumerate(b) if t)
            if x[p]%b[p]: return False
            q=x[p]//b[p]; x=[xi-q*bi for xi,bi in zip(x,b)]
        return not any(x)
    Ci=np.linalg.inv(C3); logdet=np.linalg.slogdet(C3)[1]; r=len(F.order)
    Ehyb=np.zeros(len(A))
    for i in range(len(A)):
        y=-PA[i]
        if in_coset(y):
            dz=y-mu3; Ehyb[i]=math.exp(-0.5*dz@Ci@dz-0.5*logdet-0.5*r*math.log(2*math.pi))*covS*nL**3
    # --- spectral free-triple sample
    Th=np.zeros((Ns,e)); CH=50000
    for lo in range(0,Ns,CH):
        B_=min(CH,Ns-lo)
        Th[lo:lo+B_]=sum(np.abs(F.spec(F.s_small(B_)))**2 for _ in range(3))
    m=e-1                                    # use modes 0..e-2
    X=Th[:,:m]; mu=X.mean(0); S=np.cov(X.T)
    # gamma marginals (moment match)
    kap=mu**2/X.var(0); th=X.var(0)/mu
    # KS check
    ks=max(stats.kstest(X[:20000,j],'gamma',args=(kap[j],0,th[j])).statistic for j in range(m))
    Z=np.empty_like(X)
    for j in range(m): Z[:,j]=stats.norm.ppf(np.clip(stats.gamma.cdf(X[:,j],kap[j],scale=th[j]),1e-12,1-1e-12))
    R=np.corrcoef(Z.T); Ri=np.linalg.inv(R); ldR=np.linalg.slogdet(R)[1]
    Si=np.linalg.inv(S); ldS=np.linalg.slogdet(S)[1]
    # targets
    FA=F.spec(UA); T=4*v-np.abs(FA)**2          # (|A|, e)
    Tm=T[:,:m]
    ratio=np.zeros(len(A)); alive=(Tm>0).all(1)
    for i in np.where(alive)[0]:
        t=Tm[i]
        lg=np.sum(stats.gamma.logpdf(t,kap,scale=th))
        z=stats.norm.ppf(np.clip(stats.gamma.cdf(t,kap,scale=th),1e-300,1-1e-16))
        logc=-0.5*ldR-0.5*z@(Ri-np.eye(m))@z
        lf_cop=lg+logc
        dz=t-mu; lf_g=-0.5*dz@Si@dz-0.5*ldS-0.5*m*math.log(2*math.pi)
        ratio[i]=math.exp(lf_cop-lf_g)
    Ecop=Ehyb*ratio
    out=dict(v=v,k=k,e=e,E_gauss_hyb=Ehyb.sum(),E_cop=Ecop.sum(),dead_by_positivity=int((~alive).sum()),A=len(A),ks=float(ks))
    if verbose: print(out)
    if fn:
        d2,sols=load_sols(fn,v,k); cnt=collections.Counter(s[0] for s in sols)
        ex=np.array([cnt.get(b,0) for b in A],float)
        print('exact total',ex.sum(),' E_hyb/exact %.3f  E_cop/exact %.3f'%(Ehyb.sum()/ex.sum(),Ecop.sum()/ex.sum()))
        idx=np.argsort(-Ecop); nb=8; sz=max(1,len(A)//nb)
        for q in range(0,len(A),sz):
            sel=idx[q:q+sz]; print('  cop-rank %4d..: mean E_cop %.1f  E_hyb %.1f  exact %.1f'%(q,Ecop[sel].mean(),Ehyb[sel].mean(),ex[sel].mean()))
        sel=(ex>0)&(Ecop>0)
        if sel.sum()>2: print('  log-corr(E_cop,exact) %.3f ; mean log10(exact/E_cop) %.3f ; A used but E_cop=0: %d'%(np.corrcoef(np.log(Ecop[sel]),np.log(ex[sel]))[0,1],np.mean(np.log10(ex[sel]/Ecop[sel])),int(((ex>0)&(Ecop==0)).sum())))
    return out
if __name__=="__main__":
    v=int(sys.argv[1]); k=int(sys.argv[2]); fn=sys.argv[3] if len(sys.argv)>3 else None
    run(v,k,fn)
