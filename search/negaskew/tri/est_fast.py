"""Gaussian estimate (HEURISTIC; search/NOTES.md section 92) without the lattice step: covolume 2^(e-1) (PROOFS.md P7: equal to the sampled lattice for e <= 24).  Prints log10 E for every prime v = 1 + k(k+1)/2, k = 3 mod 4, k <= KMAX.
usage: python3 est_fast.py KMAX        (the table of section 92 is python3 est_fast.py 100)"""
import sys, math, time, numpy as np
from math import comb
from sympy import isprime, primitive_root
I4=np.array([1,1j,-1,-1j])
def est(v,k,M=40000,seed=1):
    e=(v-1)//k; h=e//2; g=primitive_root(v); nc=e//2; MODE=3 if k%8==3 else 7
    rng=np.random.default_rng(seed)
    orb=np.zeros(v,dtype=np.int64); hh=pow(g,e,v); H=[pow(hh,i,v) for i in range(k)]
    for j in range(e):
        r=pow(g,j,v)
        for x in H: orb[(r*x)%v]=j+1
    Nt=np.zeros((nc,e+1,e+1))
    for t in range(nc):
        s=pow(g,t,v); np.add.at(Nt[t],(orb,orb[(np.arange(v)+s)%v]),1)
    def paf(U):
        out=np.empty((U.shape[0],2*nc),dtype=np.int64); Uc=np.conj(U)
        for t in range(nc):
            P=((U@Nt[t])*Uc).sum(1); out[:,2*t]=np.rint(P.real); out[:,2*t+1]=np.rint(P.imag)
        return out
    def perm_rows(base): return np.take_along_axis(base,np.argsort(rng.random(base.shape),axis=1),axis=1)
    def s_small(B):
        lw=np.array([2*(math.lgamma(h+1)-math.lgamma(a+1)-math.lgamma(h-a+1)) for a in range(h+1)]); w=np.exp(lw-lw.max()); w/=w.sum()
        a=rng.choice(h+1,size=B,p=w); b=h-a; idx=np.arange(e)[None,:]
        base=np.where(idx<a[:,None],0,np.where(idx<2*a[:,None],2,np.where(idx<2*a[:,None]+b[:,None],1,3)))
        ex=np.zeros((B,e+1),dtype=np.int64); ex[:,1:]=perm_rows(base); return I4[ex]
    def s_askew(B,big):
        u0=rng.integers(0,2,size=B)
        npl=((h+np.where(u0==0,1,-1))//2) if big else np.full(B,h//2)
        idx=np.arange(h)[None,:]; plus=perm_rows(idx<npl[:,None])
        c=np.where(plus,rng.integers(0,2,size=(B,h)),rng.integers(2,4,size=(B,h)))
        ex=np.zeros((B,e+1),dtype=np.int64); ex[:,0]=np.where(u0==0,0,3); ex[:,1:h+1]=c; ex[:,h+1:]=(1-c)%4; return I4[ex]
    def s_big1(B):
        lws=np.array([math.lgamma(e+1)-math.lgamma(a+2)-math.lgamma(a+1)-math.lgamma(h-a+1)-math.lgamma(h-a) for a in range(h)]); ws=np.exp(lws-lws.max()); ws/=ws.sum()
        a=rng.choice(h,size=B,p=ws); b=h-1-a; idx=np.arange(e)[None,:]
        base=np.where(idx<(a+1)[:,None],0,np.where(idx<(2*a+1)[:,None],2,np.where(idx<(2*a+1+b+1)[:,None],1,3)))
        ex=np.zeros((B,e+1),dtype=np.int64); ex[:,1:]=perm_rows(base); return I4[ex]
    def samp(f):
        CH=min(5000,M); return np.concatenate([paf(f(CH)) for _ in range(M//CH)])
    order=[]
    for t in range(nc):
        if t<nc-1: order.append(2*t)
        order.append(2*t+1)
    order=np.array(order)
    lc=lambda a,b: math.lgamma(a+1)-math.lgamma(b+1)-math.lgamma(a-b+1)
    if MODE==3:
        types=[(lambda B:s_askew(B,True),1),(s_small,3)]
        logN=(math.log(2)+lc(h,(h+1)//2)+h*math.log(2)+3*2*lc(e,h))/math.log(10); mult=1
    else:
        types=[(lambda B:s_askew(B,False),1),(s_big1,1),(s_small,2)]
        logN=(math.log(2)+lc(h,h//2)+h*math.log(2)+lc(e,h+1)+lc(e,h)+2*2*lc(e,h))/math.log(10); mult=6
    Ps=[samp(f)[:,order] for f,_ in types]; mt=[m for _,m in types]
    mu=sum(m*P.mean(0) for P,m in zip(Ps,mt)); C=sum(m*np.atleast_2d(np.cov(P.T)) for P,m in zip(Ps,mt))
    sgn,ld=np.linalg.slogdet(C); r=len(order)
    ln_dens=-0.5*r*math.log(2*math.pi)-0.5*ld-0.5*mu@np.linalg.solve(C,mu)
    logP=(ln_dens+(e-1)*math.log(2))/math.log(10)
    return dict(v=v,k=k,e=e,mode=MODE,log10N=round(logN,3),log10P=round(float(logP),3),log10E=round(float(logN+logP+math.log10(mult)),3),sdmax=round(float(np.sqrt(np.linalg.eigvalsh(C).max())),1))
if __name__=="__main__":
    kmax=int(sys.argv[1]); 
    for k in range(3,kmax+1,4):
        v=1+k*(k+1)//2
        if not isprime(v): continue
        t=time.time(); d=est(v,k); d['sec']=round(time.time()-t); print(d,flush=True)
