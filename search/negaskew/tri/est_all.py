"""Gaussian + lattice estimate (HEURISTIC; search/NOTES.md section 92) of the normalised number of solutions of the triangular family for prime v = 1 + k(k+1)/2 (est_fast.py is the faster version without the lattice step).
mode 3 (k=3 mod 8): A big skew, B,C,D small free.   mode 7 (k=7 mod 8): A small skew, one big free block (S=1+i; x3 positions x2 for S=1-i), two small free.
Prints log10 N_tot, log10 P_gauss, log10 E (the total normalised count: comparable with qskew's 'COUNT').
usage: python3 est_all.py v,k [v,k ...]        (e.g. 67,11 191,19)"""
import sys, math, time, numpy as np
from math import comb
from sympy import primitive_root, Matrix
def int_basis(vecs):
    rows=[list(map(int,x)) for x in vecs]; dim=len(rows[0]); basis=[]; col=0
    while rows and col<dim:
        while len([r for r in rows if r[col]!=0])>1:
            nz=sorted([r for r in rows if r[col]!=0],key=lambda r:abs(r[col])); p=nz[0]
            for r in nz[1:]:
                q=r[col]//p[col]
                for t in range(dim): r[t]-=q*p[t]
            rows=[r for r in rows if any(r)]
        nzr=[r for r in rows if r[col]!=0]
        if not nzr: col+=1; continue
        p=nzr[0]; basis.append(p); rows=[r for r in rows if r is not p]; col+=1
    return basis
I4=np.array([1,1j,-1,-1j])
def estimate(v,k,M=100000,seed=1,nlat=1500):
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
    def perm_rows(base):
        return np.take_along_axis(base,np.argsort(rng.random(base.shape),axis=1),axis=1)
    def s_small(B):
        w=np.array([comb(h,a)**2 for a in range(h+1)],float); w/=w.sum(); a=rng.choice(h+1,size=B,p=w); b=h-a
        idx=np.arange(e)[None,:]; base=np.where(idx<a[:,None],0,np.where(idx<2*a[:,None],2,np.where(idx<2*a[:,None]+b[:,None],1,3)))
        ex=np.zeros((B,e+1),dtype=np.int64); ex[:,1:]=perm_rows(base); return I4[ex]
    def s_askew(B,big):
        u0=rng.integers(0,2,size=B)
        if big: m=np.where(u0==0,1,-1); npl=(h+m)//2
        else: npl=np.full(B,h//2)
        idx=np.arange(h)[None,:]; base=(idx<npl[:,None]); plus=perm_rows(base)
        c=np.where(plus,rng.integers(0,2,size=(B,h)),rng.integers(2,4,size=(B,h)))
        ex=np.zeros((B,e+1),dtype=np.int64); ex[:,0]=np.where(u0==0,0,3); ex[:,1:h+1]=c; ex[:,h+1:]=(1-c)%4; return I4[ex]
    def s_big1(B):
        ws=np.array([math.exp(math.lgamma(e+1)-math.lgamma(a+2)-math.lgamma(a+1)-math.lgamma(h-a+1)-math.lgamma(h-a)) for a in range(h)]); ws/=ws.sum()
        a=rng.choice(h,size=B,p=ws); b=h-1-a; idx=np.arange(e)[None,:]
        base=np.where(idx<(a+1)[:,None],0,np.where(idx<(2*a+1)[:,None],2,np.where(idx<(2*a+1+b+1)[:,None],1,3)))
        ex=np.zeros((B,e+1),dtype=np.int64); ex[:,1:]=perm_rows(base); return I4[ex]
    def samp(f,*a):
        return np.concatenate([paf(f(20000,*a)) for _ in range(M//20000)])
    order=[]
    for t in range(nc):
        if t<nc-1: order.append(2*t)
        order.append(2*t+1)
    order=np.array(order)
    if MODE==3:
        types=[(lambda B:s_askew(B,True),1),(s_small,3)]
        logN=math.log10(2*comb(h,(h+1)//2)*2**h)+3*math.log10(comb(e,h)**2)
        mult=1
    else:
        types=[(lambda B:s_askew(B,False),1),(s_big1,1),(s_small,2)]
        logN=math.log10(2*comb(h,h//2)*2**h)+math.log10(comb(e,h+1)*comb(e,h))+2*math.log10(comb(e,h)**2)
        mult=6
    Ps=[samp(f)[:,order] for f,_ in types]; mult_t=[m for _,m in types]
    mu=sum(m*P.mean(0) for P,m in zip(Ps,mult_t)); C=sum(m*np.atleast_2d(np.cov(P.T)) for P,m in zip(Ps,mult_t))
    r=len(order)
    dl=[]
    for P in Ps: dl+= [(P[i]-P[0]) for i in range(1,nlat)]
    B=int_basis(dl); rk=len(B)
    Bm=Matrix(B); cov=math.sqrt(float((Bm*Bm.T).det()))
    x=list(map(int,-sum(m*P[0] for P,m in zip(Ps,mult_t)))); ok=True
    for b in B:
        p=next(i for i,t in enumerate(b) if t)
        if x[p]%b[p]: ok=False;break
        q=x[p]//b[p]; x=[xi-q*bi for xi,bi in zip(x,b)]
    ok=ok and not any(x)
    w,V=np.linalg.eigh(C); keep=w>1e-6*w.max(); rr=int(keep.sum())
    Vk=V[:,keep]; mk=Vk.T@mu
    ln_dens=-0.5*rr*math.log(2*math.pi)-0.5*np.log(w[keep]).sum()-0.5*(mk**2/w[keep]).sum()
    logP=(ln_dens+math.log(cov))/math.log(10)
    return dict(v=v,k=k,e=e,mode=MODE,log10N=logN,rank=rk,dens_rank=rr,log10covol=math.log10(cov),coset_ok=ok,log10P=float(logP),log10E=float(logN+logP+math.log10(mult)),sdmin=float(np.sqrt(w[keep].min())),sdmax=float(np.sqrt(w.max())))
if __name__=="__main__":
    for a in sys.argv[1:]:
        v,k=map(int,a.split(',')); t=time.time(); d=estimate(v,k); d['sec']=round(time.time()-t); print(d,flush=True)
