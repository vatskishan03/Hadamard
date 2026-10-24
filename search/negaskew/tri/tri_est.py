"""Gaussian+lattice estimate of the number of normalised triangular-pattern solutions (A big, B,C,D small with S=0, u0=1)
for v = 1 + k(k+1)/2 prime, k = 3 mod 8, H of order k.  HEURISTIC."""
import sys, math, numpy as np
from math import comb
from sympy import primitive_root, Matrix
from sympy.matrices.normalforms import hermite_normal_form
rng=np.random.default_rng(1)
def setup(v,k):
    g=primitive_root(v); e=(v-1)//k
    h=pow(g,e,v); H=[pow(h,i,v) for i in range(k)]
    orb=[0]*v; reps=[0]
    seen={0}
    idx={}
    for j in range(e):
        r=pow(g,j,v); reps.append(r)
        for x in H: orb[(r*x)%v]=j+1
    n=e+1
    # N_c(a,b) for c = 1..e (shift s = reps[c])
    N=np.zeros((n,n,n),dtype=np.int64)
    orb=np.array(orb)
    for c in range(1,n):
        s=reps[c]; np.add.at(N[c],(orb,orb[(np.arange(v)+s)%v]),1)
    neg=[0]+[((j+e//2)%e)+1 for j in range(e)]
    return dict(v=v,k=k,e=e,n=n,N=N,neg=neg,orb=orb)
def pafvec(d,ex):
    n=d['n']; N=d['N']; out=[]
    for c in range(1,n):
        if d['neg'][c]<c: continue
        Dm=(ex[:,None]-ex[None,:])%4
        M=[int(N[c][Dm==t].sum()) for t in range(4)]
        out+= [M[0]-M[2], M[1]-M[3]]
    return np.array(out,dtype=np.int64)
def sample_small(d,M):
    e=d['e']; out=[]
    while len(out)<M:
        # uniform member with sum 0: counts a,a,b,b with weight; draw random permutation of a multiset
        # choose a with prob C(e/2,a)^2 / C(e,e/2)
        h=e//2; w=np.array([comb(h,a)**2 for a in range(h+1)],float); a=rng.choice(h+1,p=w/w.sum()); b=h-a
        vals=np.array([0]*a+[2]*a+[1]*b+[3]*b); rng.shuffle(vals)
        out.append(np.concatenate([[0],vals]))
    return out
def sample_big(d,M):
    e=d['e']; h=e//2; out=[]
    for _ in range(M):
        u0=rng.integers(2); m=1 if u0==0 else -1
        # pairs j, j+h: c_j in {0,1,2,3}, partner 1-c; contribution c + i conj c -> +(1+i) for c in {0,1}, -(1+i) for {2,3}
        npl=(h+m)//2   # #plus - #minus = m
        signs=np.array([1]*npl+[-1]*(h-npl)); rng.shuffle(signs)
        c=np.where(signs>0, rng.integers(0,2,h), rng.integers(2,4,h))
        ex=np.zeros(e+1,dtype=np.int64); ex[0]=0 if u0==0 else 3
        ex[1:h+1]=c; ex[h+1:]=(1-c)%4
        out.append(ex)
    return out
def int_basis(vecs):
    # integer row echelon (Hermite-like) basis of the lattice spanned by vecs
    rows=[list(map(int,x)) for x in vecs]; dim=len(rows[0]); basis=[]
    col=0
    while rows and col<dim:
        nz=[r for r in rows if r[col]!=0]
        if not nz: col+=1; continue
        while len([r for r in rows if r[col]!=0])>1:
            nz=sorted([r for r in rows if r[col]!=0],key=lambda r:abs(r[col]))
            p=nz[0]
            for r in nz[1:]:
                q=r[col]//p[col]
                for t in range(dim): r[t]-=q*p[t]
            rows=[r for r in rows if any(r)]
        p=[r for r in rows if r[col]!=0][0]; basis.append(p); rows=[r for r in rows if r is not p]; col+=1
    return basis
def lattice_covol(diffs):
    B=Matrix(int_basis(diffs))
    G=B*B.T
    return B.shape[0], math.sqrt(float(G.det()))
def estimate(v,k,M=4000):
    d=setup(v,k); e=d['e']; h=e//2
    nA=2*comb(h,(h+1)//2)*2**h
    nS=comb(e,h)**2
    PA=np.array([pafvec(d,x) for x in sample_big(d,M)]); PS=np.array([pafvec(d,x) for x in sample_small(d,M)])
    mu=PA.mean(0)+3*PS.mean(0); C=np.cov(PA.T)+3*np.cov(PS.T)
    w,V=np.linalg.eigh(C); keep=w>1e-6*w.max(); r=int(keep.sum())
    Vk=V[:,keep]; mk=Vk.T@mu; null_resid=float(np.abs(V[:,~keep].T@mu).max()) if (~keep).any() else 0.0
    dens=-0.5*r*math.log(2*math.pi)-0.5*np.log(w[keep]).sum()-0.5*(mk**2/w[keep]).sum()
    diffs=[list(map(int,x-PA[0])) for x in PA[1:300]]+[list(map(int,x-PS[0])) for x in PS[1:300]]
    rk,cov=lattice_covol(diffs)
    lnE=math.log(nA)+3*math.log(nS)+dens+math.log(cov)
    return dict(v=v,k=k,e=e,orbits=e+1,log10_nA=math.log10(nA),log10_nS=math.log10(nS),rank=r,lat_rank=rk,null_resid=null_resid,covol=cov,log10E=lnE/math.log(10))
if __name__=="__main__":
    for v,k in [(67,11),(191,19)]+[(int(a),int(b)) for a,b in zip(sys.argv[1::2],sys.argv[2::2])]:
        print(estimate(v,k))
