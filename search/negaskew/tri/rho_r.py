"""Collision sampler of search/NOTES.md section 92: measure rho(r) = P_emp(sum of 4 PAF vectors = 0 on first r coordinates) / P_Gauss(r)  for the triangular pattern (k = 3 mod 8; for k = 7 mod 8 the pattern of tri2 mode 7):
A twisted-skew big (|R|^2=4v-3), B,C,D free with u0 = 1 and coset sum 0.   Empirical probability by hash collisions of N1 pair-sums (A,B) against N2 pair-sums (C,D).
HEURISTIC/VALIDATED: statistics of a random sample, not an exact count.
usage: python3 rho_r.py v k N1 N2 seed [rmax [perm_seed]]        (N1, N2 = numbers of pair sums; rmax = the largest number of coordinates)"""
import sys, math, time, numpy as np
from math import comb
from sympy import primitive_root
v=int(sys.argv[1]); k=int(sys.argv[2]); N1=int(float(sys.argv[3])); N2=int(float(sys.argv[4])); seed=int(sys.argv[5]); rmax=int(sys.argv[6]) if len(sys.argv)>6 else 99
perm_seed=int(sys.argv[7]) if len(sys.argv)>7 else 0
e=(v-1)//k; h=e//2; g=primitive_root(v); nc=e//2
MODE=3 if k%8==3 else 7
assert k%8 in (3,7)
rng=np.random.default_rng(seed)
orb=np.zeros(v,dtype=np.int64)
hh=pow(g,e,v); H=[pow(hh,i,v) for i in range(k)]
for j in range(e):
    r=pow(g,j,v)
    for x in H: orb[(r*x)%v]=j+1
Nt=np.zeros((nc,e+1,e+1))
for t in range(nc):
    s=pow(g,t,v)
    np.add.at(Nt[t],(orb,orb[(np.arange(v)+s)%v]),1)
I4=np.array([1,1j,-1,-1j])
def paf(U):   # U: (B,e+1) complex -> (B, 2*nc) ints (re_0,im_0,re_1,...)
    out=np.empty((U.shape[0],2*nc),dtype=np.int64)
    Uc=np.conj(U)
    for t in range(nc):
        P=((U@Nt[t])*Uc).sum(1)
        out[:,2*t]=np.rint(P.real); out[:,2*t+1]=np.rint(P.imag)
    return out
def sample_small(B):
    w=np.array([comb(h,a)**2 for a in range(h+1)],float); w/=w.sum()
    a=rng.choice(h+1,size=B,p=w); b=h-a
    # build multiset arrays: a copies of 0,2 ; b copies of 1,3
    idx=np.arange(e)[None,:]
    base=np.empty((B,e),dtype=np.int64)
    # positions 0..a-1 ->0 ; a..2a-1 ->2 ; 2a..2a+b-1 ->1 ; rest ->3
    base=np.where(idx<a[:,None],0,np.where(idx<2*a[:,None],2,np.where(idx<2*a[:,None]+b[:,None],1,3)))
    perm=np.argsort(rng.random((B,e)),axis=1)
    c=np.take_along_axis(base,perm,axis=1)
    ex=np.zeros((B,e+1),dtype=np.int64); ex[:,1:]=c
    return I4[ex]
def sample_big(B):
    u0=rng.integers(0,2,size=B); m=np.where(u0==0,1,-1)
    npl=(h+m)//2
    idx=np.arange(h)[None,:]
    base=(idx<npl[:,None])      # True -> plus pair
    perm=np.argsort(rng.random((B,h)),axis=1)
    plus=np.take_along_axis(base,perm,axis=1)
    c=np.where(plus,rng.integers(0,2,size=(B,h)),rng.integers(2,4,size=(B,h)))
    ex=np.zeros((B,e+1),dtype=np.int64); ex[:,0]=np.where(u0==0,0,3)
    ex[:,1:h+1]=c; ex[:,h+1:]=(1-c)%4
    return I4[ex]

def sample_askew_small(B):
    u0=rng.integers(0,2,size=B)
    npl=h//2
    idx=np.arange(h)[None,:]
    base=(idx<npl)
    perm=np.argsort(rng.random((B,h)),axis=1)
    plus=np.take_along_axis(np.broadcast_to(base,(B,h)),perm,axis=1)
    c=np.where(plus,rng.integers(0,2,size=(B,h)),rng.integers(2,4,size=(B,h)))
    ex=np.zeros((B,e+1),dtype=np.int64); ex[:,0]=np.where(u0==0,0,3)
    ex[:,1:h+1]=c; ex[:,h+1:]=(1-c)%4
    return I4[ex]
def sample_big1(B):
    # u0=1, S=1+i : n0=a+1,n2=a,n1=b+1,n3=b, a+b=h-1
    ws=np.array([math.exp(math.lgamma(e+1)-math.lgamma(a+2)-math.lgamma(a+1)-math.lgamma(h-a+1)-math.lgamma(h-a)) for a in range(h)])
    ws/=ws.sum(); a=rng.choice(h,size=B,p=ws); b=h-1-a
    idx=np.arange(e)[None,:]
    base=np.where(idx<(a+1)[:,None],0,np.where(idx<(2*a+1)[:,None],2,np.where(idx<(2*a+1+b+1)[:,None],1,3)))
    perm=np.argsort(rng.random((B,e)),axis=1)
    c=np.take_along_axis(base,perm,axis=1)
    ex=np.zeros((B,e+1),dtype=np.int64); ex[:,1:]=c
    return I4[ex]

# coordinate order and the dependent one
order=[]
for t in range(nc):
    if t<nc-1: order.append(2*t)
    order.append(2*t+1)
# order has 2*nc-1 = e-1 entries: re_0,im_0,...,re_{nc-2},im_{nc-2},im_{nc-1}
order=np.array(order)
if perm_seed:
    r2=np.random.default_rng(perm_seed); order=order[r2.permutation(len(order))]
nco=len(order); rmax=min(rmax,nco)
# --- Gaussian reference
M=200000
if MODE==3:
    T1,T2=sample_big,sample_small; nT1,nT2=1,3
else:
    T1,T2=sample_askew_small,sample_small
PA=np.concatenate([paf(T1(20000)) for _ in range(M//20000)])[:,order]
PS=np.concatenate([paf(T2(20000)) for _ in range(M//20000)])[:,order]
if MODE==3:
    mu=PA.mean(0)+3*PS.mean(0); C=np.cov(PA.T)+3*np.cov(PS.T)
else:
    PB=np.concatenate([paf(sample_big1(20000)) for _ in range(M//20000)])[:,order]
    mu=PA.mean(0)+PB.mean(0)+2*PS.mean(0); C=np.cov(PA.T)+np.cov(PB.T)+2*np.cov(PS.T)
# lattice of the sum projected on first r coords: HNF covolume via differences
from sympy import Matrix
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
def lattice_info(r):
    dA=[(PA[i]-PA[0])[:r] for i in range(1,3000)]; dS=[(PS[i]-PS[0])[:r] for i in range(1,3000)]
    dB=[(PB[i]-PB[0])[:r] for i in range(1,3000)] if MODE==7 else []
    B=int_basis(dA+dS+dB)
    Bm=Matrix(B); cov=math.sqrt(float((Bm*Bm.T).det())) if len(B)==r else float('nan')
    # coset test: is -(PA0+3PS0) in the lattice (+ base)?
    base0=(PA[0]+3*PS[0]) if MODE==3 else (PA[0]+PB[0]+2*PS[0]); x=list(map(int,-base0[:r])); ok=True
    for b in B:
        p=next(i for i,t in enumerate(b) if t)
        if x[p]%b[p]: ok=False;break
        q=x[p]//b[p]; x=[xi-q*bi for xi,bi in zip(x,b)]
    ok=ok and not any(x)
    return len(B),cov,ok
def gauss_P(r):
    Cr=C[:r,:r]; mr=mu[:r]; rk,cov,ok=lattice_info(r)
    if rk<r: return None,rk,cov,ok
    dens=math.exp(-0.5*mr@np.linalg.solve(Cr,mr))/((2*math.pi)**(r/2)*math.sqrt(np.linalg.det(Cr)))
    return dens*cov if ok else 0.0,rk,cov,ok
# --- sampling
def pairsum(B,T1,T2):
    return paf(T1(B))[:,order]+paf(T2(B))[:,order]
Rk=np.array([int(x)|1 for x in rng.integers(1,2**62,size=nco)],dtype=np.int64)
def stash(N,T1,T2):
    out=np.empty((N,nco),dtype=np.int16); CH=200000
    for lo in range(0,N,CH):
        B=min(CH,N-lo); out[lo:lo+B]=pairsum(B,T1,T2)
    return out
t0=time.time()
X1=stash(N1,sample_big if MODE==3 else sample_askew_small, sample_small if MODE==3 else sample_big1); X2=stash(N2,sample_small,sample_small)
print('sampled',N1,N2,'in',round(time.time()-t0),'s',flush=True)
def keyr(X,r):
    out=np.empty(X.shape[0],dtype=np.uint64); CH=2000000
    for lo in range(0,X.shape[0],CH):
        out[lo:lo+CH]=(X[lo:lo+CH,:r].astype(np.int64)*Rk[None,:r]).sum(1).view(np.uint64)
    return out
res=[]
for r in range(1,rmax+1):
    a=np.sort(keyr(X1,r)); b=(-keyr(X2,r).view(np.int64)).view(np.uint64)
    cnt=(np.searchsorted(a,b,side='right')-np.searchsorted(a,b,side='left')).sum()
    P=cnt/(N1*N2); G,rk,cov,ok=gauss_P(r)
    se=(1/math.sqrt(cnt) if cnt>0 else float('nan'))
    print(f'r={r:2d} coll={int(cnt):12d} P_emp={P:.4e} P_gauss={G:.4e} rho={P/G if G else float("nan"):.4f} (+-{se:.3f} rel) covol={cov:.1f} lat_ok={ok}',flush=True)
