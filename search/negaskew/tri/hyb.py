"""Hybrid estimate (k = 3 mod 8, search/NOTES.md section 92; HEURISTIC): exact A-list, Gaussian+lattice density of the free triple.  Per-A prediction vs exact per-A counts (when a solution file is given), and the share of the estimated mass carried by the A blocks with a large spectral peak.
usage: python3 hyb.py v k [FILE]        (FILE: the SOL lines of ../qskew or of tri2 of the exact set)"""
import sys, math, itertools, collections, numpy as np
from math import comb
import os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from trilib import *
import est_all as EA
from sympy import Matrix
v=int(sys.argv[1]); k=int(sys.argv[2]); fn=sys.argv[3] if len(sys.argv)>3 else None
d=setup(v,k); e=d['e']; h=e//2; g=d['g']; nc=e//2
rng=np.random.default_rng(5)
orb=np.zeros(v,dtype=np.int64)
for x in range(1,v): orb[x]=d['cos'][x]+1
Nt=np.zeros((nc,e+1,e+1))
for t in range(nc):
    s=pow(g,t,v); np.add.at(Nt[t],(orb,orb[(np.arange(v)+s)%v]),1)
def paf(U):
    out=np.empty((U.shape[0],2*nc),dtype=np.int64); Uc=np.conj(U)
    for t in range(nc):
        P=((U@Nt[t])*Uc).sum(1); out[:,2*t]=np.rint(P.real); out[:,2*t+1]=np.rint(P.imag)
    return out
order=[]
for t in range(nc):
    if t<nc-1: order.append(2*t)
    order.append(2*t+1)
order=np.array(order)
# all A in list
allA=[]
for u0 in (0,3):
    m=1 if u0==0 else -1; npl=(h+m)//2
    for sgn in itertools.product((1,-1),repeat=h):
        if sum(1 for s in sgn if s>0)!=npl: continue
        for ch in itertools.product((0,1),repeat=h):
            c=[0]*e
            for i in range(h):
                ci=ch[i] if sgn[i]>0 else 2+ch[i]; c[i]=ci; c[i+h]=(1-ci)%4
            allA.append((u0,tuple(c)))
def blk2U(b):
    ex=np.zeros(e+1,dtype=np.int64); ex[0]=b[0]; ex[1:]=b[1]; return I4[ex]
UA=np.array([blk2U(b) for b in allA]); PA=paf(UA)[:,order]
# small sampler
def s_small(B):
    w=np.array([comb(h,a)**2 for a in range(h+1)],float); w/=w.sum(); a=rng.choice(h+1,size=B,p=w); b=h-a
    idx=np.arange(e)[None,:]; base=np.where(idx<a[:,None],0,np.where(idx<2*a[:,None],2,np.where(idx<2*a[:,None]+b[:,None],1,3)))
    base=np.take_along_axis(base,np.argsort(rng.random(base.shape),axis=1),axis=1)
    ex=np.zeros((B,e+1),dtype=np.int64); ex[:,1:]=base; return I4[ex]
M=200000
PS=np.concatenate([paf(s_small(20000)) for _ in range(M//20000)])[:,order]
muS=PS.mean(0); CS=np.cov(PS.T)
mu3=3*muS; C3=3*CS
# lattice of small blocks
B=EA.int_basis([(PS[i]-PS[0]) for i in range(1,3000)]); rk=len(B); Bm=Matrix(B); covS=math.sqrt(float((Bm*Bm.T).det()))
base3=3*PS[0]
def in_coset(y):
    x=list(map(int,y-base3))
    for b in B:
        p=next(i for i,t in enumerate(b) if t)
        if x[p]%b[p]: return False
        q=x[p]//b[p]; x=[xi-q*bi for xi,bi in zip(x,b)]
    return not any(x)
Ci=np.linalg.inv(C3); logdet=np.linalg.slogdet(C3)[1]; r=len(order)
nL=comb(e,h)**2
pred=np.zeros(len(allA)); okc=np.zeros(len(allA),bool)
for i in range(len(allA)):
    y=-PA[i]
    okc[i]=in_coset(y)
    dz=y-mu3
    pred[i]=math.exp(-0.5*dz@Ci@dz-0.5*logdet-0.5*r*math.log(2*math.pi))*covS*nL**3 if okc[i] else 0.0
print('v',v,'k',k,'|A|',len(allA),'A passing lattice coset test:',int(okc.sum()),' E_hybrid = %.1f (log10 %.3f)'%(pred.sum(),math.log10(pred.sum())))
if fn:
    d2,sols=load_sols(fn,v,k)
    cnt=collections.Counter(s[0] for s in sols)
    ex=np.array([cnt.get(b,0) for b in allA],float)
    print('exact total',ex.sum(),' ratio hyb/exact',pred.sum()/ex.sum())
    print('A used by exact solutions:',int((ex>0).sum()),' lattice-compatible A:',int(okc.sum()),' used & compat:',int(((ex>0)&okc).sum()))
    # per-A comparison by decile of predicted
    idx=np.argsort(-pred); 
    for q in range(0,len(allA),max(1,len(allA)//8)):
        sel=idx[q:q+max(1,len(allA)//8)]
        print(' pred-rank block %4d..: mean pred %.1f mean exact %.1f'%(q,pred[sel].mean(),ex[sel].mean()))
    # correlation in log
    sel=(ex>0)&(pred>0); print('log-corr(pred,exact) over A with both >0: %.3f; mean log10(exact/pred) %.3f'%(np.corrcoef(np.log(pred[sel]),np.log(ex[sel]))[0,1],np.mean(np.log10(ex[sel]/pred[sel]))))
# --- decomposition of the hybrid mass by spectral peak of A
def peak(b):
    u=expand(d,b); F=np.fft.fft(u); g_=d['g']
    P=np.array([abs(F[pow(g_,m,v)])**2 for m in range(e)])
    return P.max()/v
pk=np.array([peak(b) for b in allA])
print('peak/v quantiles of the A list:',np.round(np.quantile(pk,[0,0.1,0.25,0.5,0.75,0.9,1]),2))
for thr in (2.5,2.8,3.0,3.2,3.5):
    print(' share of hybrid mass from A with peak > %.1f v: %.3f  (share of A list %.3f)'%(thr,pred[pk>thr].sum()/pred.sum(),(pk>thr).mean()))
if fn:
    print(' exact: max peak/v among used A: %.3f ; min among unused: %.3f'%(pk[ex>0].max(),pk[ex==0].min()))
for thr in (3.6,3.8,4.0):
    print(' share of hybrid mass from A with peak > %.1f v: %.3f (share of A list %.3f)'%(thr,pred[pk>thr].sum()/pred.sum(),(pk>thr).mean()))
