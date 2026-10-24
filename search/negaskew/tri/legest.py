"""Gaussian estimate (HEURISTIC; search/NOTES.md section 92) of the number of solutions containing a Legendre-type flat free block (k = 3 mod 8, v = 3 mod 4):  A + C + D must have PAF sum +1 (re) / 0 (im).
usage: python3 legest.py        (members 67, 191, 379, 631, 947)"""
import sys, math, numpy as np
from math import comb
import os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import est_fast as EF
from tricore import Fam
def run(v,k,M=60000):
    F=Fam(v,k,seed=3); e=F.e; h=e//2; nc=e//2
    order=F.order
    def samp(f,M):
        return np.concatenate([F.paf(f(5000))[:,order] for _ in range(M//5000)])
    # A sampler
    def s_askew(B):
        rng=F.rng; u0=rng.integers(0,2,size=B); m=np.where(u0==0,1,-1); npl=(h+m)//2
        idx=np.arange(h)[None,:]; base=(idx<npl[:,None]); plus=np.take_along_axis(base,np.argsort(rng.random(base.shape),axis=1),axis=1)
        c=np.where(plus,rng.integers(0,2,size=(B,h)),rng.integers(2,4,size=(B,h)))
        ex=np.zeros((B,e+1),dtype=np.int64); ex[:,0]=np.where(u0==0,0,3); ex[:,1:h+1]=c; ex[:,h+1:]=(1-c)%4; return EF.I4[ex]
    PA=samp(s_askew,M); PS=samp(F.s_small,M)
    mu=PA.mean(0)+2*PS.mean(0); C=np.cov(PA.T)+2*np.cov(PS.T)
    # target: re=+1, im=0 on the kept coordinates (order lists re_t (t<nc-1), im_t)
    tgt=np.array([1 if (c%2==0) else 0 for c in order],dtype=float)
    dz=tgt-mu; r=len(order)
    ln=-0.5*dz@np.linalg.solve(C,dz)-0.5*np.linalg.slogdet(C)[1]-0.5*r*math.log(2*math.pi)+(e-1)*math.log(2)
    NA=2*comb(h,(h+1)//2)*2**h; NL=comb(e,h)**2
    E_one=NA*NL*NL*math.exp(ln)        # ordered (A, C, D) with the Legendre block fixed at B
    return e,E_one*3*2                  # 3 positions, 2 Legendre-type blocks
for v,k in [(67,11),(191,19),(379,27),(631,35),(947,43)]:
    e,E=run(v,k); print(v,'e',e,'Gaussian expected number of solutions containing a Legendre-type flat block: %.3g'%E)
