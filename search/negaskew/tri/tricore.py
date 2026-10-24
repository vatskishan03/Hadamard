"""The family sampler used by hyb.py, cop.py and legest.py (HEURISTIC estimates, search/NOTES.md section 92): class Fam(v, k, seed) holds the coset structure of Z_v (v = 1 + k(k+1)/2 prime, e = (v-1)/k cosets of the subgroup of order k),
the tensors Nt of the PAF of coset-constant blocks, the periods eta_j, the samplers of small free blocks (u0 = 1, coset sum 0), the list allA() of all twisted-skew A (k = 3 mod 8, big) and the spectra F(m) = u0 + sum_j c_j eta_{j+m}."""
import sys, os, math, itertools, collections, numpy as np
from math import comb
from sympy import Matrix, primitive_root
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from trilib import setup, I4
import est_all as EA
class Fam:
    def __init__(self,v,k,seed=5):
        self.v=v; self.k=k; d=setup(v,k); self.d=d; e=d['e']; self.e=e; h=e//2; self.h=h; self.g=d['g']; self.nc=e//2
        self.rng=np.random.default_rng(seed)
        orb=np.zeros(v,dtype=np.int64)
        for x in range(1,v): orb[x]=d['cos'][x]+1
        self.orb=orb
        Nt=np.zeros((self.nc,e+1,e+1))
        for t in range(self.nc):
            s=pow(self.g,t,v); np.add.at(Nt[t],(orb,orb[(np.arange(v)+s)%v]),1)
        self.Nt=Nt
        order=[]
        for t in range(self.nc):
            if t<self.nc-1: order.append(2*t)
            order.append(2*t+1)
        self.order=np.array(order)
        # periods eta_j = sum_{x in C_j} exp(2 pi i x / v)
        eta=np.zeros(e,dtype=complex)
        for x in range(1,v): eta[d['cos'][x]]+=np.exp(2j*np.pi*x/v)
        self.eta=eta
        self.E=np.array([[eta[(j+m)%e] for j in range(e)] for m in range(e)])   # F(m)=u0+ sum_j c_j eta_{j+m} = u0 + c @ E[m,:]
    def paf(self,U):
        nc=self.nc; out=np.empty((U.shape[0],2*nc),dtype=np.int64); Uc=np.conj(U)
        for t in range(nc):
            P=((U@self.Nt[t])*Uc).sum(1); out[:,2*t]=np.rint(P.real); out[:,2*t+1]=np.rint(P.imag)
        return out
    def spec(self,U):
        """U: (B,e+1) complex (u0, c_0..c_{e-1}); returns F (B,e): F[:,m]"""
        return U[:,[0]]+U[:,1:]@self.E.T
    def s_small(self,B):
        e=self.e;h=self.h
        w=np.array([comb(h,a)**2 for a in range(h+1)],float); w/=w.sum(); a=self.rng.choice(h+1,size=B,p=w); b=h-a
        idx=np.arange(e)[None,:]; base=np.where(idx<a[:,None],0,np.where(idx<2*a[:,None],2,np.where(idx<2*a[:,None]+b[:,None],1,3)))
        base=np.take_along_axis(base,np.argsort(self.rng.random(base.shape),axis=1),axis=1)
        ex=np.zeros((B,e+1),dtype=np.int64); ex[:,1:]=base; return I4[ex]
    def allA(self):
        e=self.e;h=self.h; out=[]
        for u0 in (0,3):
            m=1 if u0==0 else -1; npl=(h+m)//2
            for sgn in itertools.product((1,-1),repeat=h):
                if sum(1 for s in sgn if s>0)!=npl: continue
                for ch in itertools.product((0,1),repeat=h):
                    c=[0]*e
                    for i in range(h):
                        ci=ch[i] if sgn[i]>0 else 2+ch[i]; c[i]=ci; c[i+h]=(1-ci)%4
                    out.append((u0,tuple(c)))
        return out
    def blk2U(self,b):
        ex=np.zeros(self.e+1,dtype=np.int64); ex[0]=b[0]; ex[1:]=b[1]; return I4[ex]
