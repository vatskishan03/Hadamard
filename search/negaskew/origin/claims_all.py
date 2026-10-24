"""Definitive pass over ALL solutions of a file: exact checks of (J_k),(d_k),(Q), energies, constancy of Z1, mode z, pattern.
 python3 claims_all.py m q "[hK]" file"""
import sys, os
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
from common import *
from collections import Counter
m=int(sys.argv[1]); q=int(sys.argv[2]); hK=eval(sys.argv[3]); fn=sys.argv[4]
S=Shape(m,hK); orb=np.array(S.orb); nu=S.nu
sols=load_solutions(fn)
Uarr=np.array(U)
idx=(np.arange(m)[:,None]+np.arange(m)[None,:])%m
def paf(x): return (x[:,None]*np.conj(x[idx])).sum(axis=0)    # sum_j x(j) conj x(j+k)
def cc(a,b): return (a[:,None]*np.conj(b[idx])).sum(axis=0)
def typ(b):
    d=(b[2][0]-b[1][0])%4
    e=(U[b[1][0]]+U[b[2][0]])/2; x=U[b[0][0]]
    c=round((q-1)*abs(e)**2+2*(x*e.conjugate()).real)
    return {2:'E',0:'F',1:'M',3:'M'}[d]+str(c)
mx=0; modes=Counter(); pats=Counter(); offdef=0; en=Counter()
for sol in sols:
    Ps=np.zeros(m,complex); Pt=np.zeros(m,complex); Z1=np.zeros(m,complex); Q2=np.zeros(m,complex); Pe0=0
    ne=0
    for b in sol:
        x0=Uarr[np.array(b[0])[orb]]; x1=Uarr[np.array(b[1])[orb]]; x2=Uarr[np.array(b[2])[orb]]
        e=(x1+x2)/2; s=x0-e; t=(x1-x2)/2
        if np.any(np.abs(e[1:])>1e-9): offdef+=1
        Ps+=paf(s); Pt+=paf(t); Z1+=cc(e,s)+cc(s,e); Q2+=cc(t,s)-cc(s,t); ne+=np.sum(np.abs(e)**2)
    mx=max(mx,np.max(np.abs(Z1[1:]-Pt[1:])),np.max(np.abs(Ps[1:]+q*Z1[1:])),np.max(np.abs(Q2)))
    z1=Z1[1:]
    const=np.max(np.abs(z1-z1[0]))<1e-9
    modes[(const, round(-z1[0].real,3), round(z1[0].imag,3))]+=1
    pats[(typ(sol[0]),tuple(sorted(typ(b) for b in sol[1:])))]+=1
    en[(round(Ps[0].real,3),round(Pt[0].real,3),round(ne,3))]+=1
print('solutions',len(sols),' max violation of (J_k),(d_k),(Q):',mx,' blocks with defects off the origin:',offdef)
print('Z1 constant?, z, imag :',dict(modes))
print('energies (c_s, c_t, n/2):',dict(en))
for k,v in sorted(pats.items()): print(k,v)
