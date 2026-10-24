import sys, os, itertools, json
sys.path.insert(0,'.')
from sstage import *
from common import load_solutions
import numpy as np
m=13; hK=[3]
epsA=0.5+0.5j; xA=1; epsB=-1; epsD=-1
sh=kshape([m],hK); nu=sh['nu']; neg=sh['neg']; N=np.array(sh['N'])
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
        q_=np.einsum('na,ab,nb->n',u[:,1:],M,np.conj(u[:,1:]))+S*np.conj(u[:,j])+np.conj(S)*u[:,neg[j]]
        out[:,j-1]=q_
    return out
PA,PB,PD,PC=pafs(uA,SA),pafs(uB,SB),pafs(uD,SD),pafs(uC,SC)
tot=PA+PB+PD
def key(P): return (np.rint(2*P.real).astype(np.int64), np.rint(2*P.imag).astype(np.int64))
Rt,It=key(-tot); Rc,Ic=key(PC)
from collections import defaultdict
d=defaultdict(list)
for i in range(len(allC)): d[(tuple(Rc[i]),tuple(Ic[i]))].append(i)
S_sols=[]
for i in range(ncomb):
    k=(tuple(Rt[i]),tuple(It[i]))
    for j in d.get(k,[]): S_sols.append((i,j))
print('S-solutions',len(S_sols))
# engine level-1 instances for this role pattern
sols=load_solutions('runs/299_all.out')
q=23
def eorig(b): return (U[b[1][0]]+U[b[2][0]])/2
eng=set()
for s in sols:
    ok = (s[0][0][0]==0 and abs(eorig(s[0])-(0.5+0.5j))<1e-9 and all(s[X][0][0]==0 for X in (1,2,3)) and abs(eorig(s[1])+1)<1e-9 and abs(eorig(s[2]))<1e-9 and abs(eorig(s[3])+1)<1e-9)
    if ok: eng.add(tuple(tuple(s[X][0][1:]) for X in range(4)))
print('engine level-1 instances for role (A12,F20,E,F20, x0=1):',len(eng))
def expo(u): return tuple(int(round(np.angle(x)/(np.pi/2)))%4 for x in u[1:])
Sset=set()
feat=[]
for (i,j) in S_sols:
    key_=(expo(uA[i]),expo(uB[i]),expo(uC if False else np.concatenate([[0],Uarr[allC[j]]])),expo(uD[i]))
    Sset.add(key_)
print('distinct (S) quadruples',len(Sset),' of which in engine set',len(Sset&eng),' engine not in S',len(eng-Sset))
pickle=__import__('pickle'); pickle.dump((Sset,eng),open('runs/sstage299.pkl','wb'))

# features: row sums of s
from collections import Counter
cf=Counter(); cn=Counter()
ell=sh['ell']
def rowsum(u,S): return S+sum(ell[o]*u[o] for o in range(1,nu))
for (i,j) in S_sols:
    uc=np.concatenate([[0],Uarr[allC[j]]])
    key_=(expo(uA[i]),expo(uB[i]),expo(uc),expo(uD[i]))
    rs=(rowsum(uA[i],SA),rowsum(uB[i],SB),rowsum(uc,SC),rowsum(uD[i],SD))
    f=tuple(round(abs(r)**2,2) for r in rs)
    f=(f[0],tuple(sorted([f[1],f[3]])),f[2])
    (cf if key_ in eng else cn)[f]+=1
print('completing:'); [print(k,v) for k,v in sorted(cf.items())]
print('failing:'); [print(k,v) for k,v in sorted(cn.items())]
