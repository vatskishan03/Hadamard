"""list the DC patterns (R_A,R_B,R_C,R_D) classes of a DC-feasible shape: python3 dcpat.py v '[gens]'"""
import sys, json, os, numpy as np
from math import isqrt
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
from dcfeas import orbit_data, subset_sums
v=int(sys.argv[1]); gens=json.loads(sys.argv[2])
H,orb,reps,sz,neg=orbit_data(v,gens); n=len(sz); T=16*v; B=isqrt(T)+1
Z=subset_sums(sz[1:],start=(1,-1),M=v+5); Zc=Z[np.abs(Z)<=B]
r=np.array(np.meshgrid(Zc,Zc,Zc,Zc,indexing='ij')).reshape(4,-1).T
P=(r**2).sum(1); Q=r[:,0]*r[:,1]+r[:,1]*r[:,2]+r[:,2]*r[:,3]-r[:,0]*r[:,3]
ok=(P<=T)&(P*P-2*Q*Q>=0)&(P+np.sqrt(2)*Q>=-1e-9)&(P-np.sqrt(2)*Q>=-1e-9)
L=r[ok]; PL=P[ok]; QL=Q[ok]
pairs=[a for a in range(1,n) if a<neg[a]]; ps=[sz[a] for a in pairs]
R1={1,-1};R2={1,-1}
for s in ps: R1={x+t*s for x in R1 for t in (-2,0,2)}; R2={x+t*s for x in R2 for t in (-2,2)}
SA=[(1,a,b,a) for a in R1 for b in R2 if 1+2*a*a+b*b<=T]
from collections import defaultdict
byPQ=defaultdict(list)
for R,p,q in zip(L.tolist(),PL.tolist(),QL.tolist()): byPQ[(p,q)].append(tuple(R))
keys=list(byPQ)
res=[]
for A in SA:
    Pa=1+2*A[1]**2+A[2]**2; Qa=2*A[1]*A[2]
    for k1 in keys:
        for k2 in keys:
            if k2<k1: continue
            for k3 in keys:
                if k3<k2: continue
                if Pa+k1[0]+k2[0]+k3[0]==T and Qa+k1[1]+k2[1]+k3[1]==0: res.append((A,k1,k2,k3))
print('shape v=%d n=%d sizes=%s'%(v,n,sorted(set(sz))),' |Z|=%d |SL|=%d |SA|=%d'%(len(Zc),len(keys),len(SA)))
print('DC class patterns (R_A; (P,Q) of R_B<=R_C<=R_D):',len(res))
for A,k1,k2,k3 in res[:12]: print(' R_A=',A,' (P,Q):',k1,k2,k3,' e.g. R_B in',byPQ[k1][:2])
