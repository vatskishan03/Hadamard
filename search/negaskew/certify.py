import json,sys,hashlib,os
import numpy as np
from qdrv import *
v=int(sys.argv[1]); gens=json.loads(sys.argv[2]); sol=sys.argv[3].split()
d=shape(v,gens)
E=[expand(d,[int(c) for c in p]) for p in sol]
def paf_all(e):
    u=np.array([1j**x for x in e])
    return np.array([np.sum(u[(np.arange(v)+s)%v]*np.conj(u)) for s in range(v)])
tot=sum(paf_all(e) for e in E)
print('G-level sum PAF: peak',tot[0],' max|off-peak|',np.max(np.abs(tot[1:])))
uA=np.array([1j**x for x in E[0]])
print('A twisted-skew:',all(abs(uA[(-g)%v]-1j*np.conj(uA[g]))<1e-9 for g in range(1,v)),' u(0)=',uA[0])
H=real_matrix(v,E); print('ring-realisation real matrix: order',H.shape[0],check_real(H))
# negacyclic form
n=2*v
def to_nega(e):
    a=[0]*n
    for b in range(v):
        w=(1+1j)*[1,1j,-1,-1j][e[b]%4]; wp=w*[1,-1j,-1,1j][b%4]; wp=complex(round(wp.real),round(wp.imag))
        iv=[1,1j,-1,-1j][v%4]   # i^v exactly
        # a_b + a_{b+v} i^v = w'  ; i^v = +-i
        if iv==-1j: a[b]=int(wp.real); a[b+v]=int(-wp.imag)
        else: a[b]=int(wp.real); a[b+v]=int(wp.imag)
    return a
S=[to_nega(e) for e in E]
assert all(set(s)<={1,-1} for s in S)
a=S[0]; print('A palindromic',all(a[j]==a[n-j] for j in range(1,n)),'a0',a[0])
cert=dict(order=8*v,construction="goethals-seidel-negacyclic",skew=True,group=f"negacyclic Z[x]/(x^{n}+1) = Z[i][Z_{v}]",v=v,n=n,
          rowsums=[int(sum(s)) for s in S],sequences=["".join('+' if x>0 else '-' for x in s) for s in S],
          quaternary_orbit_exponents=sol,multiplier_generators=gens)
os.makedirs('out',exist_ok=True)
json.dump(cert,open(f'out/cert_{8*v}.json','w'))
print('rowsums',cert['rowsums'])
