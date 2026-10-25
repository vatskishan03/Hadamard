"""Anatomy of the order-850 core of hadamard/6796: each sequence as a function of (ind x, ind z) for x in GF(25)
(coordinates (x0, x1) -> x0 + x1 t, t^2 = 2t + 1, as the certificate's symmetry M) and z in Z_17 (primitive root 3)."""
import json, os, numpy as np, itertools, sys
CERT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'hadamard', '6796', 'certificate.json')
c=json.load(open(CERT))['core']
S=[np.array([0 if ch=='0' else (1 if ch=='+' else -1) for ch in s]) for s in c['sequences']]
p=5; q2=17
# GF(25) = GF(5)[t]/(t^2-2t-1); element (a,b) = a + b t
def mul(u,v):
    a,b=u; c_,d=v
    # (a+bt)(c+dt) = ac + (ad+bc)t + bd t^2, t^2 = 2t+1
    return ((a*c_+b*d)%p, (a*d+b*c_+2*b*d)%p)
def order(u):
    x=u;k=1
    while x!=(1,0): x=mul(x,u);k+=1
    return k
gen=None
for a in range(5):
    for b in range(5):
        if (a,b)!=(0,0) and order((a,b))==24: gen=(a,b);break
    if gen: break
ind25={}
x=(1,0)
for k in range(24): ind25[x]=k; x=mul(x,gen)
# t = (0,1): index of M
print('gen',gen,'ind(t)=',ind25[(0,1)], 'ind(-1)=',ind25[(4,0)], 'ind(2)=',ind25[(2,0)])
g17=3; ind17={pow(g17,k,17):k for k in range(16)}
def idx(x0,x1,z): return x0*85+x1*17+z
names='ABCD'
for f,name in zip(S,names):
    print('==',name,'sum',f.sum())
    ax1=[(ind25[(x0,x1)], f[idx(x0,x1,0)]) for x0 in range(5) for x1 in range(5) if (x0,x1)!=(0,0)]
    ax1.sort(); print(' axis GF25 (ind mod 24):', ''.join('+' if v>0 else ('-' if v<0 else '0') for _,v in ax1))
    ax2=sorted((ind17[z], f[idx(0,0,z)]) for z in range(1,17))
    print(' axis Z17 (ind mod 16):', ''.join('+' if v>0 else ('-' if v<0 else '0') for _,v in ax2))
    T=np.zeros((24,16),int)
    for x0 in range(5):
        for x1 in range(5):
            if (x0,x1)==(0,0): continue
            for z in range(1,17):
                T[ind25[(x0,x1)], ind17[z]] = f[idx(x0,x1,z)]
    # check dependence on (i mod a, j mod b) and on i - k j combos
    for a in (2,4,8,12,24):
        for b in (2,4,8,16):
            ok=all(T[i,j]==T[i%a,j%b] for i in range(24) for j in range(16))
            if ok: print(' interior depends only on (ind25 mod %d, ind17 mod %d)'%(a,b)); break
        else: continue
        break
    # dependence on linear combos s*i + t*j mod m
    for m in (2,4,8):
        for s in range(m):
            for t in range(m):
                d={}
                ok=True
                for i in range(24):
                    for j in range(16):
                        key=(s*i+t*j)%m
                        if d.setdefault(key,T[i,j])!=T[i,j]: ok=False;break
                    if not ok: break
                if ok and (s or t): print(' interior = f((%d*i+%d*j) mod %d):'%(s,t,m), [d[k] for k in sorted(d)])
