"""anatomy of the real skew solution of order 2156 (K = Z7^2, q = 11, H_K = <(2,2)>): defect support and costs per block, the K-level identities.  usage: python3 real2156.py [HIT.json]   (default data/c539_hit.json, the HIT file of the run that found skew/2156; a replay of search/skewplug/data/repro539 writes the same file to runs/)"""
import sys, json, numpy as np
d=json.load(open(sys.argv[1] if len(sys.argv) > 1 else 'data/c539_hit.json'))
print(d.keys())
B=np.array(d['blocks']); print(B.shape)
q=11
names='ABCD'
def idx(x,y): return 7*x+y
orb={}
for x in range(7):
    for y in range(7):
        o=frozenset([(x,y),((2*x)%7,(2*y)%7),((4*x)%7,(4*y)%7)])
        orb[(x,y)]=o
orbs=sorted(set(orb.values()),key=lambda o:min(o))
print('orbits',len(orbs))
for X in range(4):
    s,t,e=B[X]
    print(names[X],'sum s^2+q t^2 =',None)
    pos=[(k//7,k%7) for k in range(49) if e[k]!=0]
    vals=[e[7*x+y] for (x,y) in pos]
    print(' e support',len(pos),'origin?',(0,0) in pos,'e values',sorted(set(vals)),'n+',vals.count(1),'n-',vals.count(-1))
    # orbits with defect
    for o in orbs:
        ev=[e[7*x+y] for (x,y) in o]
        if any(ev): print('   orbit',sorted(o)[:3],'size',len(o),'e',ev[0],' s',[s[7*x+y] for (x,y) in sorted(o)][:3],'t',[t[7*x+y] for (x,y) in sorted(o)][:3])
print('--- costs')
from collections import Counter
cc=Counter(); tot=0
for X in range(4):
    s,t,e=B[X]
    for k in range(49):
        if e[k]!=0:
            x0=s[k]+e[k]; c=(q-1)*e[k]**2+2*x0*e[k]; cc[(names[X],int(c))]+=1; tot+=c
print(cc, 'total',tot)
# K-level real identity check: level-1 sum w w* = 4 m q delta with w = s + q e ; all lags
import numpy as np
m=49
def idx2(k): return (k//7,k%7)
def paf2(w):
    W=w.reshape(7,7); F=np.fft.fft2(W); return np.fft.ifft2(np.abs(F)**2).real
tot2=np.zeros((7,7))
for X in range(4):
    s,t,e=B[X]; w=(s+q*e).astype(float); tot2+=paf2(w)
print('level-1 total PAF (should be 4mq=2156 at 0 and 0 elsewhere):', round(tot2[0,0],6), np.abs(tot2.flatten()[1:]).max())
tot3=np.zeros((7,7))
for X in range(4):
    s,t,e=B[X]; tot3+=paf2(s.astype(float))+q*paf2(t.astype(float))
print('(d) sum(ss*+q tt*):', round(tot3[0,0],6), np.abs(tot3.flatten()[1:]).max())
