import sys
sys.path.insert(0,'.')
from common import *
from collections import Counter
m=int(sys.argv[1]); q=int(sys.argv[2]); fn=sys.argv[3]
S=Shape(m,[1]); orb=np.array(S.orb)
sols=load_solutions(fn)
QR=set((k*k)%m for k in range(1,m))
def two_valued_on_classes(x, shift):
    # x: array over K (len m), test constancy on QR+shift and on NQR+shift (excluding 'shift' point itself)
    vals={0:set(),1:set()}
    for k in range(m):
        if k==shift: continue
        d=(k-shift)%m
        vals[0 if d in QR else 1].add(round(x[k].real*2)/2+1j*round(x[k].imag*2)/2)
    return len(vals[0])==1 and len(vals[1])==1
cnt=Counter(); tot=0
pafflat=Counter()
for sol in sols:
    tot+=1
    for X,b in enumerate(sol):
        x0=np.array(U)[np.array(b[0])[orb]]; x1=np.array(U)[np.array(b[1])[orb]]; x2=np.array(U)[np.array(b[2])[orb]]
        e=(x1+x2)/2; s=x0-e; t=(x1-x2)/2
        for nm,f in (('s',s),('t',t),('x0',x0)):
            for sh in range(m):
                if two_valued_on_classes(f,sh): cnt[(X,nm)]+=1; break
    # PAF value histograms for E blocks C,D (blocks with e=0 at origin?) -- use s
print(tot,'solutions;  number with block-function 2-valued on QR/NQR (some shift):')
for X in range(4): print('ABCD'[X],{nm:cnt[(X,nm)] for nm in ('s','t','x0')})
