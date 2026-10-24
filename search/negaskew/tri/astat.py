"""Which twisted-skew blocks A of the list of PROOFS.md P3 (k = 3 mod 8, A big) occur in solutions: the deviation score sum_m (|F_A(m)|^2 - (v-3))^2 and the extreme values of |F_A(m)|^2 over the e cosets, used against unused A, and the number of solutions per A by score decile (HEURISTIC statistics of EXACT solution sets; search/NOTES.md section 92).
usage: python3 astat.py v k FILE        (FILE: SOL lines of ../qskew or of tri2; k = 3 mod 8)"""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from trilib import *
from math import comb
v=int(sys.argv[1]); k=int(sys.argv[2]); fn=sys.argv[3]
d,sols=load_sols(fn,v,k); e=d['e']; h=e//2
cntA=collections.Counter(s[0] for s in sols)
# all A blocks of the list (big, k=3 mod 8)
allA=[]
import itertools
for u0 in (0,3):
    m=1 if u0==0 else -1; npl=(h+m)//2
    for signs in itertools.product((1,-1),repeat=h):
        if sum(1 for s in signs if s>0)!=npl: continue
        for ch in itertools.product((0,1),repeat=h):
            c=[0]*e
            for i in range(h):
                ci=ch[i] if signs[i]>0 else 2+ch[i]
                c[i]=ci; c[i+h]=(1-ci)%4
            allA.append((u0,tuple(c)))
print('A list size',len(allA))
def score(b):
    R,F=spectrum(d,b); P=np.abs(F)**2
    return float(((P-(v-3))**2).sum()), float(P.max()), float(P.min())
sc={b:score(b) for b in allA}
used=[sc[b][0] for b in cntA]; unused=[sc[b][0] for b in allA if b not in cntA]
print('used A: n',len(used),'mean dev2',np.mean(used),'min',min(used),'max',max(used))
print('unused A: n',len(unused),'mean dev2',np.mean(unused),'min',min(unused),'max',max(unused))
# solutions per A as function of score
rows=sorted(((sc[b][0],cntA.get(b,0)) for b in allA))
import math
bins=np.quantile([r[0] for r in rows],np.linspace(0,1,11))
for i in range(10):
    sel=[r for r in rows if bins[i]<=r[0]<=bins[i+1]]
    print('score decile',i,'range %.0f-%.0f'%(bins[i],bins[i+1]),'mean sols per A %.1f'%np.mean([r[1] for r in sel]),'frac used %.2f'%np.mean([r[1]>0 for r in sel]))
print('--- max/min of |F_A|^2 over m')
mx_u=sorted(set(round(sc[b][1],3) for b in cntA)); mx_n=sorted(set(round(sc[b][1],3) for b in allA if b not in cntA))
print('used: max|F|^2 values (min..max)',mx_u[:3],mx_u[-3:],' unused:',mx_n[:3],mx_n[-3:])
mn_u=sorted(set(round(sc[b][2],3) for b in cntA)); mn_n=sorted(set(round(sc[b][2],3) for b in allA if b not in cntA))
print('used: min|F|^2 values',mn_u[:3],mn_u[-3:],' unused:',mn_n[:3],mn_n[-3:])
