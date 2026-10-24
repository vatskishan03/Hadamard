"""The sign vectors (sigma, sigma') of the twisted-skew blocks A (PROOFS.md P3: c_j + i conj(c_j) = sigma_j (1+i), c_j - i conj(c_j) = sigma'_j (1-i)) and two invariants, <sigma, sigma'> and the parity of the number of minus signs of sigma', counted for used and unused A (k = 3 mod 8; search/NOTES.md section 92).
usage: python3 astruct.py v k FILE        (FILE: SOL lines of ../qskew or of tri2)"""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from trilib import *
import itertools
v=int(sys.argv[1]); k=int(sys.argv[2]); fn=sys.argv[3]
d,sols=load_sols(fn,v,k); e=d['e']; h=e//2
cntA=collections.Counter(s[0] for s in sols)
def sig(b):
    # sigma_j = sign of c_j + i conj(c_j) = +-(1+i);  sigma'_j = sign of c_j - i conj(c_j) = +-(1-i)
    sg=[];sp=[]
    for j in range(h):
        c=I4[b[1][j]]
        z=c+1j*np.conj(c); sg.append(int(np.sign(z.real)))
        z2=c-1j*np.conj(c); sp.append(int(np.sign(z2.real)))
    return b[0],tuple(sg),tuple(sp)
rows=[]
for A,cnt in cntA.items(): pass
# enumerate all A of the list
allA=[]
for u0 in (0,3):
    m=1 if u0==0 else -1; npl=(h+m)//2
    for sgn in itertools.product((1,-1),repeat=h):
        if sum(1 for s in sgn if s>0)!=npl: continue
        for ch in itertools.product((0,1),repeat=h):
            c=[0]*e
            for i in range(h):
                ci=ch[i] if sgn[i]>0 else 2+ch[i]; c[i]=ci; c[i+h]=(1-ci)%4
            allA.append((u0,tuple(c)))
print('v',v,'A list',len(allA),'used',len(cntA))
# invariants: dot(sigma,sigma'), parity of #minus in sigma'
inv=collections.defaultdict(lambda:[0,0])
for b in allA:
    u0,sg,sp=sig(b)
    key=(u0,sum(a*bb for a,bb in zip(sg,sp)))
    inv[key][0 if b in cntA else 1]+=1
for kk in sorted(inv): print(' (u0, <sigma,sigma_prime>)',kk,'used/unused',inv[kk])
inv2=collections.defaultdict(lambda:[0,0])
for b in allA:
    u0,sg,sp=sig(b)
    key=(u0,sum(1 for x in sp if x<0)%2)
    inv2[key][0 if b in cntA else 1]+=1
for kk in sorted(inv2): print(' (u0, #minus in sigma_prime mod 2)',kk,'used/unused',inv2[kk])
