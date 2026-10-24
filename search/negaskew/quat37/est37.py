import json, sys
import numpy as np
from classical import is_pp
sh=json.load(open(sys.argv[1]))
def symbols(q):
    # (count, c) per orbit for a free block: type0: 4 symbols c=0 ; type1: (a,b): c=(q-1)+2*Re(i^(a-b)): a-b=0 -> q+1, 2 -> q-3, 1,3 -> q-1 ; type2: 4*4 symbols: c=(q-1)/2 + Re(x0 conj e): e in {(1+i)/2...}: Re(x e*) = +-1/2 -> c=(q-1)/2 +- 1/2 ... 
    cs={}
    cs[0]=4
    for a in range(4):
        for b in range(4):
            d=(a-b)&3; c=(q-1)+2*(1 if d==0 else -1 if d==2 else 0)
            cs[c]=cs.get(c,0)+1
    # type 2: x0=i^a, e=(UR+UI i)/2 for m: ER=1,EI=1 etc -> e = (+-1 +- i)/2 ; Re(x0 conj e) = (+-1)/2 or 0? x0 conj e real part = UR*ER+UI*EI over 2... they use c=(q-1)/2+UR[a]*ER[m]+UI[a]*EI[m] (ER,EI = +-1) => c = (q-1)/2 + {1,-1,0,0}...
    UR=[1,0,-1,0]; UI=[0,1,0,-1]; ER=[1,-1,-1,1]; EI=[1,1,-1,-1]
    for a in range(4):
        for m in range(4):
            c=(q-1)//2+UR[a]*ER[m]+UI[a]*EI[m]; cs[c]=cs.get(c,0)+1
    return cs
rows=[]
for s in sh:
    if not s['ident']: continue
    q=s['q']; cs=symbols(q); cap=4*s['K']
    poly=np.zeros(cap+1,dtype=np.float64); poly[0]=1
    # free block incl. orbit 0 (x0 fixed=0 normalised: type0 count 1 etc. -> approx ignore)
    for l in s['ells']:
        new=np.zeros_like(poly)
        for c,cnt in cs.items():
            w=l*c
            if w>cap: continue
            new[w:]+=cnt*poly[:cap+1-w]
        poly=new
    nL=poly.sum()
    rows.append((s['nu'],s['v'],s['K'],s['q'],s['gens'],s['ells'][:3],nL,bool(is_pp(8*s['v']-1))))
for r in sorted(rows,key=lambda r:(r[0],r[6])): print(r[0],r[1],'K',r[2],'q',r[3],r[4],r[5],'~|L|<=%.2e'%r[6],'PALEY8' if r[7] else '')
