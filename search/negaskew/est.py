"""Gaussian estimate of the number of quaternary skew solutions of an orbit shape (HEURISTIC; the Gaussian density at 0 of
the 2 (n-1) real PAF equations against the 2^(7 (n-1) + 1) normalised quadruples).  est(v, gens) returns
(bits of freedom, log2 of the equation probability, their sum = log2 of the expected number of normalised solutions).
The estimate is within about 3 bits of the exact counts when a row-sum pattern exists (calibration: v = 21, 7, 11) and
is useless for Paley shapes (NOTES.md sections 88 and 89).
usage: python3 est.py     the calibration shapes v = 21 <16>, 33, 7 <2>, 11 <3>, 191 <180> (needs qdrv.py)."""
import numpy as np, math, sys, json
from qdrv import shape
def est(v,gens,mode=1):
    d=shape(v,gens); n=d['n']; N=d['N']
    lg=0.0
    for c in range(1,n):
        if d['neg'][c]<c: continue  # pair (c, neg c): Re and Im equations
        M=N[c].astype(float)
        sre=sum((M[a,b]+M[b,a])**2/2 for a in range(n) for b in range(a+1,n))
        sim=sum((M[a,b]-M[b,a])**2/2 for a in range(n) for b in range(a+1,n))
        for s2 in (sre,sim):
            tot=4*s2
            if tot<=0: continue
            p=1.0/math.sqrt(2*math.pi*tot)
            lg+=math.log2(p)
    bits=7*(n-1)+1   # B,C,D: 2(n-1) each normalised; A: n-1 +1
    return bits,lg,bits+lg
if __name__=="__main__":
    for v,g in [(21,[16]),(33,[ [h for h in range(2,33) if h%3==1 and h%11==3][0] ]),(7,[2]),(11,[3]),(191,[180])]:
        print(v,est(v,g))
