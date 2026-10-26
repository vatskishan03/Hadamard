"""Control: the whole join machinery on a small analogue, by brute force.

usage: python3 small_e2e.py P [U] [--out DIR]
G = Z2^3 x Z_P, H = <(tau, U)> (U a unit of Z_P of order 1 or 7 modulo P), Q = Z2^3, T = the 7 nontrivial characters, N = |G| + 7,
core order 2N. Lists ALL H-invariant A-type (symmetric, zero at 0) and C-type sequences (exhaustively, 2^n sign vectors), forms
every row-sum pattern, joins the four lists on the PAF classes (hash of pair sums, every hit re-checked on all classes) and
prints the number of solutions per pattern; the solutions are saved to OUT/small_P_U.pkl. The same classes, targets, orbit
tables and class vectors (zs.ZS) as in the n = 955 run are used, so this tests them against exhaustive enumeration.
"""
import argparse, itertools, os, numpy as np, sys, subprocess, time
import tbio
from zq import *
_ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
_ap.add_argument("p", type=int); _ap.add_argument("u", type=int, nargs="?", default=1)
_ap.add_argument("--out", default=None, help="output directory (default: current directory)")
_a = _ap.parse_args()
_out = tbio.outdir(_a.out)
tau3 = mat([[0,0,1],[1,0,1],[0,1,0]])
Tn = [k for k in itertools.product(range(2),repeat=3) if any(k)]
p=_a.p; u=_a.u
H = block_diag(tau3, mat([[u]]))
st = ZS([2,2,2,p],[H],phi=lambda c:c[:,:3]%2,T=Tn,qmods=[2,2,2],name='small')
M0=st.m+7
print(st.info(), "M0",M0,"core order",2*M0)
pats=st.patterns(sumsq=4*M0-2)
print("patterns",pats)
def allvec(kind):
    orb,sz,free=st.kind_orbits(kind)
    n=len(free)
    S=np.array(list(itertools.product((1,-1),repeat=n)),dtype=np.int8)
    return S, st.seq(kind,S)
SA,XA=allvec("L"); SC,XC=allvec("A")
print("A list",len(XA),"C list",len(XC))
keep=[c for c in range(st.C) if c!=st.zero_class]
VA=st.paf_classes(XA)[:,keep]; VC=st.paf_classes(XC)[:,keep]
rA=XA.sum(1); rC=XC.sum(1)
T=st.tvec[keep]
found=[]
for pat in pats:
    a,b,c,d=pat
    ia=np.flatnonzero(rA==a); ib=np.flatnonzero(rA==b); ic=np.flatnonzero(rC==c); idd=np.flatnonzero(rC==d)
    print(pat,len(ia),len(ib),len(ic),len(idd))
    # hash join
    rng=np.random.default_rng(1); co=rng.integers(1,2**40,size=VA.shape[1]).astype(np.int64)
    hA=(VA[ia]@co); hB=(VA[ib]@co); hC=(VC[ic]@co); hD=(VC[idd]@co); hT=int(T@co)
    left={}
    for i,x in enumerate(hA):
        for j,y in enumerate(hC):
            left.setdefault(int(x+y),[]).append((i,j))
    cnt=0
    for k,x in enumerate(hB):
        for l,y in enumerate(hD):
            key=hT-int(x)-int(y)
            if key in left:
                for (i,j) in left[key]:
                    tot=VA[ia[i]]+VA[ib[k]]+VC[ic[j]]+VC[idd[l]]
                    if (tot==T).all():
                        cnt+=1; found.append((pat,XA[ia[i]],XA[ib[k]],XC[ic[j]],XC[idd[l]]))
    print(" solutions:",cnt)
import pickle; pickle.dump((st.mods,found),open(os.path.join(_out,f'small_{p}_{u}.pkl'),'wb'))
