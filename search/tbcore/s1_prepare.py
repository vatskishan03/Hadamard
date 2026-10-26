"""Lists of one tuple of OUT/tuples.json: the exhaustive A and B lists and sampled C and D lists (library + small driver).

Library (used by s1_tstats.py, s1_prod.py): C_combos / sample_C (sample C- or D-type sequences with prescribed Z29-class
z and Q-compression q, uniformly over the sign choices consistent with the class), classvec (the PAF class values of
sequences, the zero class dropped), Tvec (the target on the same classes).

Driver: python3 s1_prepare.py TUPLE_INDEX LC SEED OUT_NPZ [--out DIR]
  builds the four lists of tuple TUPLE_INDEX (|C|, |D| = LC samples, seeded rng) and saves their class vectors and sequences
  to OUT_NPZ (relative to --out); s1_export.py turns such a file into the binary input of join4.
"""
import argparse, sys, time, itertools, os, numpy as np
from math import comb
import tbio
from s1_lists import *
Svals=[-7,-5,-3,-1,1,3,5,7]
bits_by_sum={}
allbits=np.array(list(itertools.product((1,-1),repeat=7)),dtype=np.int8)
for S in Svals:
    bits_by_sum[S]=allbits[allbits.sum(1)==S]
def C_combos(z, q):
    out=[]
    for c0 in (1,-1):
      for cE in (1,-1):
        if c0+7*cE!=z[0]: continue
        for ct in itertools.product((1,-1),repeat=4):
            S=tuple(z[t+1]-ct[t] for t in range(4))
            if any(abs(s)>7 or s%2==0 for s in S): continue
            if c0+7*sum(ct)!=q[0] or cE+sum(S)!=q[1]: continue
            w=1
            for s in S: w*=comb(7,(7+s)//2)
            out.append((c0,cE,ct,S,w))
    return out
def sample_C(z,q,L,rng):
    combos=C_combos(z,q)
    ws=np.array([c[4] for c in combos],dtype=float); tot=ws.sum()
    pick=rng.choice(len(combos),size=L,p=ws/tot)
    f0=np.zeros((L,p),dtype=np.int8); f1=np.zeros((L,p),dtype=np.int8)
    for ci,(c0,cE,ct,S,w) in enumerate(combos):
        sel=np.flatnonzero(pick==ci); n=len(sel)
        if n==0: continue
        f0[np.ix_(sel,[0])]=c0
        for t in range(4):
            for zz in O[t]: f0[sel,zz]=ct[t]
        f1[sel,0]=cE
        for t in range(4):
            bl=bits_by_sum[S[t]]
            ch=bl[rng.integers(0,len(bl),size=n)]
            for i,zz in enumerate(O[t]): f1[sel,zz]=ch[:,i]
    return f0,f1,tot
def classvec(X):
    V=st.paf_classes(X)
    keep=[c for c in range(st.C) if c!=st.zero_class]
    return V[:,keep].astype(np.int16)
keep=[c for c in range(st.C) if c!=st.zero_class]
Tvec=st.tvec[keep].astype(np.int32)
if __name__=="__main__":
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("ti", type=int); ap.add_argument("LC", type=int); ap.add_argument("seed", type=int); ap.add_argument("outp")
    ap.add_argument("--out", default=None, help="output directory (default: current directory)")
    args = ap.parse_args()
    out = tbio.outdir(args.out)
    ti, LC, seed, outp = args.ti, args.LC, args.seed, os.path.join(out, args.outp)
    tuples=tbio.load_tuples(out)
    zA,zB,zC,zD,w=tuples[ti]
    rng=np.random.default_rng(seed)
    t0=time.time()
    f0,f1=enum_A(); XA=build(f0,f1)
    z5,q0,q1=zq_compress(XA)
    def Alist(zv):
        m=(q0==0)&(q1==3)&(z5==np.array(zv)).all(1)
        return XA[m]
    LA=Alist(zA); LB=Alist(zB)
    f0c,f1c,totC=sample_C(zC,(-1,1),LC,rng); XC=build(f0c,f1c)
    f0d,f1d,totD=sample_C(zD,(-1,1),LC,rng); XD=build(f0d,f1d)
    assert check_inv(XC[:1000]) and check_inv(XD[:1000])
    print("tuple",ti,"weight %.3g"%w,"|A|=%d |B|=%d fibre|C|=%.3g fibre|D|=%.3g sampled %d"%(len(LA),len(LB),totC,totD,LC), "time %.1f"%(time.time()-t0)); sys.stdout.flush()
    VA,VB,VC,VD=[classvec(x) for x in (LA,LB,XC,XD)]
    # sanity: Q-compress of C samples
    zc,qc0,qc1=zq_compress(XC); assert (qc0==-1).all() and (qc1==1).all() and (zc==np.array(zC)).all()
    np.savez(outp,VA=VA,VB=VB,VC=VC,VD=VD,T=Tvec,XA=LA,XB=LB,XC=XC,XD=XD)
    print("saved",outp,"time %.1f"%(time.time()-t0))
