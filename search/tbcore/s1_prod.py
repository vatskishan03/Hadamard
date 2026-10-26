"""Step 3: the production sampler of the n = 955 core (exact join of exhaustive A, B lists with sampled C, D lists).

usage: python3 s1_prod.py LC MAXIT SEED LCH WALL [--out DIR] [--replay IT,A,B,C,D]

Reads OUT/tuples.json and OUT/trows.npy (s1_tuples.py, s1_tstats.py). One iteration: a tuple is drawn uniformly from the
'balanced' tuples with |A| = |B| = 882 whose C and D fibres hold at least 3 LC sequences; the A and B lists (882 sequences each) are
the exhaustive lists of that tuple, the C and D lists are LC uniform samples of the fibres; the four lists of PAF class values are
joined exactly by join4 (hash of (A, C) sums against (B, D) sums in 2^LCH key chunks, every hit re-checked on all classes; the C
program is compiled into OUT on first use). Every hit is verified (sum of the four PAF classes = target) and written to
OUT/SOLUTION_<iteration>_<a>_<b>_<c>_<d>.json; the run stops after the first iteration with a hit. All randomness
comes from numpy.random.default_rng(SEED): tuple choice and the samples, in that order, so a run is determined by
(LC, SEED) up to the NumPy version (PCG64 stream; numpy 2.5.3 here). WALL is the wall-clock budget in seconds
(checked between iterations). --replay IT,A,B,C,D draws the lists of iterations 0..IT without joining and takes the indices
A, B, C, D as the solution of iteration IT (re-creates a known solution in about a minute without the join).
Original production run: python3 s1_prod.py 400000 4 101 5 900  (iteration 0 = tuple 60, join 190 s, two hits).
"""
import argparse, sys, time, json, subprocess, resource, os, numpy as np
import tbio
from s1_prepare import *

ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
ap.add_argument("LC", type=int, help="size of the sampled C and D lists")
ap.add_argument("maxit", type=int, help="maximal number of iterations")
ap.add_argument("seed0", type=int, help="seed of numpy.random.default_rng")
ap.add_argument("lch", type=int, help="log2 of the number of key chunks of the join (memory)")
ap.add_argument("wall", type=float, help="wall-clock budget (s)")
ap.add_argument("--out", default=None, help="output directory (default: current directory)")
ap.add_argument("--replay", default=None, help="IT,A,B,C,D: skip the join, take these indices as the hit of iteration IT")
args = ap.parse_args()
out = tbio.outdir(args.out)
LC, maxit, seed0, lch, wall = args.LC, args.maxit, args.seed0, args.lch, args.wall
replay = tuple(int(x) for x in args.replay.split(",")) if args.replay else None

tuples=tbio.load_tuples(out)
rows=np.load(tbio.trows_path(out))
# choose tuples with |A|=|B|=882 (balanced small A,B) and both C fibres >= LC*3
ok=[i for i in range(len(tuples)) if rows[i,1]==882 and rows[i,2]==882 and rows[i,3]>=3*LC and rows[i,4]>=3*LC]
ok2=[i for i in range(len(tuples)) if rows[i,1]==1960 and rows[i,2]==1960 and rows[i,3]>=3*(LC*882//1960) and rows[i,4]>=3*(LC*882//1960)]
print("usable tuples: 882-type",len(ok),"1960-type",len(ok2)); sys.stdout.flush()
rng=np.random.default_rng(seed0)
f0,f1=enum_A(); XA=build(f0,f1); z5,q0,q1=zq_compress(XA)
mQ=(q0==0)&(q1==3)
def Alist(zv):
    m=mQ&(z5==np.array(zv)).all(1)
    return XA[m]
T0=time.time(); total_cmb=0.0
def lim():
    """child process: lowest priority and at most 3 GB of address space"""
    os.nice(19)
    soft, hard = resource.getrlimit(resource.RLIMIT_AS)
    cap = 3_000_000_000 if hard == resource.RLIM_INFINITY else min(3_000_000_000, hard)
    resource.setrlimit(resource.RLIMIT_AS,(cap,cap))
exe = None if replay else tbio.build_join4(out)
fn=os.path.join(out,'prod_tmp.bin')
for it in range(maxit):
    if time.time()-T0>wall: break
    if replay and it>replay[0]: break
    ti=int(rng.choice(ok))
    zA,zB,zC,zD,w=tuples[ti]
    LA=Alist(zA); LB=Alist(zB)
    f0c,f1c,totC=sample_C(zC,(-1,1),LC,rng); XC=build(f0c,f1c)
    f0d,f1d,totD=sample_C(zD,(-1,1),LC,rng); XD=build(f0d,f1d)
    if replay and it<replay[0]:
        print(f"it {it} tuple {ti} (replay: lists drawn, join skipped)"); continue
    if replay:
        sols=[["SOL"]+[str(x) for x in replay[1:]]]
        print(f"it {it} tuple {ti} |A|={len(LA)} |B|={len(LB)} L={LC} (replay) sols {len(sols)}"); sys.stdout.flush()
    else:
        VA,VB,VC,VD=[classvec(x) for x in (LA,LB,XC,XD)]
        tbio.write_binary_lists(fn,VA,VB,VC,VD,Tvec)
        t1=time.time()
        r=subprocess.run([exe,fn,str(lch)],capture_output=True,text=True,preexec_fn=lim)
        cmb=float(len(VA))*len(VB)*len(VC)*len(VD); total_cmb+=cmb
        sols=[l.split() for l in r.stdout.splitlines() if l.startswith('SOL')]
        print(f"it {it} tuple {ti} |A|={len(LA)} |B|={len(LB)} L={LC} combos {cmb:.2e} total {total_cmb:.2e} join {time.time()-t1:.0f}s sols {len(sols)} elapsed {time.time()-T0:.0f}s"); sys.stdout.flush()
        if r.returncode!=0: print("join error",r.returncode,r.stderr[-300:]); break
    if sols:
        for s in sols:
            a,b,c,d=map(int,s[1:5])
            Xs=[LA[a],LB[b],XC[c],XD[d]]
            # exact verification of all PAF classes
            tot=np.zeros(st.C,dtype=np.int64)
            for X in Xs: tot+=st.paf_classes(X[None,:])[0]
            print("verify: PAF sum classes equal target:",bool((tot==st.tvec).all()))
            sol={"tuple":ti,"A":Xs[0].tolist(),"B":Xs[1].tolist(),"C":Xs[2].tolist(),"D":Xs[3].tolist(),"ok":bool((tot==st.tvec).all())}
            with open(os.path.join(out,f"SOLUTION_{it}_{a}_{b}_{c}_{d}.json"),"w") as fo: json.dump(sol,fo)
        break
if os.path.exists(fn): os.remove(fn)
print("finished",time.time()-T0,"total combos %.3e"%total_cmb)
