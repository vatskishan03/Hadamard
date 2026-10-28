#!/usr/bin/env python3
"""Planted positive control of tower.c at real size (and an unplanted resolution of the same (R,w)-tuple for timing).

Takes the heaviest R-tuple of the pattern (18, 0, 14, 14), runs stage 1 with NS sampled C, D members, takes the (R,w)-tuple of
the given priority rank, draws one member from each of its four fibres (A, B from the complete lists, C and D by the
rejection sampler), and sets the targets (the 35 PAF columns and the line-1 and line-2 targets) to the sums of the four
members' values.  The tower must then find the planted quadruple in the right child and grandchild, and its output must
equal the planted members and pass the exact PAF check.  Both runs print the tower statistics (children, grandchildren,
leaves, covered quadruples, time).

    python3 test_tower.py [NS=3e4] [RANK=0] [SEED=99]
Needs ./pjoin and ./tower.
"""
import sys, os, time, tempfile
import numpy as np
from engine import *
from lines import paf_line

ns = int(float(sys.argv[1])) if len(sys.argv) > 1 else 30000
rank = int(sys.argv[2]) if len(sys.argv) > 2 else 0
seed = int(sys.argv[3]) if len(sys.argv) > 3 else 99
ctx = Ctx(); S, F = ctx.S, ctx.F
t0 = time.time()
pat = (18, 0, 14, 14)
L, tl, wts, uniq = ctx.rtuples(pat)
j0 = int(np.argmax(wts))
R = RTuple(ctx, pat, j0, L, tl)
wd = tempfile.mkdtemp(prefix="towertest_")
quads, cl = stage1(ctx, R, ns, 5, wd)
pr = tuple_priority(ctx, R, cl, quads); o = np.argsort(-pr)
quad = quads[o[rank]]
print("stage 1: %d (R,w)-tuples (%.0fs); (R,w)-tuple of rank %d: class indices %s" % (len(quads), time.time() - t0, rank, quad.tolist()), flush=True)
prepare_members(ctx, R)
a, b, c, d = [int(x) for x in quad]
spec, out = os.path.join(wd, "spec.bin"), os.path.join(wd, "out.txt")
# unplanted
t1 = time.time(); ia, ib = write_spec(ctx, R, spec, a, b, R.r[2], cl[0][c], R.r[3], cl[1][d])
hits, stats, l1 = run_tower(spec, out)
print("UNPLANTED:", l1, "|", stats, "| wall %.1fs" % (time.time() - t1), flush=True)
# planted
rng = np.random.default_rng(seed)
XA = R.XA[ia]; XB = R.XB[ib]
pa, pb = int(rng.integers(len(XA))), int(rng.integers(len(XB)))
XC, _ = F.sample_CD(R.r[2], cl[0][c], 1, rng, max_tries=5000); XD, _ = F.sample_CD(R.r[3], cl[1][d], 1, rng, max_tries=5000)
cells = [F.cells_A(XA[pa:pa + 1]), F.cells_A(XB[pb:pb + 1]), F.cells_C(XC), F.cells_C(XD)]
seqs = [S.expand(cc, k)[0] for cc, k in zip(cells, "AACC")]
Tp = sum(S.paf_classes(sq[None])[:, S.cols][0] for sq in seqs).astype(np.int32)
Xs = [XA[pa:pa + 1], XB[pb:pb + 1], XC, XD]
def labs(X, sym, row):
    l = np.einsum("ncz,c->nz", X.astype(np.int64), row)
    return l[:, ctx.symidx] if sym else l
tau1p = sum(paf_line(labs(Xs[q], q < 2, ctx.row1))[0] for q in range(4)); tau2p = sum(paf_line(labs(Xs[q], q < 2, ctx.row2))[0] for q in range(4))
t1 = time.time(); write_spec(ctx, R, spec, a, b, R.r[2], cl[0][c], R.r[3], cl[1][d], T=Tp, tau1=tau1p, tau2=tau2p)
hits, stats, l1 = run_tower(spec, out)
print("PLANTED:", l1, "|", stats, "| wall %.1fs" % (time.time() - t1))
found = False
for h in hits:
    sq = hit_to_sequences(ctx, dict(A=XA[h[0]].tolist(), B=XB[h[1]].tolist(), C=h[2], D=h[3]))
    exact = bool((sum(S.paf_direct(x) for x in sq) == sum(S.paf_direct(x) for x in seqs)).all())
    same = all(np.array_equal(sq[i], seqs[i]) for i in range(4))
    print(" hit: exact PAF check %s; equals the planted quadruple: %s" % (exact, same))
    found = found or (exact and same)
for f in os.listdir(wd): os.remove(os.path.join(wd, f))
os.rmdir(wd)
print("PLANTED CONTROL PASSED" if found else "PLANTED CONTROL FAILED")
sys.exit(0 if found else 1)
