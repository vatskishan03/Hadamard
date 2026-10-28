#!/usr/bin/env python3
"""Positive end-to-end control on the small analogue G = Z2 x Z3 (m = 6, N = 9, order 18) of the Q = Z3 structure of QH(498).

All cores by brute force (A, B symmetric, C, D free, sum_k PAF_k = 24 at 0, -12 on {b = 0} \\ {0}, 0 elsewhere), then for the
first five cores that have a border: the Q = Z3 border (../border.py), K by ../assemble.py (structure "18") and the exact check
K K* = 18 I.  The same border and assembly code produces the certificate of QH(498).

    python3 small_analogue.py [--out DIR]          with --out, writes DIR/complex-hadamard/C18.json for the first core
"""
import itertools, os, sys
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from assemble import assemble

out = sys.argv[sys.argv.index("--out") + 1] if "--out" in sys.argv else None
m = 6
coords = np.array(list(itertools.product(range(2), range(3))))
idx = lambda a, b: (a % 2) * 3 + (b % 3)
NEG = np.array([idx(-a, -b) for a, b in coords])
ADD = np.array([[idx(coords[g][0] + coords[h][0], coords[g][1] + coords[h][1]) for h in range(m)] for g in range(m)])
phi = np.array([b for a, b in coords])
def paf(X): return np.array([int((X * X[ADD[g]]).sum()) for g in range(m)])
target = np.array([24 if g == 0 else (-12 if phi[g] == 0 else 0) for g in range(m)])
sym = [np.array(s) for s in itertools.product((1, -1), repeat=6) if all(s[g] == s[NEG[g]] for g in range(m))]
allseq = [np.array(s) for s in itertools.product((1, -1), repeat=6)]
PA = {tuple(s): paf(s) for s in allseq}
sols = []
for A in sym:
    for B in sym:
        pab = PA[tuple(A)] + PA[tuple(B)]
        for C in allseq:
            pc = pab + PA[tuple(C)]
            for D in allseq:
                if (pc + PA[tuple(D)] == target).all():
                    sols.append((A, B, C, D))
print("cores found:", len(sols))
done = 0
for (A, B, C, D) in sols:
    seqs = ["".join("+" if v == 1 else "-" for v in X) for X in (A, B, C, D)]
    try:
        cert = assemble("18", seqs, out if done == 0 else None)
    except AssertionError as ex:
        if "no border" in str(ex): continue
        raise
    done += 1
    print("core %d: border found; K K* = 18 I exact; sha256 %s" % (done, cert["sha256"][:16]))
    if done >= 5: break
assert done >= 1
print("SMALL ANALOGUE PASSED")
