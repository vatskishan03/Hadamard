#!/usr/bin/env python3
"""Level 0 of the tower (Q = Z2^2): compressions, pools and the exact join of the four blocks.

The compression y(a) = sum of a sequence over phi^-1(a) has the "spike" form y = (y0, y1, y1, y1) because H permutes the
three nonzero a.  `qdist` is the exact distribution of y over all sequences of one kind (dynamic program over the cells,
so the pool of every compression is a count, not an estimate); `tuples` joins the four blocks exactly on the Q-level PAF
equations and returns every ordered Q-tuple.  A Q-tuple is named by its row sums (|r_A|, |r_B|, |r_C|, |r_D|), r = y0 + 3 y1
(block signs are free: X -> -X maps y -> -y).  For structure "i" the pattern used for QH(358) is (18, 0, 14, 14).

    python3 qlevel.py [i|i3]          prints the patterns, their tuples and pool sizes
"""
import sys, json, math
from collections import defaultdict
import numpy as np
from structure import Struct

def cell_q(S, kind):
    orb = S.orbA if kind == "A" else S.orbC
    n = int(orb.max()) + 1
    W = np.zeros((n, 4), dtype=np.int64)
    np.add.at(W, (orb, S.q), 1)
    assert (W[:, 1] == W[:, 2]).all() and (W[:, 1] == W[:, 3]).all()
    return W[:, 0], W[:, 1]

def qdist(S, kind):
    w0, w1 = cell_q(S, kind)
    d = {(0, 0): 1}
    for a, b in zip(w0.tolist(), w1.tolist()):
        nd = defaultdict(int)
        for (y0, y1), n in d.items():
            nd[(y0 + a, y1 + b)] += n; nd[(y0 - a, y1 - b)] += n
        d = nd
    return dict(d)

def qtarget(S):
    qT = np.zeros(4, dtype=np.int64); np.add.at(qT, S.q, S.target)
    return int(qT[0]), int(qT[1])

def paf_spike(y0, y1):
    return (y0 * y0 + 3 * y1 * y1, 2 * y0 * y1 + 2 * y1 * y1)

def tuples(S):
    dA, dC = qdist(S, "A"), qdist(S, "C")
    tg = qtarget(S)
    def pairs(d):
        out = defaultdict(list)
        ks = list(d)
        for a in ks:
            pa = paf_spike(*a)
            for b in ks:
                pb = paf_spike(*b)
                out[(pa[0] + pb[0], pa[1] + pb[1])].append((a, b))
        return out
    AB, CD = pairs(dA), pairs(dC)
    sols = []
    for key, lst in AB.items():
        k2 = (tg[0] - key[0], tg[1] - key[1])
        for (a, b) in lst:
            for (c, d_) in CD.get(k2, ()):
                sols.append((a, b, c, d_))
    return sols, dA, dC, tg

rs = lambda y: y[0] + 3 * y[1]

if __name__ == "__main__":
    nm = sys.argv[1] if len(sys.argv) > 1 else "i"
    want = tuple(int(x) for x in sys.argv[2].split(",")) if len(sys.argv) > 2 else None
    S = Struct(nm)
    sols, dA, dC, tg = tuples(S)
    print(nm, "Q target", tg, "|A comps|", len(dA), "|C comps|", len(dC), "ordered Q tuples", len(sols), "totals 2^%d 2^%d" % (round(math.log2(sum(dA.values()))), round(math.log2(sum(dC.values())))))
    # normalised patterns: block negation (X -> -X) maps y -> -y ; use |row sums|
    pat = defaultdict(list)
    for t in sols:
        key = tuple(abs(rs(y)) for y in t)
        pat[key].append(t)
    print("ordered patterns (|rA|,|rB|,|rC|,|rD|):", len(pat))
    for k in sorted(pat):
        if want is not None and k != want: continue
        n = len(pat[k]); t = pat[k][0]
        print(k, "ordered tuples (incl. signs)", n, "e.g.", t, "pools", ["%.3g" % dd[y] for dd, y in zip((dA, dA, dC, dC), t)])
    # unordered-A/B, unordered-C/D classes and sign-normalised
    norm = set()
    for k in pat:
        a, b, c, d = k
        norm.add((tuple(sorted((a, b), reverse=True)), tuple(sorted((c, d), reverse=True))))
    print("patterns up to A<->B, C<->D and signs:", len(norm))
