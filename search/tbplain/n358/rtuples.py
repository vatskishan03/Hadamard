#!/usr/bin/env python3
"""R-level tuples (exact, explicit) and the exact sampler of sequences with prescribed R-compression.

`explicit_tuples(RL, tup)` joins the four R-vector lists of a Q-tuple exactly (dictionary join on the PAF vectors over E)
and returns the lists and every R-tuple (i, j, k, l) with sum_k PAF_{r_k} = target on E, as indices into the lists
(`RLevel.vectors` order).  `RSampler` draws uniform sign vectors of a kind (A or C) with prescribed class values r (a
dynamic program per E-class, exact weights); stage 1 of the tower uses it for the C and D members.  Library module.
"""
from collections import defaultdict
import numpy as np

def explicit_tuples(RL, tup):
    kinds = ["A", "A", "C", "C"]
    L = [RL.vectors(k, y) for k, y in zip(kinds, tup)]
    groups = []
    for (R, W, P) in L:
        d = defaultdict(list)
        for i, p in enumerate(map(tuple, P.tolist())): d[p].append(i)
        groups.append(d)
    tE = tuple(int(x) for x in RL.tE)
    ab = defaultdict(list)                       # needed (C+D) PAF vector -> list of (pa, pb)
    for pa in groups[0]:
        for pb in groups[1]:
            ab[tuple(t - a - b for t, a, b in zip(tE, pa, pb))].append((pa, pb))
    out = []
    for pc, ic in groups[2].items():
        for pd, idd in groups[3].items():
            key = tuple(a + b for a, b in zip(pc, pd))
            if key in ab:
                for (pa, pb) in ab[key]:
                    for i in groups[0][pa]:
                        for j in groups[1][pb]:
                            for k in ic:
                                for l in idd:
                                    out.append((i, j, k, l))
    return L, out

class RSampler:
    """uniform sign vectors of a kind with prescribed class values r (per E-class independent DP)"""
    def __init__(self, S, RL):
        self.S, self.RL = S, RL
        self.d = {}
        for kind in "AC":
            orb = S.orbA if kind == "A" else S.orbC
            ncell = int(orb.max()) + 1; cellsize = np.bincount(orb)
            cl = RL.cls[kind]; clsize = np.bincount(cl)
            cellclass = np.array([cl[RL.piE[orb == c]][0] for c in range(ncell)])
            per = []
            for k in range(int(cl.max()) + 1):
                cells = np.flatnonzero(cellclass == k)
                contrib = [int(cellsize[c] // clsize[k]) for c in cells]
                f = [None] * (len(cells) + 1); f[len(cells)] = {0: 1.0}
                for i in range(len(cells) - 1, -1, -1):
                    dd = {}
                    for v, c in f[i + 1].items():
                        for s in (1, -1): dd[v + s * contrib[i]] = dd.get(v + s * contrib[i], 0.0) + c
                    f[i] = dd
                per.append((cells, contrib, f))
            self.d[kind] = (ncell, per)
    def sample(self, kind, r, n, rng):
        ncell, per = self.d[kind]
        s = np.zeros((n, ncell), dtype=np.int8)
        for k, (cells, contrib, f) in enumerate(per):
            rem = np.full(n, int(r[k]), dtype=np.int64)
            for i, c in enumerate(cells):
                nxt = f[i + 1]
                vals = np.unique(rem)
                lut_p = {int(v): nxt.get(int(v) - contrib[i], 0.0) for v in vals}; lut_m = {int(v): nxt.get(int(v) + contrib[i], 0.0) for v in vals}
                cp = np.array([lut_p[int(v)] for v in rem]); cm = np.array([lut_m[int(v)] for v in rem])
                up = rng.random(n) < cp / (cp + cm)
                s[:, c] = np.where(up, 1, -1); rem = rem - s[:, c] * contrib[i]
            assert (rem == 0).all()
        return s
