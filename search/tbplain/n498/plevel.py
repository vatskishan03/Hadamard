#!/usr/bin/env python3
"""P-level (Z123 level) objects of QH(498): balanced ternary slices and the class vectors of slice triples.

The Z2 coordinate splits a sequence x over G = Z2 x Z123 (Z123 = Z3 x Z41, elements (b, z)) into p = (x(0, .) + x(1, .))/2 and
m = (x(0, .) - x(1, .))/2, both constant on the orbits of H.  p takes values in {-1, 0, 1} and is the "P part"; m is a sign on the
complement of the support of p (the "M part", completed in band.c / mcomp.py).  On the nine Z41-orbits (z -> 10 z, orbit 0 = {0})
p is a triple of slices u_b in {-1, 0, 1}^9, b in Z3.  For the compressions of the production pattern, y_A = y_C = y_D = (2, 0, 0)
and y_B = 0 (row sums 2, 0, 2, 2), the fibre sums force u_b[0] = 1 for b = 0 of A, C, D and u_b[0] = 0 for every other slice, and
u_b[1..8] to sum to zero ("balanced": 1107 vectors).  So slices are indexed 0..1106 (centre 0) and 1107..2213 (centre 1), NU = 2214.
A- and B-type objects (x symmetric) are a symmetric slice for b = 0 (19 per centre) and a free slice y with u_2 = reflection of u_1
(object id s * 1107 + y); C- and D-type objects are three free slices (id (c0 * 1107 + c1) * 1107 + c2, c0 the centre-1 slice).
`pcv` is the 14-column P class vector (5 columns with b' = 0, 9 with b' = 1) of a slice triple, `PTARGET` the target of the sum of the
four.  Library module (verify_sol.py, mcomp.py); running the file self-tests the tables.

    python3 plevel.py
"""
import itertools, sys
import numpy as np
from structure import *

# ---- slice universe U: 9-vectors (center c in {0,1}, entries o=1..8 ternary balanced)
def _c_order():
    out = []
    for code in range(6561):
        c = code; v = []
        for i in range(8):
            v.append(c % 3 - 1); c //= 3
        if sum(v) == 0:
            out.append(v)
    return out
V8 = np.array(_c_order(), dtype=np.int64)     # same enumeration order as ptgen.c
assert len(V8) == 1107
NV = len(V8)
U = np.zeros((2 * NV, 9), dtype=np.int64)
U[:NV, 1:] = V8
U[NV:, 1:] = V8
U[NV:, 0] = 1
NU = 2 * NV
# reflection o -> NEGO[o]: (refl v)_o = v_{NEGO[o]}
RF = np.zeros(NU, dtype=np.int64)
key2idx = {tuple(u): i for i, u in enumerate(U)}
for i, u in enumerate(U):
    RF[i] = key2idx[tuple(u[NEGO])]
# G0[u][o'] = corr(u,u)[o'];  G1[u][w][o'] = corr(u,w)[o']
G1 = np.zeros((NU, NU, 9), dtype=np.int8)
for op in range(9):
    G1[:, :, op] = (U @ NN[op] @ U.T).astype(np.int8)
G0 = np.array([G1[i, i] for i in range(NU)])
SYMIDX = {c: [i for i in range(NU) if (U[i] == U[i][NEGO]).all() and U[i][0] == c] for c in (0, 1)}
assert len(SYMIDX[0]) == 19 and len(SYMIDX[1]) == 19
CLS0 = [0, 1, 2, 3, 6]          # (0,o') class reps: o'=0 (weight), pairs (1,4),(2,5),(3,7),(6,8)
PTARGET = np.array([243, -6, -6, -6, -6] + [0] * 9)

def pcv(u0, u1, u2):
    """P class vector (14) of slice triples given as index arrays (broadcast)"""
    g0 = G0[u0].astype(np.int32) + G0[u1] + G0[u2]
    g1 = G1[u0, u1].astype(np.int32) + G1[u1, u2] + G1[u2, u0]
    return np.concatenate([g0[..., CLS0], g1], axis=-1)

if __name__ == "__main__":
    # consistency: symmetry of the (0,o') columns
    rng = np.random.default_rng(1)
    u0 = rng.integers(NV, NU, 1000); u1 = rng.integers(0, NV, 1000); u2 = rng.integers(0, NV, 1000)
    g0 = G0[u0].astype(np.int32) + G0[u1] + G0[u2]
    for o in range(9):
        assert (g0[:, o] == g0[:, NEGO[o]]).all()
    print("P-level tables OK; NU", NU, "G1 bytes", G1.nbytes)
