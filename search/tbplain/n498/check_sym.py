#!/usr/bin/env python3
"""Check that the canonical (symmetry-reduced) set of P-tuples of band.c meets every orbit of the full set of P-tuples.

The group of symmetries of the problem acting on P-tuples has order 32: the multipliers z -> 6^k z of Z41 acting on all four blocks
(k = 0..3: the quotient of Z41* by <H, -1> has order 4; -1 itself acts on A and B trivially), the independent reflections of C and
of D (a reflection does not change the PAF), and the swap C <-> D.  band.c keeps A canonical under the multipliers, B canonical under
the stabiliser of A, and of C, D the representatives with id(C) <= id(reflected C), id(D) <= id(reflected D), id(C) <= id(D).
Given the PTUPLE lines of two runs of the same band (DUMPT=1) -- FULL with NOCANA=1 NOCANB=1 NOREFL=1 NODEDUP=1 (no reduction) and
KEEP (reductions on) -- this script computes the orbit of every tuple of FULL under that group (k = 0..7 here, which contains it) and
checks that KEEP is a subset of FULL and that every orbit of FULL contains a tuple of KEEP, so that no solution class is lost.
It prints the number of orbits and the histogram of |orbit ∩ KEEP|.

    python3 check_sym.py FULL.sol KEEP.sol
"""
import sys, itertools, numpy as np
from structure import *
import plevel as PL
NV = PL.NV
U = PL.U
key2idx = {tuple(u): i for i, u in enumerate(U)}
g = 6
def sigma_map(k):
    gk = pow(g, k, 41)
    sg = [int(orb_of[(gk * OREP[o]) % 41]) for o in range(9)]
    m = np.zeros(PL.NU, dtype=int)
    for u in range(PL.NU):
        v2 = [0] * 9
        for o in range(9): v2[sg[o]] = U[u][o]
        m[u] = key2idx[tuple(v2)]
    return m
SIG = [sigma_map(k) for k in range(8)]
symlist = {1: PL.SYMIDX[1], 0: PL.SYMIDX[0]}
symidx = {c: {u: i for i, u in enumerate(symlist[c])} for c in (0, 1)}
def sig_sym(k, cen, oid):
    s, y = divmod(oid, NV)
    u0 = symlist[cen][s]
    return symidx[cen][SIG[k][u0]] * NV + SIG[k][y]
def sig_c(k, cid):
    c2 = cid % NV; c1 = (cid // NV) % NV; c0 = cid // (NV * NV)
    return ((SIG[k][NV + c0] - NV) * NV + SIG[k][c1]) * NV + SIG[k][c2]
def refl_c(cid):
    c2 = cid % NV; c1 = (cid // NV) % NV; c0 = cid // (NV * NV)
    r0 = PL.RF[NV + c0] - NV; r1 = PL.RF[c2]; r2 = PL.RF[c1]
    return int((r0 * NV + r1) * NV + r2)
def images(t):
    a, b, c, d = t
    out = set()
    for k in range(8):
        a2, b2, c2, d2 = sig_sym(k, 1, a), sig_sym(k, 0, b), sig_c(k, c), sig_c(k, d)
        for rc in (0, 1):
            for rd in (0, 1):
                c3 = refl_c(c2) if rc else c2; d3 = refl_c(d2) if rd else d2
                out.add((a2, b2, c3, d3)); out.add((a2, b2, d3, c3))
    return out
def load(fn):
    return {tuple(map(int, l.split()[1:])) for l in open(fn) if l.startswith("PTUPLE")}
F = load(sys.argv[1]); K = load(sys.argv[2])
print("full", len(F), "canonical", len(K))
print("K subset of F:", K <= F)
orb_of_t = {}
meet_all = True; counts = []
for t in F:
    im = images(t)
    inter = im & K
    counts.append(len(inter))
    if not inter: meet_all = False
print("every orbit of F meets K:", meet_all, " intersection size histogram:", {c: counts.count(c) for c in sorted(set(counts))})
# distinct orbits in F
orbits = set()
for t in F: orbits.add(frozenset(images(t)))
print("distinct orbits in F:", len(orbits), " orbits hit by K:", sum(1 for o in orbits if o & K))
