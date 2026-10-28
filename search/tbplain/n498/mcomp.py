#!/usr/bin/env python3
"""M-completion of a P-tuple in numpy (reference implementation of the M part of band.c): variables, quadratic forms, lists.

Given the P part of one sequence (a slice triple p), the M part m is a sign on every position (b, orbit) where p vanishes, with
the reflection constraint for symmetric A and B (x(-g) = x(g): m(b, o) = m(-b, -o)); the free signs are the "variables" of the
object.  A class value of the PAF at the 13 M columns (b' = 0 with o' in {1, 2, 3, 6}, b' = 1 with o' = 0..8) is a quadratic form
in the signs: PAF_X(0, g') = 2 (PAF_p + PAF_m)(g'), PAF_X(1, g') = 2 (PAF_p - PAF_m)(g'), so the target for the M columns is 0 in
the sum over the four blocks.  `variables`, `qform`, `mlist` (all 2^n sign vectors with their column vectors), `mslices` (signs ->
slice triple m), `to_words` (p, m -> the six 9-bit row words of the sequence, X(0, b, o) = p + m, X(1, b, o) = p - m).
Running the file validates the quadratic forms against the exact PAF over G on random objects.

    python3 mcomp.py
"""
import numpy as np, itertools
from structure import *
import plevel as PL

# M columns: (b'=0, o' in [1,2,3,6]), (b'=1, o'=0..8)
MCOLS = [(0, 1), (0, 2), (0, 3), (0, 6)] + [(1, o) for o in range(9)]
NMC = 13
WQ = np.zeros((NMC, 27, 27), dtype=np.int64)    # WQ[c][P1][P2], P = 9b+o
for c, (bp, op) in enumerate(MCOLS):
    for b in range(3):
        for o in range(9):
            for o2 in range(9):
                WQ[c, 9 * b + o, 9 * ((b + bp) % 3) + o2] += NN[op, o, o2]

def pslices(typ, u0, u1, u2):
    return np.array([PL.U[u0], PL.U[u1], PL.U[u2]])      # 3 x 9

def variables(typ, p):
    """list of variables, each a list of positions P=9b+o; typ 'S' (symmetric A/B) or 'C'"""
    vs = []
    if typ == 'C':
        for b in range(3):
            for o in range(9):
                if p[b, o] == 0:
                    vs.append([9 * b + o])
    else:
        for o in range(9):
            q = int(NEGO[o])
            if o <= q and p[0, o] == 0:
                vs.append(sorted({9 * 0 + o, 9 * 0 + q}))
        for o in range(9):
            if p[1, o] == 0:
                vs.append([9 * 1 + o, 9 * 2 + int(NEGO[o])])
    return vs

def qform(vs):
    n = len(vs)
    Q = np.zeros((NMC, n, n), dtype=np.int64)
    for i, vi in enumerate(vs):
        for j, vj in enumerate(vs):
            Q[:, i, j] = WQ[:, vi][:, :, vj].sum(axis=(1, 2))
    return Q

def signs(n):
    m = np.arange(1 << n)
    return (1 - 2 * ((m[:, None] >> np.arange(n)[None, :]) & 1)).astype(np.int64)

def mlist(typ, p, vs=None, Q=None):
    vs = variables(typ, p) if vs is None else vs
    Q = qform(vs) if Q is None else Q
    n = len(vs)
    S = signs(n)
    V = np.einsum('mi,cij,mj->mc', S, Q, S)
    return vs, Q, S, V

def mslices(vs, s):
    m = np.zeros((3, 9), dtype=np.int64)
    for i, v in enumerate(vs):
        for P in v:
            m[P // 9, P % 9] = s[i]
    return m

def to_words(p, m):
    """p, m: 3x9 slices -> rows (6 words): X(0,b,o)=p+m, X(1,b,o)=p-m"""
    R = np.zeros(6, dtype=np.int64)
    for a in range(2):
        for b in range(3):
            x = p[b] + m[b] if a == 0 else p[b] - m[b]
            assert set(np.abs(x).tolist()) == {1}, (p[b], m[b])
            R[3 * a + b] = sum(1 << o for o in range(9) if x[o] == -1)
    return R

def moments(Q):
    """mean vector (diag part) and covariance of the class vector of uniform random signs"""
    n = Q.shape[1]
    d = np.einsum('cii->c', Q).astype(np.float64)
    K = Q + Q.transpose(0, 2, 1)
    iu = np.triu_indices(n, 1)
    Kf = K[:, iu[0], iu[1]].astype(np.float64)          # c x npairs
    return d, Kf @ Kf.T

if __name__ == "__main__":
    rng = np.random.default_rng(3)
    # validation against exact PAF: random C-type object, random m signs
    for typ in ('C', 'S'):
        for trial in range(4):
            if typ == 'C':
                u0 = rng.integers(PL.NV, PL.NU); u1 = rng.integers(0, PL.NV); u2 = rng.integers(0, PL.NV)
            else:
                u0 = PL.SYMIDX[1][rng.integers(19)]; u1 = rng.integers(0, PL.NV); u2 = PL.RF[u1]
            p = pslices(typ, u0, u1, u2)
            vs, Q, S, V = mlist(typ, p)
            k = rng.integers(len(S))
            m = mslices(vs, S[k])
            Rw = to_words(p, m)
            X = words_to_full(Rw)
            if typ == 'S':
                assert is_symmetric(X)
            pf = paf_full(X)
            cv = classvec(Rw[None, :])[0]
            # M columns: PAF_X(a=0,g') = 2(PAF_p + PAF_m)(g'), PAF_X(a=1,g') = 2(PAF_p - PAF_m)
            pp = PL.pcv(np.array(u0), np.array(u1), np.array(u2))            # 14 P-columns
            for c, (bp, op) in enumerate(MCOLS):
                g0 = gidx(0, bp, OREP[op]); g1 = gidx(1, bp, OREP[op])
                # P value for (bp, op): column index in P vector
                if bp == 0:
                    pcol = PL.CLS0.index(op)
                else:
                    pcol = 5 + op
                assert pf[g0] == 2 * (pp[pcol] + V[k, c]), (typ, c)
                assert pf[g1] == 2 * (pp[pcol] - V[k, c]), (typ, c)
            print(typ, "n vars", len(vs), "OK")
