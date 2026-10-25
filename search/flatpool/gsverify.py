#!/usr/bin/env python3
"""Exact verification of a Goethals-Seidel quadruple from orbit signs; builds the GS array and checks H H^T = 4v I."""
import numpy as np, sys, json
import gsflat

def seq_from_mask(S, f0, mask):
    m = S['m']
    s = [(-1 if (mask >> j) & 1 else 1) for j in range(m)]
    return gsflat.full_seq(S, f0, s), s

def circ(a):
    v = len(a)
    idx = (np.arange(v)[None, :] - np.arange(v)[:, None]) % v
    return a[idx]            # C[i,j] = a[(j-i) mod v]

def gs_matrix(A, B, C, D):
    v = len(A)
    a, b, c, d = circ(A), circ(B), circ(C), circ(D)
    R = np.eye(v, dtype=np.int64)[::-1]
    br, cr, dr = b @ R, c @ R, d @ R
    H = np.block([[a, br, cr, dr],
                  [-br, a, d.T @ R, -(c.T @ R)],
                  [-cr, -(d.T @ R), a, b.T @ R],
                  [-dr, c.T @ R, -(b.T @ R), a]])
    return H

def verify(seqs, check_matrix=True):
    v = len(seqs[0])
    tot = np.zeros(v, dtype=np.int64)
    for s in seqs:
        assert set(np.unique(s).tolist()) <= {-1, 1}
        tot += gsflat.paf_exact(s)
    paf_ok = bool(np.all(tot[1:] == 0)) and tot[0] == 4 * v
    res = dict(v=v, paf_sum_zero=paf_ok)
    if check_matrix:
        H = gs_matrix(*[np.asarray(s, dtype=np.int64) for s in seqs])
        Hf = H.astype(np.float32)
        G = Hf @ Hf.T
        n = 4 * v
        res['order'] = n
        res['HHt_equals_nI'] = bool(np.array_equal(G.astype(np.int64), n * np.eye(n, dtype=np.int64)))
        res['entries_pm1'] = bool(np.all(np.abs(H) == 1))
    return res
