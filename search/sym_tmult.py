#!/usr/bin/env python3
"""Symmetric multiplication of the propus array by an odd order t (search/NOTES.md section 49).

For symmetric pairwise commuting A, B, D of order q with A^2 + 2B^2 + D^2 = 4qI, put X_k = sum_j M_kj (x) Y_j
(k, j in {a, b, d}, Y = (A, B, D)) with M_kj symmetric (0, +-1) elements of the group ring of an abelian group G of
odd order t, the supports of M_k. partitioning G for each k. If

    sum_k c_k M_kj M_ki = t c_j delta_ij,   c = (1, 2, 1),

then X_a, X_b, X_d are symmetric, commute and X_a^2 + 2 X_b^2 + X_d^2 = 4qt I, so the propus array gives a symmetric
Hadamard matrix of order 4qt. At the trivial character the row sums give an integer matrix with
M^T diag(1, 2, 1) M = t diag(1, 2, 1), so det(M)^2 = t^3 and t must be a square. This program searches all such
M exhaustively (meet in the middle with 64-bit hashes, then an exact check):

    python3 search/sym_tmult.py 1      control: the signed permutations
    python3 search/sym_tmult.py 3      none (and 5: none; t not a square)
    python3 search/sym_tmult.py 9      none over Z_9
    python3 search/sym_tmult.py 3x3    none over Z_3 x Z_3
"""
import itertools, sys

import numpy as np


def search(dims, limit=20):
    t = int(np.prod(dims))
    els = list(itertools.product(*[range(d) for d in dims]))
    idx = {e: i for i, e in enumerate(els)}
    orbits, seen = [], set()
    for e in els:
        if e not in seen:
            o = sorted({e, tuple((-x) % d for x, d in zip(e, dims))})
            seen |= set(o)
            orbits.append([idx[x] for x in o])
    c = np.array([1, 2, 1])
    rows = list(itertools.product(itertools.product(range(3), (1, -1)), repeat=len(orbits)))
    F = np.zeros((len(rows), 3, t), dtype=np.int64)
    for r, row in enumerate(rows):
        for o, (j, s) in zip(orbits, row):
            F[r, j, o] = s
    axes = tuple(range(2, 2 + len(dims)))
    Fh = np.fft.fftn(F.reshape(len(rows), 3, *dims), axes=axes)
    G = np.einsum("rj...,ri...->rji...", Fh, Fh)
    G = np.rint(np.fft.ifftn(G, axes=tuple(range(3, 3 + len(dims)))).real).astype(np.int64).reshape(len(rows), -1)
    target = np.zeros((3, 3, t), dtype=np.int64)
    for j in range(3):
        target[j, j, 0] = t * c[j]
    target = target.reshape(-1)
    w = np.random.default_rng(1).integers(1, 2 ** 61, size=G.shape[1], dtype=np.int64)
    with np.errstate(over="ignore"):
        hG, hb, ht = G @ w, np.sort((2 * G) @ w), target @ w
        found = []
        for ia in range(len(rows)):                 # rows a and d have weight 1, row b weight 2
            rest = ht - hG[ia] - hG
            pos = np.searchsorted(hb, rest)
            pos[pos >= len(hb)] = 0
            for idd in np.flatnonzero(hb[pos] == rest):
                need = target - G[ia] - G[idd]
                for b in np.flatnonzero(((2 * G) == need).all(1)):
                    found.append((ia, int(idd), int(b)))
            if len(found) >= limit:
                break
    return len(orbits), len(rows), [[rows[k] for k in s] for s in found]


if __name__ == "__main__":
    dims = [int(x) for x in sys.argv[1].split("x")]
    n_orbits, n_rows, sols = search(dims)
    print(f"group {sys.argv[1]}: {n_orbits} orbits of x -> -x, {n_rows} rows per k, {len(sols)} solutions"
          + (" (first 20)" if len(sols) >= 20 else ""))
    for s in sols[:3]:
        print("  a:", s[0], " d:", s[1], " b:", s[2])
