"""Skew Goethals-Seidel structures over G = Z_p x Z_p (p prime) with odd-order automorphism groups
H <= GL(2, p): subgroups of a Singer cycle, of the split torus, and the unipotent/Borel groups. An
element of odd order never has eigenvalue -1, so every such H admits a skew first block.

  python3 skew_noncyclic.py rank p [MINO MAXO]
  python3 skew_noncyclic.py search p NAME [maxpool]      # NAME as printed by rank
"""
import math, sys, itertools
from math import gcd, log10
import numpy as np
from pathlib import Path as _P; sys.path[:0] = [str(_P(__file__).resolve().parent), str(_P(__file__).resolve().parent.parent / 'verify')]
import workdir as W
import spence_lib as L
import skew_rank as R


def mat_order(M, p):
    I = np.eye(2, dtype=np.int64); X = M % p; k = 1
    while not np.array_equal(X, I):
        X = X @ M % p; k += 1
        if k > p * p: return None
    return k


def odd_part(n):
    while n % 2 == 0: n //= 2
    return n


def singer(p):
    """2x2 matrix of multiplication by a generator of GF(p^2)^*."""
    c = next(x for x in range(2, p) if pow(x, (p - 1) // 2, p) == p - 1)       # non-residue
    for a in range(p):
        for b in range(1, p):
            M = np.array([[a, c * b], [b, a]], dtype=np.int64) % p
            if mat_order(M, p) == p * p - 1: return M
    raise ValueError


def mpow(M, e, p):
    R_ = np.eye(2, dtype=np.int64); B = M % p
    while e:
        if e & 1: R_ = R_ @ B % p
        B = B @ B % p; e >>= 1
    return R_


def candidates(p):
    out = []
    S = singer(p); N = p * p - 1
    for d in range(3, odd_part(N) + 1, 2):
        if odd_part(N) % d == 0: out.append((f"singer{d}", [mpow(S, N // d, p)]))
    g = next(x for x in range(2, p) if all(pow(x, (p - 1) // q, p) != 1 for q in range(2, p) if (p - 1) % q == 0 and all(q % r for r in range(2, q))))
    m = odd_part(p - 1)
    if m > 1:
        a = pow(g, (p - 1) // m, p)                      # order m
        for d1 in [d for d in range(1, m + 1) if m % d == 0]:
            for k in range(m):
                x = pow(a, m // d1, p) if d1 > 1 else 1
                y = pow(a, k, p)
                if d1 == 1 and k == 0: continue
                out.append((f"diag{d1}_{k}", [np.array([[x, 0], [0, y]], dtype=np.int64)]))
        out.append(("torus", [np.array([[a, 0], [0, 1]]), np.array([[1, 0], [0, a]])]))
    U = np.array([[1, 1], [0, 1]], dtype=np.int64)
    out.append(("unip", [U]))
    if m > 1:
        for d in [d for d in range(3, m + 1) if m % d == 0]:
            x = pow(a, m // d, p)
            for k in range(d):
                out.append((f"borel{d}_{k}", [U, np.array([[x, 0], [0, pow(x, k, p)]], dtype=np.int64)]))
    return out


def structure(p, gens):
    G = L.AbGroup([p, p])
    return G, L.Structure(G, [G.perm_from_matrix(M) for M in gens])


if __name__ == '__main__':
    p = int(sys.argv[2])
    if sys.argv[1] == 'rank':
        mino, maxo = (int(sys.argv[3]), int(sys.argv[4])) if len(sys.argv) > 4 else (6, 90)
        rng = np.random.default_rng(1)
        seen = set()
        for name, gens in candidates(p):
            G, st = structure(p, gens)
            key = tuple(sorted(np.bincount(st.orb).tolist())) + (st.C,)
            if not (mino <= st.O <= maxo): continue
            # rank() expects a cyclic group; reuse its body through a shim
            res = R.rank_structure(st, p * p, 4e6, rng=rng)
            for lE2, lE, q, lN in sorted(res, reverse=True)[:1]:
                print(f"{lE2:7.2f} p={p} {name} orbits={st.O} classes={st.C} rows=(1,{q[0]},{q[1]},{q[2]}) log10E={lE:.2f} "
                      f"log10N=({','.join(f'{x:.1f}' for x in lN)})", flush=True)
