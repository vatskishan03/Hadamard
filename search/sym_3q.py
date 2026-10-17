"""Symmetric analogue of skew_3q.py (skew/PALEY_SKEW.md): is there a placement of the T-matrices of order 3 with a
symmetric A_c and two of the other A_j with equal Gram matrices? Then the propus array (A_c, A_j, A_j, A_l) would be a
symmetric Hadamard matrix of order 3(q + 1). Exhaustive over the 192 placements and the base points of orbits 1-3.

    python3 search/sym_3q.py 7 11 19 23 31 43      # hits only for q = 11
"""
import sys, itertools
import numpy as np
sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
from skew_3q import F_table

def blocks(q, r):
    s = (q + 1) // 4
    F, z = F_table(q)
    Fx = lambda d: np.where((d // (q + 1)) % 2 == 1, -1, 1) * F[d % (q + 1)]
    y = np.arange(s)
    N = [[(-1) ** k * Fx(j - k + 4 * (y + r[j] - r[k])) + ((k == j) & (y == 0)) for j in range(4)] for k in range(4)]
    wrap = np.where(y[None, :] >= y[:, None], 1, -1)
    nega = lambda row: wrap * row[(y[None, :] - y[:, None]) % s]
    return [[nega(N[k][j]) for j in range(4)] for k in range(4)], s

def run(q):
    hits = []
    Z = np.roll(np.eye(3, dtype=np.int64), 1, axis=1)
    Zp = [np.eye(3, dtype=np.int64), Z, Z @ Z]
    for r in itertools.product(range(2), repeat=3):          # base-point shifts of orbits 1..3 (0 or -1)
        rr = [0] + [-x for x in r]
        Nm, s = blocks(q, rr)
        for zero in range(4):
            others = [k for k in range(4) if k != zero]
            for perm in itertools.permutations(range(3)):
                for signs in itertools.product((1, -1), repeat=3):
                    T = [np.zeros((3, 3), np.int64)] * 4
                    T = list(T)
                    for k, p, sg in zip(others, perm, signs):
                        T[k] = sg * Zp[p]
                    A = [sum(np.kron(T[k], Nm[k][j]) for k in range(4)) for j in range(4)]
                    sym = [j for j in range(4) if np.array_equal(A[j], A[j].T)]
                    if not sym:
                        continue
                    G = [A[j] @ A[j].T for j in range(4)]
                    for c in sym:
                        rest = [j for j in range(4) if j != c]
                        for a, b in itertools.combinations(rest, 2):
                            if np.array_equal(G[a], G[b]):
                                hits.append((tuple(rr), zero, perm, signs, c, a, b))
    return hits

for q in [int(x) for x in sys.argv[1:]]:
    h = run(q)
    print(q, len(h), h[:2], flush=True)
