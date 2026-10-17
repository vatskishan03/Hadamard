#!/usr/bin/env python3
"""Seberry-Yamada's Lemma 26 with a repeated block (U_3 = U_4, V_3 = V_4; Seberry-Balonin's Theorem 3) with every U_k,
V_k in the cyclotomic scheme of GF(p) of index e, -1 in the subgroup (search/NOTES.md section 49). The cells are the
identity and the e classes; each k in {1, 2, 3} puts every cell in U_k or in V_k, with a sign. The conditions, on the
eigenvalues j = 0 .. e (Gauss periods):
    row sums u_1(0) = 1, u_2(0) = u_3(0) = 0,
    u_1^2 + u_2^2 + 2 u_3^2 = 2p + 1 - 2p [j = 0],   v_1^2 + v_2^2 + 2 v_3^2 = 2p + 1.
A solution gives symmetric Williamson-type matrices of order v = 2p + 1 with W_3 = W_4, hence a symmetric Hadamard
matrix of order 4v (propus array). Meet in the middle over the assignments of k = 1, 2 against k = 3, then an exact
float check.

    python3 search/l26_scheme.py lo hi emax     every prime p = 1 (mod 4) in [lo, hi] and every even e <= emax
                                                dividing (p - 1)/2

Results: 5 400 8 finds only p = 5 (e = 2, v = 11) and p = 13 (e = 6, v = 27); 200 1250 4 finds nothing.
"""
import itertools, sys, cmath, math
import numpy as np
from sympy import isprime, primitive_root

def scan(p, e, out):
    g = primitive_root(p)
    f = (p - 1) // e
    cls = np.zeros(p, dtype=np.int64) - 1
    x = 1
    for i in range(p - 1):
        cls[x] = i % e
        x = x * g % p
    # Gauss periods eta_i = sum_{x in C_i} cos(2 pi x / p)
    eta = np.zeros(e)
    for xx in range(1, p):
        eta[cls[xx]] += math.cos(2 * math.pi * xx / p)
    # lam[c][j]: eigenvalue of cell c on eigenspace j (j = 0 trivial, j = 1..e: characters a in C_{j-1})
    lam = np.zeros((e + 1, e + 1))
    lam[0, :] = 1
    for i in range(e):
        lam[i + 1, 0] = f
        for j in range(e):
            lam[i + 1, j + 1] = eta[(i + j) % e]
    assigns = list(itertools.product((0, 1), (1, -1), repeat=e + 1))   # flattened pairs (UorV, sign) per cell
    U = np.zeros((len(assigns), e + 1)); Vv = np.zeros((len(assigns), e + 1))
    for a, asg in enumerate(assigns):
        for c in range(e + 1):
            side, s = asg[2 * c], asg[2 * c + 1]
            (U if side == 0 else Vv)[a] += s * lam[c]
    rs = np.rint(U[:, 0]).astype(int)
    tU = np.full(e + 1, 2 * p + 1.0); tU[0] = 1.0
    tV = np.full(e + 1, 2 * p + 1.0)
    F = np.hstack([U ** 2, Vv ** 2])
    T = np.concatenate([tU, tV])
    i1 = np.flatnonzero(rs == 1); i0 = np.flatnonzero(rs == 0)
    if len(i1) == 0 or len(i0) == 0:
        return 0
    key = lambda M: np.round(M * 1000).astype(np.int64)
    K3 = key(2 * F[i0])
    table = {}
    for idx, row in zip(i0, K3):
        table.setdefault(row.tobytes(), []).append(idx)
    found = 0
    for a1 in i1:
        rest = T - F[a1] - F[i0]                  # over a2
        R = key(rest)
        for a2, row in zip(i0, R):
            for a3 in table.get(row.tobytes(), []):
                # exact recheck in floats
                if np.allclose(F[a1] + F[a2] + 2 * F[a3], T, atol=1e-6):
                    found += 1
                    if found <= 3:
                        out.write(f"p={p} e={e} a1={assigns[a1]} a2={assigns[a2]} a3={assigns[a3]}\n")
    return found

if __name__ == "__main__":
    lo, hi, emax = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
    for p in range(lo, hi + 1):
        if not isprime(p) or p % 4 != 1:
            continue
        for e in range(2, emax + 1, 2):
            if ((p - 1) // 2) % e:
                continue
            n = scan(p, e, sys.stdout)
            if n:
                print(f"p={p} e={e}: {n} solutions  (v = {2 * p + 1})", flush=True)
