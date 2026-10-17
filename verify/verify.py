#!/usr/bin/env python3
"""Rebuild and check the Hadamard matrices in this repository.

Every order directory holds a certificate.json. For each certificate this script
  1. expands the construction data,
  2. checks the defining identities of the ingredients (complementary sequences,
     T-matrices, Williamson matrices, conference matrices, skew inputs),
  3. assembles the matrix H and checks H H^T = nI, and
  4. compares the SHA-256 of the canonical CSV encoding with the recorded value.

The canonical encoding is one row per line, entries 1 / -1 separated by commas, each
line terminated by "\\n" (numpy.savetxt with fmt="%d", delimiter=",").

Products are formed in float32, one block of rows at a time. The matrices multiplied
have entries 0 and +-1 and fewer than 2^24 columns, so every partial sum is an integer
of absolute value below 2^24 and the computation is exact.

The matrices of 14708, 14932 and of the orders above 20000 (20084, 21748, 24124, 24308, 24412, 27604, 29404,
31316, 31412, 35636, 38468, 39260, 51580) are not stored. Rebuilding and checking one of them takes from 30 seconds
(six cores) to 40 minutes (51580, about 3.5 GB of memory); the orders 24308 to 39260 from near-Williamson matrices
are assembled in int8 (0.6 to 1.5 GB).

Yamada certificates accept a prime q or a prime power q = p^k (fields "p" and "poly" give GF(q)).

    python3 verify/verify.py              check every order in the repository
    python3 verify/verify.py 2092 4492    check the named orders (hadamard/2092, hadamard/4492)
    python3 verify/verify.py skew/820     check one skew (or symmetric) order
    python3 verify/verify.py --write      also write H<n>.csv next to each certificate
"""
import argparse
import hashlib
import json
import itertools
import math
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
SIGN = {"+": 1, "-": -1, "0": 0}


def circulant(row):
    """Circulant matrix with first row `row`: A[i, j] = row[(j - i) mod n]."""
    row = np.asarray(row, dtype=np.int64)
    n = len(row)
    idx = np.arange(n)
    return row[(idx[None, :] - idx[:, None]) % n]


def circulant8(row):
    """The circulant of a +-1 row as an int8 matrix, without int64 index arrays (for large orders)."""
    row = np.asarray(row).astype(np.int8)
    n = len(row)
    w = np.concatenate([row, row])
    view = np.lib.stride_tricks.as_strided(w[n:], shape=(n, n), strides=(-w.strides[0], w.strides[0]))
    return np.ascontiguousarray(view)


def back_identity(n):
    R = np.zeros((n, n), dtype=np.int64)
    R[np.arange(n), n - 1 - np.arange(n)] = 1
    return R


def goethals_seidel(A, B, C, D):
    rev = lambda X: X[:, ::-1]                              # X R, with R = back_identity(n)
    return np.block([[A, rev(B), rev(C), rev(D)],
                     [-rev(B), A, rev(D.T), -rev(C.T)],
                     [-rev(C), -rev(D.T), A, rev(B.T)],
                     [-rev(D), rev(C.T), -rev(B.T), A]])


def abelian_group(spec):
    """G = Z_{n1}^{k1} x ... from a spec such as "3^2x11": (order, addition table, negation), elements indexed in
    mixed radix over the coordinates (the convention of search/gsgroup.py)."""
    mods = []
    for part in spec.split('x'):
        n, k = (part.split('^') + ['1'])[:2]
        mods += [int(n)] * int(k)
    elts = np.array(list(itertools.product(*[range(m) for m in mods])), dtype=np.int64).reshape(-1, len(mods))
    mods_a = np.array(mods, dtype=np.int64)
    radix = np.array([int(np.prod(mods[i + 1:])) for i in range(len(mods))], dtype=np.int64)
    add = ((elts[:, None, :] + elts[None, :, :]) % mods_a) @ radix
    neg = ((-elts) % mods_a) @ radix
    return len(elts), add, neg


def group_matrix(add, neg, x):
    """X[g, h] = x(h - g)"""
    return np.asarray(x, dtype=np.int64)[add[:, neg].T]


def goethals_seidel_group(add, neg, A, B, C, D):
    """the Goethals-Seidel array over an abelian group, with R[g, h] = [g + h = 0]"""
    t = len(neg)
    R = np.zeros((t, t), dtype=np.int64)
    R[np.arange(t), neg] = 1
    return np.block([[A, B @ R, C @ R, D @ R],
                     [-B @ R, A, D.T @ R, -C.T @ R],
                     [-C @ R, -D.T @ R, A, B.T @ R],
                     [-D @ R, C.T @ R, -B.T @ R, A]])


def paf_sum(rows):
    """Sum over the given sequences of their periodic autocorrelation functions
    (exact integer arithmetic)."""
    n = len(rows[0])
    total = np.zeros(n, dtype=np.int64)
    for x in rows:
        x = np.asarray(x, dtype=np.int64)
        total += np.array([x @ np.roll(x, -s) for s in range(n)], dtype=np.int64)
    return total


def is_pm1(X, block=2048):
    """True if every entry of X is 1 or -1 (checked a block of rows at a time)."""
    return all((np.abs(X[i:i + block]) == 1).all() for i in range(0, X.shape[0], block))


def gram_is_scalar(X, c, block=2048):
    """True if X X^T = cI, for X with entries 0, +-1. The product is formed in float32
    one pair of row blocks at a time (exact, see the module docstring); only blocks on
    and above the diagonal are needed."""
    n, k = X.shape
    if k >= 2 ** 24:
        raise ValueError("too many columns for an exact float32 product")
    for i in range(0, n, block):
        A = X[i:i + block].astype(np.float32)
        for j in range(i, n, block):
            G = A @ (A if j == i else X[j:j + block].astype(np.float32)).T
            if j == i:
                G[np.diag_indices_from(G)] -= c
            if G.any():
                return False
    return True


def build_goethals_seidel(cert):
    v = cert["v"]
    gens = cert["multipliers"] if "multipliers" in cert else [cert["multiplier"]]
    H = {1}
    while True:  # subgroup of Z_v^* generated by the multipliers
        grown = H | {h * g % v for h in H for g in gens}
        if grown == H:
            break
        H = grown
    H = sorted(H)
    if len(H) != cert["subgroup_order"]:
        raise ValueError("multipliers do not generate a subgroup of the stated order")
    seqs = []
    for reps in cert["blocks"]:
        x = np.ones(v, dtype=np.int64)
        for r in reps:
            x[[r * h % v for h in H]] = -1
        seqs.append(x)
    if [int(x.sum()) for x in seqs] != cert["rowsums"]:
        raise ValueError("row sums differ from the certificate")
    mu = cert.get("spin_multiplier")
    if mu is not None:
        t = np.arange(v)
        if not all(np.array_equal(seqs[k + 1], seqs[1][t * pow(mu, k, v) % v]) for k in range(3)):
            raise ValueError("the last three blocks are not multiplier images of each other")
    paf = paf_sum(seqs)
    if paf[0] != 4 * v or np.any(paf[1:] != 0):
        raise ValueError("sequences are not complementary")
    return goethals_seidel(*map(circulant, seqs)), seqs


def build_propus(cert):
    """Symmetric Hadamard matrix of order 4v from circulant +-1 matrices A, B, D of order v with
    A symmetric and A A^T + 2 B B^T + D D^T = 4v I (a propus-type quadruple (A, B, B, D)):

        [ A    B    B    D  ]
        [ B^T  DR  -A   -BR ]
        [ B^T -A   -DR   BR ]
        [ D^T -BR   BR  -A  ]

    (R the back-diagonal identity; XR is symmetric for circulant X and A is symmetric, so H = H^T.) The three
    sequences are given like Goethals-Seidel blocks: constant on the orbits of the multiplier
    subgroup, by one representative of each orbit on which the sequence is -1."""
    v = cert["v"]
    gens = cert["multipliers"]
    H = {1}
    while True:
        grown = H | {h * g % v for h in H for g in gens}
        if grown == H:
            break
        H = grown
    H = sorted(H)
    if len(H) != cert["subgroup_order"]:
        raise ValueError("multipliers do not generate a subgroup of the stated order")
    seqs = []
    for reps in cert["blocks"]:
        x = np.ones(v, dtype=np.int64)
        for r in reps:
            x[[r * h % v for h in H]] = -1
        seqs.append(x)
    if [int(x.sum()) for x in seqs] != cert["rowsums"]:
        raise ValueError("row sums differ from the certificate")
    A, B, D = seqs
    if not np.array_equal(A, A[(-np.arange(v)) % v]):
        raise ValueError("the first sequence is not symmetric")
    paf = paf_sum([A, B, D]) + paf_sum([B])
    if paf[0] != 4 * v or np.any(paf[1:] != 0):
        raise ValueError("PAF_A + 2 PAF_B + PAF_D is not 4v at 0 and 0 elsewhere")
    cA, cB, cD = map(circulant, seqs)
    R = back_identity(v)
    BR, DR = cB @ R, cD @ R
    return np.block([[cA, cB, cB, cD], [cB.T, DR, -cA, -BR], [cB.T, -cA, -DR, BR], [cD.T, -BR, BR, -cA]]), seqs


def build_cooper_wallis(cert):
    """T-matrices of order t (Cooper-Wallis) give OD(4t; t, t, t, t); substituting
    Williamson matrices of order w for its four variables gives H of order 4tw."""
    T = [np.array([SIGN[c] for c in s], dtype=np.int64) for s in cert["t_matrices"]]
    t = len(T[0])
    if not np.all(np.count_nonzero(np.array(T), axis=0) == 1):
        raise ValueError("T-matrix supports do not partition Z_t")
    paf = paf_sum(T)
    if paf[0] != t or np.any(paf[1:] != 0):
        raise ValueError("T-matrix condition fails")
    W = [circulant([SIGN[c] for c in s]) for s in cert["williamson"]]
    w = W[0].shape[0]
    if not all(np.array_equal(m, m.T) for m in W) or \
            not np.array_equal(sum(m @ m for m in W), 4 * w * np.eye(w, dtype=np.int64)):
        raise ValueError("Williamson condition fails")
    X = [circulant(x) for x in T]
    coeff = {"A": (X[0], X[1], X[2], X[3]), "B": (X[1], -X[0], -X[3], X[2]),
             "C": (X[2], X[3], -X[0], -X[1]), "D": (X[3], -X[2], X[1], -X[0])}
    H = 0
    for k in range(4):
        H = H + np.kron(goethals_seidel(*(coeff[c][k] for c in "ABCD")), W[k])
    return H, None


def build_cooper_wallis_near_williamson(cert):
    """Cooper-Wallis with Williamson-type matrices that are not Williamson matrices.

    As in `build_cooper_wallis`, T-matrices of order t give an orthogonal design OD(4t; t, t, t, t) and
    four +-1 matrices W_1, ..., W_4 of order w with W_k W_l^T = W_l W_k^T (k != l) and
    W_1 W_1^T + ... + W_4 W_4^T = 4w I substituted for its variables give a Hadamard matrix of order 4tw.
    The Williamson-type matrices come from near-Williamson matrices (Kharaghani, Mohammadian and
    Tayfeh-Rezaie, arXiv:2605.08661): circulant +-1 matrices A, B, C, D of order w, with B, C, D symmetric
    and A arbitrary, such that A A^T + B B^T + C C^T + D D^T = 4w I. Then (A R, B, C, D), with R the
    back-identity, are Williamson-type matrices: A R is symmetric and commutes with the symmetric circulants.
    The certificate gives the first rows of A, B, C, D under "near_williamson". Unlike `build_cooper_wallis`
    the matrix is assembled block by block in int8 (a Kronecker product in int64 would need 8 bytes per
    entry, 12 GB at order 38468)."""
    T = [np.array([SIGN[c] for c in s], dtype=np.int64) for s in cert["t_matrices"]]
    t = len(T[0])
    if len(T) != 4 or not np.all(np.count_nonzero(np.array(T), axis=0) == 1):
        raise ValueError("T-matrix supports do not partition Z_t")
    paf = paf_sum(T)
    if paf[0] != t or np.any(paf[1:] != 0):
        raise ValueError("T-matrix condition fails")
    seq = [np.array([SIGN[c] for c in s], dtype=np.int64) for s in cert["near_williamson"]]
    w = len(seq[0])
    if len(seq) != 4 or any(len(s) != w or not np.all(np.abs(s) == 1) for s in seq):
        raise ValueError("near-Williamson sequences must be four +-1 sequences of equal length")
    if any(not np.array_equal(s, s[-np.arange(w) % w]) for s in seq[1:]):
        raise ValueError("B, C, D must be symmetric circulants")
    paf = paf_sum(seq)
    if paf[0] != 4 * w or np.any(paf[1:] != 0):
        raise ValueError("near-Williamson condition fails")
    R = back_identity(w)
    W = [circulant(seq[0]) @ R] + [circulant(s) for s in seq[1:]]
    if not all(np.array_equal(W[k] @ W[l].T, W[l] @ W[k].T) for k in range(4) for l in range(k + 1, 4)) or \
            not np.array_equal(sum(m @ m.T for m in W), 4 * w * np.eye(w, dtype=np.int64)):
        raise ValueError("Williamson-type condition fails")
    X = [circulant(x) for x in T]
    coeff = {"A": (X[0], X[1], X[2], X[3]), "B": (X[1], -X[0], -X[3], X[2]),
             "C": (X[2], X[3], -X[0], -X[1]), "D": (X[3], -X[2], X[1], -X[0])}
    E = [goethals_seidel(*(coeff[c][k] for c in "ABCD")) for k in range(4)]
    I4t = np.eye(4 * t, dtype=np.int64)
    if not all(np.array_equal(E[k] @ E[k].T, t * I4t) for k in range(4)) or \
            any((E[k] @ E[l].T + E[l] @ E[k].T).any() for k in range(4) for l in range(k + 1, 4)):
        raise ValueError("T-matrices do not give an OD(4t; t, t, t, t)")
    if not np.all(sum(np.abs(e) for e in E) == 1):
        raise ValueError("an entry of the orthogonal design carries more or less than one variable")
    which, sign = np.argmax(np.stack([np.abs(e) for e in E]), axis=0), sum(E).astype(np.int8)
    W8 = np.stack(W).astype(np.int8)
    N = 4 * t * w
    H = np.empty((N, N), dtype=np.int8)
    for i in range(4 * t):                                   # block row i: blocks sign[i, j] * W[which[i, j]]
        H[i * w:(i + 1) * w] = (W8[which[i]] * sign[i][:, None, None]).transpose(1, 0, 2).reshape(w, N)
    return H, None


def bordered_williamson_type(seq):
    """Williamson-type matrices of order N = m + 1 from four +-1 sequences x_1, ..., x_4 on Z_m (search/bnw.py,
    williamson-type/README.md): x_2, x_3, x_4 symmetric and PAF_1 + ... + PAF_4 = -4 off 0. Then the row sums are
    +-2 at one index j and 0 elsewhere; with M_k the circulant of x_k (M_1 replaced by M_1 R) and s = r_j / 2,
    W_k = [[1, b_k e^T], [c_k e, M_k]] with c_k = 1, b_k = s for k != j and c_j = -1, b_j = -s. Checked here:
    W_k W_l^T = W_l W_k^T and sum W_k W_k^T = 4N I."""
    m = len(seq[0])
    if len(seq) != 4 or any(len(s) != m or not np.all(np.abs(s) == 1) for s in seq):
        raise ValueError("expected four +-1 sequences of equal length")
    if any(not np.array_equal(s, s[-np.arange(m) % m]) for s in seq[1:]):
        raise ValueError("x_2, x_3, x_4 must be symmetric")
    paf = paf_sum(seq)
    if paf[0] != 4 * m or np.any(paf[1:] != -4):
        raise ValueError("the PAFs do not sum to -4 off 0")
    r = [int(s.sum()) for s in seq]
    if sorted(abs(x) for x in r) != [0, 0, 0, 2]:
        raise ValueError(f"row sums {r}")
    j = [abs(x) for x in r].index(2)
    s, N = r[j] // 2, m + 1
    Ms = [circulant(seq[0]) @ back_identity(m)] + [circulant(x) for x in seq[1:]]
    W = []
    for k in range(4):
        Wk = np.empty((N, N), dtype=np.int64)
        Wk[0, 0] = 1
        Wk[0, 1:] = -s if k == j else s
        Wk[1:, 0] = -1 if k == j else 1
        Wk[1:, 1:] = Ms[k]
        W.append(Wk)
    if not all(np.array_equal(W[k] @ W[l].T, W[l] @ W[k].T) for k in range(4) for l in range(k + 1, 4)) or \
            not np.array_equal(sum(x @ x.T for x in W), 4 * N * np.eye(N, dtype=np.int64)):
        raise ValueError("Williamson-type condition fails")
    return W


def assemble_cooper_wallis(T, W):
    """OD(4t; t, t, t, t) from T-matrices T_1..T_4 of order t (checked), with the four Williamson-type matrices W
    substituted for its variables; assembled block by block in int8, as in build_cooper_wallis_near_williamson."""
    t = len(T[0])
    if len(T) != 4 or not np.all(np.count_nonzero(np.array(T), axis=0) == 1):
        raise ValueError("T-matrix supports do not partition Z_t")
    paf = paf_sum(T)
    if paf[0] != t or np.any(paf[1:] != 0):
        raise ValueError("T-matrix condition fails")
    X = [circulant(x) for x in T]
    coeff = {"A": (X[0], X[1], X[2], X[3]), "B": (X[1], -X[0], -X[3], X[2]),
             "C": (X[2], X[3], -X[0], -X[1]), "D": (X[3], -X[2], X[1], -X[0])}
    E = [goethals_seidel(*(coeff[c][k] for c in "ABCD")) for k in range(4)]
    I4t = np.eye(4 * t, dtype=np.int64)
    if not all(np.array_equal(E[k] @ E[k].T, t * I4t) for k in range(4)) or \
            any((E[k] @ E[l].T + E[l] @ E[k].T).any() for k in range(4) for l in range(k + 1, 4)):
        raise ValueError("T-matrices do not give an OD(4t; t, t, t, t)")
    if not np.all(sum(np.abs(e) for e in E) == 1):
        raise ValueError("an entry of the orthogonal design carries more or less than one variable")
    which, sign = np.argmax(np.stack([np.abs(e) for e in E]), axis=0), sum(E).astype(np.int8)
    W8 = np.stack(W).astype(np.int8)
    w = W8.shape[1]
    N = 4 * t * w
    H = np.empty((N, N), dtype=np.int8)
    for i in range(4 * t):
        H[i * w:(i + 1) * w] = (W8[which[i]] * sign[i][:, None, None]).transpose(1, 0, 2).reshape(w, N)
    return H


def build_miyamoto_blocks(cert):
    """Miyamoto's Hadamard matrix of order 4q (`miyamoto`) as sixteen bordered blocks, times T-matrices of order t:
    a Hadamard matrix of order 4qt (search/miyamoto_blocks.py, NOTES.md section 23).

    q = 1 (mod 4) is a prime and g ("g") a primitive root mod q; m = (q - 1)/2. Paley's conference matrix is indexed by
    infinity, 0, g^(2j), g^(2j+1) (j < m). K = [[A, B], [-B^T, A^T]] with A, B the negacyclic matrices of the pair
    "ng" = (a, b), checked here to be a negaperiodic Golay pair of length m (NAF_a + NAF_b = 2m at 0, 0 elsewhere).
    If the certificate has "pg" instead, (a, b) must be a periodic Golay pair (PAF_a + PAF_b = 2m delta) and A, B are
    circulant. With N_kj the blocks of order q of Miyamoto's matrix, T_k the circulant T-matrices ("T", checked),
    A_j = sum_k T_k (x) N_kj, and R = R_t (x) R_B, where R_t is x -> -x on Z_t and R_B fixes the border and sends
    (b, j) to (b, -j) if j = 0 or the pair is periodic, to (1 - b, -j) otherwise, the matrix is the Goethals-Seidel array
    [[A1, A2 R, A3 R, A4 R], [-A2 R, A1, A4^T R, -A3^T R], [-A3 R, -A4^T R, A1, A2^T R], [-A4 R, A3^T R, -A2^T R, A1]].
    Its H H^T = 4qt I is checked densely by build_checked."""
    q, g = cert["q"], cert["g"]
    if q % 4 != 1 or q < 5 or any(q % d == 0 for d in range(2, math.isqrt(q) + 1)):
        raise ValueError("q must be a prime = 1 (mod 4)")
    if g % q == 0 or any(pow(g, (q - 1) // p, q) == 1 for p in range(2, q) if (q - 1) % p == 0 and
                         all(p % d for d in range(2, math.isqrt(p) + 1))):
        raise ValueError("g is not a primitive root mod q")
    m = (q - 1) // 2
    periodic = "pg" in cert
    a, b = (np.array([SIGN[c] for c in s], dtype=np.int64) for s in cert["pg" if periodic else "ng"])
    if len(a) != m or len(b) != m or not (np.all(np.abs(a) == 1) and np.all(np.abs(b) == 1)):
        raise ValueError("the pair must be two +-1 sequences of length (q - 1)/2")
    sgn = 1 if periodic else -1
    afs = [int(a[:m - s] @ a[s:] + sgn * (a[m - s:] @ a[:s]) + b[:m - s] @ b[s:] + sgn * (b[m - s:] @ b[:s]))
           for s in range(m)]
    if afs[0] != 2 * m or any(afs[1:]):
        raise ValueError("the pair is not a (nega)periodic Golay pair")
    def twisted(x):
        idx = np.arange(m)[None, :] - np.arange(m)[:, None]
        return np.where(idx >= 0, 1, sgn) * x[idx % m]
    A, B = twisted(a), twisted(b)
    K = np.block([[A, B], [-B.T, A.T]]).astype(np.int8)
    chi = -np.ones(q, dtype=np.int64)
    chi[0] = 0
    chi[(np.arange(1, q, dtype=np.int64) ** 2) % q] = 1
    pts = np.array([0] + [pow(g, 2 * j, q) for j in range(m)] + [pow(g, 2 * j + 1, q) for j in range(m)], dtype=np.int64)
    C = np.ones((q + 1, q + 1), dtype=np.int8)
    C[0, 0] = 0
    C[1:, 1:] = chi[(pts[None, :] - pts[:, None]) % q]
    H = miyamoto(C, K)
    T = [np.array([SIGN[c] for c in s], dtype=np.int64) for s in cert["T"]]
    t = len(T[0])
    if len(T) != 4 or t != cert["t"] or 4 * q * t != cert["order"] or not np.all(np.count_nonzero(np.array(T), axis=0) == 1):
        raise ValueError("T-matrices: four rows of length t with supports partitioning Z_t, and order 4qt, required")
    paf = paf_sum(T)
    if paf[0] != t or np.any(paf[1:] != 0):
        raise ValueError("T-matrix condition fails")
    Rb = np.zeros(q, dtype=np.int64)
    for half in range(2):
        for j in range(m):
            Rb[1 + half * m + j] = 1 + (half if (j == 0 or periodic) else 1 - half) * m + (-j) % m
    R = (((-np.arange(t)) % t)[:, None] * q + Rb[None, :]).reshape(-1)
    n = t * q
    Aj = []
    for j in range(4):
        X = np.zeros((n, n), dtype=np.int8)
        for k in range(4):
            Nkj = H[k * q:(k + 1) * q, j * q:(j + 1) * q]
            for u in range(t):
                for w in range(t):
                    if T[k][(w - u) % t]:
                        X[u * q:(u + 1) * q, w * q:(w + 1) * q] = T[k][(w - u) % t] * Nkj
        Aj.append(X)
    AR, ATR = [x[:, R] for x in Aj], [x.T[:, R] for x in Aj]
    rows = [[Aj[0], AR[1], AR[2], AR[3]], [-AR[1], Aj[0], ATR[3], -ATR[2]],
            [-AR[2], -ATR[3], Aj[0], ATR[1]], [-AR[3], ATR[2], -ATR[1], Aj[0]]]
    G = np.empty((4 * n, 4 * n), dtype=np.int8)
    for i in range(4):
        for j in range(4):
            G[i * n:(i + 1) * n, j * n:(j + 1) * n] = rows[i][j]
    return G, None


def build_cooper_wallis_miyamoto(cert):
    """Cooper-Wallis with Williamson-type matrices of a prime order q = 1 (mod 8) from Seberry-Yamada's Corollary
    8.8(ii) (JCMCC 7 (1990), Corollary 29(ii); Seberry-Yamada's own, through Miyamoto's lemma), with a symmetric
    Hadamard matrix M of order
    n = (q - 1)/2 as input. "input" names the certificate of M (rebuilt and checked here, its matrix pinned by
    "input_sha256"); "g" is a primitive root mod q. With s_i = g^(2i), t_j = g^(1 - 2j) (i, j < n) and chi the
    Legendre symbol, U_1 = -[chi(s_i - s_j)], U_2 = [chi(s_i - t_j)], U_3 = U_4 = 0, V = (I, 0, M, M), and

        W_1 = [[1, -e], [-e^T, S_1]],   W_k = [[1, e], [e^T, S_k]] (k > 1),   S_k = U_k (x) J_2 + V_k (x) [[1, -1], [-1, 1]].

    Checked here: the W_k are +-1, symmetric, pairwise amicable and sum W_k W_k^T = 4q I (exact float64 products).
    They are substituted into the OD(4t; t, t, t, t) of the T-matrices as in `assemble_cooper_wallis`."""
    q, g = cert["q"], cert["g"]
    if q % 8 != 1 or any(q % d == 0 for d in range(2, math.isqrt(q) + 1)):
        raise ValueError("q must be a prime = 1 (mod 8)")
    if g % q == 0 or any(pow(g, (q - 1) // p, q) == 1 for p in range(2, q) if (q - 1) % p == 0 and
                         all(p % d for d in range(2, math.isqrt(p) + 1))):
        raise ValueError("g is not a primitive root mod q")
    inp = json.loads((ROOT / cert["input"]).read_text())
    M, _, digest = build_checked(inp)
    n = (q - 1) // 2
    if digest != cert["input_sha256"] or M.shape[0] != n or not np.array_equal(M, M.T):
        raise ValueError("the input is not the stated symmetric Hadamard matrix of order (q - 1)/2")
    M = M.astype(np.int64)
    chi = -np.ones(q, dtype=np.int64)
    chi[0] = 0
    chi[(np.arange(1, q, dtype=np.int64) ** 2) % q] = 1
    s = np.array([pow(g, 2 * i, q) for i in range(n)], dtype=np.int64)
    t = np.array([pow(g, (1 - 2 * j) % (q - 1), q) for j in range(n)], dtype=np.int64)
    Z, I = np.zeros((n, n), dtype=np.int64), np.eye(n, dtype=np.int64)
    U = [-chi[(s[:, None] - s[None, :]) % q], chi[(s[:, None] - t[None, :]) % q], Z, Z]
    V = [I, Z, M, M]
    J2, D2 = np.ones((2, 2), dtype=np.int64), np.array([[1, -1], [-1, 1]], dtype=np.int64)
    W = []
    for k in range(4):
        b = -1 if k == 0 else 1
        Wk = np.empty((q, q), dtype=np.int64)
        Wk[0, 0], Wk[0, 1:], Wk[1:, 0] = 1, b, b
        Wk[1:, 1:] = np.kron(U[k], J2) + np.kron(V[k], D2)
        W.append(Wk)
    F = [x.astype(np.float64) for x in W]
    if not all(np.all(np.abs(x) == 1) and np.array_equal(x, x.T) for x in W) or \
            not all(np.array_equal(F[k] @ F[l].T, F[l] @ F[k].T) for k in range(4) for l in range(k + 1, 4)) or \
            not np.array_equal(sum(x @ x.T for x in F), 4.0 * q * np.eye(q)):
        raise ValueError("Williamson-type condition fails")
    T = [np.array([SIGN[c] for c in r], dtype=np.int64) for r in cert["t_matrices"]]
    if len(T[0]) != cert["t"] or 4 * cert["t"] * q != cert["order"]:
        raise ValueError("order is not 4tq")
    return assemble_cooper_wallis(T, W), None


def build_cooper_wallis_bordered(cert):
    """Cooper-Wallis with bordered Williamson-type matrices (bordered_williamson_type): T-matrices of order t under
    "t_matrices" and the four sequences on Z_(w-1) under "bordered_williamson" give a Hadamard matrix of order 4tw."""
    T = [np.array([SIGN[c] for c in s], dtype=np.int64) for s in cert["t_matrices"]]
    seq = [np.array([SIGN[c] for c in s], dtype=np.int64) for s in cert["bordered_williamson"]]
    return assemble_cooper_wallis(T, bordered_williamson_type(seq)), None


def build_spence(cert):
    """Spence's bordered construction. Four +-1 sequences on an abelian group G of order v
    (elements in mixed-radix order of `group`) whose autocorrelations sum to -4 off the
    identity, with row sums (2, 0, 0, 0), give a Hadamard matrix of order 4v + 4:

        H = [ K        Q (x) 1^T ]
            [ P (x) 1  M         ]

    where M is the group-developed Goethals-Seidel array, P and Q are the 4 x 4 Hadamard
    matrices in the certificate, and K = -QP/2."""
    mods = cert["group"]
    v = int(np.prod(mods))
    coords = np.array(np.unravel_index(np.arange(v), mods)).T
    index = lambda c: np.ravel_multi_index(tuple(np.moveaxis(c % mods, -1, 0)), mods)
    seqs = [np.array([SIGN[c] for c in s], dtype=np.int64) for s in cert["sequences"]]
    if [int(x.sum()) for x in seqs] != [2, 0, 0, 0]:
        raise ValueError("row sums must be (2, 0, 0, 0)")
    total = np.zeros(v, dtype=np.int64)
    for s in range(v):
        shift = index(coords + coords[s])
        total[s] = sum(int(x @ x[shift]) for x in seqs)
    if total[0] != 4 * v or np.any(total[1:] != -4):
        raise ValueError("sequences do not satisfy the Spence condition")
    sub = index(coords[None, :, :] - coords[:, None, :])          # X[g, h] = x(h - g)
    R = np.zeros((v, v), dtype=np.int64)
    R[np.arange(v), index(-coords)] = 1
    A, B, C, D = (x[sub] for x in seqs)
    M = np.block([[A, B @ R, C @ R, D @ R],
                  [-B @ R, A, D.T @ R, -C.T @ R],
                  [-C @ R, -D.T @ R, A, B.T @ R],
                  [-D @ R, C.T @ R, -B.T @ R, A]])
    P, Q = np.array(cert["P"], dtype=np.int64), np.array(cert["Q"], dtype=np.int64)
    one = np.ones((1, v), dtype=np.int64)
    K = -(Q @ P) // 2
    return np.vstack([np.hstack([K, np.kron(Q, one)]), np.hstack([np.kron(P, one.T), M])]), None


def build_goethals_seidel_group(cert):
    """Goethals-Seidel array over a non-cyclic abelian group G = Z_m1 x ... x Z_mk.

    The four sequences are listed over G in mixed-radix order. Block X is the group-developed matrix
    X[g, h] = x(h - g), and R[g, h] = 1 when g + h = 0. Each sequence must be invariant under the stated
    automorphisms (integer matrices acting on coordinates).
    """
    mods = np.array(cert["group"], dtype=np.int64)
    v = int(np.prod(mods))
    coords = np.array(np.unravel_index(np.arange(v), tuple(int(m) for m in mods))).T.astype(np.int64)
    index = lambda C: np.ravel_multi_index(tuple((C % mods).T), tuple(int(m) for m in mods))
    seqs = [np.array([SIGN[c] for c in s], dtype=np.int64) for s in cert["sequences"]]
    if any(len(x) != v for x in seqs):
        raise ValueError("sequence length differs from the group order")
    if [int(x.sum()) for x in seqs] != cert["rowsums"]:
        raise ValueError("row sums differ from the certificate")
    for M in cert.get("automorphisms", []):
        p = index(coords @ np.array(M, dtype=np.int64).T)
        if len(np.unique(p)) != v:
            raise ValueError("an automorphism is not a bijection")
        if not all(np.array_equal(x[p], x) for x in seqs):
            raise ValueError("a sequence is not invariant under the stated automorphisms")
    diff = index(coords[None, :, :] - coords[:, None, :])
    A, B, C, D = (x[diff] for x in seqs)
    R = np.zeros((v, v), dtype=np.int64)
    R[np.arange(v), index(-coords)] = 1
    H = np.block([[A, B @ R, C @ R, D @ R],
                  [-B @ R, A, D.T @ R, -C.T @ R],
                  [-C @ R, -D.T @ R, A, B.T @ R],
                  [-D @ R, C.T @ R, -B.T @ R, A]])
    return H, None



def build_propus_bordered(cert):
    """Bordered propus array: a symmetric Hadamard matrix of order 4(m + 1) over an abelian group
    G = Z_m1 x ... x Z_mk of order m.

    A, B, D are +-1 sequences over G (mixed-radix order) with A(-g) = A(g) and
    PAF_A + 2 PAF_B + PAF_D = -4 at every g != 0 (4m at 0); their row sums satisfy a^2 + 2b^2 + d^2 = 4.
    With group-developed blocks X[g, h] = x(h - g), R[g, h] = [g + h = 0] and the propus array
        P = [[A, B, B, D], [B^T, DR, -A, -BR], [B^T, -A, -DR, BR], [D^T, -BR, BR, -A]],
    the matrix H = [[L, K (x) e^T], [K^T (x) e, P]] with K a Hadamard matrix of order 4 (certificate field "K") and
    L = -K S K^T / 4 (S the same array on the row sums) is symmetric with H H^T = 4(m + 1) I: P P^T = 4(m + 1) I - 4 J
    blockwise, K^T K = 4 I, L K + K S = 0 and L^2 = 4 I."""
    mods = np.array(cert["group"], dtype=np.int64)
    m = int(np.prod(mods))
    coords = np.array(np.unravel_index(np.arange(m), tuple(int(x) for x in mods))).T.astype(np.int64)
    index = lambda C: np.ravel_multi_index(tuple((C % mods).T), tuple(int(x) for x in mods))
    seqs = [np.array([SIGN[c] for c in s], dtype=np.int64) for s in cert["sequences"]]
    if len(seqs) != 3 or any(len(x) != m for x in seqs):
        raise ValueError("need three sequences of length |G|")
    if [int(x.sum()) for x in seqs] != cert["rowsums"]:
        raise ValueError("row sums differ from the certificate")
    a, b, d = cert["rowsums"]
    if a * a + 2 * b * b + d * d != 4:
        raise ValueError("row sums do not satisfy a^2 + 2b^2 + d^2 = 4")
    for M in cert.get("automorphisms", []):
        p = index(coords @ np.array(M, dtype=np.int64).T)
        if len(np.unique(p)) != m:
            raise ValueError("an automorphism is not a bijection")
        if not all(np.array_equal(x[p], x) for x in seqs):
            raise ValueError("a sequence is not invariant under the stated automorphisms")
    A, B, D = seqs
    neg = index(-coords)
    if not np.array_equal(A, A[neg]):
        raise ValueError("the first sequence is not symmetric")
    paf = lambda x: np.array([int(np.dot(x, x[index(coords + coords[g])])) for g in range(m)], dtype=np.int64)
    tot = paf(A) + 2 * paf(B) + paf(D)
    if tot[0] != 4 * m or np.any(tot[1:] != -4):
        raise ValueError("PAF_A + 2 PAF_B + PAF_D is not 4m at 0 and -4 elsewhere")
    K = np.array(cert["K"], dtype=np.int64)
    if K.shape != (4, 4) or not np.all(np.abs(K) == 1) or not np.array_equal(K.T @ K, 4 * np.eye(4, dtype=np.int64)):
        raise ValueError("K is not a Hadamard matrix of order 4")
    S = np.array([[a, b, b, d], [b, d, -a, -b], [b, -a, -d, b], [d, -b, b, -a]], dtype=np.int64)
    L4 = -(K @ S @ K.T)
    if not np.all(np.abs(L4) == 4):
        raise ValueError("-K S K^T / 4 is not a +-1 matrix")
    diff = index(coords[None, :, :] - coords[:, None, :])
    cA, cB, cD = (x[diff] for x in seqs)
    R = np.zeros((m, m), dtype=np.int64)
    R[np.arange(m), neg] = 1
    BR, DR = cB @ R, cD @ R
    P = np.block([[cA, cB, cB, cD], [cB.T, DR, -cA, -BR], [cB.T, -cA, -DR, BR], [cD.T, -BR, BR, -cA]])
    Mb = np.kron(K, np.ones((1, m), dtype=np.int64))
    H = np.block([[L4 // 4, Mb], [Mb.T, P]])
    if not np.array_equal(H, H.T):
        raise ValueError("the bordered array is not symmetric")
    return H, None


def build_propus_qbordered(cert):
    """Z_q-bordered propus array: a symmetric Hadamard matrix of order 4(m + q) over an abelian group
    G = Z_m1 x ... x Z_mk of order m with a surjective homomorphism phi: G -> Z_q, phi(g) = sum_i c_i g_i (mod q)
    (certificate fields "q" and "quotient" = [c_1, ..., c_k]; its kernel K has index q).

    A, B, D are +-1 sequences over G (mixed-radix order) with A(-g) = A(g) and
    PAF_A + 2 PAF_B + PAF_D = 4m at 0, -4q on K \\ {0} and 0 off K; their row sums satisfy a^2 + 2b^2 + d^2 = 4q.
    With group-developed blocks X[g, h] = x(h - g), R[g, h] = [g + h = 0] and the propus array
        P = [[A, B, B, D], [B^T, DR, -A, -BR], [B^T, -A, -DR, BR], [D^T, -BR, BR, -A]],
    let Mc be the same array over Z_q on the compressions y(c) = sum_{phi(g) = c} x(g), Rh a Hadamard matrix of
    order 4q (field "R") with Rh Mc Rh^T = -4q X for a +-1 matrix X, and R the 4q x 4m matrix with
    R[:, (i, g)] = Rh[:, (i, phi(g))]. Then H = [[X, R], [R^T, P]] is symmetric with H H^T = 4(m + q) I:
    P^2 + R^T R = 4(m + q) I blockwise (R^T R = I_4 (x) 4q [phi(g) = phi(h)]), R P = -X R (compression commutes with
    the group-developed blocks), R R^T = 4m I and X^2 = 4q I."""
    mods = np.array(cert["group"], dtype=np.int64)
    m = int(np.prod(mods))
    q = int(cert["q"])
    cs = np.array(cert["quotient"], dtype=np.int64)
    if len(cs) != len(mods) or any((int(c) * int(n)) % q for c, n in zip(cs, mods)):
        raise ValueError("phi is not a homomorphism G -> Z_q")
    if math.gcd(q, *[int(c) for c in cs]) != 1:
        raise ValueError("phi is not surjective")
    coords = np.array(np.unravel_index(np.arange(m), tuple(int(x) for x in mods))).T.astype(np.int64)
    index = lambda C: np.ravel_multi_index(tuple((C % mods).T), tuple(int(x) for x in mods))
    phi = (coords @ cs) % q
    seqs = [np.array([SIGN[c] for c in s], dtype=np.int64) for s in cert["sequences"]]
    if len(seqs) != 3 or any(len(x) != m for x in seqs):
        raise ValueError("need three sequences of length |G|")
    if [int(x.sum()) for x in seqs] != cert["rowsums"]:
        raise ValueError("row sums differ from the certificate")
    a, b, d = cert["rowsums"]
    if a * a + 2 * b * b + d * d != 4 * q:
        raise ValueError("row sums do not satisfy a^2 + 2b^2 + d^2 = 4q")
    A, B, D = seqs
    neg = index(-coords)
    if not np.array_equal(A, A[neg]):
        raise ValueError("the first sequence is not symmetric")
    paf = lambda x: np.array([int(np.dot(x, x[index(coords + coords[g])])) for g in range(m)], dtype=np.int64)
    tot = paf(A) + 2 * paf(B) + paf(D)
    want = np.where(phi == 0, -4 * q, 0); want[0] = 4 * m
    if not np.array_equal(tot, want):
        raise ValueError("PAF_A + 2 PAF_B + PAF_D is not 4m at 0, -4q on the kernel of phi and 0 off it")
    def propus_array(xs, diff, refl):
        cA, cB, cD = (x[diff] for x in xs)
        BR, DR = cB @ refl, cD @ refl
        return np.block([[cA, cB, cB, cD], [cB.T, DR, -cA, -BR], [cB.T, -cA, -DR, BR], [cD.T, -BR, BR, -cA]])
    ys = [np.array([int(x[phi == c].sum()) for c in range(q)], dtype=np.int64) for x in seqs]
    zq = np.arange(q)
    Rq = np.zeros((q, q), dtype=np.int64); Rq[zq, (-zq) % q] = 1
    Mc = propus_array(ys, (zq[None, :] - zq[:, None]) % q, Rq)
    Rh = np.array(cert["R"], dtype=np.int64)
    if Rh.shape != (4 * q, 4 * q) or not np.all(np.abs(Rh) == 1) or not np.array_equal(Rh @ Rh.T, 4 * q * np.eye(4 * q, dtype=np.int64)):
        raise ValueError("R is not a Hadamard matrix of order 4q")
    L = -(Rh @ Mc @ Rh.T)
    if not np.all(np.abs(L) == 4 * q):
        raise ValueError("-R Mc R^T / 4q is not a +-1 matrix")
    C3 = (coords[None, :, :] - coords[:, None, :]) % mods                 # C3[g, h] = h - g
    diff = np.ravel_multi_index(tuple(np.moveaxis(C3, -1, 0)), tuple(int(x) for x in mods))
    Rm = np.zeros((m, m), dtype=np.int64); Rm[np.arange(m), neg] = 1
    P = propus_array(seqs, diff, Rm)                 # X[g, h] = x(h - g), as Mc over Z_q
    Rb = Rh[:, [i * q + int(phi[g]) for i in range(4) for g in range(m)]]
    H = np.block([[L // (4 * q), Rb], [Rb.T, P]])
    if not np.array_equal(H, H.T):
        raise ValueError("the bordered array is not symmetric")
    return H, None


def build_bordered_gs(cert):
    """Bordered Goethals-Seidel construction over an abelian group G of order v with a
    subgroup K of index q. Four +-1 sequences on G whose autocorrelations sum to
    `autocorrelation["subgroup"]` on K \\ {0} and to `autocorrelation["elsewhere"]` off K
    give, with a suitable border of size m (K = G and m = 4 is Spence's construction),

        H = [ X  R ]      R[x, (i, g)] = Rhat[x, i, g + K]
            [ C  M ]      C[(i, g), x] = Chat[i, g + K, x]

    where M is the Goethals-Seidel array developed over G:

        [ A    B2    C2    D2  ]      A[g, h]  = a(h - g)
        [-B2   A     D2'  -C2' ]      Y2[g, h] = y(g + h)
        [-C2  -D2'   A     B2' ]      Y2'[g, h] = y(-g - h)
        [-D2   C2'  -B2'   A   ]

    Elements of G are numbered in mixed-radix order of `group` (first factor slowest),
    and the cosets of K in order of their smallest element."""
    mods = cert["group"]
    v = int(np.prod(mods))
    coords = np.array(np.unravel_index(np.arange(v), mods)).T
    index = lambda c: np.ravel_multi_index(tuple(np.moveaxis(c % mods, -1, 0)), mods)
    K, frontier = {0}, [0]
    gens = [index(np.array(g)) for g in cert["subgroup_generators"]]
    while frontier:
        x = frontier.pop()
        for g in gens:
            y = int(index(coords[x] + coords[g]))
            if y not in K:
                K.add(y)
                frontier.append(y)
    K = np.array(sorted(K))
    cos = -np.ones(v, dtype=np.int64)
    q = 0
    for g in range(v):
        if cos[g] < 0:
            cos[index(coords[g] + coords[K])] = q
            q += 1
    seqs = [np.array([SIGN[c] for c in s], dtype=np.int64) for s in cert["sequences"]]
    if [int(x.sum()) for x in seqs] != cert["rowsums"]:
        raise ValueError("row sums differ from the certificate")
    for M_ in cert.get("symmetries", []):
        img = index((coords @ np.array(M_).T))
        if not all(np.array_equal(x[img], x) for x in seqs):
            raise ValueError("sequences are not invariant under the stated symmetries")
    target = np.where(cos == 0, cert["autocorrelation"]["subgroup"], cert["autocorrelation"]["elsewhere"])
    for s in range(1, v):
        shift = index(coords + coords[s])
        if sum(int(x @ x[shift]) for x in seqs) != target[s]:
            raise ValueError("sequences do not satisfy the autocorrelation condition")
    k = cert["border"]
    X = np.array([[SIGN[c] for c in r] for r in k["X"]], dtype=np.int64)
    m = X.shape[0]
    Rhat = np.array([[SIGN[c] for c in r] for r in k["Rhat"]], dtype=np.int64).reshape(m, 4, q)
    Chat = np.array([[SIGN[c] for c in r] for r in k["Chat"]], dtype=np.int64).reshape(4, q, m)
    one = index(coords[:, None, :] - coords[None, :, :]).T          # one[g, h] = h - g
    two = index(coords[:, None, :] + coords[None, :, :])            # two[g, h] = g + h
    neg = index(-coords)
    A = seqs[0][one]
    B, C, D = (x[two] for x in seqs[1:])
    Br, Cr, Dr = (x[neg][two] for x in seqs[1:])
    M = np.block([[A, B, C, D], [-B, A, Dr, -Cr], [-C, -Dr, A, Br], [-D, Cr, -Br, A]])
    R = Rhat[:, :, cos].reshape(m, 4 * v)
    Cb = Chat[:, cos, :].reshape(4 * v, m)
    return np.block([[X, R], [Cb, M]]), seqs


def paley_conference(q):
    """Paley's symmetric conference matrix of order q + 1 for a prime q = 1 (mod 4).
    Index 0 is the point at infinity and index 1 + x is x in Z_q: C[0, 0] = 0,
    C[0, j] = C[j, 0] = 1 for j > 0, and C[1 + x, 1 + y] = chi(y - x), where chi is the
    quadratic character (chi(0) = 0)."""
    chi = -np.ones(q, dtype=np.int64)
    chi[np.arange(1, q) ** 2 % q] = 1
    chi[0] = 0
    C = np.ones((q + 1, q + 1), dtype=np.int8)
    C[0, 0] = 0
    C[1:, 1:] = circulant(chi)
    return C


def miyamoto(C, K):
    """Miyamoto's array, in the block layout of SageMath's
    hadamard_matrix_miyamoto_construction.

    C is a symmetric conference matrix of order q + 1 whose row and column 0 are
    (0, 1, ..., 1), K is a Hadamard matrix of order q - 1, and m = (q - 1)/2. First, as
    in SageMath, the k-th index j in 2..m+1 with C[1, j] = -1 is exchanged with the k-th
    index j in m+2..q with C[1, j] = 1 (rows and columns), so that row 1 becomes
    (1, 0, 1, ..., 1, -1, ..., -1). In blocks of order m write

        C[2:, 2:] = [[-C1, C2], [C2^T, C4]],        K = [[K1, K2], [-K3, K4]],

        U = [[C1,   C2, 0,    0 ],          V = [[I,    0,    K1, K2],
             [C2^T, C4, 0,    0 ],               [0,    I,    K3, K4],
             [0,    0,  C1,   C2],               [K1^T, K3^T, I,  0 ],
             [0,    0,  C2^T, C4]],              [K2^T, K4^T, 0,  I ]],

    and T_rc = [[U_rc + V_rc, U_rc - V_rc], [U_rc - V_rc, U_rc + V_rc]], 0 <= r, c < 4.
    H is the 4 x 4 array of blocks of order q

        H_rc = S_rc [[1, s e], [s e^T, T_rc]],   s = -1 if r = c, else 1,

    with e the all-ones row of length q - 1 and S the skew Hadamard matrix below."""
    q = C.shape[0] - 1
    m = (q - 1) // 2
    neg = [j for j in range(2, m + 2) if C[1, j] == -1]
    pos = [j for j in range(m + 2, q + 1) if C[1, j] == 1]
    p = np.arange(q + 1)
    p[neg], p[pos] = p[pos], p[neg]
    C = C[np.ix_(p, p)]
    C1, C2, C4 = -C[2:m + 2, 2:m + 2], C[2:m + 2, m + 2:], C[m + 2:, m + 2:]
    K1, K2, K3, K4 = K[:m, :m], K[:m, m:], -K[m:, :m], K[m:, m:]
    I, O = np.eye(m, dtype=np.int8), np.zeros((m, m), dtype=np.int8)
    U = [[C1, C2, O, O], [C2.T, C4, O, O], [O, O, C1, C2], [O, O, C2.T, C4]]
    V = [[I, O, K1, K2], [O, I, K3, K4], [K1.T, K3.T, I, O], [K2.T, K4.T, O, I]]
    S = np.array([[1, 1, 1, 1], [-1, 1, 1, -1], [-1, -1, 1, 1], [-1, 1, -1, 1]], dtype=np.int8)
    H = np.empty((4 * q, 4 * q), dtype=np.int8)
    for r in range(4):
        for c in range(4):
            s, B = S[r, c], H[r * q:(r + 1) * q, c * q:(c + 1) * q]
            B[0, 0] = s
            B[0, 1:] = B[1:, 0] = -s if r == c else s
            B[1:m + 1, 1:m + 1] = B[m + 1:, m + 1:] = s * (U[r][c] + V[r][c])
            B[1:m + 1, m + 1:] = B[m + 1:, 1:m + 1] = s * (U[r][c] - V[r][c])
    return H


def build_miyamoto(cert):
    """Miyamoto's construction (M. Miyamoto, J. Combin. Theory Ser. A 57 (1991) 86-108):
    for q = 1 (mod 4), a symmetric conference matrix of order q + 1 and a Hadamard matrix
    of order q - 1 give a Hadamard matrix of order 4q, assembled as in `miyamoto`.

    The certificate gives q, the conference matrix ("paley": Paley's matrix, which needs
    q prime here) and, under "K", the complete certificate of the Hadamard matrix of
    order q - 1. K is rebuilt with its own builder and checked like any other
    certificate, including its SHA-256."""
    q = cert["q"]
    if q % 4 != 1 or q < 5 or any(q % d == 0 for d in range(2, math.isqrt(q) + 1)):
        raise ValueError("q must be a prime = 1 (mod 4)")
    if cert["conference_matrix"] != "paley":
        raise ValueError(f"unknown conference matrix {cert['conference_matrix']!r}")
    if cert["K"]["order"] != q - 1:
        raise ValueError("K must have order q - 1")
    try:
        K = build_checked(cert["K"])[0].astype(np.int8)
    except ValueError as exc:
        raise ValueError(f"K: {exc}") from None
    C = paley_conference(q)
    if not np.array_equal(C, C.T) or np.diagonal(C).any() or not gram_is_scalar(C, q):
        raise ValueError("C is not a symmetric conference matrix")
    return miyamoto(C, K), None


def is_skew(X, block=2048):
    """True if X + X^T = 2I (checked a block of rows at a time)."""
    n = X.shape[0]
    for i in range(0, n, block):
        E = X[i:i + block].astype(np.int16) + X[:, i:i + block].T
        E[np.arange(E.shape[0]), i + np.arange(E.shape[0])] -= 2
        if E.any():
            return False
    return True


def build_skew_doubling(cert):
    """[[S, S], [-S^T, S^T]] is a skew-Hadamard matrix of order 2m when S is one of order m.
    The certificate gives, under "S", the complete certificate of S, which is rebuilt and
    checked like any other certificate and then checked to satisfy S + S^T = 2I."""
    if 2 * cert["S"]["order"] != cert["order"]:
        raise ValueError("S must have half the order")
    try:
        S = build_checked(cert["S"])[0].astype(np.int8)
    except ValueError as exc:
        raise ValueError(f"S: {exc}") from None
    if not is_skew(S):
        raise ValueError("S is not skew (S + S^T != 2I)")
    return np.block([[S, S], [-S.T, S.T]]), None


def build_yamada(cert):
    """Yamada's construction (M. Yamada, J. Austral. Math. Soc. Ser. A 46 (1989) 371-383,
    Theorem 4, which he credits to Z. Kiyasu): a prime q = 5 (mod 8) and a skew-Hadamard
    matrix Q of order (q + 3)/2 give a Hadamard matrix of order 4(q + 2).

    Let n0 = (q + 1)/2, F = GF(q)[w]/(w^2 - r) with r a non-square, xi a primitive element
    of F, psi the quadratic character of GF(q) with psi(0) = 0 and Tr(a + bw) = 2a. The
    circulants alpha, beta of order n0 with first rows psi(Tr xi^(4m)) and
    psi(Tr xi^(4m + n0)), m = 0, ..., n0 - 1, are symmetric with
    alpha alpha^T + beta beta^T = qI (Yamada's Lemma 3). Normalise Q to [[1, e^T], [-e, S + I]].
    With P = [[1, -1], [-1, 1]] and J2 the all-one matrix of order 2,

        A = P (x) beta + J2 (x) I,   B = P (x) alpha,   C = D = J2 (x) S + P (x) I

    satisfy A A^T + B B^T + C C^T + D D^T = 4(q + 2)I - 4J, and Theorem 1 gives

        H = [ L          K4 (x) e^T ]      H0 = [  A     B     C     D   ]
            [ M^T (x) e  H0         ]           [ -B^T   A^T  -D^T   C^T ]
                                                [ -C^T   D     A^T  -B   ]
                                                [ -D^T  -C     B^T   A   ]

    with L = 2I - J4, M as below and K4 = -LM/2. The certificate gives q, r, xi and, under
    "Q", the complete certificate of Q."""
    q, r, (xa, xb) = cert["q"], cert["r"], cert["xi"]
    if "poly" in cert:                                          # q = p^k: GF(q) = GF(p)[x]/(poly)
        return _build_yamada_field(cert)
    if q % 8 != 5 or any(q % d == 0 for d in range(2, math.isqrt(q) + 1)):
        raise ValueError("q must be a prime = 5 (mod 8)")
    if cert["Q"]["order"] != (q + 3) // 2:
        raise ValueError("Q must have order (q + 3)/2")
    psi = -np.ones(q, dtype=np.int64)
    psi[np.arange(1, q) ** 2 % q] = 1
    psi[0] = 0
    if psi[r % q] != -1:
        raise ValueError("r is a square")
    mul = lambda x, y: ((x[0] * y[0] + r * x[1] * y[1]) % q, (x[0] * y[1] + x[1] * y[0]) % q)

    def power(x, e):
        y = (1, 0)
        while e:
            if e & 1:
                y = mul(y, x)
            x, e = mul(x, x), e >> 1
        return y
    N, primes, m, d = q * q - 1, [], q * q - 1, 2
    while d * d <= m:                                           # prime divisors of |F*|
        if m % d == 0:
            primes.append(d)
            while m % d == 0:
                m //= d
        d += 1
    if m > 1:
        primes.append(m)
    if any(power((xa, xb), N // p) == (1, 0) for p in primes):
        raise ValueError("xi is not primitive")
    n0 = (q + 1) // 2
    alpha, beta = np.empty(n0, dtype=np.int64), np.empty(n0, dtype=np.int64)
    a, b, x4 = (1, 0), power((xa, xb), n0), power((xa, xb), 4)
    for m in range(n0):
        alpha[m], beta[m] = psi[2 * a[0] % q], psi[2 * b[0] % q]
        a, b = mul(a, x4), mul(b, x4)
    return _yamada_matrix(cert, q, alpha, beta)


def _build_yamada_field(cert):
    """build_yamada for a prime power q = p^k = 5 (mod 8). Elements of GF(q) are integers 0 .. q - 1 whose base-p
    digits (least significant first) are the coefficients of a polynomial in x modulo cert["poly"] (monic, given by
    its coefficients c_0 .. c_{k-1}: x^k = -(c_0 + c_1 x + ... + c_{k-1} x^{k-1})). r is a non-square of GF(q) and
    xi = (xa, xb) = xa + xb w a primitive element of GF(q^2) = GF(q)[w]/(w^2 - r), in the same encoding."""
    q, r, (xa, xb), poly = cert["q"], cert["r"], cert["xi"], cert["poly"]
    p, k = cert["p"], len(poly)
    if p ** k != q or q % 8 != 5 or any(p % d == 0 for d in range(2, math.isqrt(p) + 1)):
        raise ValueError("q must be a prime power p^k = 5 (mod 8)")
    if cert["Q"]["order"] != (q + 3) // 2:
        raise ValueError("Q must have order (q + 3)/2")
    digits = np.array([[(e // p ** j) % p for j in range(k)] for e in range(q)], dtype=np.int64)
    weights = p ** np.arange(k)
    enc = lambda d: int((np.asarray(d) % p) @ weights)
    def xmul(e):                                                # e * x
        d = np.zeros(k + 1, dtype=np.int64); d[1:] = digits[e]
        top = d[k]; d = d[:k] - top * np.asarray(poly, dtype=np.int64)
        return enc(d)
    # log tables from the element x; check that x generates GF(q)* (so poly is primitive, hence irreducible)
    exp = np.empty(q - 1, dtype=np.int64); log = -np.ones(q, dtype=np.int64)
    e = 1
    for i in range(q - 1):
        if log[e] >= 0:
            raise ValueError("poly is not primitive")
        exp[i], log[e] = e, i
        e = xmul(e)
    if e != 1:
        raise ValueError("poly is not primitive")
    add = lambda u, v: enc(digits[u] + digits[v])
    fm = lambda u, v: 0 if u == 0 or v == 0 else int(exp[(log[u] + log[v]) % (q - 1)])
    psi = np.where(log % 2 == 0, 1, -1).astype(np.int64); psi[0] = 0
    if psi[r] != -1:
        raise ValueError("r is a square")
    mul = lambda x, y: (add(fm(x[0], y[0]), fm(r, fm(x[1], y[1]))), add(fm(x[0], y[1]), fm(x[1], y[0])))
    def power(x, e):
        y = (1, 0)
        while e:
            if e & 1:
                y = mul(y, x)
            x, e = mul(x, x), e >> 1
        return y
    N, primes, m, d = q * q - 1, [], q * q - 1, 2
    while d * d <= m:
        if m % d == 0:
            primes.append(d)
            while m % d == 0:
                m //= d
        d += 1
    if m > 1:
        primes.append(m)
    if any(power((xa, xb), N // f) == (1, 0) for f in primes):
        raise ValueError("xi is not primitive")
    n0 = (q + 1) // 2
    alpha, beta = np.empty(n0, dtype=np.int64), np.empty(n0, dtype=np.int64)
    a, b, x4 = (1, 0), power((xa, xb), n0), power((xa, xb), 4)
    two = lambda u: add(u, u)                                   # Tr(a + b w) = 2a
    for m in range(n0):
        alpha[m], beta[m] = psi[two(a[0])], psi[two(b[0])]
        a, b = mul(a, x4), mul(b, x4)
    return _yamada_matrix(cert, q, alpha, beta)


def _yamada_matrix(cert, q, alpha, beta):
    n0 = (q + 1) // 2
    s = np.arange(n0)
    if not (np.array_equal(alpha, alpha[-s % n0]) and np.array_equal(beta, beta[-s % n0])):
        raise ValueError("alpha or beta is not symmetric")
    paf = np.array([alpha @ np.roll(alpha, -k) + beta @ np.roll(beta, -k) for k in range(n0)])
    if paf[0] != q or np.any(paf[1:] != 0):
        raise ValueError("alpha alpha^T + beta beta^T != qI")
    try:
        Q = build_checked(cert["Q"])[0].astype(np.int8)
    except ValueError as exc:
        raise ValueError(f"Q: {exc}") from None
    if not is_skew(Q):
        raise ValueError("Q is not skew (Q + Q^T != 2I)")
    d = Q[0].copy()
    d[0] = 1
    S = (Q * d[:, None] * d[None, :])[1:, 1:] - np.eye(n0, dtype=np.int8)
    P, J2, I = np.array([[1, -1], [-1, 1]], dtype=np.int8), np.ones((2, 2), dtype=np.int8), np.eye(n0, dtype=np.int8)
    A = np.kron(P, circulant(beta).astype(np.int8)) + np.kron(J2, I)
    B = np.kron(P, circulant(alpha).astype(np.int8))
    C = np.kron(J2, S) + np.kron(P, I)
    D = C
    L = 2 * np.eye(4, dtype=np.int64) - np.ones((4, 4), dtype=np.int64)
    M = np.array([[1, 1, 1, 1], [1, -1, -1, 1], [1, 1, -1, -1], [1, -1, 1, -1]], dtype=np.int64)
    K4 = -(L @ M) // 2
    n = q + 1
    H = np.empty((4 * n + 4, 4 * n + 4), dtype=np.int8)
    H[:4, :4] = L
    H[:4, 4:] = np.kron(K4, np.ones((1, n), dtype=np.int64))
    H[4:, :4] = np.kron(M.T, np.ones((n, 1), dtype=np.int64))
    for i, row in enumerate([[A, B, C, D], [-B.T, A.T, -D.T, C.T], [-C.T, D, A.T, -B], [-D.T, -C, B.T, A]]):
        for j, X in enumerate(row):
            H[4 + i * n:4 + (i + 1) * n, 4 + j * n:4 + (j + 1) * n] = X
    return H, None


def _skew_propus_od(n, blocks, rowsums):
    """The symmetric OD x P + y N of order 4n from a skew-propus quadruple (see build_skew_propus_plug); checked."""
    m = 4 * n
    seqs = []
    for reps in blocks:
        x = np.ones(n, dtype=np.int64)
        x[reps] = -1
        seqs.append(x)
    A3, A1, A0 = seqs
    if [int(x.sum()) for x in seqs] != rowsums:
        raise ValueError("row sums differ from the certificate")
    idx = (-np.arange(n)) % n
    if not np.array_equal(A3, A3[idx]):
        raise ValueError("A3 is not symmetric")
    if A0[0] != 1 or not np.array_equal(A0[1:], -A0[idx][1:]):
        raise ValueError("A0 is not of skew type")
    paf = paf_sum([A3, A1, A1, A0])
    if paf[0] != m or np.any(paf[1:] != 0):
        raise ValueError("PAF_A3 + 2 PAF_A1 + PAF_A0 is not 4n at 0 and 0 elsewhere")
    R = back_identity(n)
    A, B = circulant(A3), circulant(A1)
    D0 = circulant(A0) - np.eye(n, dtype=np.int64)
    Z = np.zeros((n, n), dtype=np.int64)
    BR, RB = B @ R, R @ B
    N = np.block([[-A, BR, BR, D0 @ R], [BR, R @ D0, A, -RB], [BR, A, -R @ D0, RB], [D0 @ R, -RB, RB, A]])
    P = np.block([[Z, Z, Z, R], [Z, R, Z, Z], [Z, Z, -R, Z], [R, Z, Z, Z]])
    _check_sym_od(P, N)
    return P, N, seqs


def build_skew_propus_plug(cert):
    """Symmetric Hadamard matrix of order m(m - 1), m = 4n, from a skew-propus quadruple of odd length n.

    A3 (symmetric), A1 and A0 (of skew type: a0[0] = 1, a0[-g] = -a0[g]) are +-1 sequences on Z_n with
    PAF_A3 + 2 PAF_A1 + PAF_A0 = 4n delta. With A = A3, B = A1, D0 = A0 - I (circulants) and R the back identity,

        N = [ -A     BR     BR     D0 R ]        P = [ 0  0   0  R ]
            [ BR     R D0   A     -RB   ]            [ 0  R   0  0 ]
            [ BR     A     -R D0   RB   ]            [ 0  0  -R  0 ]
            [ D0 R  -RB     RB     A    ]            [ R  0   0  0 ]

    gives the symmetric orthogonal design x P + y N of type (1, m - 1) (Dokovic, arXiv:2508.17141, Theorem 1,
    Balonin array): P and N are symmetric with complementary supports, P^2 = I, N^2 = (m - 1) I, PN + NP = 0;
    all of this is checked. M = P + N is a symmetric Hadamard matrix of order m. Moving its first diagonal entry
    -1 to position 0 and normalising row and column 0 to (-1, 1, ..., 1) by a diagonal +-1 similarity leaves a
    symmetric D of order m - 1 with D J = J and D^2 = m I - J (checked), and

        H = P (x) J_{m-1} + N (x) D

    is symmetric with H H^T = I (x) ((m - 1) J) + (PN + NP) (x) J + N^2 (x) D^2 = m (m - 1) I."""
    n = cert["n"]
    m = 4 * n
    P, N, seqs = _skew_propus_od(n, cert["blocks"], cert["rowsums"])
    D = _core(P + N)
    J = np.ones((m - 1, m - 1), dtype=np.int64)
    H = np.kron(P.astype(np.int8), J.astype(np.int8)) + np.kron(N.astype(np.int8), D.astype(np.int8))
    return H, seqs


def build_skew_amicable_plug(cert):
    """Skew-Hadamard matrix of order m(m - 1), m = 4n, by Seberry's theorem (Ars Combin. 1978; amicable Hadamard
    matrices of order m and a skew-type OD(m; 1, 1, m - 2)), with the amicable pair taken from a skew-propus
    quadruple of length n.

    1. The skew-propus quadruple ("blocks") gives the symmetric OD x P + y N of type (1, m - 1) of
       build_skew_propus_plug. S = I + PN is skew-Hadamard, M = P + N symmetric Hadamard, and S M^T = M S^T.
    2. Let r be the first index with P[r, r] = +1 and e = S[r]. Conjugating S and M by diag(e) makes row r of both
       all ones; with r moved first, S = [[1, e], [-e, I + W]] and M = [[1, e], [e, D]], where W is skew, W e = 0,
       (I + W)(I + W)^T = mI - J, D is symmetric, D e = -e, D^2 = mI - J and W D + D W = 0 (all checked).
    3. "g": four +-1 sequences of length n, k0 and k1 of skew type and g2, g3 symmetric, with
       PAF_k0 + PAF_k1 + PAF_g2 + PAF_g3 = 4n delta. The Goethals-Seidel array of the blocks a I + c(K0 - I),
       b I + c(K1 - I), c G2, c G3 is a skew-type orthogonal design a I + b B + c C of type (1, 1, m - 2)
       (B, C skew, checked).
    4. H = I (x) (I + W) + B (x) J + C (x) D. Then H + H^T = 2I, and
       H H^T = I (x) [(I + W)(I + W)^T + J^2 + (m - 2) D^2] + (B + B^T) (x) J - C (x) (W D + D W)
               + (B C^T + C B^T) (x) J D = I (x) [mI - J + (m - 1) J + (m - 2)(mI - J)] = m (m - 1) I.
    Instead of the skew-propus quadruple, "paley": q (a prime = 3 (mod 4)) takes the amicable pair from Paley's
    symmetric OD(q + 1; 1, q) (the classical amicable order q + 1).
    General form (Seberry's Theorem 7): instead of "g", "od" = {"type": [1, k, N - k - 1], "blocks": four
    length-t lists of codes 1 (a, only at the first entry of the first block), +-2 (b), +-3 (c)} gives in the
    Goethals-Seidel array a skew-type OD aI + bB + cC of order N = 4t = k m (checked), and then
    H H^T = I (x) [mI - J + k(m - 1) J + (N - k - 1)(mI - J)] = N (m - 1) I."""
    if "paley" in cert:
        # amicable pair from Paley's symmetric OD(q + 1; 1, q), q = 3 (mod 4) (the classical amicable order q + 1):
        # "paley": q for a prime q, or {"p": p, "poly": [...]} for q = p^k
        pal = cert["paley"]
        P, N = _paley_sym_od(pal) if isinstance(pal, int) else _paley_sym_od_field(pal["p"], pal["poly"])
        _check_sym_od(P, N)
        m = P.shape[0]
        n, seqs = None, []
    else:
        n = cert["n"]
        m = 4 * n
        P, N, seqs = _skew_propus_od(n, cert["blocks"], cert["rowsums"])
    I = np.eye(m, dtype=np.int64)
    S = I + P @ N
    M = P + N
    if not (np.array_equal(S + S.T, 2 * I) and np.array_equal(S @ M.T, M @ S.T)):
        raise ValueError("S = I + PN and M = P + N are not an amicable pair with S skew")
    r = int(np.flatnonzero(np.diag(P) == 1)[0])
    e = S[r].copy()
    S1, M1 = S * e[:, None] * e[None, :], M * e[:, None] * e[None, :]
    order = [r] + [i for i in range(m) if i != r]
    S1, M1 = S1[np.ix_(order, order)], M1[np.ix_(order, order)]
    if not ((S1[0] == 1).all() and (M1[0] == 1).all()):
        raise ValueError("normalisation failed")
    W = S1[1:, 1:] - np.eye(m - 1, dtype=np.int64)
    D = M1[1:, 1:]
    I1, J1, e1 = np.eye(m - 1, dtype=np.int64), np.ones((m - 1, m - 1), dtype=np.int64), np.ones(m - 1, dtype=np.int64)
    if not (np.array_equal(W.T, -W) and not (W @ e1).any() and np.array_equal((I1 + W) @ (I1 + W).T, m * I1 - J1)
            and np.array_equal(D, D.T) and np.array_equal(D @ e1, -e1) and np.array_equal(D @ D, m * I1 - J1)
            and not (W @ D + D @ W).any()):
        raise ValueError("the cores W, D are not amicable cores")
    if "od" in cert:
        typ, blocks = cert["od"]["type"], cert["od"]["blocks"]
        t = len(blocks[0])
        No = 4 * t
        if (len(blocks) != 4 or any(len(b) != t for b in blocks) or typ[0] != 1 or sum(typ) != No
                or typ[1] * m != No):
            raise ValueError("od: need four blocks of one length t, type (1, k, 4t - k - 1) and 4t = k m")
        if any(c not in (1, 2, -2, 3, -3) for b in blocks for c in b) or [c == 1 for b in blocks for c in b].count(True) != 1 \
                or blocks[0][0] != 1:
            raise ValueError("od: codes must be 1 (once, first entry of block 1), +-2 or +-3")
        parts = {}
        if "group" in cert["od"]:
            # blocks over an abelian group (mixed-radix element order, see abelian_group), R = negation
            tg, add, neg = abelian_group(cert["od"]["group"])
            if tg != t:
                raise ValueError("od: the group order differs from the block length")
            for v in (1, 2, 3):
                parts[v] = goethals_seidel_group(add, neg, *[group_matrix(add, neg, [int(np.sign(c)) if abs(c) == v
                                                                                    else 0 for c in b]) for b in blocks])
        else:
            for v in (1, 2, 3):
                parts[v] = goethals_seidel(*[circulant(np.array([int(np.sign(c)) if abs(c) == v else 0 for c in b],
                                                                dtype=np.int64)) for b in blocks])
        Am, Bm, Cm = parts[1], parts[2], parts[3]
        Io = np.eye(No, dtype=np.int64)
        k, c = typ[1], typ[2]
        if not (np.array_equal(Am, Io) and np.array_equal(Bm.T, -Bm) and np.array_equal(Cm.T, -Cm)
                and np.array_equal(Bm @ Bm.T, k * Io) and np.array_equal(Cm @ Cm.T, c * Io)
                and not (Bm @ Cm.T + Cm @ Bm.T).any()
                and np.array_equal(np.abs(Am) + np.abs(Bm) + np.abs(Cm), np.ones((No, No), dtype=np.int64))):
            raise ValueError(f"a I + b B + c C is not a skew-type OD({No}; 1, {k}, {c})")
        g = []
    else:
        if n is None:
            raise ValueError("G-matrices need the skew-propus form of the amicable pair")
        g = []
        for reps in cert["g"]:
            x = np.ones(n, dtype=np.int64)
            x[reps] = -1
            g.append(x)
        k0, k1, g2, g3 = g
        idx = (-np.arange(n)) % n
        for k in (k0, k1):
            if k[0] != 1 or not np.array_equal(k[1:], -k[idx][1:]):
                raise ValueError("k0 or k1 is not of skew type")
        for x in (g2, g3):
            if not np.array_equal(x, x[idx]):
                raise ValueError("g2 or g3 is not symmetric")
        paf = paf_sum(g)
        if paf[0] != m or np.any(paf[1:] != 0):
            raise ValueError("the four sequences of g do not have PAF sum 4n delta")
        In = np.eye(n, dtype=np.int64)
        Zn = np.zeros((n, n), dtype=np.int64)
        Bm = goethals_seidel(Zn, In, Zn, Zn)
        Cm = goethals_seidel(circulant(k0) - In, circulant(k1) - In, circulant(g2), circulant(g3))
        if not (np.array_equal(Bm.T, -Bm) and np.array_equal(Cm.T, -Cm) and np.array_equal(Bm @ Bm.T, I)
                and np.array_equal(Cm @ Cm.T, (m - 2) * I) and not (Bm @ Cm.T + Cm @ Bm.T).any()
                and np.array_equal(np.abs(Bm) + np.abs(Cm) + I, np.ones((m, m), dtype=np.int64))):
            raise ValueError("a I + b B + c C is not a skew-type OD(m; 1, 1, m - 2)")
    H = (np.kron(np.eye(Bm.shape[0], dtype=np.int8), (I1 + W).astype(np.int8))
         + np.kron(Bm.astype(np.int8), J1.astype(np.int8)) + np.kron(Cm.astype(np.int8), D.astype(np.int8)))
    return H, seqs + g


def _is_prime(q):
    return q > 1 and all(q % p for p in range(2, int(math.isqrt(q)) + 1))


def _check_sym_od(P, N):
    """x P + y N is a symmetric orthogonal design of type (1, m - 1)."""
    m = P.shape[0]
    I = np.eye(m, dtype=np.int64)
    if not (np.array_equal(P, P.T) and np.array_equal(N, N.T)):
        raise ValueError("P or N is not symmetric")
    if not np.array_equal(np.abs(P) + np.abs(N), np.ones((m, m), dtype=np.int64)):
        raise ValueError("P and N do not have complementary supports")
    if not (np.array_equal(P @ P, I) and np.array_equal(N @ N, (m - 1) * I) and not (P @ N + N @ P).any()):
        raise ValueError("x P + y N is not an orthogonal design of type (1, m - 1)")


def _core(M):
    """D with D symmetric, D J = J, D^2 = m I - J from a symmetric Hadamard matrix M of order m: move the first
    diagonal entry -1 to position 0 and normalise row and column 0 to (-1, 1, ..., 1) by a diagonal +-1 similarity."""
    m = M.shape[0]
    r = int(np.flatnonzero(np.diag(M) == -1)[0])
    order = [r] + [i for i in range(m) if i != r]
    M = M[np.ix_(order, order)]
    e = M[0].copy()
    e[0] = 1
    M = M * e[:, None] * e[None, :]
    D = M[1:, 1:]
    J = np.ones((m - 1, m - 1), dtype=np.int64)
    if not (np.array_equal(D, D.T) and np.array_equal(D @ J, J) and np.array_equal(D @ D, m * np.eye(m - 1, dtype=np.int64) - J)):
        raise ValueError("the core D is not symmetric with D J = J and D^2 = m I - J")
    return D


def _paley_sym_od(q):
    """Paley's skew-Hadamard matrix I + C of order q + 1 (q = 3 mod 4 prime; index 0 is infinity, 1 + x is x) and
    P = (-1) + (x -> -x): P C P = -C, so x P + y C P is a symmetric OD(q + 1; 1, q)."""
    if not (_is_prime(q) and q % 4 == 3):
        raise ValueError("q must be a prime = 3 (mod 4)")
    chi = -np.ones(q, dtype=np.int64)
    chi[np.arange(1, q) ** 2 % q] = 1
    chi[0] = 0
    m = q + 1
    C = np.zeros((m, m), dtype=np.int64)
    C[0, 1:] = 1
    C[1:, 0] = -1
    C[1:, 1:] = circulant(chi)
    P = np.zeros((m, m), dtype=np.int64)
    P[0, 0] = -1
    P[1 + np.arange(q), 1 + (-np.arange(q)) % q] = 1
    return P, C @ P


def _paley_sym_od_field(p, poly):
    """_paley_sym_od for a prime power q = p^k = 3 (mod 4). Elements of GF(q) are integers 0 .. q - 1 whose base-p
    digits (least significant first) are the coefficients of a polynomial in x modulo the monic poly given by its
    coefficients c_0 .. c_{k-1} (x^k = -(c_0 + ... + c_{k-1} x^{k-1})); poly must be primitive (checked)."""
    k = len(poly)
    q = p ** k
    if not (_is_prime(p) and q % 4 == 3):
        raise ValueError("q must be a prime power = 3 (mod 4)")
    digits = np.array([[(e // p ** j) % p for j in range(k)] for e in range(q)], dtype=np.int64)
    weights = p ** np.arange(k)
    enc = lambda d: int((np.asarray(d) % p) @ weights)
    def xmul(e):
        d = np.zeros(k + 1, dtype=np.int64); d[1:] = digits[e]
        top = d[k]; d = d[:k] - top * np.asarray(poly, dtype=np.int64)
        return enc(d)
    log = -np.ones(q, dtype=np.int64)
    e = 1
    for i in range(q - 1):
        if log[e] >= 0:
            raise ValueError("poly is not primitive")
        log[e] = i
        e = xmul(e)
    if e != 1:
        raise ValueError("poly is not primitive")
    chi = np.where(log % 2 == 0, 1, -1).astype(np.int64)
    chi[0] = 0
    neg = np.array([enc(-digits[u]) for u in range(q)], dtype=np.int64)
    sub = np.array([[enc(digits[y] - digits[x]) for y in range(q)] for x in range(q)], dtype=np.int64)   # y - x
    m = q + 1
    C = np.zeros((m, m), dtype=np.int64)
    C[0, 1:] = 1
    C[1:, 0] = -1
    C[1:, 1:] = chi[sub]
    P = np.zeros((m, m), dtype=np.int64)
    P[0, 0] = -1
    P[1 + np.arange(q), 1 + neg] = 1
    return P, C @ P


def _paley2_doubled(q, k):
    """Symmetric Hadamard matrix of order 2^(k+1) (q + 1): Paley II [[C + I, C - I], [C - I, -C - I]] for a prime
    q = 1 (mod 4), then k Kronecker products with [[1, 1], [1, -1]]."""
    if not (_is_prime(q) and q % 4 == 1):
        raise ValueError("q must be a prime = 1 (mod 4)")
    C = paley_conference(q).astype(np.int64)
    I = np.eye(q + 1, dtype=np.int64)
    M = np.block([[C + I, C - I], [C - I, -C - I]])
    for _ in range(k):
        M = np.kron(np.array([[1, 1], [1, -1]], dtype=np.int64), M)
    return M


def build_sym_od_plug(cert):
    """Symmetric Hadamard matrix from a symmetric orthogonal design x P + y N of type (1, m - 1) (P, N symmetric with
    complementary supports, P^2 = I, N^2 = (m - 1) I, PN + NP = 0) by H = P (x) X + N (x) D, where D (order r) is the
    core of a symmetric Hadamard matrix (D symmetric, D J = J, D^2 = (r + 1) I - J) and X is symmetric, commutes with
    D, and X^2 + (m - 1) D^2 = m r I:

      plug "J":       X = J,      r = m - 1, D from P + N itself          (order m (m - 1));
      plug "J-2I":    X = J - 2I, r = m + 3, D from a symmetric Hadamard matrix of order m + 4   (order m (m + 3)),

    since (J - 2I)^2 = 4I + (r - 4) J and (m - 1)((r + 1) I - J) add up to m r I exactly when r = m + 3. Then
    H H^T = P^2 (x) X^2 + (PN + NP) (x) XD + N^2 (x) D^2 = I (x) (X^2 + (m - 1) D^2) = m r I, and H = H^T.

    "design": {"type": "paley", "q": q} (q = 3 mod 4 prime). "hadamard" (plug "J-2I"): {"type": "paley2", "q": q,
    "doublings": k}. Everything (the design identities, the core identities, H H^T, H = H^T) is checked."""
    des = cert["design"]
    if des["type"] == "paley":
        P, N = _paley_sym_od(des["q"])
    else:
        raise ValueError("unknown design type")
    _check_sym_od(P, N)
    m = P.shape[0]
    if cert["plug"] == "J":
        D = _core(P + N)
        X = np.ones((m - 1, m - 1), dtype=np.int64)
    elif cert["plug"] == "J-2I":
        h = cert["hadamard"]
        if h["type"] != "paley2":
            raise ValueError("unknown Hadamard source")
        M = _paley2_doubled(h["q"], h["doublings"])
        if M.shape[0] != m + 4 or not np.array_equal(M, M.T) or not np.array_equal(M @ M, (m + 4) * np.eye(m + 4, dtype=np.int64)):
            raise ValueError("the matrix of order m + 4 is not a symmetric Hadamard matrix")
        D = _core(M)
        X = np.ones((m + 3, m + 3), dtype=np.int64) - 2 * np.eye(m + 3, dtype=np.int64)
    else:
        raise ValueError("unknown plug")
    r = D.shape[0]
    if not (np.array_equal(X @ D, D @ X) and np.array_equal(X @ X + (m - 1) * (D @ D), m * r * np.eye(r, dtype=np.int64))):
        raise ValueError("X and D do not commute or X^2 + (m - 1) D^2 != m r I")
    return np.kron(P.astype(np.int8), X.astype(np.int8)) + np.kron(N.astype(np.int8), D.astype(np.int8)), None


def build_two_block(cert):
    """Bordered two-block Hadamard matrices of order 2(w + 1) over an abelian group G of odd order w (search/szek2.py).

    "kind": "skew": a of skew type (a(0) = 1, a(-x) = -a(x)), b any +-1 function, both with row sum 1,
        PAF_a + PAF_b = -2 off 0, and with A, B group-developed (X[g, h] = x(h - g))
            [[1, 1, e, e], [-1, 1, e, -e], [-e^T, -e^T, A, B], [-e^T, e^T, -B^T, A^T]]
        is skew-Hadamard (H + H^T = 2I uses only A + A^T = 2I; b need not be symmetric).
    "kind": "symmetric": a symmetric, b any, both with row sum -1, PAF_a + PAF_b = -2 off 0, and
            [[X, E], [E^T, M]],  X = [[1, 1], [1, -1]], E = [[e, e], [e, -e]], M = [[A, B], [B^T, -A]]
        is a symmetric Hadamard matrix.
    "group": spec of G as in abelian_group (default: cyclic of order len(a)); "a", "b": the sequences as lists of
    +-1 in the mixed-radix element order. Every hypothesis and the result are checked."""
    a = np.array(cert["a"], dtype=np.int64)
    b = np.array(cert["b"], dtype=np.int64)
    w = len(a)
    tg, add, neg = abelian_group(cert.get("group", str(w)))
    if tg != w or len(b) != w or w % 2 == 0:
        raise ValueError("two-block: a and b must have the (odd) group order")
    if not (np.abs(a) == 1).all() or not (np.abs(b) == 1).all():
        raise ValueError("two-block: entries must be +-1")
    pa = np.array([int(a @ a[add[:, s]]) for s in range(w)])
    pb = np.array([int(b @ b[add[:, s]]) for s in range(w)])
    if pa[0] != w or pb[0] != w or np.any((pa + pb)[1:] != -2):
        raise ValueError("two-block: PAF_a + PAF_b is not -2 off 0")
    A, B = group_matrix(add, neg, a), group_matrix(add, neg, b)
    e = np.ones((1, w), dtype=np.int64)
    if cert["kind"] == "skew":
        if a[0] != 1 or np.any(a[neg][1:] != -a[1:]) or a.sum() != 1 or b.sum() != 1:
            raise ValueError("two-block skew: a must be of skew type and both row sums 1")
        one = np.ones((1, 1), dtype=np.int64)
        H = np.vstack([np.hstack([one, one, e, e]), np.hstack([-one, one, e, -e]),
                       np.hstack([-e.T, -e.T, A, B]), np.hstack([-e.T, e.T, -B.T, A.T])])
    elif cert["kind"] == "symmetric":
        if np.any(a[neg] != a) or a.sum() != -1 or b.sum() != -1:
            raise ValueError("two-block symmetric: a must be symmetric and both row sums -1")
        X = np.array([[1, 1], [1, -1]], dtype=np.int64)
        E = np.vstack([np.hstack([e, e]), np.hstack([e, -e])])
        H = np.block([[X, E], [E.T, np.block([[A, B], [B.T, -A]])]])
    else:
        raise ValueError("two-block: unknown kind")
    return H, [a, b]


def build_paley_blocks(cert):
    """Paley's matrix of order q + 1 = 4s (q = 3 mod 8; a prime, or a prime power with cert["field"], see
    _build_paley_blocks_field) as a 4 x 4 array of circulant s x s blocks, times
    T-matrices of order t (Seberry-Yamada 1992, Theorem 3.10), as a Goethals-Seidel quadruple of length n = st.
    Points of PG(1, q) are 0..q-1 and q (infinity), with vectors u_x = (1, x) and u_inf = (0, -1); Paley's matrix is
    P[a, b] = [a = b] + chi(det(u_a, u_b)).  The torus element g is (x -> -1/(x + tau))^power; P(ga, gb) = e(a) e(b)
    P(a, b); E(g a) = c e(a) E(a) along each g-cycle (cycles in order of their least point, each from its least point)
    makes g an automorphism of E P E, whose blocks along the cycles are circulant.  A_j(x) = sum_k T_k(x mod t)
    N_kj(x mod s); the certificate's sequences must equal the rebuilt ones."""
    q, t, tau, power, c = cert["q"], cert["t"], cert["tau"], cert["power"], cert["c"]
    if q % 8 != 3:
        raise ValueError("q must be = 3 (mod 8)")
    if not _is_prime(q):
        return _build_paley_blocks_field(cert)
    s = (q + 1) // 4
    grouped = math.gcd(s, t) != 1
    if s * t != cert["v"] or grouped != (cert.get("group") == f"Z{t} x Z{s}"):
        raise ValueError("v must be s t, with group Z_t x Z_s recorded exactly when gcd(s, t) > 1")
    chi = np.full(q, -1, dtype=np.int8)
    chi[0] = 0
    chi[(np.arange(1, q, dtype=np.int64) ** 2) % q] = 1
    u0 = np.array([1] * q + [0], dtype=np.int64)
    u1 = np.array(list(range(q)) + [q - 1], dtype=np.int64)
    P = np.empty((q + 1, q + 1), dtype=np.int8)
    for a in range(0, q + 1, 1024):                         # P[a, b] = [a = b] + chi(det(u_a, u_b)), by row blocks
        r = slice(a, min(a + 1024, q + 1))
        P[r] = chi[(u0[r, None] * u1[None, :] - u1[r, None] * u0[None, :]) % q]
    P[np.arange(q + 1), np.arange(q + 1)] += 1
    def step(x):                                           # x -> -1/(x + tau) on PG(1, q)
        if x == q:
            return 0
        den = (x + tau) % q
        return q if den == 0 else (q - 1) * pow(den, q - 2, q) % q
    g = list(range(q + 1))
    for _ in range(power):
        g = [step(x) for x in g]
    return _paley_blocks_finish(cert, P, np.array(g))


def _paley_blocks_finish(cert, P, g):
    """the common part of build_paley_blocks: cycles of g, signs, circulant blocks, T-matrices, sequences, array"""
    q, t, c = cert["q"], cert["t"], cert["c"]
    s = (q + 1) // 4
    grouped = math.gcd(s, t) != 1
    if s * t != cert["v"] or grouped != (cert.get("group") == f"Z{t} x Z{s}"):
        raise ValueError("v must be s t, with group Z_t x Z_s recorded exactly when gcd(s, t) > 1")
    cyc, seen = [], set()
    for a in range(q + 1):
        if a not in seen:
            cy = [a]
            while g[cy[-1]] != a:
                cy.append(int(g[cy[-1]]))
            seen.update(cy)
            cyc.append(cy)
    if sorted(len(cy) for cy in cyc) != [s] * 4:
        raise ValueError("g does not have four cycles of length s")
    e = (P[0] * P[g[0], g]).astype(np.int8)
    for a in range(0, q + 1, 1024):                         # P(ga, gb) = e(a) e(b) P(a, b), by row blocks
        r = slice(a, min(a + 1024, q + 1))
        if not np.array_equal(e[r, None] * e[None, :] * P[r], P[np.ix_(g[r], g)]):
            raise ValueError("g is not a signed automorphism of Paley's matrix")
    E = np.zeros(q + 1, dtype=np.int64)
    for cy in cyc:
        E[cy[0]] = 1
        for a in range(1, s):
            E[cy[a]] = c * e[cy[a - 1]] * E[cy[a - 1]]
        if c * e[cy[-1]] * E[cy[-1]] != 1:
            raise ValueError("sign normalisation does not close")
    order = np.array([a for cy in cyc for a in cy])
    E8 = E.astype(np.int8)
    B = (E8[order, None] * P[np.ix_(order, order)] * E8[None, order]).astype(np.int8)
    N = [[B[k * s, j * s:(j + 1) * s].astype(np.int64) for j in range(4)] for k in range(4)]
    for k in range(4):
        for j in range(4):
            if not np.array_equal(B[k * s:(k + 1) * s, j * s:(j + 1) * s], circulant(N[k][j])):
                raise ValueError("blocks are not circulant")
    T = [np.array([SIGN[ch] for ch in row], dtype=np.int64) for row in cert["T"]]
    if len(T) != 4 or any(len(r) != t for r in T) or not (np.abs(np.array(T)).sum(0) == 1).all():
        raise ValueError("T-matrices must be four rows of length t with disjoint supports covering Z_t")
    if not np.array_equal(sum(circulant(r) @ circulant(r).T for r in T), t * np.eye(t, dtype=np.int64)):
        raise ValueError("not T-matrices")
    if grouped:                                            # Z_t x Z_s, element (i, a) at index i s + a
        seq2 = [sum(np.outer(T[k], N[k][j]) for k in range(4)) for j in range(4)]
        seqs = [a.reshape(-1) for a in seq2]
    else:
        x = np.arange(s * t)
        seqs = [sum(T[k][x % t] * N[k][j][x % s] for k in range(4)) for j in range(4)]
    if [("".join("+" if v > 0 else "-" for v in a)) for a in seqs] != cert["sequences"]:
        raise ValueError("the sequences differ from the rebuilt ones")
    if grouped:
        paf = sum(np.array([[int((a * np.roll(np.roll(a, -di, 0), -da, 1)).sum()) for da in range(s)]
                            for di in range(t)]) for a in seq2).reshape(-1)
    else:
        paf = paf_sum(seqs)
    if paf[0] != 4 * s * t or np.any(paf[1:] != 0):
        raise ValueError("sequences are not complementary")
    if not grouped:
        return goethals_seidel(*map(circulant8, seqs)), seqs
    it, ia = np.arange(t), np.arange(s)
    dt, ds = (it[None, :] - it[:, None]) % t, (ia[None, :] - ia[:, None]) % s
    gm = [a.astype(np.int8)[dt[:, None, :, None], ds[None, :, None, :]].reshape(s * t, s * t) for a in seq2]
    neg = (((-it) % t)[:, None] * s + ((-ia) % s)[None, :]).reshape(-1)    # R[g, h] = [g + h = 0]: X R = X[:, neg]
    A, B, C, D = gm
    rv = lambda X: X[:, neg]
    return np.block([[A, rv(B), rv(C), rv(D)],
                     [-rv(B), A, rv(D.T), -rv(C.T)],
                     [-rv(C), -rv(D.T), A, rv(B.T)],
                     [-rv(D), rv(C.T), -rv(B.T), A]]), seqs


def _build_paley_blocks_field(cert):
    """build_paley_blocks for q = p^k, k > 1. cert["field"] gives p, k and f = (f_0, ..., f_(k-1)); an element is the
    integer sum d_i p^i of its coordinates in the basis 1, x, ..., x^(k-1) of GF(p)[x]/(x^k + f_(k-1) x^(k-1) + ... + f_0).
    x must have multiplicative order q - 1, which also proves that the polynomial is irreducible, so the quotient is
    the field GF(q). The rest is as for prime q: P[a, b] = [a = b] + chi(det(u_a, u_b)) on PG(1, q), the torus element
    (x -> -1/(x + tau))^power, the signs E, the circulant blocks and the sequences."""
    q, t, tau, power, c = cert["q"], cert["t"], cert["tau"], cert["power"], cert["c"]
    p, k, f = cert["field"]["p"], cert["field"]["k"], list(cert["field"]["poly"])
    if not (_is_prime(p) and k > 1 and p ** k == q and len(f) == k and all(0 <= v < p for v in f)):
        raise ValueError("field must be GF(p^k) = GF(q) with k coefficients mod p")
    pw = np.array([p ** i for i in range(k)], dtype=np.int64)
    digits = (np.arange(q, dtype=np.int64)[:, None] // pw[None, :]) % p
    enc = lambda d: int(sum(int(v) * p ** i for i, v in enumerate(d)))
    exp, d = [1], [1] + [0] * (k - 1)
    for _ in range(q - 1):                                 # powers of x
        top, d = d[-1], [0] + d[:-1]
        d = [(d[i] - top * f[i]) % p for i in range(k)]
        exp.append(enc(d))
    if exp[q - 1] != 1 or len(set(exp[:q - 1])) != q - 1 or 0 in exp:
        raise ValueError("x is not primitive modulo the given polynomial")
    log = np.full(q, -1, dtype=np.int64)
    log[np.array(exp[:q - 1])] = np.arange(q - 1)
    sub = lambda a, b: (((digits[a] - digits[b]) % p) * pw).sum(-1)
    mul = lambda a, b: 0 if a == 0 or b == 0 else exp[(log[a] + log[b]) % (q - 1)]
    inv = lambda a: exp[(-log[a]) % (q - 1)]
    chi = np.where(log % 2 == 0, 1, -1).astype(np.int8)
    chi[0] = 0
    one = 1
    mone = int(sub(0, one))
    # u_x = (1, x), u_inf = (0, -1); det(u_a, u_b) = a0 b1 - a1 b0
    u0 = [one] * q + [0]
    u1 = list(range(q)) + [mone]
    P = np.empty((q + 1, q + 1), dtype=np.int8)
    for a in range(q + 1):
        row = [int(sub(mul(u0[a], u1[b]), mul(u1[a], u0[b]))) for b in range(q + 1)]
        P[a] = chi[row]
    P[np.arange(q + 1), np.arange(q + 1)] += 1
    if not np.array_equal(P.astype(np.int64) @ P.T.astype(np.int64), (q + 1) * np.eye(q + 1, dtype=np.int64)):
        raise ValueError("Paley's matrix over GF(q) is not Hadamard")
    def step(x):                                           # x -> -1/(x + tau) on PG(1, q)
        if x == q:
            return 0
        den = int(sub(x, sub(0, tau)))                     # x + tau
        return q if den == 0 else int(sub(0, inv(den)))
    g = list(range(q + 1))
    for _ in range(power):
        g = [step(x) for x in g]
    return _paley_blocks_finish(cert, P, np.array(g))


BUILDERS = {"goethals-seidel": build_goethals_seidel, "goethals-seidel-group": build_goethals_seidel_group,
            "cooper-wallis": build_cooper_wallis,
            "cooper-wallis-near-williamson": build_cooper_wallis_near_williamson,
            "cooper-wallis-bordered": build_cooper_wallis_bordered,
            "cooper-wallis-miyamoto": build_cooper_wallis_miyamoto, "miyamoto-blocks": build_miyamoto_blocks,
            "spence": build_spence, "miyamoto": build_miyamoto,
            "bordered-goethals-seidel": build_bordered_gs, "skew-doubling": build_skew_doubling,
            "yamada": build_yamada, "propus": build_propus, "propus-bordered": build_propus_bordered,
            "propus-qbordered": build_propus_qbordered, "skew-propus-plug": build_skew_propus_plug,
            "sym-od-plug": build_sym_od_plug, "skew-amicable-plug": build_skew_amicable_plug,
            "two-block": build_two_block, "paley-blocks": build_paley_blocks}


def csv_chunks(H, rows=512):
    """The canonical CSV encoding of a +-1 matrix, a block of rows at a time. The bytes
    are those of numpy.savetxt(f, H, fmt="%d", delimiter=","), but the text of the whole
    matrix is never held in memory."""
    for i in range(0, H.shape[0], rows):
        X = H[i:i + rows]
        if not (np.abs(X) == 1).all():
            raise ValueError("entries are not all +-1")
        cell = np.empty(X.shape + (3,), dtype=np.uint8)    # "-1," or "1,"; 0 = no byte
        cell[..., 0] = np.where(X < 0, ord("-"), 0)
        cell[..., 1] = ord("1")
        cell[..., 2] = ord(",")
        cell[:, -1, 2] = ord("\n")
        cell = cell.ravel()
        yield cell[cell != 0].tobytes()


def canonical_csv(H):
    return b"".join(csv_chunks(H))


def build_checked(cert):
    """Build the matrix of a certificate and check its order, its entries, H H^T = nI and
    the SHA-256 of its canonical CSV. Returns H, the sequences (or None) and the digest."""
    H, seqs = BUILDERS[cert["construction"]](cert)
    n = H.shape[0]
    if n != cert["order"]:
        raise ValueError(f"built order {n}, certificate says {cert['order']}")
    if not is_pm1(H):
        raise ValueError("entries are not all +-1")
    if not gram_is_scalar(H, n):
        raise ValueError("H H^T != nI")
    if cert.get("skew") and not np.array_equal(H + H.T, 2 * np.eye(n, dtype=H.dtype)):
        raise ValueError("certificate claims a skew-Hadamard matrix, but H + H^T != 2I")
    if cert.get("symmetric") and not np.array_equal(H, H.T):
        raise ValueError("certificate claims a symmetric Hadamard matrix, but H != H^T")
    sha = hashlib.sha256()
    for chunk in csv_chunks(H):
        sha.update(chunk)
    digest = sha.hexdigest()
    if digest != cert["sha256"]:
        raise ValueError(f"SHA-256 mismatch: {digest}")
    return H, seqs, digest


def check(directory, write=False):
    cert = json.loads((directory / "certificate.json").read_text())
    H, seqs, digest = build_checked(cert)
    seq_file = directory / "sequences.txt"
    if seqs is not None and seq_file.exists():
        rows = [ln.strip() for ln in seq_file.read_text().splitlines() if ln.strip()]
        if rows != ["".join("+" if v > 0 else "-" for v in x) for x in seqs]:
            raise ValueError("sequences.txt does not match the certificate")
    if write:
        with open(directory / f"H{H.shape[0]}.csv", "wb") as f:
            for chunk in csv_chunks(H):
                f.write(chunk)
    return H.shape[0], digest


def check_t_matrices(path):
    T = [[SIGN[c] for c in ln.strip()] for ln in path.read_text().splitlines() if ln.strip()]
    t = len(T[0])
    if len(T) != 4 or any(len(r) != t for r in T):
        raise ValueError("expected four rows of equal length")
    if not np.all(np.count_nonzero(np.array(T), axis=0) == 1):
        raise ValueError("supports do not partition Z_t")
    paf = paf_sum(T)
    if paf[0] != t or np.any(paf[1:] != 0):
        raise ValueError("autocorrelation condition fails")
    return t


def check_near_williamson(path):
    """Four rows over + and -: A circulant, B, C, D symmetric circulants and
    A A^T + B B^T + C C^T + D D^T = 4w I (checked through the periodic autocorrelations)."""
    S = [[SIGN[c] for c in ln.strip()] for ln in path.read_text().splitlines() if ln.strip()]
    w = len(S[0])
    if len(S) != 4 or any(len(r) != w or 0 in r for r in S):
        raise ValueError("expected four +-1 rows of equal length")
    if any(r != [r[-k % w] for k in range(w)] for r in S[1:]):
        raise ValueError("rows 2 to 4 must be symmetric")
    paf = paf_sum(S)
    if paf[0] != 4 * w or np.any(paf[1:] != 0):
        raise ValueError("near-Williamson condition fails")
    return w


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("orders", nargs="*", help="orders (2092) or directories (skew/820) to check (default: all)")
    ap.add_argument("--write", action="store_true", help="write H<n>.csv next to each certificate")
    args = ap.parse_args()
    def resolve(o):
        """An order (2092: hadamard/2092) or a directory (skew/820, symmetric/532)."""
        return ROOT / "hadamard" / o if o.isdigit() else ROOT / o

    dirs = [resolve(o) for o in args.orders] if args.orders else \
        [p.parent for family in ("hadamard", "skew", "symmetric")
         for p in sorted(ROOT.glob(f"{family}/*/certificate.json"), key=lambda p: int(p.parent.name))]
    failed = 0
    for d in dirs:
        try:
            n, digest = check(d, args.write)
            c = json.loads((d / "certificate.json").read_text())
            kind = "  (skew: H + H^T = 2I)" if c.get("skew") else "  (symmetric: H = H^T)" if c.get("symmetric") else ""
            print(f"{d.relative_to(ROOT).as_posix():<16}  ok    H H^T = {n} I   sha256 {digest}{kind}")
        except Exception as exc:  # report and continue with the remaining orders
            failed += 1
            print(f"{d.relative_to(ROOT).as_posix():<16}  FAIL  {exc}")
    if not args.orders:
        for f in sorted((ROOT / "t-matrices").glob("T*.txt"), key=lambda p: int(p.stem[1:])):
            try:
                print(f"{f.relative_to(ROOT).as_posix():<16}  ok    T-matrices of order {check_t_matrices(f)}")
            except Exception as exc:
                failed += 1
                print(f"{f.relative_to(ROOT).as_posix():<16}  FAIL  {exc}")
        for f in sorted((ROOT / "near-williamson").glob("NW*.txt"), key=lambda p: int(p.stem[2:])):
            try:
                print(f"{f.relative_to(ROOT).as_posix():<16}  ok    near-Williamson matrices of order {check_near_williamson(f)}")
            except Exception as exc:
                failed += 1
                print(f"{f.relative_to(ROOT).as_posix():<16}  FAIL  {exc}")
    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()
