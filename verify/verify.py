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

The matrices of most orders above 10000 are not stored (each README says so). Rebuilding and checking one of them
takes from seconds to 40 minutes (51580, about 3.5 GB of memory); the Cooper-Wallis and near-Williamson orders are
assembled in int8 (0.6 to 1.5 GB).

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
    return assemble_cooper_wallis(T, W), None                # sum_k E_k (x) W_k, assembled in int8


def build_cooper_wallis_near_williamson(cert):
    """Cooper-Wallis with Williamson-type matrices that are not Williamson matrices.

    As in `build_cooper_wallis`, T-matrices of order t give an orthogonal design OD(4t; t, t, t, t) and
    four +-1 matrices W_1, ..., W_4 of order w with W_k W_l^T = W_l W_k^T (k != l) and
    W_1 W_1^T + ... + W_4 W_4^T = 4w I substituted for its variables give a Hadamard matrix of order 4tw.
    The Williamson-type matrices come from near-Williamson matrices (Kharaghani, Mohammadian and
    Tayfeh-Rezaie, arXiv:2605.08661): circulant +-1 matrices A, B, C, D of order w, with B, C, D symmetric
    and A arbitrary, such that A A^T + B B^T + C C^T + D D^T = 4w I. Then (A R, B, C, D), with R the
    back-identity, are Williamson-type matrices: A R is symmetric and commutes with the symmetric circulants.
    The certificate gives the first rows of A, B, C, D under "near_williamson". As in `build_cooper_wallis`
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
    +-2 at one index j and 0 elsewhere; with M_k the circulant of x_k (M_1 replaced by M_1 R, R the back-diagonal
    permutation i -> m - 1 - i) and s = r_j / 2,
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


def miyamoto_nw_k(rows, m):
    """K of order 2m = 4w for build_miyamoto_blocks from near-Williamson matrices (a, b, c, d) of order w (checked):
    miyamoto_wt_k of (A R_w, B, C, D) with R_w the back identity."""
    a, b, c, d = (np.array([SIGN[ch] for ch in r], dtype=np.int64) for r in rows)
    w = len(a)
    if 2 * w != m or any(len(x) != w for x in (b, c, d)):
        raise ValueError("near-Williamson matrices of order (q - 1)/4 required")
    if any(not np.array_equal(x, np.roll(x[::-1], 1)) for x in (b, c, d)):
        raise ValueError("b, c, d must be symmetric")
    paf = paf_sum([a, b, c, d])
    if paf[0] != 4 * w or np.any(paf[1:] != 0):
        raise ValueError("the periodic autocorrelations of a, b, c, d must sum to 4w delta")
    return miyamoto_wt_k([circulant(a) @ back_identity(w), circulant(b), circulant(c), circulant(d)], m)


def unpack_symmetric(s, w):
    """the symmetric +-1 matrix of order w whose upper triangle (row by row, diagonal included) is the bit string
    base64 `s` (np.packbits order, 1 for +1)"""
    import base64
    k = w * (w + 1) // 2
    bits = np.unpackbits(np.frombuffer(base64.b64decode(s), dtype=np.uint8))
    if len(bits) != 8 * ((k + 7) // 8) or bits[k:].any():
        raise ValueError("packed symmetric matrix of the wrong length")
    X = np.zeros((w, w), dtype=np.int64)
    X[np.triu_indices(w)] = 2 * bits[:k].astype(np.int64) - 1
    return X + np.triu(X, 1).T


def miyamoto_wt_k(X, m):
    """K of order 2m = 4w for build_miyamoto_blocks from symmetric Williamson-type matrices X_1..X_4 of order w,
    checked here: +-1 entries, symmetric, pairwise commuting, X_1^2 + X_2^2 + X_3^2 + X_4^2 = 4w I. K is their
    Williamson array, relabelled so that the half swap of every quadrant becomes j -> 1 - j."""
    w = X[0].shape[0]
    if 2 * w != m or len(X) != 4 or any(x.shape != (w, w) for x in X):
        raise ValueError("four matrices of order (q - 1)/4 required")
    if any(not np.all(np.abs(x) == 1) or not np.array_equal(x, x.T) for x in X):
        raise ValueError("Williamson-type matrices: symmetric +-1 matrices required")
    if any(not np.array_equal(X[i] @ X[j], X[j] @ X[i]) for i in range(4) for j in range(i + 1, 4)):
        raise ValueError("Williamson-type matrices must commute")
    if not np.array_equal(sum(x @ x for x in X), 4 * w * np.eye(w, dtype=np.int64)):
        raise ValueError("Williamson-type matrices: the sum of squares must be 4w I")
    X1, X2, X3, X4 = X
    K0 = np.block([[X1, X2, X3, X4], [-X2, X1, -X4, X3], [-X3, X4, X1, -X2], [-X4, -X3, X2, X1]])
    P = np.zeros((m, m), dtype=np.int64)
    used, y = set(), 0
    for j in range(m):
        if j in used:
            continue
        P[j, y] = P[(1 - j) % m, w + y] = 1
        used |= {j, (1 - j) % m}
        y += 1
    Pi = np.kron(np.eye(2, dtype=np.int64), P)
    return (Pi @ K0 @ Pi.T).astype(np.int8)


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
    If the certificate has "wt" instead, these are symmetric Williamson-type matrices X_1..X_4 of order w = (q - 1)/4
    (each a base64 bit string of its upper triangle, row by row with the diagonal, 1 for +1), checked here to be
    symmetric, pairwise commuting and with X_1^2 + ... + X_4^2 = 4w I. Then the quadrants [[X, Y], [-Y, X]] of their
    Williamson array K0 (order 4w = 2m) lie in a commutative algebra and the half swap transposes them;
    K = (I_2 (x) P) K0 (I_2 (x) P)^T, where P relabels the indices (h, y) of each quadrant (h = 0, 1; y < w) to Z_m: the
    pairs {j, 1 - j} of Z_m in order of their smaller element j become (0, y) -> j and (1, y) -> 1 - j for y = 0, 1, ...
    R_B fixes the border and sends (b, j) to (b, 1 - j). If the certificate has "nw" = (a, b, c, d), these must be
    near-Williamson matrices of order w (circulant, b, c, d symmetric, PAF sum 4w delta), and X = (A R_w, B, C, D) with
    R_w the back identity are used in the same way.
    Its H H^T = 4qt I is checked densely by build_checked."""
    q, g = cert["q"], cert["g"]
    if q % 4 != 1 or q < 5 or any(q % d == 0 for d in range(2, math.isqrt(q) + 1)):
        raise ValueError("q must be a prime = 1 (mod 4)")
    if g % q == 0 or any(pow(g, (q - 1) // p, q) == 1 for p in range(2, q) if (q - 1) % p == 0 and
                         all(p % d for d in range(2, math.isqrt(p) + 1))):
        raise ValueError("g is not a primitive root mod q")
    m = (q - 1) // 2
    periodic = "pg" in cert
    if "wt" in cert:
        K = miyamoto_wt_k([unpack_symmetric(x, m // 2) for x in cert["wt"]], m)
    elif "nw" in cert:
        K = miyamoto_nw_k(cert["nw"], m)
    else:
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
            if "nw" in cert or "wt" in cert:
                Rb[1 + half * m + j] = 1 + half * m + (1 - j) % m
            else:
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


def _symmetric_input(spec):
    """A symmetric Hadamard matrix from a short description: "paley-amicable" (the symmetric member P + N of Paley's
    amicable pair of order q + 1, q = 3 (mod 4) a prime, or {"p", "poly"} for a prime power), "paley-ii" (Paley's
    [[C + I, C - I], [C - I, -C - I]] of order 2(q + 1), q = 1 (mod 4) a prime), "sylvester" ([[1, 1], [1, -1]] to the
    power k), "kron" (Kronecker product of "factors"), "propus" (a propus quadruple (A, B, B, D) over Z_v given like a
    "propus" certificate, built and checked by build_propus) or "certificate" (a certificate in this repository, "path",
    with the SHA-256 of its matrix pinned by "sha256"). Checked: +-1, symmetric, M M = nI."""
    kind = spec["kind"]
    if kind == "paley-amicable":
        q = spec["q"]
        P, N = _paley_sym_od(q) if isinstance(q, int) else _paley_sym_od_field(q["p"], q["poly"])
        M = P + N
    elif kind == "paley-ii":
        M = _paley2_doubled(spec["q"], 0)
    elif kind == "sylvester":
        M = np.ones((1, 1), dtype=np.int64)
        for _ in range(spec["k"]):
            M = np.kron(np.array([[1, 1], [1, -1]], dtype=np.int64), M)
    elif kind == "kron":
        M = np.ones((1, 1), dtype=np.int64)
        for f in spec["factors"]:
            M = np.kron(M, _symmetric_input(f))
    elif kind == "propus":
        M = build_propus(spec)[0]
    elif kind == "certificate":
        M, _, digest = build_checked(json.loads((ROOT / spec["path"]).read_text()))
        if digest != spec["sha256"]:
            raise ValueError("the input certificate does not give the pinned matrix")
        M = M.astype(np.int64)
    else:
        raise ValueError(f"unknown input kind {kind}")
    n = M.shape[0]
    if not (np.all(np.abs(M) == 1) and np.array_equal(M, M.T) and np.array_equal(M @ M, n * np.eye(n, dtype=np.int64))):
        raise ValueError(f"the {kind} input is not a symmetric Hadamard matrix")
    return M


def build_symmetric_propus_miyamoto(cert):
    """Symmetric Hadamard matrix of order 4q, q = 1 (mod 8) a prime power, from a symmetric Hadamard matrix M of order
    n = (q - 1)/2. Seberry and Yamada's Corollary 29(ii) (JCMCC 7 (1990); their 1992 Corollary 8.8(ii), through
    Miyamoto's lemma) builds, from matrices with the properties of Yamada's (their Remark 19), here
    U_1 = -[chi(s_i - s_j)] and U_2 = [chi(s_i - t_j)] over GF(q) (s_i = g^(2i), t_j = g^(1 - 2j), g primitive, chi the
    quadratic character), with U_3 = U_4 = 0 and V = (I, 0, M, M), the matrices

        W_1 = [[1, -e], [-e^T, S_1]],   W_k = [[1, e], [e^T, S_k]] (k > 1),   S_k = U_k (x) J_2 + V_k (x) [[1, -1], [-1, 1]],

    which are symmetric, pairwise amicable, with sum W_k W_k^T = 4q I, and W_3 = W_4 because U_3 = U_4 and V_3 = V_4.
    That repeated block is hypothesis (vi) of Seberry and Balonin's propus variation of Miyamoto's lemma (Australas. J.
    Combin. 69 (2017), Theorem 3; Lemma 3 of arXiv:1512.01732), with their indices 2, 3 for our 3, 4, and the propus
    array

        H = [[W_1, W_3, W_3, W_2], [W_3, W_2, -W_1, -W_3], [W_3, -W_1, -W_2, W_3], [W_2, -W_3, W_3, -W_1]]

    is a symmetric Hadamard matrix of order 4q. "q" is a prime (with "g" a primitive root) or {"p", "poly"} for
    q = p^k (elements are integers whose base-p digits are polynomial coefficients; the class of x is primitive,
    checked). "input" describes M (_symmetric_input). Everything is checked here."""
    n, U1, U2 = _miyamoto_upart(cert)
    M = _symmetric_input(cert["input"])
    if M.shape[0] != n:
        raise ValueError("the input must have order (q - 1)/2")
    Z, I = np.zeros((n, n), dtype=np.int64), np.eye(n, dtype=np.int64)
    return _lemma26_propus([U1, U2, Z, Z], [I, Z, M, M]), None


def _miyamoto_upart(cert):
    """n = (q - 1)/2, U_1 = -[chi(s_i - s_j)] and U_2 = [chi(s_i - t_j)] (s_i = g^(2i), t_j = g^(1 - 2j), i, j < n) over
    GF(q), q = 1 (mod 8), from cert["q"] (and cert["g"]) as in build_symmetric_propus_miyamoto."""
    qs = cert["q"]
    if isinstance(qs, int):
        q = qs
        if q % 8 != 1 or not _is_prime(q):
            raise ValueError("q must be a prime = 1 (mod 8)")
        g = cert["g"]
        exp = np.empty(q - 1, dtype=np.int64)
        e = 1
        for i in range(q - 1):
            exp[i] = e
            e = e * g % q
        sub = lambda a, b: (a - b) % q
    else:
        p, poly = qs["p"], qs["poly"]
        k = len(poly)
        q = p ** k
        if q % 8 != 1 or not _is_prime(p):
            raise ValueError("q must be a prime power = 1 (mod 8)")
        digits = np.array([[(x // p ** j) % p for j in range(k)] for x in range(q)], dtype=np.int64)
        weights = p ** np.arange(k)
        exp = np.empty(q - 1, dtype=np.int64)
        e = 1
        for i in range(q - 1):
            exp[i] = e
            d = np.zeros(k + 1, dtype=np.int64)
            d[1:] = digits[e]
            e = int(((d[:k] - d[k] * np.asarray(poly, dtype=np.int64)) % p) @ weights)
        sub = lambda a, b: ((digits[a] - digits[b]) % p) @ weights
    if len(set(exp.tolist())) != q - 1 or 0 in set(exp.tolist()):
        raise ValueError("the generator is not primitive")
    log = np.full(q, -1, dtype=np.int64)
    log[exp] = np.arange(q - 1)
    chi = np.where(log % 2 == 0, 1, -1).astype(np.int64)
    chi[0] = 0
    n = (q - 1) // 2
    s = exp[2 * np.arange(n)]
    t = exp[(1 - 2 * np.arange(n)) % (q - 1)]
    return n, -chi[sub(s[:, None], s[None, :])], chi[sub(s[:, None], t[None, :])]


def build_symmetric_propus_weighing(cert):
    """Symmetric Hadamard matrix of order 4q, q = 1 (mod 8) a prime power, from ANY Hadamard matrix K of order
    n = (q - 1)/2 (NOTES.md §30). The U-part is that of build_symmetric_propus_miyamoto: H_0 is the propus array of
    W_1 = [[1, -e], [-e^T, U_1 (x) J_2]], W_2 = [[1, e], [e^T, U_2 (x) J_2]] and W_3 = W_4 = [[1, e], [e^T, 0]]
    (_miyamoto_upart). Its zeros are, in the cores of the four block rows and columns (m = q, core index
    k m + 1 + 2i + a, a = 0, 1), the 2 x 2 blocks at the nonzeros of the support pattern

        P = [[I, J, J, 0], [J, 0, I, J], [J, I, 0, J], [0, J, J, I]]

    (U_1 has zero diagonal and no other zero, U_2 has no zero). For any symmetric (0, +-1) matrix V of order 4n with
    support P and V^2 = qI, put V[k n + i, l n + j] [[1, -1], [-1, 1]] at rows k m + 1 + 2i + (0, 1), columns
    l m + 1 + 2j + (0, 1): this L has H_0 L = L H_0 = 0 (J_2 D_2 = 0, e D_2 = 0) and the same L^2 as the lift of
    Seberry-Yamada's V = (I, 0, M, M), so H = H_0 + L is a symmetric Hadamard matrix. With h = n/2, K = [A | B]
    (A the first h columns), E = diag(1^h, (-1)^h), X = [A^T; A^T] and Y = [B, B],

        V = [[E, X, -E X, 0], [X^T, 0, I, Y], [-X^T E, I, 0, -Y E], [0, Y^T, -E Y^T, E]]

    and V^2 = qI because the columns of K are orthogonal. "K" is the complete certificate of K, rebuilt and checked
    with its own builder; "q" (and "g") as in build_symmetric_propus_miyamoto. Everything is checked here."""
    n, U1, U2 = _miyamoto_upart(cert)
    if cert["K"]["order"] != n:
        raise ValueError("K must have order (q - 1)/2")
    try:
        K = build_checked(cert["K"])[0].astype(np.int8)
    except ValueError as exc:
        raise ValueError(f"K: {exc}") from None
    h, m = n // 2, 2 * n + 1
    A, B = K[:, :h], K[:, h:]
    X, Y = np.vstack([A.T, A.T]), np.hstack([B, B])
    E = np.diag(np.repeat(np.array([1, -1], dtype=np.int8), h))
    I, Z = np.eye(n, dtype=np.int8), np.zeros((n, n), dtype=np.int8)
    V = np.block([[E, X, -E @ X, Z], [X.T, Z, I, Y], [-X.T @ E, I, Z, -Y @ E], [Z, Y.T, -E @ Y.T, E]])
    P = np.block([[I, 1 - Z, 1 - Z, Z], [1 - Z, Z, I, 1 - Z], [1 - Z, I, Z, 1 - Z], [Z, 1 - Z, 1 - Z, I]])
    if not np.array_equal(V, V.T) or not np.array_equal(np.abs(V), P) or not gram_is_scalar(V, 2 * n + 1):
        raise ValueError("V is not a symmetric weighing matrix of weight q with support P")
    J2 = np.ones((2, 2), dtype=np.int8)
    W = []
    for k, Uk in enumerate([U1, U2, None, None]):
        Wk = np.zeros((m, m), dtype=np.int8)
        Wk[0, 0], Wk[0, 1:], Wk[1:, 0] = 1, -1 if k == 0 else 1, -1 if k == 0 else 1
        if Uk is not None:
            Wk[1:, 1:] = np.kron(Uk.astype(np.int8), J2)
        W.append(Wk)
    Wa, Wb, Wd = W[0], W[2], W[1]
    H = np.block([[Wa, Wb, Wb, Wd], [Wb, Wd, -Wa, -Wb], [Wb, -Wa, -Wd, Wb], [Wd, -Wb, Wb, -Wa]])
    idx = (m * np.arange(4)[:, None] + 1 + np.arange(2 * n)[None, :]).ravel()
    core = H[np.ix_(idx, idx)]
    if core[np.kron(P, J2) == 1].any():
        raise ValueError("the U-part is not zero on the support of V")
    H[np.ix_(idx, idx)] = core + np.kron(V, np.array([[1, -1], [-1, 1]], dtype=np.int8))
    if not np.array_equal(H, H.T):
        raise ValueError("H is not symmetric")
    return H, None


def _lemma26_propus(U, V):
    """Seberry-Yamada's Lemma 26 (Seberry-Balonin's propus variation): from (0, +-1) matrices U_k, V_k of order n, the
    bordered matrices W_1 = [[1, -e], [-e^T, S_1]], W_k = [[1, e], [e^T, S_k]] (k > 1) of order m = 2n + 1, with
    S_k = U_k (x) J_2 + V_k (x) [[1, -1], [-1, 1]], checked here to be +-1, symmetric, pairwise amicable, with
    sum W_k W_k^T = 4m I and W_3 = W_4; returns the propus array of A = W_1, B = W_3, D = W_2."""
    n = U[0].shape[0]
    m = 2 * n + 1
    J2, D2 = np.ones((2, 2), dtype=np.int64), np.array([[1, -1], [-1, 1]], dtype=np.int64)
    W = []
    for k in range(4):
        b = -1 if k == 0 else 1
        Wk = np.empty((m, m), dtype=np.int64)
        Wk[0, 0], Wk[0, 1:], Wk[1:, 0] = 1, b, b
        Wk[1:, 1:] = np.kron(U[k], J2) + np.kron(V[k], D2)
        W.append(Wk)
    F = [x.astype(np.float64) for x in W]
    if not all(np.all(np.abs(x) == 1) and np.array_equal(x, x.T) for x in W) or not np.array_equal(W[2], W[3]) or \
            not all(np.array_equal(F[k] @ F[l].T, F[l] @ F[k].T) for k in range(4) for l in range(k + 1, 4)) or \
            not np.array_equal(sum(x @ x.T for x in F), 4.0 * m * np.eye(m)):
        raise ValueError("the W_k are not symmetric pairwise amicable Williamson-type matrices with W_3 = W_4")
    A, B, D = (x.astype(np.int8) for x in (W[0], W[2], W[1]))
    return np.block([[A, B, B, D], [B, D, -A, -B], [B, -A, -D, B], [D, -B, B, -A]])


def unpack_conference_core(s, n):
    """the symmetric (0, +-1) matrix of order n with zero diagonal whose strict upper triangle (row by row) is the bit
    string base64 `s` (np.packbits order, 1 for +1)"""
    import base64
    k = n * (n - 1) // 2
    bits = np.unpackbits(np.frombuffer(base64.b64decode(s), dtype=np.uint8))
    if len(bits) != 8 * ((k + 7) // 8) or bits[k:].any():
        raise ValueError("packed core of the wrong length")
    X = np.zeros((n, n), dtype=np.int64)
    X[np.triu_indices(n, 1)] = 2 * bits[:k].astype(np.int64) - 1
    return X + X.T


def build_symmetric_propus_sb(cert):
    """Symmetric Hadamard matrix of order 4(q + 2) by Seberry and Balonin's Corollary 1 (Australas. J. Combin. 69
    (2017), Theorem 3 with Corollary 1; the census's Seberry-Balonin series), for a prime q = 1 (mod 4) and the core Q
    of a symmetric conference matrix of order n + 1, n = (q + 1)/2: Lemma 26 with U = (I, 0, Q, Q) and
    V = (M, N, I, I), where M, N are Yamada's symmetric circulants of order n with M^2 + N^2 = qI and M of zero
    diagonal (first rows beta, alpha of build_yamada, from "r" and "xi"), gives symmetric Williamson-type matrices of
    order q + 2 with W_3 = W_4, and the propus array. Q is given under "core" as its packed strict upper triangle
    (unpack_conference_core) with the SHA-256 of its rows; checked: symmetric, zero diagonal, Q e = 0, Q^2 = nI - J."""
    q = cert["q"]
    if q % 4 != 1 or not _is_prime(q):
        raise ValueError("q must be a prime = 1 (mod 4)")
    n = (q + 1) // 2
    core = cert["core"]
    if core["order"] != n:
        raise ValueError("the core must have order (q + 1)/2")
    Q = unpack_conference_core(core["upper"], n)
    if hashlib.sha256(Q.astype(np.int8).tobytes()).hexdigest() != core["sha256"]:
        raise ValueError("the core does not match its SHA-256")
    I = np.eye(n, dtype=np.int64)
    if (Q @ np.ones(n, dtype=np.int64)).any() or not np.array_equal(Q @ Q, n * I - 1):
        raise ValueError("Q is not the core of a symmetric conference matrix")
    alpha, beta = _yamada_alpha_beta(q, cert["r"], tuple(cert["xi"]))
    M, N = circulant(beta), circulant(alpha)
    if M[0, 0] != 0 or not (np.array_equal(M, M.T) and np.array_equal(N, N.T)) or \
            not np.array_equal(M @ M + N @ N, q * I):
        raise ValueError("Yamada's circulants are not symmetric with M^2 + N^2 = qI")
    Z = np.zeros((n, n), dtype=np.int64)
    return _lemma26_propus([I, Z, Q, Q], [M, N, I, I]), None


def build_symmetric_propus_skewcore(cert):
    """Symmetric Hadamard matrix of order 4(q + 2) from a skew-propus quadruple of length l = (q + 3)/8 (NOTES.md
    §29): Seberry-Yamada's Corollary 27(i) with the core of a skew-Hadamard matrix that has an involutive
    anti-automorphism in place of Paley's skew core. "skew_propus" gives A3 (symmetric), A1, A0 (skew type) as +-
    strings with PAF_A3 + 2 PAF_A1 + PAF_A0 = 4l delta. With A, B, D0 the circulants of A3, A1 and A0 - I and R the
    back-diagonal permutation, N = [[-A, BR, BR, D0R], [BR, RD0, A, -RB], [BR, A, -RD0, RB], [D0R, -RB, RB, A]] and
    P = [[0, 0, 0, R], [0, R, 0, 0], [0, 0, -R, 0], [R, 0, 0, 0]] are symmetric with N^2 = (4l - 1)I, P^2 = I and
    PN = -NP, so S = I + NP is skew-Hadamard with PSP = S^T. S is normalised at a fixed point of P, and its core
    B (order n = 4l - 1, skew, BJ = 0) is relabelled so that P acts as x -> -x on Z_n (RBR = B^T, R now that
    reversal). Lemma 26 with U = (I, 0, BR, BR), V = (M, N, R, R) and Yamada's circulants M, N (q = 2n - 1, "r",
    "xi") gives the W_k (_lemma26_propus). Everything is checked here."""
    s = [np.array([SIGN[c] for c in x], dtype=np.int64) for x in cert["skew_propus"]]
    l = len(s[0])
    a3, a1, a0 = s
    t = (-np.arange(l)) % l
    if any(len(x) != l for x in s) or not np.array_equal(a3, a3[t]) or a0[0] != 1 or \
            not np.array_equal(a0[t][1:], -a0[1:]):
        raise ValueError("A3 must be symmetric and A0 of skew type, all of one length")
    paf = paf_sum([a3, a1, a0]) + paf_sum([a1])
    if paf[0] != 4 * l or np.any(paf[1:] != 0):
        raise ValueError("PAF_A3 + 2 PAF_A1 + PAF_A0 is not 4l delta")
    A, B, D0 = circulant(a3), circulant(a1), circulant(a0) - np.eye(l, dtype=np.int64)
    Rl = back_identity(l)
    Z = np.zeros((l, l), dtype=np.int64)
    N = np.block([[-A, B @ Rl, B @ Rl, D0 @ Rl], [B @ Rl, Rl @ D0, A, -Rl @ B], [B @ Rl, A, -Rl @ D0, Rl @ B],
                  [D0 @ Rl, -Rl @ B, Rl @ B, A]])
    P = np.block([[Z, Z, Z, Rl], [Z, Rl, Z, Z], [Z, Z, -Rl, Z], [Rl, Z, Z, Z]])
    m4 = 4 * l
    I4 = np.eye(m4, dtype=np.int64)
    if not (np.array_equal(N, N.T) and np.array_equal(P, P.T) and np.array_equal(P @ P, I4) and
            np.array_equal(N @ N, (m4 - 1) * I4) and not (P @ N + N @ P).any()):
        raise ValueError("xP + yN is not a symmetric orthogonal design of type (1, 4l - 1)")
    S = I4 + N @ P
    if not (np.array_equal(S + S.T, 2 * I4) and np.array_equal(S @ S.T, m4 * I4) and np.array_equal(P @ S @ P, S.T)):
        raise ValueError("S = I + NP is not a skew-Hadamard matrix with anti-automorphism P")
    perm = np.argmax(np.abs(P), axis=0)
    n = m4 - 1
    for i0 in [i for i in range(m4) if perm[i] == i]:
        d = np.where(np.arange(m4) == i0, 1, S[i0])
        order = [i0] + [j for j in range(m4) if j != i0]
        Sn = ((d[:, None] * S) * d[None, :])[np.ix_(order, order)]
        W = Sn[1:, 1:] - np.eye(n, dtype=np.int64)
        pos = {v: k for k, v in enumerate(order)}
        sig = [pos[perm[order[k + 1]]] - 1 for k in range(n)]
        fixc = [k for k in range(n) if sig[k] == k]
        Sg = np.zeros((n, n), dtype=np.int64)
        Sg[sig, np.arange(n)] = 1
        if len(fixc) != 1 or not np.array_equal(Sg @ W @ Sg.T, W.T):
            continue
        lab, x = [None] * n, 1
        lab[fixc[0]] = 0
        for k in range(n):
            if lab[k] is None:
                lab[k], lab[sig[k]], x = x, (-x) % n, x + 1
        Qp = np.zeros((n, n), dtype=np.int64)
        Qp[lab, np.arange(n)] = 1
        Bc = Qp @ W @ Qp.T
        break
    else:
        raise ValueError("no fixed point of P gives a core with the reversal as anti-automorphism")
    I = np.eye(n, dtype=np.int64)
    R = I[(-np.arange(n)) % n]
    if not (np.array_equal(R @ Bc @ R, Bc.T) and np.array_equal(Bc + Bc.T, 0 * I) and not (Bc @ np.ones(n)).any()
            and np.array_equal(Bc @ Bc.T, n * I - 1)):
        raise ValueError("the relabelled core is not a skew core with the reversal as anti-automorphism")
    q = cert["q"]
    if q != 2 * n - 1 or q % 4 != 1 or not _is_prime(q):
        raise ValueError("q must be the prime 2n - 1 = 8l - 3, = 1 (mod 4)")
    alpha, beta = _yamada_alpha_beta(q, cert["r"], tuple(cert["xi"]))
    M, Nc = circulant(beta), circulant(alpha)
    if M[0, 0] != 0 or not (np.array_equal(M, M.T) and np.array_equal(Nc, Nc.T)) or \
            not np.array_equal(M @ M + Nc @ Nc, q * I):
        raise ValueError("Yamada's circulants are not symmetric with M^2 + N^2 = qI")
    Zn = np.zeros((n, n), dtype=np.int64)
    return _lemma26_propus([I, Zn, Bc @ R, Bc @ R], [M, Nc, R, R]), None


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
    alpha, beta = _yamada_alpha_beta(q, r, (xa, xb))
    return _yamada_matrix(cert, q, alpha, beta)


def _yamada_alpha_beta(q, r, xi):
    """Yamada's circulant first rows alpha, beta of order (q + 1)/2 for a prime q (see build_yamada); r must be a
    non-square and xi = xa + xb w primitive in GF(q)[w]/(w^2 - r), both checked."""
    xa, xb = xi
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
    return alpha, beta


def _build_yamada_field(cert):
    """build_yamada for a prime power q = p^k = 5 (mod 8) (_yamada_field_alpha_beta)."""
    if cert["Q"]["order"] != (cert["q"] + 3) // 2:
        raise ValueError("Q must have order (q + 3)/2")
    q, alpha, beta = _yamada_field_alpha_beta(cert)
    return _yamada_matrix(cert, q, alpha, beta)


def _yamada_field_alpha_beta(cert):
    """Yamada's alpha, beta (build_yamada) for a prime power q = p^k = 5 (mod 8). Elements of GF(q) are integers
    0 .. q - 1 whose base-p digits (least significant first) are the coefficients of a polynomial in x modulo
    cert["poly"] (monic, given by its coefficients c_0 .. c_{k-1}: x^k = -(c_0 + c_1 x + ... + c_{k-1} x^{k-1})). r is
    a non-square of GF(q) and xi = (xa, xb) = xa + xb w a primitive element of GF(q^2) = GF(q)[w]/(w^2 - r), in the
    same encoding."""
    q, r, (xa, xb), poly = cert["q"], cert["r"], cert["xi"], cert["poly"]
    p, k = cert["p"], len(poly)
    if p ** k != q or q % 8 != 5 or any(p % d == 0 for d in range(2, math.isqrt(p) + 1)):
        raise ValueError("q must be a prime power p^k = 5 (mod 8)")
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
    return q, alpha, beta


def _yamada_inputs(cert, q, alpha, beta):
    """Checks alpha, beta (symmetric, alpha alpha^T + beta beta^T = qI) and returns the skew core S of the skew-Hadamard
    matrix of cert["Q"], normalised to [[1, e^T], [-e, S + I]] (S^T = -S, S J = 0, S S^T = nI - J, n = (q + 1)/2)."""
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
    return (Q * d[:, None] * d[None, :])[1:, 1:] - np.eye(n0, dtype=np.int8)


def _yamada_matrix(cert, q, alpha, beta):
    n0 = (q + 1) // 2
    S = _yamada_inputs(cert, q, alpha, beta)
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


def build_symmetric_yamada_kiyasu(cert):
    """Symmetric Hadamard matrix of order 4v, v = q + 2, for a prime power q = 5 (mod 8) and ANY skew-Hadamard matrix of
    order (q + 3)/2 (NOTES.md §31). n = (q + 1)/2; S is the normalised skew core and alpha, beta are Yamada's circulants,
    as in build_yamada (same fields "q", "r", "xi", "Q", and "p", "poly" for a prime power); M = circ(beta) (zero
    diagonal), N = circ(alpha), M^2 + N^2 = qI, G = ((M + N) + i(M - N))/2, P = [[1, -1], [-1, 1]], J_2 the all-one
    matrix of order 2, e the all-one row of length 2n. The quaternary matrix of order 2v (blocks of orders 1, 2n, 1, 2n)

        L = [[-1,    -e,                           -1,     i e                          ],
             [-e^T,  (1 + i)/2 I (x) J_2 + G (x) P,  -i e^T, -i I (x) P - S (x) J_2         ],
             [-1,    -i e,                          1,     e                            ],
             [i e^T, -i I (x) P + S (x) J_2,         e^T,   (i - 1)/2 I (x) J_2 + G* (x) P]]

    is symmetric with entries +-1, +-i, and L L* = 2vI: on the vectors x (x) (1, -1) it acts as
    2 [[G, -iI], [-iI, G*]] (G G* = qI/2 since MN = NM), on the borders and the x (x) (1, 1) as a matrix of order 2n + 2
    whose square norm needs only S^T = -S, S J = 0, S^2 = J - nI, and the two parts do not meet (J_2 P = 0, e P = 0).
    It is Yamada's Lemma 26 data U = (I, 0, S, -S), V = (M, N, I, I) (Theorem 4, after Kiyasu) put in quaternary form.
    With L = X + iY (X, Y symmetric, of disjoint supports, X^2 + Y^2 = 2vI, XY = YX),
    H = [[X + Y, X - Y], [X - Y, -X - Y]] is a symmetric Hadamard matrix. Everything is checked here."""
    if "poly" in cert:
        q, alpha, beta = _yamada_field_alpha_beta(cert)
    else:
        q = cert["q"]
        if q % 8 != 5 or not _is_prime(q):
            raise ValueError("q must be a prime = 5 (mod 8)")
        alpha, beta = _yamada_alpha_beta(q, cert["r"], tuple(cert["xi"]))
    if cert["Q"]["order"] != (q + 3) // 2:
        raise ValueError("Q must have order (q + 3)/2")
    S = _yamada_inputs(cert, q, alpha, beta).astype(np.int16)
    n = (q + 1) // 2
    M, N = circulant(beta).astype(np.int16), circulant(alpha).astype(np.int16)
    I, J2 = np.eye(n, dtype=np.int16), np.ones((2, 2), dtype=np.int16)
    P = np.array([[1, -1], [-1, 1]], dtype=np.int16)
    k = np.kron
    e, o, z = np.full((1, 2 * n), 2, np.int16), np.zeros((1, 1), np.int16), np.zeros((1, 2 * n), np.int16)
    c = lambda x: np.full((1, 1), x, np.int16)
    X2 = np.block([[c(-2), -e, c(-2), z],                     # 2 Re L
                   [-e.T, k(I, J2) + k(M + N, P), z.T, -2 * k(S, J2)],
                   [c(-2), z, c(2), e],
                   [z.T, 2 * k(S, J2), e.T, -k(I, J2) + k(M + N, P)]])
    Y2 = np.block([[o, z, o, e],                              # 2 Im L
                   [z.T, k(I, J2) + k(M - N, P), -e.T, -2 * k(I, P)],
                   [o, -e, o, z],
                   [e.T, -2 * k(I, P), z.T, k(I, J2) - k(M - N, P)]])
    if (X2 % 2).any() or (Y2 % 2).any():
        raise ValueError("L has non-integral parts")
    X, Y = (X2 // 2).astype(np.int8), (Y2 // 2).astype(np.int8)
    del X2, Y2
    if not (np.abs(X) + np.abs(Y) == 1).all() or not (np.array_equal(X, X.T) and np.array_equal(Y, Y.T)):
        raise ValueError("L is not a symmetric matrix with entries +-1, +-i")
    H = np.block([[X + Y, X - Y], [X - Y, -X - Y]])
    if not np.array_equal(H, H.T):
        raise ValueError("H is not symmetric")
    return H, None


def build_paley_skew(cert):
    """Paley's skew-Hadamard matrix I + C of order q + 1, q = 3 (mod 4) a prime ("q") or a prime power ({"p", "poly"}
    under "q", as in _paley_sym_od_field); C[0, x] = 1, C[x, 0] = -1, C[x, y] = chi(y - x)."""
    qs = cert["q"]
    P, CP = _paley_sym_od(qs) if isinstance(qs, int) else _paley_sym_od_field(qs["p"], qs["poly"])
    C = CP @ P
    m = C.shape[0]
    H = (np.eye(m, dtype=np.int64) + C).astype(np.int8)
    if not is_skew(H):
        raise ValueError("I + C is not skew")
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


def build_skew_amicable_product(cert):
    """Skew-Hadamard matrix of order h m by J. Wallis's product (J. Combin. Theory Ser. A 21 (1976); Seberry-Yamada
    1992, SXI). If I + S is a skew-Hadamard matrix of order h and (A, M) an amicable pair of order m (A skew-Hadamard,
    M symmetric Hadamard, A M^T = M A^T), then K = I (x) A + S (x) M is a skew-Hadamard matrix of order h m:
    K + K^T = I (x) (A + A^T) + (S + S^T) (x) M = 2I, and
    K K^T = I (x) A A^T + (S + S^T) (x) A M^T + S S^T (x) M M^T = m I + (h - 1) m I.
    "factor": the certificate fields of build_goethals_seidel for a skew-Hadamard matrix of order h = 4v (first
    block of skew type; checked). "amicable": {"n", "blocks", "rowsums"}, a skew-propus quadruple of length n; the
    pair A = I + PN, M = P + N of order m = 4n is that of build_skew_amicable_plug (checked)."""
    F, fseqs = build_goethals_seidel(cert["factor"])
    h = F.shape[0]
    Ih = np.eye(h, dtype=np.int64)
    if not (np.array_equal(F + F.T, 2 * Ih) and np.array_equal(F @ F.T, h * Ih)):
        raise ValueError("the factor is not a skew-Hadamard matrix")
    am = cert["amicable"]
    P, N, aseqs = _skew_propus_od(am["n"], am["blocks"], am["rowsums"])
    m = P.shape[0]
    I = np.eye(m, dtype=np.int64)
    A, M = I + P @ N, P + N
    if not (np.array_equal(A + A.T, 2 * I) and np.array_equal(M, M.T) and np.array_equal(A @ A.T, m * I)
            and np.array_equal(M @ M, m * I) and np.array_equal(A @ M.T, M @ A.T)):
        raise ValueError("A = I + PN and M = P + N are not an amicable pair")
    K = (np.kron(np.eye(h, dtype=np.int8), A.astype(np.int8))
         + np.kron((F - Ih).astype(np.int8), M.astype(np.int8)))
    return K, fseqs + aseqs


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


def _pgp_s_form(cert, n):
    """The Hadamard matrix K = [[P, Q], [Q, P]] of order 2n of a periodic Golay pair (a, b) over an abelian group
    Gamma of order n with an element t of order 2: K' = [[A, B], [-B^T, A^T]] (A, B group-developed) is Hadamard, the
    translation by t (on rows and columns of both halves) is a fixed-point-free involutory automorphism, and ordering
    the indices as (coset c, half r, transversal element u), with element u + c t, makes it the block swap."""
    a = np.array(cert["a"], dtype=np.int64)
    b = np.array(cert["b"], dtype=np.int64)
    tg, add, neg = abelian_group(cert.get("group", str(len(a))))
    t = int(cert["t"])
    if tg != n or len(a) != n or len(b) != n or t == 0 or add[t, t] != 0:
        raise ValueError("periodic Golay pair: wrong length or t is not of order 2")
    if not (np.abs(a) == 1).all() or not (np.abs(b) == 1).all():
        raise ValueError("periodic Golay pair: entries must be +-1")
    paf = np.array([int(a @ a[add[:, s]] + b @ b[add[:, s]]) for s in range(n)])
    if paf[0] != 2 * n or np.any(paf[1:] != 0):
        raise ValueError("(a, b) is not a periodic Golay pair")
    A, B = group_matrix(add, neg, a), group_matrix(add, neg, b)
    K = np.block([[A, B], [-B.T, A.T]])
    U = [g for g in range(n) if g < add[g, t]]
    order = [r * n + (u if c == 0 else int(add[u, t])) for c in (0, 1) for r in (0, 1) for u in U]
    K = K[np.ix_(order, order)]
    if not np.array_equal(K[n:, n:], K[:n, :n]) or not np.array_equal(K[n:, :n], K[:n, n:]):
        raise ValueError("K is not of the form [[P, Q], [Q, P]]")
    if not gram_is_scalar(K, 2 * n):
        raise ValueError("K is not a Hadamard matrix")
    return K, [a, b]


def build_symmetric_involution(cert):
    """Symmetric Hadamard matrix of order 4q, q = 1 (mod 4) a prime, from a Hadamard matrix K = [[P, Q], [Q, P]] of
    order q - 1 (symmetric/THEOREM_F.md), here from a periodic Golay pair of length (q - 1)/2 (_pgp_s_form).

    Points of the projective line: x = 0 .. q - 1 and infinity; chi is the quadratic character. On the q - 1 points
    x != 0 put A[x, y] = chi(x + y) (zero exactly at y = -x) and Z[x, y] = [x + y = 0]. Row r < (q - 1)/2 of K is the
    point r + 1 and row (q - 1)/2 + r the point -(r + 1), so that x -> -x acts on the rows of K as the block swap; Kp
    is K with its rows so relabelled. Then

        F = [[A + Z, A - Z], [A - Z, A + Z]],   B = [ (1; 1) | (chi; chi) | (k_j; -k_j) ],   P = -B^T F B / (4(q - 1)),

    where k_j are the columns of Kp, satisfy B^T B = 2(q - 1) I, F B = -2 B P and F^2 + 2 B B^T = 4q I, and P is a
    symmetric signed permutation matrix without fixed points. For a non-square nu, V0[x, y] = chi(xy - nu),
    V0[x, inf] = chi(x), V0[inf, inf] = 1 is a symmetric weighing matrix of order q + 1 and weight q whose zeros are
    the fixed-point-free involution x -> nu/x (0 <-> inf); V is V0 relabelled so that its zeros are those of P. Then

        H = [[F, B, B], [B^T, P + V, P - V], [B^T, P - V, P + V]]

    is symmetric with entries +-1 and H^2 = 4q I. "q", "nu" (optional, default the least non-square) and the pair
    ("a", "b", "group", "t") as in _pgp_s_form. Every step is checked here."""
    q = int(cert["q"])
    if not _is_prime(q) or q % 4 != 1:
        raise ValueError("symmetric-involution: q must be a prime = 1 (mod 4)")
    m, n = q - 1, (q - 1) // 2
    K, seqs = _pgp_s_form(cert, n)
    chi = -np.ones(q, dtype=np.int64)
    chi[np.arange(1, q) ** 2 % q] = 1
    chi[0] = 0
    nu = int(cert.get("nu", next(x for x in range(2, q) if chi[x] == -1)))
    if chi[nu % q] != -1:
        raise ValueError("nu must be a non-square")
    x = np.arange(1, q)
    A = chi[(x[:, None] + x[None, :]) % q]
    Z = ((x[:, None] + x[None, :]) % q == 0).astype(np.int64)
    point = np.concatenate([np.arange(1, n + 1), q - np.arange(1, n + 1)])
    Kp = np.zeros_like(K)
    Kp[point - 1] = K
    F = np.block([[A + Z, A - Z], [A - Z, A + Z]])
    B = np.hstack([np.ones((2 * m, 1), dtype=np.int64), np.concatenate([chi[x], chi[x]])[:, None],
                   np.vstack([Kp, -Kp])])
    if not (np.abs(F) == 1).all() or not np.array_equal(F, F.T):
        raise ValueError("F is not a symmetric +-1 matrix")
    if not np.array_equal(B.T @ B, 2 * m * np.eye(q + 1, dtype=np.int64)):
        raise ValueError("B^T B != 2(q - 1) I")
    G = -(B.T @ F @ B)
    if np.any(G % (4 * m)):
        raise ValueError("B^T F B is not divisible by 4(q - 1)")
    P = G // (4 * m)
    if (not np.array_equal(P, P.T) or not np.array_equal(np.abs(P).sum(axis=1), np.ones(q + 1, dtype=np.int64))
            or np.any(np.diag(P) != 0)):
        raise ValueError("P is not a fixed-point-free symmetric signed permutation matrix")
    if not np.array_equal(F @ B, -2 * B @ P) or not np.array_equal(F @ F + 2 * B @ B.T,
                                                                    4 * q * np.eye(2 * m, dtype=np.int64)):
        raise ValueError("F B != -2 B P or F^2 + 2 B B^T != 4q I")
    V0 = np.zeros((q + 1, q + 1), dtype=np.int64)
    f = np.arange(q)
    V0[:q, :q] = chi[(f[:, None] * f[None, :] - nu) % q]
    V0[:q, q] = V0[q, :q] = chi[f]
    V0[q, q] = 1
    if not np.array_equal(V0, V0.T) or not np.array_equal(V0 @ V0, q * np.eye(q + 1, dtype=np.int64)):
        raise ValueError("V0 is not a symmetric weighing matrix of weight q")
    def pairs(M):
        z = [int(np.flatnonzero(row)[0]) for row in M]
        if any(z[z[i]] != i or z[i] == i for i in range(len(z))):
            raise ValueError("the zero pattern is not a fixed-point-free involution")
        return sorted((i, z[i]) for i in range(len(z)) if i < z[i])
    tau = np.zeros(q + 1, dtype=np.int64)
    for (i, j), (u, w) in zip(pairs(np.abs(P)), pairs((V0 == 0).astype(np.int64))):
        tau[i], tau[j] = u, w
    V = V0[np.ix_(tau, tau)]
    if not np.array_equal(V == 0, P != 0):
        raise ValueError("the zeros of V are not the support of P")
    H = np.block([[F, B, B], [B.T, P + V, P - V], [B.T, P - V, P + V]])
    return H, seqs


def build_paley_blocks(cert):
    """Paley's matrix of order q + 1 = 4s (q = 3 mod 4, q > 3; a prime, or a prime power with cert["field"], see
    _build_paley_blocks_field) as a 4 x 4 array of circulant s x s blocks, times
    T-matrices of order t (Seberry-Yamada 1992, Theorem 3.10), as a Goethals-Seidel quadruple of length n = st.
    Points of PG(1, q) are 0..q-1 and q (infinity), with vectors u_x = (1, x) and u_inf = (0, -1); Paley's matrix is
    P[a, b] = [a = b] + chi(det(u_a, u_b)).  The torus element g is (x -> -1/(x + tau))^power; P(ga, gb) = e(a) e(b)
    P(a, b); E(g a) = c e(a) E(a) along each g-cycle (cycles in order of their least point, each from its least point)
    makes g an automorphism of E P E, whose blocks along the cycles are circulant.  A_j(x) = sum_k T_k(x mod t)
    N_kj(x mod s); the certificate's sequences must equal the rebuilt ones.  Optional fields: "shifts" r moves the
    base point of cycle k by r_k (N_kj(y) -> N_kj(y + r_j - r_k)), and "columns" lists the block columns in the order
    of the sequences; skew certificates (3(q + 1), see search/paley_blocks.py) use both."""
    q, t, tau, power, c = cert["q"], cert["t"], cert["tau"], cert["power"], cert["c"]
    if q % 4 != 3 or q < 7:
        raise ValueError("q must be = 3 (mod 4) and at least 7")
    if not _is_prime(q):
        return _build_paley_blocks_field(cert)
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
    """the common part of build_paley_blocks: cycles of g, signs, circulant blocks, T-matrices, sequences, array.
    For q = 7 (mod 8), s is even and g^s acts as -I: every cycle closes with the sign -1, the blocks are negacyclic
    (N_kj[i, j] = N~_kj(j - i) with N~(y + s) = -N~(y)), the certificate says "negacyclic": true and
    "group": "Z_t x N_s", and the Goethals-Seidel array uses R: (i, a) -> (-i, s - 1 - a), so that R X R = X^T."""
    q, t, c = cert["q"], cert["t"], cert["c"]
    s = (q + 1) // 4
    nega = s % 2 == 0
    if bool(cert.get("negacyclic")) != nega:
        raise ValueError("the blocks are negacyclic exactly when s is even")
    grouped = nega or math.gcd(s, t) != 1
    group = f"Z{t} x N{s}" if nega else f"Z{t} x Z{s}"
    if s * t != cert["v"] or grouped != (cert.get("group") == group):
        raise ValueError("v must be s t, with the group (Z_t x Z_s, or Z_t x N_s if s is even) recorded exactly when "
                         "gcd(s, t) > 1 or s is even")
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
        if c * e[cy[-1]] * E[cy[-1]] != (-1 if nega else 1):
            raise ValueError("sign normalisation does not close")
    order = np.array([a for cy in cyc for a in cy])
    E8 = E.astype(np.int8)
    B = (E8[order, None] * P[np.ix_(order, order)] * E8[None, order]).astype(np.int8)
    N = [[B[k * s, j * s:(j + 1) * s].astype(np.int64) for j in range(4)] for k in range(4)]
    ia = np.arange(s)
    wrap = np.where(ia[None, :] >= ia[:, None], 1, -1) if nega else 1      # negacyclic: sign -1 below the diagonal
    ext = lambda row, y: (np.where((y // s) % 2 == 1, -1, 1) if nega else 1) * row[y % s]   # N~(y), y any integer
    for k in range(4):
        for j in range(4):
            blk = wrap * N[k][j][(ia[None, :] - ia[:, None]) % s]
            if not np.array_equal(B[k * s:(k + 1) * s, j * s:(j + 1) * s], blk):
                raise ValueError("blocks are not circulant" + (" (negacyclic)" if nega else ""))
    r, cols = cert.get("shifts", [0, 0, 0, 0]), cert.get("columns", [0, 1, 2, 3])
    if len(r) != 4 or not all(0 <= x < (2 * s if nega else s) for x in r) or sorted(cols) != [0, 1, 2, 3]:
        raise ValueError("shifts must be four residues mod s (2s if negacyclic) and columns a permutation of 0..3")
    N = [[ext(N[k][j], ia + r[j] - r[k]) for j in cols] for k in range(4)]   # new base points: N_kj(y + r_j - r_k)
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
    if grouped:                                            # sum_x a(x) a~(x + d), a~ antiperiodic in a if negacyclic
        sh = lambda a, da: np.concatenate([a, -a if nega else a], axis=1)[:, da:da + s]
        paf = sum(np.array([[int((a * np.roll(sh(a, da), -di, 0)).sum()) for da in range(s)]
                            for di in range(t)]) for a in seq2).reshape(-1)
    else:
        paf = paf_sum(seqs)
    if paf[0] != 4 * s * t or np.any(paf[1:] != 0):
        raise ValueError("sequences are not complementary")
    if not grouped:
        return goethals_seidel(*map(circulant8, seqs)), seqs
    it, ia = np.arange(t), np.arange(s)
    dt, ds = (it[None, :] - it[:, None]) % t, (ia[None, :] - ia[:, None]) % s
    sg = (wrap if nega else np.ones((s, s), dtype=np.int64)).astype(np.int8)[None, :, None, :]
    gm = [(a.astype(np.int8)[dt[:, None, :, None], ds[None, :, None, :]] * sg).reshape(s * t, s * t) for a in seq2]
    neg = (((-it) % t)[:, None] * s + ((s - 1 - ia) if nega else (-ia) % s)[None, :]).reshape(-1)   # X R = X[:, neg]
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


def _paley_conference_field(qs):
    """Paley's symmetric conference matrix of order q + 1 for a prime power q = 1 (mod 4): "q" is a prime, or
    {"p", "poly"} for q = p^k with elements of GF(q) encoded as in _paley_sym_od_field (poly primitive, checked).
    Index 0 is the point at infinity and index 1 + x is x: C[0, 0] = 0, C[0, j] = C[j, 0] = 1, C[1 + x, 1 + y] =
    chi(y - x)."""
    if isinstance(qs, int):
        if not (_is_prime(qs) and qs % 4 == 1):
            raise ValueError("q must be a prime = 1 (mod 4)")
        return paley_conference(qs).astype(np.int64)
    p, poly = qs["p"], qs["poly"]
    k = len(poly)
    q = p ** k
    if not (_is_prime(p) and q % 4 == 1):
        raise ValueError("q must be a prime power = 1 (mod 4)")
    digits = np.array([[(e // p ** j) % p for j in range(k)] for e in range(q)], dtype=np.int64)
    weights = p ** np.arange(k)
    enc = lambda d: int((np.asarray(d) % p) @ weights)
    log = -np.ones(q, dtype=np.int64)
    e = 1
    for i in range(q - 1):
        if log[e] >= 0:
            raise ValueError("poly is not primitive")
        log[e] = i
        d = np.zeros(k + 1, dtype=np.int64); d[1:] = digits[e]
        e = enc(d[:k] - d[k] * np.asarray(poly, dtype=np.int64))
    if e != 1:
        raise ValueError("poly is not primitive")
    chi = np.where(log % 2 == 0, 1, -1).astype(np.int64)
    chi[0] = 0
    sub = np.array([[enc(digits[y] - digits[x]) for y in range(q)] for x in range(q)], dtype=np.int64)   # y - x
    C = np.ones((q + 1, q + 1), dtype=np.int64)
    C[0, 0] = 0
    C[1:, 1:] = chi[sub]
    return C


def build_symmetric_conference_plug(cert):
    """Symmetric Hadamard matrix of order 2n(n - 1) from a symmetric conference matrix C of order n
    (symmetric/CONFERENCE_PLUG.md). C = [[0, e^T], [e, Q]] is Paley's matrix of order n = q + 1 for the prime power
    q = 1 (mod 4) under "q" (_paley_conference_field). With m = n - 1 and J = J_m, X = [[J, J], [J, -J]] and
    Y = [[Q + I, Q - I], [Q - I, -Q - I]] are symmetric and anticommute (Q 1 = 0), and H = I_n (x) X + C (x) Y.
    Checked here: C = C^T, zero diagonal, C^2 = mI, Q 1 = 0 and XY + YX = 0."""
    C = _paley_conference_field(cert["q"])
    n = C.shape[0]
    m = n - 1
    if (not np.array_equal(C, C.T) or np.any(np.diag(C) != 0) or np.any(C[0, 1:] != 1)
            or not np.array_equal(C @ C, m * np.eye(n, dtype=np.int64))):
        raise ValueError("C is not a normalised symmetric conference matrix")
    Q = C[1:, 1:]
    if np.any(Q.sum(axis=1) != 0):
        raise ValueError("Q 1 != 0")
    I = np.eye(m, dtype=np.int64)
    J = np.ones((m, m), dtype=np.int64)
    X = np.block([[J, J], [J, -J]])
    Y = np.block([[Q + I, Q - I], [Q - I, -Q - I]])
    if np.any(X @ Y + Y @ X != 0):
        raise ValueError("X and Y do not anticommute")
    H = np.kron(np.eye(n, dtype=np.int8), X.astype(np.int8)) + np.kron(C.astype(np.int8), Y.astype(np.int8))
    return H, None


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
            "skew-amicable-product": build_skew_amicable_product,
            "symmetric-propus-miyamoto": build_symmetric_propus_miyamoto,
            "symmetric-propus-sb": build_symmetric_propus_sb,
            "symmetric-propus-skewcore": build_symmetric_propus_skewcore,
            "symmetric-propus-weighing": build_symmetric_propus_weighing,
            "symmetric-yamada-kiyasu": build_symmetric_yamada_kiyasu, "paley-skew": build_paley_skew,
            "two-block": build_two_block, "paley-blocks": build_paley_blocks,
            "symmetric-involution": build_symmetric_involution,
            "symmetric-conference-plug": build_symmetric_conference_plug}


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
