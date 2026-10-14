"""Symmetric Hadamard matrices of order m (m - 1) from symmetric orthogonal designs OD(m; 1, m - 1).

A symmetric OD(m; 1, m - 1) is Y = x P + y N where P is a symmetric signed permutation matrix (P^2 = I), N is
symmetric with N[i, j] = 0 exactly where P[i, j] != 0, N^2 = (m - 1) I and P N + N P = 0. Equivalently,
S = I + N P is a skew-Hadamard matrix and P an involutive signed permutation with P C P = -C for C = N P
(an anti-automorphism of S of order 2). Then M = P + N is a symmetric Hadamard matrix of order m. Normalise M
symmetrically (E M E with E diagonal +-1, after moving a diagonal entry -1 to position 0) so that row and
column 0 are (-1, 1, ..., 1); the rest, D, is symmetric with D J = J and D^2 = m I - J, and

    H = P (x) J_{m-1} + N (x) D

is a symmetric Hadamard matrix of order m (m - 1), because J and D commute:
    H H^T = P^2 (x) J^2 + (P N + N P) (x) J + N^2 (x) D^2 = (m - 1) I (x) J + (m - 1) I (x) (m I - J) = m (m - 1) I.
For prime powers q = 3 (mod 8) and the design built from Xia-Xia-Seberry-Wu difference families this is
Dokovic's theorem (arXiv:2508.17141, Theorem 2, where D comes from the Paley core); any symmetric Hadamard matrix
of order m can supply D. The design exists e.g. for

  * Paley: m = q + 1, q = 3 (mod 4) a prime (prime power), P = (-1) + (x -> -x);
  * skew-propus quadruples of odd length n (m = 4n): A3 symmetric, A1 used twice, A0 of skew type with
    PAF_A3 + 2 PAF_A1 + PAF_A0 = 4n delta, put into the Balonin array with D = x I + y (A0 - I) (Dokovic's Y).
    Since the row sum of A0 is 1, a3^2 + 2 a1^2 = 4n - 1 is necessary (it fails for n = 23: 91 = 7 * 13).

  python3 search/sym_mm1.py test                    (Paley q = 3, 7, 11, 19, 23 and every skew-propus length <= 15)
  python3 search/sym_mm1.py propus n [--all]        (exhaustive skew-propus search of length n; first solution
                                                     or all, up to multipliers, shifts of A1 and A0 -> A0^T)
  python3 search/sym_mm1.py mpropus n gens [--all]  (the same, constant on the orbits of the multiplier group <gens>)
  python3 search/sym_mm1.py mszek v gens [--all]    (Legendre pairs (skew type, symmetric) of length v on the orbits
                                                     of <gens>: the Szekeres-type design of order 2v + 2)
  python3 search/sym_mm1.py build -- n A3 A1 A0 [out]    (build and check the matrix of order 4n (4n - 1))
  python3 search/sym_mm1.py certify -- n A3 A1 A0        (write symmetric/<4n(4n-1)>/: certificate.json
                                                          (construction "skew-propus-plug"), sequences.txt, H.csv.gz)
"""
import argparse, gzip, hashlib, itertools, json, math, sys, time
import numpy as np
from pathlib import Path as _P
ROOT = _P(__file__).resolve().parent.parent
sys.path[:0] = [str(ROOT / 'verify'), str(ROOT / 'search')]
import verify

SIGN = {'+': 1, '-': -1}


def check_od(P, N):
    """True if x P + y N is a symmetric OD(m; 1, m - 1)."""
    m = P.shape[0]
    I = np.eye(m, dtype=np.int64)
    P, N = P.astype(np.int64), N.astype(np.int64)
    if not (np.array_equal(P, P.T) and np.array_equal(N, N.T)):
        return False
    if not (np.array_equal(np.abs(P).sum(1), np.ones(m, np.int64)) and np.array_equal(P @ P, I)):
        return False
    if not np.array_equal(np.abs(P) + np.abs(N), np.ones((m, m), np.int64)):
        return False
    return np.array_equal(N @ N, (m - 1) * I) and not (P @ N + N @ P).any()


def core(M):
    """D of order m - 1 with D symmetric, D J = J, D^2 = m I - J, from a symmetric Hadamard matrix M of order m."""
    M = M.astype(np.int64)
    r = int(np.flatnonzero(np.diag(M) == -1)[0]) if (np.diag(M) == -1).any() else None
    if r is None:
        M, r = -M, 0
    order = [r] + [i for i in range(M.shape[0]) if i != r]
    M = M[np.ix_(order, order)]
    e = M[0].copy()
    e[0] = 1
    M = M * e[:, None] * e[None, :]
    assert M[0, 0] == -1 and (M[0, 1:] == 1).all()
    return M[1:, 1:]


def plug(P, N, D=None):
    """H = P (x) J + N (x) D, order m (m - 1)."""
    m = P.shape[0]
    if D is None:
        D = core(P + N)
    J = np.ones((m - 1, m - 1), dtype=np.int8)
    return (np.kron(P.astype(np.int8), J) + np.kron(N.astype(np.int8), D.astype(np.int8))).astype(np.int8)


def plug3(P, N, M):
    """H = P (x) (J - 2I) + N (x) D, D the core of a symmetric Hadamard matrix M of order m + 4: order m (m + 3).
    (J - 2I)^2 = 4I + (m - 1) J and (m - 1) D^2 = (m - 1)(m + 4) I - (m - 1) J add up to m (m + 3) I."""
    m = P.shape[0]
    D = core(M)
    r = D.shape[0]
    assert r == m + 3
    X = np.ones((r, r), dtype=np.int8) - 2 * np.eye(r, dtype=np.int8)
    return (np.kron(P.astype(np.int8), X) + np.kron(N.astype(np.int8), D.astype(np.int8))).astype(np.int8)


def paley2(q, doublings=0):
    """Symmetric Hadamard matrix [[C + I, C - I], [C - I, -C - I]] of order 2(q + 1) (q = 1 mod 4 prime), doubled."""
    C = verify.paley_conference(q).astype(np.int64)
    I = np.eye(q + 1, dtype=np.int64)
    M = np.block([[C + I, C - I], [C - I, -C - I]])
    for _ in range(doublings):
        M = np.kron(np.array([[1, 1], [1, -1]], dtype=np.int64), M)
    return M


def is_symmetric_hadamard(H):
    n = H.shape[0]
    return np.array_equal(H, H.T) and verify.is_pm1(H) and verify.gram_is_scalar(H, n)


def paley(q):
    """Symmetric OD(q + 1; 1, q) from Paley's skew-Hadamard matrix, q = 3 (mod 4) prime."""
    chi = -np.ones(q, dtype=np.int64)
    chi[np.arange(1, q) ** 2 % q] = 1
    chi[0] = 0
    m = q + 1
    C = np.zeros((m, m), dtype=np.int64)
    C[0, 1:] = 1
    C[1:, 0] = -1
    C[1:, 1:] = verify.circulant(chi)
    P = np.zeros((m, m), dtype=np.int64)
    P[0, 0] = -1
    for x in range(q):
        P[1 + x, 1 + (-x) % q] = 1
    return P, C @ P


def balonin(A3, A1, A0):
    """Dokovic's symmetric OD(4n; 1, 4n - 1) in the Balonin array:
        [ -A   BR   BR   DR ]
        [ BR   RD   A   -RB ]
        [ BR   A   -RD   RB ]
        [ DR  -RB   RB   A  ]       with A = y A3, B = y A1, D = x I + y (A0 - I)."""
    n = len(A3)
    R = verify.back_identity(n)
    A, B = verify.circulant(A3), verify.circulant(A1)
    D0 = verify.circulant(A0) - np.eye(n, dtype=np.int64)
    Z = np.zeros((n, n), dtype=np.int64)
    BR, RB = B @ R, R @ B
    N = np.block([[-A, BR, BR, D0 @ R], [BR, R @ D0, A, -RB], [BR, A, -R @ D0, RB], [D0 @ R, -RB, RB, A]])
    P = np.block([[Z, Z, Z, R], [Z, R, Z, Z], [Z, Z, -R, Z], [R, Z, Z, Z]])
    return P, N


def szekeres(a, b):
    """Symmetric OD(2v + 2; 1, 2v + 1) from a Legendre pair on Z_v with a of skew type and b symmetric (row sums 1):
    S = [[1, 1, e, -e], [-1, 1, -e, -e], [-e, e, A, B], [e, e, -B, A^T]] is skew-Hadamard and
    P = [[0, 1], [1, 0]] + diag(R, -R) (R: g -> -g) satisfies P C P = -C; returns (P, N) with N = C P."""
    v = len(a)
    A, B = verify.circulant(a), verify.circulant(b)
    e = np.ones((v, 1), dtype=np.int64)
    m = 2 * v + 2
    S = np.zeros((m, m), dtype=np.int64)
    S[:2, :2] = [[1, 1], [-1, 1]]
    S[0, 2:2 + v], S[0, 2 + v:] = 1, -1
    S[1, 2:2 + v], S[1, 2 + v:] = -1, -1
    S[2:2 + v, :2] = np.hstack([-e, e])
    S[2 + v:, :2] = np.hstack([e, e])
    S[2:2 + v, 2:2 + v], S[2:2 + v, 2 + v:] = A, B
    S[2 + v:, 2:2 + v], S[2 + v:, 2 + v:] = -B, A.T
    C = S - np.eye(m, dtype=np.int64)
    R = np.zeros((v, v), dtype=np.int64)
    R[np.arange(v), (-np.arange(v)) % v] = 1
    P = np.zeros((m, m), dtype=np.int64)
    P[0, 1] = P[1, 0] = 1
    P[2:2 + v, 2:2 + v] = R
    P[2 + v:, 2 + v:] = -R
    return P, C @ P


def amicable_cores(P, N):
    """Amicable cores from a symmetric OD xP + yN of order m: S = I + PN (skew-Hadamard) and M = P + N (symmetric
    Hadamard) are amicable (S M^T = M S^T). Normalise both at the first fixed point r of P with P[r, r] = +1 by
    E = diag(S[r]): row r of E S E and of E M E is all ones; then E S E = [[1, e], [-e, I + W]] and
    E M E = [[1, e], [e, D]] (r moved first) with W skew, W e = 0, (I + W)(I + W)^T = mI - J, D symmetric, D e = -e,
    D^2 = mI - J and W D + D W = 0."""
    m = P.shape[0]
    S_ = np.eye(m, dtype=np.int64) + P @ N
    M = P + N
    r = next(i for i in range(m) if P[i, i] == 1)
    e = S_[r].copy()
    S1 = S_ * e[:, None] * e[None, :]
    M1 = M * e[:, None] * e[None, :]
    order = [r] + [i for i in range(m) if i != r]
    S1, M1 = S1[np.ix_(order, order)], M1[np.ix_(order, order)]
    assert (S1[0] == 1).all() and (M1[0] == 1).all()
    W = S1[1:, 1:] - np.eye(m - 1, dtype=np.int64)
    D = M1[1:, 1:]
    return W, D


def g_matrices(n, log=print):
    """Four +-1 sequences of odd length n, k0 and k1 of skew type, g2 and g3 symmetric, with
    PAF_k0 + PAF_k1 + PAF_g2 + PAF_g3 = 4n delta (they give a skew-type OD(4n; 1, 1, 4n - 2) in the
    Goethals-Seidel array, with a and b on the diagonals of k0 and k1). First solution found, or None."""
    h = (n - 1) // 2
    bits = np.array(list(itertools.product([1, -1], repeat=h)), dtype=np.int8)
    K = np.ones((len(bits), n), dtype=np.int8)
    K[:, 1:h + 1] = bits
    K[:, n - h:] = -bits[:, ::-1]
    sb = np.array(list(itertools.product([1, -1], repeat=h + 1)), dtype=np.int8)
    G = np.empty((len(sb), n), dtype=np.int8)
    G[:, 0] = sb[:, 0]
    G[:, 1:h + 1] = sb[:, 1:]
    G[:, n - h:] = sb[:, 1:][:, ::-1]
    PK, PG = pafs(K, h), pafs(G, h)
    rs = G.sum(1, dtype=np.int64)
    pair = {}
    for i in range(len(K)):
        for j in range(i, len(K)):
            pair.setdefault((PK[i] + PK[j]).tobytes(), (i, j))
    for s2 in range(1, 2 * n, 2):
        for s3 in range(s2, 2 * n, 2):
            if 2 + s2 * s2 + s3 * s3 != 4 * n:
                continue
            I2, I3 = np.flatnonzero(rs == s2), np.flatnonzero(rs == s3)
            log(f'g-matrices n = {n}: row sums 1, 1, {s2}, {s3}: {len(I2)} x {len(I3)} symmetric pairs')
            for a in I2:
                for b in I3:
                    key = (-(PG[a] + PG[b])).astype(PK.dtype).tobytes()
                    if key in pair:
                        i, j = pair[key]
                        return K[i].astype(np.int64), K[j].astype(np.int64), G[a].astype(np.int64), G[b].astype(np.int64)
    return None


def skew_od_from_g(k0, k1, g2, g3):
    """Coefficient matrices (B, C) of the skew-type OD aI + bB + cC of order 4n built from G-matrices in the
    Goethals-Seidel array (blocks aI + c(K0 - I), bI + c(K1 - I), c G2, c G3)."""
    n = len(k0)
    I = np.eye(n, dtype=np.int64)
    R = verify.back_identity(n)
    Z = np.zeros((n, n), dtype=np.int64)
    def gs(A, B, C, D):
        return np.block([[A, B @ R, C @ R, D @ R], [-B @ R, A, D.T @ R, -C.T @ R],
                         [-C @ R, -D.T @ R, A, B.T @ R], [-D @ R, C.T @ R, -B.T @ R, A]])
    Bm = gs(Z, I, Z, Z)
    Cm = gs(verify.circulant(k0) - I, verify.circulant(k1) - I, verify.circulant(g2), verify.circulant(g3))
    return Bm, Cm


def skew_od_from_codes(blocks):
    """Coefficient matrices (A, B, C) of the orthogonal design xA + yB + zC of order 4t given by four circulant
    blocks in the Goethals-Seidel array, each block a length-t list of codes 0, +-1 (x), +-2 (y), +-3 (z)."""
    t = len(blocks[0])
    R = verify.back_identity(t)
    def gs(A, B, C, D):
        return np.block([[A, B @ R, C @ R, D @ R], [-B @ R, A, D.T @ R, -C.T @ R],
                         [-C @ R, -D.T @ R, A, B.T @ R], [-D @ R, C.T @ R, -B.T @ R, A]])
    out = []
    for v in (1, 2, 3):
        parts = [verify.circulant([int(np.sign(c)) if abs(c) == v else 0 for c in b]) for b in blocks]
        out.append(gs(*parts))
    return tuple(out)


def check_skew_od(A, B, C, weights):
    """xI + yB + zC is a skew-type OD of type (1, b, c): A = I, B and C skew, BB^T = bI, CC^T = cI,
    BC^T + CB^T = 0, supports partitioning the matrix."""
    n = A.shape[0]
    b, c = weights
    I = np.eye(n, dtype=np.int64)
    ok = (np.array_equal(A, I) and np.array_equal(B.T, -B) and np.array_equal(C.T, -C)
          and np.array_equal(B @ B.T, b * I) and np.array_equal(C @ C.T, c * I) and not (B @ C.T + C @ B.T).any()
          and np.array_equal(np.abs(A) + np.abs(B) + np.abs(C), np.ones((n, n), dtype=np.int64)))
    return bool(ok)


def skew_seberry(P, N, Bm, Cm):
    """Seberry's 1978 plug-in (Ars Combin. 6, Theorem 7): H = I (x) (I + W) + B (x) J + C (x) D, a skew-Hadamard
    matrix of order N(m' - 1) from the amicable cores (W, D) of order m' - 1 (m' = order of P) and a skew-type
    OD(N; 1, k, N - k - 1) aI + bB + cC with N = k m' (k = 1: G-matrices, OD(m'; 1, 1, m' - 2))."""
    W, D = amicable_cores(P, N)
    r = P.shape[0] - 1
    I = np.eye(r, dtype=np.int8)
    J = np.ones((r, r), dtype=np.int8)
    return (np.kron(np.eye(Bm.shape[0], dtype=np.int8), I + W.astype(np.int8)) + np.kron(Bm.astype(np.int8), J)
            + np.kron(Cm.astype(np.int8), D.astype(np.int8))).astype(np.int8)


def pafs(X, half):
    """Periodic autocorrelations at shifts 1..half of the rows of X (int8 array, one sequence per row)."""
    X = X.astype(np.int32)
    return np.stack([(X * np.roll(X, -s, axis=1)).sum(1) for s in range(1, half + 1)], axis=1)


def skew_propus(n, find_all=False, log=print):
    """All (A3, A1, A0): A3 symmetric, A0 of skew type (a0[0] = 1, a0[-j] = -a0[j]), PAF_A3 + 2 PAF_A1 + PAF_A0 = 4n
    delta; row sums of A3 and A1 taken positive.  Solutions are reduced modulo multipliers (applied to all three),
    shifts of A1 and A0 -> A0^T."""
    h = (n - 1) // 2
    sols, seen = [], set()
    reps = [(a3, a1) for a3 in range(1, 2 * n, 2) for a1 in range(1, 2 * n, 2) if a3 * a3 + 2 * a1 * a1 == 4 * n - 1]
    log(f'n = {n}: row sums (a3, a1) with a3^2 + 2 a1^2 = {4 * n - 1}: {reps}')
    # skew-type A0: free signs at 1..h
    bits = np.array(list(itertools.product([1, -1], repeat=h)), dtype=np.int8)
    A0s = np.ones((len(bits), n), dtype=np.int8)
    A0s[:, 1:h + 1] = bits
    A0s[:, n - h:] = -bits[:, ::-1]
    P0 = pafs(A0s, h)
    # symmetric A3: a3[0] and signs at 1..h
    sb = np.array(list(itertools.product([1, -1], repeat=h + 1)), dtype=np.int8)
    A3all = np.empty((len(sb), n), dtype=np.int8)
    A3all[:, 0] = sb[:, 0]
    A3all[:, 1:h + 1] = sb[:, 1:]
    A3all[:, n - h:] = sb[:, 1:][:, ::-1]
    units = [u for u in range(1, n) if math.gcd(u, n) == 1]
    for a3, a1 in reps:
        A3s = A3all[A3all.sum(1, dtype=np.int64) == a3]
        P3 = pafs(A3s, h)
        # A1 with row sum a1: (n - a1) / 2 entries -1
        k = (n - a1) // 2
        t = time.time()
        table = {}
        for comb in itertools.combinations(range(1, n), k - 1) if k else [()]:
            # A1 up to shifts: put a -1 at position 0 (every A1 with k >= 1 has a shift with a1[0] = -1)
            x = np.ones(n, dtype=np.int8)
            if k:
                x[0] = -1
                x[list(comb)] = -1
            key = tuple(int((x.astype(np.int32) * np.roll(x, -s)).sum()) for s in range(1, h + 1))
            table.setdefault(key, []).append(x)
        log(f'  (a3, a1) = ({a3}, {a1}): |A3| {len(A3s)}, |A0| {len(A0s)}, A1 classes {len(table)} [{time.time() - t:.1f}s]')
        for i3 in range(len(A3s)):
            T = -(P3[i3][None, :] + P0)
            ok = (T % 2 == 0).all(1)
            for i0 in np.flatnonzero(ok):
                key = tuple(int(v) for v in T[i0] // 2)
                for x in table.get(key, ()):
                    A3, A1, A0 = A3s[i3], x, A0s[i0]
                    tot = verify.paf_sum([A3, A1, A1, A0])
                    assert tot[0] == 4 * n and not tot[1:].any()
                    canon = min(_canon(A3, A1, A0, u, n) for u in units)
                    if canon in seen:
                        continue
                    seen.add(canon)
                    sols.append((A3.copy(), A1.copy(), A0.copy()))
                    if not find_all:
                        return sols
    return sols


def subgroup(gens, n):
    H = {1}
    while True:
        G = H | {h * g % n for h in H for g in gens}
        if G == H:
            return sorted(H)
        H = G


def skew_propus_mult(n, gens, find_all=False, log=print):
    """Skew-propus quadruples of length n constant on the orbits of the multiplier group H = <gens> (-1 not in H):
    A1 and A0 on H-orbits (A0 of skew type), A3 on <H, -1>-orbits. Exhaustive over all orbit sign vectors."""
    H = subgroup(gens, n)
    if n - 1 in H:
        raise ValueError('-1 in H: no skew-type sequence is H-invariant')
    orbs, seen = [], set()
    for g in range(n):
        if g not in seen:
            o = sorted({g * h % n for h in H})
            seen |= set(o)
            orbs.append(o)
    idx = {g: i for i, o in enumerate(orbs) for g in o}
    neg = [idx[(-o[0]) % n] for o in orbs]
    O = len(orbs)
    h = (n - 1) // 2
    # A1: any sign per orbit
    def expand(bits):
        x = np.empty(n, dtype=np.int8)
        for i, o in enumerate(orbs):
            x[o] = bits[i]
        return x
    reps = [(a3, a1) for a3 in range(1, 2 * n, 2) for a1 in range(1, 2 * n, 2) if a3 * a3 + 2 * a1 * a1 == 4 * n - 1]
    pairs = [i for i in range(O) if i < neg[i]]                  # orbit pairs {o, -o} (o != 0)
    zero = idx[0]
    sols = []
    A0s = []
    for bits in itertools.product([1, -1], repeat=len(pairs)):
        b = [0] * O
        b[zero] = 1
        for i, s_ in zip(pairs, bits):
            b[i], b[neg[i]] = s_, -s_
        A0s.append(expand(b))
    A0s = np.array(A0s)
    A3s = []
    for bits in itertools.product([1, -1], repeat=len(pairs) + 1):
        b = [0] * O
        b[zero] = bits[0]
        for i, s_ in zip(pairs, bits[1:]):
            b[i] = b[neg[i]] = s_
        A3s.append(expand(b))
    A3s = np.array(A3s)
    A1all = np.array([expand(bits) for bits in itertools.product([1, -1], repeat=O)])
    P0 = pafs(A0s, h)
    log(f'n = {n}, H = <{",".join(map(str, gens))}> of order {len(H)}: {O} orbits; |A0| {len(A0s)}, |A3| {len(A3s)}, |A1| {len(A1all)}; rows {reps}')
    for a3, a1 in reps:
        X3 = A3s[A3s.sum(1, dtype=np.int64) == a3]
        X1 = A1all[np.abs(A1all.sum(1, dtype=np.int64)) == a1]
        P1 = pafs(X1, h)
        table = {}
        for i in range(len(X1)):
            table.setdefault(P1[i].tobytes(), []).append(i)
        P3 = pafs(X3, h)
        for i3 in range(len(X3)):
            T = -(P3[i3][None, :] + P0)
            ok = (T % 2 == 0).all(1)
            for i0 in np.flatnonzero(ok):
                key = (T[i0] // 2).astype(P1.dtype).tobytes()
                for i1 in table.get(key, ()):
                    A3, A1, A0 = X3[i3], X1[i1], A0s[i0]
                    tot = verify.paf_sum([A3, A1, A1, A0])
                    assert tot[0] == 4 * n and not tot[1:].any()
                    sols.append((A3.copy(), A1.copy(), A0.copy()))
                    if not find_all:
                        return sols
    return sols


def szekeres_mult(v, gens, find_all=False, log=print):
    """Legendre pairs on Z_v (v odd) with a of skew type and b symmetric (row sums 1), both constant on the orbits
    of the multiplier group H = <gens> (-1 not in H; b then on the <H, -1>-orbits): PAF_a + PAF_b = -2 off 0."""
    H = subgroup(gens, v)
    if v - 1 in H:
        raise ValueError('-1 in H')
    orbs, seen = [], set()
    for g in range(v):
        if g not in seen:
            o = sorted({g * h % v for h in H})
            seen |= set(o)
            orbs.append(o)
    idx = {g: i for i, o in enumerate(orbs) for g in o}
    neg = [idx[(-o[0]) % v] for o in orbs]
    pairs = [i for i in range(len(orbs)) if i < neg[i]]
    zero = idx[0]
    h = (v - 1) // 2
    def expand(b):
        x = np.empty(v, dtype=np.int8)
        for i, o in enumerate(orbs):
            x[o] = b[i]
        return x
    As, Bs = [], []
    for bits in itertools.product([1, -1], repeat=len(pairs)):
        b = [0] * len(orbs)
        b[zero] = 1
        for i, s_ in zip(pairs, bits):
            b[i], b[neg[i]] = s_, -s_
        As.append(expand(b))
        b0 = 1 if h % 2 == 0 else -1
        b2 = [0] * len(orbs)
        b2[zero] = b0
        for i, s_ in zip(pairs, bits):
            b2[i] = b2[neg[i]] = s_
        x = expand(b2)
        if x.sum() == 1:
            Bs.append(x)
    As, Bs = np.array(As), np.array(Bs)
    log(f'v = {v}, H = <{",".join(map(str, gens))}> of order {len(H)}: {len(orbs)} orbits, |a| {len(As)}, |b| {len(Bs)}')
    if not len(Bs):
        return []
    PB = pafs(Bs, h)
    table = {}
    for i in range(len(Bs)):
        table.setdefault(PB[i].tobytes(), []).append(i)
    PA = pafs(As, h)
    sols = []
    for i in range(len(As)):
        key = (-2 - PA[i]).astype(PB.dtype).tobytes()
        for j in table.get(key, ()):
            sols.append((As[i].copy(), Bs[j].copy()))
            if not find_all:
                return sols
    return sols


def _canon(A3, A1, A0, u, n):
    idx = (u * np.arange(n)) % n
    a3, a1, a0 = A3[idx], A1[idx], A0[idx]
    s1 = min(tuple(np.roll(a1, -s)) for s in range(n))
    a0t = a0[(-np.arange(n)) % n]
    s0 = min(tuple(a0), tuple(a0t))
    return (tuple(a3), s1, s0)


def fmt(x):
    return ''.join('+' if v > 0 else '-' for v in x)


def parse(s):
    return np.array([SIGN[c] for c in s], dtype=np.int64)


def test_szekeres():
    import subprocess
    for v in (5, 9, 13, 21, 29, 33, 41):
        out = subprocess.run([str(_P('/tmp/hadamard-work/bin/szek_pair')), str(v), '1'], capture_output=True, text=True).stdout
        line = [l for l in out.splitlines() if l.startswith('SOL')][0].split()
        a, b = parse(line[1]), parse(line[2])
        P, N = szekeres(a, b)
        S = np.eye(len(P), dtype=np.int64) + N @ P
        sk = np.array_equal(S + S.T, 2 * np.eye(len(P), dtype=np.int64)) and np.array_equal(S @ S.T, len(P) * np.eye(len(P), dtype=np.int64))
        ok = check_od(P, N)
        H = plug(P, N)
        print(f'Szekeres-type v = {v}: skew-Hadamard {sk}, symmetric OD({2 * v + 2}; 1, {2 * v + 1}) {ok}, order {H.shape[0]}: {is_symmetric_hadamard(H)}')


def test():
    for q in (3, 7, 11, 19, 23):
        P, N = paley(q)
        assert check_od(P, N), q
        H = plug(P, N)
        print(f'Paley q = {q}: symmetric OD({q + 1}; 1, {q}) ok, symmetric Hadamard of order {H.shape[0]}: {is_symmetric_hadamard(H)}')
    H2 = np.array([[1, 1], [1, -1]], dtype=np.int64)
    sylvester16 = np.kron(np.kron(H2, H2), np.kron(H2, H2))
    for q, M in ((7, paley2(5)), (11, sylvester16), (19, paley2(5, 1)), (23, paley2(13))):
        P, N = paley(q)
        H = plug3(P, N, M)
        print(f'Paley q = {q} with J - 2I and a symmetric Hadamard matrix of order {M.shape[0]}: order {H.shape[0]}: {is_symmetric_hadamard(H)}')
    for n in range(3, 16, 2):
        sols = skew_propus(n, log=lambda *a: None)
        if not sols:
            print(f'skew-propus n = {n}: none')
            continue
        A3, A1, A0 = sols[0]
        P, N = balonin(A3, A1, A0)
        assert check_od(P, N), n
        H = plug(P, N)
        print(f'skew-propus n = {n}: A3 {fmt(A3)} A1 {fmt(A1)} A0 {fmt(A0)}; order {H.shape[0]}: {is_symmetric_hadamard(H)}')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('cmd', choices=['test', 'propus', 'mpropus', 'mszek', 'build', 'certify'])
    ap.add_argument('args', nargs='*')
    ap.add_argument('--all', action='store_true')
    a = ap.parse_args()
    if a.cmd == 'test':
        test()
        test_szekeres()
    elif a.cmd == 'propus':
        n = int(a.args[0])
        t = time.time()
        sols = skew_propus(n, a.all)
        for A3, A1, A0 in sols:
            print(f'SOL n={n} A3 {fmt(A3)} A1 {fmt(A1)} A0 {fmt(A0)}')
        print(f'{len(sols)} solution(s) [{time.time() - t:.1f}s]')
    elif a.cmd == 'mpropus':
        n, gens = int(a.args[0]), [int(x) for x in a.args[1].split(',')]
        t = time.time()
        sols = skew_propus_mult(n, gens, a.all)
        for A3, A1, A0 in sols:
            print(f'SOL n={n} A3 {fmt(A3)} A1 {fmt(A1)} A0 {fmt(A0)}')
        print(f'{len(sols)} solution(s) [{time.time() - t:.1f}s]')
    elif a.cmd == 'mszek':
        v, gens = int(a.args[0]), [int(x) for x in a.args[1].split(',')]
        t = time.time()
        sols = szekeres_mult(v, gens, a.all)
        for x, y in sols:
            print(f'SOL v={v} a {fmt(x)} b {fmt(y)}')
        print(f'{len(sols)} solution(s) [{time.time() - t:.1f}s]')
    elif a.cmd == 'certify':
        from certify import write_certificate
        n = int(a.args[0])
        seqs = [parse(s) for s in a.args[1:4]]
        m = 4 * n
        order = m * (m - 1)
        out = ROOT / 'symmetric' / str(order)
        out.mkdir(parents=True, exist_ok=True)
        cert = {"order": order, "construction": "skew-propus-plug", "symmetric": True, "n": n, "m": m,
                "rowsums": [int(x.sum()) for x in seqs],
                "blocks": [np.flatnonzero(x < 0).tolist() for x in seqs]}
        H, _ = verify.build_skew_propus_plug(cert)
        sha = hashlib.sha256()
        with gzip.open(out / f'H{order}.csv.gz', 'wb', compresslevel=9) as f:
            for chunk in verify.csv_chunks(H):
                sha.update(chunk)
                f.write(chunk)
        cert["sha256"] = sha.hexdigest()
        write_certificate(cert, out / 'certificate.json')
        (out / 'sequences.txt').write_text(''.join(fmt(x) + '\n' for x in seqs))
        o, digest = verify.check(out)
        print(f'symmetric order {o}: sha256 {digest}')
    else:
        n = int(a.args[0])
        A3, A1, A0 = (parse(s) for s in a.args[1:4])
        P, N = balonin(A3, A1, A0)
        print('symmetric OD:', check_od(P, N))
        H = plug(P, N)
        print(f'order {H.shape[0]}: symmetric Hadamard {is_symmetric_hadamard(H)}')
        if len(a.args) > 4:
            np.save(a.args[4], H)


if __name__ == '__main__':
    main()
