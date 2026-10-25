#!/usr/bin/env python3
"""Symmetric Williamson-type matrices: four symmetric +-1 matrices X_1..X_4 of order w that commute pairwise (for
symmetric matrices that is the same as pairwise amicable) with X_1^2 + X_2^2 + X_3^2 + X_4^2 = 4w I. Their Williamson
array is the input K of the Miyamoto blocks of search/miyamoto_blocks.py (hadamard/MIYAMOTO_BLOCKS.md), for
q = 4w + 1 a prime.

Builders, all from Seberry and Yamada (J. Combin. Math. Combin. Comput. 7 (1990); 1992 survey, Section 8):

    circ(w)              symmetric circulant Williamson matrices of order w (meet in the middle; w <= 13)
    circm(w, g1, ...)    the same, constant on the orbits of <g1, ..., -1> (search/miyamoto_blocks.py williamson_mult)
    nw(w)                (A R, B, C, D) from the near-Williamson matrices near-williamson/NW<w>.txt
    turyn(q)             Turyn's Williamson matrices of order (q + 1)/2 for a prime power q = 1 (mod 4) (J. Combin. Theory
                         A 12 (1972)): (I + B, I - B, A, A) with Yamada's symmetric circulants A, B (B of zero diagonal),
                         A^2 + B^2 = q I (yamada_pair)
    one()                order 1: (1, 1, 1, 1)
    cor25(U, V)          order 2uv from orders u and v (Corollary 25 / 8.5): E, F = (A +- B)/2, G, H = (C +- D)/2,
                         T = (E (x) I2, F (x) I2, G (x) s, H (x) s) with s = [[0, 1], [1, 0]], and
                         L = T1 X + T2 Y + T3 Z + T4 W, M = -T1 Y + T2 X + T3 W - T4 Z,
                         N = -T1 Z - T2 W + T3 X + T4 Y, P = -T1 W + T2 Z - T3 Y + T4 X (Kronecker products)
    cor27i(p)            order p + 2 for a prime p = 5 (mod 8) with n = (p + 1)/2 a prime (Corollary 27(i)): Lemma 26 with
                         U = (I, 0, B R, B R), V = (M, N, R, R), B the Paley core of order n, R the reversal, M, N
                         Yamada's symmetric circulants of order n with M^2 + N^2 = p I (M of zero diagonal)
    cor31(q)             order 2q + 1 for a prime q = 1 (mod 4) with 2q - 1 = r a prime power (Corollary 31): Lemma 26
                         with U = (I, 0, Q, Q), V = (M, N, I, I), Q Paley's symmetric core of order q and M, N symmetric
                         circulants of order q with M^2 + N^2 = r I, M of zero diagonal: Yamada's pair (yamada_pair)
    cor27c(p, NAME)      order p + 2 for a prime p = 1 (mod 4) from a symmetric conference matrix of order (p + 3)/2 in
                         search/data/NAME (rows of 0, +, -, first row and column +; Corollary 27): Lemma 26 with
                         U = (I, 0, Q, Q), V = (M, N, I, I), Q its core and M, N Yamada's symmetric circulants of order
                         (p + 1)/2 with M^2 + N^2 = p I (yamada_pair), e.g.
                         cor27c(233, "conference118_5fe566972f21790b.txt")
    cor29c(q)            order q for a prime power q = 1 (mod 4) with (q - 1)/2 - 1 = p a prime power = 1 (mod 4) (Corollary
                         29(ii) = 8.8(ii) with Paley's symmetric conference matrix C of order (q - 1)/2: M3, M4 = C +- I)
    cor29r(q, W)         the same with M3 = [[Y, X], [X, -Z]], M4 = [[Z, X], [X, -Y]] from symmetric Williamson-type
                         matrices W = (Y, X, X, Z) of order (q - 1)/4 with a repeated block (NOTES.md section 7), e.g.
                         cor29r(773, cor29h(193, ...)): cor29h, cor27i and cor31 give a repeated block
    cor29h(q, H)         the same with a symmetric Hadamard matrix H of order (q - 1)/2 (M3 = M4 = H), e.g.
                         kron(paley2(p), h2(), ...): Paley II, [[C + I, C - I], [C - I, -C - I]] of order 2(p + 1) for a
                         prime power p = 1 (mod 4), times Sylvester's h2() = [[1, 1], [1, -1]]
    cor29q(q, f1, f2, ...) the same with (M3, M4) any symmetric, commuting +-1 pair of order (q - 1)/2 with
                         M3^2 + M4^2 = (q - 1) I, that is a symmetric quaternary Hadamard matrix W = (M3 + M4)/2 +
                         i (M4 - M3)/2 (NOTES.md section 19), built as the Kronecker product of factors: qconf(p) =
                         (C + I, C - I) with Paley's conference matrix of order p + 1 (Turyn), qh(H) = (H, H) with a
                         symmetric Hadamard matrix, qwt(W) = the repeated-block arrangement of cor29r. It contains
                         cor29c, cor29h and cor29r as the cases with one factor; e.g. cor29q(2953, qconf(17), qconf(81)).

    python3 search/symwt.py "cor25(circ(7), circ(17))"      # builds and checks; prints the order
"""
import itertools, math, sys
from pathlib import Path
import numpy as np
from sympy import isprime

ROOT = Path(__file__).resolve().parent.parent
sys.path[:0] = [str(ROOT / "search"), str(ROOT / "verify")]
import miyamoto_blocks as MB
import miyamoto_wt as MW
import verify as V


def check(X):
    """symmetric +-1, pairwise commuting, sum of squares 4w I (products in float64, exact for these sizes)"""
    w = X[0].shape[0]
    assert len(X) == 4 and all(x.shape == (w, w) for x in X)
    assert all(set(np.unique(x)) <= {-1, 1} for x in X), "entries"
    assert all(np.array_equal(x, x.T) for x in X), "symmetric"
    F = [x.astype(np.float64) for x in X]
    assert all(np.array_equal(F[i] @ F[j], F[j] @ F[i]) for i in range(4) for j in range(i + 1, 4)), "commute"
    assert np.array_equal(sum(f @ f for f in F), 4.0 * w * np.eye(w)), "sum of squares"
    return X


def circ(w):
    return check([V.circulant(r) for r in MB.williamson(w)])


def circm(w, *gens):
    W = MB.williamson_mult(w, list(gens))
    assert W is not None, "no Williamson matrices with this multiplier group"
    return check([V.circulant(r) for r in W])


def nw(w):
    a, b, c, d = MB.load_nw(f"near-williamson/NW{w}.txt")
    return check([V.circulant(a) @ V.back_identity(w), V.circulant(b), V.circulant(c), V.circulant(d)])


def turyn(q):
    al, be = yamada_pair(q)
    A, B = V.circulant(al), V.circulant(be)
    I = np.eye(len(al), dtype=np.int64)
    return check([I + B, I - B, A, A])


def one():
    return [np.ones((1, 1), np.int64)] * 4


def cor25(U, Vv):
    A, B, C, D = U
    E, F, G, H = (A + B) // 2, (A - B) // 2, (C + D) // 2, (C - D) // 2
    I2, s = np.eye(2, dtype=np.int64), np.array([[0, 1], [1, 0]], np.int64)
    T1, T2, T3, T4 = np.kron(E, I2), np.kron(F, I2), np.kron(G, s), np.kron(H, s)
    X, Y, Z, W = Vv
    k = np.kron
    return check([k(T1, X) + k(T2, Y) + k(T3, Z) + k(T4, W), -k(T1, Y) + k(T2, X) + k(T3, W) - k(T4, Z),
                  -k(T1, Z) - k(T2, W) + k(T3, X) + k(T4, Y), -k(T1, W) + k(T2, Z) - k(T3, Y) + k(T4, X)])


def lemma26(U, Vv):
    """Seberry-Yamada Lemma 26: symmetric Williamson-type matrices of order 2n + 1 from U_i, V_i of order n (border
    signs found by trying all 16)"""
    n = U[0].shape[0]
    J2, P = np.ones((2, 2), np.int64), np.array([[1, -1], [-1, 1]], np.int64)
    S = [np.kron(U[i], J2) + np.kron(Vv[i], P) for i in range(4)]
    m = 2 * n + 1
    for signs in itertools.product((1, -1), repeat=4):
        X = []
        for i in range(4):
            Xi = np.zeros((m, m), np.int64)
            Xi[0, 0] = 1
            Xi[0, 1:] = Xi[1:, 0] = signs[i]
            Xi[1:, 1:] = S[i]
            X.append(Xi)
        try:
            return check(X)
        except AssertionError:
            continue
    raise ValueError("Lemma 26: no border signs")


def legendre(p):
    chi = -np.ones(p, np.int64); chi[np.arange(1, p) ** 2 % p] = 1; chi[0] = 0
    return chi


def reversal(n):
    R = np.zeros((n, n), np.int64); R[np.arange(n), (-np.arange(n)) % n] = 1
    return R


def cor27i(p):
    import yamada_blocks as Y
    n = (p + 1) // 2
    assert p % 8 == 5 and isprime(p) and isprime(n)
    al, be = Y.alpha_beta(p)
    B = V.circulant(legendre(n))                 # Paley core, skew since n = 3 (mod 4)
    R = reversal(n)
    I, Z = np.eye(n, dtype=np.int64), np.zeros((n, n), np.int64)
    return lemma26([I, Z, B @ R, B @ R], [V.circulant(be), V.circulant(al), R, R])


def mn_pair(q, r):
    """symmetric circulants M (zero diagonal) and N of order q with M^2 + N^2 = r I, by a meet in the middle on the
    periodic autocorrelations, with the row sums 0 and sqrt(r) (r an odd square)"""
    h, s = (q - 1) // 2, math.isqrt(r)
    assert s * s == r and s % 2 == 1
    def rows(d0, plus):
        out = []
        for P in itertools.combinations(range(h), plus):
            x = -np.ones(q, np.int64); x[0] = d0
            for k in P:
                x[k + 1] = x[q - k - 1] = 1
            out.append(x)
        return np.array(out)
    paf = lambda X: np.rint(np.fft.irfft(np.abs(np.fft.rfft(X, axis=1)) ** 2, n=q, axis=1)).astype(np.int64)[:, 1:h + 1]
    Ms = rows(0, h // 2) if h % 2 == 0 else None
    Ns = rows(1, (s - 1 + 2 * h) // 4) if (s - 1 + 2 * h) % 4 == 0 else rows(-1, (s + 1 + 2 * h) // 4)
    left = {}
    for x, pa in zip(Ms, paf(Ms)):
        left.setdefault(pa.tobytes(), x)
    for y, pb in zip(Ns, paf(Ns)):
        x = left.get((-pb).tobytes())
        if x is not None:
            M, N = V.circulant(x), V.circulant(y)
            assert np.array_equal(M @ M + N @ N, r * np.eye(q, dtype=np.int64))
            return M, N
    raise ValueError("no pair")


def cor31(q):
    r = 2 * q - 1
    Q = V.circulant(legendre(q))                 # symmetric core, q = 1 (mod 4)
    al, be = yamada_pair(r)
    M, N = V.circulant(be), V.circulant(al)
    I, Z = np.eye(q, dtype=np.int64), np.zeros((q, q), np.int64)
    return lemma26([I, Z, Q, Q], [M, N, I, I])


def cor27c(p, name):
    rows = (ROOT / "search" / "data" / name).read_text().split()
    C = np.array([[{"0": 0, "+": 1, "-": -1}[c] for c in r] for r in rows], np.int64)
    n = C.shape[0] - 1
    assert 2 * n == p + 1 and (C[0, 1:] == 1).all() and (C[1:, 0] == 1).all()
    Q = C[1:, 1:]                                # symmetric core, Q^2 = nI - J
    assert np.array_equal(Q, Q.T) and np.array_equal(Q @ Q, n * np.eye(n, dtype=np.int64) - 1)
    al, be = yamada_pair(p)
    M, N = V.circulant(be), V.circulant(al)
    I, Z = np.eye(n, dtype=np.int64), np.zeros((n, n), np.int64)
    return lemma26([I, Z, Q, Q], [M, N, I, I])


def gf_chi(q):
    """GF(q), q = p^k, as integers 0..q-1 (base-p digits = coefficients): subtraction table and quadratic character"""
    D, chi, _, _ = gf_tables(q)
    return D, chi


def gf_tables(q):
    """GF(q) as integers 0..q-1: subtraction table D[x, y] = x - y, quadratic character, addition and multiplication
    tables"""
    from sympy import factorint
    (p, k), = factorint(q).items()
    digits = lambda x: [(x // p ** i) % p for i in range(k)]
    num = lambda d: sum(c * p ** i for i, c in enumerate(d))
    def mulmod(a, b, f):
        r = [0] * (2 * k - 1)
        for i, x in enumerate(a):
            for j, y in enumerate(b):
                r[i + j] = (r[i + j] + x * y) % p
        for d in range(2 * k - 2, k - 1, -1):            # reduce by the monic f of degree k
            c = r[d]
            if c:
                for i in range(k + 1):
                    r[d - k + i] = (r[d - k + i] - c * f[i]) % p
        return r[:k]
    f = [0, 1] if k == 1 else next(c + [1] for c in (digits(x) for x in range(q)) if _irreducible(c + [1], p))
    dg = [digits(x) for x in range(q)]
    M = np.array([[num(mulmod(dg[x], dg[y], f)) for y in range(q)] for x in range(q)])
    sq = set(M[np.arange(1, q), np.arange(1, q)].tolist())
    chi = np.array([0] + [1 if x in sq else -1 for x in range(1, q)], np.int64)
    D = np.array([[num([(a - b) % p for a, b in zip(dg[x], dg[y])]) for y in range(q)] for x in range(q)])
    A = np.array([[num([(a + b) % p for a, b in zip(dg[x], dg[y])]) for y in range(q)] for x in range(q)])
    return D, chi, A, M


def _irreducible(f, p):
    """f (coefficients, lowest first, monic) has no monic factor over GF(p) of degree <= deg f / 2"""
    k = len(f) - 1
    def divides(g):
        r = list(f)
        for d in range(k, len(g) - 2, -1):
            c = r[d]
            if c:
                for i in range(len(g)):
                    r[d - len(g) + 1 + i] = (r[d - len(g) + 1 + i] - c * g[i]) % p
        return not any(r[:len(g) - 1])
    return not any(divides(list(g) + [1]) for e in range(1, k // 2 + 1) for g in itertools.product(range(p), repeat=e))


def yamada_pair(q):
    """Yamada's symmetric circulants alpha, beta of order (q + 1)/2 with alpha^2 + beta^2 = q I, beta of zero diagonal,
    for a prime power q = 1 (mod 4): psi(Tr xi^(4m)), psi(Tr xi^(4m + (q + 1)/2)) over GF(q^2) = GF(q)[w]/(w^2 - r)
    (search/yamada_blocks.py alpha_beta for primes)"""
    from sympy import primefactors
    _, psi, A, M = gf_tables(q)
    r = next(x for x in range(2, q) if psi[x] == -1)
    mul = lambda x, y: (A[M[x[0], y[0]], M[r, M[x[1], y[1]]]], A[M[x[0], y[1]], M[x[1], y[0]]])
    def power(x, e):
        y = (1, 0)
        while e:
            if e & 1: y = mul(y, x)
            x, e = mul(x, x), e >> 1
        return y
    N = q * q - 1
    xi = next((a, b) for a in range(q) for b in range(1, q) if all(power((a, b), N // f) != (1, 0) for f in primefactors(N)))
    n0 = (q + 1) // 2
    al, be = np.empty(n0, np.int64), np.empty(n0, np.int64)
    a, b, x4 = (1, 0), power(xi, n0), power(xi, 4)
    for m in range(n0):
        al[m], be[m] = psi[A[a[0], a[0]]], psi[A[b[0], b[0]]]
        a, b = mul(a, x4), mul(b, x4)
    Ca, Cb = V.circulant(al), V.circulant(be)
    assert be[0] == 0 and np.array_equal(Ca, Ca.T) and np.array_equal(Cb, Cb.T)
    assert np.array_equal(Ca @ Ca + Cb @ Cb, q * np.eye(n0, dtype=np.int64))
    return al, be


def conference(q):
    """Paley's symmetric conference matrix of order q + 1 (q = 1 mod 4 a prime power)"""
    assert q % 4 == 1
    D, chi = gf_chi(q)
    C = np.ones((q + 1, q + 1), np.int64); C[0, 0] = 0
    C[1:, 1:] = chi[D.T]                                # C[x, y] = chi(y - x)
    assert np.array_equal(C, C.T) and np.array_equal(C @ C, q * np.eye(q + 1, dtype=np.int64))
    return C


def paley2(p):
    C = conference(p)
    I = np.eye(p + 1, dtype=np.int64)
    H = np.block([[C + I, C - I], [C - I, -C - I]])
    assert np.array_equal(H, H.T) and np.array_equal(H @ H, len(H) * np.eye(len(H), dtype=np.int64))
    return H


def h2():
    return np.array([[1, 1], [1, -1]], np.int64)


def kron(*Hs):
    H = np.ones((1, 1), np.int64)
    for x in Hs:
        H = np.kron(H, x)
    return H


def miyamoto_lemma(q, M3, M4):
    """Seberry-Yamada's Lemma 26 in Miyamoto's form (search/miyamoto_wt.py williamson_type) for a prime power q = 1
    (mod 4): U_1 = -[chi(g^(2i) - g^(2j))], U_2 = [chi(g^(2i) - g^(1 - 2j))], V = (I, 0, M3, M4), g primitive in GF(q)"""
    if isprime(q):
        return MW.williamson_type(q, M3, M4)
    from sympy import primefactors
    D, chi, _, M = gf_tables(q)
    def power(x, e):
        y = 1
        while e:
            if e & 1: y = M[y, x]
            x, e = M[x, x], e >> 1
        return y
    g = next(x for x in range(2, q) if all(power(x, (q - 1) // f) != 1 for f in primefactors(q - 1)))
    n = (q - 1) // 2
    S = [power(g, 2 * j) for j in range(n)]
    Nr = [power(g, (1 - 2 * j) % (q - 1)) for j in range(n)]
    QSS, QSN = chi[D[np.ix_(S, S)]], chi[D[np.ix_(S, Nr)]]
    I, Z = np.eye(n, dtype=np.int64), np.zeros((n, n), np.int64)
    U, Vv = [-QSS, QSN, Z, Z], [I, Z, M3, M4]
    J2, D2 = np.ones((2, 2), np.int64), np.array([[1, -1], [-1, 1]], np.int64)
    e = np.ones((1, 2 * n), np.int64)
    X = []
    for i in range(4):
        sg = -1 if i == 0 else 1
        Si = np.kron(U[i], J2) + np.kron(Vv[i], D2)
        X.append(np.block([[np.ones((1, 1), np.int64), sg * e], [sg * e.T, Si]]))
    return X


def cor29c(q):
    n = (q - 1) // 2
    C = conference(n - 1)
    I = np.eye(n, dtype=np.int64)
    return check(miyamoto_lemma(q, C + I, C - I))


def cor29h(q, H):
    assert np.array_equal(H, H.T) and np.array_equal(H @ H, len(H) * np.eye(len(H), dtype=np.int64))
    return check(miyamoto_lemma(q, H, H))


def cor29r(q, W):
    w = W[0].shape[0]
    assert 4 * w + 1 == q
    i, j = next((i, j) for i in range(4) for j in range(i + 1, 4) if np.array_equal(W[i], W[j]))
    X = W[i]
    Y, Z = (W[k] for k in range(4) if k not in (i, j))
    return check(miyamoto_lemma(q, np.block([[Y, X], [X, -Z]]), np.block([[Z, X], [X, -Y]])))


def qconf(p):
    C = conference(p)
    I = np.eye(p + 1, dtype=np.int64)
    return C + I, C - I


def qh(H):
    assert np.array_equal(H, H.T) and np.array_equal(H @ H, len(H) * np.eye(len(H), dtype=np.int64))
    return H, H


def qwt(W):
    i, j = next((i, j) for i in range(4) for j in range(i + 1, 4) if np.array_equal(W[i], W[j]))
    X = W[i]
    Y, Z = (W[k] for k in range(4) if k not in (i, j))
    return np.block([[Y, X], [X, -Z]]), np.block([[Z, X], [X, -Y]])


def cor29q(q, *factors):
    """Lemma 26 with M3 + i M4 = (1 + i) (W_1 (x) W_2 (x) ...), W_k = ((A + B) + i (B - A))/2 for the factor (A, B)"""
    re, im = np.ones((1, 1), np.int64), np.zeros((1, 1), np.int64)
    for A, B in factors:
        a, b = (A + B) // 2, (B - A) // 2
        re, im = np.kron(re, a) - np.kron(im, b), np.kron(re, b) + np.kron(im, a)
    M3, M4 = re - im, re + im
    n = (q - 1) // 2
    assert M3.shape == (n, n) and np.array_equal(M3, M3.T) and np.array_equal(M4, M4.T)
    F3, F4 = M3.astype(np.float64), M4.astype(np.float64)
    assert np.array_equal(F3 @ F4, F4 @ F3) and np.array_equal(F3 @ F3 + F4 @ F4, 2.0 * n * np.eye(n))
    return check(miyamoto_lemma(q, M3, M4))


def build(recipe):
    env = {f: globals()[f] for f in ("circ", "circm", "nw", "turyn", "one", "cor25", "cor27i", "cor27c", "cor31",
                                     "cor29c", "cor29h", "cor29r", "cor29q", "qconf", "qh", "qwt", "paley2", "h2", "kron")}
    return check(eval(recipe, {"__builtins__": {}}, env))


if __name__ == "__main__":
    X = build(sys.argv[1])
    print(f"{sys.argv[1]}: symmetric Williamson-type matrices of order {X[0].shape[0]}, checked")
