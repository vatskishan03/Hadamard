"""Yamada's Hadamard matrix of order 4(q + 2) (after Kiyasu; search/NOTES.md section 13) in sixteen commuting bordered
blocks, times T-matrices, and the Seberry-Yamada Lemma 26 that makes the route classical (NOTES.md section 24).

Yamada: q = 5 (mod 8), symmetric circulants alpha, beta of order n0 = (q + 1)/2 with alpha^2 + beta^2 = qI, and the skew
core S of a skew-Hadamard matrix of order n0 + 1. The inner blocks A = P x beta + J2 x I, B = P x alpha,
C = D = J2 x S + P x I lie in P x circ(n0) + J2 x alg(S, J), which is commutative (P J2 = 0). With the border columns
reassigned to the block columns as (0, 3, 1, 2) and signed (1, -1, -1, -1), the sixteen bordered blocks commute and
diag(1, I2 x R) transposes them all whenever R S R = S^T for the reversal R; Paley cores (n0 a prime) have this. Then
A_j = sum_k T_k x N_kj go into the Goethals-Seidel array (Seberry-Yamada 1990, Theorem 4).

    python3 search/yamada_blocks.py check q t      prime q = 5 mod 8 with (q + 1)/2 prime; T-matrices of order t from
                                                   base sequences; checks the blocks and H H^T = 4(q + 2)t I
    python3 search/yamada_blocks.py lemma26 w      prime w = s^2 + 4: bordered two-circulant skew matrix of order
                                                   2(w + 1) with an involutive anti-automorphism fixing two points, and
                                                   Lemma 26 (Corollary 27(i) form): symmetric Williamson-type matrices
                                                   of order 4w + 3 when 4w + 1 is a prime
"""
import itertools, sys
from pathlib import Path
import numpy as np
from sympy import isprime, primefactors, primitive_root

sys.path.insert(0, str(Path(__file__).resolve().parent))
import paley_blocks as PB


def circ(x):
    n = len(x); t = np.arange(n)
    return x[(t[None, :] - t[:, None]) % n]                  # X[i, j] = x(j - i)


def legendre(p):
    chi = -np.ones(p, np.int64); chi[np.arange(1, p) ** 2 % p] = 1; chi[0] = 0
    return chi


def reversal(n):
    R = np.zeros((n, n), np.int64); R[np.arange(n), (-np.arange(n)) % n] = 1
    return R


def alpha_beta(q):
    """Yamada's Lemma 3: first rows psi(Tr xi^(4m)), psi(Tr xi^(4m + n0)) over GF(q^2) = GF(q)[w]/(w^2 - r); beta(0) = 0"""
    psi = legendre(q)
    r = next(x for x in range(2, q) if psi[x] == -1)
    mul = lambda x, y: ((x[0] * y[0] + r * x[1] * y[1]) % q, (x[0] * y[1] + x[1] * y[0]) % q)

    def power(x, e):
        y = (1, 0)
        while e:
            if e & 1: y = mul(y, x)
            x, e = mul(x, x), e >> 1
        return y
    N = q * q - 1
    xi = next((a, b) for a in range(q) for b in range(1, q)
              if all(power((a, b), N // f) != (1, 0) for f in primefactors(N)))
    n0 = (q + 1) // 2
    al, be = np.empty(n0, np.int64), np.empty(n0, np.int64)
    a, b, x4 = (1, 0), power(xi, n0), power(xi, 4)
    for m in range(n0):
        al[m], be[m] = psi[2 * a[0] % q], psi[2 * b[0] % q]
        a, b = mul(a, x4), mul(b, x4)
    assert np.array_equal(circ(al) @ circ(al) + circ(be) @ circ(be), q * np.eye(n0, dtype=np.int64))
    return al, be


def blocks(q):
    """the sixteen bordered blocks (order q + 2) of Yamada's matrix with the Paley core of order (q + 1)/2"""
    n0 = (q + 1) // 2
    al, be = alpha_beta(q)
    S = circ(legendre(n0))
    P, J2, I = np.array([[1, -1], [-1, 1]]), np.ones((2, 2), np.int64), np.eye(n0, dtype=np.int64)
    A, B = np.kron(P, circ(be)) + np.kron(J2, I), np.kron(P, circ(al))
    C = np.kron(J2, S) + np.kron(P, I); D = C
    L = 2 * np.eye(4, dtype=np.int64) - np.ones((4, 4), np.int64)
    M = np.array([[1, 1, 1, 1], [1, -1, -1, 1], [1, 1, -1, -1], [1, -1, 1, -1]])
    K4 = -(L @ M) // 2
    H0 = [[A, B, C, D], [-B.T, A.T, -D.T, C.T], [-C.T, D, A.T, -B], [-D.T, -C, B.T, A]]
    rho, dl = (0, 3, 1, 2), (1, -1, -1, -1)
    m = 2 * n0
    N = [[None] * 4 for _ in range(4)]
    for i in range(4):
        for j in range(4):
            X = np.zeros((m + 1, m + 1), np.int64)
            X[0, 0], X[0, 1:], X[1:, 0] = dl[j] * L[i, rho[j]], K4[i, j], dl[j] * M[rho[j], i]
            X[1:, 1:] = H0[i][j]
            N[i][j] = X
    RB = np.zeros((m + 1, m + 1), np.int64); RB[0, 0] = 1; RB[1:, 1:] = np.kron(np.eye(2, dtype=np.int64), reversal(n0))
    return N, RB


def check(q, t):
    assert q % 8 == 5 and isprime(q) and isprime((q + 1) // 2)
    N, RB = blocks(q)
    m = q + 2
    H = np.block(N)
    assert np.array_equal(H @ H.T, 4 * m * np.eye(4 * m, dtype=np.int64)), "Yamada's matrix"
    Ns = [x for row in N for x in row]
    assert all(np.array_equal(x @ y, y @ x) for x in Ns for y in Ns), "blocks do not commute"
    assert all(np.array_equal(RB @ x @ RB, x.T) for x in Ns), "reversal does not transpose"
    T = PB.t_matrices(t)
    Ct = [circ(r) for r in T]
    Aj = [sum(np.kron(Ct[k], N[k][j]) for k in range(4)) for j in range(4)]
    R = np.kron(reversal(t), RB)
    AR = [x @ R for x in Aj]; ATR = [x.T @ R for x in Aj]
    G = np.block([[Aj[0], AR[1], AR[2], AR[3]], [-AR[1], Aj[0], ATR[3], -ATR[2]],
                  [-AR[2], -ATR[3], Aj[0], ATR[1]], [-AR[3], ATR[2], -ATR[1], Aj[0]]])
    n = 4 * m * t
    assert set(np.unique(G)) <= {-1, 1} and np.array_equal(G @ G.T, n * np.eye(n, dtype=np.int64))
    print(f"q = {q}, t = {t}: sixteen commuting bordered blocks of order {m}; H H^T = {n} I")


def skew_two_circulant(w):
    """skew S of order 2(w + 1) from U = C0 u C1 (skew) and V = squares (symmetric), PAF_U + PAF_V = -2 (w = s^2 + 4),
    and the signed involution sigma with sigma S sigma = -S: borders swapped, x -> -x on both orbits (signs +, -)"""
    g = primitive_root(w)
    ind = {pow(g, k, w): k for k in range(w - 1)}
    U, V = np.zeros(w, np.int64), np.zeros(w, np.int64)
    for x in range(1, w):
        U[x] = 1 if ind[x] % 4 in (0, 1) else -1
        V[x] = 1 if ind[x] % 4 in (0, 2) else -1
    paf = lambda v: np.array([int(v @ np.roll(v, -d)) for d in range(w)])
    for V[0] in (1, -1):
        if ((paf(U) + paf(V))[1:] == -2).all():
            break
    else:
        raise ValueError("w is not of the form s^2 + 4")
    Uc, Vc, N = circ(U), circ(V), 2 * w + 2
    for u1, u2, s, c in itertools.product((1, -1), repeat=4):
        S = np.zeros((N, N), np.int64)
        S[0, 1], S[1, 0] = s, -s
        S[0, 2:2 + w], S[0, 2 + w:] = u1, u2
        S[1, 2:2 + w], S[1, 2 + w:] = -c * u1, c * u2
        S[2:, :2] = -S[:2, 2:].T
        S[2:2 + w, 2:2 + w], S[2:2 + w, 2 + w:] = Uc, Vc
        S[2 + w:, 2:2 + w], S[2 + w:, 2 + w:] = -Vc.T, -Uc
        sig = np.zeros((N, N), np.int64); sig[0, 1] = sig[1, 0] = c
        for x in range(w):
            sig[2 + x, 2 + (-x) % w] = 1; sig[2 + w + x, 2 + w + (-x) % w] = -1
        if np.array_equal(S.T, -S) and np.array_equal(S @ S.T, (N - 1) * np.eye(N, dtype=np.int64)) \
                and np.array_equal(sig @ S @ sig, -S):
            return S, sig
    raise ValueError("no border")


def lemma26(w):
    S, sig = skew_two_circulant(w)
    N = S.shape[0]
    r0 = next(i for i in range(N) if sig[i, i] == 1)              # the fixed point with sign +
    d = S[r0].copy(); d[r0] = 1
    T = S * d[:, None] * d[None, :]
    rest = [i for i in range(N) if i != r0]
    B = T[np.ix_(rest, rest)]                                     # skew core of order n0 = 2w + 1
    tau = (np.diag(d) @ sig @ np.diag(d))[np.ix_(rest, rest)]
    tau = -tau if (tau.sum(1) == -1).all() else tau
    assert (tau >= 0).all() and np.array_equal(tau @ B @ tau, -B)
    n0 = B.shape[0]; q = 2 * n0 - 1
    assert q % 8 == 5 and isprime(q), "4w + 1 must be a prime = 5 (mod 8)"
    perm, rev = tau.argmax(1), (-np.arange(n0)) % n0               # relabel so that tau becomes the reversal
    def pairs(pm):
        fx, pr, seen = [], [], set()
        for i in range(n0):
            if i not in seen:
                seen |= {i, pm[i]}; (fx if pm[i] == i else pr).append((i, pm[i]))
        return fx, pr
    (fT, pT), (fR, pR) = pairs(perm), pairs(rev)
    phi = np.empty(n0, np.int64); phi[fR[0][0]] = fT[0][0]
    for (a, b), (x, y) in zip(pR, pT): phi[a], phi[b] = x, y
    B = B[np.ix_(phi, phi)]
    R = reversal(n0)
    assert np.array_equal(R @ B @ R, -B)
    al, be = alpha_beta(q)
    I, Z = np.eye(n0, dtype=np.int64), np.zeros((n0, n0), np.int64)
    U, V = [I, Z, B @ R, B @ R], [circ(be), circ(al), R, R]
    J2, P = np.ones((2, 2), np.int64), np.array([[1, -1], [-1, 1]])
    m = 2 * n0 + 1
    for signs in itertools.product((1, -1), repeat=4):
        X = []
        for i in range(4):
            Xi = np.zeros((m, m), np.int64); Xi[0, 0] = 1
            Xi[0, 1:] = Xi[1:, 0] = signs[i]
            Xi[1:, 1:] = np.kron(U[i], J2) + np.kron(V[i], P)
            X.append(Xi)
        if np.array_equal(sum(x @ x.T for x in X), 4 * m * np.eye(m, dtype=np.int64)) and \
                all(np.array_equal(X[i] @ X[j].T, X[j] @ X[i].T) for i in range(4) for j in range(4)):
            assert all(np.array_equal(x, x.T) for x in X)
            print(f"w = {w}: skew core of order {n0} with an involutive anti-automorphism; symmetric Williamson-type "
                  f"matrices of order {m} (border signs {signs})")
            return X
    raise ValueError("no border signs")


if __name__ == "__main__":
    if sys.argv[1] == "check":
        check(int(sys.argv[2]), int(sys.argv[3]))
    elif sys.argv[1] == "lemma26":
        lemma26(int(sys.argv[2]))
