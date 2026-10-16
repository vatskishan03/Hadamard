#!/usr/bin/env python3
"""Paley's matrix of order q + 1 (q = 3 mod 4 prime) is of quaternion type but not of Williamson type (NOTES.md §21).

Index the points of PG(1, q) by vectors of F_q^2 = F_q[i] and let H[u, v] = chi(lambda_u lambda_v) (delta + chi(det)) on
representatives; SL(2, q) acts monomially with -1 acting as -I.  Q8 = <x -> i x, x -> c conj(x)> (N(c) = -1) lies in
the dicyclic group of order 2(q + 1) and acts semiregularly, so in a basis adapted to its N = (q + 1)/4 orbits
H is the Williamson array of four blocks X_1, X_i, X_j, X_k.  H H^T = 4N I then needs only [X_i, X_1] = -[X_k, X_j] and
its cyclic versions ([X, Y] = X Y^T - Y X^T), not pairwise amicability.

    python3 search/quaternion_type.py pattern 19     # block pattern and amicability of X_1, X_i, X_j, X_k
    python3 search/quaternion_type.py rebase 19      # every choice of base point per orbit (rows and columns)
    python3 search/quaternion_type.py plug 19 3      # substitute into the Cooper-Wallis array of T-matrices of order t
"""
import itertools
import sys
import numpy as np


class Paley:
    def __init__(self, q):
        assert q % 4 == 3
        self.q, self.N = q, (q + 1) // 4
        self.chi = np.array([0] + [1 if pow(x, (q - 1) // 2, q) == 1 else -1 for x in range(1, q)], np.int64)
        c = next((a, b) for a in range(q) for b in range(q) if (a * a + b * b + 1) % q == 0)
        al = lambda u: (-u[1] % q, u[0] % q)                                        # x -> i x
        be = lambda u: ((c[0] * u[0] + c[1] * u[1]) % q, (c[1] * u[0] - c[0] * u[1]) % q)   # x -> c conj(x)
        self.G = [lambda u: u, al, be, lambda u: al(be(u))]                         # 1, i, j, k (mod -1)
        seen, self.reps = set(), []
        for p in [(1, y) for y in range(q)] + [(0, 1)]:
            if p in seen:
                continue
            self.reps.append(p)
            seen |= {self.norm(g(p))[0] for g in self.G}
        assert len(self.reps) == self.N

    def norm(self, u):
        a, b = u[0] % self.q, u[1] % self.q
        return ((1, b * pow(a, -1, self.q) % self.q), a) if a else ((0, 1), b)

    def entry(self, u, v):
        (ru, lu), (rv, lv) = self.norm(u), self.norm(v)
        d = (ru[0] * rv[1] - ru[1] * rv[0]) % self.q
        return self.chi[lu * lv % self.q] * ((1 if ru == rv else 0) + self.chi[d])

    def matrix(self):
        vec = [g(r) for g in self.G for r in self.reps]
        return np.array([[self.entry(u, v) for v in vec] for u in vec], np.int64)


def blocks(H, N):
    return [H[:N, h * N:(h + 1) * N] for h in range(4)]


def amicable(X):
    return all(np.array_equal(X[a] @ X[b].T, X[b] @ X[a].T) for a in range(4) for b in range(a + 1, 4))


def pattern(q):
    P = Paley(q)
    H, N = P.matrix(), P.N
    assert np.array_equal(H @ H.T, 4 * N * np.eye(4 * N, dtype=np.int64))
    X = blocks(H, N)
    pat = [[next(f"{'+-'[s < 0]}X{'1ijk'[w]}" for w in range(4) for s in (1, -1)
                 if np.array_equal(H[g * N:(g + 1) * N, h * N:(h + 1) * N], s * X[w])) for h in range(4)] for g in range(4)]
    br = lambda a, b: X[a] @ X[b].T - X[b] @ X[a].T
    print(f"q = {q}, N = {N}: Hadamard; block pattern {pat}")
    rel = [(a, b, c, d) for a, b, c, d in ((1, 0, 3, 2), (2, 0, 1, 3), (3, 0, 2, 1))]
    signs = [('+' if np.array_equal(br(a, b), br(c, d)) else '-' if np.array_equal(br(a, b), -br(c, d)) else '?')
             for a, b, c, d in rel]
    print(f"  pairwise amicable: {amicable(X)};  [X_i, X_1] = {signs[0]}[X_k, X_j], [X_j, X_1] = {signs[1]}[X_i, X_k], "
          f"[X_k, X_1] = {signs[2]}[X_j, X_i]; [X_i, X_1] = 0: {not br(1, 0).any()}")


def rebase(q):
    P = Paley(q)
    N = P.N
    E = np.zeros((4, N, 4, 4, N), np.int64)
    for g, k, w, h, l in itertools.product(range(4), range(N), range(4), range(4), range(N)):
        E[g, k, w, h, l] = P.entry(P.G[g](P.reps[k]), P.G[w](P.G[h](P.reps[l])))
    found = tried = 0
    for gs in itertools.product(range(4), repeat=N):
        for hs in itertools.product(range(4), repeat=N):
            tried += 1
            X = [np.array([[E[gs[k], k, w, hs[l], l] for l in range(N)] for k in range(N)]) for w in range(4)]
            found += amicable(X)
    print(f"q = {q}, N = {N}: {found} of {tried} base-point choices give pairwise amicable blocks")


def plug(q, t):
    P = Paley(q)
    N = P.N
    X = blocks(P.matrix(), N)
    circ = lambda v: np.array([[v[(j - i) % len(v)] for j in range(len(v))] for i in range(len(v))], np.int64)
    T = None
    for assign in itertools.product(range(4), repeat=t):          # first circulant T-matrices of order t
        for signs in itertools.product((1, -1), repeat=t):
            rows = [np.zeros(t, np.int64) for _ in range(4)]
            for p, (a, s) in enumerate(zip(assign, signs)):
                rows[a][p] = s
            M = [circ(r) for r in rows]
            if np.array_equal(sum(m @ m.T for m in M), t * np.eye(t, dtype=np.int64)):
                T = M
                break
        if T:
            break
    R = np.eye(t, dtype=np.int64)[::-1]

    def bh(a, b, c, d):
        A1 = a * T[0] + b * T[1] + c * T[2] + d * T[3]
        A2 = -b * T[0] + a * T[1] + d * T[2] - c * T[3]
        A3 = -c * T[0] - d * T[1] + a * T[2] + b * T[3]
        A4 = -d * T[0] + c * T[1] - b * T[2] + a * T[3]
        return np.block([[A1, A2 @ R, A3 @ R, A4 @ R], [-A2 @ R, A1, A4.T @ R, -A3.T @ R],
                         [-A3 @ R, -A4.T @ R, A1, A2.T @ R], [-A4 @ R, A3.T @ R, -A2.T @ R, A1]])
    Mu = [bh(*[int(i == u) for i in range(4)]) for u in range(4)]
    n = 4 * t * N
    good = sum(np.array_equal(Hb @ Hb.T, n * np.eye(n, dtype=np.int64))
               for perm in itertools.permutations(range(4)) for sg in itertools.product((1, -1), repeat=4)
               for Hb in [sum(np.kron(Mu[u], sg[u] * X[perm[u]]) for u in range(4))])
    print(f"q = {q}, N = {N}, t = {t}: {good} of 384 substitutions give a Hadamard matrix of order {n}")


if __name__ == "__main__":
    cmd = sys.argv[1] if len(sys.argv) > 1 else ""
    if cmd == "pattern":
        pattern(int(sys.argv[2]))
    elif cmd == "rebase":
        rebase(int(sys.argv[2]))
    elif cmd == "plug":
        plug(int(sys.argv[2]), int(sys.argv[3]))
    else:
        print(__doc__)
