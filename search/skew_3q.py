#!/usr/bin/env python3
"""Skew-Hadamard matrices of order 3(q + 1) for primes q = 3 (mod 4), from a closed formula (independent of
search/paley_blocks.py, which builds the same family through PG(1, q) and also covers prime powers).

GF(q^2) = GF(q)[i] with i^2 = -1 (q = 3 mod 4), z a primitive element, a = z^(q+1) (a generator of GF(q)*, a
non-square).  The points of PG(1, q) are x_m = z^m, m = 0..q, and det(u, v) = Im(conj(u) v) with conj(u) = u^q, so
Paley's skew-Hadamard matrix of order q + 1 is

    P(x_m, x_n) = [m = n] + (-1)^m F(n - m),    F(d) = chi(Im z^d)   (F(0) = 0),

with F(d + q + 1) = -F(d) and F(-d) = (-1)^(d+1) F(d).  Multiplication by z^4 has the four orbits m = k (mod 4) of
length s = (q + 1)/4, so P is a 4 x 4 array of negacyclic s x s blocks N_kj(y) = [k = j, y = 0] + (-1)^k F(j - k + 4y)
(N(y + s) = -N(y)).  Frobenius (z -> z^q) is an anti-automorphism of P that fixes the orbits 0 and 2 and swaps 1 and 3;
it gives N_30(y) = N_10(1 - y) for every y.  With the base point of orbit 3 moved back by one step (N'_3j(y) =
N_3j(y + 1), N'_k3(y) = N_k3(y - 1)) and the T-matrices of order 3

    T_0 = I,  T_1 = Z,  T_2 = 0,  T_3 = -Z^2        (Z the cyclic shift of order 3),

A_j = sum_k T_k (x) N'_kj are +-1 matrices on Z_3 x N_s (circulant in the first factor, negacyclic in the second),
sum_j A_j A_j^T = 12 s I, and A_0 = I (x) N_00 + Z (x) N_10 - Z^2 (x) N_10^T is skew-type.  The Goethals-Seidel array
with A_0 first, R = R_3 (x) J_s (x -> -x on Z_3, the back-identity on N_s), is a skew-Hadamard matrix of order 12s.

    python3 search/skew_3q.py 11 23 331 599     # build and check H H^T = nI and H + H^T = 2I exactly

Checked for all 86 primes q = 3 (mod 4) from 7 to 991.
"""
import sys
import numpy as np


def is_prime(n):
    return n > 1 and all(n % p for p in range(2, int(n ** 0.5) + 1))


def F_table(q):
    """F(d) = chi(Im z^d) for d = 0..q (z the first primitive element a + bi in lexicographic order of (b, a))"""
    order = q * q - 1
    fac = [p for p in range(2, order + 1) if order % p == 0 and is_prime(p)]
    mul = lambda u, v: ((u[0] * v[0] - u[1] * v[1]) % q, (u[0] * v[1] + u[1] * v[0]) % q)
    def power(u, e):
        r = (1, 0)
        while e:
            if e & 1:
                r = mul(r, u)
            u, e = mul(u, u), e >> 1
        return r
    z = next((a, b) for b in range(1, q) for a in range(q) if all(power((a, b), order // p) != (1, 0) for p in fac))
    chi = np.full(q, -1, np.int64)
    chi[0] = 0
    chi[(np.arange(1, q) ** 2) % q] = 1
    F, w = np.zeros(q + 1, np.int64), (1, 0)
    for d in range(q + 1):
        F[d] = chi[w[1]]
        w = mul(w, z)
    return F, z


def build(q):
    assert is_prime(q) and q % 4 == 3 and q >= 7
    s = (q + 1) // 4
    F, z = F_table(q)
    Fx = lambda d: np.where((d // (q + 1)) % 2 == 1, -1, 1) * F[d % (q + 1)]      # antiperiodic extension
    y = np.arange(s)
    r = [0, 0, 0, -1]                                                             # base point of orbit 3 moved
    N = [[(-1) ** k * Fx(j - k + 4 * (y + r[j] - r[k])) + ((k == j) & (y == 0)) for j in range(4)] for k in range(4)]
    assert np.array_equal(N[3][0], N[1][0][(-y) % s] * np.where(y == 0, 1, -1))   # N'_30(y) = N_10(-y) (negacyclic)
    Z = np.roll(np.eye(3, dtype=np.int64), 1, axis=1)
    T = [np.eye(3, dtype=np.int64), Z, np.zeros((3, 3), np.int64), -(Z @ Z)]
    wrap = np.where(y[None, :] >= y[:, None], 1, -1)
    nega = lambda row: wrap * row[(y[None, :] - y[:, None]) % s]                  # negacyclic matrix, first row row
    A = [sum(np.kron(T[k], nega(N[k][j])) for k in range(4)) for j in range(4)]
    R = np.kron(np.eye(3, dtype=np.int64)[(-np.arange(3)) % 3], np.eye(s, dtype=np.int64)[::-1])
    A0, B, C, D = A
    H = np.block([[A0, B @ R, C @ R, D @ R],
                  [-B @ R, A0, D.T @ R, -C.T @ R],
                  [-C @ R, -D.T @ R, A0, B.T @ R],
                  [-D @ R, C.T @ R, -B.T @ R, A0]])
    return H, z


if __name__ == "__main__":
    for q in map(int, sys.argv[1:]):
        H, z = build(q)
        n = H.shape[0]
        Hf = H.astype(np.float64)                           # exact: integer sums of at most n < 2^53 terms +-1
        pm1 = bool(np.all(np.abs(H) == 1))
        had = pm1 and np.array_equal(Hf @ Hf.T, n * np.eye(n))
        skew = np.array_equal(H + H.T, 2 * np.eye(n, dtype=np.int64))
        print(f"q = {q}, z = {z[0]} + {z[1]}i: order {n} = 3(q + 1), Hadamard {had}, skew {skew}")
