#!/usr/bin/env python3
"""Hadamard matrices of order 4st from Paley's matrix of order q + 1 = 4s (q = 3 mod 8 prime) and T-matrices of order t.

1. Paley's skew-Hadamard matrix H = I + S on PG(1, q) is invariant, up to signs, under PSL(2, q):
   H(gx, gy) = e(x) e(y) H(x, y).  An element g of order s = (q + 1)/4 (odd) in a non-split torus has no fixed point,
   so PG(1, q) splits into four g-orbits of length s.  Rescaling the points by signs E with E(gx) = c e(x) E(x) (possible
   because s is odd) makes g a plain automorphism, and ordering each orbit along g turns H into a 4 x 4 array of
   circulant s x s blocks N_kj.
2. Seberry-Yamada (1992), Theorem 3.10: with T-matrices T_1, ..., T_4 of order t (circulant, disjoint supports,
   sum T_k T_k^T = t I) the matrices A_j = sum_k T_k (x) N_kj are +-1 and group-developed over Z_t x Z_s, and
   sum_j A_j A_j^T = sum_{k,l} T_k T_l^T (x) sum_j N_kj N_lj^T = 4st I by the orthogonality of the block rows of H.
   For gcd(s, t) = 1 they are circulants of length n = st, a Goethals-Seidel quadruple, and the Goethals-Seidel array
   gives a Hadamard matrix of order 4n.

    python3 search/paley_blocks.py check 1051              # the 16 circulant blocks of Paley's matrix of order 1052
    python3 search/paley_blocks.py quad 1051 3 out.json    # GS quadruple of length 789 (T-matrices of order 3)
    python3 search/paley_blocks.py certify 1051 3          # hadamard/3156/certificate.json and sequences.txt
"""
import json
import sys
import numpy as np


def is_prime(n):
    return n > 1 and all(n % p for p in range(2, int(n ** 0.5) + 1))


def paley(q):
    """Paley's skew-Hadamard matrix of order q + 1 on the points 0..q-1, q (= infinity)"""
    chi = np.array([0] + [1 if pow(x, (q - 1) // 2, q) == 1 else -1 for x in range(1, q)], np.int64)
    x = np.arange(q)
    S = np.zeros((q + 1, q + 1), np.int64)
    S[:q, :q] = chi[(x[None, :] - x[:, None]) % q]
    S[q, :q], S[:q, q] = 1, -1
    return S + np.eye(q + 1, dtype=np.int64)


def mobius(q, a, b, c, d):
    """the permutation x -> (a x + b)/(c x + d) of PG(1, q), infinity = q"""
    p = np.empty(q + 1, np.int64)
    for x in range(q):
        den = (c * x + d) % q
        p[x] = q if den == 0 else (a * x + b) * pow(den, -1, q) % q
    p[q] = q if c % q == 0 else a * pow(c, -1, q) % q
    return p


def cycles(p):
    seen, out = np.zeros(len(p), bool), []
    for i in range(len(p)):
        if not seen[i]:
            cyc, j = [], i
            while not seen[j]:
                seen[j] = True
                cyc.append(j)
                j = p[j]
            out.append(cyc)
    return out


def torus_element(q):
    """a fixed-point-free element of PSL(2, q) of order s = (q + 1)/4, from x -> -1/(x + tau)"""
    s = (q + 1) // 4
    for tau in range(1, q):
        disc = (tau * tau - 4) % q
        if pow(disc, (q - 1) // 2, q) != q - 1:          # x^2 - tau x + 1 must be irreducible
            continue
        p = mobius(q, 0, q - 1, 1, tau)
        for k in (1, 2):
            pk = p
            for _ in range(k - 1):
                pk = pk[p]
            if {len(c) for c in cycles(pk)} == {s}:
                return pk, tau, k
    raise ValueError("no torus element of order s")


def blocks(q):
    """the 4 x 4 array of circulant s x s blocks N[k][j] (first rows) and the full rearranged matrix"""
    assert is_prime(q) and q % 8 == 3
    s, n = (q + 1) // 4, q + 1
    H = paley(q)
    g, tau, kpow = torus_element(q)
    Hg = H[np.ix_(g, g)]                                  # Hg[x, y] = H(gx, gy) = e(x) e(y) H(x, y)
    e = H[0] * Hg[0]                                      # e(0) e(y), with e(0) := 1
    assert np.array_equal(e[:, None] * e[None, :] * H, Hg), "g is not a signed automorphism"
    cyc = cycles(g)
    for c in (1, -1):
        E, ok = np.zeros(n, np.int64), True
        for cy in cyc:
            E[cy[0]] = 1
            for a in range(1, s):
                E[cy[a]] = c * e[cy[a - 1]] * E[cy[a - 1]]
            ok &= c * e[cy[-1]] * E[cy[-1]] == E[cy[0]]
        if ok:
            break
    assert ok, "sign normalisation failed"
    order = [x for cy in cyc for x in cy]
    Hb = (E[:, None] * H * E[None, :])[np.ix_(order, order)]
    N = [[Hb[k * s:(k + 1) * s, j * s:(j + 1) * s] for j in range(4)] for k in range(4)]
    for k in range(4):
        for j in range(4):
            r = N[k][j][0]
            assert all(np.array_equal(N[k][j][a], np.roll(r, a)) for a in range(s)), f"block {k},{j} not circulant"
    F = Hb.astype(np.float32)                              # exact: integer sums of at most n < 2^24 terms +-1
    assert np.array_equal(F @ F.T, n * np.eye(n, dtype=np.float32))
    return [[N[k][j][0] for j in range(4)] for k in range(4)], Hb, {"tau": tau, "power": kpow, "c": c}


def base_sequences(m):
    """base sequences BS(m + 1, m): A, B of length m + 1 and C, D of length m with N_A + N_B + N_C + N_D = 0 off 0
    (aperiodic autocorrelations), by a meet in the middle; the first solution in lexicographic order"""
    import itertools
    def nac(x):
        return tuple(int(np.dot(x[:-k], x[k:])) for k in range(1, len(x)))
    seqs = lambda L: [np.array(b) for b in itertools.product((1, -1), repeat=L)]
    S1, S0 = seqs(m + 1), seqs(m)
    need = {}
    for c in S0:
        for d in S0:
            v = tuple(a + b for a, b in zip(nac(c), nac(d))) + (0,) if m > 0 else (0,)
            need.setdefault(v, (c, d))
    for a in S1:
        for b in S1:
            v = tuple(-(x + y) for x, y in zip(nac(a), nac(b)))
            key = v[:-1] + (0,) if m > 0 else (0,)
            if v[-1] == 0 and key in need:
                c, d = need[key]
                return a, b, c, d
    raise ValueError("no base sequences")


def t_matrices(t):
    """first rows of T-matrices of order t = 2m + 1 from base sequences BS(m + 1, m) (Turyn):
    ((A + B)/2; 0), ((A - B)/2; 0), (0; (C + D)/2), (0; (C - D)/2)"""
    m = (t - 1) // 2
    if m == 0:
        return np.array([[1], [0], [0], [0]], np.int64)
    a, b, c, d = base_sequences(m)
    z1, z0 = np.zeros(m + 1, np.int64), np.zeros(m, np.int64)
    T = np.array([np.concatenate([(a + b) // 2, z0]), np.concatenate([(a - b) // 2, z0]),
                  np.concatenate([z1, (c + d) // 2]), np.concatenate([z1, (c - d) // 2])], np.int64)
    circ = lambda v: np.array([np.roll(v, k) for k in range(t)])
    assert (np.abs(T).sum(0) == 1).all()
    assert np.array_equal(sum(circ(r) @ circ(r).T for r in T), t * np.eye(t, dtype=np.int64))
    return T


def quadruple(q, T):
    """GS quadruple of length n = s t (gcd(s, t) = 1): A_j(x) = sum_k T_k(x mod t) N_kj(x mod s) via CRT"""
    Nrow, _, info = blocks(q)
    s, t = (q + 1) // 4, T.shape[1]
    assert np.gcd(s, t) == 1
    n = s * t
    x = np.arange(n)
    A = [sum(T[k][x % t] * Nrow[k][j][x % s] for k in range(4)) for j in range(4)]
    paf = sum(np.array([int(a @ np.roll(a, -d)) for d in range(n)]) for a in A)
    assert set(np.unique(np.concatenate(A))) <= {-1, 1}
    assert paf[0] == 4 * n and not paf[1:].any(), "PAF sum is not zero"
    return A, info


if __name__ == "__main__":
    cmd = sys.argv[1] if len(sys.argv) > 1 else ""
    if cmd == "check":
        q = int(sys.argv[2])
        Nrow, Hb, info = blocks(q)
        print(f"q = {q}: Paley's matrix of order {q + 1} is a 4 x 4 array of circulant blocks of order {(q + 1) // 4}")
    elif cmd == "quad":
        q, t, out = int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
        T = t_matrices(t)
        A, info = quadruple(q, T)
        n = len(A[0])
        json.dump({"n": n, "q": q, "t": t, **info, "T": ["".join("0+-"[int(v)] for v in r) for r in T],
                   "sequences": ["".join("+" if v > 0 else "-" for v in a) for a in A]}, open(out, "w"))
        print(f"q = {q}, t = {t}: Goethals-Seidel quadruple of length {n} (PAF sum 0 off 0), order {4 * n}")
    elif cmd == "certify":
        import hashlib
        from pathlib import Path
        root = Path(__file__).resolve().parent.parent
        sys.path.insert(0, str(root / "verify"))
        import verify as V
        q, t = int(sys.argv[2]), int(sys.argv[3])
        T = t_matrices(t)
        A, info = quadruple(q, T)
        n = len(A[0])
        cert = {"order": 4 * n, "construction": "paley-blocks", "v": n, "q": q, "t": t, **{k: int(v) for k, v in info.items()},
                "T": ["".join("0+-"[int(v)] for v in r) for r in T],
                "sequences": ["".join("+" if v > 0 else "-" for v in a) for a in A]}
        H, _ = V.build_paley_blocks(cert)
        sha = hashlib.sha256()
        for chunk in V.csv_chunks(H):
            sha.update(chunk)
        cert["sha256"] = sha.hexdigest()
        d = root / "hadamard" / str(4 * n)
        d.mkdir(exist_ok=True)
        (d / "certificate.json").write_text(json.dumps(cert, indent=1) + "\n")
        (d / "sequences.txt").write_text("\n".join(cert["sequences"]) + "\n")
        print(f"wrote {d}/certificate.json (order {4 * n}, sha256 {cert['sha256']})")
    else:
        print(__doc__)
