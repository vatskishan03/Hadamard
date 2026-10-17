#!/usr/bin/env python3
"""Hadamard matrices of order 4st from Paley's matrix of order q + 1 = 4s (q = 3 mod 8 a prime or a prime power) and
T-matrices of order t.

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
    python3 search/paley_blocks.py certify 523 29 T29.txt  # with T-matrices read from a file (four rows over 0, +, -)
    python3 search/paley_blocks.py certify 243 151 t-matrices/T151.txt   # q = 3^5: GF(q) with a primitive polynomial

For gcd(s, t) > 1 the sequences live on Z_t x Z_s (certificate field "group").  For a prime power q = p^k the field
is GF(p)[x]/(f) with the first primitive f (certificate field "field"); points are coded as integers sum d_i p^i.

3. Skew-Hadamard matrices of order 3(q + 1).  A reflection w of the dihedral normaliser of the torus in PGL(2, q) with
   non-square determinant is an anti-automorphism of H (H(wx, wy) = f(x) f(y) H(y, x)) that inverts g, fixes two of the
   four g-orbits (c and another) and swaps the other two (b and d).  Hence N_dc(y) = eps N_bc(delta - y), with a
   constant sign eps because s is odd.  Moving the base point of orbit d by -delta (N'_kj(y) = N_kj(y + r_j - r_k),
   r_d = -delta) gives N'_dc = eps N'_bc^T, and with the T-matrices T_c = I, T_b = Z, T_d = -eps Z^2 (Z the cyclic
   shift of order 3) and 0 on the fourth orbit, A_c = I (x) N'_cc + Z (x) N'_bc - Z^2 (x) N'_bc^T is skew-type.  The
   Goethals-Seidel array with A_c first is a skew-Hadamard matrix of order 12s = 3(q + 1).  For q = 7 (mod 8), s is even,
   g^s acts as -I and the blocks are negacyclic; the same steps work over Z_3 x N_s (circulant in the Z_3 part,
   negacyclic in the other), with R: (i, a) -> (-i, s - 1 - a).  So there is a skew-Hadamard matrix of order 3(q + 1)
   for every prime power q = 3 (mod 4).

    python3 search/paley_blocks.py skew 331                # skew/996/certificate.json and sequences.txt
    python3 search/paley_blocks.py skew 27 /tmp/x          # the same in another directory (here q = 3^3, order 84)
    python3 search/paley_blocks.py skew 599                # q = 7 (mod 8): skew/1800, negacyclic blocks
"""
import json
import sys
import numpy as np


def is_prime(n):
    return n > 1 and all(n % p for p in range(2, int(n ** 0.5) + 1))


def prime_power(q):
    for p in range(2, q + 1):
        if q % p == 0:
            k, r = 0, q
            while r % p == 0:
                r //= p
                k += 1
            return (p, k) if r == 1 else None
    return None


class GF:
    """GF(p^k) for k > 1: an element is the integer sum d_i p^i of its coefficients d_i in the basis 1, x, ..., x^(k-1)
    modulo the first monic primitive polynomial x^k + f_(k-1) x^(k-1) + ... + f_0 in lexicographic order of
    (f_0, ..., f_(k-1)); exp[i] = x^i and log are the tables of the cyclic group GF(q)*."""
    def __init__(self, q):
        self.p, self.k = prime_power(q)
        self.q, p, k = q, self.p, self.k
        self.pw = p ** np.arange(k)
        import itertools
        for f in itertools.product(range(p), repeat=k):
            if f[0] == 0:
                continue
            exp, d = [1], [1] + [0] * (k - 1)
            for _ in range(q - 2):
                top = d[-1]
                d = [0] + d[:-1]
                d = [(d[i] - top * f[i]) % p for i in range(k)]
                exp.append(sum(c * p ** i for i, c in enumerate(d)))
            if len(set(exp)) == q - 1 and 1 not in exp[1:]:
                self.poly = list(f)
                break
        self.exp = np.array(exp, np.int64)
        self.log = np.full(q, -1, np.int64)
        self.log[self.exp] = np.arange(q - 1)
        self.digits = (np.arange(q)[:, None] // self.pw[None, :]) % p

    def sub(self, a, b):
        return (((self.digits[a] - self.digits[b]) % self.p) * self.pw).sum(-1)

    def add(self, a, b):
        return (((self.digits[a] + self.digits[b]) % self.p) * self.pw).sum(-1)

    def neg(self, a):
        return self.sub(0, a)

    def mul(self, a, b):
        return 0 if a == 0 or b == 0 else int(self.exp[(self.log[a] + self.log[b]) % (self.q - 1)])

    def inv(self, a):
        return int(self.exp[(-self.log[a]) % (self.q - 1)])

    def chi(self):
        c = np.where(self.log % 2 == 0, 1, -1).astype(np.int8)
        c[0] = 0
        return c


def paley(q):
    """Paley's skew-Hadamard matrix of order q + 1 on the points 0..q-1, q (= infinity), as int8"""
    if not is_prime(q):
        F = GF(q)
        chi = F.chi()
        S = np.zeros((q + 1, q + 1), np.int8)
        x = np.arange(q)
        for a in range(0, q, 256):
            b = min(a + 256, q)
            S[a:b, :q] = chi[F.sub(x[None, :], x[a:b, None])]
        S[q, :q], S[:q, q] = 1, -1
        S[np.arange(q + 1), np.arange(q + 1)] += 1
        return S
    chi = np.full(q, -1, np.int8)
    chi[0] = 0
    chi[(np.arange(1, q, dtype=np.int64) ** 2) % q] = 1
    x = np.arange(q)
    S = np.zeros((q + 1, q + 1), np.int8)
    for a in range(0, q, 1024):
        b = min(a + 1024, q)
        S[a:b, :q] = chi[(x[None, :] - x[a:b, None]) % q]
    S[q, :q], S[:q, q] = 1, -1
    S[np.arange(q + 1), np.arange(q + 1)] += 1
    return S


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
    if not is_prime(q):
        F = GF(q)
        chi, four = F.chi(), F.add(F.add(1, 1), F.add(1, 1))
        for tau in range(1, q):
            if chi[F.sub(F.mul(tau, tau), four)] != -1:          # x^2 - tau x + 1 must be irreducible
                continue
            p = np.empty(q + 1, np.int64)
            for x in range(q):
                den = F.add(x, tau)
                p[x] = q if den == 0 else F.neg(F.inv(den))
            p[q] = 0
            for k in (1, 2):
                pk = p if k == 1 else p[p]
                if {len(c) for c in cycles(pk)} == {s}:
                    return pk, tau, k
        raise ValueError("no torus element of order s")
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


def nega_ext(row, y, nega):
    """row~(y) for any integer y: the periodic extension, or the antiperiodic one (row~(y + s) = -row~(y)) if nega"""
    s = len(row)
    return (np.where((y // s) % 2 == 1, -1, 1) if nega else 1) * row[y % s]


def blocks(q):
    """the 4 x 4 array of circulant s x s blocks N[k][j] (first rows) and the full rearranged matrix.
    For q = 7 (mod 8), s is even, every lift of g to SL(2, q) has order 2s, so g^s acts as -I: the signs close with -1
    on every cycle and the blocks are negacyclic, N_kj[i, j] = N~_kj(j - i) with N~(y + s) = -N~(y)
    (info["negacyclic"])."""
    assert prime_power(q) and q % 4 == 3 and q >= 7
    s, n = (q + 1) // 4, q + 1
    nega = s % 2 == 0
    H = paley(q)
    g, tau, kpow = torus_element(q)
    e = H[0] * H[g[0], g]                                 # e(0) e(y), with e(0) := 1, from H(gx, gy) = e(x) e(y) H(x, y)
    for a in range(0, n, 1024):
        assert np.array_equal(e[a:a + 1024, None] * e[None, :] * H[a:a + 1024], H[np.ix_(g[a:a + 1024], g)]), \
            "g is not a signed automorphism"
    e = e.astype(np.int64)
    cyc = cycles(g)
    for c in (1, -1):
        E, ok = np.zeros(n, np.int64), True
        for cy in cyc:
            E[cy[0]] = 1
            for a in range(1, s):
                E[cy[a]] = c * e[cy[a - 1]] * E[cy[a - 1]]
            ok &= c * e[cy[-1]] * E[cy[-1]] == (-1 if nega else 1) * E[cy[0]]
        if ok:
            break
    assert ok, "sign normalisation failed"
    order = np.array([x for cy in cyc for x in cy])
    E8 = E.astype(np.int8)
    Hb = E8[order, None] * H[np.ix_(order, order)] * E8[None, order]
    del H
    N = [[Hb[k * s:(k + 1) * s, j * s:(j + 1) * s] for j in range(4)] for k in range(4)]
    idx = np.arange(s)
    wrap = np.where(idx[None, :] >= idx[:, None], 1, -1) if nega else 1
    for k in range(4):
        for j in range(4):
            r = N[k][j][0]
            assert np.array_equal(N[k][j], wrap * r[(idx[None, :] - idx[:, None]) % s]), f"block {k},{j} not circulant"
    F = Hb.astype(np.float32)                              # exact: integer sums of at most n < 2^24 terms +-1
    assert np.array_equal(F @ F.T, n * np.eye(n, dtype=np.float32))
    info = {"tau": tau, "power": kpow, "c": c, **({"negacyclic": True} if nega else {})}
    return [[N[k][j][0].astype(np.int64) for j in range(4)] for k in range(4)], Hb, info


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


def load_t(path):
    """T-matrices from a file of four first rows over 0, +, - (as in t-matrices/), checked"""
    rows = [l.strip() for l in open(path) if l.strip() and not l.startswith('#')]
    T = np.array([[{'0': 0, '+': 1, '-': -1}[ch] for ch in r] for r in rows], np.int64)
    t = T.shape[1]
    assert T.shape[0] == 4 and (np.abs(T).sum(0) == 1).all()
    paf = sum(np.array([int(r @ np.roll(r, -d)) for d in range(t)]) for r in T)
    assert paf[0] == t and not paf[1:].any(), "not T-matrices"
    return T


def skew_relation(Nrow, s, nega=False):
    """the first block column c (in order) with two other block rows b, d, a shift delta and a sign eps such that
    N_dc(y) = eps N_bc(delta - y) (N~ antiperiodic if nega)"""
    for c in range(4):
        for b in range(4):
            for d in range(4):
                if len({b, c, d}) < 3:
                    continue
                x, y = Nrow[d][c], Nrow[b][c]
                for delta in range(2 * s if nega else s):
                    for eps in (1, -1):
                        if np.array_equal(x, eps * nega_ext(y, delta - np.arange(s), nega)):
                            return c, b, d, delta, eps
    raise ValueError("no block column with related blocks")


def skew_data(q):
    """the re-based blocks (shifts r), the column order (c first) and the T-matrices of order 3 that make A_c
    skew-type"""
    Nrow, _, info = blocks(q)
    s, nega = (q + 1) // 4, "negacyclic" in info
    c, b, d, delta, eps = skew_relation(Nrow, s, nega)
    r = [0, 0, 0, 0]
    r[d] = (-delta) % (2 * s if nega else s)
    T = np.zeros((4, 3), np.int64)
    T[c, 0], T[b, 1], T[d, 2] = 1, 1, -eps
    cols = [c] + [j for j in range(4) if j != c]
    return Nrow, info, r, cols, T


def quadruple(q, T, Nrow=None, info=None, shifts=(0, 0, 0, 0), cols=(0, 1, 2, 3)):
    """GS quadruple of length n = s t: A_j(x) = sum_k T_k(x mod t) N_kj(x mod s) via CRT when gcd(s, t) = 1; otherwise
    A_j(i, a) = sum_k T_k(i) N_kj(a) on Z_t x Z_s, flattened as i s + a (always so for negacyclic blocks, Z_t x N_s).
    The blocks are first re-based, N_kj(y) -> N_kj(y + r_j - r_k), and the block columns taken in the order cols"""
    if Nrow is None:
        Nrow, _, info = blocks(q)
    s, t = (q + 1) // 4, T.shape[1]
    nega, r = "negacyclic" in info, shifts
    Nrow = [[nega_ext(Nrow[k][j], np.arange(s) + r[j] - r[k], nega) for j in cols] for k in range(4)]
    n = s * t
    if np.gcd(s, t) == 1 and not nega:
        x = np.arange(n)
        A = [sum(T[k][x % t] * Nrow[k][j][x % s] for k in range(4)) for j in range(4)]
        paf = sum(np.array([int(a @ np.roll(a, -d)) for d in range(n)]) for a in A)
    else:
        A2 = [sum(np.outer(T[k], Nrow[k][j]) for k in range(4)) for j in range(4)]
        sh = lambda a, da: np.concatenate([a, -a if nega else a], axis=1)[:, da:da + s]
        paf = sum(np.array([[int((a * np.roll(sh(a, da), -di, 0)).sum()) for da in range(s)] for di in range(t)])
                  for a in A2).reshape(-1)
        A = [a.reshape(-1) for a in A2]
        info = {**info, "group": f"Z{t} x {'N' if nega else 'Z'}{s}"}
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
        T = load_t(sys.argv[4]) if len(sys.argv) > 4 else t_matrices(t)
        assert T.shape[1] == t
        A, info = quadruple(q, T)
        n = len(A[0])
        cert = {"order": 4 * n, "construction": "paley-blocks", "v": n, "q": q, "t": t,
                **{k: (v if isinstance(v, str) else int(v)) for k, v in info.items()},
                "T": ["".join("0+-"[int(v)] for v in r) for r in T],
                "sequences": ["".join("+" if v > 0 else "-" for v in a) for a in A]}
        if not is_prime(q):
            F = GF(q)
            cert = {**{k: v for k, v in cert.items() if k not in ("T", "sequences")},
                    "field": {"p": F.p, "k": F.k, "poly": F.poly}, "T": cert["T"], "sequences": cert["sequences"]}
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
    elif cmd == "skew":
        import hashlib
        from pathlib import Path
        root = Path(__file__).resolve().parent.parent
        sys.path.insert(0, str(root / "verify"))
        import verify as V
        q = int(sys.argv[2])
        Nrow, info, r, cols, T = skew_data(q)
        A, info = quadruple(q, T, Nrow, info, r, cols)
        n = len(A[0])
        cert = {"order": 4 * n, "construction": "paley-blocks", "skew": True, "v": n, "q": q, "t": 3,
                **{k: (v if isinstance(v, (str, bool)) else int(v)) for k, v in info.items()},
                "shifts": [int(x) for x in r], "columns": cols,
                "T": ["".join("0+-"[int(v)] for v in row) for row in T],
                "sequences": ["".join("+" if v > 0 else "-" for v in a) for a in A]}
        if not is_prime(q):
            F = GF(q)
            cert = {**{k: v for k, v in cert.items() if k not in ("T", "sequences")},
                    "field": {"p": F.p, "k": F.k, "poly": F.poly}, "T": cert["T"], "sequences": cert["sequences"]}
        H, _ = V.build_paley_blocks(cert)
        assert V.is_skew(H), "H + H^T != 2I"
        sha = hashlib.sha256()
        for chunk in V.csv_chunks(H):
            sha.update(chunk)
        cert["sha256"] = sha.hexdigest()
        d = Path(sys.argv[3]) if len(sys.argv) > 3 else root / "skew" / str(4 * n)
        if (d / "certificate.json").exists():
            sys.exit(f"{d}/certificate.json exists; not overwritten")
        d.mkdir(parents=True, exist_ok=True)
        (d / "certificate.json").write_text(json.dumps(cert, indent=1) + "\n")
        (d / "sequences.txt").write_text("\n".join(cert["sequences"]) + "\n")
        print(f"wrote {d}/certificate.json (skew order {4 * n}, sha256 {cert['sha256']})")
    else:
        print(__doc__)
