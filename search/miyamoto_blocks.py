#!/usr/bin/env python3
"""Miyamoto's Hadamard matrix of order 4q as a 4 x 4 array of bordered blocks in a commutative algebra, multiplied by
T-matrices of order t: a Hadamard matrix of order 4qt (NOTES.md §23).

q = 1 (mod 4) is a prime, m = (q - 1)/2 and g a primitive root mod q. Paley's conference matrix C of order q + 1 is
indexed by infinity, 0, g^(2j) and g^(2j+1) (j in Z_m); its quadratic-residue blocks are then circulant in j. The
Hadamard matrix K of order q - 1 is of 2N-type, K = [[A, B], [-B^T, A^T]] with A, B the negacyclic matrices of a
negaperiodic Golay pair (a, b) of length m (Balonin and Djokovic, arXiv:1508.00640). Miyamoto's array (verify.py
`miyamoto`) then has sixteen blocks N_rc = s [[1, +-e], [+-e^T, J_2 (x) U + D_2 (x) V]] with U circulant and V negacyclic
over Z_m, D_2 = [[1, -1], [-1, 1]]. These blocks commute (the two Kronecker components are orthogonal idempotent
parts, each commutative), and the permutation R_B that fixes the border and sends (b, j) to (b, -j) for j = 0 and to
(1 - b, -j) for j != 0 satisfies R_B N R_B = N^T (on the negacyclic part it acts as diag(1, -1, ..., -1) times the
reversal). With T-matrices T_1..T_4 of order t, A_j = sum_k T_k (x) N_kj are commuting +-1 matrices of order tq with
R A_j R = A_j^T for R = R_t (x) R_B and sum_j A_j A_j^T = 4qt I, so the Goethals-Seidel array gives H(4qt).

Negaperiodic Golay pairs (NG pairs): (1, c_1..c_p), (-1, c_1..c_p) of length p + 1 from the negacyclic Paley
conference matrix of order p + 1 (p an odd prime power; their Prop. 6), (1, c_2, c_4, ...), (c_1, c_3, ...) of length
(p + 1)/2 for p = 3 (mod 4) (Prop. 8), from Williamson matrices of odd order w through their array (29) (length 2w),
and Turyn's multiplication by a Golay pair of length g (their Prop. 3, formulas (15) and (16)).

A periodic Golay pair (PG pair, PAF_a + PAF_b = 2m delta) works as well: then A, B and every V are circulant, and the
plain reversal (b, j) -> (b, -j) transposes the blocks. The PG pair of length 202 of Djokovic and Kotsireas (Numer.
Algorithms 69 (2015), arXiv:1310.5773, Section 3: subgroup H = {1, 87, 95, 137, 185} of Z_202^*) is `dk202`; the same
formulas (15), (16) multiply it by a Golay pair (Georgiou, Stylianou, Drosou and Koukouvinos 2014, Theorem 1).

    python3 search/miyamoto_blocks.py selftest
    python3 search/miyamoto_blocks.py certify q t SPEC [Tfile]   (SPEC: paley1:p | paley2:p | golay10:paley1:121 |
                                                                  golay2:dk202 | ...)
"""
import hashlib, itertools, json, sys
import numpy as np
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "verify"))
sys.path.insert(0, str(ROOT / "search"))
import verify as V


def factor(n):
    f, p = {}, 2
    while p * p <= n:
        while n % p == 0:
            f[p] = f.get(p, 0) + 1
            n //= p
        p += 1
    if n > 1:
        f[n] = f.get(n, 0) + 1
    return f


def prim_root(q):
    ps = list(factor(q - 1))
    return next(g for g in range(2, q) if all(pow(g, (q - 1) // p, q) != 1 for p in ps))


# ---------- negaperiodic Golay pairs ----------

def naf(x):
    """negaperiodic autocorrelations NAF_x(s), s = 0..v-1"""
    v, x = len(x), np.asarray(x, np.int64)
    return np.array([int(x[:v - s] @ x[s:]) - int(x[v - s:] @ x[:s]) for s in range(v)])


def is_ng(a, b):
    v = len(a)
    s = naf(a) + naf(b)
    return len(b) == v and s[0] == 2 * v and not s[1:].any()


def paf(x):
    v, x = len(x), np.asarray(x, np.int64)
    return np.array([int(x @ np.roll(x, -s)) for s in range(v)])


def is_pg(a, b):
    v = len(a)
    s = paf(a) + paf(b)
    return len(b) == v and s[0] == 2 * v and not s[1:].any()


def pg_dk202():
    """the periodic Golay pair of length 202 of Djokovic and Kotsireas (arXiv:1310.5773, Section 3), from the orbits
    of H = {1, 87, 95, 137, 185} on Z_202: -1 on the union of the orbits of the listed representatives"""
    v, H = 202, (1, 87, 95, 137, 185)
    reps = ((2, 4, 9, 11, 12, 13, 18, 20, 22, 24, 25, 26, 38, 41, 50, 51, 53, 55, 67, 76),
            (1, 3, 4, 6, 8, 9, 11, 12, 16, 17, 20, 25, 39, 41, 48, 52, 67, 76, 101))
    out = []
    for J in reps:
        x = np.ones(v, np.int64)
        x[sorted({h * j % v for h in H for j in J})] = -1
        out.append(x)
    assert is_pg(*out), "not a periodic Golay pair"
    return tuple(out)


def nega_conference_row(p):
    """first row (c_0 = 0, c_1, ..., c_p) of a negacyclic conference matrix of order p + 1, p an odd prime power:
    points of PG(1, p) as w^k F* (w primitive in E = GF(p^2), F = GF(p)), c_k = chi_F(D(1, w^k) / D(1, w)) with the
    alternating form D(x, y) = x y^p - x^p y; D(w^i, w^j) = N(w)^i D(1, w^(j-i)) and chi(N(w)) = -1 make
    [(-1)^i C_ij] negacyclic"""
    pk = factor(p)
    if len(pk) == 1 and list(pk.values())[0] == 1:
        # E = F(sqrt r): (x0, x1) = x0 + x1 sqrt r; Frobenius is conjugation, so D(1, z) = -2 z1 sqrt r
        r = next(x for x in range(2, p) if pow(x, (p - 1) // 2, p) == p - 1)
        mul = lambda x, y: ((x[0] * y[0] + r * x[1] * y[1]) % p, (x[0] * y[1] + x[1] * y[0]) % p)
        def power(x, e):
            out = (1, 0)
            while e:
                if e & 1: out = mul(out, x)
                x, e = mul(x, x), e >> 1
            return out
        order = p * p - 1
        qs = list(factor(order))
        w = next((x0, x1) for x0 in range(p) for x1 in range(1, p)
                 if all(power((x0, x1), order // l) != (1, 0) for l in qs))
        chi = lambda z: 0 if z % p == 0 else (1 if pow(z, (p - 1) // 2, p) == 1 else -1)
        z, zs = (1, 0), []
        for k in range(p + 1):
            zs.append(z[1])
            z = mul(z, w)
        inv1 = pow(zs[1], p - 2, p)
        return np.array([chi(x * inv1) for x in zs], np.int64)
    import paley_blocks as PB
    E = PB.GF(p * p)                      # w = x, the generator of E* in PB.GF's tables
    qE = p * p - 1
    def D(k):                             # log of D(1, w^k) = w^(kp) - w^k, or None if 0
        d = E.sub(int(E.exp[(k * p) % qE]), int(E.exp[k % qE]))
        return None if d == 0 else int(E.log[d])
    l1 = D(1)
    row = [0]
    for k in range(1, p + 1):
        l = (D(k) - l1) % qE              # the quotient lies in F* = <w^(p + 1)>
        assert l % (p + 1) == 0
        row.append(1 if (l // (p + 1)) % 2 == 0 else -1)
    return np.array(row, np.int64)


def negacyclic(x):
    v = len(x)
    M = np.empty((v, v), np.int64)
    for i in range(v):
        M[i] = np.concatenate([-x[v - i:], x[:v - i]])
    return M


def ng_paley1(p):
    c = nega_conference_row(p)
    C = negacyclic(c)
    assert np.array_equal(C @ C.T, p * np.eye(p + 1, dtype=np.int64)), "not a negacyclic conference matrix"
    a, b = c.copy(), c.copy()
    a[0], b[0] = 1, -1
    return a, b


def ng_paley2(p):
    assert p % 4 == 3
    c = nega_conference_row(p)
    n = p + 1
    a = np.concatenate([[1], c[2:n:2]])
    b = c[1:n:2].copy()
    return a, b


def golay_pairs():
    """Golay pairs of lengths 2, 10 and 26 (Golay 1961), checked"""
    P = {2: ("++", "+-"), 10: ("++-+-+--++", "++-+++++--"),
         26: ("++++-++--+-+-+--+-+++--+++", "++++-++--+-+++++-+---++---")}
    out = {}
    for g, (x, y) in P.items():
        a = np.array([1 if c == "+" else -1 for c in x]); b = np.array([1 if c == "+" else -1 for c in y])
        nac = lambda s, k: int(s[:len(s) - k] @ s[k:])
        assert all(nac(a, k) + nac(b, k) == 0 for k in range(1, g)), f"not a Golay pair of length {g}"
        out[g] = (a, b)
    return out


def ng_golay(gp, ng):
    """Turyn's multiplication GP_g x NG_v -> NG_gv (Balonin-Djokovic (15), (16))"""
    (a, b), (c, d) = gp, ng
    g, v = len(a), len(c)
    s, t = (a + b) // 2, (a - b) // 2
    e, f = np.zeros(g * v, np.int64), np.zeros(g * v, np.int64)
    for i in range(v):
        e[g * i:g * i + g] += s * c[i] + t * d[v - 1 - i]
        f[g * i:g * i + g] += -t * c[v - 1 - i] + s * d[i]
    return e, f


def williamson(w):
    """symmetric circulant Williamson matrices of odd order w (first rows), by a meet in the middle on PAF vectors"""
    h = (w - 1) // 2
    S = []
    for bits in itertools.product((1, -1), repeat=h + 1):
        r = np.empty(w, np.int64); r[0] = bits[0]
        for k in range(1, h + 1): r[k] = r[w - k] = bits[k]
        S.append(r)
    S = np.array(S)
    P = np.array([[int(r @ np.roll(r, -d)) for d in range(1, h + 1)] for r in S])
    rs = S.sum(1)
    for quad in sorted({tuple(sorted(abs(x) for x in c)) for c in itertools.product(range(1, w + 1, 2), repeat=4)
                        if sum(x * x for x in c) == 4 * w}):
        ra, rb, rc, rd = quad
        A, B, Cc, D = (np.where(rs == r)[0] for r in quad)
        left = {}
        for i in A:
            for j in B:
                left.setdefault((P[i] + P[j]).tobytes(), (i, j))
        for k in Cc:
            for l in D:
                key = (-(P[k] + P[l])).tobytes()
                if key in left:
                    i, j = left[key]
                    return [S[i], S[j], S[k], S[l]]
    return None


def ng_williamson(W):
    """NG pair of length 2w from Williamson (or quasi-Williamson) matrices A, B, C, D of odd order w: their array
    (29) is K = [[P, Q], [-Q^T, P^T]] with P = a + x b, Q = c + x d in Z[x]/(x^2 + 1) (x) Z[Z_w], and y -> x z is an
    isomorphism from Z[y]/(y^(2w) + 1) since w is odd: y^j = (-1)^(j // 2) x^(j mod 2) z^(j mod w). The coefficients
    of P and Q in the powers of y form the pair."""
    a, b, c, d = (np.asarray(x, np.int64) for x in W)
    w = len(a)
    j = np.arange(2 * w)
    sg = np.where((j // 2) % 2 == 0, 1, -1)
    p = sg * np.where(j % 2 == 0, a[j % w], b[j % w])
    q = sg * np.where(j % 2 == 0, c[j % w], d[j % w])
    return p, q


def ng_from_spec(spec):
    """paley1:p, paley2:p, williamson:w, or golayG:<spec> (negaperiodic)"""
    kind, _, rest = spec.partition(":")
    if kind == "paley1": return ng_paley1(int(rest))
    if kind == "paley2": return ng_paley2(int(rest))
    if kind == "williamson": return ng_williamson(williamson(int(rest)))
    if kind.startswith("golay"): return ng_golay(golay_pairs()[int(kind[5:])], ng_from_spec(rest))
    raise ValueError(spec)


def pair_from_spec(spec):
    """(a, b, periodic): dk202 and golayG:...:dk202 are periodic Golay pairs, the rest negaperiodic"""
    if spec.split(":")[-1] == "dk202":
        a, b = pg_dk202()
        for kind in reversed(spec.split(":")[:-1]):
            a, b = ng_golay(golay_pairs()[int(kind[5:])], (a, b))
        assert is_pg(a, b)
        return a, b, True
    a, b = ng_from_spec(spec)
    return a, b, False


# ---------- Miyamoto's array in sixteen blocks, times T-matrices ----------

def paley_cyclo(q, g):
    """Paley's symmetric conference matrix of order q + 1 on infinity, 0, g^(2j), g^(2j+1) (j < m), int8"""
    m = (q - 1) // 2
    chi = -np.ones(q, np.int64); chi[0] = 0; chi[(np.arange(1, q) ** 2) % q] = 1
    pts = np.array([0] + [pow(g, 2 * j, q) for j in range(m)] + [pow(g, 2 * j + 1, q) for j in range(m)], np.int64)
    C = np.ones((q + 1, q + 1), np.int8)
    C[0, 0] = 0
    C[1:, 1:] = chi[(pts[None, :] - pts[:, None]) % q]
    return C


def reversal_b(q, periodic=False):
    m = (q - 1) // 2
    R = np.arange(q)
    for b in range(2):
        for j in range(m):
            R[1 + b * m + j] = 1 + (b if (j == 0 or periodic) else 1 - b) * m + (-j) % m
    return R


def circulant(x):
    x = np.asarray(x, np.int64)
    return np.array([np.roll(x, i) for i in range(len(x))])


def build(q, g, a, b, T, check_algebra=False, periodic=False):
    """H(4qt) from q, the primitive root g, an NG pair (or, if periodic, a PG pair) (a, b) of length (q - 1)/2 and
    T-matrices T (4 x t), int8"""
    m = (q - 1) // 2
    assert len(a) == m and (is_pg(a, b) if periodic else is_ng(a, b))
    mat = circulant if periodic else negacyclic
    A, B = mat(np.asarray(a)), mat(np.asarray(b))
    K = np.block([[A, B], [-B.T, A.T]]).astype(np.int8)
    H = V.miyamoto(paley_cyclo(q, g), K)
    N = [[H[r * q:(r + 1) * q, c * q:(c + 1) * q] for c in range(4)] for r in range(4)]
    Rb = reversal_b(q, periodic)
    if check_algebra:
        Ns = [x.astype(np.int64) for row in N for x in row]
        assert all(np.array_equal(x[np.ix_(Rb, Rb)], x.T) for x in Ns), "R N R != N^T"
        assert all(np.array_equal(x @ y, y @ x) for x in Ns for y in Ns), "blocks do not commute"
    T = np.asarray(T, np.int64)
    t = T.shape[1]
    Ct = [V.circulant(r).astype(np.int8) for r in T]
    Aj = [sum(np.kron(Ct[k], N[k][j]) for k in range(4)).astype(np.int8) for j in range(4)]
    R = (((-np.arange(t)) % t)[:, None] * q + Rb[None, :]).reshape(-1)
    n = t * q
    G = np.empty((4 * n, 4 * n), np.int8)
    AR = [x[:, R] for x in Aj]
    ATR = [x.T[:, R] for x in Aj]
    rows = [[Aj[0], AR[1], AR[2], AR[3]], [-AR[1], Aj[0], ATR[3], -ATR[2]],
            [-AR[2], -ATR[3], Aj[0], ATR[1]], [-AR[3], ATR[2], -ATR[1], Aj[0]]]
    for i in range(4):
        for j in range(4):
            G[i * n:(i + 1) * n, j * n:(j + 1) * n] = rows[i][j]
    return G


def load_t(t, path=None):
    if path:
        rows = [l.strip() for l in open(path) if l.strip() and not l.startswith("#")]
        return np.array([[{"0": 0, "+": 1, "-": -1}[c] for c in r] for r in rows], np.int64)
    import paley_blocks as PB
    return PB.t_matrices(t)


def selftest():
    ok = True
    for p in (3, 5, 7, 9, 11, 13, 19, 23, 25, 27, 43, 47, 121):
        a, b = ng_paley1(p); ok &= is_ng(a, b)
        if p % 4 == 3:
            a, b = ng_paley2(p); ok &= is_ng(a, b)
    a, b = ng_williamson(williamson(9)); ok &= is_ng(a, b)
    a, b = ng_from_spec("golay10:paley2:11"); ok &= is_ng(a, b) and len(a) == 60
    a, b = ng_from_spec("golay26:paley1:3"); ok &= is_ng(a, b) and len(a) == 104
    a, b, _ = pair_from_spec("golay2:dk202"); ok &= is_pg(a, b) and len(a) == 404
    print("NG and PG constructions:", ok)
    for q, spec in ((13, "paley2:11"), (29, "paley1:13"), (37, "williamson:9"), (41, "golay2:paley2:19"),
                    (53, "williamson:13"), (61, "paley2:59")):
        g = prim_root(q)
        a, b = ng_from_spec(spec)
        for t in (3, 5, 7):
            G = build(q, g, a, b, load_t(t), check_algebra=(t == 3)).astype(np.float64)
            good = np.array_equal(G @ G.T, 4 * q * t * np.eye(4 * q * t))
            ok &= good
            print(f"q = {q} ({spec}), t = {t}: order {4 * q * t}, H H^T = nI: {good}")
    # periodic Golay pairs: Golay pairs of length 2 and 26 are periodic Golay pairs (q = 5, 53)
    for q, gl in ((5, 2), (53, 26)):
        a, b = golay_pairs()[gl]
        G = build(q, prim_root(q), a, b, load_t(3), check_algebra=True, periodic=True).astype(np.float64)
        good = np.array_equal(G @ G.T, 4 * q * 3 * np.eye(4 * q * 3))
        ok &= good
        print(f"q = {q} (periodic Golay pair of length {gl}), t = 3: order {12 * q}, H H^T = nI: {good}")
    return ok


def certify(q, t, spec, tfile=None):
    g = prim_root(q)
    a, b, periodic = pair_from_spec(spec)
    T = load_t(t, tfile)
    assert T.shape[1] == t
    G = build(q, g, a, b, T, periodic=periodic)
    sha = hashlib.sha256()
    for chunk in V.csv_chunks(G):
        sha.update(chunk)
    sgn = lambda x: "".join("+" if v > 0 else "-" for v in x)
    tsg = lambda x: "".join({1: "+", -1: "-", 0: "0"}[int(v)] for v in x)
    cert = {"order": 4 * q * t, "construction": "miyamoto-blocks", "v": q * t, "q": q, "g": g, "t": t,
            "ng_source": spec, ("pg" if periodic else "ng"): [sgn(a), sgn(b)], "T": [tsg(r) for r in T],
            "sha256": sha.hexdigest()}
    d = ROOT / "hadamard" / str(4 * q * t)
    d.mkdir(exist_ok=True)
    (d / "certificate.json").write_text(json.dumps(cert, indent=1) + "\n")
    print(f"wrote {d}/certificate.json (order {4 * q * t}, sha256 {cert['sha256']})")


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "selftest":
        raise SystemExit(0 if selftest() else 1)
    if len(sys.argv) > 4 and sys.argv[1] == "certify":
        certify(int(sys.argv[2]), int(sys.argv[3]), sys.argv[4], sys.argv[5] if len(sys.argv) > 5 else None)
        raise SystemExit(0)
    print(__doc__)
