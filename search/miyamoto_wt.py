#!/usr/bin/env python3
"""Williamson-type matrices of order q, and a complex Hadamard matrix of order 2q, from a symmetric commuting pair of
order n = (q - 1)/2 (NOTES.md §19).

Seberry and Yamada's Lemma 8.6 (after Miyamoto): if U_i, V_i (i = 1..4) are (0, +-1) matrices of order n, the U_i
pairwise amicable, the V_i pairwise amicable, U_i +- V_i (+-1) matrices, U_1 of row sum 1 and U_2, U_3, U_4 of row sum
0, sum U_i U_i^T = (2n + 1) I - 2J and sum V_i V_i^T = (2n + 1) I, then

    X_1 = [[1, -e], [-e^T, S_1]],   X_i = [[1, e], [e^T, S_i]] (i > 1),   S_i = U_i (x) [[1, 1], [1, 1]] + V_i (x) [[1, -1], [-1, 1]]

are Williamson-type matrices of order 2n + 1.  Here q = 2n + 1 = 1 (mod 4) is a prime, g a primitive root,

    U_1 = -[chi(g^(2i) - g^(2j))],   U_2 = [chi(g^(2i) - g^(1 - 2j))],   U_3 = U_4 = 0,
    V_1 = I,  V_2 = 0,  V_3 = M_3,  V_4 = M_4,

with M_3, M_4 symmetric +-1 matrices, M_3 M_4 = M_4 M_3, M_3^2 + M_4^2 = 2n I.  U_1 is a symmetric circulant and U_2 a
back-circulant in the exponents, so they are amicable; U_1^2 + U_2 U_2^T = q I - 2J because Paley's core Q of order q
satisfies Q^2 = q I - J.  The pair (M_3, M_4) is the same thing as a symmetric quaternary complex Hadamard matrix
K = (M_3 + i M_4)/(1 + i) of order n.  The Williamson-type matrices give the complex Hadamard matrix
[[X, Y], [-conj(Y), conj(X)]] of order 2q, X = ((X_1 + X_2) + i (X_1 - X_2))/2, Y = ((X_3 + X_4) + i (X_3 - X_4))/2.

A doubly symmetric propus quadruple (A, W, W, B) of order m (A, B symmetric circulants, A A^T + B B^T + 2 W W^T = 4m I)
gives such a pair of order 2m (`pair_from_propus`), hence complex Hadamard matrices of order 2(4m + 1) when 4m + 1 is
a prime.

    python3 search/miyamoto_wt.py selftest      # Paley, Sylvester and small propus inputs, q = 13 ... 61
    python3 search/miyamoto_wt.py vpatterns     # the 2 x 2 V-patterns over commuting symmetric blocks (NOTES.md §19)
    python3 search/miyamoto_wt.py complete 809 complex-hadamard/inputs/H404.txt complex-hadamard/C1618.json

`complete` is the completion of NOTES.md §19: any quaternary Hadamard matrix C of order (q - 1)/2, not necessarily
symmetric, gives a quaternary Hadamard matrix of order 2q through the same U_i, V_i, with V_3 + i V_4 = (1 + i) C and
the array [[X, Y], [-Y*, X*]].
"""
import sys
import numpy as np


def prim_root(q):
    fac = [p for p in range(2, q) if (q - 1) % p == 0 and all(p % d for d in range(2, p))]
    return next(g for g in range(2, q) if all(pow(g, (q - 1) // p, q) != 1 for p in fac))


def chi_table(q):
    c = -np.ones(q, np.int64)
    c[0] = 0
    c[(np.arange(1, q) ** 2) % q] = 1
    return c


def circ(a):
    n = len(a)
    return np.array([[a[(j - i) % n] for j in range(n)] for i in range(n)], np.int64)


def u_matrices(q):
    g, chi, n = prim_root(q), chi_table(q), (q - 1) // 2
    S = [pow(g, 2 * j, q) for j in range(n)]
    Nr = [pow(g, (1 - 2 * j) % (q - 1), q) for j in range(n)]
    QSS = np.array([[chi[(x - y) % q] for y in S] for x in S], np.int64)
    QSN = np.array([[chi[(x - y) % q] for y in Nr] for x in S], np.int64)
    Z = np.zeros((n, n), np.int64)
    return [-QSS, QSN, Z, Z]


def williamson_type(q, M3, M4):
    n = (q - 1) // 2
    I = np.eye(n, dtype=np.int64)
    U, V = u_matrices(q), [I, 0 * I, M3, M4]
    J2, D2 = np.ones((2, 2), np.int64), np.array([[1, -1], [-1, 1]], np.int64)
    e = np.ones((1, 2 * n), np.int64)
    X = []
    for i in range(4):
        s = -1 if i == 0 else 1
        S = np.kron(U[i], J2) + np.kron(V[i], D2)
        X.append(np.block([[np.ones((1, 1), np.int64), s * e], [s * e.T, S]]))
    return X


def is_williamson_type(X):
    q = X[0].shape[0]
    return (all(set(np.unique(x)) <= {-1, 1} for x in X)
            and np.array_equal(sum(x @ x.T for x in X), 4 * q * np.eye(q, dtype=np.int64))
            and all(np.array_equal(X[i] @ X[j].T, X[j] @ X[i].T) for i in range(4) for j in range(i + 1, 4)))


def complex_hadamard(X):
    """exponent matrix (entry i^e) of [[X, Y], [-conj Y, conj X]], checked exactly"""
    A, B, C, D = X
    re = np.block([[(A + B) // 2, (C + D) // 2], [-(C + D) // 2, (A + B) // 2]])
    im = np.block([[(A - B) // 2, (C - D) // 2], [(C - D) // 2, -(A - B) // 2]])
    m = re.shape[0]
    ok = (np.array_equal(re @ re.T + im @ im.T, m * np.eye(m, dtype=np.int64))
          and not (im @ re.T - re @ im.T).any() and np.array_equal(re * re + im * im, np.ones((m, m), np.int64)))
    E = np.select([re == 1, im == 1, re == -1], [0, 1, 2], 3)
    return E, ok


def complex_general(X):
    """exponent matrix of K = [[X, Y], [-Y*, X*]] (* = conjugate transpose), X = ((X_1 + X_2) + i (X_1 - X_2))/2,
    Y = ((X_3 + X_4) + i (X_3 - X_4))/2, checked exactly (entries i^e and K K* = m I); for symmetric X_i this is the
    array of complex_hadamard"""
    A, B, C, D = X
    xr, xi, yr, yi = (A + B) // 2, (A - B) // 2, (C + D) // 2, (C - D) // 2
    re = np.block([[xr, yr], [-yr.T, xr.T]])
    im = np.block([[xi, yi], [yi.T, -xi.T]])
    m = re.shape[0]
    ok = (np.array_equal(re * re + im * im, np.ones((m, m), np.int64))
          and np.array_equal(re @ re.T + im @ im.T, m * np.eye(m, dtype=np.int64)) and not (im @ re.T - re @ im.T).any())
    return np.select([re == 1, im == 1, re == -1], [0, 1, 2], 3), ok


def completion(q, C):
    """the completion of the printed proof of Seberry-Yamada's Corollary 8.8(i) (NOTES.md section 19): for a prime
    q = 1 (mod 4) and any quaternary Hadamard matrix C of order (q - 1)/2 (exponent matrix, entry i^e), write
    (1 + i) C = M_3 + i M_4 and take U = (A, BR, 0, 0), V = (I, 0, M_3, M_4) in Lemma 8.6.  The X_i are not
    Williamson-type unless M_3, M_4 are symmetric, but K = [[X, Y], [-Y*, X*]] is a quaternary Hadamard matrix of
    order 2q.  Returns the exponent matrix of K and the exact check."""
    M3 = np.choose(C, [1, -1, -1, 1]).astype(np.int64)        # (1 + i) i^e = M3 + i M4
    M4 = np.choose(C, [1, 1, -1, -1]).astype(np.int64)
    return complex_general(williamson_type(q, M3, M4))


def pair_from_paley(p):
    """M3, M4 of order p + 1 from the symmetric Legendre pair (chi + d_0, -chi + d_0), p = 1 (mod 4) prime"""
    chi = chi_table(p)
    out = []
    for s in (np.where(np.arange(p) == 0, 1, chi), np.where(np.arange(p) == 0, 1, -chi)):
        C, r = circ(s), int(s.sum())
        e = np.ones((1, p), np.int64)
        out.append(np.block([[np.array([[-r]]), r * e], [r * e.T, C]]))
    return out


def sylvester(k):
    H = np.ones((1, 1), np.int64)
    for _ in range(k):
        H = np.block([[H, H], [H, -H]])
    return H


def pair_from_propus(A, W, B):
    """M3, M4 of order 2m from a doubly symmetric propus quadruple (A, W, W, B) of circulants of order m (A, B symmetric,
    A A^T + B B^T + 2 W W^T = 4m I): M3 = [[B, W], [W^T, -A]], M4 = [[A, W], [W^T, -B]], i.e. (1 + i) K with
    K = [[X, W], [W^T, -conj X]], X = ((A + B) + i (A - B))/2"""
    A, W, B = circ(A), circ(W), circ(B)
    return [np.block([[B, W], [W.T, -A]]), np.block([[A, W], [W.T, -B]])]


def find_dsym_propus(m):
    """first doubly symmetric propus quadruple of order m (brute force, small m)"""
    import itertools
    h = (m - 1) // 2
    syms = []
    for bits in itertools.product((1, -1), repeat=h + 1):
        a = [bits[0]] + [0] * (m - 1)
        for i in range(1, h + 1):
            a[i] = a[m - i] = bits[i]
        syms.append(np.array(a))
    paf = lambda x: np.array([int(x @ np.roll(x, -s)) for s in range(m)])
    P = [paf(a) for a in syms]
    need = {}
    for bits in itertools.product((1, -1), repeat=m):
        w = np.array(bits)
        need.setdefault(tuple(-2 * paf(w)[1:]), w)
    for i, a in enumerate(syms):
        for j, b in enumerate(syms):
            k = tuple((P[i] + P[j])[1:])
            if k in need:
                return a, need[k], b
    return None


def selftest():
    cases = [(13, pair_from_paley(5)), (29, pair_from_paley(13)), (37, pair_from_paley(17)), (61, pair_from_paley(29)),
             (17, [sylvester(3)] * 2)]
    for m in (3, 7, 9):
        quad = find_dsym_propus(m)
        if quad is not None:
            cases.append((4 * m + 1, pair_from_propus(*quad)))
    ok_all = True
    for q, (M3, M4) in cases:
        n = (q - 1) // 2
        good = (np.array_equal(M3, M3.T) and np.array_equal(M4, M4.T) and np.array_equal(M3 @ M4, M4 @ M3)
                and np.array_equal(M3 @ M3 + M4 @ M4, 2 * n * np.eye(n, dtype=np.int64)))
        X = williamson_type(q, M3, M4)
        wt = is_williamson_type(X)
        E, ck = complex_hadamard(X)
        print(f"q = {q}: pair of order {n} {'ok' if good else 'BAD'}; Williamson-type of order {q} "
              f"{'ok' if wt else 'FAILED'}; complex Hadamard of order {E.shape[0]} {'ok' if ck else 'FAILED'}")
        ok_all &= good and wt and ck
    return ok_all


def vpatterns():
    """the 2 x 2 block patterns for the V's of Lemma 8.6 over four commuting symmetric blocks a1..a4: V_i =
    [[p_i, q_i], [r_i, s_i]] with entries 0, +-1, +-a_k, each row of [V_1 .. V_4] holding one +-1, a1..a4 once each and
    three zeros, sum V_i V_i^T = (1 + sum a_k^2) I and V_i V_j^T = V_j V_i^T as polynomial identities.  First rows are
    taken up to relabelling and signs of the a_k.  Prints the solutions grouped by the shape of each V_i (I = the +-1,
    a = a block, 0 = zero); the U_i are zero where V_i has an a, and a zero-diagonal class block where V_i has the 1."""
    import itertools
    from collections import Counter

    def prod(e, f):
        if e is None or f is None:
            return None
        return e[0] * f[0], tuple(sorted(v for v in (e[1], f[1]) if v))

    def dot(x, y):
        c = Counter()
        for e, f in zip(x, y):
            t = prod(e, f)
            if t:
                c[t[1]] += t[0]
        return {k: v for k, v in c.items() if v}

    rows = []
    for one in range(8):
        rest = [i for i in range(8) if i != one]
        for apos in itertools.permutations(rest, 4):
            for sg in itertools.product((1, -1), repeat=5):
                r = [None] * 8
                r[one] = (sg[0], 0)
                for k, pos in enumerate(apos):
                    r[pos] = (sg[k + 1], k + 1)
                rows.append(tuple(r))
    canon = [r for r in rows if all(e is None or e[0] == 1 for e in r) and [e[1] for e in r if e and e[1]] == [1, 2, 3, 4]]
    shapes = Counter()
    for X in canon:
        x = [(X[2 * i], X[2 * i + 1]) for i in range(4)]
        for Y in rows:
            y = [(Y[2 * i], Y[2 * i + 1]) for i in range(4)]
            if any(dot(x[i], y[j]) != dot(x[j], y[i]) for i in range(4) for j in range(i + 1, 4)):
                continue
            tr = Counter()
            for i in range(4):
                for k, v in dot(x[i], y[i]).items():
                    tr[k] += v
            if any(tr.values()):
                continue
            pat = tuple(sorted("".join("0" if e is None else ("I" if e[1] == 0 else "a")
                                       for e in (X[2 * i], X[2 * i + 1], Y[2 * i], Y[2 * i + 1])) for i in range(4)))
            shapes[pat] += 1
    print(sum(shapes.values()), "solutions")
    for k, v in shapes.most_common():
        print(v, " ".join(k))


def read_exponents(path):
    """a matrix file: rows over +, - (real) or over the digits 0..3 (entry i^e)"""
    rows = [l.strip() for l in open(path) if l.strip() and not l.startswith("#")]
    return np.array([[{"+": 0, "-": 2}.get(ch, None) if ch in "+-" else int(ch) for ch in r] for r in rows], np.int64)


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "complete":
        import hashlib, json
        from pathlib import Path
        q, inp, out = int(sys.argv[2]), sys.argv[3], sys.argv[4]
        C = read_exponents(inp)
        E, ok = completion(q, C)
        assert ok, "not a complex Hadamard matrix"
        text = "".join("".join(map(str, r)) + "\n" for r in E).encode()
        root = Path(__file__).resolve().parent.parent / "complex-hadamard"
        rel = Path(inp).resolve().relative_to(root)
        cert = {"order": 2 * q, "N": q, "construction": "miyamoto-completion", "q": q, "g": prim_root(q),
                "input": str(rel), "input_sha256": hashlib.sha256(open(inp, "rb").read()).hexdigest(),
                "sha256": hashlib.sha256(text).hexdigest()}
        Path(out).write_text(json.dumps(cert, indent=1) + "\n")
        print(f"q = {q}: quaternary Hadamard matrix of order {2 * q} from {rel}, sha256 {cert['sha256']}")
        raise SystemExit(0)
    if len(sys.argv) > 1 and sys.argv[1] == "cw":
        import hashlib, json
        from pathlib import Path
        root = Path(__file__).resolve().parent.parent
        sys.path.insert(0, str(root / "verify"))
        import verify as V
        q, inp, tfile = int(sys.argv[2]), Path(sys.argv[3]).resolve(), Path(sys.argv[4])
        M, _, digest = V.build_checked(json.loads(inp.read_text()))
        assert M.shape[0] == (q - 1) // 2 and np.array_equal(M, M.T) and q % 8 == 1, "need a symmetric H((q - 1)/2)"
        M = M.astype(np.int64)
        X = williamson_type(q, M, M)
        F = [x.astype(np.float64) for x in X]
        assert np.array_equal(sum(x @ x.T for x in F), 4.0 * q * np.eye(q)) and \
            all(np.array_equal(F[i] @ F[j].T, F[j] @ F[i].T) for i in range(4) for j in range(i + 1, 4)), "not Williamson type"
        rows = [l.strip() for l in tfile.read_text().splitlines() if l.strip()]
        T = [np.array([V.SIGN[c] for c in r], np.int64) for r in rows]
        t = len(rows[0])
        H = V.assemble_cooper_wallis(T, X)
        sha = hashlib.sha256()
        for chunk in V.csv_chunks(H):
            sha.update(chunk)
        cert = {"order": 4 * t * q, "construction": "cooper-wallis-miyamoto", "t": t, "q": q, "g": prim_root(q),
                "input": inp.relative_to(root).as_posix(), "input_sha256": digest, "t_matrices": rows,
                "sha256": sha.hexdigest()}
        d = root / "hadamard" / str(4 * t * q)
        d.mkdir(exist_ok=True)
        (d / "certificate.json").write_text(json.dumps(cert, indent=1) + "\n")
        print(f"wrote {d}/certificate.json (order {4 * t * q}, sha256 {cert['sha256']})")
        raise SystemExit(0)
    if len(sys.argv) > 1 and sys.argv[1] == "selftest":
        raise SystemExit(0 if selftest() else 1)
    if len(sys.argv) > 1 and sys.argv[1] == "vpatterns":
        vpatterns()
        raise SystemExit(0)
    print(__doc__)
