#!/usr/bin/env python3
"""Check the complex (quaternary) Hadamard matrices in complex-hadamard/.

A quaternary Hadamard matrix of order n has entries 1, -1, i, -i and K K* = n I. Each certificate is the file
complex-hadamard/C<n>.json. Certificates with construction "szollosi-ssxx" give four +-1 sequences A, B, C, D of
length N = n/2 (first rows of circulant matrices) with A and B symmetric and A A^T + B B^T + C C^T + D D^T = 4N I,
and the matrix is Szöllősi's array (arXiv:2603.09572, Theorem 4), R the back-diagonal identity:

    K = [[ (A+B)/2 + i(A-B)/2,          (C+DR)/2 + i(C-DR)/2 ],
         [ (C^T+DR)/2 - i(C^T-DR)/2,   -(A+B)/2 + i(A-B)/2  ]].

Certificates with construction "szollosi-ssxx-group" use group-developed matrices over G = Z_{m_1} x ... x Z_{m_r}
(field "group"; elements in mixed-radix order, the last factor fastest) and R the matrix of g -> -g.
Certificates with construction "szollosi-bordered" (or "-group") give four sequences on a group of order N - 1
with sum of PAFs -4 off 0 and the border (x, s, t, y, u, v) of the bordered array (search/bsz.py):
K = [[X, Y], [Y^*, -X^*]], X = [[x, s e^T], [t e, (A+B)/2 + i(A-B)/2]], Y = [[y, u e^T], [v e, (C+DR)/2 + i(C-DR)/2]].

Certificates with construction "szollosi-tbordered" give four sequences on a group G of order N - t, the fibre
labels of a map phi: G -> Q and a border (Xb, Rh, Ch) of exponent matrices (complex-hadamard/TBORDERED.md):
K = [[Xb, R], [R', Kc]] with Kc the array above over G, R[a, (i, g)] = Rh[a, i q + phi(g)], R'[(i, g), b] =
Ch[i q + phi(g), b].

Certificates with construction "miyamoto-completion" give q, a primitive root g and an input file: the completion
of Seberry-Yamada's Corollary 8.8(i) (completion_build, NOTES.md section 19).

Certificates with construction "paley" give a prime power q = p^k = 3 (mod 4) (fields "q", "p", "k", and "poly" for
k > 1) and the order (q + 1)/2: the quaternary matrix obtained from Paley's skew Hadamard matrix of order q + 1 and
the signed permutation u -> -1/u that commutes with it (paley_build, complex-hadamard/PALEY.md).

Certificates with construction "complex-yamada-kiyasu" give a prime q = 1 (mod 4), Yamada's parameters r and xi and,
under "core", a quaternary conference matrix of order (q + 3)/2 as the four sequences of a zero-diagonal Szöllősi
array over the group "group", plain or bordered (optional fields "border" {t, v}, "reflection", "rowsums",
"symmetry_generators", "symmetry_order"). K of order 2(q + 2) is built by verify/verify.py (complex_yamada_kiyasu,
NOTES.md §31, §34), which also builds the Hadamard matrices of order 4(q + 2) in hadamard/ from the same data; for
these certificates only the final check below is independent of verify.py. All other constructions are built here.

The script checks that the built matrix has the order stated in the certificate and in the file name (and that the
field "N" is n/2), that every entry is 1, -1, i or -i, and that K K* = n I, and compares the SHA-256 of the canonical
text of K: n lines of n digits, digit e for the entry i^e. K K* is formed from the real and imaginary parts as float64
matrix products, which are exact: every entry is 0 or +-1 and every partial sum is an integer of absolute value at
most n, far below 2^53.

    python3 verify/verify_complex.py [n ...]          all certificates when no order is given
    python3 verify/verify_complex.py --write n         also write complex-hadamard/C<n>.txt
"""
import hashlib, importlib.util, json, sys
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parent.parent / "complex-hadamard"
_VERIFY = None


def _real():
    """verify/verify.py, loaded once by path: it holds the builder of the complex-yamada-kiyasu certificates."""
    global _VERIFY
    if _VERIFY is None:
        path = Path(__file__).resolve().with_name("verify.py")
        spec = importlib.util.spec_from_file_location("hadamard_verify", path)
        _VERIFY = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(_VERIFY)
    return _VERIFY


def circ(s):
    x = np.array([1 if c == "+" else -1 for c in s], dtype=np.int64)
    n = len(x); t = np.arange(n)
    return x[(t[None, :] - t[:, None]) % n]


def group_developed(s, mods):
    """X[g, h] = x(h - g) over G = Z_{m_1} x ... x Z_{m_r}, elements in mixed-radix order (last factor fastest)."""
    x = np.array([1 if c == "+" else -1 for c in s], dtype=np.int64)
    mods = [int(m) for m in mods]
    v = int(np.prod(mods))
    coords = np.array(np.unravel_index(np.arange(v), tuple(mods))).T
    d = (coords[None, :, :] - coords[:, None, :]) % np.array(mods)
    return x[np.ravel_multi_index(tuple(np.moveaxis(d, 2, 0)), tuple(mods))], coords, mods


def bordered_build(cert):
    """Bordered Szöllősi array: sequences on G of order N - 1 (cyclic Z_(N-1) or the mixed-radix group "group"),
    A and B symmetric, sum of PAFs -4 off 0; X = [[x, s e^T], [t e, (A+B)/2 + i(A-B)/2]], Y = [[y, u e^T],
    [v e, (C+DR)/2 + i(C-DR)/2]], K = [[X, Y], [Y^*, -X^*]], border (x, s, t, y, u, v) from the certificate."""
    seqs = cert["sequences"]
    m = len(seqs[0]); N = m + 1
    mods = cert.get("group", [m])
    mats = [group_developed(s, mods) for s in seqs]
    A, B, C, D = (M[0] for M in mats)
    coords, mods = mats[0][1], mats[0][2]
    if 2 * N != cert["order"]:
        raise ValueError("lengths")
    if not (np.array_equal(A, A.T) and np.array_equal(B, B.T)):
        raise ValueError("A and B must be symmetric")
    I = np.eye(m, dtype=np.int64)
    if not np.array_equal(A @ A.T + B @ B.T + C @ C.T + D @ D.T, 4 * m * I - 4 * (np.ones((m, m), np.int64) - I)):
        raise ValueError("A A^T + B B^T + C C^T + D D^T != 4(N-1) I - 4(J - I)")
    neg = np.ravel_multi_index(tuple(((-coords) % np.array(mods)).T), tuple(mods))
    R = np.zeros((m, m), np.int64); R[np.arange(m), neg] = 1
    DR = D @ R
    (xr, xi), (sr, si), (tr, ti), (yr, yi), (ur, ui), (vr, vi) = cert["border"]

    def bordered(c, top, left, Mr, Mi):
        Kr = np.empty((N, N), np.int64); Ki = np.empty((N, N), np.int64)
        Kr[0, 0], Ki[0, 0] = c; Kr[0, 1:], Ki[0, 1:] = top; Kr[1:, 0], Ki[1:, 0] = left
        Kr[1:, 1:], Ki[1:, 1:] = Mr, Mi
        return Kr, Ki

    XR, XI = bordered((xr, xi), (sr, si), (tr, ti), (A + B) // 2, (A - B) // 2)
    YR, YI = bordered((yr, yi), (ur, ui), (vr, vi), (C + DR) // 2, (C - DR) // 2)
    return np.block([[XR, YR], [YR.T, -XR.T]]), np.block([[XI, YI], [-YI.T, XI.T]])


def tbordered_build(cert):
    """Character-bordered Szöllősi array (complex-hadamard/TBORDERED.md): sequences on G = "group" (order m), A and B
    symmetric; phi: G -> Q given as the list of fibre labels ("phi", q = max + 1); border Xb (2t x 2t), Rh (2t x 2q),
    Ch (2q x 2t) as exponent matrices (entry i^e). K = [[Xb, R], [R', Kc]], Kc = [[X, Y], [Y^*, -X^*]] as in
    szollosi-ssxx-group, R[a, (i, g)] = Rh[a, i q + phi(g)], R'[(i, g), b] = Ch[i q + phi(g), b]."""
    mods = cert["group"]
    mats = [group_developed(s, mods) for s in cert["sequences"]]
    A, B, C, D = (M[0] for M in mats)
    coords, mods = mats[0][1], mats[0][2]
    m = A.shape[0]
    phi = np.array(cert["phi"], dtype=np.int64)
    q = int(phi.max()) + 1
    Xb, Rh, Ch = (np.array(cert["border"][k], dtype=np.int64) for k in ("Xb", "Rh", "Ch"))
    t = Xb.shape[0] // 2
    if 2 * (m + t) != cert["order"] or len(phi) != m or Rh.shape != (2 * t, 2 * q) or Ch.shape != (2 * q, 2 * t):
        raise ValueError("sizes")
    if not (np.array_equal(A, A.T) and np.array_equal(B, B.T)):
        raise ValueError("A and B must be symmetric")
    neg = np.ravel_multi_index(tuple(((-coords) % np.array(mods)).T), tuple(mods))
    R = np.zeros((m, m), np.int64); R[np.arange(m), neg] = 1
    DR = D @ R
    Kcr = np.block([[(A + B) // 2, (C + DR) // 2], [(C.T + DR) // 2, -(A + B) // 2]])
    Kci = np.block([[(A - B) // 2, (C - DR) // 2], [-(C.T - DR) // 2, (A - B) // 2]])
    re, im = np.array([1, 0, -1, 0]), np.array([0, 1, 0, -1])
    cols = np.concatenate([phi, q + phi])
    Kr = np.block([[re[Xb], re[Rh][:, cols]], [re[Ch][cols, :], Kcr]])
    Ki = np.block([[im[Xb], im[Rh][:, cols]], [im[Ch][cols, :], Kci]])
    return Kr, Ki


def paley_build(cert):
    """Paley's skew Hadamard matrix H = I + S of order q + 1 (q = p^k = 3 mod 4) and the signed permutation P of
    u -> -1/u, which commutes with H and has P^2 = -I (complex-hadamard/PALEY.md). Points: x = 0..q-1 (base-p digits
    of the polynomial; "poly" = c_0, ..., c_(k-1), 1 is the modulus for k > 1) and infinity = q; S[u, v] =
    chi(det(rep u, rep v)) with reps (x, 1), (1, 0). K[k, l] = (H[f_k, f_l] + i s_l H[f_k, g(f_l)])/(1 + i), f the
    smaller point of each pair {u, g(u)}, P e_u = s_u e_g(u). Field arithmetic here is by polynomial
    multiplication (squares and inverses found directly), independent of search/paley_complex.py."""
    q, p, k = int(cert["q"]), int(cert["p"]), int(cert["k"])
    if p ** k != q or q % 4 != 3 or 2 * cert["order"] != q + 1:
        raise ValueError("parameters")
    f = [int(c) for c in cert.get("poly", [0, 1])]
    if len(f) != k + 1 or f[-1] != 1:
        raise ValueError("modulus")
    digits = [[(x // p ** i) % p for i in range(k)] for x in range(q)]
    enc = lambda v: sum((c % p) * p ** i for i, c in enumerate(v))

    def mul(x, y):
        a, b = digits[x], digits[y]
        pr = [0] * (2 * k - 1)
        for i in range(k):
            for j in range(k):
                pr[i + j] += a[i] * b[j]
        for d in range(2 * k - 2, k - 1, -1):          # t^k = -(c_0 + ... + c_(k-1) t^(k-1))
            c, pr[d] = pr[d], 0
            for i in range(k):
                pr[d - k + i] -= c * f[i]
        return enc(pr[:k])

    sq = {mul(y, y) for y in range(1, q)}
    if len(sq) != (q - 1) // 2:
        raise ValueError("the modulus is not irreducible")
    chi = np.array([0] + [1 if x in sq else -1 for x in range(1, q)], np.int64)
    neg = np.array([enc([-c for c in digits[x]]) for x in range(q)], np.int64)
    D = np.array(digits, np.int64)
    w = p ** np.arange(k)
    sub = ((D[:, None, :] - D[None, :, :]) % p) @ w          # sub[x, y] = x - y
    n = q + 1
    H = np.empty((n, n), np.int64)
    H[:q, :q] = chi[sub]
    H[q, :q] = chi[1]                                        # det((1, 0), (y, 1)) = 1
    H[:q, q] = chi[neg[1]]                                   # det((x, 1), (1, 0)) = -1
    H[np.arange(n), np.arange(n)] = 1
    # g(x, y) = (y, -x): (x, 1) -> (1, -x) = (-x)(-1/x, 1) for x != 0; (0, 1) -> (1, 0); (1, 0) -> (0, -1) = -(0, 1)
    inv = {}
    for x in range(1, q):
        e, r, b = q - 2, 1, x                                # x^(q - 2)
        while e:
            if e & 1:
                r = mul(r, b)
            b, e = mul(b, b), e >> 1
        if mul(r, x) != 1:
            raise ValueError("inverse")
        inv[x] = r
    img = np.array([q] + [int(neg[inv[x]]) for x in range(1, q)] + [0], np.int64)
    s = np.array([chi[1]] + [int(chi[neg[x]]) for x in range(1, q)] + [int(chi[neg[1]])], np.int64)
    first = np.array([u for u in range(n) if u < img[u]], np.int64)
    if len(first) != n // 2:
        raise ValueError("pairs")
    A = H[np.ix_(first, first)]
    B = s[first][None, :] * H[np.ix_(first, img[first])]
    # (A + iB)/(1 + i) = ((A + B) + i (B - A))/2
    return (A + B) // 2, (B - A) // 2


def completion_build(cert):
    """The completion of the printed proof of Seberry-Yamada's Corollary 8.8(i) (NOTES.md section 19). q = 1 (mod 4)
    a prime, g a primitive root, n = (q - 1)/2, chi the quadratic character, s_i = g^(2i) and t_j = g^(1 - 2j).
    U_1 = -[chi(s_i - s_j)], U_2 = [chi(s_i - t_j)], U_3 = U_4 = 0; V_1 = I, V_2 = 0, V_3 + i V_4 = (1 + i) C for the
    input quaternary Hadamard matrix C of order n (complex-hadamard/<input>, rows over +, - or the digits e of i^e).
    S_k = U_k (x) [[1, 1], [1, 1]] + V_k (x) [[1, -1], [-1, 1]], X_1 = [[1, -e], [-e^T, S_1]], X_k = [[1, e], [e^T,
    S_k]] (k > 1), X = ((X_1 + X_2) + i (X_1 - X_2))/2, Y = ((X_3 + X_4) + i (X_3 - X_4))/2 and
    K = [[X, Y], [-Y*, X*]], * the conjugate transpose. Written from the statement, independently of
    search/miyamoto_wt.py; the input is checked exactly (C C* = n I) and its SHA-256 compared."""
    q, g = int(cert["q"]), int(cert["g"])
    if q % 4 != 1 or any(q % d == 0 for d in range(2, int(q ** 0.5) + 1)) or cert["order"] != 2 * q:
        raise ValueError("q must be a prime = 1 (mod 4) and the order 2q")
    fac = [d for d in range(2, q) if (q - 1) % d == 0 and all(d % e for e in range(2, int(d ** 0.5) + 1))]
    if any(pow(g, (q - 1) // d, q) == 1 for d in fac):
        raise ValueError("g is not a primitive root")
    raw = (ROOT / cert["input"]).read_bytes()
    if hashlib.sha256(raw).hexdigest() != cert["input_sha256"]:
        raise ValueError("input file changed")
    rows = [l.strip() for l in raw.decode().splitlines() if l.strip()]
    E = np.array([[{"+": 0, "-": 2, "0": 0, "1": 1, "2": 2, "3": 3}[ch] for ch in r] for r in rows], np.int64)
    n = (q - 1) // 2
    if E.shape != (n, n):
        raise ValueError("input must have order (q - 1)/2")
    Cr = np.choose(E, [1, 0, -1, 0]); Ci = np.choose(E, [0, 1, 0, -1])
    if not (np.array_equal(Cr @ Cr.T + Ci @ Ci.T, n * np.eye(n, dtype=np.int64)) and not (Ci @ Cr.T - Cr @ Ci.T).any()):
        raise ValueError("the input is not a quaternary Hadamard matrix")
    V3, V4 = Cr - Ci, Cr + Ci                               # (1 + i)(Cr + i Ci)
    chi = -np.ones(q, np.int64); chi[0] = 0; chi[(np.arange(1, q) ** 2) % q] = 1
    s_ = np.array([pow(g, 2 * i, q) for i in range(n)]); t_ = np.array([pow(g, (1 - 2 * j) % (q - 1), q) for j in range(n)])
    U1 = -chi[(s_[:, None] - s_[None, :]) % q]; U2 = chi[(s_[:, None] - t_[None, :]) % q]
    Z, I = np.zeros((n, n), np.int64), np.eye(n, dtype=np.int64)
    J2, D2 = np.ones((2, 2), np.int64), np.array([[1, -1], [-1, 1]], np.int64)
    e = np.ones((1, 2 * n), np.int64)
    X = []
    for k, (U, V) in enumerate(((U1, I), (U2, Z), (Z, V3), (Z, V4))):
        S = np.kron(U, J2) + np.kron(V, D2)
        sg = -1 if k == 0 else 1
        X.append(np.block([[np.ones((1, 1), np.int64), sg * e], [sg * e.T, S]]))
    xr, xi = (X[0] + X[1]) // 2, (X[0] - X[1]) // 2
    yr, yi = (X[2] + X[3]) // 2, (X[2] - X[3]) // 2
    return np.block([[xr, yr], [-yr.T, xr.T]]), np.block([[xi, yi], [yi.T, -xi.T]])


def build(cert):
    if cert["construction"] == "complex-yamada-kiyasu":
        Kr, Ki = _real().complex_yamada_kiyasu(cert)
        return Kr.astype(np.int64), Ki.astype(np.int64)
    if cert["construction"] == "miyamoto-completion":
        return completion_build(cert)
    if cert["construction"] == "paley":
        return paley_build(cert)
    if cert["construction"] == "szollosi-tbordered":
        return tbordered_build(cert)
    if cert["construction"] in ("szollosi-bordered", "szollosi-bordered-group"):
        return bordered_build(cert)
    if cert["construction"] == "szollosi-ssxx-group":
        mats = [group_developed(s, cert["group"]) for s in cert["sequences"]]
        A, B, C, D = (m[0] for m in mats)
        coords, mods = mats[0][1], mats[0][2]
        N = A.shape[0]
        if 2 * N != cert["order"]:
            raise ValueError("lengths")
        if not (np.array_equal(A, A.T) and np.array_equal(B, B.T)):
            raise ValueError("A and B must be symmetric")
        I = np.eye(N, dtype=np.int64)
        if not np.array_equal(A @ A.T + B @ B.T + C @ C.T + D @ D.T, 4 * N * I):
            raise ValueError("A A^T + B B^T + C C^T + D D^T != 4N I")
        neg = np.ravel_multi_index(tuple(((-coords) % np.array(mods)).T), tuple(mods))
        R = np.zeros((N, N), np.int64); R[np.arange(N), neg] = 1
        DR = D @ R
        Kr = np.block([[(A + B) // 2, (C + DR) // 2], [(C.T + DR) // 2, -(A + B) // 2]])
        Ki = np.block([[(A - B) // 2, (C - DR) // 2], [-(C.T - DR) // 2, (A - B) // 2]])
        return Kr, Ki
    if cert["construction"] != "szollosi-ssxx":
        raise ValueError("unknown construction")
    A, B, C, D = (circ(s) for s in cert["sequences"])
    N = A.shape[0]
    if any(X.shape != (N, N) for X in (B, C, D)) or 2 * N != cert["order"]:
        raise ValueError("lengths")
    if not (np.array_equal(A, A.T) and np.array_equal(B, B.T)):
        raise ValueError("A and B must be symmetric")
    I = np.eye(N, dtype=np.int64)
    if not np.array_equal(A @ A.T + B @ B.T + C @ C.T + D @ D.T, 4 * N * I):
        raise ValueError("A A^T + B B^T + C C^T + D D^T != 4N I")
    t = np.arange(N)
    R = np.zeros((N, N), np.int64); R[t, N - 1 - t] = 1
    DR = D @ R
    Kr = np.block([[(A + B) // 2, (C + DR) // 2], [(C.T + DR) // 2, -(A + B) // 2]])
    Ki = np.block([[(A - B) // 2, (C - DR) // 2], [-(C.T - DR) // 2, (A - B) // 2]])
    return Kr, Ki


def check(Kr, Ki):
    n = Kr.shape[0]
    if not np.all(np.abs(Kr) + np.abs(Ki) == 1):
        raise ValueError("an entry is not in {1, -1, i, -i}")
    I = np.eye(n, dtype=np.int64)
    # K K* = (Kr + i Ki)(Kr^T - i Ki^T); float64 products are exact: entries 0, +-1 and every partial sum an integer
    # of absolute value <= n < 2^53
    Fr, Fi = Kr.astype(np.float64), Ki.astype(np.float64)
    re, im = (Fr @ Fr.T + Fi @ Fi.T).astype(np.int64), (Fi @ Fr.T - Fr @ Fi.T).astype(np.int64)
    if not (np.array_equal(re, n * I) and not im.any()):
        raise ValueError("K K* != n I")


def text(Kr, Ki):
    e = np.where(Kr == 1, 0, np.where(Ki == 1, 1, np.where(Kr == -1, 2, 3)))
    return "".join("".join(map(str, row)) + "\n" for row in e)


def main():
    args = sys.argv[1:]
    write = "--write" in args
    orders = [int(a) for a in args if a.isdigit()]
    files = [ROOT / f"C{n}.json" for n in orders] if orders else sorted(ROOT.glob("C*.json"),
                                                                        key=lambda p: int(p.stem[1:]))
    bad = 0
    for f in files:
        stem = int(f.stem[1:])
        try:
            cert = json.loads(f.read_text())
            n = cert["order"]
            if "N" in cert and 2 * cert["N"] != n:
                raise ValueError(f'"N" is {cert["N"]}, the order is {n}')
            Kr, Ki = build(cert)
            if Kr.shape[0] != n or n != stem:
                raise ValueError(f"built order {Kr.shape[0]}, certificate {n}, file {f.name}")
            check(Kr, Ki)
            txt = text(Kr, Ki)
            sha = hashlib.sha256(txt.encode()).hexdigest()
            ok = sha == cert.get("sha256")
            print(f"{n:6d}  {'ok' if ok else 'HASH MISMATCH'}  K K* = {n} I   sha256 {sha}")
            bad += not ok
            if write:
                (ROOT / f"C{n}.txt").write_text(txt)
        except Exception as ex:      # report and continue with the remaining orders
            why = (f"missing certificate field or unrecognised key {ex}" if isinstance(ex, KeyError)
                   else str(ex) if isinstance(ex, ValueError) else f"{type(ex).__name__}: {ex}")
            print(f"{stem:6d}  FAIL  {why}"); bad += 1
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
