#!/usr/bin/env python3
"""Check the complex (quaternary) Hadamard matrices in complex-hadamard/.

A quaternary Hadamard matrix of order n has entries 1, -1, i, -i and K K* = n I. Each certificate
complex-hadamard/C<n>.json gives four +-1 sequences A, B, C, D of length N = n/2 (first rows of circulant matrices)
with A and B symmetric and A A^T + B B^T + C C^T + D D^T = 4N I, and the matrix is Szöllősi's array
(arXiv:2603.09572, Theorem 4), R the back-diagonal identity:

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

This script rebuilds K from the sequences alone, checks every entry and K K* = n I exactly in integer arithmetic
(real and imaginary parts separately), and compares the SHA-256 of the canonical text of K: n lines of n digits,
digit e for the entry i^e.

    python3 verify/verify_complex.py [n ...]          all certificates when no order is given
    python3 verify/verify_complex.py --write n         also write complex-hadamard/C<n>.txt
"""
import hashlib, json, sys
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parent.parent / "complex-hadamard"


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


def build(cert):
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
    re = Kr @ Kr.T + Ki @ Ki.T                 # K K* = (Kr + i Ki)(Kr^T - i Ki^T)
    im = Ki @ Kr.T - Kr @ Ki.T
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
        cert = json.loads(f.read_text())
        n = cert["order"]
        try:
            Kr, Ki = build(cert)
            check(Kr, Ki)
            txt = text(Kr, Ki)
            sha = hashlib.sha256(txt.encode()).hexdigest()
            ok = sha == cert.get("sha256")
            print(f"{n:6d}  {'ok' if ok else 'HASH MISMATCH'}  K K* = {n} I   sha256 {sha}")
            bad += not ok
            if write:
                (ROOT / f"C{n}.txt").write_text(txt)
        except ValueError as ex:
            print(f"{n:6d}  FAIL  {ex}"); bad += 1
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
