#!/usr/bin/env python3
"""Check the complex (quaternary) Hadamard matrices in complex-hadamard/.

A quaternary Hadamard matrix of order n has entries 1, -1, i, -i and K K* = n I. Each certificate
complex-hadamard/C<n>.json gives four +-1 sequences A, B, C, D of length N = n/2 (first rows of circulant matrices)
with A and B symmetric and A A^T + B B^T + C C^T + D D^T = 4N I, and the matrix is Szöllősi's array
(arXiv:2603.09572, Theorem 4), R the back-diagonal identity:

    K = [[ (A+B)/2 + i(A-B)/2,          (C+DR)/2 + i(C-DR)/2 ],
         [ (C^T+DR)/2 - i(C^T-DR)/2,   -(A+B)/2 + i(A-B)/2  ]].

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


def build(cert):
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
