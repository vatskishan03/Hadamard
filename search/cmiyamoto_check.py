#!/usr/bin/env python3
"""Exact checks of Miyamoto's Theorem 7 in quaternary form (hadamard/COMPLEX_MIYAMOTO.md) for the small q.

For each q the quaternary Hadamard matrix K of order q - 1 is built, verify.complex_miyamoto(q, Kr, Ki) gives
H = Hr + i Hi of order 4q, and the script checks that every entry is a unit and that H H* = 4qI exactly (float64
products of integers far below 2^53).

- q = 7, 19, 31, 43: K = C + iI with C Paley's symmetric conference matrix of order q - 1 (q - 2 a prime = 1 mod 4).
- q = 223, 263, 331, 431: K from complex-hadamard/C222.json, C262.json, C330.json and C430.json.
- with --big also q = 1559 and 1619 (complex-hadamard/C1558.json and C1618.json, the matrices of orders 12472 and
  12952), checked through Turyn's double in float32 (exact below 2^24); about a minute each.
- with --general the theorem with a skew conference matrix C of order q + 1 other than Paley's, for q = 15, 39, 63
  (not prime powers): C is Paley's skew conference matrix of order 8, 20 or 8 doubled once, once or three times
  (skew Hadamard H -> [[H, H], [-H^T, H^T]]), under a random simultaneous signed permutation; K as for Paley's q.

    python3 search/cmiyamoto_check.py [--big] [--general] [--q Q ...]
"""
import argparse, json, sys, time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "verify"))
import verify, verify_complex

PALEY_Q = [7, 19, 31, 43]
FILE_Q = {223: "C222.json", 263: "C262.json", 331: "C330.json", 431: "C430.json",
          1559: "C1558.json", 1619: "C1618.json"}


def is_prime(n):
    return n > 1 and all(n % p for p in range(2, int(n ** 0.5) + 1))


def paley_K(q):
    """C + iI of order q - 1, C Paley's symmetric conference matrix of order p + 1, p = q - 2 a prime = 1 (mod 4)"""
    p = q - 2
    if not (is_prime(p) and p % 4 == 1):
        raise ValueError("q - 2 must be a prime = 1 (mod 4)")
    sq = {x * x % p for x in range(1, p)}
    chi = np.array([0] + [1 if x in sq else -1 for x in range(1, p)])
    C = np.zeros((p + 1, p + 1), dtype=np.int64)
    C[0, 1:] = C[1:, 0] = 1
    idx = np.arange(p)
    C[1:, 1:] = chi[(idx[None, :] - idx[:, None]) % p]
    return C, np.eye(p + 1, dtype=np.int64)


def paley_skew(p):
    """skew conference matrix of order p + 1, p = 3 (mod 4) a prime"""
    sq = {x * x % p for x in range(1, p)}
    chi = np.array([0] + [1 if x in sq else -1 for x in range(1, p)])
    C = np.zeros((p + 1, p + 1), dtype=np.int64)
    C[0, 1:], C[1:, 0] = 1, -1
    idx = np.arange(p)
    C[1:, 1:] = chi[(idx[None, :] - idx[:, None]) % p]
    return C


def double_skew(C):
    """skew conference matrix of order 2n from one of order n, through the skew Hadamard matrix [[H, H], [-H^T, H^T]]"""
    E = np.eye(C.shape[0], dtype=np.int64)
    H = E + C
    return np.block([[H, H], [-H.T, H.T]]) - np.eye(2 * C.shape[0], dtype=np.int64)


GENERAL = {15: (7, 1), 39: (19, 1), 63: (7, 3)}     # q: (p, number of doublings)


def general_C(q, seed=1):
    p, k = GENERAL[q]
    C = paley_skew(p)
    for _ in range(k):
        C = double_skew(C)
    rng = np.random.default_rng(seed)
    perm, s = rng.permutation(q + 1), rng.choice([-1, 1], q + 1)
    return (s[:, None] * C * s[None, :])[np.ix_(perm, perm)]


def file_K(name):
    cert = json.loads((ROOT / "complex-hadamard" / name).read_text())
    Kr, Ki = verify_complex.build(cert)
    return np.asarray(Kr, dtype=np.int64), np.asarray(Ki, dtype=np.int64)


def check(q, Kr, Ki, big=False, C=None):
    n = q - 1
    G = Kr @ Kr.T + Ki @ Ki.T
    assert np.array_equal(G, n * np.eye(n, dtype=np.int64)) and not np.any(Ki @ Kr.T - Kr @ Ki.T), "K is not Hadamard"
    Hr, Hi = verify.complex_miyamoto(q, Kr, Ki, C)
    Hr, Hi = np.asarray(Hr), np.asarray(Hi)
    assert Hr.shape == (4 * q, 4 * q)
    assert np.all(np.abs(Hr.astype(np.int64)) + np.abs(Hi.astype(np.int64)) == 1), "entries are not units"
    if not big:
        R, I = Hr.astype(np.float64), Hi.astype(np.float64)
        ok = (np.array_equal(R @ R.T + I @ I.T, 4 * q * np.eye(4 * q)) and not np.any(I @ R.T - R @ I.T))
    else:
        X, Y = Hr + Hi, Hr - Hi
        D = np.block([[X, Y], [-Y, X]]).astype(np.float32)
        del Hr, Hi, X, Y
        N, ok = D.shape[0], True
        for s in range(0, N, 2048):
            B = D[s:s + 2048] @ D.T
            B[np.arange(B.shape[0]), s + np.arange(B.shape[0])] -= N
            ok = ok and not np.any(B)
    return ok


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--big", action="store_true", help="also q = 1559 and 1619")
    ap.add_argument("--general", action="store_true", help="also q = 15, 39, 63 with non-Paley skew conference matrices")
    ap.add_argument("--q", type=int, nargs="*", help="only these q")
    a = ap.parse_args()
    qs = a.q or PALEY_Q + [223, 263, 331, 431] + ([1559, 1619] if a.big else [])
    bad = 0
    for q in qs:
        t = time.time()
        Kr, Ki = paley_K(q) if q in PALEY_Q or q not in FILE_Q else file_K(FILE_Q[q])
        ok = check(q, Kr, Ki, big=q > 1000)
        bad += not ok
        print(f"q = {q}: H H* = {4 * q}I {'holds' if ok else 'FAILS'} ({time.time() - t:.1f} s)", flush=True)
    if a.general:
        for q in GENERAL:
            t = time.time()
            Kr, Ki = paley_K(q)
            ok = check(q, Kr, Ki, C=general_C(q))
            bad += not ok
            print(f"q = {q} (skew conference matrix of order {q + 1}, not Paley's): H H* = {4 * q}I "
                  f"{'holds' if ok else 'FAILS'} ({time.time() - t:.1f} s)", flush=True)
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
