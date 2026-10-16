#!/usr/bin/env python3
"""Exact check of a quaternary Hadamard matrix in the text format of Zaitsev's page (one row per line, characters
1, -, i, j for 1, -1, i, -i), and of its real double.

Turyn's doubling: if K = A + iB has entries 1, -1, i, -i and K K* = m I, then A Aᵀ + B Bᵀ = m I and B Aᵀ = A Bᵀ, so
with X = A + B and Y = A - B (both ±1 matrices) H = [[X, Y], [-Y, X]] satisfies H Hᵀ = 2m I.

    python3 -I search/quaternary_check.py QH1142.txt            # K K* = 1142 I
    python3 -I search/quaternary_check.py QH1142.txt --double   # and H Hᵀ = 2284 I
"""
import sys
import numpy as np

ENTRY = {'1': (1, 0), '-': (-1, 0), 'i': (0, 1), 'j': (0, -1)}


def load(path):
    rows = [line.strip() for line in open(path) if line.strip()]
    if any(len(r) != len(rows) or set(r) - set(ENTRY) for r in rows):
        raise SystemExit(f"{path}: not a square matrix over 1, -, i, j")
    A = np.array([[ENTRY[c][0] for c in r] for r in rows], np.int64)
    B = np.array([[ENTRY[c][1] for c in r] for r in rows], np.int64)
    return A, B


def main():
    path = sys.argv[1]
    A, B = load(path)
    m = A.shape[0]
    eye = np.eye(m, dtype=np.int64)
    ok = np.array_equal(A @ A.T + B @ B.T, m * eye) and np.array_equal(B @ A.T, A @ B.T)
    print(f"{path}: quaternary Hadamard of order {m}: {'ok' if ok else 'FAILED'}")
    if ok and "--double" in sys.argv:
        X, Y = A + B, A - B
        H = np.block([[X, Y], [-Y, X]])
        good = np.array_equal(H @ H.T, 2 * m * np.eye(2 * m, dtype=np.int64))
        print(f"real double of order {2 * m}: {'ok' if good else 'FAILED'}")
        ok = good
    raise SystemExit(0 if ok else 1)


if __name__ == "__main__":
    main()
