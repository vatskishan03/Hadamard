#!/usr/bin/env python3
"""Rebuild and check the Hadamard matrices in this repository.

Every order directory holds a certificate.json. For each certificate this script
  1. expands the construction data,
  2. checks the defining identity (periodic autocorrelations of the four sequences
     sum to zero off the origin),
  3. assembles the matrix H and checks H H^T = nI, and
  4. compares the SHA-256 of the canonical CSV encoding with the recorded value.

The canonical encoding is one row per line, entries 1 / -1 separated by commas, each
line terminated by "\\n" (numpy.savetxt with fmt="%d", delimiter=",").

Products are formed in float64. Every entry of H H^T is an integer of absolute value
at most n < 2^53, so the computation is exact.

    python3 verify/verify.py              check every order in the repository
    python3 verify/verify.py 2092 4492    check the named orders
    python3 verify/verify.py --write      also write H<n>.csv next to each certificate
"""
import argparse
import hashlib
import io
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
SIGN = {"+": 1, "-": -1, "0": 0}


def circulant(row):
    """Circulant matrix with first row `row`: A[i, j] = row[(j - i) mod n]."""
    row = np.asarray(row, dtype=np.int64)
    n = len(row)
    idx = np.arange(n)
    return row[(idx[None, :] - idx[:, None]) % n]


def back_identity(n):
    R = np.zeros((n, n), dtype=np.int64)
    R[np.arange(n), n - 1 - np.arange(n)] = 1
    return R


def goethals_seidel(A, B, C, D):
    R = back_identity(A.shape[0])
    return np.block([[A, B @ R, C @ R, D @ R],
                     [-B @ R, A, D.T @ R, -C.T @ R],
                     [-C @ R, -D.T @ R, A, B.T @ R],
                     [-D @ R, C.T @ R, -B.T @ R, A]])


def paf_sum(rows):
    """Sum over the given sequences of their periodic autocorrelation functions
    (exact integer arithmetic)."""
    n = len(rows[0])
    total = np.zeros(n, dtype=np.int64)
    for x in rows:
        x = np.asarray(x, dtype=np.int64)
        total += np.array([x @ np.roll(x, -s) for s in range(n)], dtype=np.int64)
    return total


def build_goethals_seidel(cert):
    v, g = cert["v"], cert["multiplier"]
    H = sorted({pow(g, i, v) for i in range(v)})
    if len(H) != cert["subgroup_order"]:
        raise ValueError("multiplier does not generate a subgroup of the stated order")
    seqs = []
    for reps in cert["blocks"]:
        x = np.ones(v, dtype=np.int64)
        for r in reps:
            x[[r * h % v for h in H]] = -1
        seqs.append(x)
    if [int(x.sum()) for x in seqs] != cert["rowsums"]:
        raise ValueError("row sums differ from the certificate")
    mu = cert.get("spin_multiplier")
    if mu is not None:
        t = np.arange(v)
        if not all(np.array_equal(seqs[k + 1], seqs[1][t * pow(mu, k, v) % v]) for k in range(3)):
            raise ValueError("the last three blocks are not multiplier images of each other")
    paf = paf_sum(seqs)
    if paf[0] != 4 * v or np.any(paf[1:] != 0):
        raise ValueError("sequences are not complementary")
    return goethals_seidel(*map(circulant, seqs)), seqs


def build_cooper_wallis(cert):
    """T-matrices of order t (Cooper-Wallis) give OD(4t; t, t, t, t); substituting
    Williamson matrices of order w for its four variables gives H of order 4tw."""
    T = [np.array([SIGN[c] for c in s], dtype=np.int64) for s in cert["t_matrices"]]
    t = len(T[0])
    if not np.all(np.count_nonzero(np.array(T), axis=0) == 1):
        raise ValueError("T-matrix supports do not partition Z_t")
    paf = paf_sum(T)
    if paf[0] != t or np.any(paf[1:] != 0):
        raise ValueError("T-matrix condition fails")
    W = [circulant([SIGN[c] for c in s]) for s in cert["williamson"]]
    w = W[0].shape[0]
    if not all(np.array_equal(m, m.T) for m in W) or \
            not np.array_equal(sum(m @ m for m in W), 4 * w * np.eye(w, dtype=np.int64)):
        raise ValueError("Williamson condition fails")
    X = [circulant(x) for x in T]
    coeff = {"A": (X[0], X[1], X[2], X[3]), "B": (X[1], -X[0], -X[3], X[2]),
             "C": (X[2], X[3], -X[0], -X[1]), "D": (X[3], -X[2], X[1], -X[0])}
    H = 0
    for k in range(4):
        H = H + np.kron(goethals_seidel(*(coeff[c][k] for c in "ABCD")), W[k])
    return H, None


def build_spence(cert):
    """Spence's bordered construction. Four +-1 sequences on an abelian group G of order v
    (elements in mixed-radix order of `group`) whose autocorrelations sum to -4 off the
    identity, with row sums (2, 0, 0, 0), give a Hadamard matrix of order 4v + 4:

        H = [ K        Q (x) 1^T ]
            [ P (x) 1  M         ]

    where M is the group-developed Goethals-Seidel array, P and Q are the 4 x 4 Hadamard
    matrices in the certificate, and K = -QP/2."""
    mods = cert["group"]
    v = int(np.prod(mods))
    coords = np.array(np.unravel_index(np.arange(v), mods)).T
    index = lambda c: np.ravel_multi_index(tuple(np.moveaxis(c % mods, -1, 0)), mods)
    seqs = [np.array([SIGN[c] for c in s], dtype=np.int64) for s in cert["sequences"]]
    if [int(x.sum()) for x in seqs] != [2, 0, 0, 0]:
        raise ValueError("row sums must be (2, 0, 0, 0)")
    total = np.zeros(v, dtype=np.int64)
    for s in range(v):
        shift = index(coords + coords[s])
        total[s] = sum(int(x @ x[shift]) for x in seqs)
    if total[0] != 4 * v or np.any(total[1:] != -4):
        raise ValueError("sequences do not satisfy the Spence condition")
    sub = index(coords[None, :, :] - coords[:, None, :])          # X[g, h] = x(h - g)
    R = np.zeros((v, v), dtype=np.int64)
    R[np.arange(v), index(-coords)] = 1
    A, B, C, D = (x[sub] for x in seqs)
    M = np.block([[A, B @ R, C @ R, D @ R],
                  [-B @ R, A, D.T @ R, -C.T @ R],
                  [-C @ R, -D.T @ R, A, B.T @ R],
                  [-D @ R, C.T @ R, -B.T @ R, A]])
    P, Q = np.array(cert["P"], dtype=np.int64), np.array(cert["Q"], dtype=np.int64)
    one = np.ones((1, v), dtype=np.int64)
    K = -(Q @ P) // 2
    return np.vstack([np.hstack([K, np.kron(Q, one)]), np.hstack([np.kron(P, one.T), M])]), None


BUILDERS = {"goethals-seidel": build_goethals_seidel, "cooper-wallis": build_cooper_wallis,
            "spence": build_spence}


def canonical_csv(H):
    buf = io.BytesIO()
    np.savetxt(buf, H, fmt="%d", delimiter=",")
    return buf.getvalue()


def check(directory, write=False):
    cert = json.loads((directory / "certificate.json").read_text())
    H, seqs = BUILDERS[cert["construction"]](cert)
    n = H.shape[0]
    if n != cert["order"]:
        raise ValueError(f"built order {n}, certificate says {cert['order']}")
    if not np.isin(H, (-1, 1)).all():
        raise ValueError("entries are not all +-1")
    Hf = H.astype(np.float64)
    if not np.array_equal(Hf @ Hf.T, n * np.eye(n)):
        raise ValueError("H H^T != nI")
    seq_file = directory / "sequences.txt"
    if seqs is not None and seq_file.exists():
        rows = [ln.strip() for ln in seq_file.read_text().splitlines() if ln.strip()]
        if rows != ["".join("+" if v > 0 else "-" for v in x) for x in seqs]:
            raise ValueError("sequences.txt does not match the certificate")
    data = canonical_csv(H)
    digest = hashlib.sha256(data).hexdigest()
    if digest != cert["sha256"]:
        raise ValueError(f"SHA-256 mismatch: {digest}")
    if write:
        (directory / f"H{n}.csv").write_bytes(data)
    return n, digest


def check_t_matrices(path):
    T = [[SIGN[c] for c in ln.strip()] for ln in path.read_text().splitlines() if ln.strip()]
    t = len(T[0])
    if len(T) != 4 or any(len(r) != t for r in T):
        raise ValueError("expected four rows of equal length")
    if not np.all(np.count_nonzero(np.array(T), axis=0) == 1):
        raise ValueError("supports do not partition Z_t")
    paf = paf_sum(T)
    if paf[0] != t or np.any(paf[1:] != 0):
        raise ValueError("autocorrelation condition fails")
    return t


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("orders", nargs="*", help="order directories to check (default: all)")
    ap.add_argument("--write", action="store_true", help="write H<n>.csv next to each certificate")
    args = ap.parse_args()
    dirs = [ROOT / o for o in args.orders] if args.orders else \
        sorted((p.parent for p in ROOT.glob("*/certificate.json")), key=lambda p: int(p.name))
    failed = 0
    for d in dirs:
        try:
            n, digest = check(d, args.write)
            print(f"{d.name:>6}  ok    H H^T = {n} I   sha256 {digest}")
        except Exception as exc:  # report and continue with the remaining orders
            failed += 1
            print(f"{d.name:>6}  FAIL  {exc}")
    if not args.orders:
        for f in sorted((ROOT / "t-matrices").glob("T*.txt"), key=lambda p: int(p.stem[1:])):
            try:
                print(f"{f.stem:>6}  ok    T-matrices of order {check_t_matrices(f)}")
            except Exception as exc:
                failed += 1
                print(f"{f.stem:>6}  FAIL  {exc}")
    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()
