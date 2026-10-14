#!/usr/bin/env python3
"""Certificates for Cooper-Wallis Hadamard matrices from a near-Williamson quadruple (search/nw_exh.c HIT line).

A HIT line gives circulant +-1 sequences A, B, C, D of order w with B, C, D symmetric and sum of PAFs = 4w delta.
Then (A R, B, C, D) are Williamson-type matrices (verify.py, build_cooper_wallis_near_williamson), and with
T-matrices of order t (from base sequences BS(m + 1, m), t = 2m + 1, found here by meet-in-the-middle) the
Cooper-Wallis construction gives a Hadamard matrix of order 4tw.

    python3 search/nw_certify.py HITFILE [--line K] --t 5 9 15 --out DIR [--write-csv]

writes DIR/<order>/certificate.json (construction "cooper-wallis-near-williamson"), checks the matrix with
verify.py (entries, H H^T = nI) and records the SHA-256 of its canonical CSV; with --write-csv also H<order>.csv.gz.
"""
import argparse, gzip, hashlib, itertools, json, sys
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "verify"))
import verify


def npaf(x):
    n = len(x)
    return np.array([int(np.dot(x[:n - s], x[s:])) for s in range(1, n)], dtype=np.int64)


def base_sequences(m):
    """BS(m + 1, m) by meet in the middle over (A, B) and (C, D), first entries +1."""
    L1 = {}
    for a in itertools.product((1, -1), repeat=m):
        a = np.array((1,) + a)
        for b in itertools.product((1, -1), repeat=m):
            b = np.array((1,) + b)
            L1.setdefault(tuple(npaf(a) + npaf(b)), (a, b))
    for c in itertools.product((1, -1), repeat=m - 1):
        c = np.array((1,) + c)
        for d in itertools.product((1, -1), repeat=m - 1):
            d = np.array((1,) + d)
            s = npaf(c) + npaf(d) if m > 1 else np.zeros(0, np.int64)
            key = tuple(-np.concatenate([s, [0]]))
            if key in L1:
                return L1[key] + (c, d)
    raise ValueError(f"no base sequences BS({m + 1}, {m})")


def t_matrices(t):
    m = (t - 1) // 2
    a, b, c, d = base_sequences(m)
    T = np.zeros((4, t), np.int64)
    T[0, :m + 1] = (a + b) // 2; T[1, :m + 1] = (a - b) // 2
    T[2, m + 1:] = (c + d) // 2; T[3, m + 1:] = (c - d) // 2
    assert np.all(np.count_nonzero(T, axis=0) == 1)
    paf = verify.paf_sum(list(T))
    assert paf[0] == t and not paf[1:].any()
    return ["".join({1: "+", -1: "-", 0: "0"}[int(x)] for x in row) for row in T]


def read_hit(path, k):
    hits = [ln for ln in Path(path).read_text().splitlines() if ln.startswith("HIT")]
    d = dict(kv.split("=", 1) for kv in hits[k].split()[1:])
    return [d[x] for x in "ABCD"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("hits"); ap.add_argument("--line", type=int, default=0)
    ap.add_argument("--t", type=int, nargs="+", required=True)
    ap.add_argument("--out", required=True); ap.add_argument("--write-csv", action="store_true")
    a = ap.parse_args()
    nw = read_hit(a.hits, a.line)
    w = len(nw[0])
    seq = [np.array([1 if ch == "+" else -1 for ch in s], np.int64) for s in nw]
    paf = verify.paf_sum(seq)
    assert paf[0] == 4 * w and not paf[1:].any(), "not a near-Williamson quadruple"
    assert all(np.array_equal(s, s[-np.arange(w) % w]) for s in seq[1:]), "B, C, D must be symmetric"
    print(f"near-Williamson quadruple of order {w}: row sums {[int(s.sum()) for s in seq]}")
    for t in a.t:
        cert = {"order": 4 * t * w, "construction": "cooper-wallis-near-williamson", "t": t, "w": w,
                "t_matrices": t_matrices(t), "near_williamson": nw}
        H, _ = verify.BUILDERS[cert["construction"]](cert)
        n = H.shape[0]
        assert n == 4 * t * w and verify.is_pm1(H) and verify.gram_is_scalar(H, n), "check failed"
        sha = hashlib.sha256()
        for chunk in verify.csv_chunks(H):
            sha.update(chunk)
        cert["sha256"] = sha.hexdigest()
        d = Path(a.out) / str(n); d.mkdir(parents=True, exist_ok=True)
        (d / "certificate.json").write_text(json.dumps(cert, indent=1) + "\n")
        if a.write_csv:
            with gzip.open(d / f"H{n}.csv.gz", "wb", compresslevel=9) as f:
                for chunk in verify.csv_chunks(H):
                    f.write(chunk)
        print(f"order {n} = 4 * {t} * {w}: H H^T = {n} I checked, sha256 {cert['sha256']}")


if __name__ == "__main__":
    main()
