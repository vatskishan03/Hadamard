#!/usr/bin/env python3
"""Gaussian estimate with the full covariance of the class values (ssxx.gauss uses only the diagonal).

    python3 search/gauss_full.py N g [g2] --rows a,b,c,d [--bordered] [--K 4000]
"""
import argparse, math, sys
from math import gcd, log10, pi
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
from bnw import Struct
import ssxx

KINDS = ["L", "L", "A", "A"]


def estimates(st, rows, bordered, K=4000, seed=1):
    rng = np.random.default_rng(seed)
    tg = -4 if bordered else 0
    lN, mus, covs, lats = 0.0, [], [], []
    for kind, r in zip(KINDS, rows):
        n = st.counts(kind).get(r, 0)
        if bordered and r == 0:
            n //= 2
        if n == 0:
            return None
        lN += log10(n)
        V = st.values(kind, st.sample(kind, r, K, rng))
        mus.append(V.mean(0)); covs.append(np.cov(V.T, bias=True))
        Vi = np.rint(V).astype(np.int64)
        lat = []
        for c in range(st.C):
            u = np.unique(Vi[:, c])
            g = int(np.gcd.reduce(np.diff(u))) if len(u) > 1 else 0
            lat.append((g, int(u[0] % g) if g else int(u[0])))
        lats.append(lat)
    mu = sum(mus); Sg = sum(covs); C = st.C
    var = np.diag(Sg).copy()
    keep = [c for c in range(C) if var[c] > 1e-9]
    lg = 0.0
    for c in range(C):
        g, res = 0, 0
        for lat in lats:
            g = gcd(g, lat[c][0]); res += lat[c][1]
        if var[c] <= 1e-9:
            if abs(tg - mu[c]) > 1e-6:
                return None
            continue
        if g and (tg - res) % g:
            return None
        lg += log10(max(g, 1))
    # drop the first kept class (linear dependency through the row sums), as ssxx.gauss does
    first = keep[0]
    c0 = keep[1:]
    g_first = 0
    for lat in lats:
        g_first = gcd(g_first, lat[first][0])
    lg -= log10(max(g_first, 1))
    d = (tg - mu)[c0]
    Sk = Sg[np.ix_(c0, c0)]
    sign, logdet = np.linalg.slogdet(Sk)
    quad = float(d @ np.linalg.solve(Sk, d))
    full = lN + lg + (-0.5 * len(c0) * math.log(2 * pi) - 0.5 * logdet - 0.5 * quad) / math.log(10)
    diag = lN + lg + sum(-0.5 * math.log(2 * pi * Sg[c, c]) - (tg - mu[c]) ** 2 / (2 * Sg[c, c]) for c in c0) / math.log(10)
    # correlation summary
    D = np.sqrt(np.diag(Sk)); R = Sk / np.outer(D, D)
    off = R[~np.eye(len(c0), dtype=bool)]
    return full, diag, float(np.abs(off).max()), float(off.mean())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("N", type=int); ap.add_argument("gens", type=int, nargs="+")
    ap.add_argument("--rows", required=True); ap.add_argument("--bordered", action="store_true")
    ap.add_argument("--K", type=int, default=4000)
    a = ap.parse_args()
    st = Struct(a.N, a.gens, plain=not a.bordered)
    rows = tuple(int(x) for x in a.rows.split(","))
    r = estimates(st, rows, a.bordered, a.K)
    if r is None:
        print("impossible"); return
    full, diag, mx, mean = r
    print(f"N = {a.N} H = <{a.gens}> rows {rows}: classes {st.C}, log10 E full-cov {full:.2f}, diagonal {diag:.2f}"
          f" (max |corr| {mx:.2f}, mean corr {mean:.3f})")


if __name__ == "__main__":
    main()
