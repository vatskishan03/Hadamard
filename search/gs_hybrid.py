#!/usr/bin/env python3
"""Goethals-Seidel quadruples (real Hadamard matrices of order 4N) in multiplier structures: all four blocks constant
on the H-orbits of Z_N, sum of PAFs 0 off zero. Diagonal Gaussian vs hybrid (gauss_mc) estimates per row pattern.

    python3 search/gs_hybrid.py N [--mino 20] [--maxo 60] [--maxpool 1e11]
"""
import argparse, math, sys
from math import gcd, log10, pi
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
from bnw import Struct, subgroups
import gauss_mc

K4 = ["A", "A", "A", "A"]


def gs_patterns(st):
    N = st.N; lim = math.isqrt(4 * N)
    ok = sorted(r for r in st.counts("A") if 0 < r <= lim)
    out = []
    for a in ok:
        for b in ok:
            if b < a: continue
            for c in ok:
                if c < b: continue
                d2 = 4 * N - a * a - b * b - c * c
                if d2 < c * c: continue
                d = math.isqrt(d2)
                if d * d == d2 and d in ok:
                    out.append((a, b, c, d))
    return out


def diag(st, rows, K=3000, seed=1):
    rng = np.random.default_rng(seed); lN = 0.0; mus = []; vs = []; lats = []
    for r in rows:
        n = st.counts("A").get(r, 0)
        if n == 0: return -math.inf
        lN += log10(n)
        V = st.values("A", st.sample("A", r, K, rng)); mus.append(V.mean(0)); vs.append(V.var(0))
        Vi = np.rint(V).astype(np.int64); lat = []
        for c in range(st.C):
            u = np.unique(Vi[:, c]); g = int(np.gcd.reduce(np.diff(u))) if len(u) > 1 else 0
            lat.append((g, int(u[0] % g) if g else int(u[0])))
        lats.append(lat)
    lp, dropped = 0.0, False
    for c in range(st.C):
        mu = sum(x[c] for x in mus); var = sum(x[c] for x in vs); g = 0; res = 0
        for lat in lats: g = gcd(g, lat[c][0]); res += lat[c][1]
        if var < 1e-9:
            if abs(mu) > 1e-6: return -math.inf
            continue
        if g and (-res) % g: return -math.inf
        if not dropped: dropped = True; continue
        lp += -0.5 * log10(2 * pi * var) - mu ** 2 / (2 * var) / math.log(10) + log10(max(g, 1))
    return lN + lp


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("Ns", type=int, nargs="+"); ap.add_argument("--mino", type=int, default=20)
    ap.add_argument("--maxo", type=int, default=60); ap.add_argument("--maxpool", type=float, default=1e11)
    ap.add_argument("--top", type=int, default=4)
    a = ap.parse_args()
    for N in a.Ns:
        for H, gens in subgroups(N).items():
            if (N - 1) in H:          # GS blocks need not be symmetric; -1 in H forces symmetric blocks (Williamson)
                continue
            orb = sum(gcd(h - 1, N) for h in H) // len(H)
            if not (a.mino <= orb <= a.maxo): continue
            st = Struct(N, gens, plain=True)
            pats = gs_patterns(st)
            res = []
            for pat in pats:
                pools = [st.counts("A").get(r, 0) for r in pat]
                if max(pools) > a.maxpool: continue
                d = diag(st, pat, K=1500)
                res.append((d, pat, pools))
            res.sort(key=lambda x: -x[0])
            for d, pat, pools in res[:a.top]:
                h = gauss_mc.estimate(st, pat, False, Kb=8000, K=600000, kinds=K4)
                hv = None if (h is None or isinstance(h[0], str)) else round(float(h[0]), 2)
                print(f"N = {N} H = <{gens}> |H| = {len(st.H)} orbits {len(st.sizesA)} classes {st.C} rows {pat}: "
                      f"diag {d:.2f} hybrid {hv} pools {[f'{p:.2g}' for p in pools]}", flush=True)


if __name__ == "__main__":
    main()
