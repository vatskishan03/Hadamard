#!/usr/bin/env python3
"""Frontier map for the low-energy k-tree: every multiplier structure (cyclic Z_N and bordered Z_(N-1)) of an open
order whose orbit count lies in [MINO, MAXO], with the Gaussian estimate E of every row pattern and the pool sizes.
Score = log10 E - (log10 P_C + log10 P_D) / 2 + 7.5, the rough log10 of the solutions inside low-energy lists of
3e7 members (retention ~ (L/P)^0.5 per list, fitted on N = 171, 183, 213, 341). Per-(kind, row) statistics are cached.

    python3 search/frontier_scan.py N [N ...] [--mino 20] [--maxo 75] [--maxpool 1e13] [--out FILE]
"""
import argparse, json, math, sys
from math import gcd, log10, pi
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
from bnw import Struct, subgroups
import ssxx, bsz

KINDS = ["L", "L", "A", "A"]


def stats(st, kind, r, K, rng, cache):
    key = (kind, r)
    if key not in cache:
        n = st.counts(kind).get(r, 0)
        if n == 0:
            cache[key] = None
        else:
            V = st.values(kind, st.sample(kind, r, K, rng))
            Vi = np.rint(V).astype(np.int64)
            lat = []
            for c in range(st.C):
                u = np.unique(Vi[:, c])
                g = int(np.gcd.reduce(np.diff(u))) if len(u) > 1 else 0
                lat.append((g, int(u[0] % g) if g else int(u[0])))
            cache[key] = (n, V.mean(0), V.var(0), lat)
    return cache[key]


def gauss(st, rows, bordered, K, rng, cache):
    tg = -4 if bordered else 0
    lN, mus, vs, lats = 0.0, [], [], []
    for kind, r in zip(KINDS, rows):
        s = stats(st, kind, r, K, rng, cache)
        if s is None:
            return -math.inf
        n, mu, var, lat = s
        if bordered and r == 0:
            n //= 2
        lN += log10(n); mus.append(mu); vs.append(var); lats.append(lat)
    lp, dropped = 0.0, False
    for c in range(st.C):
        mu = sum(x[c] for x in mus); var = sum(x[c] for x in vs)
        g, res = 0, 0
        for lat in lats:
            g = gcd(g, lat[c][0]); res += lat[c][1]
        if var < 1e-9:
            if abs(tg - mu) > 1e-6:
                return -math.inf
            continue
        if g and (tg - res) % g:
            return -math.inf
        if not dropped:
            dropped = True; continue
        lp += -0.5 * log10(2 * pi * var) - (tg - mu) ** 2 / (2 * var) / math.log(10) + log10(max(g, 1))
    return lN + lp


def scan(N, mino, maxo, maxpool, K, out):
    res = []
    for bordered in (False, True):
        m = N - 1 if bordered else N
        for H, gens in subgroups(m).items():
            if (m - 1) in H:
                continue
            orb = sum(gcd(h - 1, m) for h in H) // len(H)
            if not (mino <= orb <= maxo):
                continue
            st = Struct(N, gens, plain=not bordered)
            pats = [p for _, p in bsz.patterns(st)] if bordered else ssxx.patterns(st)
            rng = np.random.default_rng(1); cache = {}
            for pat in pats:
                sizes = [st.counts(k).get(r, 0) for k, r in zip(KINDS, pat)]
                if min(sizes) == 0 or max(sizes[2:]) > maxpool or max(sizes[:2]) > 1e8:
                    continue
                e = gauss(st, pat, bordered, K, rng, cache)
                if e == -math.inf:
                    continue
                work = min(sizes[o[0]] * sizes[o[1]] + sizes[o[2]] * sizes[o[3]]
                           for o in ([0, 1, 2, 3], [0, 2, 1, 3], [0, 3, 1, 2]))
                score = e - 0.5 * (log10(sizes[2]) + log10(sizes[3])) + 7.5
                r = {"N": N, "bordered": bordered, "gens": list(gens), "H": len(st.H), "OA": len(st.sizesA),
                     "OL": len(st.sizesL), "C": st.C, "rows": list(map(int, pat)), "logE": round(float(e), 2),
                     "score": round(float(score), 2), "work": float(work), "sizes": [float(s) for s in sizes]}
                res.append(r)
                if out:
                    out.write(json.dumps(r) + "\n"); out.flush()
    res.sort(key=lambda r: -r["score"])
    best = [(r["score"], r["logE"], "b" if r["bordered"] else "c", r["gens"], r["OA"], r["rows"],
             "%.2g/%.2g" % (r["sizes"][2], r["sizes"][3]), "exact" if r["work"] <= 2e10 else "") for r in res[:3]]
    print(N, best, flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("Ns", type=int, nargs="+")
    ap.add_argument("--mino", type=int, default=20)
    ap.add_argument("--maxo", type=int, default=75)
    ap.add_argument("--maxpool", type=float, default=1e13)
    ap.add_argument("--K", type=int, default=300)
    ap.add_argument("--out")
    a = ap.parse_args()
    out = open(a.out, "a") if a.out else None
    for N in a.Ns:
        scan(N, a.mino, a.maxo, a.maxpool, a.K, out)


if __name__ == "__main__":
    main()
