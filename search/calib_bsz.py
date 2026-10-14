#!/usr/bin/env python3
"""Calibration of the solution-count estimates (diagonal Gaussian bsz.gauss/ssxx.gauss vs hybrid gauss_mc) against
exact counts (complete lists + join4x) on small structures.

    python3 search/calib_bsz.py [--bordered] [--maxwork 2e9] [--mino 10] [--maxo 30] N [N ...]
"""
import argparse, math, sys, time
from math import gcd
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
from bnw import Struct, subgroups, enum_block, join4
import bsz, ssxx, gauss_mc

KINDS = ["L", "L", "A", "A"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("Ns", type=int, nargs="+"); ap.add_argument("--bordered", action="store_true")
    ap.add_argument("--maxwork", type=float, default=2e9)
    ap.add_argument("--mino", type=int, default=10); ap.add_argument("--maxo", type=int, default=30)
    a = ap.parse_args()
    for N in a.Ns:
        m = N - 1 if a.bordered else N
        for H, gens in subgroups(m).items():
            if (m - 1) in H:
                continue
            orb = sum(gcd(h - 1, m) for h in H) // len(H)
            if not (a.mino <= orb <= a.maxo):
                continue
            st = Struct(N, gens, plain=not a.bordered)
            pats = [p for _, p in bsz.patterns(st)] if a.bordered else ssxx.patterns(st)
            for pat in pats:
                sizes = [st.counts(k).get(r, 0) for k, r in zip(KINDS, pat)]
                if a.bordered:
                    sizes = [s // 2 if r == 0 else s for s, r in zip(sizes, pat)]
                if min(sizes) == 0:
                    continue
                order = min(([0, 1, 2, 3], [0, 2, 1, 3], [0, 3, 1, 2]),
                            key=lambda o: sizes[o[0]] * sizes[o[1]] + sizes[o[2]] * sizes[o[3]])
                work = sizes[order[0]] * sizes[order[1]] + sizes[order[2]] * sizes[order[3]]
                if work > a.maxwork:
                    continue
                d = bsz.gauss(st, pat, K=2000) if a.bordered else ssxx.gauss(st, pat, K=2000)
                h = gauss_mc.estimate(st, pat, a.bordered, Kb=5000, K=500000)
                hv = h[0] if (h is not None and not isinstance(h[0], str)) else (-99 if h else -999)
                t0 = time.time()
                lists = []
                for k, (kind, r) in enumerate(zip(KINDS, pat)):
                    V, M = enum_block(st, kind, r, 1 if (a.bordered and r == 0) else 0, f"cal{k}")
                    lists.append(V)
                sols, total = join4([lists[k] for k in order], st.C, 1, None, target=-4 if a.bordered else 0)
                print(f"{N:4d} {'b' if a.bordered else 'c'} H={str(gens):10s} |H|={len(st.H):3d} O={len(st.sizesA):3d} "
                      f"C={st.C:3d} rows {str(pat):16s} pools {[f'{s:.2g}' for s in sizes]} exact {total} "
                      f"(log {math.log10(total) if total else float('-inf'):.2f}) diag {d:.2f} hybrid {hv:.2f} "
                      f"[{time.time() - t0:.0f}s]", flush=True)


if __name__ == "__main__":
    main()
