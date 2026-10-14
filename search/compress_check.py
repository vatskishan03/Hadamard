#!/usr/bin/env python3
"""Necessary conditions from the compressions: for every prime-power divisor d of the group order, join the four
blocks' compressions to Z_d exactly (search/cchain4.py, level 1). An empty join rules the structure/pattern out.
A level whose state space exceeds maxstates is skipped (None); --dmax skips the levels above it (the large prime
levels take minutes and gigabytes, the small ones milliseconds).

    python3 search/compress_check.py ssxx N g[,g2] a,b,c,d [--dmax D]     (cyclic Szöllősi)
    python3 search/compress_check.py bsz N g[,g2] r1,r2,r3,r4 [--dmax D]  (bordered, rows like 0,0,2,0)
"""
import argparse, sys, time
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
import cchain4 as CC
import chain as CH


def check(mode, N, gens, pat, maxstates=2_000_000, dmax=None, stop=False):
    """{d: number of level-1 solutions, 0 (ruled out), None (state space too large) or "large" (join too large)};
    stop: return at the first empty level"""
    kinds, T, bordered = CC.MODES[mode]
    v = N - 1 if bordered else N
    ps = [p for p in range(2, v) if v % p == 0 and all(p % r for r in range(2, int(p ** 0.5) + 1))]
    cand = sorted({p ** e for p in ps for e in range(1, 8) if v % (p ** e) == 0 and p ** e < v})
    out = {}
    for d in cand:
        if dmax and d > dmax:
            continue
        c1 = CC.KChain(v, gens, [d], T, kinds)
        Sd = [c1.states(k, r, maxstates) for k, r in zip(kinds, pat)]
        if any(x is None for x in Sd):
            out[d] = None; continue
        U = [np.unique(K_, axis=0).astype(np.int64) for K_, _ in Sd]
        try:
            out[d] = 0 if any(len(u) == 0 for u in U) else len(CH.join4([c1.paf(0, u) for u in U], c1.target(0)))
        except CH.JoinTooLarge:
            out[d] = "large"
        if stop and out[d] == 0:
            break
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=["ssxx", "bsz"]); ap.add_argument("N", type=int)
    ap.add_argument("gens"); ap.add_argument("pat")
    ap.add_argument("--dmax", type=int); ap.add_argument("--maxstates", type=float, default=2e6)
    a = ap.parse_args()
    mode, N, gens, pat = a.mode, a.N, [int(x) for x in a.gens.split(",")], tuple(int(x) for x in a.pat.split(","))
    t0 = time.time()
    r = check(mode, N, gens, pat, int(a.maxstates), a.dmax)
    dead = [d for d, n in r.items() if n == 0]
    print(f"{mode} N = {N} H = <{gens}> rows {pat}: level-1 joins {r} -> "
          f"{'RULED OUT at d = ' + str(dead[0]) if dead else 'not excluded'} [{time.time() - t0:.0f}s]")


if __name__ == "__main__":
    main()
