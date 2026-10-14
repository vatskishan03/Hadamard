#!/usr/bin/env python3
"""Compression-chain sweep (attack/cchain4.py) over the cyclic Szöllősi structures that the exact sweep
(attack/ssxx_sweep.py) skips: every H <= Z_N^* with -1 not in H, at most 46 orbits and pair work above the exact
budget, fewest orbits first, every row-sum pattern, with time budgets per pattern and per N.

    python3 attack/cchain_sweep.py N [N ...] [--minwork 2e10] [--maxo 46] [--patbudget 300] [--budget 3600] [--log FILE]
"""
import argparse, contextlib, io, sys, time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import cchain4, ssxx, bsz
from bnw import Struct, subgroups


def structures(N, minwork, maxo, mode="ssxx"):
    out = []
    v = N if mode == "ssxx" else N - 1
    for H, gens in subgroups(v).items():
        if mode == "ssxx" and (N - 1) in H:
            continue
        st = Struct(N, gens, plain=(mode == "ssxx"))
        OA, OL = len(st.sizesA), len(st.sizesL)
        if OA > maxo:
            continue
        pats = []
        allp = ssxx.patterns(st) if mode == "ssxx" else [p for _, p in bsz.patterns(st) if bsz.z2_ok(st, p)]
        for p in allp:
            sizes = [st.counts(k).get(r, 0) for k, r in zip(ssxx.KINDS, p)]
            if min(sizes) == 0:
                continue
            w = ssxx.pair_cost(sizes, ssxx.best_pairing(sizes))
            if w > minwork:
                pats.append((w, p))
        if pats:
            pats.sort()
            out.append((OA, len(H), gens, [p for _, p in pats]))
    out.sort()
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("N", type=int, nargs="+")
    ap.add_argument("--minwork", type=float, default=2e10)
    ap.add_argument("--maxo", type=int, default=46)
    ap.add_argument("--patbudget", type=float, default=300)
    ap.add_argument("--budget", type=float, default=3600)
    ap.add_argument("--log")
    ap.add_argument("--mode", default="ssxx", choices=["ssxx", "bsz"])
    a = ap.parse_args()
    for N in a.N:
        t0 = time.time()
        sts = structures(N, a.minwork, a.maxo, a.mode)
        print(f"N = {N}: {len(sts)} structures for the chain", flush=True)
        done = False
        for OA, h, gens, pats in sts:
            for pat in pats:
                if time.time() - t0 > a.budget:
                    break
                buf = io.StringIO()
                try:
                    with contextlib.redirect_stdout(buf):
                        r = cchain4.run(a.mode, N, list(gens), tuple(int(x) for x in pat), budget=a.patbudget,
                                        log=a.log)
                except Exception as ex:
                    print(f"  H = <{gens}> rows {pat}: error {ex!r}", flush=True); continue
                res = [ln for ln in buf.getvalue().splitlines() if "RESULT" in ln or "too large" in ln]
                print(f"  |H| = {h} <{','.join(map(str, gens))}> O_A = {OA} rows {tuple(pat)}: "
                      f"{res[-1].split(': ', 1)[-1] if res else buf.getvalue().strip()[-120:]}", flush=True)
                if r.get("found"):
                    print(buf.getvalue(), flush=True)
                    print(f"N = {N}: FOUND", flush=True)
                    done = True; break
            if done or time.time() - t0 > a.budget:
                break
        if not done:
            print(f"N = {N}: nothing ({time.time() - t0:.0f}s)", flush=True)


if __name__ == "__main__":
    main()
