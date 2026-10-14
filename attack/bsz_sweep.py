#!/usr/bin/env python3
"""Exhaustive sweep of bordered Szöllősi structures (attack/bsz.py, attack/bsz_group.py) for the open complex
Hadamard orders 2N: every multiplier subgroup H <= Z_(N-1)^* and, for N = 1 (mod 4), every spin structure on
Z_2^2 x Z_((N-1)/4), whose four complete lists can be joined within the pair-work budget, cheapest first.

    python3 attack/bsz_sweep.py N [N ...] [--maxpairs 2e10] [--budget SECONDS] [--log FILE]
"""
import argparse, math, sys, time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import bsz, bsz_group
from bnw import Struct, subgroups


def plan(N, maxpairs):
    jobs = []
    for H, gens in subgroups(N - 1).items():
        st = Struct(N, gens, plain=False)
        for name, pat in bsz.patterns(st):
            pools = [st.counts(k).get(r, 0) // (2 if r == 0 else 1) for k, r in zip(bsz.KINDS, pat)]
            if min(pools) == 0:
                continue
            s = sorted(pools)
            work = min(pools[o[0]] * pools[o[1]] + pools[o[2]] * pools[o[3]]
                       for o in ([0, 1, 2, 3], [0, 2, 1, 3], [0, 3, 1, 2]))
            enum = sum(2.0 ** len(st.sizes(k)) for k in bsz.KINDS)
            if work <= maxpairs and max(pools) <= 3e8 and enum <= 2 ** 31 and bsz.z2_ok(st, pat):
                jobs.append((work + enum, "cyc", gens, name, None))
    if (N - 1) % 4 == 0:
        for mods, mats, tag, O in bsz_group.spin_structures(N, maxo=64):
            bs = bsz_group.BStruct(N, mods, mats, tag)
            for name, pat in bs.patterns():
                pools = [bs.count(k, r) // (2 if r == 0 else 1) for k, r in zip(bsz.KINDS, pat)]
                if min(pools) == 0:
                    continue
                work = min(pools[o[0]] * pools[o[1]] + pools[o[2]] * pools[o[3]]
                           for o in ([0, 1, 2, 3], [0, 2, 1, 3], [0, 3, 1, 2]))
                enum = 2 * 2.0 ** bs.stL.O + 2 * 2.0 ** bs.stA.O
                if work <= maxpairs and max(pools) <= 3e8 and enum <= 2 ** 31:
                    jobs.append((work + enum, "spin", (mods, mats, tag), name, None))
    jobs.sort(key=lambda j: j[0])
    return jobs


def run(N, maxpairs, budget, log):
    t0 = time.time()
    jobs = plan(N, maxpairs)
    print(f"N = {N}: {len(jobs)} bordered structure/pattern jobs within {maxpairs:.2g} pair work", flush=True)
    for work, kind, spec, name, _ in jobs:
        if time.time() - t0 > budget:
            print(f"N = {N}: budget exhausted", flush=True); break
        try:
            if kind == "cyc":
                hit = bsz.search(N, list(spec), log=log, maxpairs=maxpairs, pats=[name])
            else:
                mods, mats, tag = spec
                bs = bsz_group.BStruct(N, mods, mats, tag)
                print(f"N = {N}, spin {tag}, 2 in {name}: O_A = {bs.stA.O}, O_L = {bs.stL.O}, classes = {bs.C}",
                      flush=True)
                hit = bsz_group.search(bs, log=log, maxpairs=maxpairs)
        except Exception as ex:
            print(f"  error: {ex!r}", flush=True); continue
        if hit:
            print(f"N = {N}: FOUND", flush=True)
            return True
    print(f"N = {N}: nothing", flush=True)
    return False


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("N", type=int, nargs="+")
    ap.add_argument("--maxpairs", type=float, default=2e10)
    ap.add_argument("--budget", type=float, default=1800)
    ap.add_argument("--log")
    a = ap.parse_args()
    for N in a.N:
        run(N, a.maxpairs, a.budget, a.log)


if __name__ == "__main__":
    main()
