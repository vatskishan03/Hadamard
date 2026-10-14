#!/usr/bin/env python3
"""Streaming k-tree search for Szöllősi structures on Z_N whose arbitrary blocks C, D are too large to list
(attack/kstream.c): A and B (symmetric, on the <H, -1>-orbits) are listed completely; C and D (on the H-orbits) are
enumerated by depth-first search, joined to A and to B on k classes, and matched on the remaining classes by
fingerprints. Every hit is rebuilt as Szöllősi's complex Hadamard matrix and checked exactly.

The pass-1 table holds the pairs (A, C) with a + c in the chosen target set on the k classes K1; the expected
number of solutions found is E times the probability that a solution's a + c falls in the target set.

    python3 attack/kstream.py N g a,b,c,d [--k 6] [--spread 2] [--wide 3] [--maxm1 2e8] [--threads 4] [--log FILE]
"""
import argparse, itertools, json, math, os, struct, subprocess, sys, time
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import workdir as W
import ssxx
from bnw import Struct, mask_to_y


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("N", type=int); ap.add_argument("gens"); ap.add_argument("pat")
    ap.add_argument("--k", type=int, default=6)
    ap.add_argument("--wide", type=int, default=3, help="classes of K1 that get two target values")
    ap.add_argument("--maxm1", type=float, default=2e8)
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--log")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--fracC", type=float, default=1.0)
    ap.add_argument("--fracD", type=float, default=1.0)
    ap.add_argument("--vals", type=int, default=2, help="target values per wide class")
    a = ap.parse_args()
    N = a.N; gens = [int(x) for x in a.gens.split(",")]; pat = [int(x) for x in a.pat.split(",")]
    st = Struct(N, gens, plain=True)
    C = st.C
    VA, MA = ssxx.enum(st, "L", pat[0]); VB, MB = ssxx.enum(st, "L", pat[1])
    OA = len(st.sizesA)
    print(f"N = {N}, H = <{a.gens}> |H| = {len(st.H)}: lists A {len(VA)}, B {len(VB)}; C, D on {OA} orbits "
          f"(pools {st.counts('A').get(pat[2], 0):.3g}, {st.counts('A').get(pat[3], 0):.3g}); {C} classes", flush=True)
    rng = np.random.default_rng(a.seed)
    VC = st.values("A", st.sample("A", pat[2], 4000, rng))
    mu = VA.mean(0) + VC.mean(0)
    res = (int(VA[0, 0]) + int(np.rint(VC[0, 0]))) % 4
    K1 = list(range(a.k))
    choices = []
    for q, c in enumerate(K1):
        # lattice values = res (mod 4) nearest the mean of a + c
        base = math.floor((mu[c] - res) / 4) * 4 + res
        cand = sorted([base + 4 * j for j in range(-3, 5)], key=lambda t: abs(t - mu[c]))
        choices.append(cand[:a.vals] if q < a.wide else cand[:1])
    targets = list(itertools.product(*choices))
    print(f"K1 = {K1}, {len(targets)} targets, e.g. {targets[0]}", flush=True)
    path = W.WORK / f"kstream_{N}_{os.getpid()}.in"
    with open(path, "wb") as f:
        f.write(struct.pack("8i", OA, C, a.k, len(targets), pat[2], pat[3], len(VA), len(VB)))
        f.write(np.asarray(st.sizesA, np.int32).tobytes())
        f.write(np.asarray(st.S["A"], np.int32).tobytes())
        f.write(np.asarray(K1, np.int32).tobytes())
        f.write(np.asarray(targets, np.int32).tobytes())
        f.write(np.zeros(C, np.int32).tobytes())
        f.write(np.ascontiguousarray(VA, np.int16).tobytes()); f.write(np.ascontiguousarray(VB, np.int16).tobytes())
    t0 = time.time()
    proc = subprocess.Popen([W.bin("kstream"), str(path), str(int(a.maxm1)), str(a.threads), str(a.fracC),
                             str(a.fracD)], stdout=subprocess.PIPE, text=True)
    found = 0
    for line in proc.stdout:
        line = line.strip()
        if line.startswith("STATS") or line.startswith("OVERFLOW"):
            print(line, f"[{time.time() - t0:.0f}s]", flush=True)
        if not line.startswith("SOL"):
            continue
        ai, bi, cc, dc = (int(x) for x in line.split()[1:])
        ya = mask_to_y(MA[ai], len(st.sizesL)); yb = mask_to_y(MB[bi], len(st.sizesL))
        yc = np.array([1 if (cc >> o) & 1 else -1 for o in range(OA)]); yd = np.array([1 if (dc >> o) & 1 else -1 for o in range(OA)])
        seqs = [st.expand("L", ya), st.expand("L", yb), st.expand("A", yc), st.expand("A", yd)]
        ssxx.complex_hadamard(seqs)
        found += 1
        rec = {"N": N, "gens": gens, "subgroup_order": len(st.H), "rowsums": [int(x.sum()) for x in seqs],
               "method": "kstream", "sequences": ["".join("+" if v > 0 else "-" for v in x) for x in seqs]}
        if a.log:
            with open(a.log, "a") as f:
                f.write(json.dumps(rec) + "\n")
        print(f"verified: complex Hadamard matrix of order {2 * N} (rows {rec['rowsums']}) [{time.time() - t0:.0f}s]",
              flush=True)
    proc.wait()
    path.unlink()
    print(f"done: {found} verified solutions [{time.time() - t0:.0f}s]", flush=True)


if __name__ == "__main__":
    main()
