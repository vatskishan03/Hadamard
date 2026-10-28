#!/usr/bin/env python3
"""Input and output of pjoin: write the four lists and the target, run pjoin, read its histogram and stored quadruples.

File written by `write_lists(fn, [V0, V1, V2, V3], T)`: int32 m, n0, n1, n2, n3, then the four int16 matrices (n_k x m, row major)
and the int32 target T[m].  `run_pjoin(fn, k0, nbits, seed, part, nparts, lstore, out)` runs ./pjoin (same directory) with these
arguments (see pjoin.c); `read_out(out, k0, m)` returns (histogram {prefix length: count}, [(l, a, b, c, d), ...] of the quadruples
with prefix length >= lstore, the STATS line).  Library module (engine.py, test_pjoin.py, small_analogue.py).
"""
import numpy as np, subprocess, os, struct

def write_lists(fn, V, T):
    m = V[0].shape[1]
    with open(fn, "wb") as f:
        f.write(struct.pack("5i", m, *[len(v) for v in V]))
        for v in V:
            f.write(np.ascontiguousarray(v, dtype=np.int16).tobytes())
        f.write(np.asarray(T, dtype=np.int32).tobytes())

def run_pjoin(fn, k0, nbits, seed, part, nparts, lstore, out, exe=None):
    exe = exe or os.path.join(os.path.dirname(os.path.abspath(__file__)), "pjoin")
    subprocess.run([exe, fn, str(k0), str(nbits), str(seed), str(part), str(nparts), str(lstore), out], check=True)

def read_out(out, k0, m):
    hist = {}; quads = []; stats = None
    for l in open(out):
        if l.startswith("HIST"):
            vals = list(map(int, l.split()[1:]))
            for i, v in enumerate(vals): hist[k0 + i] = v
        elif l.startswith("Q"):
            quads.append(tuple(map(int, l.split()[1:])))
        elif l.startswith("# STATS"):
            stats = l.strip()
    return hist, quads, stats
