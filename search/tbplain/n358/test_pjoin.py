#!/usr/bin/env python3
"""pjoin against a numpy brute force: 12 random small instances (4 lists, 5 to 8 columns, planted and noisy targets), all quadruples.

For every instance and for 1, 8, 16 (in 3 parts) and 64 (in 2 parts) key buckets (bucket ranges split between processes)
the histogram of the agreeing-prefix length that pjoin reports must equal the brute-force histogram, for every prefix
length from K0 to m.  Needs ./pjoin (gcc -O2 -march=native -o pjoin pjoin.c).

    python3 test_pjoin.py
"""
import os, sys, tempfile
import numpy as np
from pio import write_lists, run_pjoin, read_out

rng = np.random.default_rng(5)
tmp = tempfile.mkdtemp(prefix="pjtest_")
fn, out = os.path.join(tmp, "pj.bin"), os.path.join(tmp, "pj.out")
ok = True
for trial in range(12):
    m = int(rng.integers(5, 9)); k0 = int(rng.integers(2, 5)); n = [int(rng.integers(25, 45)) for _ in range(4)]
    V = [rng.integers(-2, 3, size=(nn, m)) * (1 if rng.random() < 0.5 else 2) for nn in n]
    # target near the middle, sometimes planted
    pl = [int(rng.integers(nn)) for nn in n]
    T = sum(V[k][pl[k]] for k in range(4)) if trial % 3 == 0 else np.rint(sum(v.mean(0) for v in V)).astype(int)
    # brute force: histogram of the agreeing prefix length over all quadruples
    S01 = (V[0][:, None, :] + V[1][None, :, :]).reshape(-1, m)
    S23 = (V[2][:, None, :] + V[3][None, :, :]).reshape(-1, m)
    hist = np.zeros(m + 2, dtype=np.int64)
    for s in S01:
        d = (S23 + s[None, :]) == T[None, :]
        pref = np.cumprod(d, axis=1).sum(1)           # prefix length
        hist += np.bincount(pref, minlength=m + 2)
    for nbits, nparts in ((0, 1), (3, 1), (4, 3), (6, 2)):
        write_lists(fn, V, T)
        tot = {}
        for part in range(nparts):
            run_pjoin(fn, k0, nbits, trial + 1, part, nparts, m, out)
            h, q, st = read_out(out, k0, m)
            for k, v in h.items(): tot[k] = tot.get(k, 0) + v
        for l in range(k0, m + 1):
            if tot.get(l, 0) != hist[l]:
                ok = False; print("MISMATCH trial", trial, "nbits", nbits, "l", l, tot.get(l), hist[l])
for f in (fn, out): os.path.exists(f) and os.remove(f)
os.rmdir(tmp)
print("ALL OK" if ok else "FAILED")
sys.exit(0 if ok else 1)
