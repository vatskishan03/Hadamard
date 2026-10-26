"""Control of join4.c: the exact four-list join against a numpy brute force on random lattice-structured instances.

usage: python3 test_join4.py [--out DIR]
60 random instances (class vectors of 1..8 columns with steps, offsets, constant columns, negative values; list sizes 5..150;
two thirds of the targets are the sum of one random vector of each list, i.e. a planted solution, the others such a sum plus noise in -2..2): the set of
solutions (a, b, c, d) printed by ./join4 (at 1, 4 and 16 key chunks) must equal the brute force set (all pairs
(A, B) against all pairs (C, D) by exact comparison through a numpy join on a random linear key, every hit re-checked on the full vector).
join4 is compiled into --out (default: current directory).
"""
import argparse, os, subprocess, sys
import numpy as np
import tbio

ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
ap.add_argument("--out", default=None, help="output directory (default: current directory)")
args = ap.parse_args()
out = tbio.outdir(args.out)
exe = tbio.build_join4(out)
fn = os.path.join(out, "test_join4.bin")

def brute(Fs, T):
    n = [len(f) for f in Fs]
    key = np.random.default_rng(1).integers(1, 2**40, size=Fs[0].shape[1])
    K = [f.astype(np.int64) @ key for f in Fs]
    L = np.add.outer(K[0], K[2]).ravel(); R = np.add.outer(K[1], K[3]).ravel()   # (A,C) against (B,D), as join4
    tk = int(np.asarray(T, dtype=np.int64) @ key)
    order = np.argsort(R, kind="stable"); Rs = R[order]
    lo = np.searchsorted(Rs, tk - L, "left"); hi = np.searchsorted(Rs, tk - L, "right")
    sols = set()
    for h in np.flatnonzero(hi > lo):
        a, c = divmod(int(h), n[2])
        for p in range(lo[h], hi[h]):
            b, d = divmod(int(order[p]), n[3])
            if np.array_equal(Fs[0][a] + Fs[1][b] + Fs[2][c] + Fs[3][d], T): sols.add((a, b, c, d))
    return sols

def run_join(Fs, T, lch):
    tbio.write_binary_lists(fn, Fs[0], Fs[1], Fs[2], Fs[3], np.asarray(T, dtype=np.int32))
    r = subprocess.run([exe, fn, str(lch)], capture_output=True, text=True, check=True)
    sols = {tuple(int(x) for x in l.split()[1:5]) for l in r.stdout.splitlines() if l.startswith("SOL")}
    done = [l for l in r.stdout.splitlines() if l.startswith("DONE")][0].split()
    assert int(done[2]) == len(sols), "solution count differs from the SOL lines"
    return sols

rng = np.random.default_rng(7)
ok = True; nsol_total = 0; nplanted = 0
for trial in range(60):
    while True:
        C = int(rng.integers(1, 9)); m = int(rng.choice([5, 30, 80, 150]))
        step = rng.integers(1, 5, size=C); off = rng.integers(-20, 20, size=C)
        spread = rng.integers(2, 12, size=C); spread[1:][rng.random(C - 1) < 0.2] = 0     # some constant columns (not all)
        sizes = [m + int(rng.integers(0, 30)) for _ in range(4)]
        est = np.prod([float(x) for x in sizes]) / np.prod(4.0 * spread + 1)
        if est < 2e5: break
    Fs = [(rng.integers(0, spread[None, :] + 1, size=(sizes[k], C)) * step[None, :] + off[None, :]).astype(np.int16) for k in range(4)]
    if trial % 3:
        T = sum(F[rng.integers(len(F))].astype(np.int64) for F in Fs); nplanted += 1
    else:
        T = sum(F[rng.integers(len(F))].astype(np.int64) for F in Fs) + rng.integers(-2, 3, size=C)
    b = brute(Fs, T)
    nsol_total += len(b)
    for lch in (0, 2, 4):
        e = run_join(Fs, T, lch)
        if e != b:
            ok = False; print("MISMATCH trial", trial, "C", C, "sizes", sizes, "lch", lch, len(e), len(b))
os.remove(fn)
print("join4 random tests: 60 instances (%d planted), %d solutions in total, 3 chunk settings each:" % (nplanted, nsol_total), "ALL OK" if ok else "FAILURES")
sys.exit(0 if ok else 1)
