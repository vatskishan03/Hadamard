"""join4 against a numpy brute force on random integer lists (with lattice structure, negative values, constant columns, planted solutions)"""
import sys, os, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
import numpy as np
import hlib

def brute(Fs, T):
    n = [len(f) for f in Fs]
    key = np.random.default_rng(1).integers(1, 2**40, size=Fs[0].shape[1])
    K = [f.astype(np.int64) @ key for f in Fs]
    L = np.add.outer(K[0], K[1]).ravel(); R = np.add.outer(K[2], K[3]).ravel()
    tk = int(np.asarray(T, dtype=np.int64) @ key)
    order = np.argsort(R, kind='stable'); Rs = R[order]
    lo = np.searchsorted(Rs, tk - L, 'left'); hi = np.searchsorted(Rs, tk - L, 'right')
    sols = set()
    for h in np.flatnonzero(hi > lo):
        a, b = divmod(int(h), n[1])
        for p in range(lo[h], hi[h]):
            c, d = divmod(int(order[p]), n[3])
            if np.array_equal(Fs[0][a] + Fs[1][b] + Fs[2][c] + Fs[3][d], T): sols.add((a, b, c, d))
    return sols

rng = np.random.default_rng(7)
ok = True
for trial in range(60):
    while True:
        C = int(rng.integers(1, 9)); m = int(rng.choice([5, 30, 120, 400]))
        step = rng.integers(1, 5, size=C); off = rng.integers(-20, 20, size=C)
        spread = rng.integers(2, 12, size=C); spread[1:][rng.random(C - 1) < 0.2] = 0     # some constant columns (not all)
        sizes = [m + int(rng.integers(0, 30)) for _ in range(4)]
        est = np.prod([float(x) for x in sizes]) / np.prod(4.0 * spread + 1)
        if est < 3e5: break
    Fs = [(rng.integers(0, spread[None, :] + 1, size=(sizes[k], C)) * step[None, :] + off[None, :]).astype(np.int32) for k in range(4)]
    # target: a planted solution (or a random one)
    if trial % 3:
        T = sum(F[rng.integers(len(F))].astype(np.int64) for F in Fs)
    else:
        T = sum(F[rng.integers(len(F))].astype(np.int64) for F in Fs) + rng.integers(-2, 3, size=C)
    b = brute(Fs, T)
    for budget in (50, 100000):
        tot, out, stats = hlib.Problem.join4(Fs, T, budget=budget)
        e = set(map(tuple, out.tolist()))
        if not (tot == len(b) and e == b):
            ok = False; print('MISMATCH trial', trial, 'C', C, 'm', m, 'budget', budget, tot, len(b))
print('join4 random tests:', 'ALL OK' if ok else 'FAILURES')
