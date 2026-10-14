#!/usr/bin/env python3
"""RRR (Elser) search for Goethals-Seidel quadruples over Z_v constant on the orbits of H = <g>.

A: entries +-1;  B: sum_i |X_i(k)|^2 = 4v at every frequency k.  Both projections preserve
H-invariance (X(hk) = X(k) for H-invariant x), so iterates stay in the symmetric subspace;
we re-average over orbits each step to remove round-off drift."""
import sys, time
import numpy as np

def orbits(v, gens):
    orb = -np.ones(v, dtype=np.int64); k = 0
    for t in range(v):
        if orb[t] < 0:
            S, todo = {t}, [t]
            while todo:
                y = todo.pop()
                for g in gens:
                    z = y * g % v
                    if z not in S: S.add(z); todo.append(z)
            orb[list(S)] = k; k += 1
    return orb, np.bincount(orb)

def make(v, gens):
    orb, sz = orbits(v, gens)
    def sym(x):
        m = np.stack([np.bincount(orb, weights=r) for r in x]) / sz
        return m[:, orb]
    def PB(x):
        X = np.fft.rfft(x, axis=1)
        nrm = np.sqrt((np.abs(X) ** 2).sum(0))
        sc = 2 * np.sqrt(v) / np.maximum(nrm, 1e-12)
        return np.fft.irfft(X * sc, n=v, axis=1)
    return orb, sz, sym, PB

def is_sol(s, v):
    X = np.fft.rfft(s, axis=1)
    return np.abs((np.abs(X) ** 2).sum(0) - 4 * v).max() < 1e-6

def run(v, gens, beta=0.5, iters=10**6, seed=0, report=0):
    orb, sz, sym, PB = make(v, gens)
    rng = np.random.default_rng(seed)
    x = sym(rng.standard_normal((4, v)))
    t = time.time()
    for it in range(iters):
        pb = PB(x)
        r = 2 * pb - x
        pa = np.where(r >= 0, 1.0, -1.0)
        x = sym(x + beta * (pa - pb))
        if it % 5 == 0:
            s = np.where(pb >= 0, 1.0, -1.0)
            if is_sol(s, v): return it, s, time.time() - t
            if is_sol(pa, v): return it, pa, time.time() - t
        if report and it % report == 0:
            print(f"  it={it} |pa-pb|^2={((pa - pb) ** 2).sum():.1f}", flush=True)
    return None, None, time.time() - t

if __name__ == '__main__':
    v = int(sys.argv[1]); gens = [int(g) for g in sys.argv[2].split(',')]
    beta = float(sys.argv[3]); iters = int(float(sys.argv[4])); seeds = int(sys.argv[5])
    rep = int(sys.argv[6]) if len(sys.argv) > 6 else 0
    orb, sz = orbits(v, gens)
    print(f"v={v} gens={gens} orbits={len(sz)}", flush=True)
    for sd in range(seeds):
        it, s, dt = run(v, gens, beta, iters, sd, rep)
        msg = f"SOLVED it={it} rowsums={sorted(abs(int(r)) for r in s.sum(1))}" if it is not None else "no solution"
        print(f"seed={sd}: {msg} ({dt:.1f}s)", flush=True)
        if it is not None:
            np.save(f"sol_{v}_{sd}.npy", s)
