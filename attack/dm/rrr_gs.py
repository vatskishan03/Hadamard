#!/usr/bin/env python3
"""Relaxed-reflect-reflect (Elser) search for Goethals-Seidel quadruples of length v.

Constraint A: entries +-1.
Constraint B: for every frequency k, sum_i |X_i(k)|^2 = 4v   (X_i = DFT of row i).
A intersect B = GS quadruples (sum of periodic autocorrelations = 4v delta).
RRR:  x <- x + beta * (P_A(2 P_B(x) - x) - P_B(x))."""
import sys, time
import numpy as np

def PB(x, v):
    X = np.fft.rfft(x, axis=1)                    # (4, v//2+1)
    nrm = np.sqrt((np.abs(X) ** 2).sum(0))        # per frequency
    scale = np.where(nrm > 1e-12, 2 * np.sqrt(v) / np.maximum(nrm, 1e-12), 1.0)
    return np.fft.irfft(X * scale, n=v, axis=1)

def PA(x):
    return np.where(x >= 0, 1.0, -1.0)

def check(s, v):
    X = np.fft.rfft(s, axis=1)
    return np.allclose((np.abs(X) ** 2).sum(0), 4 * v, atol=1e-6)

def run(v, beta=0.5, iters=10**6, seed=0, report=0):
    rng = np.random.default_rng(seed)
    x = rng.standard_normal((4, v))
    t = time.time()
    for it in range(iters):
        pb = PB(x, v)
        pa = PA(2 * pb - x)
        x = x + beta * (pa - pb)
        if it % 10 == 0:
            s = PA(pb)
            if check(s, v): return it, s, time.time() - t
            if check(pa, v): return it, pa, time.time() - t
        if report and it % report == 0:
            err = np.abs(pa - pb).sum()
            print(f"  it={it} err={err:.2f}", flush=True)
    return None, None, time.time() - t

if __name__ == '__main__':
    v = int(sys.argv[1]); beta = float(sys.argv[2]) if len(sys.argv) > 2 else 0.5
    iters = int(float(sys.argv[3])) if len(sys.argv) > 3 else 10**5
    seeds = int(sys.argv[4]) if len(sys.argv) > 4 else 3
    for sd in range(seeds):
        it, s, dt = run(v, beta, iters, sd)
        print(f"v={v} beta={beta} seed={sd}: {'SOLVED at it=%d' % it if it is not None else 'no solution'} ({dt:.1f}s)", flush=True)
