#!/usr/bin/env python3
"""RRR search for GS quadruples over Z_v constant on orbits of H, in orbit coordinates.

y: (4, k) orbit values; X_i(f_j) = (y_i @ P)[j], P[O, j] = sum_{x in O} zeta^(f_j x).
P_B scales each frequency-orbit 4-vector to norm 2 sqrt(v); y' = Re(X' @ Pinv).
Euclidean metric of x-space (orbit weights |O|) is preserved by working with linear updates."""
import sys, time
import numpy as np

def orbits(v, gens):
    orb = -np.ones(v, dtype=np.int64); reps = []
    for t in range(v):
        if orb[t] < 0:
            S, todo = {t}, [t]
            while todo:
                y = todo.pop()
                for g in gens:
                    z = y * g % v
                    if z not in S: S.add(z); todo.append(z)
            orb[list(S)] = len(reps); reps.append(t)
    return orb, np.bincount(orb), np.array(reps)

class GS:
    def __init__(self, v, gens):
        self.v = v
        self.orb, self.sz, self.reps = orbits(v, gens)
        k = len(self.sz)
        x = np.arange(v)
        Z = np.exp(2j * np.pi * np.outer(x, self.reps) / v)        # (v, k): zeta^(f_j x)
        P = np.zeros((k, k), dtype=complex)
        np.add.at(P, self.orb, Z)
        self.P = P; self.Pinv = np.linalg.inv(P)
        self.k = k
    def PB(self, y):
        X = y @ self.P
        nrm = np.sqrt((np.abs(X) ** 2).sum(0))
        X = X * (2 * np.sqrt(self.v) / np.maximum(nrm, 1e-12))
        return (X @ self.Pinv).real
    def expand(self, y):
        return y[:, self.orb]
    def is_sol(self, ys):
        X = ys @ self.P
        return np.abs((np.abs(X) ** 2).sum(0) - 4 * self.v).max() < 1e-6

def run(g, beta=0.5, iters=10**6, seed=0, report=0, check=1):
    rng = np.random.default_rng(seed)
    y = rng.standard_normal((4, g.k))
    t = time.time()
    for it in range(iters):
        pb = g.PB(y)
        pa = np.where(2 * pb - y >= 0, 1.0, -1.0)
        y = y + beta * (pa - pb)
        if it % check == 0:
            s = np.where(pb >= 0, 1.0, -1.0)
            if g.is_sol(s): return it, s, time.time() - t
            if g.is_sol(pa): return it, pa, time.time() - t
        if report and it % report == 0:
            print(f"  it={it} |pa-pb|^2={(((pa - pb) ** 2) * g.sz).sum():.1f}", flush=True)
    return None, None, time.time() - t

if __name__ == '__main__':
    v = int(sys.argv[1]); gens = [int(x) for x in sys.argv[2].split(',')]
    beta = float(sys.argv[3]); iters = int(float(sys.argv[4])); seeds = int(sys.argv[5])
    rep = int(sys.argv[6]) if len(sys.argv) > 6 else 0
    g = GS(v, gens)
    print(f"v={v} gens={gens} orbits={g.k}", flush=True)
    for sd in range(seeds):
        it, s, dt = run(g, beta, iters, sd, rep)
        if it is not None:
            full = g.expand(s)
            tot = sum(np.array([full[i] @ np.roll(full[i], -t) for t in range(v)]) for i in range(4))
            ok = tot[0] == 4 * v and not np.any(tot[1:])
            print(f"seed={sd}: SOLVED it={it} rowsums={sorted(abs(int(r)) for r in full.sum(1))} PAFcheck={ok} ({dt:.1f}s)", flush=True)
            np.save(f"sol_{v}_{gens[0]}_{sd}.npy", full)
        else:
            print(f"seed={sd}: no solution ({dt:.1f}s)", flush=True)
