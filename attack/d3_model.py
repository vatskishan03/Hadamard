"""Gaussian model for the d=3 twisted structure (X0, X1, mu X1, mu^2 X1) over Z_n."""
import sys, math
from math import log10, pi
import numpy as np
sys.path.insert(0, '/home/user/Hadamard/attack')
import spence_lib as L

def setup(n, h_gen, mu):
    G = L.AbGroup([n])
    stH = L.Structure(G, [G.perm_from_matrix(np.array([[h_gen]]))])
    stK = L.Structure(G, [G.perm_from_matrix(np.array([[h_gen]])), G.perm_from_matrix(np.array([[mu]]))])
    # K-classes reps; for each, the H-class index of s, mu s, mu^2 s
    fold = np.array([[stH.cls[(pow(mu, j, n) * s) % n] for j in range(3)] for s in stK.reps])
    return G, stH, stK, fold

def values(stH, stK, fold, Y0, Y1):
    V0 = stK.class_values(Y0)                       # (b, CK)
    V1h = stH.class_values(Y1)                      # (b, CH)
    V1 = V1h[:, fold].sum(2)                        # (b, CK)
    return V0, V1

if __name__ == '__main__':
    n, hg, mu = map(int, sys.argv[1:4])
    G, stH, stK, fold = setup(n, hg, mu)
    print(f"n={n} H-orbits={stH.O} K-orbits={stK.O} K-classes={stK.C} H-classes={stH.C}")
    rng = np.random.default_rng(1)
    for r0 in range(1, int(2 * math.sqrt(n)) + 2, 2):
        rem = 4 * n - r0 * r0
        if rem % 3: continue
        r1 = math.isqrt(rem // 3)
        if 3 * r1 * r1 != rem or r1 % 2 == 0: continue
        N0, N1 = stK.pool(r0), stH.pool(r1)
        print(f"pattern r0={r0} r1={r1}: pools N0={N0:.3g} N1={N1:.3g}")
        if N0 == 0 or N1 == 0: continue
        Y0 = stK.sample(r0, 3000, rng); Y1 = stH.sample(r1, 3000, rng)
        V0, V1 = values(stH, stK, fold, Y0, Y1)
        V = V0 + V1
        lp = []; dead = False
        for c in range(stK.C):
            u0 = np.unique(V0[:, c]); u1 = np.unique(V1[:, c])
            g = int(np.gcd.reduce(np.concatenate([np.diff(u0), np.diff(u1)]).astype(np.int64)))
            res = (int(u0[0]) + int(u1[0])) % g if g else int(u0[0] + u1[0])
            if g and res % g != 0: dead = True; break
            var = V[:, c].var(); mu_c = V[:, c].mean()
            lp.append((stK.csize[c] * math.sqrt(var), -0.5 * log10(2 * pi * var) - mu_c ** 2 / (2 * var) / math.log(10) + log10(max(g, 1))))
        if dead: print("   lattice excludes 0"); continue
        lp.sort(); lP = sum(t for _, t in lp[:-1])
        print(f"   log10 E = {log10(N0) + log10(N1) + lP:.2f}  (log10 P = {lP:.2f}); sd X0 {np.sqrt(V0.var(0)).mean():.1f} sd X1 {np.sqrt(V1.var(0)).mean():.1f}")
