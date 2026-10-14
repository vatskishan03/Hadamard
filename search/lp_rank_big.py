"""Gaussian E for plain Legendre-pair structures (A, B both H-invariant, row sums 1) with MINO..MAXO orbits."""
import os
import sys, re, math
from math import log10, pi
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lp_search as LP, tlp_model as TM
mino, maxo = int(sys.argv[1]), int(sys.argv[2])
ns = [int(a) for a in sys.argv[3:]] or [int(m.group(2)) for line in open(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'FRONTIER.md')) for m in [re.match(r'\|\s*(\d+)\s*\|\s*(\d+)\s*\|', line)] if m]
rng = np.random.default_rng(3)
for n in ns:
    l = 2 * n - 1
    for H, gens in LP.subgroups(l).items():
        lab, k = TM.orbit_labels(l, H)
        if not (mino <= k <= maxo): continue
        sizes = np.bincount(lab, minlength=k)
        cnt = {0: 1}
        for s in sizes:
            nxt = {}
            for t, c in cnt.items():
                nxt[t + s] = nxt.get(t + s, 0) + c; nxt[t - s] = nxt.get(t - s, 0) + c
            cnt = nxt
        N = cnt.get(1, 0)
        if N == 0: continue
        L = set(H) | {(-h) % l for h in H}
        cls = -np.ones(l, dtype=np.int64); reps = []
        for s in range(1, l):
            if cls[s] < 0:
                for x in L: cls[s * x % l] = len(reps)
                reps.append(s)
        Y = []
        while sum(len(y) for y in Y) < 3000:
            Z = np.where(rng.random((20000, k)) < 0.5, 1, -1); Y.append(Z[(Z @ sizes) == 1])
        Y = np.concatenate(Y)[:3000]
        A = Y[:, lab]; x = np.arange(l)
        V = np.stack([(A * A[:, (x + s) % l]).sum(1) for s in reps], 1).astype(float)
        terms = []; dead = False
        for c in range(len(reps)):
            u = np.unique(V[:, c]); g = int(np.gcd.reduce(np.diff(u).astype(np.int64))) if len(u) > 1 else 0
            m, var = V[:, c].mean(), 2 * V[:, c].var()          # sum of two independent sequences
            if g and (int(-2 - 2 * u[0])) % g: dead = True; break
            terms.append((math.sqrt(var), -0.5 * log10(2 * pi * var) - (2 * m + 2) ** 2 / (2 * var) / math.log(10) + log10(max(g, 1))))
        if dead: continue
        terms.sort(); lE = 2 * log10(N) + sum(t for _, t in terms[:-1])
        print(f"{lE:7.2f} n={n} l={l} |H|={len(H)} gens={gens} orbits={k} classes={len(reps)} pool={N:.3g}", flush=True)
