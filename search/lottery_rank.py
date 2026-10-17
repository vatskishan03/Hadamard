"""Rank small Goethals-Seidel multiplier structures of all open orders by E (Gaussian) and exact pair work."""
import os
import re, sys, math
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import spence_lib as L, gs_lib as GS
from gs_sweep_subgroups import subgroups

fr = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'FRONTIER.md')).read()
ns = [int(n) for _, n in re.findall(r'^\|\s*(\d+)\s*\|\s*(\d+)\s*\|', fr, re.M)]
MIN_O, MAX_O = int(sys.argv[1]), int(sys.argv[2])
rows = []
for n in ns:
    for H, gens in subgroups(n).items():
        if (n - 1) in H: continue
        seen = set(); O = 0
        for t in range(n):
            if t not in seen:
                seen |= {t * h % n for h in H}; O += 1
        if not (MIN_O <= O <= MAX_O): continue
        G = L.AbGroup([n])
        st = L.Structure(G, [G.perm_from_matrix(np.array([[g]])) for g in gens])
        pats = [q for q in GS.patterns(n) if all(st.pool(r) > 0 for r in q)]
        if not pats: continue
        work = 0.0
        for q in pats:
            s = sorted(st.pool(r) for r in q)
            work += s[3] * s[0] + s[2] * s[1]
        best = GS.gs_estimate(st, samples=2000)
        lE = best['lE'] if best else -math.inf
        rows.append((lE, n, len(H), gens, O, st.C, work, best['pattern'] if best else None))
        print(f"{lE:7.2f}  n={n:5d} |H|={len(H):4d} gens={gens} O={O:3d} C={st.C:3d} work={work:.2e} {best['pattern'] if best else ''}", flush=True)
