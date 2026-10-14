"""Rank slice-admissible Spence structures Z_2 x Z_m (n = 3 mod 4) by Gaussian E and exhaustive work."""
import os
import re, sys, math
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import spence_lib as L, spence_slice_check as SC
fr = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'FRONTIER.md')).read()
ns = [int(n) for _, n in re.findall(r'^\|\s*(\d+)\s*\|\s*(\d+)\s*\|', fr, re.M)]
MAXO = int(sys.argv[1])
for n in ns:
    if n % 4 != 3: continue
    m = (n - 1) // 2
    for H in SC.subgroups(m):
        sizes = SC.orbit_sizes(m, H)
        if 2 * len(sizes) > MAXO: continue
        ok, wit = SC.admissible(n, sizes)
        if not ok: continue
        G = L.AbGroup([2, m])
        st = L.Structure(G, [G.perm_from_matrix(np.array([[1, 0], [0, h]])) for h in H])
        est = st.estimate(samples=2000)
        nA, nB = st.pool(2), st.pool(0)
        print(f"n={n} |H|={len(H)} orbits={st.O} C={st.C} slice={wit} pools=({nA:.2e},{nB:.2e}) work={nA*nB+nB*nB:.2e} lE={est['lE'] if est else None}", flush=True)
