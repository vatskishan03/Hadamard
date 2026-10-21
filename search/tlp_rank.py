"""Rank the twisted Legendre pair structures of tlp_model.py for the open orders of FRONTIER.md.

    python3 search/tlp_rank.py MAXO

For every n in the second column of the table in FRONTIER.md put l = 2n - 1 (a Legendre pair of length l gives a
Hadamard matrix of order 2(l + 1) = 4n). Every structure (H, mu) of tlp_model.twisted_structures(l, MAXO) with 6 to
MAXO orbits of H on Z_l is scored by tlp_model.model (Gaussian estimate lE = log10 of the expected number of
solutions, with the pool size); one line is printed for each structure with lE >= -6.
"""
import os
import sys, re, math
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import tlp_model as TM
maxo = int(sys.argv[1])
ns = [int(m.group(2)) for line in open(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'FRONTIER.md')) for m in [re.match(r'\|\s*(\d+)\s*\|\s*(\d+)\s*\|', line)] if m]
rng = np.random.default_rng(2)
for n in ns:
    l = 2 * n - 1
    for H, gens, mu, lab, k in TM.twisted_structures(l, maxo):
        if k < 6: continue
        m = TM.model(l, H, mu, lab, k, samples=2000, rng=rng)
        if m is None or not math.isfinite(m['lE']) or m['lE'] < -6: continue
        print(f"{m['lE']:7.2f} n={n} l={l} |H|={len(H)} gens={gens} mu={mu} orbits={k} classes={m['C']} pool={m['N']:.3g}", flush=True)
