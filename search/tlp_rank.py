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
