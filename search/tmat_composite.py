"""Cooper-Wallis routes with composite t: T-matrices over Z_t (t = n/w composite, w Williamson order).
T-matrices of order t <=> GS quadruples of length t with A*B*C*D = 1 pointwise.  For each multiplier structure
on Z_t (orbits <= MAXO, -1 not in H) we take the Gaussian GS estimate E_GS and the heuristic
E_T = E_GS * 2^-(k-1) (k orbits; the orbit {0} is forced by A(0)B(0)C(0)D(0) = 1 up to a global sign choice)."""
import os
import re, sys, math
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import spence_lib as L, gs_lib as GS
from gs_sweep_subgroups import subgroups
from sympy import isprime
WILL = list(range(1, 18, 2))
MAXO = int(sys.argv[1])
ns = [int(m.group(2)) for line in open(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'FRONTIER.md')) for m in [re.match(r'\|\s*(\d+)\s*\|\s*(\d+)\s*\|', line)] if m]
for n in ns:
    for w in WILL:
        if n % w: continue
        t = n // w
        if isprime(t) or t < 100: continue
        for H, gens in subgroups(t).items():
            if (t - 1) in H: continue
            G = L.AbGroup([t]); st = L.Structure(G, [G.perm_from_matrix(np.array([[g]])) for g in gens])
            if not (6 <= st.O <= MAXO): continue
            best = GS.gs_estimate(st, samples=1500)
            if not best or not math.isfinite(best['lE']): continue
            lET = best['lE'] - (st.O - 1) * math.log10(2)
            print(f"{lET:7.2f} n={n} w={w} t={t} |H|={len(H)} gens={gens} orbits={st.O} classes={st.C} pattern={best['pattern']} lE_GS={best['lE']:.2f}", flush=True)
