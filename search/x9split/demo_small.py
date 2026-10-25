"""The small cases of NOTES.md §97: circulant splits of Williamson order 3 and of three Goethals-Seidel quadruples of
order 11 (each checked through the plug-in, H of order 108 and 396), and the infeasible propus splits for w = 5..11."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from csplit import solve, verify
from gen_small import gs_quads, propus_quads
W3 = [np.array([1, 1, 1]), np.array([1, -1, -1]), np.array([1, -1, -1]), np.array([1, -1, -1])]
for Q in [W3] + gs_quads(11, count=3):
    st, t, lab = solve('gs', Q, T=120)
    print(f"GS w={len(Q[0])}: {st} {t:.1f}s", "(parts ok, symmetric, H Hadamard, order) =", verify('gs', Q, lab) if lab else None, flush=True)
for w in (5, 7, 9, 11):
    for P in propus_quads(w, count=3):
        st, t, lab = solve('propus', P, T=120)
        print(f"propus w={w}: {st} {t:.1f}s", flush=True)
