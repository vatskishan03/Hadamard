"""Validation of the real plug-in engine (gs_exh) against the independent G-level brute force (gbf.count_G) for SPENCE shapes with ODD |G| (never covered by validate.py,
which skipped them): spence mode, phi = -4, row sums (+-1)^4."""
import sys, os, json, time
import numpy as np
from gsk import *
from gs_exh import exhaust
from gbf import count_G
from pk import signed_subgroups
out = sys.argv[1]; tmax = float(sys.argv[2]) if len(sys.argv) > 2 else 60
fo = open(out, 'a')
Ks = [[3], [5], [7], [9], [3, 3], [11], [13], [15], [5, 3]]
jobs = []
for r in (3, 5, 7, 9):
    for K in Ks:
        try: sub = signed_subgroups(K, maxgens=2)
        except Exception: sub = [([(1,) * len(K)], [])]
        for (G_, gens) in sub:
            hg = [(tuple(u), int(s)) for (u, s) in gens]
            jobs.append((int(np.prod(K)) * r, K, r, hg))
jobs.sort(key=lambda j: (j[0], len(j[3])))
nm = nnz = 0
for (Gsz, K, r, hg) in jobs:
    if Gsz % 2 == 0: continue
    rec = dict(K=K, r=r, hg=[[list(u), p] for u, p in hg], Gsz=Gsz)
    try: sh = GShape(K, r, hg, 'spence', None)
    except AssertionError as ex: continue
    rec['nu'] = sh.nu
    if sh.nu > 7: continue
    t0 = time.time()
    try:
        gc, info = count_G(K, r, hg, 'spence', None)
    except Exception as ex:
        continue
    try:
        tot, sols = exhaust(sh, tmax=tmax, max_inst=1_500_000, tl2=30.0)
    except Exception as ex:
        rec['status'] = 'error ' + repr(ex)[:100]; fo.write(json.dumps(rec) + '\n'); fo.flush(); continue
    rec.update(G_count=gc, K_count=tot['E_G'], exhaustive=tot['exhaustive'])
    rec['status'] = ('MATCH' if tot['E_G'] == gc else 'MISMATCH') if tot['exhaustive'] else 'partial'
    nm += rec['status'] == 'MATCH'; nnz += (rec['status'] == 'MATCH' and gc > 0)
    fo.write(json.dumps(rec) + '\n'); fo.flush()
    print(rec['status'], K, r, hg, gc, tot['E_G'], '%.1fs' % (time.time() - t0), flush=True)
print('matches', nm, 'nonzero', nnz)
