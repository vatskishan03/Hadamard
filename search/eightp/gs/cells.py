"""cells.py -- size of the exact search of a tau-alive real plug-in shape: for every allocation token (n.m x 4) and every exact row pattern ((sigma, eps) x 4) the pattern-restricted list sizes
(enum_pattern count_only) and the cost proxy of the class-restricted join  max(L0*L2, L1*L3) (blocks sorted by size).   usage: cells.py Kmods r hg [maxtok]"""
import sys, os, json, math, collections, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gsk import *
from prod_tau import tokens_for
Kmods = json.loads(sys.argv[1]); r = int(sys.argv[2]); hg = [(tuple(u), ph) for u, ph in json.loads(sys.argv[3])]
maxtok = int(sys.argv[4]) if len(sys.argv) > 4 else None
sh = GShape(Kmods, r, hg, 'spence', None)
toks = tokens_for(sh)
print('shape', Kmods, r, hg, 'orbits', sh.nu, 'tokens', len(toks), flush=True)
if maxtok: toks = toks[:maxtok]
cells = 0; tot_cost = 0.0; maxcost = 0.0; tot_inst = 0.0; costs = []; maxL = 0
t0 = time.time()
cache = {}
def cnt(n, m, sg, ep):
    k = (n, m, sg, ep)
    if k not in cache:
        cans = (sg == 0 and ep == 0)
        c = enum_pattern(sh, n, m, sg, ep, cans, count_only=True)
        cache[k] = c
    return cache[k]
for tk in toks:
    al = tuple(n for n, m in tk); ms = tuple(m for n, m in tk)
    P = row_patterns(sh, al, ms)
    if P is None: print('pattern enumeration too big for', tk); continue
    for pat in P:
        Ls = sorted(cnt(al[b], ms[b], pat[b][0], pat[b][1]) for b in range(4))
        cost = max(float(Ls[0]) * Ls[2], float(Ls[1]) * Ls[3]) if Ls[3] > 0 else 0.0
        cells += 1; tot_cost += cost; maxcost = max(maxcost, cost); costs.append(cost); maxL = max(maxL, Ls[3])
costs = np.array(costs) if costs else np.zeros(1)
print(f'cells (token, pattern) {cells}   sum join cost {tot_cost:.3g}   max cell {maxcost:.3g}   max list {maxL:.3g}   cells with cost > 1e9: {(costs > 1e9).sum()}  > 1e11: {(costs > 1e11).sum()}   [{time.time()-t0:.0f}s]')
