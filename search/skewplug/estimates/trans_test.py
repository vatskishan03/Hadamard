"""Translation invariance of level 2: for composite K the points of K fixed by H_K (Fix) act on level-1 instances by translating B, C or D separately (the K-level system is a sum of per-block
autocorrelation vectors).  For instances of the slot (0,3,12,18)/(0,1,12,18) of K = Z35, q = 7, H_K = <11> (Fix of size 5) every instance is translated by random points of Fix in B, C, D; the fibre
sizes, status and completion count of both level-2 engines (skf and sl2) must be identical among the translates ('with differing verdict/fibres among translates: 0').
usage (from search/skewplug):  python3 estimates/trans_test.py"""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, '..')); sys.path.insert(0, HERE)
DATA = os.path.join(HERE, '..', 'data'); RUNS = os.path.join(HERE, '..', 'runs')
import json, random, time, warnings, itertools
warnings.filterwarnings('ignore')
import numpy as np
from sk import SShape
from sk_samp import level1
from gl1 import L1
from skf_io import skf_header, run_skf
from sl2io import header_text, run_sl2
sh = SShape([35], 7, [[11]]); l1 = L1(sh); N = sh.N
al = (0, 3, 12, 18); ms = (0, 1, 12, 18); pat = ((20, -1), (15, -6), (5, -2))
cnt, S, E, st = level1(sh, l1, al, ms, pat, 2000, 6000, 6000, seed=3, A_list=sh.enum_A(al[0], ms[0]), maxout=300, hseed=1)
print('instances', cnt, len(S))
fix = [int(P[0]) for P in sh.kCD.pos if len(P) == 1]; print('fix', fix)
rng = random.Random(1)
insts = []; tags = []
for i in range(0, min(len(S), 200), 10):
    base = [(S[i, b], E[i, b]) for b in range(4)]
    for tb, tc, td in [(0, 0, 0)] + [(rng.choice(fix), rng.choice(fix), rng.choice(fix)) for _ in range(3)]:
        out = [base[0]]
        for X, t in zip((1, 2, 3), (tb, tc, td)):
            idx = np.asarray(sh.shift[t]); out.append((base[X][0][idx], base[X][1][idx]))
        insts.append(out); tags.append((i, (tb, tc, td)))
r, ss = run_skf(sh, insts, maxsol=1, tlimit=60.0, hdr=skf_header(sh))
r2, ss2 = run_sl2(sh, insts, maxsol=1, tlimit=60.0, maxfib=1500000, stop_first=False, hdr=header_text(sh))
by = {}
for k, (i, t) in enumerate(tags): by.setdefault(i, []).append((r[k]['compl'], r[k].get('fib'), r[k]['status'], r2[k]['compl'], r2[k]['status']))
bad = 0
for i, L in by.items():
    if len(set(L)) != 1: bad += 1; print('DIFF', i, L)
print('base instances', len(by), 'with differing verdict/fibres among translates:', bad, 'example', list(by.values())[0][:2])
