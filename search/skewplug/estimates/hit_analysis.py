"""Anatomy of the hit that found skew/1892: re-runs the tenth round (index 9) of the seeded loop of prod_ks.py (K = Z43, q = 11, H_K = <6>, slots of data/cost473.jsonl, seed 1) and gives ALL level-1
instances of that round to level 2 (skf) with all completions: which instances complete, how many completions each (expected: 8 instances, 2 completions each = one t-solution up to t -> -t), and how the
completing instances share blocks (the 2 x 2 x 2 combinations of two A, two B, two C and one D: ONE solution cluster, not eight independent solutions).
usage (from search/skewplug):  python3 estimates/hit_analysis.py        (about 1 minute; writes runs/hit_analysis.json)"""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, '..')); sys.path.insert(0, HERE)
DATA = os.path.join(HERE, '..', 'data'); RUNS = os.path.join(HERE, '..', 'runs')
import json, math, random, time
import numpy as np
from sk import SShape
from sk_pat import row_patterns_skew
from sk_samp import level1
from gl1 import L1
from sl2io import header_text, run_sl2
from skf_io import skf_header, run_skf
sh = SShape([43], 11, [6]); l1 = L1(sh); hdr = header_text(sh)
slots = []
for l in open(os.path.join(DATA, 'cost473.jsonl')):
    r = json.loads(l)
    if 't_l2_inst' in r and r.get('n1_per_round', 0) > 0:
        c1 = r['t_l1_round'] / r['n1_per_round']; eff = 10 ** r['logrho'] / (c1 + r['t_l2_inst'])
        slots.append((math.log10(eff), tuple(r['al']), tuple(r['ms']), r.get('logN1', 0), r.get('logrho', 0)))
slots.sort(key=lambda x: -x[0]); slots = slots[:40]
mx = slots[0][0]; wts = [10 ** (s[0] - mx) for s in slots]
rnd = 9; rng = random.Random(1 * 1000003 + rnd)
(logE, al, ms, lN1, lrho) = rng.choices(slots, weights=wts, k=1)[0]
pats = row_patterns_skew(sh, al, ms); A_list = sh.enum_A(al[0], ms[0]); pat = rng.choice(pats)
print(al, ms, pat)
cnt, S, E, st = level1(sh, l1, al, ms, pat, 10000, 60000, 60000, seed=rng.randrange(10 ** 6), A_list=A_list, maxout=20000, hseed=rnd + 1)
print('level-1 instances', cnt)
res_all = {}; sols_all = []
t0 = time.time()
for c0 in range(0, len(S), 200):
    insts = [[(S[i, b], E[i, b]) for b in range(4)] for i in range(c0, min(len(S), c0 + 200))]
    res, ss = run_skf(sh, insts, maxsol=1000, tlimit=60.0, hdr=skf_header(sh))
    for i, d in res.items(): res_all[c0 + i] = d['compl']
    for (i, ts) in ss: sols_all.append((c0 + i, ts))
print('level 2 time %.0f s' % (time.time() - t0))
comp = {i: c for i, c in res_all.items() if c > 0}
print('instances with completions', len(comp), 'completions (x2 gauge)', sum(comp.values()), comp)
def key(i, blocks): return tuple(tuple(S[i, b].tolist()) + tuple(E[i, b].tolist()) for b in blocks)
for i in comp: print(i, 'A-key hash', hash(key(i, [0])) % 10**6, 'B', hash(key(i, [1])) % 10**6, 'C', hash(key(i, [2])) % 10**6, 'D', hash(key(i, [3])) % 10**6)
print('distinct solutions printed', len(sols_all))
os.makedirs(RUNS, exist_ok=True)
json.dump(dict(al=al, ms=ms, pat=pat, n1=int(cnt), comp=comp), open(os.path.join(RUNS, 'hit_analysis.json'), 'w'))
