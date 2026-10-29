"""Validation of the fast level 2 (skf2 = pl2f2.c) against the exact level 2 of gl2.c (sl2) on level-1 instances: per instance the completion count, the fibre sizes and status, and the SETS of printed
solutions must be identical (all completions, no time limit).  Instances are sampled by the level 1 of sk_samp for one slot of a shape.
usage (from search/skewplug):  python3 val_skf.py 'K' q 'hunits' 'alloc' 'ms' LB LC LD NINST SEED
Last line: RESULT ... N instances, W with completions, C completions, count/status mismatches 0, solution-set mismatches 0, time sl2 .. skf .. (speedup ..x).
val_skf.sh runs the battery of 22 slots of seven shapes (7500 instances in the recorded run: 0 mismatches; skf2 1.8 to 6.7 times faster than sl2)."""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import json, time, random
import numpy as np
from sk import SShape
from sk_pat import allocations, m_splits, row_patterns_skew
from sk_samp import level1
from gl1 import L1
from sl2io import header_text, run_sl2
from skf_io import skf_header, run_skf
K = json.loads(sys.argv[1]); q = int(sys.argv[2]); hu = json.loads(sys.argv[3]); al = tuple(json.loads(sys.argv[4])); ms = tuple(json.loads(sys.argv[5]))
L = [int(float(x)) for x in sys.argv[6:9]]; ninst = int(sys.argv[9]); seed = int(sys.argv[10])
sh = SShape(K, q, hu); l1 = L1(sh); h1 = header_text(sh); h2 = skf_header(sh)
pats = row_patterns_skew(sh, al, ms); A = sh.enum_A(al[0], ms[0])
rng = random.Random(seed); S_all = []; E_all = []; tot = 0
for r in range(200):
    pat = rng.choice(pats) if pats else None
    cnt, S, E, st = level1(sh, l1, al, ms, pat, L[0], L[1], L[2], seed=rng.randrange(10**6), A_list=A, maxout=5000, hseed=r + 1)
    if cnt: S_all.append(S); E_all.append(E); tot += len(S)
    if tot >= ninst: break
if not S_all: print('NOINST'); sys.exit(0)
S = np.concatenate(S_all)[:ninst]; E = np.concatenate(E_all)[:ninst]
insts = [[(S[i, b], E[i, b]) for b in range(4)] for i in range(len(S))]
t0 = time.time(); r1, s1 = run_sl2(sh, insts, maxsol=100000, tlimit=1e9, maxfib=50_000_000, stop_first=False, hdr=h1); t1 = time.time() - t0
t0 = time.time(); r2, s2 = run_skf(sh, insts, maxsol=100000, tlimit=1e9, hdr=h2); t2 = time.time() - t0
bad = 0; ncompl = 0; nwith = 0
for i in range(len(insts)):
    a = r1.get(i); b = r2.get(i)
    if a is None or b is None: bad += 1; print('MISSING', i, a is None, b is None); continue
    if a['compl'] != b['compl'] or a['status'] != b['status']:
        bad += 1; print('MISMATCH', i, a['raw'], '|||', b['raw'])
    ncompl += a['compl']; nwith += a['compl'] > 0
def key(sols):
    d = {}
    for (i, ts) in sols: d.setdefault(i, set()).add(tuple(tuple(int(x) for x in t) for t in ts))
    return d
k1 = key(s1); k2 = key(s2)
setbad = sum(1 for i in range(len(insts)) if k1.get(i, set()) != k2.get(i, set()))
print(f'RESULT K={K} q={q} hu={hu} al={al} ms={ms}: {len(insts)} instances, {nwith} with completions, {ncompl} completions, count/status mismatches {bad}, solution-set mismatches {setbad}, time sl2 {t1:.1f}s skf {t2:.1f}s (speedup {t1/max(t2,1e-9):.2f}x)')
