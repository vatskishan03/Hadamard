"""Measures, per slot of a shape, the level-1 yield and the level-2 cost per instance (engine sl2 or skf): for each slot of an estimate file (records with logE, best first) three level-1 rounds with a random
row-sum pattern and level 2 on up to 60 of the instances; writes n1_per_round, t_l1_round, t_l2_inst (seconds, machine dependent).  A file of such records is a 'cost file' of prod_ks.py / prod_dedup.py
(weight = estimated completions per second); data/cost473.jsonl is one for K = Z43, q = 11, H_K = <6>.
usage (from search/skewplug):  python3 estimates/slot_cost.py SLOTFILE OUT.jsonl 'K' q 'hunits' LB LC LD ENGINE [TOPN]
e.g.  python3 estimates/slot_cost.py data/est473.jsonl runs/cost473.jsonl "[43]" 11 "[[6]]" 10000 60000 60000 sl2 40"""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, '..')); sys.path.insert(0, HERE)
DATA = os.path.join(HERE, '..', 'data'); RUNS = os.path.join(HERE, '..', 'runs')
import json, time, random
import numpy as np
from sk import SShape
from sk_pat import row_patterns_skew
from sk_samp import level1
from gl1 import L1
from sl2io import header_text, run_sl2
from skf_io import skf_header, run_skf
slotfile, out = sys.argv[1], sys.argv[2]; K = json.loads(sys.argv[3]); q = int(sys.argv[4]); hu = json.loads(sys.argv[5]); LB, LC, LD = [int(float(x)) for x in sys.argv[6:9]]; engine = sys.argv[9]; topn = int(sys.argv[10]) if len(sys.argv) > 10 else 100000
sh = SShape(K, q, hu); l1 = L1(sh); h1 = header_text(sh) if engine == 'sl2' else skf_header(sh)
slots = []
for l in open(slotfile):
    r = json.loads(l)
    if 'logE' in r and r['logE'] > -50: slots.append(r)
slots.sort(key=lambda r: -r['logE']); slots = slots[:topn]
done = set()
if os.path.exists(out):
    for l in open(out): r = json.loads(l); done.add((tuple(r['al']), tuple(r['ms'])))
rng = random.Random(11)
for r in slots:
    al = tuple(r['al']); ms = tuple(r['ms'])
    if (al, ms) in done: continue
    try:
        pats = row_patterns_skew(sh, al, ms); A = sh.enum_A(al[0], ms[0]) if pats else None
        if not pats or A is None: rec = dict(al=al, ms=ms, logE=r['logE'], logrho=r.get('logrho'), logN1=r.get('logN1'), nopat=1)
        else:
            tt = 0; n1 = 0; S_all = []; E_all = []
            for k in range(3):
                pat = rng.choice(pats); t0 = time.time()
                c, S, E, st = level1(sh, l1, al, ms, pat, LB, LC, LD, seed=rng.randrange(10**6), A_list=A, maxout=3000, hseed=k + 1); tt += time.time() - t0; n1 += c
                if c: S_all.append(S[:20]); E_all.append(E[:20])
            rec = dict(al=al, ms=ms, logE=r['logE'], logrho=r.get('logrho'), logN1=r.get('logN1'), npat=len(pats), n1_per_round=n1 / 3, t_l1_round=tt / 3)
            if S_all:
                S = np.concatenate(S_all); E = np.concatenate(E_all); insts = [[(S[i, b], E[i, b]) for b in range(4)] for i in range(len(S))]
                t0 = time.time()
                if engine == 'sl2': res, _ = run_sl2(sh, insts, maxsol=1, tlimit=30.0, maxfib=1500000, stop_first=True, hdr=h1)
                else: res, _ = run_skf(sh, insts, maxsol=1, tlimit=30.0, hdr=h1)
                rec['t_l2_inst'] = (time.time() - t0) / len(insts); rec['n_l2'] = len(insts)
                rec['status'] = sorted(set(d['status'] for d in res.values()))
    except Exception as ex:
        rec = dict(al=al, ms=ms, err=str(ex))
    open(out, 'a').write(json.dumps(rec) + '\n')
