"""Gaussian estimates (est_skew.estimate: logN1, logE, logrho) of a random sample of the (allocation, m-split) slots of a shape, appended to a jsonl file (resumable; records with an error are kept with an 'err' field).
usage (from search/skewplug):  python3 estimates/est_slots.py 'K' q 'hunits' SEED NSLOTS OUT.jsonl
data/est473.jsonl is the 58 slots with a finite estimate of:  python3 estimates/est_slots.py "[43]" 11 "[[6]]" 1 400 runs/est473.jsonl    (the slot file of the loop prod_dedup.py)"""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, '..')); sys.path.insert(0, HERE)
DATA = os.path.join(HERE, '..', 'data'); RUNS = os.path.join(HERE, '..', 'runs')
import json, math, warnings, random, time, os
warnings.filterwarnings('ignore')
from sk import SShape
from sk_pat import allocations, m_splits
import est_skew as es
K = json.loads(sys.argv[1]); q = int(sys.argv[2]); hu = json.loads(sys.argv[3]); seed = int(sys.argv[4]); nslots = int(sys.argv[5]); out = sys.argv[6]
sh = SShape(K, q, hu)
pool = [(a, ms) for (a, em) in allocations(sh) for ms in m_splits(a, em)]
rng = random.Random(seed); rng.shuffle(pool)
done = set()
if os.path.exists(out):
    for l in open(out):
        r = json.loads(l); done.add((tuple(r['al']), tuple(r['ms'])))
print('pool', len(pool), 'done', len(done), flush=True)
cnt = 0
for (a, ms) in pool:
    if cnt >= nslots: break
    if (a, ms) in done: continue
    cnt += 1
    try: r = es.estimate(sh, a, ms, M=250)
    except Exception as ex: r = dict(err=str(ex))
    rec = dict(al=a, ms=ms, **{k: (round(v, 2) if isinstance(v, float) else v) for k, v in r.items()})
    open(out, 'a').write(json.dumps(rec) + '\n')
