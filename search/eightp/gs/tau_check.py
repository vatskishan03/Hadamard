"""tau_check.py -- re-derives the verdicts of the exact tau filter (tau_filter.tau_alive, search/NOTES.md section 92) for the shapes of the recorded sweep ../data/gs/tau_sweep2.jsonl with at most OMAX orbits and
compares them with the recorded ones (alive / dead / undecided, number of block states and of distinct feature vectors).  The whole sweep took 588 CPU seconds (the 18 and 20 orbit shapes 440 of them).
usage: tau_check.py [OMAX [OMIN]]      (default 12, 3)       prints 'MISMATCH ...' lines and a summary"""
import sys, os, json, time
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from gsk import *
from tau_filter import tau_alive
omax = int(sys.argv[1]) if len(sys.argv) > 1 else 12; omin = int(sys.argv[2]) if len(sys.argv) > 2 else 3
rows = [json.loads(l) for l in open(os.path.join(HERE, '..', 'data', 'gs', 'tau_sweep2.jsonl'))]
n = bad = 0; alive = dead = 0; t0 = time.time()
for r in rows:
    if not (omin <= r['orb'] <= omax) or r['alive'] is None: continue
    hg = [(tuple(u), ph) for u, ph in r['hg']]
    sh = GShape(r['K'], r['r'], hg, 'spence', None)
    a, info = tau_alive(sh)
    n += 1; alive += (a is True); dead += (a is False)
    ref = r['info']
    if a != r['alive'] or (isinstance(info, dict) and (info['blockstates'] != ref['blockstates'] or info['features'] != ref['features'])):
        bad += 1; print('MISMATCH', r['n'], r['K'], r['r'], r['hg'], 'now', a, info, 'recorded', r['alive'], ref, flush=True)
print('tau_check: %d shapes with %d..%d orbits re-derived (%d alive, %d dead), %d mismatches with the record  [%.0f s]' % (n, omin, omax, alive, dead, bad, time.time() - t0))
