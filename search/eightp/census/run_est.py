"""Gaussian/lattice estimate (est4.estimate_shape, HEURISTIC) of every live shape of a census file (../data/census/cand*.jsonl): usage run_est.py in.jsonl out.jsonl M nfull topk  (resumable)"""
import sys, json, time, signal, math
import os
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, '..', 'gs'))
from est4 import *
class TO(Exception): pass
def _al(a, b): raise TO()
signal.signal(signal.SIGALRM, _al)
inp, out = sys.argv[1], sys.argv[2]
M = int(sys.argv[3]); nfull = int(sys.argv[4]); topk = int(sys.argv[5])
done = set()
try:
    for l in open(out):
        d = json.loads(l); done.add(d['key'])
except FileNotFoundError: pass
f = open(out, 'a')
rows = [json.loads(l) for l in open(inp)]
rows = [r for r in rows if r.get('status') == 'live']
t00 = time.time(); k = 0
for d in rows:
    key = json.dumps([d['mode'], d['n'], d['r'], d['K'], d['hg'], d['qspec']])
    if key in done: continue
    hg = [(tuple(u), ph) for u, ph in d['hg']]
    qs = d['qspec']
    if qs: qs = (list(qs[0]), [tuple(a) for a in qs[1]])
    t0 = time.time()
    res = dict(d); res['key'] = key
    try:
        sh = Shape4(d['K'], d['r'], hg, qs)
        signal.alarm(100)
        try:
            e = estimate_shape(sh, M=M, nfull=nfull, topk=topk, ntop=3)
        finally:
            signal.alarm(0)
        res.update(est=e)
        res['status2'] = 'ok'
    except TO:
        res['status2'] = 'timeout'
    except Exception as ex:
        res['status2'] = 'error:' + str(ex)[:100]
    res['secs'] = time.time() - t0
    f.write(json.dumps(res, default=lambda o: (float(o) if hasattr(o, 'dtype') else str(o))) + '\n'); f.flush(); k += 1
    if k % 50 == 0: print(k, len(rows), '%.0fs' % (time.time() - t00), flush=True)
print('done', k, '%.0fs' % (time.time() - t00))
