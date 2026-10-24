"""apply tau_filter.tau_alive (the exact tau filter of search/NOTES.md section 92) to every live Spence shape of the census files (../data/census/cand*.jsonl);
usage: tau_sweep.py OUT.jsonl CAND1 [CAND2 ...]  (resumable; shapes in ascending order of the orbit count)"""
import sys, os, json, time, signal
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gsk import *
from tau_filter import tau_alive
class TO(Exception): pass
def _al(a, b): raise TO()
signal.signal(signal.SIGALRM, _al)
out = sys.argv[1]; files = sys.argv[2:]
done = set()
if os.path.exists(out):
    for l in open(out): done.add(json.loads(l)['skey'])
fo = open(out, 'a')
rows = []
for f in files: rows += [json.loads(l) for l in open(f)]
rows = [r for r in rows if r['status'] == 'live' and r['mode'] == 'b']
rows.sort(key=lambda r: (r['orb'], r['n']))
for r in rows:
    sk = json.dumps([r['n'], r['r'], r['K'], r['hg']])
    if sk in done: continue
    done.add(sk)
    hg = [(tuple(u), ph) for u, ph in r['hg']]
    rec = dict(skey=sk, n=r['n'], K=r['K'], r=r['r'], hg=r['hg'], orb=r['orb'])
    t0 = time.time()
    try:
        signal.alarm(150)
        try:
            sh = GShape(r['K'], r['r'], hg, 'spence', None)
            a, info = tau_alive(sh)
        finally:
            signal.alarm(0)
        rec['alive'] = a; rec['info'] = info if isinstance(info, str) else {k: int(v) for k, v in info.items()}
    except TO:
        rec['alive'] = None; rec['info'] = 'timeout'
    except Exception as ex:
        rec['alive'] = None; rec['info'] = 'error ' + repr(ex)[:100]
    rec['secs'] = round(time.time() - t0, 1)
    fo.write(json.dumps(rec) + '\n'); fo.flush()
print('done')
