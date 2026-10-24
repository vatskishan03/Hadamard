"""exh_sweep.py -- exhaustive two-level verdicts (gs_exh.exhaust, validated against G-level brute force) for the live Spence shapes of the census of
n = 2p (../data/census/cand8p.jsonl; the shapes of the verdicts of search/NOTES.md section 92 were also taken from ../data/gs/cand_alive.jsonl with EST=../data/census/est_all8p.jsonl).  usage: exh_sweep.py OUT.jsonl OMIN OMAX TMAX [N]   (resumable: shapes already in OUT are skipped)"""
import sys, os, json, time
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from gsk import *
from gs_exh import exhaust
out = sys.argv[1]; omin = int(sys.argv[2]); omax = int(sys.argv[3]); tmax = float(sys.argv[4])
nsel = int(sys.argv[5]) if len(sys.argv) > 5 else None
CAND = os.environ.get('CAND', os.path.join(HERE, '..', 'data', 'census', 'cand8p.jsonl')); EST = os.environ.get('EST', os.path.join(HERE, '..', 'data', 'census', 'est8p.jsonl'))
rows = [json.loads(l) for l in open(CAND)]
est = {}
for l in open(EST):
    d = json.loads(l); est[d['key']] = d
live = [r for r in rows if r['status'] == 'live' and omin <= r['orb'] <= omax and (nsel is None or r['n'] == nsel)]
def skey(r): return json.dumps([r['mode'], r['n'], r['r'], r['K'], r['hg']])
done = set()
if os.path.exists(out):
    for l in open(out): done.add(json.loads(l)['skey'])
fo = open(out, 'a')
def lE(r):
    k = json.dumps([r['mode'], r['n'], r['r'], r['K'], r['hg'], r['qspec']])
    e = est.get(k, {}).get('est', {})
    return e.get('logE', -99)
live.sort(key=lambda r: (r['orb'], -lE(r)))
for r in live:
    sk = skey(r)
    if sk in done: continue
    hg = [(tuple(u), ph) for u, ph in r['hg']]
    rec = dict(skey=sk, n=r['n'], K=r['K'], r=r['r'], hg=r['hg'], orb=r['orb'], estlogE=lE(r))
    t0 = time.time()
    try:
        sh = GShape(r['K'], r['r'], hg, 'spence', None)
        rec['nu'] = sh.nu; rec['feasible_E'] = sh.feasible_totals()[:12]
        tot, sols = exhaust(sh, tmax=tmax, max_inst=3_000_000, tl2=60.0)
        rec.update(E_G=tot['E_G'], n1=tot['n1'], exhaustive=tot['exhaustive'], allocs=tot['allocs'], skipped=[(list(a), w) for a, w in tot['skipped']][:6], flags=tot['flags'][:3], nsol=len(sols))
        if sols:
            import pickle
            pickle.dump([[(a.tolist(), b.tolist(), c.tolist()) for (a, b, c) in s] for s in sols], open('found_%d_%d.pkl' % (r['n'], abs(hash(sk)) % 10**6), 'wb'))
    except Exception as ex:
        rec['error'] = repr(ex)[:200]
    rec['secs'] = round(time.time() - t0, 1)
    fo.write(json.dumps(rec) + '\n'); fo.flush()
    print(rec['n'], rec['K'], rec['r'], rec['hg'], 'orb', rec['orb'], 'E_G', rec.get('E_G'), 'exh', rec.get('exhaustive'), rec.get('error', ''), rec['secs'], flush=True)
