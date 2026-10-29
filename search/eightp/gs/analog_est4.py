"""analog_est4.py -- est4 against the EXACT two-level count (gs_exh.exhaust) on a small analogue Z_m x F_r, H_K = <h>, of the 5912 shape (Z211 x F7, <19>, 15 K-orbits; README, "Order 5912").
Appends one JSON record to out.jsonl: the exact K-level count (K_count, exact_logE; exh False = a lower bound), the number of level-1 instances N1, and est4's log10 E and log10 N1.
usage: analog_est4.py m r h tmax out.jsonl      e.g. analog_est4.py 31 7 2 600 ../data/gs/p5912/analog.jsonl"""
import sys, os, json, math, time
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from gsk import GShape
from gs_exh import exhaust
import k4, est4

m, r, h = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]); tmax = float(sys.argv[4]); out = sys.argv[5]
hg = [((h,), 1)]
rec = dict(m=m, r=r, h=h, G=m * r, korbits=len(set(frozenset(pow(h, j, m) * x % m for j in range(m)) for x in range(m))))
t0 = time.time()
tot, _ = exhaust(GShape([m], r, hg, 'spence', None), tmax=tmax, max_inst=3_000_000, tl2=30.0)
rec.update(exh=tot['exhaustive'], K_count=tot['compl'], E_G=tot['E_G'], N1=tot['n1'], skipped=len(tot['skipped']), secs=round(time.time() - t0, 1), nalloc=tot['allocs'])
try:
    e = est4.estimate_shape(k4.Shape4([m], r, hg, ([], [()])), M=600, nfull=200, topk=6, ntop=1)
    rec['est_logE'] = e.get('logE'); rec['est_logN1'] = e.get('logN1')
except Exception as ex:
    rec['est_err'] = repr(ex)[:100]
rec['exact_logE'] = math.log10(tot['compl']) if tot['compl'] > 0 else None
rec['exact_logN1'] = math.log10(tot['n1']) if tot['n1'] > 0 else None
print(json.dumps(rec), flush=True)
open(out, 'a').write(json.dumps(rec) + '\n')
