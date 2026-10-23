"""Validation of prime-power q (q = 27 = 3^3, F_27 = GF(27), bf.py with GF tables): G-level brute force (bf.gcount) vs K-level brute force (bf.kcount) vs the two-level pipeline (sk_exh) on ten shapes
over K in {Z7 (with <2> and with H_K = 1), Z3, Z13 <3>, Z5, Z9 <4>, Z11 <3>, Z19 <7>, Z3^2 unipotent, Z7^2 with the order-21 group}.  Status MATCH when the counts agree.  Index 0 (Z7 <2> x GF(27)) is
nonzero: 6,144 ordered solutions; index 9 (Z7^2 x GF(27), order 21: v = 1323) is EXHAUSTIVELY empty (G-level brute force over G = 1323, about 6 minutes on one core).
usage (from search/skewplug):  python3 patterns/val_q27.py [i,j,...]      (default all; one JSON line per shape is appended to runs/val_q27.jsonl)"""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, '..')); sys.path.insert(0, HERE)
import json, time, warnings
warnings.filterwarnings('ignore')
from sk import SShape
from sk_exh import exhaust
from bf import gcount, kcount
SH = [([7], 27, [[2]]), ([7], 27, []), ([3], 27, []), ([13], 27, [[3]]), ([5], 27, []), ([9], 27, [[4]]), ([11], 27, [[3]]), ([19], 27, [[7]]),
      ([3, 3], 27, [[[1, 1], [0, 1]]]), ([7, 7], 27, [[[2, 0], [0, 2]], [[1, 1], [0, 1]]])]
os.makedirs(os.path.join(HERE, '..', 'runs'), exist_ok=True)
sel = [int(x) for x in sys.argv[1].split(',')] if len(sys.argv) > 1 else range(len(SH))
for (K, q, hu) in [SH[i] for i in sel]:
    sh = SShape(K, q, hu); t0 = time.time()
    rec = dict(K=K, q=q, hu=hu, G=sh.G, nuK=sh.kCD.nu)
    try: rec['gcount'] = gcount(K, q, hu, maxorb=21)[0]
    except Exception as ex: rec['gcount'] = None; rec['gerr'] = str(ex)[:60]
    try: rec['kcount'] = kcount(sh, maxu=7)[0]
    except Exception as ex: rec['kcount'] = None; rec['kerr'] = str(ex)[:60]
    try:
        tot, sols = exhaust(sh, tmax=100); rec['pipe'] = tot['E_G']; rec['exh'] = tot['exhaustive']
    except Exception as ex: rec['pipe'] = None; rec['perr'] = str(ex)[:60]
    vals = [v for v in (rec['gcount'], rec['kcount'], rec.get('pipe')) if v is not None]
    rec['status'] = 'MATCH' if len(set(vals)) == 1 and len(vals) >= 2 else ('MISMATCH' if len(set(vals)) > 1 else 'PARTIAL')
    rec['secs'] = round(time.time() - t0, 1)
    open(os.path.join(HERE, '..', 'runs', 'val_q27.jsonl'), 'a').write(json.dumps(rec) + '\n'); print(rec, flush=True)
