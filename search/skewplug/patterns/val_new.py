"""Validation of the orbit structures with matrix-valued H_K and non-cyclic K: G-level brute force (bf.gcount) vs K-level brute force (bf.kcount) vs the two-level pipeline (sk_exh): exact counts of
ordered solutions on 38 small shapes; status MATCH when all available counts agree.  Shape indices: 0-4 Z3^2 with the unipotent matrix of order 3 (q = 3, 7, 11, 19, 23); 5-7 and 18-22 Z3^3 with
a Singer matrix of order 13 (q = 3, 7, 11 and 19, 23, 31, 43, 47); 8-15 K = Z3 x Z9, Z3 x Z15, Z3 x Z21 with unit tuples; 16-17 Z5^2 unipotent (q = 3, 7); 23-28 Z7^2 with the order-21 group <2 I, unipotent>
(q = 3, 7, 11, 19, 23, 27; index 25 is v = 539); 29-33 Z5^3 with a Singer matrix of order 31 (q = 3 ... 23); 34-37 Z7^3 with the Singer matrix of order 171 (q = 3, 7, 11, 19; index 34 is v = 1029, 35 is v = 2401).
Only Z3^3 x F7 (index 6, 6,144 ordered solutions) is nonzero among the non-cyclic matrix shapes; the zero counts back the EXHAUSTIVELY empty verdicts of 539, 1029 and 2401.
usage (from search/skewplug):  python3 patterns/val_new.py [i,j,...]      (default all; one JSON line per shape is appended to runs/val_new.jsonl)"""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, '..')); sys.path.insert(0, HERE)
REPO = os.path.abspath(os.path.join(HERE, '..', '..', '..'))
import json, time, warnings
warnings.filterwarnings('ignore')
import numpy as np
from sk import SShape
from sk_exh import exhaust
from bf import gcount, kcount
U3 = [[1, 1], [0, 1]]
S13 = [[0, 0, 1], [0, 1, 1], [2, 1, 1]]      # an element of order 13 of GL_3(F_3) (Singer cycle)
SH = [([3, 3], 3, [U3]), ([3, 3], 7, [U3]), ([3, 3], 11, [U3]), ([3, 3], 19, [U3]), ([3, 3], 23, [U3]),
      ([3, 3, 3], 3, [S13]), ([3, 3, 3], 7, [S13]), ([3, 3, 3], 11, [S13]),
      ([3, 9], 3, [[1, 4]]), ([3, 9], 7, [[1, 4]]), ([3, 9], 11, [[1, 4]]), ([3, 9], 3, [[1, 7]]),
      ([3, 15], 3, [[1, 4]]), ([3, 15], 7, [[1, 4]]), ([3, 21], 3, [[1, 4]]), ([3, 21], 7, [[1, 2]]),
      ([5, 5], 3, [[[1, 1], [0, 1]]]), ([5, 5], 7, [[[1, 1], [0, 1]]])]
SC = [[2, 0], [0, 2]]
SH += [([3, 3, 3], 19, [S13]), ([3, 3, 3], 23, [S13]), ([3, 3, 3], 31, [S13]), ([3, 3, 3], 43, [S13]), ([3, 3, 3], 47, [S13]),
       ([7, 7], 3, [SC, [[1, 1], [0, 1]]]), ([7, 7], 7, [SC, [[1, 1], [0, 1]]]), ([7, 7], 11, [SC, [[1, 1], [0, 1]]]), ([7, 7], 19, [SC, [[1, 1], [0, 1]]]), ([7, 7], 23, [SC, [[1, 1], [0, 1]]]), ([7, 7], 27, [SC, [[1, 1], [0, 1]]])]
S31 = [[1, 2, 3], [3, 1, 2], [4, 1, 3]]
S171 = [[0, 3, 5], [0, 0, 3], [1, 4, 2]]
SH += [([5, 5, 5], 3, [S31]), ([5, 5, 5], 7, [S31]), ([5, 5, 5], 11, [S31]), ([5, 5, 5], 19, [S31]), ([5, 5, 5], 23, [S31]),
       ([7, 7, 7], 3, [S171]), ([7, 7, 7], 7, [S171]), ([7, 7, 7], 11, [S171]), ([7, 7, 7], 19, [S171])]
os.makedirs(os.path.join(HERE, '..', 'runs'), exist_ok=True)
sel = [int(x) for x in sys.argv[1].split(',')] if len(sys.argv) > 1 else range(len(SH))
for i in sel:
    K, q, hu = SH[i]
    t0 = time.time(); sh = SShape(K, q, hu)
    rec = dict(i=i, K=K, q=q, hu=hu, G=sh.G, nuK=sh.kCD.nu)
    try: rec['gcount'] = gcount(K, q, hu)[0]
    except Exception as ex: rec['gcount'] = None; rec['gerr'] = str(ex)[:80]
    try: rec['kcount'] = kcount(sh, maxu=7)[0]
    except Exception as ex: rec['kcount'] = None; rec['kerr'] = str(ex)[:80]
    big = max([v for v in (rec['gcount'], rec['kcount']) if v is not None] or [0])
    if big > 3e8:
        rec['pipe'] = None
    else:
        try:
            tot, sols = exhaust(sh, tmax=200)
            rec['pipe'] = tot['E_G']; rec['pipe_exhaustive'] = tot['exhaustive']
            if sols: rec['sol_check_G'] = bool(sh.g_target_skew_ok(sh.glevel_blocks(sols[0])))
        except Exception as ex:
            rec['pipe'] = None; rec['perr'] = str(ex)[:100]
    vals = [v for v in (rec['gcount'], rec['kcount'], rec.get('pipe')) if v is not None]
    rec['status'] = 'MATCH' if len(set(vals)) == 1 and len(vals) >= 2 else ('MISMATCH' if len(set(vals)) > 1 else 'PARTIAL')
    rec['secs'] = round(time.time() - t0, 1)
    open(os.path.join(HERE, '..', 'runs', 'val_new.jsonl'), 'a').write(json.dumps(rec) + '\n'); print(rec, flush=True)
