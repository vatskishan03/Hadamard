"""Exhaustive verdict for a skew plug-in shape by the allocation pipeline (../sk_exh.py: explicit lists for every allocation and m-split, exact join, exact level 2).
usage (from search/skewplug):  python3 patterns/exh_shape.py OUT.jsonl 'K' q 'hunits' tmax
e.g.  python3 patterns/exh_shape.py runs/exh_979.jsonl "[89]" 11 "[[2]]" 600
Environment: MAXINST (level-1 instances above which an allocation is skipped, default 3e6), PAIRS_HI (join size above which it is skipped, default 3e9), PAIRS_LO (skip the joins up to that size).
The RESULT line gives the number of allocations, E_G (ordered level-2 completions, = the number of K-level solutions) and the skipped allocations;
the shape is EXHAUSTIVELY empty iff exhaustive is true and nskipped = 0 and found = 0."""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, '..')); sys.path.insert(0, HERE)
import json, time, warnings
warnings.filterwarnings('ignore')
from sk import SShape
from sk_exh import exhaust
out = sys.argv[1]; K = json.loads(sys.argv[2]); q = int(sys.argv[3]); hu = json.loads(sys.argv[4]); tmax = float(sys.argv[5])
os.makedirs(os.path.dirname(out) or '.', exist_ok=True)
sh = SShape(K, q, hu)
t0 = time.time()
tot, sols = exhaust(sh, tmax=tmax, verbose=True, max_inst=int(float(os.environ.get('MAXINST', 3e6))), pairs_lo=float(os.environ.get('PAIRS_LO', 0)), pairs_hi=float(os.environ.get('PAIRS_HI', 3e9)))
rec = dict(K=K, q=q, hu=hu, nu=sh.kCD.nu, allocs=tot['allocs'], exhaustive=tot['exhaustive'], E_G=tot['E_G'], N1_G=tot['N1_G'], skipped=[(list(a), str(r)) for a, r in tot['skipped']][:20], nskipped=len(tot['skipped']), flags=len(tot['flags']), secs=round(time.time() - t0, 1), found=len(sols))
print('RESULT', json.dumps(rec), flush=True)
open(out, 'a').write(json.dumps(rec) + '\n')
if sols:
    from sk_cert import certify
    print(certify(sh, sols[0], tag=f'x{sh.G}'))
