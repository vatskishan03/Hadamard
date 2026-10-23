"""EXACT row-sum (DC) feasibility of the slots of a shape.  A slot is an (allocation, m-split) pair; ../sk_pat.row_patterns_skew lists the exact row-sum patterns of the four blocks
(sum of squares 4|G|, parity, defect contributions).  A slot without a pattern carries no solution; a shape in which NO slot has a pattern is EXACTLY empty (a sample of slots proves nothing).
usage (from search/skewplug):  python3 patterns/rowscan.py 'K' q 'hunits' NSAMPLE [EREL [SECS]]
  NSAMPLE slots are drawn at random from the pool (all E levels, or only the top EREL levels below the maximum E); give NSAMPLE >= the pool size to test every slot.
Output: pool size, number tested, number with at least one row-sum pattern, and the time.  With NSAMPLE = 100000:
  "[109]" 11 "[[16]]" -> pool 3056, 0 with patterns (1199);  "[109]" 19 "[[16]]" -> pool 5, 0 (2071);  "[73]" 11 "[[2]]" -> pool 317, 0 (803);  "[127]" 11 "[[22]]" -> pool 7870, 0 (1397, the pool takes longer than 240 s: SECS = 3000)."""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, '..')); sys.path.insert(0, HERE)
import json, random, time, warnings
warnings.filterwarnings('ignore')
from sk import SShape
from sk_pat import allocations, m_splits, row_patterns_skew
K = json.loads(sys.argv[1]); q = int(sys.argv[2]); hu = json.loads(sys.argv[3]); ns = int(sys.argv[4]); erel = int(sys.argv[5]) if len(sys.argv) > 5 else 99
sh = SShape(K, q, hu); al = allocations(sh); emax = max(sum(a) for a, em in al)
pool = [(a, ms) for (a, em) in al if sum(a) >= emax - erel for ms in m_splits(a, em)]
rng = random.Random(1); rng.shuffle(pool); t0 = time.time(); good = 0; tested = 0; np_ = 0
for (a, ms) in pool[:ns]:
    if time.time() - t0 > float(sys.argv[6] if len(sys.argv) > 6 else 240): break
    tested += 1
    p = row_patterns_skew(sh, a, ms)
    if p: good += 1; np_ += len(p)
print(f'{K} q={q} hu={hu}: pool {len(pool)} slots, tested {tested}, with row-sum patterns {good}, patterns {np_}, {time.time()-t0:.0f}s')
