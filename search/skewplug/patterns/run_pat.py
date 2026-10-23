"""Exhaustive verdict for a skew plug-in shape by the pattern pipeline: big-big filter (bbpat) -> exact row-sum filter -> symmetry reduction (patsym) ->
for every remaining pattern class the exhaustive level 1 + level 2 of pat_exh.  Resumable: classes already in OUT.jsonl are skipped.  A hit is certified through ../sk_cert.py.
usage (from search/skewplug):  python3 patterns/run_pat.py OUT.jsonl 'K' q 'hunits' tlimit
e.g.  python3 patterns/run_pat.py runs/pat_589.jsonl "[31]" 19 "[[5]]" 600
Result lines: 'alive patterns A, mod symmetry B, done C' then one line per class with n1 (level-1 instances), compl (level-2 completions), status.
The shape is EXHAUSTIVELY empty iff every class has compl = 0 and status ok (and the run ends without 'timeout')."""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, '..')); sys.path.insert(0, HERE)
import json, time, warnings
warnings.filterwarnings('ignore')
import numpy as np
from sk import SShape
from bbpat import enum_patterns, rowfilter
from patsym import Sym
from pat_exh import solve_pattern
from gl1 import L1
from sl2io import header_text
out = sys.argv[1]; K = json.loads(sys.argv[2]); q = int(sys.argv[3]); hu = json.loads(sys.argv[4]); tlim = float(sys.argv[5])
os.makedirs(os.path.dirname(out) or '.', exist_ok=True)
sh = SShape(K, q, hu); t0 = time.time()
r = enum_patterns(sh, tlimit=600)
assert not r['timed_out']
r['alive'], nrm = rowfilter(sh, r['alive'])
sym = Sym(sh); red = sym.reduce(r['alive'])
keys = sorted(red.keys())
done = {}
if os.path.exists(out):
    for l in open(out):
        d = json.loads(l)
        if d.get('K') == K and d.get('q') == q and d.get('hu') == hu and 'pat' in d: done[json.dumps(d['pat'])] = d
print(f'{K} q={q} hu={hu}: alive patterns {len(r["alive"])}, mod symmetry {len(keys)}, done {len(done)}', flush=True)
l1 = L1(sh); hdr = header_text(sh)
allok = True
for i, key in enumerate(keys):
    A, free = key[0], key[1]
    pk = json.dumps([A, free])
    if pk in done: continue
    if time.time() - t0 > tlim: allok = False; break
    rec, sols = solve_pattern(sh, l1, hdr, A, free[0], free[1], free[2])
    rec.update(K=K, q=q, hu=hu, pat=[A, free], idx=i)
    open(out, 'a').write(json.dumps(rec) + '\n')
    print(i, len(keys), {k: v for k, v in rec.items() if k not in ('K', 'hu')}, flush=True)
    if sols:
        for s in sols: assert sh.check_solution(s)
        from sk_cert import certify
        print('HIT', certify(sh, sols[0], tag=f'x{sh.G}'), flush=True)
# summary over all classes of this shape (including those decided in earlier runs of the same OUT file)
cls = {}
for l in open(out) if os.path.exists(out) else []:
    d = json.loads(l)
    if d.get('K') == K and d.get('q') == q and d.get('hu') == hu and 'pat' in d: cls[json.dumps(d['pat'])] = d
mine = [cls[json.dumps([k[0], k[1]])] for k in keys if json.dumps([k[0], k[1]]) in cls]
n1 = sum(d['n1'] for d in mine); compl = sum(d['compl'] for d in mine); bad = sum(1 for d in mine if d['status'] != 'ok')
verdict = ('EXACTLY empty (no pattern survives the filters)' if not keys else 'SOLUTIONS FOUND' if compl else
           'EXHAUSTIVELY empty' if len(mine) == len(keys) and not bad else 'NOT DECIDED (time limit, too many instances or level-2 timeouts)')
print(f'RESULT {K} q={q} hu={hu}: classes {len(keys)}, done {len(mine)}, level-1 instances {n1}, completions {compl}, not-ok classes {bad} -> {verdict}', flush=True)
