"""production sampler for the skew plug-in (one core): rounds of [random allocation / m-split / exact row-sum pattern -> level 1 (exact A list x sampled B, C, D) -> exact prefilter + score ranking -> level 2 (sl2, stop at first completion)].
usage: python3 sk_worker.py NAME SEEDBASE SECS 'K' q ALLOCS_JSON [LB LC LD n2 topfrac]      (K = '[17]', ALLOCS_JSON = '[[nA,nB,nC,nD],...]' or 'E9,10' for all allocations of total weight 9 or 10)
Writes runs/NAME.jsonl (one record per round), on a hit runs/NAME.HIT_<time>.json and a certificate (out/skew_<order>_NAME.json); stops on runs/NAME.STOP or after a hit."""
import sys, os, time, json, random, warnings
warnings.filterwarnings('ignore')
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from sk import SShape
from sk_pat import allocations, m_splits, row_patterns_skew
from sk_samp import level1, level2
from sk_filter import prefilter
from gl1 import L1
from sl2io import header_text
from sk_cert import certify

name, seedbase, secs = sys.argv[1], int(sys.argv[2]), float(sys.argv[3])
Kmods = json.loads(sys.argv[4]); q = int(sys.argv[5]); spec = sys.argv[6]
LB, LC, LD, n2, topfrac = 2000, 40000, 40000, 300, 0.01
if len(sys.argv) > 7: LB, LC, LD, n2 = int(float(sys.argv[7])), int(float(sys.argv[8])), int(float(sys.argv[9])), int(float(sys.argv[10]))
if len(sys.argv) > 11: topfrac = float(sys.argv[11])
os.makedirs('runs', exist_ok=True)
sh = SShape(Kmods, q); l1 = L1(sh); hdr = header_text(sh)
allall = allocations(sh)
if spec.startswith('E'): Es = {int(x) for x in spec[1:].split(',')}; allocs = [(a, em) for (a, em) in allall if sum(a) in Es]
else:
    want = {tuple(a) for a in json.loads(spec)}; allocs = [(a, em) for (a, em) in allall if tuple(a) in want]
pool = []
for (a, em) in allocs:
    for ms in m_splits(a, em): pool.append((a, ms))
log = f'runs/{name}.jsonl'
rnd0 = sum(1 for _ in open(log)) if os.path.exists(log) else 0
t_start = time.time(); rnd = rnd0
cache = {}
print(f'{name}: {len(pool)} (allocation, m-split) pairs, resuming at round {rnd0}', flush=True)
while time.time() - t_start < secs and not os.path.exists(f'runs/{name}.STOP'):
    rng = random.Random(seedbase * 1000003 + rnd)
    al, ms = rng.choice(pool)
    key = (al, ms)
    if key not in cache:
        pats = row_patterns_skew(sh, al, ms); A_list = sh.enum_A(al[0], ms[0]) if pats else None
        cache[key] = (pats, A_list)
    pats, A_list = cache[key]
    rec = dict(round=rnd, al=al, ms=ms, t0=time.time())
    if not pats or A_list is None:
        rec['skip'] = 'no pattern / no A list'; open(log, 'a').write(json.dumps(rec) + '\n'); rnd += 1; continue
    pat = rng.choice(pats)
    rec['pat'] = pat
    t0 = time.time()
    cnt, S, E, st = level1(sh, l1, al, ms, pat, LB, LC, LD, seed=rng.randrange(10 ** 6), A_list=A_list, maxout=100000, hseed=rnd + 1)
    rec['n1'] = int(cnt); rec['t_l1'] = round(time.time() - t0, 2)
    if cnt:
        ok, score = prefilter(sh, S, E)
        rec['n_pass'] = int(ok.sum())
        idx = np.flatnonzero(ok)
        if len(idx):
            idx = idx[np.argsort(score[idx])]
            idx = idx[:max(1, min(n2, int(len(idx) * topfrac) + 1 if topfrac < 1 else n2))]
            t2 = time.time()
            r, sols = level2(sh, S[idx], E[idx], hdr, chunk=25, stop_first=True, tl2=60.0, want=3)
            rec.update(n_l2=int(len(idx)), t_l2=round(time.time() - t2, 2), nok=int(r['nok']), stat=r['stat'])
            if sols:
                for s in sols:
                    assert sh.check_solution(s), 'K-level check failed'
                res = certify(sh, sols[0], tag=name)
                rec['HIT'] = res
                json.dump(dict(blocks=[[x.tolist() for x in b] for b in sols[0]], al=al, ms=ms, cert=res), open(f'runs/{name}.HIT_{int(time.time())}.json', 'w'))
                open(log, 'a').write(json.dumps(rec) + '\n')
                print('HIT', res, flush=True)
                break
    open(log, 'a').write(json.dumps(rec) + '\n')
    rnd += 1
print(f'{name}: done after {rnd - rnd0} rounds', flush=True)
