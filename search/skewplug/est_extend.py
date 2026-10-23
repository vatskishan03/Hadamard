"""extend est/NAME.jsonl: every slot (allocation, m-split) of the top E window, exact row-sum feasibility, full Gaussian logE (est_skew.estimate, M samples).
usage: est_extend.py NAME 'K' q 'hu' M SECS [EREL]   -> est/NAME.ext.jsonl (append, resume-safe); merge with est_merge.py"""
import sys, json, random, time, warnings, os, math
warnings.filterwarnings('ignore')
from sk import SShape
from sk_pat import allocations, m_splits, row_patterns_skew
import est_skew as es
name = sys.argv[1]; K = json.loads(sys.argv[2]); q = int(sys.argv[3]); hu = json.loads(sys.argv[4]); M = int(sys.argv[5]); secs = float(sys.argv[6]); erel = int(sys.argv[7]) if len(sys.argv) > 7 else 5
sh = SShape(K, q, hu); al = allocations(sh); emax = max(sum(a) for a, em in al)
pool = [(a, ms) for (a, em) in al if sum(a) >= emax - erel for ms in m_splits(a, em)]
rng = random.Random(int(os.environ.get('EST_SEED', '7'))); rng.shuffle(pool)
out = f'est/{name}.ext.jsonl' if os.environ.get('EST_SEED', '7') == '7' else f'est/{name}.ext{os.environ["EST_SEED"]}.jsonl'; done = set()
if os.path.exists(out):
    for l in open(out):
        try: r = json.loads(l); done.add((tuple(r['al']), tuple(r['ms'])))
        except Exception: pass
t0 = time.time(); n = 0; nf = 0
for (a, ms) in pool:
    if time.time() - t0 > secs: break
    key = (tuple(a), tuple(ms))
    if key in done: continue
    n += 1
    rec = dict(al=list(a), ms=list(ms))
    try:
        p = row_patterns_skew(sh, a, ms)
    except Exception as ex:
        p = None
    if not p:
        rec['nopat'] = 1
    else:
        nf += 1
        try:
            e = es.estimate(sh, tuple(a), tuple(ms), M=M)
            rec.update(logN1=round(e['logN1'], 2), logE=round(e['logE'], 2), logrho=round(e['logrho'], 2), npat=len(p))
        except Exception as ex:
            rec['err'] = str(ex)[:80]
    open(out, 'a').write(json.dumps(rec) + '\n')
    if n % 20 == 0: print(name, n, nf, round(time.time() - t0), flush=True)
print(name, 'done', n, nf, round(time.time() - t0), 's of', len(pool), 'slots', flush=True)
