"""production sampler for a K-symmetric (or small) skew plug-in shape (one core, resumable, stops at the first hit):
round = [slot drawn with weight ~ 10^logE (Gaussian per-slot estimate file) -> exact row-sum pattern -> level 1 (gl1: exact A list x sampled B, C, D lists) -> level 2 on ALL level-1 instances (skf = fast level 2, or sl2) -> exact K-level check -> certificate].
usage: prod_ks.py NAME SEEDBASE SECS 'K' q 'hunits' SLOTFILE TOPN LB LC LD [engine]
Writes runs/NAME.jsonl (one record per round), runs/NAME.HIT_*.json and out/skew_<order>_NAME.json on a hit; stops on runs/NAME.STOP."""
import sys, os, time, json, random, warnings, math
warnings.filterwarnings('ignore')
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from sk import SShape
from sk_pat import allocations, m_splits, row_patterns_skew
from sk_samp import level1
from gl1 import L1
from sl2io import header_text, run_sl2
from skf_io import skf_header, run_skf
name, seedbase, secs = sys.argv[1], int(sys.argv[2]), float(sys.argv[3])
Kmods = json.loads(sys.argv[4]); q = int(sys.argv[5]); hu = json.loads(sys.argv[6]); slotfile = sys.argv[7]; topn = int(sys.argv[8])
LB, LC, LD = int(float(sys.argv[9])), int(float(sys.argv[10])), int(float(sys.argv[11])); engine = sys.argv[12] if len(sys.argv) > 12 else 'skf'; COSTONLY = ('cost' in os.path.basename(slotfile))
os.makedirs('runs', exist_ok=True)
sh = SShape(Kmods, q, hu); l1 = L1(sh); hdr_sl2 = header_text(sh); hdr_skf = skf_header(sh)
slots = []
ALLMODE = (slotfile == 'ALL')
if ALLMODE:
    for (a, em) in allocations(sh):
        for ms_ in m_splits(a, em): slots.append((0.0, tuple(a), tuple(ms_), 0, 0))
    emin = int(os.environ.get('PROD_EMIN', '0'))
    if os.environ.get('PROD_EREL') and slots: emin = max(emin, max(sum(x[1]) for x in slots) - int(os.environ['PROD_EREL']))
    slots = [x for x in slots if sum(x[1]) >= emin]
    if len(slots) > 4000: random.Random(seedbase).shuffle(slots); slots = slots[:4000]
for l in (open(slotfile) if not ALLMODE else []):
    r = json.loads(l)
    if 't_l2_inst' in r and r.get('n1_per_round', 0) > 0:      # cost file: weight = Gaussian rho / (level-1 + level-2 seconds per instance)  (HEURISTIC efficiency)
        c1 = r['t_l1_round'] / r['n1_per_round']; eff = 10 ** r['logrho'] / (c1 + r['t_l2_inst'])
        slots.append((math.log10(eff), tuple(r['al']), tuple(r['ms']), r.get('logN1', 0), r.get('logrho', 0)))
    elif 'logE' in r and 't_l2_inst' not in r and 'nopat' not in r and not COSTONLY: slots.append((r['logE'], tuple(r['al']), tuple(r['ms']), r.get('logN1', 0), r.get('logrho', 0)))
slots.sort(key=lambda x: -x[0]); slots = slots[:topn]
mx = slots[0][0]; wts = [10 ** (s[0] - mx) for s in slots]
yield_ = {}
print(f'{name}: {len(slots)} slots, best logE {mx:.2f}, worst kept {slots[-1][0]:.2f}, engine {engine}', flush=True)
log = f'runs/{name}.jsonl'
rnd0 = sum(1 for _ in open(log)) if os.path.exists(log) else 0
t_start = time.time(); rnd = rnd0; cache = {}
tot = dict(rounds=0, n1=0, l2=0, nok=0, t1=0.0, t2=0.0)
while time.time() - t_start < secs and not os.path.exists(f'runs/{name}.STOP'):
    rng = random.Random(seedbase * 1000003 + rnd)
    if sum(wts) <= 0: print(name + ': no feasible slot left (all slots skipped: no row-sum pattern / A list too large)', flush=True); break
    sidx = rng.choices(range(len(slots)), weights=wts, k=1)[0]; (logE, al, ms, lN1, lrho) = slots[sidx]
    key = (al, ms)
    if key not in cache:
        try:
            pats = row_patterns_skew(sh, al, ms); A_list = sh.enum_A(al[0], ms[0], cap=3_000_000) if pats else None      # too large an exact A list (> 3M) = slot skipped (a sampled A is not implemented)
        except Exception as ex:
            pats, A_list = None, None
        if len(cache) > 60: cache.pop(next(iter(cache)))      # bounded memory
        cache[key] = (pats, A_list)
    pats, A_list = cache[key]
    rec = dict(round=rnd, al=al, ms=ms)
    if not pats or A_list is None:
        wts[sidx] = 0.0 if len(slots) > 1 else 1.0; rec['skip'] = 1; open(log, 'a').write(json.dumps(rec) + '\n'); rnd += 1; continue
    pat = rng.choice(pats); rec['pat'] = pat
    try: engine = open('engine.txt').read().strip() or engine
    except Exception: pass
    try: LB, LC, LD = [int(x) for x in open('lists.txt').read().split()[:3]]      # optional override of the list sizes (higher yield per pair for rare level-1 instances)
    except Exception: pass
    if sh.N > 64: engine = 'sl2'      # skf (pl2f lineage) is limited to n <= 64
    t0 = time.time()
    cnt, S, E, st = level1(sh, l1, al, ms, pat, LB, LC, LD, seed=rng.randrange(10 ** 6), A_list=A_list, maxout=20000, hseed=rnd + 1)
    t1 = time.time(); rec.update(n1=int(cnt), t_l1=round(t1 - t0, 2)); tot['n1'] += int(cnt); tot['t1'] += t1 - t0
    if ALLMODE: wts[sidx] = max(wts[sidx] * 0.5, 0.02) if cnt == 0 else min(wts[sidx] * 1.5, 60.0)      # exploit slots that produce level-1 instances
    hit = None
    if cnt:
        nok = 0; ndone = 0; stat = {}
        for c0 in range(0, len(S), 100):
            insts = [[(S[i, b], E[i, b]) for b in range(4)] for i in range(c0, min(len(S), c0 + 100))]
            if engine in ('skf', 'skf2'): res, ss = run_skf(sh, insts, maxsol=3, tlimit=60.0, hdr=hdr_skf)
            else: res, ss = run_sl2(sh, insts, maxsol=3, tlimit=60.0, maxfib=1500000, stop_first=True, hdr=hdr_sl2)
            ndone += len(insts)
            for i, d in res.items():
                stat[d['status']] = stat.get(d['status'], 0) + 1; nok += d['compl'] > 0
            if ss:
                i, ts = ss[0]
                sol = [(insts[i][X][0].astype(np.int64), ts[X].astype(np.int64), insts[i][X][1].astype(np.int64)) for X in range(4)]
                assert sh.check_solution(sol), 'K-level check failed'
                from sk_cert import certify
                res_c = certify(sh, sol, tag=(name + '_r' + str(rnd)) if os.environ.get('PROD_CONTINUE') else name); hit = res_c
                json.dump(dict(blocks=[[x.tolist() for x in b] for b in sol], al=al, ms=ms, cert=res_c), open(f'runs/{name}.HIT_{int(time.time())}.json', 'w'))
                break
        t2 = time.time(); rec.update(n_l2=ndone, t_l2=round(t2 - t1, 2), nok=int(nok), stat=stat); tot['l2'] += ndone; tot['nok'] += nok; tot['t2'] += t2 - t1
    tot['rounds'] += 1
    if hit: rec['HIT'] = hit
    open(log, 'a').write(json.dumps(rec) + '\n')
    rnd += 1
    if hit:
        print('HIT', hit, flush=True)
        if not os.environ.get('PROD_CONTINUE'): break
    if tot['rounds'] % 20 == 0: print(f'{name}: rounds {tot["rounds"]} n1 {tot["n1"]} level-2 runs {tot["l2"]} completions {tot["nok"]} t1 {tot["t1"]:.0f}s t2 {tot["t2"]:.0f}s', flush=True)
print(f'{name}: done: ' + json.dumps(tot), flush=True)
