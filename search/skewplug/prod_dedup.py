"""Production loop for K-symmetric shapes WITH SYMMETRY DEDUPLICATION of the level-1 instances (prod_ks.py plus per-worker engine/list overrides, adaptive slot weights, a pattern bandit and
bounded level 2 per round).  Rounds are seeded by (SEEDBASE, round index); the log runs/NAME.jsonl is replayed at start, so the loop is resumable and a log can be copied to continue or to replay a run.
usage (from this directory):  python3 -u prod_dedup.py NAME SEEDBASE SECS 'K' q 'hunits' SLOTS TOPN LB LC LD ENGINE
  SLOTS   EST    = slots from est/NAME.jsonl (Gaussian estimates est_skew.estimate of a random sample of the top allocations, built on first use: PROD_NEST slots, E window PROD_EREL, PROD_ESTM samples)
          FILE   = an estimate file (a name containing 'est', e.g. data/est473.jsonl: records with logE), or a cost file (a name containing 'cost', e.g. data/cost473.jsonl: measured yield and level-2 cost)
          ALL    = every (allocation, m-split) of the E window (at most 4000 at random; weights x1.5 after a round with new instances, x0.5 after an empty round)
  TOPN    number of best slots kept;  LB LC LD  list sizes of B, C, D;  ENGINE  skf (fast level 2, cyclic K with |K| <= 64) or sl2
What it adds to prod_ks.py:  (a) every level-1 instance is reduced to a canonical 64-bit key under the EXACT symmetries of the K-level system (each block separately: reflection k -> -k; B, C, D also negated
and translated by any point of K fixed by H_K; runs/NAME.seen), so only NEW classes go to level 2 -- the raw level-1 counts of the sampler are dominated by symmetric copies (a group of order up to 2^7 |Fix|^3);
(b) a pattern bandit: a few row-sum patterns of a slot carry nearly all its level-1 instances, 75% of the rounds re-use patterns that produced new classes; (c) slot weights multiplied by the yield of new classes
per wall-clock second (EST mode); (d) at most PROD_L2CAP (default 3000) level-2 instances per round; (e) automatic certification of a hit (sk_cert.certify -> out/skew_<order>_<NAME>.json).
Environment: PROD_EMIN / PROD_EREL (E window), PROD_CONTINUE=1 (go on after a hit), PROD_MAXA (cap of the exact A list, default 3e6), PROD_MAXOUT (stored level-1 instances per round, default 20000),
PROD_L2CAP, PROD_NEST, PROD_ESTM.  Per-worker overrides: runs/NAME.lists ('LB LC LD'), runs/NAME.engine; runs/NAME.STOP ends the loop.
Because the slot weights of EST mode use measured wall-clock times, a FRESH run is not bit-reproducible; a replay of a logged state is (see README: the sixth matrix of order 1892)."""

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
MAXA = int(float(os.environ.get('PROD_MAXA', '3e6')))
os.makedirs('runs', exist_ok=True)
sh = SShape(Kmods, q, hu); l1 = L1(sh); hdr_sl2 = header_text(sh); hdr_skf = skf_header(sh) if (len(Kmods) == 1 and sh.N <= 64) else None

def build_est(sh, path, nslots=800, erel=5, M=100, seed=1):
    """est file: full Gaussian estimate (logN1, logE, logrho) of a random sample of (allocation, m-split) slots among the allocations with E within erel of the maximum"""
    import est_skew as es
    al_ = allocations(sh); emax = max(sum(a) for a, em in al_); emin = emax - erel
    pool = [(a, em) for (a, em) in al_ if sum(a) >= emin]
    rg = random.Random(seed); rg.shuffle(pool); cnt = 0; t0 = time.time()
    with open(path + '.tmp', 'w') as fo:
        for (a, em) in pool:
            if cnt >= nslots: break
            mss = m_splits(a, em); rg.shuffle(mss)
            for ms in mss[:3]:
                cnt += 1
                try:
                    e = es.estimate(sh, tuple(a), tuple(ms), M=M)
                    if e['logE'] > -50: fo.write(json.dumps(dict(al=a, ms=ms, logN1=round(e['logN1'], 2), logE=round(e['logE'], 2), logrho=round(e['logrho'], 2))) + '\n')
                except Exception as ex:
                    pass
    os.replace(path + '.tmp', path)
    print(f'{name}: est built {cnt} slots in {time.time() - t0:.0f}s', flush=True)

# ---- symmetry dedup (EXACT symmetries of the K-level system: each block X separately may be translated by a point fixed by H_K, reflected (k -> -k) and, for B, C, D, negated;
#      A may be reflected; the K-level equations are sums of per-block autocorrelation-type vectors).  A level-1 instance and its images have the same level-2 verdict.
_N = sh.N
_Wh = np.random.default_rng(12345).integers(1, 2 ** 62, size=_N, dtype=np.uint64) * np.uint64(2) + np.uint64(1)
_fix = [int(P[0]) for P in sh.kCD.pos if len(P) == 1]
_neg = sh.negperm
_permsB = [np.asarray(sh.shift[t], dtype=np.int64)[_neg] for t in _fix] + [np.asarray(sh.shift[t], dtype=np.int64) for t in _fix]
_permsA = [np.arange(_N), np.asarray(_neg, dtype=np.int64)]
_MUL = [np.uint64(0x9E3779B97F4A7C15), np.uint64(0xC2B2AE3D27D4EB4F), np.uint64(0x165667B19E3779F9), np.uint64(0x27D4EB2F165667C5)]
def _canon(Sx, Ex, perms, negs):
    best = None
    for sg in ((1, -1) if negs else (1,)):
        v = ((sg * Sx.astype(np.int64) + 2) * 3 + (sg * Ex.astype(np.int64) + 1) + 1).astype(np.uint64)
        for P in perms:
            h = (v[:, P] * _Wh).sum(axis=1, dtype=np.uint64)
            best = h if best is None else np.minimum(best, h)
    return best
def inst_keys(S, E):
    key = np.zeros(len(S), dtype=np.uint64)
    for X in range(4):
        key = key + _canon(S[:, X, :], E[:, X, :], _permsA if X == 0 else _permsB, X > 0) * _MUL[X]
    return key
seen_path = None
slots = []
if slotfile == 'EST':
    os.makedirs('est', exist_ok=True); slotfile = f'est/{name}.jsonl'
    if not os.path.exists(slotfile): build_est(sh, slotfile, nslots=int(os.environ.get('PROD_NEST', '300')), erel=int(os.environ.get('PROD_EREL', '5')), M=int(os.environ.get('PROD_ESTM', '100')))
ALLMODE = (slotfile == 'ALL'); ESTMODE = (slotfile.startswith('est/') or 'est' in os.path.basename(slotfile)) and not COSTONLY
if ALLMODE:
    allocs = allocations(sh)
    emin = int(os.environ.get('PROD_EMIN', '0'))
    if os.environ.get('PROD_EREL') and allocs: emin = max(emin, max(sum(a) for a, em in allocs) - int(os.environ['PROD_EREL']))
    for (a, em) in allocs:
        if sum(a) < emin: continue
        for ms_ in m_splits(a, em): slots.append((0.0, tuple(a), tuple(ms_), 0, 0))
    if len(slots) > 4000: random.Random(seedbase).shuffle(slots); slots = slots[:4000]
for l in (open(slotfile) if not ALLMODE else []):
    r = json.loads(l)
    if 't_l2_inst' in r and r.get('n1_per_round', 0) > 0:
        c1 = r['t_l1_round'] / r['n1_per_round']; eff = 10 ** r['logrho'] / (c1 + r['t_l2_inst'])
        slots.append((math.log10(eff), tuple(r['al']), tuple(r['ms']), r.get('logN1', 0), r.get('logrho', 0)))
    elif 'logE' in r and 't_l2_inst' not in r and 'nopat' not in r and not COSTONLY: slots.append((r['logE'], tuple(r['al']), tuple(r['ms']), r.get('logN1', 0), r.get('logrho', 0)))
if not slots and not ALLMODE:      # est produced no usable slot: fall back to ALL mode (random slots of the E window)
    print(f'{name}: no est slots, falling back to ALL mode', flush=True)
    ALLMODE = True; ESTMODE = False
    allocs = allocations(sh); emin = int(os.environ.get('PROD_EMIN', '0'))
    if os.environ.get('PROD_EREL') and allocs: emin = max(emin, max(sum(a) for a, em in allocs) - int(os.environ['PROD_EREL']))
    for (a, em) in allocs:
        if sum(a) < emin: continue
        for ms_ in m_splits(a, em): slots.append((0.0, tuple(a), tuple(ms_), 0, 0))
    if len(slots) > 4000: random.Random(seedbase).shuffle(slots); slots = slots[:4000]
slots.sort(key=lambda x: -x[0]); slots = slots[:topn]
mx = slots[0][0]; wts = [10 ** (s[0] - mx) for s in slots]; w0 = list(wts); ema = [10.0] * len(wts)      # slot yield EMA: new classes per wall second (optimistic start); the slot weight is w0 * (0.5 + ema)
sidx_of = {(s[1], s[2]): i for i, s in enumerate(slots)}
log = f'runs/{name}.jsonl'
rnd0 = 0
pstat = {}      # (al, ms) -> {json(pat): [rounds, n1]}: pattern bandit (a few row-sum patterns carry all level-1 instances of a slot)
if os.path.exists(log):
    for l in open(log):
        rnd0 += 1
        try: r = json.loads(l)
        except Exception: continue
        k = (tuple(r.get('al', ())), tuple(r.get('ms', ())))
        if 'pat' in r and 'n1' in r:
            e_ = pstat.setdefault(k, {}).setdefault(json.dumps(r['pat']), [0, 0, 0.0]); e_[0] += 1; e_[1] += r.get('n1_new', r['n1']); e_[2] = r.get('n1_new', r['n1']) if e_[0] == 1 else 0.5 * e_[2] + 0.5 * r.get('n1_new', r['n1'])
        if ALLMODE and k in sidx_of and 'n1' in r:      # replay the adaptive weights (ALL mode: multiplicative)
            i = sidx_of[k]; wts[i] = max(wts[i] * 0.5, 0.02 * w0[i]) if r.get('n1_new', r['n1']) == 0 else min(wts[i] * 1.5, 60.0 * max(w0[i], 0.05))
        if ESTMODE and k in sidx_of and 'n1' in r and r.get('t_wall'):
            i = sidx_of[k]; ema[i] = 0.6 * ema[i] + 0.4 * (r.get('n1_new', r['n1']) / max(r['t_wall'], 0.3))
        if r.get('skip') and k in sidx_of: wts[sidx_of[k]] = 0.0 if len(slots) > 1 else 1.0
seen_path = f'runs/{name}.seen'
seen = set(np.fromfile(seen_path, dtype=np.uint64).tolist()) if os.path.exists(seen_path) else set()
print(f'{name}: {len(slots)} slots, best {mx:.2f}, engine {engine}, N {sh.N} nu {sh.kCD.nu}, resume at round {rnd0}', flush=True)
t_start = time.time(); rnd = rnd0; cache = {}
tot = dict(rounds=0, n1=0, l2=0, nok=0, t1=0.0, t2=0.0)
while time.time() - t_start < secs and not os.path.exists(f'runs/{name}.STOP'):
    t_round0 = time.time()
    rng = random.Random(seedbase * 1000003 + rnd)
    if sum(wts) <= 0: print(name + ': no feasible slot left', flush=True); break
    sidx = rng.choices(range(len(slots)), weights=([(w * (0.5 + ema[i]) if w > 0 else 0.0) for i, w in enumerate(wts)] if ESTMODE else wts), k=1)[0]; (logE, al, ms, lN1, lrho) = slots[sidx]
    key = (al, ms)
    if key not in cache:
        try:
            pats = row_patterns_skew(sh, al, ms); A_list = sh.enum_A(al[0], ms[0], cap=MAXA) if pats else None
        except Exception as ex:
            pats, A_list = None, None
        if len(cache) > 60: cache.pop(next(iter(cache)))
        cache[key] = (pats, A_list)
    pats, A_list = cache[key]
    rec = dict(round=rnd, al=al, ms=ms)
    if not pats or A_list is None:
        wts[sidx] = 0.0 if len(slots) > 1 else 1.0; rec['skip'] = 1; rec['t_wall'] = round(time.time() - t_round0, 3); open(log, 'a').write(json.dumps(rec) + '\n'); rnd += 1; continue
    ps_ = pstat.get(key, {})
    prod_ = [(p_, e_) for p_, e_ in ps_.items() if e_[2] >= 1.0]
    if prod_ and rng.random() < 0.6:
        tot_ = sum(e_[2] for p_, e_ in prod_); x_ = rng.random() * tot_; pj = prod_[-1][0]
        for p_, e_ in prod_:
            x_ -= e_[2]
            if x_ <= 0: pj = p_; break
        pat = json.loads(pj); pat = tuple(tuple(x) for x in pat)
    else:
        pat = rng.choice(pats)
    rec['pat'] = pat
    try: LB, LC, LD = [int(x) for x in open(f'runs/{name}.lists').read().split()[:3]]
    except Exception: pass
    try: engine = open(f'runs/{name}.engine').read().strip() or engine
    except Exception: pass
    eng = engine
    if hdr_skf is None: eng = 'sl2'
    t0 = time.time()
    cnt, S, E, st = level1(sh, l1, al, ms, pat, LB, LC, LD, seed=rng.randrange(10 ** 6), A_list=A_list, maxout=int(os.environ.get('PROD_MAXOUT', '20000')), hseed=rnd + 1)
    nnew = 0; newkeys = np.zeros(0, dtype=np.uint64)
    if cnt and len(S):
        keys_ = inst_keys(S, E); uk_, first_ = np.unique(keys_, return_index=True)
        msk_ = np.fromiter((int(k_) not in seen for k_ in uk_), dtype=bool, count=len(uk_))
        sidx_ = first_[msk_]; ord_ = np.argsort(sidx_); sel_ = sidx_[ord_]; newkeys = uk_[msk_][ord_]; S = S[sel_]; E = E[sel_]; nnew = len(sel_)
    t1 = time.time(); rec.update(n1=int(cnt), n1_new=int(nnew), t_l1=round(t1 - t0, 2)); tot['n1'] += int(cnt); tot['t1'] += t1 - t0
    e_ = pstat.setdefault(key, {}).setdefault(json.dumps(pat), [0, 0, 0.0]); e_[0] += 1; e_[1] += int(nnew); e_[2] = float(nnew) if e_[0] == 1 else 0.5 * e_[2] + 0.5 * float(nnew)
    if ALLMODE: wts[sidx] = max(wts[sidx] * 0.5, 0.02 * w0[sidx]) if nnew == 0 else min(wts[sidx] * 1.5, 60.0 * max(w0[sidx], 0.05))
    hit = None
    if nnew:
        nok = 0; ndone = 0; stat = {}
        L2CAP = int(os.environ.get('PROD_L2CAP', '3000'))
        for c0 in range(0, len(S), 100):
            if ndone >= L2CAP or time.time() - t_start > secs + 45: break      # bounded level-2 work per round (the unprocessed instances are not marked as seen)
            insts = [[(S[i, b], E[i, b]) for b in range(4)] for i in range(c0, min(len(S), c0 + 100))]
            if eng in ('skf', 'skf2'): res, ss = run_skf(sh, insts, maxsol=3, tlimit=60.0, hdr=hdr_skf)
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
                json.dump(dict(blocks=[[x.tolist() for x in b] for b in sol], al=al, ms=ms, cert=res_c, Kmods=Kmods, q=q, hu=hu), open(f'runs/{name}.HIT_{int(time.time())}.json', 'w'))
                break
        t2 = time.time(); rec.update(n_l2=ndone, t_l2=round(t2 - t1, 2), nok=int(nok), stat=stat); tot['l2'] += ndone; tot['nok'] += nok; tot['t2'] += t2 - t1
        if not hit:
            seen.update(int(k_) for k_ in newkeys[:ndone]); open(seen_path, 'ab').write(newkeys[:ndone].tobytes())
    tot['rounds'] += 1
    if hit: rec['HIT'] = hit
    rec['t_wall'] = round(time.time() - t_round0, 3)
    if ESTMODE: ema[sidx] = 0.6 * ema[sidx] + 0.4 * (nnew / max(rec['t_wall'], 0.3))
    open(log, 'a').write(json.dumps(rec) + '\n')
    rnd += 1
    if hit:
        print('HIT', hit, flush=True)
        if not os.environ.get('PROD_CONTINUE'): break
    if tot['rounds'] % 20 == 0: print(f'{name}: rounds {tot["rounds"]} n1 {tot["n1"]} level-2 runs {tot["l2"]} completions {tot["nok"]} t1 {tot["t1"]:.0f}s t2 {tot["t2"]:.0f}s', flush=True)
print(f'{name}: done: ' + json.dumps(tot), flush=True)
