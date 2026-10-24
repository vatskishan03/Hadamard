"""Scan worker: est_extend-style scan of the slots of the top energy window (EREL) of a shape: exact row-sum feasibility (row_patterns_skew) + Gaussian logE/logrho (est_skew.estimate, M samples).
Time-budgeted (cumulative over invocations), resumable.  Pools up to POOL_FULL canonical slots are scanned in a fixed shuffled order (est_extend's), larger pools by uniform random sampling of
(allocation, split) (coverage = processed / ordered-tuple upper bound of the pool).
usage (in search/skewplug): python3 scan38/scan38.py SHAPES.json IDXLIST SECS [EREL=5] [M=60]
writes est/NAME.ext.jsonl (one record per slot: nopat | logN1, logE, logrho, npat) and est/NAME.hdr.json."""
import sys, os, json, time, random, warnings, signal
warnings.filterwarnings('ignore')
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, '..'))      # scan38/ and the modules of search/skewplug
DATA = os.path.join(HERE, 'data')
import numpy as np
from pshape import shape_from_row
from sk_pat import m_splits, row_patterns_skew
from rowfast2 import row_count_fast
from awin import alloc_window, AllocSampler, AllocSamplerE, valid_totals
import est_skew as es
POOL_FULL = int(os.environ.get('POOL_FULL', '60000'))

def shape_name(i, r, wmode='top'):
    return 'S%04d_v%d_q%d_K%s_n%d_o%d' % (i, r['v'], r['q'], 'x'.join(map(str, r['kmods'])), r['nu'], r['order']) + ('' if wmode == 'top' else '_' + wmode)

def pool_ub(win):
    tot = 0
    for (a, em) in win:
        poly = np.array([1.0])
        for j, n in enumerate(a):
            p = np.zeros(n + 1)
            for mm in range(0, n + 1, 2 if j == 0 else 1): p[mm] = 1
            poly = np.convolve(poly, p)
        tot += float(poly[em]) if em < len(poly) else 0.0
    return tot

def sample_slot(win, rng):
    """uniform-ish random (allocation, split): allocation uniform, split by rejection; canonical (B, C, D equal-n sorted as in m_splits)"""
    for _ in range(200):
        a, em = win[rng.randrange(len(win))]
        nA, nB, nC, nD = a
        for _ in range(200):
            mA = 2 * rng.randint(0, nA // 2); mB = rng.randint(0, nB); mC = rng.randint(0, nC); mD = em - mA - mB - mC
            if 0 <= mD <= nD:
                t = sorted([(nB, mB), (nC, mC), (nD, mD)], key=lambda z: (z[0], z[1]))
                return (tuple(a), (mA,) + tuple(m for (_, m) in t))
    return None

def run_shape(i, r, secs, erel=5, M=60, outdir='est', seed=7, wmode='top'):
    name = shape_name(i, r, wmode)
    out = f'{outdir}/{name}.ext.jsonl'; hdr = f'{outdir}/{name}.hdr.json'
    t0 = time.time()
    sh = shape_from_row(r)
    lazy = None
    if wmode == 'lo':      # q = 3: the identity gives no upper bound on E; the top window is the degenerate all-defect corner.  window = E in [N, 1.6 N] (calibrated on the known-solvable 309 = Z103 x F3 <13>: feasible slots and the maximal logE are at E = 1.1 N .. 1.6 N); sampled lazily
        lazy = AllocSampler(sh, sh.N, int(1.6 * sh.N)); win = [None] if lazy.valid else []; emin, emax = (lazy.valid[0], lazy.valid[-1]) if lazy.valid else (None, None)
    elif sh.N > int(os.environ.get('LAZY_N', '150')):      # huge windows (millions of allocations): lazy sampling of allocations in the top window
        vt, emin, emax = valid_totals(sh)
        lazy = AllocSamplerE(sh, [E for E in vt if E >= emax - erel]) if vt else None
        win = [None] if (lazy is not None and lazy.valid) else []
    else:
        win, emin, emax = alloc_window(sh, erel)
    info = dict(window=wmode, name=name, v=r['v'], q=r['q'], kmods=r['kmods'], nu=sh.kCD.nu, order=r['order'], nwin=len(win), erel=erel, M=M)
    if not win:
        info.update(pool=0, status='EXACT_DEAD_identity', secs=round(time.time() - t0, 1), scan_secs=0.0)
        json.dump(info, open(hdr, 'w')); return info
    ub = pool_ub(win) if lazy is None else float('nan')
    info.update(Emin=emin, Emax=emax, pool_ub=ub)
    prior = 0.0
    if os.path.exists(hdr):
        try: prior = float(json.load(open(hdr)).get('scan_secs', 0.0))
        except Exception: prior = 0.0
    budget = max(0.0, secs - prior)
    done = set()
    if os.path.exists(out):
        for l in open(out):
            try: rr = json.loads(l); done.add((tuple(rr['al']), tuple(rr['ms'])))
            except Exception: pass
    rng = random.Random(seed + len(done))
    if lazy is None and ub <= POOL_FULL:
        pool = [(a, ms) for (a, em) in win for ms in m_splits(a, em)]
        random.Random(seed).shuffle(pool)
        mode = 'full'; npool = len(pool)
        info.update(pool=npool, setup_secs=round(time.time() - t0, 1), mode=mode)
        gen = iter(pool)
    else:
        mode = 'sample'; npool = ub
        info.update(pool=npool, setup_secs=round(time.time() - t0, 1), mode=mode)
        def g():
            while True:
                if lazy is not None:
                    al_ = lazy.draw(rng)
                    if al_ is None: return
                    a, em = al_; nA, nB, nC, nD = a
                    for _ in range(200):
                        mA = 2 * rng.randint(0, nA // 2); mB = rng.randint(0, nB); mC = rng.randint(0, nC); mD = em - mA - mB - mC
                        if 0 <= mD <= nD:
                            tt = sorted([(nB, mB), (nC, mC), (nD, mD)], key=lambda z: (z[0], z[1]))
                            yield (tuple(a), (mA,) + tuple(m for (_, m) in tt)); break
                    continue
                s = sample_slot(win, rng)
                if s is None: return
                yield s
        gen = g()
    t1 = time.time(); n = 0; nf = 0; misses = 0
    tot_budget = secs          # cumulative budget (seconds)
    def elapsed(): return prior + (time.time() - t1)
    class _TMO(Exception): pass
    def _alarm(sig, frm): raise _TMO()
    signal.signal(signal.SIGALRM, _alarm)
    TMO_ROW = float(os.environ.get('TMO_ROW', '12')); TMO_EST = float(os.environ.get('TMO_EST', '60'))
    def eval_slot(fo, a, ms):
        rec = dict(al=list(a), ms=list(ms))
        unchecked = False
        try:
            signal.setitimer(signal.ITIMER_REAL, TMO_ROW)
            p = row_count_fast(sh, a, ms)
        except _TMO:
            p = -1; unchecked = True          # block_classes too slow (huge orbit weights): the row-sum filter is skipped, logE is then an UPPER bound
        except Exception as ex:
            p = 0; rec['perr'] = str(ex)[:60]
        finally:
            signal.setitimer(signal.ITIMER_REAL, 0)
        if not p:
            rec['nopat'] = 1
        else:
            try:
                signal.setitimer(signal.ITIMER_REAL, TMO_EST)
                e = es.estimate(sh, tuple(a), tuple(ms), M=M)
                rec.update(logN1=round(e['logN1'], 2), logE=round(e['logE'], 2), logrho=round(e['logrho'], 2), npat=int(p))
                if unchecked: rec['unchecked'] = 1
            except _TMO:
                rec['tmo_est'] = 1
            except Exception as ex:
                rec['err'] = str(ex)[:80]
            finally:
                signal.setitimer(signal.ITIMER_REAL, 0)
        fo.write(json.dumps(rec) + '\n'); fo.flush()
        return rec
    bestal = {}
    def load_best():
        bestal.clear()
        if os.path.exists(out):
            for l in open(out):
                try: rr = json.loads(l)
                except Exception: continue
                if 'logE' in rr:
                    k = tuple(rr['al'])
                    if rr['logE'] > bestal.get(k, -1e9): bestal[k] = rr['logE']
    with open(out, 'a') as fo:
        phaseA_end = 0.6 * tot_budget if mode == 'sample' else tot_budget
        for (a, ms) in gen:
            if elapsed() > phaseA_end: break
            key = (tuple(a), tuple(ms))
            if key in done:
                misses += 1
                if mode == 'sample' and misses > 20000: break
                continue
            done.add(key); n += 1
            eval_slot(fo, a, ms)
        if mode == 'sample' and elapsed() < tot_budget:
            # phase B: all splits of the best allocations found so far (best first)
            load_best()
            for al_ in sorted(bestal, key=lambda k: -bestal[k])[:30]:
                em = lazy.em_of(sum(al_)) if lazy is not None else next(e_ for (a_, e_) in win if tuple(a_) == al_)
                for ms in m_splits(list(al_), em):
                    if elapsed() > tot_budget: break
                    key = (al_, tuple(ms))
                    if key in done: continue
                    done.add(key); n += 1
                    eval_slot(fo, al_, ms)
                if elapsed() > tot_budget: break
    cnt = len(done)
    complete = (mode == 'full' and cnt >= npool)
    info.update(processed=cnt, scan_secs=round(prior + (time.time() - t1), 1), secs=round(time.time() - t0, 1), coverage=(round(cnt / max(1, npool), 5) if npool == npool else None), status='COMPLETE' if complete else 'PARTIAL')
    json.dump(info, open(hdr, 'w'))
    return info

if __name__ == '__main__':
    rows = json.load(open(sys.argv[1]))
    idxs = [int(x) for x in sys.argv[2].split(',')]
    secs = float(sys.argv[3]); erel = int(sys.argv[4]) if len(sys.argv) > 4 else 5; M = int(sys.argv[5]) if len(sys.argv) > 5 else 60
    os.makedirs('est', exist_ok=True)
    for i in idxs:
        try:
            info = run_shape(i, rows[i], secs, erel, M)
        except Exception as ex:
            info = dict(idx=i, error=str(ex)[:200])
        print(json.dumps(info), flush=True)
