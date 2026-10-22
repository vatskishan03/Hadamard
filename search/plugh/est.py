"""Gaussian / lattice local-limit estimator for orbit-constant plug-in systems (ksym.Shape).
E   = sum over defect-weight allocations (nA, nB, nC, nD) of  #blocks * P(sum of four independent block vectors = T)      (full system)
N1  = same with the U-columns only and (s,e)-blocks only (distinct level-1 instances);  rho = E / N1 = completions per level-1 instance.
Also returns the typical number of t-sign bits n_t (level-2 MITM cost ~ 2^(n_t/2)).
All log10 of ORDERED counts in the ensemble of ksym (block sign gauge included).  HEURISTIC (Gaussian local limit; see search/NOTES.md 80-82)."""
import numpy as np, math, itertools, sys, time
from pk import *

LN10 = math.log(10)


def logsumexp(a):
    a = np.asarray(a, float)
    if len(a) == 0: return -np.inf
    m = a.max()
    return m + math.log(np.exp(a - m).sum())


class Ens:
    """statistics of one (kind, weight n) ensemble for mode 'all' (full vector, uniform over blocks) or 'u' (U columns, uniform over (s,e)-configs)"""
    pass


def kind_dpu(K):
    """DP over free units for the (s,e)-count: factor per unit: nondefect 2 (1 for A/B zero unit), defect 4 (2 for zero unit)"""
    nu = K.nu
    fi = K.fi
    maxw = K.shape.N
    dp = [[0.0] * (maxw + 1) for _ in range(len(fi) + 1)]
    dp[len(fi)][0] = 1.0
    for ii in range(len(fi) - 1, -1, -1):
        i = fi[ii]; wi = int(K.w[i]); zu = (i == K.zero_unit)
        fn = 1.0 if zu else 2.0; fd = 2.0 if zu else 4.0
        for s_ in range(maxw + 1):
            dp[ii][s_] = fn * dp[ii + 1][s_] + (fd * dp[ii + 1][s_ - wi] if s_ >= wi else 0.0)
    # forced units contribute defect factor only
    ff = 1.0
    for i in range(nu):
        if K.forced[i]: ff *= (2.0 if i == K.zero_unit else 4.0)
    return dp, ff


def sample_mode(K, n, M, rng, mode):
    """mode 'u': defect subsets drawn with probability proportional to 2^(#nondef) 4^(#def) (distinct (s,e) configs)"""
    if mode == 'all': return K.sample(n, M, rng)
    if not hasattr(K, '_dpu'): K._dpu = kind_dpu(K)
    dp, ff = K._dpu
    fn = np.array([1.0 if i == K.zero_unit else 2.0 for i in range(K.nu)]); fd = np.array([2.0 if i == K.zero_unit else 4.0 for i in range(K.nu)])
    return K.sample(n, M, rng, dp=dp, fn=fn, fd=fd)


def lognblocks(K, n, mode):
    if mode == 'all': return K.lognblocks(n)
    if not hasattr(K, '_dpu'): K._dpu = kind_dpu(K)
    dp, ff = K._dpu
    c = dp[0][n - K.fw]
    if c <= 0: return None
    return math.log10(c) + math.log10(ff)


def build_ens(sh, K, n, mode, M, Mlat, rng):
    cols = sh.nU if mode == 'u' else sh.R
    s, t, e = sample_mode(K, n, max(M, Mlat), rng, mode)
    V = sh.eqvec(s, t, e, 'u' if mode == 'u' else 'all').astype(np.float64)
    out = Ens()
    out.mean = V[:M].mean(0); out.cov = np.cov(V[:M].T) if cols > 1 else np.array([[V[:M, 0].var()]])
    Vi = V[:Mlat].astype(np.int64)
    lat = Lat(cols); lat.add_many(list(np.unique(Vi[1:] - Vi[0], axis=0)))
    out.lat = lat; out.v0 = Vi[0]; out.logn = lognblocks(K, n, mode)
    # typical number of t-sign bits: nondefect non-forced units that carry a tau bit
    # (defect indicator not stored; recompute from the sample: t != 0 on the unit's first position)
    firstpos = np.array([int(p[0]) for p in K.pos])
    tb = (t[:M][:, firstpos] != 0).sum(1)
    out.tbits = float(tb.mean())
    # number of free x0/s bits
    return out


def estimate(sh, mode='all', M=6000, Mlat=1200, seed=1, topk=25, verbose=False, Ewin=None, maxalloc=400000):
    rng = np.random.default_rng(seed)
    cols = sh.nU if mode == 'u' else sh.R
    T = sh.T[:cols]
    Kp = sh.N - (1 if sh.bordered else 0)
    q = sh.q
    lo = -(-4 * Kp // (q + 1)); hi = (4 * Kp) // (q - 3) if q > 3 else 4 * sh.N
    if Ewin is not None: lo, hi = Ewin
    WA = [n for n in sh.kAB.weights() if n <= hi]; WC = [n for n in sh.kCD.weights() if n <= hi]
    # prune weights so that they can complete the window
    ens = {}
    t0 = time.time()
    needA = set(); needC = set()
    allocs = []
    for i, nA in enumerate(WA):
        for nB in WA[i:]:
            for j, nC in enumerate(WC):
                for nD in WC[j:]:
                    E = nA + nB + nC + nD
                    if E < lo or E > hi: continue
                    allocs.append((nA, nB, nC, nD))
    if len(allocs) > maxalloc:
        raise RuntimeError(f'too many allocations {len(allocs)}')
    for a in allocs:
        needA.update([a[0], a[1]]); needC.update([a[2], a[3]])
    for n in sorted(needA): ens[('A', n)] = build_ens(sh, sh.kAB, n, mode, M, Mlat, rng)
    for n in sorted(needC): ens[('C', n)] = build_ens(sh, sh.kCD, n, mode, M, Mlat, rng)
    if verbose: print(f'  sampled {len(ens)} ensembles, {len(allocs)} allocations [{time.time()-t0:.0f}s]', flush=True)
    # Gaussian-only scoring
    scores = []
    for a in allocs:
        keys = [('A', a[0]), ('A', a[1]), ('C', a[2]), ('C', a[3])]
        mu = sum(ens[k].mean for k in keys); S = sum(ens[k].cov for k in keys)
        w, U = np.linalg.eigh(S)
        keep = w > 1e-7 * max(w.max(), 1e-12)
        d = U[:, keep].T @ (T - mu)
        lg = -0.5 * float((d ** 2 / w[keep]).sum()) - 0.5 * float(np.log(2 * np.pi * w[keep]).sum())
        mult = (1 if a[0] == a[1] else 2) * (1 if a[2] == a[3] else 2)
        logcnt = sum(ens[k].logn for k in keys) * LN10 + math.log(mult)
        scores.append(lg + logcnt)
    order = np.argsort(scores)[::-1]
    # full evaluation of the top allocations
    full = {}
    for idx in order[:topk]:
        a = allocs[idx]
        keys = [('A', a[0]), ('A', a[1]), ('C', a[2]), ('C', a[3])]
        v0 = sum(ens[k].v0 for k in keys)
        lat = ens[keys[0]].lat.copy()
        for k in keys[1:]:
            lat.add_many(ens[k].lat.rows())
        inlat = lat.contains((T - v0).round().astype(np.int64).tolist())
        if not inlat:
            full[idx] = None; continue
        mu = sum(ens[k].mean for k in keys); S = sum(ens[k].cov for k in keys)
        Bm = np.array(lat.rows(), np.float64); r = len(Bm)
        Uu, ss, Vt = np.linalg.svd(Bm, full_matrices=False); Qb = Vt[:r].T
        Cm = Qb.T @ S @ Qb; dm = Qb.T @ (T - mu)
        w, Uc = np.linalg.eigh(Cm)
        keep = w > 1e-7 * max(w.max(), 1e-12)
        dm2 = Uc[:, keep].T @ dm
        logphi = -0.5 * float((dm2 ** 2 / w[keep]).sum()) - 0.5 * float(np.log(2 * np.pi * w[keep]).sum())
        lv = lat.logvol()
        mult = (1 if a[0] == a[1] else 2) * (1 if a[2] == a[3] else 2)
        logcnt = sum(ens[k].logn for k in keys) * LN10 + math.log(mult)
        full[idx] = (logcnt + logphi + lv, logphi + lv - scores[idx] + logcnt, r)
    # lattice correction for the rest: median of (full - gaussian) over evaluated nonempty ones; empty (coset) fraction
    ev = [v for v in full.values() if v is not None]
    ne = sum(1 for v in full.values() if v is None)
    tot_terms = []
    for idx in order[:topk]:
        v = full[idx]
        if v is not None: tot_terms.append(v[0])
    if ev:
        corr = float(np.median([v[1] for v in ev]))        # (lattice vol + gaussian-restriction) correction relative to the gaussian-only score
        frac_alive = len(ev) / (len(ev) + ne)
    else:
        corr = 0.0; frac_alive = 0.0
    rest = [scores[i] + corr + math.log(max(frac_alive, 1e-300)) for i in order[topk:]]
    allterms = tot_terms + rest
    logE = logsumexp(allterms) / LN10 if allterms else -np.inf
    # typical t-bits at the dominant allocation(s)
    tb = None
    if ev:
        bestidx = max((i for i in order[:topk] if full[i] is not None), key=lambda i: full[i][0])
        a = allocs[bestidx]
        keys = [('A', a[0]), ('A', a[1]), ('C', a[2]), ('C', a[3])]
        tb = sum(ens[k].tbits for k in keys); besta = a
    else:
        besta = None
    info = dict(logE=logE, n_alloc=len(allocs), alive=(len(ev), ne), best_alloc=besta, tbits=tb, corr=corr / LN10,
                top=[(allocs[i], (full[i][0] / LN10 if full[i] else None)) for i in order[:5]])
    return info


def shape_stats(sh):
    nA = sh.kAB.nu; nC = sh.kCD.nu
    return nA, nC


if __name__ == '__main__':
    N = int(sys.argv[1]); q = int(sys.argv[2]); bord = sys.argv[3] == 'b'
    hg = [tuple(int(x) for x in a.split(':')) if ':' in a else (int(a), 1) for a in sys.argv[4].split(',')] if len(sys.argv) > 4 and sys.argv[4] else []
    sh = Shape(N, q, hg, bord)
    t0 = time.time()
    r = estimate(sh, 'all', verbose=True)
    print('E:', r, f'[{time.time()-t0:.0f}s]')
    r1 = estimate(sh, 'u', verbose=False)
    print('N1:', r1['logE'])
    print('rho = 10^', r['logE'] - r1['logE'])
