"""estimate2: Monte-Carlo-over-allocations version of est.estimate with lazy lattices (for big shapes with wide defect windows)."""
import numpy as np, math, itertools, sys, time
from pk import *
from est import kind_dpu, sample_mode, lognblocks, logsumexp, LN10


class EnsLazy:
    def __init__(self, sh, K, n, mode, M, Mlat, rng):
        self.sh, self.K, self.n, self.mode, self.M, self.Mlat = sh, K, n, mode, M, Mlat
        cols = sh.nU if mode == 'u' else (sh.nU + 1 if mode == 'u0' else sh.R)
        self.cols = cols
        s, t, e = sample_mode(K, n, max(M, Mlat), rng, 'u' if mode == 'u0' else mode)
        V = sh.eqvec(s, t, e, mode if mode in ('u', 'u0') else 'all').astype(np.float64)
        self.mean = V[:M].mean(0)
        self.cov = np.cov(V[:M].T) if cols > 1 else np.array([[V[:M, 0].var()]])
        self.V0 = V[:Mlat].astype(np.int64)
        self.v0 = self.V0[0]
        self.logn = lognblocks(K, n, mode)
        firstpos = np.array([int(p[0]) for p in K.pos])
        self.tbits = float((t[:M][:, firstpos] != 0).sum(1).mean())
        self._lat = None

    @property
    def lat(self):
        if self._lat is None:
            lat = Lat(self.cols)
            lat.add_many(list(np.unique(self.V0[1:] - self.V0[0], axis=0)))
            self._lat = lat
            self.V0 = None
        return self._lat


def alloc_tables(W, wmax):
    """counts of ordered-by-size pairs (n1 <= n2) by sum, and lists"""
    pairs = {}
    for i, a in enumerate(W):
        for b in W[i:]:
            pairs.setdefault(a + b, []).append((a, b))
    return pairs


def estimate2(sh, mode='all', Malloc=2500, M=3000, Mlat=500, topk=24, seed=1, verbose=False, Ewin=None, allocs_override=None):
    rng = np.random.default_rng(seed)
    cols = sh.nU if mode == 'u' else (sh.nU + 1 if mode == 'u0' else sh.R)
    T = sh.T[:cols] if mode != 'u0' else np.append(sh.T[:sh.nU], 0.0)
    Kp = sh.N - (1 if sh.bordered else 0)
    q = sh.q
    lo = -(-4 * Kp // (q + 1)); hi = (4 * Kp) // (q - 3) if q > 3 else 4 * sh.N
    if Ewin is not None: lo, hi = Ewin
    WA = [n for n in sh.kAB.weights() if n <= hi]; WC = [n for n in sh.kCD.weights() if n <= hi]
    pA = alloc_tables(WA, hi); pC = alloc_tables(WC, hi)
    cA = {s: len(v) for s, v in pA.items()}; cC = {s: len(v) for s, v in pC.items()}
    combos = [(sA, sC, cA[sA] * cC[sC]) for sA in cA for sC in cC if lo <= sA + sC <= hi]
    total = float(sum(c for _, _, c in combos))
    if total == 0: return dict(logE=-np.inf, n_alloc=0, alive=(0, 0), best_alloc=None, tbits=None, corr=0.0, top=[])
    # sample allocations
    if allocs_override is not None:
        allocs = list(allocs_override); total = float(len(allocs)); wmul = 1.0
    elif total <= Malloc:
        allocs = []
        for sA, sC, c in combos:
            for (a, b) in pA[sA]:
                for (cc, d) in pC[sC]:
                    allocs.append((a, b, cc, d))
        wmul = 1.0
    else:
        probs = np.array([c for _, _, c in combos], float); probs /= probs.sum()
        idx = rng.choice(len(combos), size=Malloc, p=probs)
        allocs = []
        for i in idx:
            sA, sC, c = combos[i]
            a, b = pA[sA][rng.integers(len(pA[sA]))]; cc, d = pC[sC][rng.integers(len(pC[sC]))]
            allocs.append((a, b, cc, d))
        wmul = total / Malloc
    ens = {}
    def get(kind, n):
        key = (kind, n)
        if key not in ens:
            K = sh.kAB if kind == 'A' else sh.kCD
            ens[key] = EnsLazy(sh, K, n, mode, M, Mlat, rng)
        return ens[key]
    t0 = time.time()
    scores = []
    for a in allocs:
        keys = [get('A', a[0]), get('A', a[1]), get('C', a[2]), get('C', a[3])]
        mu = sum(k.mean for k in keys); S = sum(k.cov for k in keys)
        w, U = np.linalg.eigh(S)
        keep = w > 1e-7 * max(w.max(), 1e-12)
        d = U[:, keep].T @ (T - mu)
        lg = -0.5 * float((d ** 2 / w[keep]).sum()) - 0.5 * float(np.log(2 * np.pi * w[keep]).sum())
        mult = (1 if a[0] == a[1] else 2) * (1 if a[2] == a[3] else 2)
        logcnt = sum(k.logn for k in keys) * LN10 + math.log(mult)
        scores.append(lg + logcnt)
    if verbose: print(f'  {len(allocs)} sampled allocations of {total:.3g}, {len(ens)} ensembles [{time.time()-t0:.0f}s]', flush=True)
    order = np.argsort(scores)[::-1]
    full = {}
    for idx_ in order[:topk]:
        a = allocs[idx_]
        keys = [get('A', a[0]), get('A', a[1]), get('C', a[2]), get('C', a[3])]
        v0 = sum(k.v0 for k in keys)
        lat = keys[0].lat.copy()
        for k in keys[1:]: lat.add_many(k.lat.rows())
        if not lat.contains((T - v0).round().astype(np.int64).tolist()):
            full[idx_] = None; continue
        mu = sum(k.mean for k in keys); S = sum(k.cov for k in keys)
        Bm = np.array(lat.rows(), np.float64); r = len(Bm)
        Uu, ss, Vt = np.linalg.svd(Bm, full_matrices=False); Qb = Vt[:r].T
        Cm = Qb.T @ S @ Qb; dm = Qb.T @ (T - mu)
        w, Uc = np.linalg.eigh(Cm)
        keep = w > 1e-7 * max(w.max(), 1e-12)
        dm2 = Uc[:, keep].T @ dm
        logphi = -0.5 * float((dm2 ** 2 / w[keep]).sum()) - 0.5 * float(np.log(2 * np.pi * w[keep]).sum())
        lv = lat.logvol()
        mult = (1 if a[0] == a[1] else 2) * (1 if a[2] == a[3] else 2)
        logcnt = sum(k.logn for k in keys) * LN10 + math.log(mult)
        full[idx_] = (logcnt + logphi + lv, logphi + lv - scores[idx_] + logcnt, r)
    ev = [v for v in full.values() if v is not None]
    ne = sum(1 for v in full.values() if v is None)
    terms = [full[i][0] for i in order[:topk] if full[i] is not None]
    if ev:
        corr = float(np.median([v[1] for v in ev])); frac_alive = len(ev) / (len(ev) + ne)
    else:
        corr = 0.0; frac_alive = 0.0
    rest = [scores[i] + corr + math.log(max(frac_alive, 1e-300)) for i in order[topk:]]
    allterms = terms + rest
    logsum = logsumexp(allterms) if allterms else -np.inf
    # scale: sampled allocations stand for `total` allocations
    if total > Malloc: logsum += math.log(total / len(allocs))
    logE = logsum / LN10
    tb = None; besta = None
    if ev:
        bi = max((i for i in order[:topk] if full[i] is not None), key=lambda i: full[i][0])
        a = allocs[bi]; besta = a
        tb = sum(get(k, n).tbits for k, n in (('A', a[0]), ('A', a[1]), ('C', a[2]), ('C', a[3])))
    return dict(logE=logE, n_alloc=total, alive=(len(ev), ne), best_alloc=besta, tbits=tb, corr=corr / LN10,
                top=[(allocs[i], (full[i][0] / LN10 if full[i] else None)) for i in order[:3]])


if __name__ == '__main__':
    N = int(sys.argv[1]); q = int(sys.argv[2]); bord = sys.argv[3] == 'b'
    hg = [tuple(int(x) for x in a.split(':')) if ':' in a else (int(a), 1) for a in sys.argv[4].split(',')] if len(sys.argv) > 4 and sys.argv[4] else []
    sh = Shape(N, q, hg, bord)
    t0 = time.time()
    r = estimate2(sh, 'all', verbose=True)
    print('E:', r, f'[{time.time()-t0:.0f}s]')
    r1 = estimate2(sh, 'u')
    print('N1:', r1['logE'], ' rho = 10^', r['logE'] - r1['logE'], ' tbits', r['tbits'])
