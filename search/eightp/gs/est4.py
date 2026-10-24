"""est4.py -- Gaussian/lattice estimator for the REAL plug-in (four free blocks, all of kind CD), adapted from est2.py (search/NOTES.md sections 82, 85).
Per allocation a = (n1 <= n2 <= n3 <= n4) of defect weights: E_a = ordered number of K-level solutions (block sign gauge included), N1_a = number of level-1 (s,e) instances, rho_a = E_a / N1_a.
HEURISTIC (Gaussian local limit + lattice index); calibrated on the exact counts of odd-|G| Spence shapes in calib_odd.py."""
import numpy as np, math, itertools, sys, time
from k4 import *
from est import kind_dpu, sample_mode, lognblocks, logsumexp, LN10
from est2 import EnsLazy


def perm_mult(a):
    from collections import Counter
    m = 24
    for v in Counter(a).values(): m //= math.factorial(v)
    return m


def defect_window(sh):
    r = sh.q; D = sh.Did
    if D < 0: return 1, 0
    lo = -(-D // (r + 1)); hi = D // (r - 3) if r > 3 else 4 * sh.N
    return lo, min(hi, 4 * sh.N)


def enum_allocs(W, lo, hi):
    """all sorted 4-tuples of the sorted array W with lo <= sum <= hi, as an (n,4) int array"""
    W = np.array(sorted(W)); nW = len(W); out = []
    for i in range(nW):
        for j in range(i, nW):
            s2 = W[i] + W[j]
            if s2 + 2 * W[j] > hi + 10 ** 9: pass
            for k in range(j, nW):
                s3 = s2 + W[k]
                # l >= k and lo <= s3 + W[l] <= hi
                l0 = max(k, int(np.searchsorted(W, lo - s3, 'left'))); l1 = int(np.searchsorted(W, hi - s3, 'right'))
                if l1 > l0:
                    ls = np.arange(l0, l1)
                    out.append(np.stack([np.full(len(ls), W[i]), np.full(len(ls), W[j]), np.full(len(ls), W[k]), W[ls]], 1))
    return np.concatenate(out) if out else np.zeros((0, 4), int)


class Est4:
    def __init__(self, sh, M=1000, Mlat=None, seed=1):
        self.sh = sh; self.K = sh.kCD; self.M = M
        self.Mlat = Mlat if Mlat else max(sh.R + 40, 120)
        self.rng = np.random.default_rng(seed)
        self.cache = {}
        self.lo, self.hi = defect_window(sh)
        W = [n for n in self.K.weights() if n <= self.hi]
        self.step = 1
        if len(W) > 36:
            self.step = -(-len(W) // 36)
            W = W[::self.step]
        self.W = W
        self.logscale = 4 * math.log(self.step)

    def ens(self, mode, n):
        key = (mode, n)
        if key not in self.cache:
            self.cache[key] = EnsLazy(self.sh, self.K, n, mode, self.M, self.Mlat, self.rng)
        return self.cache[key]

    def cols(self, mode):
        return self.sh.nU if mode == 'u' else self.sh.R

    def target(self, mode):
        return self.sh.T[:self.cols(mode)]

    def single(self, mode, a):
        """lattice-corrected log10 E for the allocation a (sorted 4-tuple); returns (logE, gauss_only_logE, tbits) ; logE = -inf if T not in the lattice coset"""
        T = self.target(mode)
        keys = [self.ens(mode, n) for n in a]
        v0 = sum(k.v0 for k in keys)
        lat = keys[0].lat.copy()
        for k in keys[1:]: lat.add_many(k.lat.rows())
        mult = perm_mult(a)
        logcnt = sum(k.logn for k in keys) * LN10 + math.log(mult)
        mu = sum(k.mean for k in keys); S = sum(k.cov for k in keys)
        w, U = np.linalg.eigh(S); keep = w > 1e-7 * max(w.max(), 1e-12)
        d = U[:, keep].T @ (T - mu)
        lg = -0.5 * float((d ** 2 / w[keep]).sum()) - 0.5 * float(np.log(2 * np.pi * w[keep]).sum())
        tb = sum(k.tbits for k in keys)
        if not lat.contains((T - v0).round().astype(np.int64).tolist()):
            return -np.inf, (lg + logcnt) / LN10, tb
        Bm = np.array(lat.rows(), np.float64); r = len(Bm)
        Uu, ss, Vt = np.linalg.svd(Bm, full_matrices=False); Qb = Vt[:r].T
        Cm = Qb.T @ S @ Qb; dm = Qb.T @ (T - mu)
        w, Uc = np.linalg.eigh(Cm); keep = w > 1e-7 * max(w.max(), 1e-12)
        dm2 = Uc[:, keep].T @ dm
        logphi = -0.5 * float((dm2 ** 2 / w[keep]).sum()) - 0.5 * float(np.log(2 * np.pi * w[keep]).sum())
        lv = lat.logvol()
        return (logcnt + logphi + lv) / LN10, (lg + logcnt) / LN10, tb

    def total(self, mode='all', nfull=250, topk=8, verbose=False):
        sh = self.sh
        A = enum_allocs(self.W, self.lo, self.hi)
        if len(A) == 0: return dict(logE=-np.inf, nalloc=0, best=None, tbits=None, top=[])
        wset = sorted(set(A.ravel().tolist()))
        T = self.target(mode)
        c2 = [0, sh.nU] if mode != 'u' else [0]
        mu2 = {n: self.ens(mode, n).mean[c2] for n in wset}; C2 = {n: self.ens(mode, n).cov[np.ix_(c2, c2)] if len(c2) > 1 else np.atleast_2d(self.ens(mode, n).cov[0, 0]) for n in wset}
        logn = {n: self.ens(mode, n).logn for n in wset}
        idx = {n: i for i, n in enumerate(wset)}
        MU = np.array([mu2[n] for n in wset]); CC = np.array([C2[n] for n in wset]); LN = np.array([logn[n] for n in wset])
        ia = np.vectorize(idx.get)(A)
        mu = MU[ia].sum(1); S = CC[ia].sum(1) + 1e-9 * np.eye(len(c2))
        d = (T[c2][None, :] - mu)
        sol = np.linalg.solve(S, d[:, :, None])[:, :, 0]
        quad = (d * sol).sum(1)
        logdet = np.linalg.slogdet(S)[1]
        mult = np.array([perm_mult(tuple(a)) for a in A]) if len(A) < 200000 else np.full(len(A), 12.0)
        sc1 = -0.5 * quad - 0.5 * logdet + math.log(2 * math.pi) * (-len(c2) / 2.0) + np.log(mult) + LN[ia].sum(1) * LN10
        order = np.argsort(sc1)[::-1]
        sel = order[:nfull]
        # stage 2: full gaussian
        full_scores = {}
        keys_all = {}
        for i in sel:
            a = tuple(int(x) for x in A[i])
            keys = [self.ens(mode, n) for n in a]
            mu_ = sum(k.mean for k in keys); S_ = sum(k.cov for k in keys)
            w, U = np.linalg.eigh(S_); keep = w > 1e-7 * max(w.max(), 1e-12)
            dd = U[:, keep].T @ (T - mu_)
            lg = -0.5 * float((dd ** 2 / w[keep]).sum()) - 0.5 * float(np.log(2 * np.pi * w[keep]).sum())
            full_scores[i] = lg + sum(k.logn for k in keys) * LN10 + math.log(perm_mult(a))
        ordf = sorted(full_scores, key=lambda i: -full_scores[i])
        # hopeless shapes: even with the largest plausible lattice gain (0.45 decades per column) the Gaussian total is far below 1
        logEg = (logsumexp([full_scores[i] for i in ordf]) + self.logscale) / LN10
        if logEg + 0.45 * self.cols(mode) < -8:
            i0 = ordf[0]; a0 = tuple(int(x) for x in A[i0])
            return dict(logE=logEg, nalloc=len(A), best=a0, tbits=sum(self.ens(mode, n).tbits for n in a0), top=[(logEg, a0)], alive=(0, 0), skipped_lattice=True)
        # stage 3: lattice on the top
        ev = {}
        for i in ordf[:topk]:
            a = tuple(int(x) for x in A[i]); le, lg, tb = self.single(mode, a)
            ev[i] = (le, tb, a)
        good = [(i, v) for i, v in ev.items() if np.isfinite(v[0])]
        if good:
            corr = float(np.median([(v[0] * LN10 - full_scores[i]) for i, v in good])); frac = len(good) / len(ev)
        else:
            corr = 0.0; frac = 0.0
        terms = [v[0] * LN10 + self.logscale for i, v in good]
        # evaluated-but-not-lattice-checked full scores, and the rest by stage-1 (corrected by stage1->full median offset)
        off12 = float(np.median([full_scores[i] - sc1[i] for i in ordf[:min(len(ordf), 60)]]))
        for i in ordf[topk:]:
            terms.append(full_scores[i] + corr + math.log(max(frac, 1e-300)) + self.logscale)
        rest = [i for i in order[nfull:]]
        if rest:
            r1 = sc1[rest] + off12 + corr + math.log(max(frac, 1e-300))
            terms.append(logsumexp(r1) + self.logscale)
        logE = logsumexp(terms) / LN10 if terms else -np.inf
        best = max(good, key=lambda iv: iv[1][0]) if good else None
        top = sorted([(v[0], v[2]) for v in ev.values()], key=lambda x: -x[0])[:3]
        return dict(logE=logE, nalloc=len(A), best=(best[1][2] if best else None), tbits=(best[1][1] if best else None), top=top, alive=(len(good), len(ev) - len(good)))


def estimate_shape(sh, M=1000, nfull=250, topk=8, seed=1, ntop=3):
    """full estimate: total log10 E, log10 N1 (level-1 instances, total), and for the best `ntop` allocations (E_a, N1_a, rho_a, n_t)"""
    es = Est4(sh, M=M, seed=seed)
    r = es.total('all', nfull=nfull, topk=topk)
    out = dict(logE=r['logE'], nalloc=r['nalloc'], tbits=r['tbits'], best=r['best'], alive=r.get('alive'))
    if not np.isfinite(r['logE']) or r['best'] is None or r.get('skipped_lattice'):
        out['hopeless'] = bool(r.get('skipped_lattice'))
        return out
    rows = []
    for (le, a) in r['top']:
        if not np.isfinite(le): continue
        n1, _, _ = es.single('u', a)
        rows.append(dict(alloc=list(a), logE=float(le), logN1=float(n1), tbits=float(es.ens('all', a[0]).tbits + sum(es.ens('all', x).tbits for x in a[1:]))))
    out['allocs'] = rows
    r1 = es.total('u', nfull=nfull, topk=topk)
    out['logN1'] = r1['logE']
    return out
