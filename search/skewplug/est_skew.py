"""Gaussian / lattice local-limit estimator of the number of solutions of a slot of the skew plug-in (HEURISTIC; sampled block statistics, full covariance, per-feature lattice
spacing from sampled sums).  For an allocation (nA,nB,nC,nD) with m-split ms:  E_alloc = T_A T_B T_C T_D * P(sum of the four block feature vectors = target),   features
  F = [I(0); I(k), J(k) (k in repsPM, k != 0); Q(k) (k in repsQ)]   (target (4|G|; 0; 0; 0)),      level 1 only (U system):  N1 = T^se_A ... * P(sum U(k) = 0, k != 0).
estimate(sh, al, ms) returns log10 of E_alloc, N1 and rho = E/N1 (ordered quadruples, no sign gauge: B, C, D and the t-gauge are in the counts).
At v = 473 (K = Z43, q = 11, H_K = <6>) the level-1 value is within a factor 1.5 of the measured one (estimates/n1_est.py) but E is 3 to 4 decades too low: the estimate ranks slots,
it cannot discard a shape.  Library module (sk_filter-style helpers for sampling: sample_A, sample_free, feats, gauss_log10, spacing_of, T_A, T_free);
usage: python3 est_skew.py '(nA,nB,nC,nD)' '(mA,mB,mC,mD)'      (K = Z17, q = 7, H_K = 1), e.g. '(2,1,4,6)' '(2,0,2,5)'"""
import numpy as np, math, itertools, random
from math import comb


def count_DM(weights, n, m):
    """number of (D, M): D subset of items with total weight n, M subset of D with weight m, and the weighted sums  sum 2^|D|, sum 2^|D| 4^(items-|D|)"""
    from collections import defaultdict
    st = {(0, 0, 0): 1}          # (wD, wM, |D|) -> count
    for w in weights:
        new = defaultdict(int)
        for (wd, wm, k), c in st.items():
            new[(wd, wm, k)] += c                       # nondefect
            if wd + w <= n:
                new[(wd + w, wm, k + 1)] += c           # defect p
                if wm + w <= m: new[(wd + w, wm + w, k + 1)] += c      # defect m
        st = new
    return {(k): c for (wd, wm, k), c in st.items() if wd == n and wm == m}


def T_free(sh, n, m, full):
    ws = [int(w) for w in sh.kCD.w]; nu = len(ws)
    cnt = count_DM(ws, n, m)
    return sum(c * (2 ** k) * ((4 if full else 2) ** (nu - k)) for k, c in cnt.items())


def T_A(sh, n, m, full):
    ws = [2 * len(sh.spos[o]) for (o, po) in sh.apairs]; npr = len(ws)
    cnt = count_DM(ws, n, m)
    tot = 0
    for k, c in cnt.items():
        tot += c * ((2 ** k) * ((4 if full else 2) ** (npr - k)))
    return tot * (2 if full else 1)


def sample_free(sh, n, m, M, rng):
    K = sh.kCD; N = sh.N; ws = [int(w) for w in K.w]
    s, t, e = K.sample(n, M, rng)                       # defect set uniform; types iid
    s = s.astype(np.int64); t = t.astype(np.int64); e = e.astype(np.int64)
    for i in range(M):
        d = [o for o in range(K.nu) if e[i, K.pos[o][0]] != 0]
        for _ in range(200):
            rng.shuffle(d) if hasattr(rng, 'shuffle') else None
            perm = list(np.random.default_rng(int(rng.integers(1 << 30))).permutation(d)); tot = 0; Ms = []
            for o in perm:
                if tot + ws[o] <= m: Ms.append(o); tot += ws[o]
            if tot == m: break
        for o in d:
            P = K.pos[o]; ev = e[i, P[0]]
            s[i, P] = -2 * ev if o in Ms else 0
    return s, t, e


def sample_A(sh, n, m, M, rng):
    N = sh.N; pairs = sh.apairs; pos = sh.spos; npr = len(pairs); ws = [2 * len(pos[o]) for (o, po) in pairs]
    S = np.zeros((M, N), np.int64); T = np.zeros((M, N), np.int64); E = np.zeros((M, N), np.int64)
    S[:, 0] = 1; T[:, 0] = rng.integers(0, 2, M) * 2 - 1
    # all (D, M) pairs by rejection: random subset of pairs with weight n (uniform over subsets by DP-free rejection when small)
    import itertools as it
    subs = []
    def rec(j, rem, cur):
        if rem == 0: subs.append(list(cur)); return
        if j == npr or rem < 0: return
        cur.append(j); rec(j + 1, rem - ws[j], cur); cur.pop(); rec(j + 1, rem, cur)
    rec(0, n, [])
    pickl = []
    for D in subs:
        for r in range(len(D) + 1):
            for Mm in it.combinations(D, r):
                if sum(ws[j] for j in Mm) == m: pickl.append((D, set(Mm)))
    for i in range(M):
        D, Mm = pickl[int(rng.integers(len(pickl)))]
        for j, (o, po) in enumerate(pairs):
            be = 1 if rng.integers(2) else -1; bt = 1 if rng.integers(2) else -1
            if j in D:
                ev = be; sv = -2 * ev if j in Mm else 0
                S[i, pos[o]] = sv; S[i, pos[po]] = -sv; E[i, pos[o]] = ev; E[i, pos[po]] = -ev
            else:
                S[i, pos[o]] = be; S[i, pos[po]] = -be
                T[i, pos[o]] = bt; T[i, pos[po]] = bt
    return S, T, E


def feats(sh, S, T, E, which):
    q, eps = sh.q, sh.eps; cols = []
    u = S + q * E
    if which == 'U':
        for k in sh.repsPM:
            if k == 0: continue
            cols.append((u * u[:, sh.shift[k]]).sum(1))
        return np.stack(cols, 1).astype(float)
    cols.append((S * S).sum(1) + q * (T * T).sum(1))      # I(0)
    for k in sh.repsPM:
        if k == 0: continue
        sh_ = sh.shift[k]
        cols.append((S * S[:, sh_]).sum(1) + q * (T * T[:, sh_]).sum(1))
    for k in sh.repsPM:
        if k == 0: continue
        sh_ = sh.shift[k]
        cols.append(q * (E * E[:, sh_]).sum(1) + (S * E[:, sh_]).sum(1) + (E * S[:, sh_]).sum(1) - (T * T[:, sh_]).sum(1))
    for k in sh.repsQ:
        sh_ = sh.shift[k]
        cols.append((S * T[:, sh_]).sum(1) + eps * (T * S[:, sh_]).sum(1))
    return np.stack(cols, 1).astype(float)


def gauss_log10(F4, target, spacing=None):
    """F4: list of 4 arrays (M, d) of independent block features; Gaussian density of the SUM at target, per-feature lattice spacing from the sampled sums"""
    mu = sum(f.mean(0) for f in F4)
    cov = sum(np.cov(f.T) for f in F4)
    d = len(mu)
    if spacing is None:
        spacing = np.ones(d)
    # robust: add tiny ridge to the diagonal (degenerate features have zero variance: they are constant -> must equal target exactly)
    var = np.diag(cov)
    const = var < 1e-9
    ok = np.abs(mu[const] - target[const]) < 1e-6
    if not ok.all(): return -np.inf
    keep = ~const
    c = cov[np.ix_(keep, keep)]; m_ = mu[keep]; tg = target[keep]; sp = spacing[keep]
    sign, logdet = np.linalg.slogdet(c + 1e-9 * np.eye(keep.sum()))
    z = tg - m_
    q = z @ np.linalg.solve(c + 1e-9 * np.eye(keep.sum()), z)
    lg = np.log(sp).sum() - 0.5 * keep.sum() * math.log(2 * math.pi) - 0.5 * logdet - 0.5 * q
    return lg / math.log(10)


def spacing_of(F4, M=400, rng=None):
    """per-feature gcd of differences of sums of four independent random blocks"""
    rng = rng or np.random.default_rng(1)
    d = F4[0].shape[1]
    sums = sum(f[rng.integers(0, len(f), M)] for f in F4)
    out = np.ones(d)
    for j in range(d):
        v = np.rint(sums[:, j]).astype(np.int64); diffs = np.abs(v[1:] - v[0]); g = 0
        for x in diffs: g = math.gcd(g, int(x))
        out[j] = g if g > 0 else 1
    return out


def estimate(sh, al, ms, M=1500, seed=1):
    rng = np.random.default_rng(seed)
    nA, nB, nC, nD = al; mA, mB, mC, mD = ms
    blocks = [sample_A(sh, nA, mA, M, rng)] + [sample_free(sh, n, m, M, rng) for n, m in ((nB, mB), (nC, mC), (nD, mD))]
    q = sh.q; G = sh.G
    out = {}
    # level 1
    FU = [feats(sh, *b, 'U') for b in blocks]
    sp = spacing_of(FU, rng=rng)
    tgU = np.zeros(FU[0].shape[1])
    l1 = gauss_log10(FU, tgU, sp)
    TA1 = T_A(sh, nA, mA, False); TF1 = [T_free(sh, n, m, False) for n, m in ((nB, mB), (nC, mC), (nD, mD))]
    logT1 = math.log10(TA1) + sum(math.log10(x) for x in TF1)
    # full system
    FF = [feats(sh, *b, 'all') for b in blocks]
    spF = spacing_of(FF, rng=rng)
    tg = np.zeros(FF[0].shape[1]); tg[0] = 4 * G
    lf = gauss_log10(FF, tg, spF)
    TAf = T_A(sh, nA, mA, True); TFf = [T_free(sh, n, m, True) for n, m in ((nB, mB), (nC, mC), (nD, mD))]
    logTf = math.log10(TAf) + sum(math.log10(x) for x in TFf)
    out = dict(logN1=logT1 + l1, logE=logTf + lf, logT1=logT1, logTf=logTf, logP1=l1, logPf=lf)
    out['logrho'] = out['logE'] - out['logN1']
    return out


if __name__ == '__main__':
    import sys, json
    from sk import SShape
    sh = SShape([17], 7)
    al = eval(sys.argv[1]); ms = eval(sys.argv[2])
    r = estimate(sh, al, ms)
    print(json.dumps({k: round(v, 2) for k, v in r.items()}))
