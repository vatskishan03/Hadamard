"""Brute-force counters for skew Goethals-Seidel quadruples over G = K x F_q (q prime = 3 mod 4), blocks invariant under H = {1} x S and optional H_K (units of K, phi = +1).
   gcount(...)  : INDEPENDENT G-level count (raw +-1 sequences on G, integer PAFs; no Paley algebra).  A(0) = 1, A(-g) = -A(g); B, C, D free.
   kcount(...)  : K-level count from the (I), (J), (Q) equations (s,t,e blocks), all configurations, no sign reduction.
Both return the exact number of ORDERED quadruples (A, B, C, D) with sum PAF = 4|G| delta."""
import numpy as np, itertools, sys, time, json


def count4b(FA, FF, target, maxad=6e7):
    """same count with the table on (A, D) sums (size dA * nF) and a loop over (b, c): memory dA*nF instead of nF^2"""
    uA, cA = np.unique(FA, axis=0, return_counts=True)
    uF, cF = np.unique(FF, axis=0, return_counts=True)
    nF = len(uF)
    if len(uA) * nF > maxad: raise RuntimeError('AD table too large')
    rng = np.random.default_rng(778); d = FA.shape[1]
    c1 = rng.integers(1, 2 ** 61, d, dtype=np.int64) | 1; c2 = rng.integers(1, 2 ** 61, d, dtype=np.int64) | 1
    hF1 = (uF * c1).sum(1); hF2 = (uF * c2).sum(1); hA1 = (uA * c1).sum(1); hA2 = (uA * c2).sum(1)
    t1 = np.int64((np.asarray(target, dtype=np.int64) * c1).sum()); t2 = np.int64((np.asarray(target, dtype=np.int64) * c2).sum())
    P1 = (hA1[:, None] + hF1[None, :]).ravel(); P2 = (hA2[:, None] + hF2[None, :]).ravel()
    W = (cA[:, None].astype(np.int64) * cF[None, :].astype(np.int64)).ravel()
    o = np.argsort(P1, kind='stable'); P1 = P1[o]; P2 = P2[o]; W = W[o]
    total = 0
    for b in range(nF):
        need1 = t1 - hF1[b] - hF1; need2 = t2 - hF2[b] - hF2
        lo = np.searchsorted(P1, need1, 'left'); hi = np.searchsorted(P1, need1, 'right')
        sel = np.flatnonzero(hi > lo)
        if not len(sel): continue
        ln = (hi - lo)[sel]
        idx = np.repeat(lo[sel], ln) + (np.arange(ln.sum()) - np.repeat(np.cumsum(ln) - ln, ln))
        csel = np.repeat(sel, ln)
        good = P2[idx] == need2[csel]
        total += int(cF[b]) * sum(int(cF[c]) * int(w) for c, w in zip(csel[good], W[idx[good]]))
    return total, dict(nA=len(FA), nF=len(FF), dA=len(uA), dF=nF)


def count4(FA, FF, target, maxpairs=2.5e7):
    """FA, FF: int64 feature matrices (rows = configs) of the skew block and of the free blocks; count (a,b,c,d) with FA[a]+FF[b]+FF[c]+FF[d] = target"""
    uA, cA = np.unique(FA, axis=0, return_counts=True)
    uF, cF = np.unique(FF, axis=0, return_counts=True)
    nF = len(uF)
    if nF * nF > maxpairs: return count4b(FA, FF, target)
    rng = np.random.default_rng(777)
    d = FA.shape[1]
    c1 = rng.integers(1, 2 ** 61, d, dtype=np.int64) | 1; c2 = rng.integers(1, 2 ** 61, d, dtype=np.int64) | 1
    hF1 = (uF * c1).sum(1); hF2 = (uF * c2).sum(1); hA1 = (uA * c1).sum(1); hA2 = (uA * c2).sum(1)
    t1 = np.int64((np.asarray(target, dtype=np.int64) * c1).sum()); t2 = np.int64((np.asarray(target, dtype=np.int64) * c2).sum())
    P1 = (hF1[:, None] + hF1[None, :]).ravel(); P2 = (hF2[:, None] + hF2[None, :]).ravel()
    W = (cF[:, None].astype(object) * cF[None, :].astype(object)).ravel() if nF * nF < 0 else (cF[:, None].astype(np.int64) * cF[None, :].astype(np.int64)).ravel()
    o = np.argsort(P1, kind='stable'); P1 = P1[o]; P2 = P2[o]; W = W[o]
    total = 0
    for a in range(len(uA)):
        need1 = t1 - hA1[a] - hF1; need2 = t2 - hA2[a] - hF2
        lo = np.searchsorted(P1, need1, 'left'); hi = np.searchsorted(P1, need1, 'right')
        sel = np.flatnonzero(hi > lo)
        if not len(sel): continue
        ln = (hi - lo)[sel]
        idx = np.repeat(lo[sel], ln) + (np.arange(ln.sum()) - np.repeat(np.cumsum(ln) - ln, ln))
        bsel = np.repeat(sel, ln)
        good = P2[idx] == need2[bsel]
        total += int(cA[a]) * sum(int(cF[b]) * int(w) for b, w in zip(bsel[good], W[idx[good]]))
    return total, dict(nA=len(FA), nF=len(FF), dA=len(uA), dF=nF)


# ------------------------------------------------------------------------------------------- G level
def build_G(Kmods, q, hunits=()):
    from pk import KGrp
    from idx import GF, factor_prime_power
    Kg = KGrp(Kmods); N = Kg.n; n = N * q
    p_, e_ = factor_prime_power(q)
    if e_ == 1:
        sq = sorted({x * x % q for x in range(1, q)})
        mulf = lambda c, a: (c * a) % q
        negf = lambda a: (-a) % q
    else:
        F = GF(p_, e_)
        sq = sorted({int(F.mul[x, x]) for x in range(1, q)})
        mulf = lambda c, a: int(F.mul[c, a])
        def negf(a): return F.idx([(-d) % p_ for d in F.dig(a)])
    maps = []
    for c in sq:
        maps.append(np.array([k * q + mulf(c, a) for k in range(N) for a in range(q)]))
    for u in hunits:
        kp = Kg.mulperm(tuple(u) if not np.isscalar(u) else u)
        maps.append(np.array([int(kp[k]) * q + a for k in range(N) for a in range(q)]))
    neg = np.array([int(Kg.neg[k]) * q + negf(a) for k in range(N) for a in range(q)])
    mods = list(Kmods) + [p_] * e_
    coords = np.array(np.unravel_index(np.arange(n), mods)).T
    m = np.array(mods); w = np.ones(len(mods), dtype=np.int64)
    for i in range(len(mods) - 2, -1, -1): w[i] = w[i + 1] * mods[i + 1]
    return Kg, N, n, maps, neg, coords, m, w


def orbit_labels(n, maps):
    parent = list(range(n))
    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]; x = parent[x]
        return x
    for pm in maps:
        for x in range(n):
            a, b = find(x), find(int(pm[x]))
            if a != b: parent[a] = b
    roots = {}; lab = np.zeros(n, dtype=np.int64)
    for x in range(n):
        r = find(x)
        if r not in roots: roots[r] = len(roots)
        lab[x] = roots[r]
    return lab, len(roots)


def paf_feats(X, coords, m, w, reps):
    X = X.astype(np.int32); n = X.shape[1]
    out = np.zeros((X.shape[0], len(reps)), dtype=np.int64)
    for j, g in enumerate(reps):
        sh = (((coords + coords[g][None, :]) % m) * w).sum(1)
        out[:, j] = (X * X[:, sh]).sum(1)
    return out


def gcount(Kmods, q, hunits=(), maxorb=21):
    Kg, N, n, maps, neg, coords, m, w = build_G(Kmods, q, hunits)
    lab, nor = orbit_labels(n, maps)
    if nor > maxorb: raise RuntimeError(f'too many orbits {nor}')
    # free blocks
    M = 1 << nor; ids = np.arange(M, dtype=np.int64)
    vals = np.zeros((M, nor), dtype=np.int8)
    for j in range(nor): vals[:, j] = 1 - 2 * ((ids >> j) & 1)
    XF = vals[:, lab]
    # skew blocks: orbits O, -O distinct (pairs); A(0) = 1
    olab_neg = {}
    for o in range(nor):
        g = int(np.flatnonzero(lab == o)[0]); olab_neg[o] = int(lab[neg[g]])
    zero = int(lab[0]); assert olab_neg[zero] == zero and (lab == zero).sum() == 1
    pairs = []; seen = {zero}
    for o in range(nor):
        if o in seen: continue
        po = olab_neg[o]; assert po != o
        pairs.append((o, po)); seen |= {o, po}
    MA = 1 << len(pairs); ida = np.arange(MA, dtype=np.int64)
    va = np.zeros((MA, nor), dtype=np.int8); va[:, zero] = 1
    for j, (o, po) in enumerate(pairs):
        b = 1 - 2 * ((ida >> j) & 1); va[:, o] = b; va[:, po] = -b
    XA = va[:, lab]
    assert (XA[:, neg] == -XA)[:, 1:].all() if True else True
    reps = [g for g in range(n) if g <= int(neg[g])]
    FA = paf_feats(XA, coords, m, w, reps); FF = paf_feats(XF, coords, m, w, reps)
    tg = np.zeros(len(reps), dtype=np.int64); tg[reps.index(0)] = 4 * n
    tot, info = count4(FA, FF, tg)
    info.update(orbits=nor, G=n)
    return tot, info


# ------------------------------------------------------------------------------------------- K level
def kcount(sh, maxu=7):
    """K-level count: all (s,t,e) blocks, features (I over repsPM, J over repsPM, Q over repsQ)."""
    from sk import SShape
    N, q = sh.N, sh.q
    K = sh.kCD
    nu = K.nu
    if nu > maxu: raise RuntimeError(f'too many K-orbits {nu}')
    # free blocks: per orbit 8 variants
    var = []
    for i in range(nu):
        v = []
        for be in (1, -1):
            for bt in (1, -1):
                v.append((be, 0, bt))            # nondefect: s = be, e = 0, t = bt
        for be in (1, -1):
            for typ in (0, 1):
                v.append((0 if typ == 0 else -2 * be, be, 0))   # defect: s, e, t = 0
        var.append(v)
    idx = np.indices([8] * nu).reshape(nu, -1).T
    M = len(idx)
    S = np.zeros((M, N), np.int64); E = np.zeros((M, N), np.int64); T = np.zeros((M, N), np.int64)
    arr = [np.array(v, dtype=np.int64) for v in var]
    for i in range(nu):
        P = K.pos[i]; sel = arr[i][idx[:, i]]
        S[:, P] = sel[:, 0:1]; E[:, P] = sel[:, 1:2]; T[:, P] = sel[:, 2:3] * K.sig[P][None, :]
    FF = feats_K(sh, S, T, E)
    # skew blocks
    pairs = sh.apairs; npair = len(pairs); pos = sh.spos
    vA = []
    for j in range(npair):
        v = []
        for be in (1, -1):
            for bt in (1, -1): v.append((be, 0, bt))
        for be in (1, -1):
            for typ in (0, 1): v.append((0 if typ == 0 else -2 * be, be, 0))
        vA.append(np.array(v, dtype=np.int64))
    idxA = np.indices([8] * npair).reshape(npair, -1).T if npair else np.zeros((1, 0), dtype=int)
    MA = len(idxA) * 2
    SA = np.zeros((MA, N), np.int64); EA = np.zeros((MA, N), np.int64); TA = np.zeros((MA, N), np.int64)
    SA[:, 0] = 1
    TA[:, 0] = np.tile([1, -1], len(idxA)) if False else np.repeat([1, -1], len(idxA))
    base = np.concatenate([idxA, idxA]) if npair else np.zeros((2, 0), dtype=int)
    for j, (o, po) in enumerate(pairs):
        sel = vA[j][base[:, j]]
        SA[:, pos[o]] = sel[:, 0:1]; EA[:, pos[o]] = sel[:, 1:2]
        SA[:, pos[po]] = -sel[:, 0:1]; EA[:, pos[po]] = -sel[:, 1:2]
        TA[:, pos[o]] = sel[:, 2:3]; TA[:, pos[po]] = sel[:, 2:3]
    FA = feats_K(sh, SA, TA, EA)
    d = FA.shape[1]
    tg = np.zeros(d, dtype=np.int64); tg[0] = 4 * sh.G       # I(0) = 4|G| (first feature: I at k = 0 is repsPM[0] = 0)
    assert sh.repsPM[0] == 0
    tot, info = count4(FA, FF, tg)
    info.update(korbits=nu)
    return tot, info


def feats_K(sh, S, T, E):
    """(I at repsPM, J at repsPM, Q at repsQ) for arrays (M, N) of s, t, e"""
    q, eps = sh.q, sh.eps
    cols = []
    for k in sh.repsPM:
        sh_ = sh.shift[k]
        cols.append((S * S[:, sh_]).sum(1) + q * (T * T[:, sh_]).sum(1))
    for k in sh.repsPM:
        sh_ = sh.shift[k]
        cols.append(q * (E * E[:, sh_]).sum(1) + (S * E[:, sh_]).sum(1) + (E * S[:, sh_]).sum(1) - (T * T[:, sh_]).sum(1))
    for k in sh.repsQ:
        sh_ = sh.shift[k]
        cols.append((S * T[:, sh_]).sum(1) + eps * (T * S[:, sh_]).sum(1))
    return np.stack(cols, 1).astype(np.int64)
