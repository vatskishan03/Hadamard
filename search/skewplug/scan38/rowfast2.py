"""vectorised exact row-sum feasibility of a slot (same conditions as sk_pat.row_patterns_skew; returns the number of raw (sigma, eps) class triples of B, C, D passing, no sign-gauge / symmetry normalisation)"""
import numpy as np
from collections import defaultdict
from gsk import block_classes
def _classes(sh, n, m):
    # The cache lives on the shape object.  Correction (search/NOTES.md section 91): an earlier version kept module-level caches keyed by id(sh); id() is reused after a shape object is freed, so a worker
    # that ran several shapes in one process used the block classes of another shape and flagged feasible slots as having no row-sum pattern (575, 627, 957, 385 over Z55 x F7, ...).
    cc = sh.__dict__.setdefault('_rf2_cc', {})
    k = (n, m)
    if k not in cc:
        if len(cc) > 4000: cc.clear()
        cc[k] = sorted(block_classes(sh, n, m))
    return cc[k]
def row_count_fast(sh, al, ms):
    r, N = sh.r, sh.N; T2 = sh.T2
    cls = [None] + [_classes(sh, al[i], ms[i]) for i in (1, 2, 3)]
    if any(len(c) == 0 for c in cls[1:]): return 0
    Wn = [N - al[i] for i in range(4)]
    taus = [list(range(-w, w + 1, 2)) for w in Wn]
    Rsq = sh.Rsq - 1
    S01 = {a * a + b * b for a in taus[0] for b in taus[1]}; S23 = {a * a + b * b for a in taus[2] for b in taus[3]}
    def tau_ok(rem):
        if rem < 0 or rem % r: return False
        v = rem // r
        return any((v - x) in S23 for x in S01)
    cB = np.array(cls[1], dtype=np.int64); cC = np.array(cls[2], dtype=np.int64); cD = np.array(cls[3], dtype=np.int64)
    RC = cC[:, 0] + r * cC[:, 1]; RD = cD[:, 0] + r * cD[:, 1]
    sC = cC[:, 0] ** 2; sD = cD[:, 0] ** 2
    R2 = (RC[:, None] ** 2 + RD[None, :] ** 2).ravel()
    S2 = (sC[:, None] + sD[None, :]).ravel()
    order = np.argsort(R2, kind='stable'); sR = R2[order]; sS = S2[order]
    cache = {}
    cnt = 0
    for (sg, ep) in cB:
        RB = sg + r * ep
        need = Rsq - RB * RB
        lo = np.searchsorted(sR, need, 'left'); hi = np.searchsorted(sR, need, 'right')
        if hi <= lo: continue
        s2v = 1 + sg * sg + sS[lo:hi]
        vals, counts = np.unique(s2v, return_counts=True)
        for u, c in zip(vals, counts):
            u = int(u)
            if u not in cache: cache[u] = tau_ok(T2 - u)
            if cache[u]: cnt += int(c)
    return cnt


def row_feasible(sh, al, ms):
    """True iff some (sigma, eps) class triple of B, C, D passes the exact row-sum conditions (same as row_count_fast(...) > 0), early exit"""
    r, N = sh.r, sh.N; T2 = sh.T2
    cls = [None] + [_classes(sh, al[i], ms[i]) for i in (1, 2, 3)]
    if any(len(c) == 0 for c in cls[1:]): return False
    Wn = [N - al[i] for i in range(4)]
    _ok = sh.__dict__.setdefault('_rf2_ok', {})
    key = tuple(Wn)
    if key not in _ok:
        taus = [list(range(-w, w + 1, 2)) for w in Wn]
        S01 = {a * a + b * b for a in taus[0] for b in taus[1]}; S23 = {a * a + b * b for a in taus[2] for b in taus[3]}
        smax = 1 + 3 * N * N + 10
        arr = np.zeros(smax + 1, dtype=bool)
        # tau_ok(T2 - s2) for every s2 that can occur; s2 = 1 + sB^2 + sC^2 + sD^2 <= 1 + 3 N^2
        for s2 in range(0, smax + 1):
            rem = T2 - s2
            if rem < 0: break
            if rem % r: continue
            v = rem // r
            arr[s2] = any((v - x) in S23 for x in S01)
        if len(_ok) > 200: _ok.clear()
        _ok[key] = arr
    arr = _ok[key]
    Rsq = sh.Rsq - 1
    cB = np.array(cls[1], dtype=np.int64); cC = np.array(cls[2], dtype=np.int64); cD = np.array(cls[3], dtype=np.int64)
    RC = cC[:, 0] + r * cC[:, 1]; RD = cD[:, 0] + r * cD[:, 1]
    R2 = (RC[:, None] ** 2 + RD[None, :] ** 2).ravel()
    S2 = (cC[:, 0][:, None] ** 2 + cD[:, 0][None, :] ** 2).ravel()
    order = np.argsort(R2, kind='stable'); sR = R2[order]; sS = S2[order]
    for (sg, ep) in cB:
        RB = sg + r * ep; need = Rsq - RB * RB
        lo = np.searchsorted(sR, need, 'left'); hi = np.searchsorted(sR, need, 'right')
        if hi <= lo: continue
        idx = 1 + sg * sg + sS[lo:hi]
        idx = idx[idx < len(arr)]
        if idx.size and arr[idx].any(): return True
    return False
