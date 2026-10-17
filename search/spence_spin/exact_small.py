"""Exact uncapped joins for the levels whose lists are enumerable (Z2^2, Z11, and A-side / B-side of larger levels)."""
import numpy as np, sys, time
from collections import defaultdict
from struct573 import build
from levels import make

def enum_list(vals, w, rowsum, maxn=60_000_000):
    """all vectors (a_R in vals[R]) with sum_R w_R a_R = rowsum.  Returns (V (n x R) int64, weight (n) float = product of multiplicities)."""
    R = len(vals)
    arrs = [np.array(list(vals[r].keys()), dtype=np.int64) for r in range(R)]
    mult = [np.array(list(vals[r].values()), dtype=np.float64) for r in range(R)]
    # suffix ranges for pruning
    lo = np.zeros(R + 1, dtype=np.int64); hi = np.zeros(R + 1, dtype=np.int64)
    for r in range(R - 1, -1, -1):
        lo[r] = lo[r + 1] + w[r] * arrs[r].min(); hi[r] = hi[r + 1] + w[r] * arrs[r].max()
    V = np.zeros((1, 0), dtype=np.int64); W = np.ones(1); partial = np.zeros(1, dtype=np.int64)
    for r in range(R):
        idx_V = np.repeat(np.arange(len(V)), len(arrs[r]))
        newv = np.tile(arrs[r], len(V)); neww = np.tile(mult[r], len(V))
        p2 = partial[idx_V] + w[r] * newv
        ok = (p2 + lo[r + 1] <= rowsum) & (p2 + hi[r + 1] >= rowsum)
        V = np.concatenate([V[idx_V][ok], newv[ok][:, None]], axis=1); W = W[idx_V][ok] * neww[ok]; partial = p2[ok]
        if len(V) > maxn: raise MemoryError('list too large: %d' % len(V))
    assert np.all(partial == rowsum)
    return V, W

def keys(V, M):
    """class quadratic forms:  key_c = V^T M_c V  (M_c need not be symmetric)"""
    return np.stack([np.einsum('ni,ij,nj->n', V, M[c], V) for c in range(M.shape[0])], axis=1)

def exact_join(L, verbose=True):
    t = time.time()
    VA, WA = enum_list(L.valsA, L.wA, 2)
    VB, WB = enum_list(L.valsB, L.wB, 0)
    KA = keys(VA, L.MA); KB = keys(VB, L.MB)
    nA = (VA * VA) @ L.wA; nB = (VB * VB) @ L.wB
    # norm equation  |a|^2 + 3|b|^2 = 4n - 4|N|   (implied by the class equations and row sums; used as a consistency check)
    tgtKB = L.target - KA                                  # needed key of B
    dic = defaultdict(list)
    for i, row in enumerate(KB): dic[row.tobytes()].append(i)
    pairs = 0; wsum = 0.0; usedA = set(); usedB = set(); ex = None
    for i, row in enumerate(tgtKB):
        h = dic.get(row.tobytes())
        if h:
            okh = [j for j in h if nA[i] + 3 * nB[j] == L.norm]
            assert len(okh) == len(h), 'norm equation not implied?'
            pairs += len(h); wsum += WA[i] * sum(WB[j] for j in h); usedA.add(i); usedB.update(h)
            if ex is None: ex = (VA[i].copy(), VB[h[0]].copy())
    if verbose:
        print('%s: |A|=%d |B|=%d  vector pairs=%d  weighted(sign-pattern) pairs=%.4g  usedA=%d usedB=%d  %.1fs' % (L.name, len(VA), len(VB), pairs, wsum, len(usedA), len(usedB), time.time() - t))
    return dict(name=L.name, nA=len(VA), nB=len(VB), pairs=pairs, weighted=wsum, usedA=len(usedA), usedB=len(usedB), example=ex)

if __name__ == '__main__':
    for u in (3, 9):
        S = build(u)
        for name in sys.argv[1:] or ['Z2^2', 'Z11']:
            L = make(S, name)
            print('u=%d' % u, end=' ')
            exact_join(L)
