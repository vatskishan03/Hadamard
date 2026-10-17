"""Exact count of Spence solutions for one structure with moderate pools (chunked hash join)."""
import numpy as np

def pool_vectors(st, r):
    """All sign vectors with row sum r (vectorized over 2^O)."""
    O = st.O
    out = []
    step = 1 << 20
    for start in range(0, 1 << O, step):
        idx = np.arange(start, min(start + step, 1 << O), dtype=np.int64)
        bits = ((idx[:, None] >> np.arange(O)) & 1).astype(np.int8)
        Y = 1 - 2 * bits
        keep = (Y.astype(np.int64) @ st.sizes) == r
        out.append(Y[keep])
    return np.concatenate(out).astype(np.float64)

def count(st, chunk=2000):
    v = st.G.v; pm = v % 4; target = (-4 - 4 * pm) // 4
    A, B = pool_vectors(st, 2), pool_vectors(st, 0)
    VA = np.rint((st.class_values(A) - pm) / 4).astype(np.int64)
    VB = np.rint((st.class_values(B) - pm) / 4).astype(np.int64)
    R = np.random.default_rng(11).integers(1, 2**62, size=st.C, dtype=np.int64).astype(np.uint64) | np.uint64(1)
    hA = (VA.astype(np.uint64) * R).sum(1, dtype=np.uint64)
    hB = (VB.astype(np.uint64) * R).sum(1, dtype=np.uint64)
    ht = (np.full(st.C, target, dtype=np.int64).astype(np.uint64) * R).sum(dtype=np.uint64)
    right = (ht - (hB[:, None] + hB[None, :])).ravel()
    right.sort()
    total = 0
    for s in range(0, len(hA), chunk):
        left = (hA[s:s + chunk, None] + hB[None, :]).ravel()
        total += int((np.searchsorted(right, left, 'right') - np.searchsorted(right, left, 'left')).sum())
    return total, len(A), len(B)
