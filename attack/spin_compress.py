"""Exact compression test for the spin structure of spin_skew.py.

A is invariant under K = <H, mu> and of skew type, B is H-invariant with row sum b, C = B o mu, D = B o mu^2.
For d | v, compress to Z_d (x_j = sum of X(t) over t = j mod d). The compressions satisfy, as integers,

    PAF(a)(k) + PAF(b)(k) + PAF(b)(mu k) + PAF(b)(mu^2 k) = 4v [k = 0]   for k in Z_d,

because PAF(B o mu)(k) = PAF(B)(mu k). The compression a of A is listed exactly (a DP over the K-orbit pairs of
Z_v). The compression b of B is constant on the H-orbits of Z_d; with beta_i its value on orbit i, every
beta_i = v/d (mod 2), |beta_i| <= v/d, sum_i |O_i| beta_i = b, and at k = 0, sum_i |O_i| beta_i^2 =
(4v - |a|^2)/3. All integer vectors beta meeting these conditions are enumerated (a superset of the attainable
ones), and the equations at k != 0 are checked. If no pair (a, beta) survives, the structure is ruled out.

  python3 attack/spin_compress.py v mu b d [H generators ...]
"""
import sys
import numpy as np


def closure(gens, v):
    H = {1}
    while True:
        H2 = H | {h * g % v for h in H for g in gens}
        if H2 == H: return sorted(H)
        H = H2


def orbits(n, G):
    oid = -np.ones(n, dtype=np.int64); reps = []
    for x in range(n):
        if oid[x] < 0:
            for g in G: oid[x * g % n] = len(reps)
            reps.append(x)
    return oid, reps


def dedup(K):
    vt = np.dtype((np.void, K.dtype.itemsize * K.shape[1]))
    _, idx = np.unique(np.ascontiguousarray(K).view(vt).ravel(), return_index=True)
    return K[idx]


def skew_compressions(v, d, K):
    oid, reps = orbits(v, K)
    O = len(reps); t = np.arange(v)
    C = np.zeros((O, d), dtype=np.int16)
    np.add.at(C, (oid, t % d), 1)
    partner = oid[(-np.array(reps)) % v]
    assert not any(partner[o] == o for o in range(1, O)), "-1 fixes an orbit"
    S = C[0][None, :].copy()
    for o in range(1, O):
        if o < partner[o]:
            it = C[o] - C[partner[o]]
            S = dedup(np.concatenate([S + it, S - it]))
    return S


def betas(sizes, m, rowsum, norm):
    """Integer vectors beta (one value per H-orbit of Z_d) with beta_i = m (mod 2), |beta_i| <= m,
    sum sizes_i beta_i = rowsum and sum sizes_i beta_i^2 = norm."""
    n = len(sizes); out = []
    vals = [x for x in range(-m, m + 1) if (x - m) % 2 == 0]
    suf_min_sq = [0] * (n + 1)
    for i in range(n - 1, -1, -1): suf_min_sq[i] = suf_min_sq[i + 1] + sizes[i] * min(x * x for x in vals)
    cur = []
    def rec(i, s, q):
        if i == n:
            if s == rowsum and q == norm: out.append(list(cur))
            return
        rest_max = sum(sizes[i:]) * m
        if abs(rowsum - s) > rest_max: return
        for x in vals:
            q2 = q + sizes[i] * x * x
            if q2 + suf_min_sq[i + 1] > norm: continue
            cur.append(x); rec(i + 1, s + sizes[i] * x, q2); cur.pop()
    rec(0, 0, 0)
    return np.array(out, dtype=np.int64).reshape(-1, n)


def pafs(X):
    d = X.shape[1]
    return np.stack([(X * np.roll(X, -k, axis=1)).sum(1) for k in range(d)], axis=1)


def test(v, mu, b, d, hgens):
    H = closure(hgens, v); K = closure(list(hgens) + [mu], v)
    SA = skew_compressions(v, d, K).astype(np.int64)
    oidd, repsd = orbits(d, sorted({h % d for h in H}))
    sizes = np.bincount(oidd)
    k = np.arange(d)
    total, hits = 0, 0
    PA = pafs(SA)
    for norm in sorted(set(((4 * v - PA[:, 0]) // 3).tolist())):
        rows = [i for i in range(len(SA)) if 4 * v - PA[i, 0] == 3 * norm]
        need = {tuple(-PA[i]) for i in rows}            # spun PAF of b must equal -PAF(a) off zero
        for sgn in ((1, -1) if b else (1,)):
            Bt = betas(list(sizes), v // d, sgn * b, norm)
            total += len(Bt)
            if not len(Bt): continue
            X = Bt[:, oidd]
            PB = pafs(X)
            spun = PB + PB[:, (mu * k) % d] + PB[:, (mu * mu * k) % d]
            spun[:, 0] -= 4 * v
            hits += sum(1 for r in spun if tuple(r) in need)
    return len(SA), total, hits


if __name__ == '__main__':
    v, mu, b, d = map(int, sys.argv[1:5]); hg = [int(x) for x in sys.argv[5:]]
    nA, nB, hits = test(v, mu, b, d, hg)
    print(f"v={v} mu={mu} b={b} d={d} H={hg}: {nA} A-compressions, {nB} candidate B-compressions, {hits} consistent pairs"
          + ("  => ruled out" if hits == 0 else ""), flush=True)
