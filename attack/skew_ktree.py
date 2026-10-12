"""Skew Goethals-Seidel search with the four-list k-tree matcher (match4.c) for structures whose pools
are too large for the exact matcher find4b. A is skew-type (A(0)=1, A(-x)=-A(x)), B, C, D are constant
on the H-orbits; every hit is checked: H H^T = 4n I and H + H^T = 2 I.

  python3 skew_ktree.py n g[,g2] [--maxtau T] [--top M] [--pattern b,c,d] [--sols K] [--maxpool L]

Pools larger than L are replaced by L distinct uniform samples (the expected yield scales by L/pool).
"""
import argparse, json, math, os, subprocess, sys, tempfile, time
import numpy as np
sys.path.insert(0, '/home/user/Hadamard/attack'); sys.path.insert(0, '/home/user/Hadamard/verify')
import spence_lib as L
import verify
HERE = '/home/user/Hadamard/attack'


def pool_exact(st, r):
    """All +-1 orbit vectors with row sum r, by splitting the orbits into two halves."""
    O = st.O; h = O // 2
    def half(lo, hi):
        m = hi - lo
        idx = np.arange(1 << m, dtype=np.int64)
        Y = (1 - 2 * ((idx[:, None] >> np.arange(m)) & 1)).astype(np.int8)
        return Y, Y.astype(np.int64) @ st.sizes[lo:hi]
    Y1, s1 = half(0, h); Y2, s2 = half(h, O)
    order = np.argsort(s2, kind='stable'); s2s = s2[order]
    out = []
    for v in np.unique(s1):
        a = np.flatnonzero(s1 == v)
        lo, hi = np.searchsorted(s2s, r - v), np.searchsorted(s2s, r - v, side='right')
        if hi == lo: continue
        b = order[lo:hi]
        out.append(np.concatenate([np.repeat(Y1[a], len(b), 0), np.tile(Y2[b], (len(a), 1))], 1))
    return np.concatenate(out) if out else np.zeros((0, O), np.int8)


def values(st, Y, pm):
    """Class values (PAF - n mod 4)/4 as int8, and the mask of rows whose values fit in int8."""
    chunk = max(1000, int(1.5e8 // (st.C * st.O)))
    V = np.empty((len(Y), st.C), np.int16)
    for s in range(0, len(Y), chunk):
        V[s:s + chunk] = np.rint((st.class_values(Y[s:s + chunk].astype(np.float64)) - pm) / 4)
    ok = np.abs(V).max(1) <= 127
    return V.astype(np.int8), ok


def pool(st, r, maxpool, rng):
    """All sign vectors with row sum r if there are at most maxpool of them (and the orbits are few enough
    to enumerate), else maxpool distinct uniform samples."""
    N = st.pool(r)
    if N <= maxpool and st.O <= 44:
        return pool_exact(st, r)
    Y = np.unique(st.sample(r, int(min(N, maxpool) * 1.02) + 16, rng).astype(np.int8), axis=0)
    return Y[rng.permutation(len(Y))[:maxpool]]


def write(path, V):
    rec = np.zeros((len(V), 8 + V.shape[1]), np.uint8)
    rec[:, :8] = np.arange(len(V), dtype=np.uint64).view(np.uint8).reshape(-1, 8); rec[:, 8:] = V.view(np.uint8)
    rec.tofile(path)


def clamp(Vl, Vr, cap=2.5e7):
    """Smallest S such that the most likely left and right pair lists per tau stay below cap. The yield per
    unit of work grows like the smaller pair list, so the clamp is as loose as memory allows."""
    S = 0
    for (X, Y) in (Vl, Vr):
        lp, s = math.log(len(X)) + math.log(len(Y)), 0
        while lp > math.log(cap) and s < X.shape[1]:
            hx = np.bincount(X[:, s].astype(np.int64) + 128, minlength=256) / len(X)
            hy = np.bincount(Y[:, s].astype(np.int64) + 128, minlength=256) / len(Y)
            lp += math.log(np.convolve(hx, hy).max()); s += 1
        S = max(S, s)
    return S


def run(n, gens, maxtau=200, top=6, only=None, sols=4, clampS=None, maxpool=4e6, seed=1, verbose=True):
    G = L.AbGroup([n]); st = L.Structure(G, [G.perm_from_matrix(np.array([[g]])) for g in gens])
    rng = np.random.default_rng(seed); maxpool = int(maxpool)
    rep = np.array([np.flatnonzero(st.orb == o)[0] for o in range(st.O)])
    partner = st.orb[G.neg()[rep]]
    assert not np.any(partner[1:] == np.arange(1, st.O)), "-1 fixes an orbit"
    pairs = [(o, partner[o]) for o in range(1, st.O) if o < partner[o]]
    if len(pairs) <= 22:
        s = 1 - 2 * ((np.arange(2 ** len(pairs))[:, None] >> np.arange(len(pairs))) & 1)
    else:
        s = np.unique(rng.choice([-1, 1], size=(maxpool, len(pairs))), axis=0)
    YA = np.zeros((len(s), st.O), dtype=np.int8); YA[:, 0] = 1
    for j, (o, q) in enumerate(pairs): YA[:, o] = s[:, j]; YA[:, q] = -s[:, j]
    pm = n % 4
    VA, ok = values(st, YA, pm); YA, VA = YA[ok], VA[ok]
    lim = int(math.isqrt(4 * n)); pats = []
    for b in range(1, lim + 1, 2):
        for c in range(b, lim + 1, 2):
            d2 = 4 * n - 1 - b * b - c * c
            if d2 < c * c: continue
            d = math.isqrt(d2)
            if d * d == d2 and all(st.pool(r) > 0 for r in (b, c, d)): pats.append((b, c, d))
    if only: pats = [p for p in pats if p == tuple(only)]
    if verbose: print(f"n={n} gens={gens} orbits={st.O} classes={st.C} skew pool={len(YA)} patterns={pats}", flush=True)
    found = []
    for q in pats:
        t0 = time.time()
        pools, vals = {}, {}
        for r in set(q):
            Y = pool(st, r, maxpool, rng); V, ok = values(st, Y, pm)
            pools[r], vals[r] = Y[ok], V[ok]
        big = max(range(3), key=lambda k: len(pools[q[k]]))
        rest = [q[k] for k in range(3) if k != big]
        roles = [(YA, VA), (pools[rest[0]], vals[rest[0]]), (pools[rest[1]], vals[rest[1]]), (pools[q[big]], vals[q[big]])]
        sizes = [len(r[0]) for r in roles]
        S = clampS or clamp((roles[0][1], roles[3][1]), (roles[1][1], roles[2][1]))
        S = max(1, min(S, st.C - 2))
        if verbose: print(f"  pattern {q}: full pools {[st.pool(r) for r in q]:}")
        if verbose: print(f"  pattern {q}: roles A,{rest[0]},{rest[1]},{q[big]} sizes {sizes} clamp S={S} [{time.time()-t0:.0f}s prep]", flush=True)
        tmp = tempfile.mkdtemp(dir=HERE)
        paths = [os.path.join(tmp, f"{k}.bin") for k in range(4)]
        for p, (Y, V) in zip(paths, roles): write(p, V)
        cmd = [os.path.join(HERE, 'match4'), '-s', str(S), '-t', str(-pm), '-n', str(maxtau), '-k', str(sols), '-m', str(top)] + paths + [str(st.C)]
        res = subprocess.run(cmd, capture_output=True, text=True)
        for p in paths: os.remove(p)
        os.rmdir(tmp)
        err = res.stderr.strip().splitlines()
        if verbose:
            for line in err[:3] + err[-2:]: print("   ", line, flush=True)
        for line in res.stdout.splitlines():
            if not line.startswith('SOL'): continue
            ids = [int(x, 16) for x in line.split()[1:]]
            seqs = [roles[k][0][ids[k]].astype(np.int64)[st.orb] for k in range(4)]
            rs = [int(x.sum()) for x in seqs]
            tot = verify.paf_sum(seqs)
            if not (tot[0] == 4 * n and not tot[1:].any()): continue
            H = verify.goethals_seidel(*map(verify.circulant, seqs)).astype(np.float64)   # BLAS; exact below 2^53
            had = np.array_equal(H @ H.T, 4 * n * np.eye(4 * n))
            skew = np.array_equal(H + H.T, 2 * np.eye(4 * n))
            found.append(dict(n=n, gens=list(gens), pattern=rs, hadamard=bool(had), skew=bool(skew),
                              seqs=[''.join('+' if x > 0 else '-' for x in sq) for sq in seqs]))
            if verbose: print(f"    VERIFIED rows={rs}: Hadamard={had} skew={skew}", flush=True)
            if len(found) >= sols: break
        if found: break
    return found


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('n', type=int); ap.add_argument('gens')
    ap.add_argument('--maxtau', type=int, default=200); ap.add_argument('--top', type=int, default=6)
    ap.add_argument('--pattern'); ap.add_argument('--sols', type=int, default=4); ap.add_argument('--clamp', type=int)
    ap.add_argument('--maxpool', type=float, default=4e6); ap.add_argument('--seed', type=int, default=1)
    a = ap.parse_args()
    found = run(a.n, [int(g) for g in a.gens.split(',')], a.maxtau, a.top,
                [int(x) for x in a.pattern.split(',')] if a.pattern else None, a.sols, a.clamp, a.maxpool, a.seed)
    if found:
        with open(os.path.join(HERE, 'skew_hits.json'), 'a') as f:
            for h in found: f.write(json.dumps(h) + '\n')
    print(f"done: {len(found)} verified", flush=True)
