"""Skew Goethals-Seidel families with a block symmetry of order 2 ("pair-spin").

Let H <= Z_v^* be such that no H-orbit of Z_v is closed under negation (a skew block can be H-invariant), and mu in
Z_v^* with mu not in H and mu^2 in H; K = <H, mu> (it may contain -1). For a block X put T_X(s) = P_X(s) + P_X(mu s);
then P_{X o mu}(s) = P_X(mu s), so a pair of blocks (X, X o mu) contributes T_X, which is invariant under K. Layouts
(A is the skew block, A(0) = 1, A(-x) = -A(x); A o mu is skew as well):

  b: (A, A o mu, C, D)        A skew and H-invariant, C and D K-invariant;           2 + c^2 + d^2 = 4v
  c: (A, A o mu, C, C o mu)   A skew and H-invariant, C H-invariant;                 2 + 2 c^2 = 4v

The condition sum P = 0 off 0 is K-invariant and symmetric, so it holds on the classes of <K, -1>.
A layout with a K-invariant skew block is impossible: a skew block needs -1 outside the image of K in Z_p^* for every
prime p | v (the orbits of the elements of order p), so every such image has odd order and so has K, while
[K : H] = 2. (mu = -1 is allowed: then A o mu = A with the signs off 0 reversed, and T_A = 2 P_A.)

  python3 search/pairspin.py rank v1,v2,..  [--maxO 64] [--samples 20000]
      Gaussian estimate (full covariance, lattice of the sampled values) of the number of solutions of every structure
      -> $HADAMARD_WORK/pairspin_rank.jsonl
  python3 search/pairspin.py sweep v1,v2,.. [--maxA 24] [--maxpairs 3e9] [--maxlist 1e8]
      exact search of every structure with at most 2^maxA skew members and small lists -> pairspin_sweep.jsonl
  python3 search/pairspin.py search v hgens mu layout rows [--maxlist 3e8]
      exact two-list join (layout b: C x D against A; c: A against C) on the classes of <K, -1>;
      every hit is rebuilt over Z_v and checked (verify.py) -> $HADAMARD_WORK/skew_hits.jsonl
"""
import argparse, itertools, json, math, os, sys, time
from math import gcd, log10
import numpy as np
from pathlib import Path as _P
sys.path[:0] = [str(_P(__file__).resolve().parent), str(_P(__file__).resolve().parent.parent / 'verify')]
import workdir as W
import spence_lib as L
from gs_sweep_subgroups import subgroups
import verify


def closure(gens, v):
    S = {1}
    while True:
        T = S | {s * g % v for s in S for g in gens}
        if T == S: return frozenset(S)
        S = T


def order2_multipliers(v, Hs):
    """One mu per coset mu H with mu not in H and mu^2 in H."""
    out, seen = [], set()
    for u in range(2, v):
        if gcd(u, v) != 1 or u in Hs or (u * u) % v not in Hs: continue
        c = frozenset(u * h % v for h in Hs)
        if c in seen: continue
        seen.add(c); out.append(u)
    return out


class PairSpin:
    def __init__(self, v, hgens, mu):
        self.v, self.mu, self.hgens = v, mu, list(hgens)
        G = L.AbGroup([v]); self.G = G
        P = lambda g: G.perm_from_matrix(np.array([[g]]))
        self.stH = L.Structure(G, [P(g) for g in hgens])
        self.stK = L.Structure(G, [P(g) for g in hgens] + [P(mu)])
        self.H, self.K = closure(self.hgens, v), closure(self.hgens + [mu], v)
        stH, stK = self.stH, self.stK
        neg = G.neg()
        def pairs(st):
            rep = np.array([np.flatnonzero(st.orb == o)[0] for o in range(st.O)])
            partner = st.orb[neg[rep]]
            assert all(partner[o] != o for o in range(1, st.O)), 'an orbit is closed under negation'
            return [(o, int(partner[o])) for o in range(1, st.O) if o < partner[o]]
        self.pairsH = pairs(stH)
        # T_X on the K-class c for an H-invariant X: P_X at s_c plus P_X at mu s_c
        self.MH = np.zeros((stK.C, stH.O, stH.O))
        for c, s in enumerate(stK.reps):
            for k in range(2):
                self.MH[c] += stH.S[stH.cls[s * pow(mu, k, v) % v]]
        self.MH = (self.MH + self.MH.transpose(0, 2, 1)) / 2
        self.SK = (stK.S + stK.S.transpose(0, 2, 1)) / 2

    # ---- members -------------------------------------------------------------------------------------------------
    def skew_vectors(self, s, st, prs):
        Y = np.zeros((len(s), st.O)); Y[:, 0] = 1
        for j, (o, q) in enumerate(prs): Y[:, o] = s[:, j]; Y[:, q] = -s[:, j]
        return Y

    def keys_K(self, Y):                        # P_X on the K-classes, X K-invariant
        return np.einsum('bo,cop,bp->bc', Y, self.SK, Y, optimize=True)

    def keys_T(self, Y):                        # T_X on the K-classes, X H-invariant
        return np.einsum('bo,cop,bp->bc', Y, self.MH, Y, optimize=True)

    def roles(self, layout, rows):
        """[(kind, row, pool size, sampler(n, rng) -> sign vectors, key function)] for the independent roles."""
        stH, stK = self.stH, self.stK
        out = []
        if True:
            out.append(('A_H', 1, 2.0 ** len(self.pairsH),
                        lambda n, rng: self.skew_vectors(np.where(rng.random((n, len(self.pairsH))) < .5, 1., -1.), stH, self.pairsH),
                        self.keys_T))
            if layout == 'b':
                c, d = rows
                out.append(('C_K', c, stK.pool(c), lambda n, rng, c=c: stK.sample(c, n, rng), self.keys_K))
                out.append(('D_K', d, stK.pool(d), lambda n, rng, d=d: stK.sample(d, n, rng), self.keys_K))
            else:
                c, = rows
                out.append(('C_H', c, stH.pool(c), lambda n, rng, c=c: stH.sample(c, n, rng), self.keys_T))
        return out

    def sequences(self, layout, members):
        """The four +-1 sequences on Z_v from the role members (orbit sign vectors)."""
        v, mu = self.v, self.mu
        x = lambda st, Y: np.asarray(Y)[st.orb].astype(int)
        perm = (np.arange(v) * mu) % v
        A = x(self.stH, members[0])
        if layout == 'b':
            return [A, A[perm], x(self.stK, members[1]), x(self.stK, members[2])]
        C = x(self.stH, members[1])
        return [A, A[perm], C, C[perm]]


def rows_for(v, layout):
    """Admissible row sums (non-negative) of the free roles."""
    out = []
    if layout == 'b':
        for c in range(1, math.isqrt(4 * v) + 1, 2):
            r = 4 * v - 2 - c * c
            if r < c * c: break
            d = math.isqrt(r)
            if d * d == r and d % 2: out.append((c, d))
    else:
        r = 2 * v - 1; c = math.isqrt(r)
        if c * c == r: out.append((c,))
    return out


def gauss(keys, sizes, nsamp):
    """log10 of the expected number of solutions: product of the pools times the density of sum(keys) at 0, from the
    sampled keys (full covariance on the span, lattice spacing from the sampled values)."""
    mus = [k.mean(0) for k in keys]; covs = [np.cov(k.T) for k in keys]
    mu = sum(mus); Sig = sum(covs)
    w, U = np.linalg.eigh(Sig)
    keep = w > 1e-6 * max(w.max(), 1e-12)
    # coordinates outside the span must already match
    if not keep.all():
        off = U[:, ~keep].T @ mu
        if np.abs(off).max() > 1e-6: return -math.inf, int(keep.sum())
    Uk, wk = U[:, keep], w[keep]
    z = Uk.T @ (0 - mu)
    logd = -0.5 * (keep.sum() * math.log(2 * math.pi) + np.log(wk).sum()) - 0.5 * float((z * z / wk).sum())
    # lattice: each class value moves in steps of g_c (gcd of the differences of the sampled values, summed over roles)
    lat = 0.0
    D = keys[0].shape[1]
    if keep.sum() == D:
        for c in range(D):
            g = 0
            for k in keys:
                u = np.unique(np.rint(k[:, c]).astype(np.int64))
                if len(u) > 1: g = gcd(g, int(np.gcd.reduce(np.diff(u))))
            if g > 1: lat += log10(g)
    return sum(log10(max(s, 1e-300)) for s in sizes) + logd / math.log(10) + lat, int(keep.sum())


def rank(vs, maxO, samples, out):
    rng = np.random.default_rng(1)
    done = set()
    if out.exists():
        for l in open(out):
            r = json.loads(l); done.add((r['v'], tuple(r['hgens']), r['mu']))
    for v in vs:
        for Hs, gens in subgroups(v).items():
            if (v - 1) in Hs: continue
            if (v - 1) // len(Hs) > 4 * maxO: continue
            for mu in order2_multipliers(v, Hs):
                if (v, tuple(gens), mu) in done: continue
                try:
                    ps = PairSpin(v, list(gens), mu)
                except (AssertionError, ValueError):
                    continue
                if ps.stH.O > maxO: continue
                base = dict(v=v, hgens=list(gens), H=len(ps.H), mu=mu, K=len(ps.K), OH=ps.stH.O, OK=ps.stK.O, classes=ps.stK.C)
                for layout in 'bc':
                    for rows in rows_for(v, layout):
                        R = ps.roles(layout, rows)
                        sizes = [r[2] for r in R]
                        if min(sizes) <= 0: continue
                        keys = [r[4](r[3](samples, rng)) for r in R]
                        lE, rk = gauss(keys, sizes, samples)
                        rec = dict(base, layout=layout, rows=list(rows), logE=round(lE, 2), rank=rk,
                                   lists=[round(log10(s), 1) for s in sizes])
                        print(json.dumps(rec), flush=True)
                        with open(out, 'a') as f: f.write(json.dumps(rec) + '\n')


# ---- exact search ---------------------------------------------------------------------------------------------------
def enumerate_row(st, r, maxn):
    """All sign vectors on the orbits of st with row sum r (orbit 0 included), as an int8 array; None if more than maxn."""
    n = st.pool(r)
    if n > maxn: return None
    O, sizes = st.O, st.sizes
    out = []
    # meet in the middle over the orbits: split into two halves, combine the partial sums
    h = O // 2
    def half(idx):
        res = {}
        for signs in itertools.product((1, -1), repeat=len(idx)):
            s = sum(sg * int(sizes[i]) for sg, i in zip(signs, idx))
            res.setdefault(s, []).append(signs)
        return res
    L1, L2 = half(list(range(h))), half(list(range(h, O)))
    for s1, a1 in L1.items():
        a2 = L2.get(r - s1)
        if not a2: continue
        A1 = np.array(a1, np.int8); A2 = np.array(a2, np.int8)
        out.append(np.concatenate([np.repeat(A1, len(A2), 0), np.tile(A2, (len(A1), 1))], 1))
    return np.concatenate(out) if out else np.zeros((0, O), np.int8)


def keys_int(fn, Y, chunk=1 << 18):
    out = np.empty((len(Y), fn.__self__.stK.C), np.int32)
    for i in range(0, len(Y), chunk):
        out[i:i + chunk] = np.rint(fn(Y[i:i + chunk].astype(np.float64))).astype(np.int32)
    return out


def hash_rows(K, salt):
    """64-bit hashes of the integer rows of K."""
    h = np.full(len(K), np.uint64(salt))
    for c in range(K.shape[1]):
        h = (h ^ (K[:, c].astype(np.int64) + 1 << 20).astype(np.uint64)) * np.uint64(0x9E3779B97F4A7C15)
        h ^= h >> np.uint64(29)
    return h


def join(left, right, D, maxhash, log, chunk=4e6):
    """Pairs (i, j) with sum(left keys)[i] + right[j] = 0, where left is one or two (keys, labels) lists combined as
    all sums and right is one list. The smaller side is hashed, the other streamed in chunks."""
    salt = 0x1234567
    nleft = int(np.prod([len(k) for k, _ in left]))
    out = []
    def check(sl, sr, il, ir):          # candidate (left sum rows, right rows)
        ok = np.all(sl + sr == 0, axis=1)
        return il[ok], ir[ok]
    def left_rows(idx):                 # left sums for flat indices idx
        if len(left) == 1: return left[0][0][idx], (idx,)
        n2 = len(left[1][0]); i1, i2 = idx // n2, idx % n2
        return left[0][0][i1] + left[1][0][i2], (i1, i2)
    KR = right[0]
    if nleft <= len(KR) and nleft <= maxhash:
        idx = np.arange(nleft, dtype=np.int64)
        S, _ = left_rows(idx)
        hs = hash_rows(S, salt); order = np.argsort(hs); hs = hs[order]
        step = max(1, int(chunk))
        for r0 in range(0, len(KR), step):
            R = KR[r0:r0 + step]; hr = hash_rows(-R, salt)
            pos = np.searchsorted(hs, hr); pos = np.minimum(pos, len(hs) - 1)
            m = np.flatnonzero(hs[pos] == hr)
            for j in m:
                p = pos[j]
                while p < len(hs) and hs[p] == hr[j]:
                    i = order[p]; Si, sub = left_rows(np.array([i]))
                    if np.array_equal(Si[0], -R[j]): out.append((int(i), int(r0 + j)))
                    p += 1
    elif len(KR) <= maxhash:
        hs = hash_rows(-KR, salt); order = np.argsort(hs); hs = hs[order]
        step = max(1, int(chunk))
        for l0 in range(0, nleft, step):
            idx = np.arange(l0, min(nleft, l0 + step), dtype=np.int64)
            S, _ = left_rows(idx); hl = hash_rows(S, salt)
            pos = np.searchsorted(hs, hl); pos = np.minimum(pos, len(hs) - 1)
            m = np.flatnonzero(hs[pos] == hl)
            for i in m:
                p = pos[i]
                while p < len(hs) and hs[p] == hl[i]:
                    j = order[p]
                    if np.array_equal(S[i], -KR[j]): out.append((int(idx[i]), int(j)))
                    p += 1
    else:
        raise RuntimeError(f'both sides too large to hash ({nleft:.3g}, {len(KR):.3g})')
    return out


def search(v, hgens, mu, layout, rows, maxlist, log, maxhash=2e8, maxpairs=2e10):
    ps = PairSpin(v, hgens, mu)
    stH, stK = ps.stH, ps.stK
    log(f'v={v} H={sorted(ps.H)[:6]}.. |H|={len(ps.H)} mu={mu} |K|={len(ps.K)} OH={stH.O} OK={stK.O} classes={stK.C} '
        f'layout {layout} rows {rows}')
    sA = np.array(list(itertools.product((1., -1.), repeat=len(ps.pairsH))), np.int8) if len(ps.pairsH) <= 26 else None
    if sA is None: raise RuntimeError('skew list too large')
    YA = ps.skew_vectors(sA.astype(float), stH, ps.pairsH).astype(np.int8)
    KA = keys_int(ps.keys_T, YA)
    if layout == 'b':
        c, d = rows
        YC = enumerate_row(stK, c, maxlist); YD = enumerate_row(stK, d, maxlist)
        if YC is None or YD is None: raise RuntimeError('list too large')
        if len(YC) * len(YD) > maxpairs: raise RuntimeError(f'{len(YC) * len(YD):.3g} pairs')
        KC = keys_int(ps.keys_K, YC); KD = keys_int(ps.keys_K, YD)
        log(f'lists A {len(KA)}, C {len(KC)}, D {len(KD)}')
        pairs = join([(KC, None), (KD, None)], (KA, None), stK.C, maxhash, log)
        nD = len(KD)
        hits = [[YA[j], YC[i // nD], YD[i % nD]] for i, j in pairs]
    else:
        c, = rows
        YC = enumerate_row(stH, c, maxlist)
        if YC is None: raise RuntimeError('list too large')
        KC = keys_int(ps.keys_T, YC)
        log(f'lists A {len(KA)}, C {len(KC)}')
        pairs = join([(KA, None)], (KC, None), stK.C, maxhash, log)
        hits = [[YA[i], YC[j]] for i, j in pairs]
    log(f'{len(hits)} key matches')
    nver = 0
    for mem in hits:
        seqs = ps.sequences(layout, mem)
        tot = verify.paf_sum(seqs)
        good = tot[0] == 4 * v and not tot[1:].any()
        if not good: continue
        Hm = verify.goethals_seidel(*map(verify.circulant, seqs)).astype(np.float64)
        had = np.array_equal(Hm @ Hm.T, 4 * v * np.eye(4 * v)); skew = np.array_equal(Hm + Hm.T, 2 * np.eye(4 * v))
        log(f'hit: complementary={good} Hadamard={had} skew={skew} rows {[int(s.sum()) for s in seqs]}')
        if had and skew:
            nver += 1
            with open(W.hits(), 'a') as f:
                f.write(json.dumps(dict(n=v, gens=sorted(hgens), pair_multiplier=mu, layout=layout,
                                        pattern=[int(s.sum()) for s in seqs], hadamard=True, skew=True,
                                        seqs=[''.join('+' if x > 0 else '-' for x in s) for s in seqs])) + '\n')
            if nver >= 3: break
    log(f'{nver} verified solutions')
    return nver, len(hits)


def sweep(vs, maxA, maxpairs, maxlist, out, log):
    """Exact search of every structure whose lists are small enough (no estimate needed)."""
    done = set()
    if out.exists():
        for l in open(out):
            r = json.loads(l); done.add((r['v'], tuple(r['hgens']), r['mu'], r['layout'], tuple(r['rows'])))
    for v in vs:
        found = False
        for Hs, gens in sorted(subgroups(v).items(), key=lambda t: -len(t[0])):
            if found: break
            if (v - 1) in Hs: continue
            if (v - 1) // len(Hs) > 4 * 64: continue
            for mu in order2_multipliers(v, Hs):
                if found: break
                try: ps = PairSpin(v, list(gens), mu)
                except (AssertionError, ValueError): continue
                if len(ps.pairsH) > maxA: continue
                for layout in 'bc':
                    for rows in rows_for(v, layout):
                        key = (v, tuple(gens), mu, layout, tuple(rows))
                        if key in done or found: continue
                        if layout == 'b':
                            nC, nD = ps.stK.pool(rows[0]), ps.stK.pool(rows[1])
                            if nC == 0 or nD == 0 or nC * nD > maxpairs or max(nC, nD) > maxlist: continue
                        else:
                            nC = ps.stH.pool(rows[0])
                            if nC == 0 or nC > maxlist: continue
                        t0 = time.time()
                        rec = dict(v=v, hgens=list(gens), H=len(ps.H), mu=mu, K=len(ps.K), OH=ps.stH.O, OK=ps.stK.O,
                                   classes=ps.stK.C, layout=layout, rows=list(rows), A=2 ** len(ps.pairsH))
                        try:
                            nver, nm = search(v, list(gens), mu, layout, rows, maxlist, lambda m: None)
                            rec.update(matches=nm, verified=nver)
                            if nver: found = True
                        except (RuntimeError, MemoryError) as ex:
                            rec.update(error=repr(ex)[:80])
                        rec['sec'] = round(time.time() - t0, 1)
                        log(json.dumps(rec))
                        with open(out, 'a') as f: f.write(json.dumps(rec) + '\n')


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest='cmd', required=True)
    r = sub.add_parser('rank'); r.add_argument('vs'); r.add_argument('--maxO', type=int, default=64)
    r.add_argument('--samples', type=int, default=20000); r.add_argument('--out', default='pairspin_rank.jsonl')
    s = sub.add_parser('search'); s.add_argument('v', type=int); s.add_argument('hgens'); s.add_argument('mu', type=int)
    s.add_argument('layout', choices=('b', 'c')); s.add_argument('rows')
    s.add_argument('--maxlist', type=float, default=3e8)
    w = sub.add_parser('sweep'); w.add_argument('vs'); w.add_argument('--maxA', type=int, default=24)
    w.add_argument('--maxpairs', type=float, default=3e9); w.add_argument('--maxlist', type=float, default=1e8)
    w.add_argument('--out', default='pairspin_sweep.jsonl')
    a = ap.parse_args()
    t0 = time.time(); log = lambda m: print(f'{m} [{time.time() - t0:.0f}s]', flush=True)
    if a.cmd == 'rank':
        rank([int(x) for x in a.vs.split(',')], a.maxO, a.samples, W.WORK / a.out)
    elif a.cmd == 'sweep':
        sweep([int(x) for x in a.vs.split(',')], a.maxA, a.maxpairs, a.maxlist, W.WORK / a.out, lambda m: print(m, flush=True))
    else:
        search(a.v, [int(x) for x in a.hgens.split(',')], a.mu, a.layout, tuple(int(x) for x in a.rows.split(',')), a.maxlist, log)


if __name__ == '__main__':
    main()
