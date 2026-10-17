"""Propus search with compressions: symmetric Hadamard matrices (A, B, B, D) over Z_v for composite v and larger orbit
structures than propus.py can join directly.

For d | v the compression x'(j) = sum_{t = j (mod d)} X(t) of each block satisfies the compressed condition
    PAF_a'(s) + 2 PAF_b'(s) + PAF_d'(s) = 0   (s != 0 in Z_d)
and is constant on the H-orbits of Z_d. Stage 1 lists, for each block and row sum, every compression that some member
attains (dynamic program over the orbits) and joins the three compression lists exactly. Stage 2 lists, for each
compatible triple, the members of A, B and D with exactly those compressions (meet in the middle over two halves of
the orbits, keyed by row sum and compression) and joins them on all autocorrelation classes (join4x), heaviest triples
first. Nothing is sampled; with no cap reached a structure and pattern that yield nothing are ruled out.

Variant: A on the orbits of <H, -1> (symmetric), B and D on the orbits of H (the form of propus.py --freeD).

  python3 search/propus_chain.py v gens [--d D] [--rows ra,rb,rd] [--maxwork W] [--maxtriples N] [--K 1]
  python3 search/propus_chain.py v --sweep [--mino 31] [--maxo 44] ...   (every subgroup structure in the range)
  ... --plain: ordinary Hadamard matrices, Goethals-Seidel quadruples (A, B, B, D) with A on the H-orbits too
"""
import argparse, json, os, struct, subprocess, sys, tempfile, time
from collections import defaultdict
import numpy as np
from pathlib import Path as _P
sys.path[:0] = [str(_P(__file__).resolve().parent), str(_P(__file__).resolve().parent.parent / 'verify')]
import workdir as W
import propus as PR
import verify


def residue_orbits(d, H):
    oid = -np.ones(d, np.int64); reps = []
    for j in range(d):
        if oid[j] < 0:
            for h in H: oid[j * h % d] = len(reps)
            reps.append(j)
    return oid, reps


class Comp:
    """Compression of a structure st on Z_v to Z_d, on the H-orbit representatives of Z_d."""
    def __init__(self, st, v, d, H):
        self.d = d
        self.roid, self.rreps = residue_orbits(d, H)
        t = np.arange(v)
        cnt = np.zeros((st.O, d), np.int64)
        np.add.at(cnt, (st.orb, t % d), 1)
        self.cnt = cnt[:, self.rreps]                      # O x k (values at the representatives)
        self.k = len(self.rreps)

    def full(self, key):
        """Compression on all of Z_d from its values at the representatives."""
        return np.asarray(key, np.int64)[self.roid]

    def dp(self, st, r, maxstates=3e7):
        """{(compression key): count} over all sign vectors with row sum r (numpy over mixed-radix codes when they
        fit in 63 bits, a dictionary otherwise). Partial row sums that can no longer reach r are dropped; None when
        the number of states exceeds maxstates."""
        tot = np.concatenate([[int(st.sizes.sum())], self.cnt.sum(0)])          # bound of each coordinate
        rad = 2 * tot + 1
        rem = np.concatenate([np.cumsum(st.sizes[::-1])[::-1][1:], [0]]).astype(np.int64)   # sizes after orbit o
        if float(np.prod(rad.astype(np.float64))) < 2 ** 62:
            mult = np.cumprod(np.concatenate([[1], rad[:-1]])).astype(np.int64)
            off = int((tot * mult).sum())                                       # code of the all-zero key
            codes = np.array([off], np.int64); counts = np.array([1.0])
            for o in range(st.O):
                dlt = int(st.sizes[o]) * int(mult[0]) + int((self.cnt[o] * mult[1:]).sum())
                codes = np.concatenate([codes + dlt, codes - dlt]); counts = np.concatenate([counts, counts])
                keep = np.abs(r - (codes % int(rad[0]) - int(tot[0]))) <= rem[o]
                codes = codes[keep]; counts = counts[keep]
                if len(codes) > 2 * maxstates: return None
                codes, inv = np.unique(codes, return_inverse=True)
                counts = np.bincount(inv, weights=counts)
                if len(codes) > maxstates: return None
            out = {}
            for c, n in zip(codes.tolist(), counts.tolist()):
                x = c; vals = []
                for j in range(len(rad)):
                    vals.append(x % int(rad[j]) - int(tot[j])); x //= int(rad[j])
                if vals[0] == r: out[tuple(vals[1:])] = int(round(n))
            return out
        states = {(0,) + (0,) * self.k: 1}
        sizes = st.sizes
        for o in range(st.O):
            c = self.cnt[o]; s = int(sizes[o]); new = defaultdict(int)
            for key, n in states.items():
                for sg in (1, -1):
                    nk = (key[0] + sg * s,) + tuple(key[1 + i] + sg * int(c[i]) for i in range(self.k))
                    if abs(r - nk[0]) <= rem[o]: new[nk] += n
                if len(new) > maxstates / 8: return None          # tuples cost about 8 times a numpy code
            states = new
        return {k[1:]: n for k, n in states.items() if k[0] == r}


def cpaf(x):
    d = len(x)
    return np.array([int(np.dot(x, np.roll(x, -s))) for s in range(1, d // 2 + 1)], np.int64)


def cpaf_batch(keys, comp, chunk=1 << 17):
    """cpaf of every compression in keys (tuples of values at the orbit representatives), vectorized: rows of
    (len(keys), d // 2)."""
    d = comp.d; out = np.empty((len(keys), d // 2), np.int64)
    for c0 in range(0, len(keys), chunk):
        X = np.asarray(keys[c0:c0 + chunk], np.int64).reshape(-1, len(comp.rreps))[:, comp.roid]
        for s in range(1, d // 2 + 1):
            out[c0:c0 + len(X), s - 1] = np.einsum('ij,ij->i', X, np.roll(X, -s, axis=1))
    return out


def join3(PA, PB, PD, target):
    """All (i, j, k) with PA[i] + PB[j] + PD[k] == target (integer row vectors): the two smaller lists are combined
    one row at a time (numpy) and looked up by a 64-bit hash in the sorted third; candidates are checked exactly."""
    lists = [np.asarray(PA, np.int64), np.asarray(PB, np.int64), np.asarray(PD, np.int64)]
    if min(len(x) for x in lists) == 0: return []
    dim = lists[0].shape[1]
    w = np.random.default_rng(11).integers(1, 2 ** 61, size=dim, dtype=np.int64)
    h = [x @ w for x in lists]
    ht = np.asarray(target, np.int64) @ w
    order = sorted(range(3), key=lambda r: len(lists[r]))
    ix, iy, iz = order                                    # z (the largest) is looked up
    oz = np.argsort(h[iz], kind='stable'); hz = h[iz][oz]
    out = []
    for i in range(len(lists[ix])):
        need = ht - h[ix][i] - h[iy]                      # one value per row of y (wrapping int64 arithmetic)
        lo = np.searchsorted(hz, need, 'left'); hi = np.searchsorted(hz, need, 'right')
        for j in np.flatnonzero(hi > lo):
            for k in oz[lo[j]:hi[j]]:
                if np.array_equal(lists[ix][i] + lists[iy][j] + lists[iz][k], target):
                    idx = [0, 0, 0]; idx[ix] = i; idx[iy] = int(j); idx[iz] = int(k)
                    out.append(tuple(idx))
    return out


_HALVES = {}


def _halves(st, comp):
    ck = (id(st), id(comp))
    if ck in _HALVES and (_HALVES[ck][0] is not st or _HALVES[ck][1] is not comp):
        del _HALVES[ck]
    if ck not in _HALVES:
        members(st, comp, 0, (0,) * comp.k, 0)                # builds and caches the halves
    return _HALVES[ck]


_QH = {}


def members_idx(st, comp, r, key, cap):
    """(I, J): the members with row sum r and compression key are concat(Y1[I], Y2[J]) (halves as in members())."""
    _, _, Y1, K1, k1, Y2, K2, k2s, o2, w = _halves(st, comp)
    target = np.array((r,) + tuple(key), np.int64)
    need = (target @ w) - k1
    lo = np.searchsorted(k2s, need, 'left'); hi = np.searchsorted(k2s, need, 'right')
    ii = np.flatnonzero(hi > lo)
    if len(ii) == 0: return np.zeros(0, np.int64), np.zeros(0, np.int64)
    cnt = hi[ii] - lo[ii]
    if cnt.sum() > 4 * cap + 10: return None
    I = np.repeat(ii, cnt)
    starts = np.repeat(lo[ii] - np.concatenate([[0], np.cumsum(cnt)[:-1]]), cnt)
    J = o2[np.arange(len(I)) + starts]
    ok = np.all(K1[I] + K2[J] == target[None, :], axis=1)        # drop hash collisions
    I, J = I[ok], J[ok]
    if len(I) > cap: return None
    return I, J


def half_vals(st, comp, I, J, chunk=1 << 16):
    """Class values of the members concat(Y1[I], Y2[J]): Y S_c Y^T = Q1[i] + Q2[j] + Y1[i] (S12 + S21^T) Y2[j]^T, with
    the half forms Q1, Q2 cached per structure and the cross term computed per chunk (about C (O - h) operations per
    member instead of C O^2)."""
    _, _, Y1, K1, k1, Y2, K2, k2s, o2, w = _halves(st, comp)
    O = st.O; h = O // 2; C = st.C
    ck = (id(st), id(comp))
    if ck not in _QH or _QH[ck][0] is not st:
        S = st.S.astype(np.float32)                        # C x O x O
        S11, S22 = S[:, :h, :h], S[:, h:, h:]
        X = S[:, :h, h:] + np.transpose(S[:, h:, :h], (0, 2, 1))   # C x h x (O - h)
        Yf1 = Y1.astype(np.float32); Yf2 = Y2.astype(np.float32)
        Q1 = np.stack([np.einsum('ij,jk,ik->i', Yf1, S11[c], Yf1) for c in range(C)], 1)
        Q2 = np.stack([np.einsum('ij,jk,ik->i', Yf2, S22[c], Yf2) for c in range(C)], 1)
        M12 = np.ascontiguousarray(np.transpose(X, (1, 0, 2)).reshape(h, C * (O - h)))   # h x (C (O - h))
        while len(_QH) >= 2: _QH.pop(next(iter(_QH)))
        _QH[ck] = (st, Q1, Q2, M12)
    _, Q1, Q2, M12 = _QH[ck]
    out = np.empty((len(I), C), np.float32)
    for c0 in range(0, len(I), chunk):
        Ii, Jj = I[c0:c0 + chunk], J[c0:c0 + chunk]
        T = (Y1[Ii].astype(np.float32) @ M12).reshape(len(Ii), C, O - h)
        out[c0:c0 + chunk] = Q1[Ii] + Q2[Jj] + np.einsum('nck,nk->nc', T, Y2[Jj].astype(np.float32))
    return out


def members(st, comp, r, key, cap):
    """All sign vectors (int8) with row sum r and compression key (meet in the middle; the two halves are cached)."""
    O = st.O; h = O // 2
    ck = (id(st), id(comp))
    if ck in _HALVES and (_HALVES[ck][0] is not st or _HALVES[ck][1] is not comp):
        del _HALVES[ck]                                   # ids of freed objects can be reused: check identity
    if ck not in _HALVES:
        def half(lo, hi):
            m = hi - lo
            idx = np.arange(1 << m, dtype=np.int64)
            Y = (1 - 2 * ((idx[:, None] >> np.arange(m)) & 1)).astype(np.int8)
            K = np.concatenate([(Y.astype(np.int64) @ st.sizes[lo:hi])[:, None], Y.astype(np.int64) @ comp.cnt[lo:hi]], 1)
            return Y, K
        Y1, K1 = half(0, h); Y2, K2 = half(h, O)
        w = np.random.default_rng(7).integers(1, 2 ** 61, size=K1.shape[1], dtype=np.int64)
        k1 = K1 @ w; k2 = K2 @ w
        o2 = np.argsort(k2, kind='stable')
        while len(_HALVES) >= 2: _HALVES.pop(next(iter(_HALVES)))      # keep the A and the B/D structures
        _HALVES[ck] = (st, comp, Y1, K1, k1, Y2, K2, k2[o2], o2, w)   # holding st and comp keeps their ids unique
    _, _, Y1, K1, k1, Y2, K2, k2s, o2, w = _HALVES[ck]
    target = np.array((r,) + tuple(key), np.int64)
    need = (target @ w) - k1                          # K1 + K2 = target  <=>  k2 = target.w - k1 (64-bit hash)
    lo = np.searchsorted(k2s, need, 'left'); hi = np.searchsorted(k2s, need, 'right')
    out = []; total = 0
    for i in np.flatnonzero(hi > lo):
        js = o2[lo[i]:hi[i]]
        js = js[np.all(K1[i][None, :] + K2[js] == target[None, :], axis=1)]
        if len(js):
            out.append(np.concatenate([np.repeat(Y1[i][None, :], len(js), 0), Y2[js]], 1)); total += len(js)
            if total > cap: return None
    return np.concatenate(out) if out else np.zeros((0, O), np.int8)


def search(v, gens, rows, d, a, log):
    st_ad, st_b = PR.structures(v, gens, False)
    if a.plain: st_ad = st_b                     # ordinary Hadamard: A on the H-orbits too, Goethals-Seidel (A, B, B, D)
    Hs = set([1])
    while True:
        n2 = Hs | {x * g % v for x in Hs for g in gens}
        if n2 == Hs: break
        Hs = n2
    ra, rb, rd = rows; pm = v % 4
    qb = a.qb                                    # Z_qb-bordered propus: weighted PAF sum -4qb on qb Z_v (off 0), 0 off it
    tgt = np.zeros(v, np.int64)                  # (qb = 1: -4 everywhere, order 4(v + 1); qb = 0: plain propus)
    if qb: tgt[::qb] = -4 * qb
    t0 = time.time()
    for d in (d if isinstance(d, (list, tuple)) else [d]):   # candidate moduli, finest first: the first whose DP fits
        Hd = sorted({h % d for h in Hs})
        Hd_sym = sorted({h % d for h in Hs} | {(-h) % d for h in Hs})
        ca = Comp(st_ad, v, d, Hd if a.plain else Hd_sym); cb = Comp(st_b, v, d, Hd)
        Da = ca.dp(st_ad, ra); Db = cb.dp(st_b, rb) if Da is not None else None
        Dd = (ca.dp(st_ad, rd) if a.symD else cb.dp(st_b, rd)) if Db is not None else None
        if Da is not None and Db is not None and Dd is not None: break
        log(f"v={v} gens={gens} d={d} rows={rows}: compression DP above its cap")
    else:
        log(f"v={v} gens={gens} rows={rows}: no modulus with a compression DP under the cap; skipped"); return []
    log(f"v={v} gens={gens} d={d} rows={rows}: orbits A {st_ad.O}, B/D {st_b.O}, classes {st_ad.C}; compressions "
        f"A {len(Da)}, B {len(Db)}, D {len(Dd)} [{time.time() - t0:.1f}s]")
    KA, KB, KD = list(Da), list(Db), list(Dd)
    st_dd, cd = (st_ad, ca) if a.symD else (st_b, cb)       # --symD: D symmetric like A (doubly symmetric propus)
    PA = cpaf_batch(KA, ca); PB = 2 * cpaf_batch(KB, cb); PD = cpaf_batch(KD, cd)
    triples = []
    tgt_d = np.array([int(tgt[r::d].sum()) for r in range(1, d // 2 + 1)], np.int64)   # compressed target
    for i, j, k in join3(PA, PB, PD, tgt_d):
        ka, kb, kd = KA[i], KB[j], KD[k]
        triples.append((Da[ka] * Dd[kd], Da[ka], Db[kb], Dd[kd], ka, kb, kd))
    triples.sort(key=lambda x: -x[0])
    if max(st_ad.O, st_b.O) > a.maxmitm:
        log(f"  stage 1: {len(triples)} compatible triples; {st_b.O} orbits is too many for the meet in the middle, skipped")
        return []
    work = sum(t[0] for t in triples)
    log(f"  stage 1: {len(triples)} compatible triples, pair work {work:.2e} (direct join: "
        f"{st_ad.pool(ra) * min(st_b.pool(rb), st_b.pool(rd)):.2e}) [{time.time() - t0:.1f}s]")
    found = []; done_work = 0.0
    if len(triples) > 100:
        # group the triples by their D key (heaviest groups first) so that the cached lists are reused
        gw = defaultdict(float)
        for t_ in triples: gw[t_[6]] += t_[0]
        triples.sort(key=lambda t_: (-gw[t_[6]], t_[6], -t_[0]))
    cache = {}; cache_n = [0]                               # (role, key) -> (I, J, values); members held, bounded
    def lst(role, st, comp, r, key, scale):
        ck_ = (role, r, key)
        if ck_ in cache:
            cache[ck_] = cache.pop(ck_); return cache[ck_]      # most recently used last
        IJ = members_idx(st, comp, r, key, a.maxlist)
        if IJ is None: return None
        I_, J_ = IJ
        V = np.empty((len(I_), st.C), np.int16)
        for c0 in range(0, len(I_), 1 << 18):
            V[c0:c0 + (1 << 18)] = scale * np.rint((half_vals(st, comp, I_[c0:c0 + (1 << 18)], J_[c0:c0 + (1 << 18)]) - pm) / 4)
        cache[ck_] = (I_.astype(np.int32), J_.astype(np.int32), V); cache_n[0] += len(I_)
        while cache_n[0] > a.cache and len(cache) > 3:
            old_ = cache.pop(next(iter(cache))); cache_n[0] -= len(old_[0])
        return cache[ck_]
    def rows_of(st, comp, I_, J_, idx):
        _, _, Y1, _, _, Y2, _, _, _, _ = _halves(st, comp)
        return np.concatenate([Y1[I_[idx]], Y2[J_[idx]]])
    for i, (w_, na, nb, nd, ka, kb, kd) in enumerate(triples[:a.maxtriples]):
        if done_work + w_ > a.maxwork:
            log(f"  work cap reached after {i} triples ({done_work:.2e})"); break
        LA = lst('A', st_ad, ca, ra, ka, 1); LB = lst('B', st_b, cb, rb, kb, 2); LD = lst('D', st_dd, cd, rd, kd, 1)
        if LA is None or LB is None or LD is None:
            log(f"  triple {i}: a list exceeds {a.maxlist:.0e}; skipped"); continue
        Va, Vb, Vd = LA[2], LB[2], LD[2]
        Z = np.zeros((1, st_ad.C), np.int16)
        roles = sorted([(len(Va), 0), (len(Vd), 1), (len(Vb), 2)]); byrole = {0: Va, 1: Vd, 2: Vb}
        lists = [byrole[roles[0][1]], byrole[roles[1][1]], byrole[roles[2][1]], Z]
        slot = {r: j for j, (_, r) in enumerate(roles)}
        fd, path = tempfile.mkstemp(suffix='.bin', dir=W.tmpdir()); os.close(fd)
        try:
            with open(path, 'wb') as f:
                f.write(struct.pack('iiiii', *[len(x) for x in lists], st_ad.C))
                f.write(np.array([(tgt[g] - 4 * pm) // 4 for g in st_ad.reps], dtype=np.int16).tobytes())
                for x in lists: f.write(np.ascontiguousarray(x).tobytes())
            out = subprocess.run([W.bin('join4x'), path, str(a.K)], capture_output=True, text=True, check=True).stdout
        finally:
            os.remove(path)
        total = int(out.strip().split('TOTAL')[-1]); done_work += w_
        if i < 5 or i % 200 == 0 or total:
            log(f"  triple {i}: lists {len(Va)}, {len(Vb)}, {len(Vd)}: {total} solutions [{time.time() - t0:.0f}s]")
        for line in out.splitlines():
            if not line.startswith('SOL'): continue
            ids = list(map(int, line.split()[1:]))
            A = rows_of(st_ad, ca, LA[0], LA[1], ids[slot[0]]).astype(np.int64)[st_ad.orb]
            B = rows_of(st_b, cb, LB[0], LB[1], ids[slot[2]]).astype(np.int64)[st_b.orb]
            D = rows_of(st_dd, cd, LD[0], LD[1], ids[slot[1]]).astype(np.int64)[st_dd.orb]
            tot = verify.paf_sum([A, B, D]) + verify.paf_sum([B])
            if not (tot[0] == 4 * v and np.array_equal(tot[1:], tgt[1:])): log('  WRONG autocorrelations'); continue
            Kb = None
            if qb:
                Hm, Kb = PR.bordered_matrix(A, B, D) if qb == 1 else PR.qbordered_matrix(A, B, D, qb)
                n4 = 4 * (v + qb)
                ok = Hm is not None and np.array_equal(Hm, Hm.T) and np.array_equal(Hm @ Hm.T, n4 * np.eye(n4, dtype=Hm.dtype))
                log(f"  hit: Z_{qb}-bordered symmetric Hadamard of order {n4} {ok}")
            elif a.plain:
                Hm = verify.goethals_seidel(*(verify.circulant(x) for x in (A, B, B, D)))
                ok = np.array_equal(Hm @ Hm.T, 4 * v * np.eye(4 * v, dtype=Hm.dtype))
                log(f"  hit: Hadamard (Goethals-Seidel, C = B) {ok}")
            else:
                Hm = PR.symmetric_matrix_v3(A, B, D)
                ok = np.array_equal(Hm, Hm.T) and np.array_equal(Hm @ Hm.T, 4 * v * np.eye(4 * v, dtype=Hm.dtype))
                log(f"  hit: symmetric Hadamard {ok}")
            if ok:
                rec = dict(n=v, gens=list(gens), symB=False, freeD=not a.symD, plain=bool(a.plain), rows=[int(A.sum()), int(B.sum()), int(D.sum())],
                           method=f'propus_chain d={d}', seqs=[''.join('+' if x > 0 else '-' for x in s) for s in (A, B, D)])
                if qb == 1: rec.update(border=True, order=4 * (v + 1), K=Kb.tolist())
                elif qb: rec.update(border=True, q=qb, order=4 * (v + qb), R=Kb.tolist())
                with open(W.hits('gs_propus_hits.jsonl' if a.plain else 'sym_hits.jsonl'), 'a') as fh: fh.write(json.dumps(rec) + '\n')
                found.append(rec)
        if found: break
    log(f"  done: {min(len(triples), i + 1) if triples else 0} of {len(triples)} triples, work {done_work:.2e}, "
        f"{len(found)} verified [{time.time() - t0:.0f}s]")
    return found


def best_d(v, alldiv=False):
    if alldiv:                                   # every proper divisor > 1 (even v: powers of 2 and their multiples too)
        return [d for d in range(2, v) if v % d == 0]
    ps = [p for p in range(3, v) if v % p == 0 and all(p % q for q in range(2, int(p ** 0.5) + 1))]
    return ps


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('v', type=int); ap.add_argument('gens', nargs='?', default='')
    ap.add_argument('--d', type=int, default=0); ap.add_argument('--rows', default='')
    ap.add_argument('--maxwork', type=float, default=3e11); ap.add_argument('--maxtriples', type=int, default=10 ** 7)
    ap.add_argument('--maxlist', type=float, default=3e7); ap.add_argument('--K', type=int, default=1)
    ap.add_argument('--sweep', action='store_true'); ap.add_argument('--mino', type=int, default=31)
    ap.add_argument('--plain', action='store_true', help='ordinary Hadamard: Goethals-Seidel quadruples (A, B, B, D) on H-orbits')
    ap.add_argument('--maxo', type=int, default=44)
    ap.add_argument('--border', action='store_true', help='bordered propus on Z_v: PAF sum -4, symmetric Hadamard of order 4(v + 1)')
    ap.add_argument('--maxmitm', type=int, default=46, help='largest orbit count for the meet in the middle (memory ~ 2^(O/2))')
    ap.add_argument('--kdmax', type=int, default=16, help='largest number of compression coordinates tried first')
    ap.add_argument('--cache', type=float, default=2e7, help='members kept in the stage-2 list cache')
    ap.add_argument('--qborder', type=int, default=0, help='Z_q-bordered propus on Z_v (q | v): PAF sum -4q on qZ_v, 0 off it; order 4(v + q)')
    ap.add_argument('--symD', action='store_true', help='D symmetric too (doubly symmetric propus: a symmetric pair of order 2v, NOTES.md §19)')
    a = ap.parse_args()
    a.qb = a.qborder or (1 if a.border else 0)
    if a.qb > 1 and a.v % a.qb: print(f'q = {a.qb} does not divide {a.v}'); return
    log = lambda s: print(s, flush=True)
    v = a.v
    ds = [a.d] if a.d else best_d(v, bool(a.qb))
    if not ds: log("v is prime: no compression"); return
    if a.sweep:
        jobs = []
        for S, gs in PR.subgroup_gens(v).items():
            Ob = len({min(t * h % v for h in S) for t in range(v)})
            if a.mino <= Ob <= a.maxo and (a.plain or a.qb or (v - 1) not in S): jobs.append((Ob, gs))   # bordered: B, D may be symmetric too
        jobs.sort()
        log(f"v={v}: {len(jobs)} structures with {a.mino}..{a.maxo} orbits")
    else:
        jobs = [(0, [int(x) for x in a.gens.split(',')])]
    for _, gs in jobs:
        pats = [tuple(int(x) for x in a.rows.split(','))] if a.rows else PR.qborder_patterns(a.qb) if a.qb else \
            [p for p in PR.patterns(v) if not a.plain or p[0] <= p[2]]
        st_ad, st_b = PR.structures(v, gs, False)
        if a.plain: st_ad = st_b
        Hs = {1}
        while True:
            n2 = Hs | {x * g % v for x in Hs for g in gs}
            if n2 == Hs: break
            Hs = n2
        kd = {p: len(residue_orbits(p, sorted({h % p for h in Hs}))[1]) for p in ds}
        # candidate moduli: the earlier choice (the largest modulus with at most 10 orbit coordinates), unless that is
        # coarse (at most 4 coordinates) while a finer modulus has at most KDMAX coordinates and an estimated compression
        # space (product over the coordinates of the number of attainable values) of at most 1e8: then the finer ones
        # first. This sends a cyclotomic product Z_r x GF(p) of index e >= 10 to Z_p (1 + e coordinates) instead of Z_r,
        # whose pair work was infeasible (v = 183 = 3 * 61 with index 12: 4e11 instead of 1e14).
        def space(p):
            c = Comp(st_b, v, p, sorted({h % p for h in Hs})).cnt
            tot = 1.0
            for r in range(c.shape[1]):
                nz = c[:, r][c[:, r] > 0]
                tot *= (len(nz) + 1) if len(set(nz.tolist())) <= 1 else min(int(nz.sum()) + 1, 2 ** len(nz))
            return tot
        small = [p for p in ds if kd[p] <= 10]
        old = max(small) if small else min(ds, key=lambda p: kd[p])    # the earlier choice
        fine = [] if kd[old] > 4 else \
            sorted([p for p in ds if kd[old] < kd[p] <= a.kdmax and space(p) <= 1e8], key=lambda p: (-kd[p], -p))
        d = (fine + [old] + [p for p in sorted(ds, key=lambda p: kd[p]) if p not in fine and p != old]) if not a.d else ds
        for rows in pats:
            if min(st_ad.pool(rows[0]), st_b.pool(rows[1]), st_b.pool(rows[2])) == 0: continue
            if search(v, gs, rows, d, a, log):
                log(f"HIT v={v} gens={gs} rows={rows}"); return


if __name__ == '__main__':
    main()
