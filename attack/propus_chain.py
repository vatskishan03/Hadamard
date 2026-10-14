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

  python3 attack/propus_chain.py v gens [--d D] [--rows ra,rb,rd] [--maxwork W] [--maxtriples N] [--K 1]
  python3 attack/propus_chain.py v --sweep [--mino 31] [--maxo 44] ...   (every subgroup structure in the range)
  ... --plain: ordinary Hadamard matrices, Goethals-Seidel quadruples (A, B, B, D) with A on the H-orbits too
"""
import argparse, json, math, os, struct, subprocess, sys, tempfile, time
from collections import defaultdict
import numpy as np
from pathlib import Path as _P
sys.path[:0] = [str(_P(__file__).resolve().parent), str(_P(__file__).resolve().parent.parent / 'verify')]
import workdir as W
import spence_lib as L
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


def members(st, comp, r, key, cap):
    """All sign vectors (int8) with row sum r and compression key (meet in the middle; the two halves are cached)."""
    O = st.O; h = O // 2
    ck = (id(st), id(comp))
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
        _HALVES[ck] = (Y1, K1, k1, Y2, K2, k2[o2], o2, w)
    Y1, K1, k1, Y2, K2, k2s, o2, w = _HALVES[ck]
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
    Hd = sorted({h % d for h in Hs})
    Hd_sym = sorted({h % d for h in Hs} | {(-h) % d for h in Hs})
    ca = Comp(st_ad, v, d, Hd if a.plain else Hd_sym); cb = Comp(st_b, v, d, Hd)
    ra, rb, rd = rows; pm = v % 4
    T = -4 if a.border else 0                    # bordered propus: weighted PAF sum -4 off 0 (order 4(v + 1))
    t0 = time.time()
    Da = ca.dp(st_ad, ra); Db = cb.dp(st_b, rb); Dd = cb.dp(st_b, rd)
    if Da is None or Db is None or Dd is None:
        log(f"v={v} gens={gens} d={d} rows={rows}: compression DP above its cap; skipped"); return []
    log(f"v={v} gens={gens} d={d} rows={rows}: orbits A {st_ad.O}, B/D {st_b.O}, classes {st_ad.C}; compressions "
        f"A {len(Da)}, B {len(Db)}, D {len(Dd)} [{time.time() - t0:.1f}s]")
    KA, KB, KD = list(Da), list(Db), list(Dd)
    PA = np.array([cpaf(ca.full(k)) for k in KA]).reshape(len(KA), -1)
    PB = np.array([2 * cpaf(cb.full(k)) for k in KB]).reshape(len(KB), -1)
    PD = np.array([cpaf(cb.full(k)) for k in KD]).reshape(len(KD), -1)
    triples = []
    for i, j, k in join3(PA, PB, PD, np.full(d // 2, T * (v // d), np.int64)):
        ka, kb, kd = KA[i], KB[j], KD[k]
        triples.append((Da[ka] * Dd[kd], Da[ka], Db[kb], Dd[kd], ka, kb, kd))
    triples.sort(key=lambda x: -x[0])
    if max(st_ad.O, st_b.O) > 46:
        log(f"  stage 1: {len(triples)} compatible triples; {st_b.O} orbits is too many for the meet in the middle, skipped")
        return []
    work = sum(t[0] for t in triples)
    log(f"  stage 1: {len(triples)} compatible triples, pair work {work:.2e} (direct join: "
        f"{st_ad.pool(ra) * min(st_b.pool(rb), st_b.pool(rd)):.2e}) [{time.time() - t0:.1f}s]")
    found = []; done_work = 0.0
    for i, (w_, na, nb, nd, ka, kb, kd) in enumerate(triples[:a.maxtriples]):
        if done_work + w_ > a.maxwork:
            log(f"  work cap reached after {i} triples ({done_work:.2e})"); break
        Ya = members(st_ad, ca, ra, ka, a.maxlist); Yb = members(st_b, cb, rb, kb, a.maxlist); Yd = members(st_b, cb, rd, kd, a.maxlist)
        if Ya is None or Yb is None or Yd is None:
            log(f"  triple {i}: a list exceeds {a.maxlist:.0e}; skipped"); continue
        def vals(st, Y, scale):
            return (scale * np.rint((st.class_values(Y) - pm) / 4)).astype(np.int16)
        Va, Vb, Vd = vals(st_ad, Ya, 1), vals(st_b, Yb, 2), vals(st_b, Yd, 1)
        Z = np.zeros((1, st_ad.C), np.int16)
        roles = sorted([(len(Va), 0), (len(Vd), 1), (len(Vb), 2)]); byrole = {0: Va, 1: Vd, 2: Vb}
        lists = [byrole[roles[0][1]], byrole[roles[1][1]], byrole[roles[2][1]], Z]
        slot = {r: j for j, (_, r) in enumerate(roles)}
        fd, path = tempfile.mkstemp(suffix='.bin', dir=W.tmpdir()); os.close(fd)
        try:
            with open(path, 'wb') as f:
                f.write(struct.pack('iiiii', *[len(x) for x in lists], st_ad.C))
                f.write(np.full(st_ad.C, (T - 4 * pm) // 4, dtype=np.int16).tobytes())
                for x in lists: f.write(np.ascontiguousarray(x).tobytes())
            out = subprocess.run([W.bin('join4x'), path, str(a.K)], capture_output=True, text=True, check=True).stdout
        finally:
            os.remove(path)
        total = int(out.strip().split('TOTAL')[-1]); done_work += w_
        if i < 5 or i % 200 == 0 or total:
            log(f"  triple {i}: lists {len(Ya)}, {len(Yb)}, {len(Yd)}: {total} solutions [{time.time() - t0:.0f}s]")
        for line in out.splitlines():
            if not line.startswith('SOL'): continue
            ids = list(map(int, line.split()[1:]))
            A = Ya[ids[slot[0]]].astype(np.int64)[st_ad.orb]; B = Yb[ids[slot[2]]].astype(np.int64)[st_b.orb]
            D = Yd[ids[slot[1]]].astype(np.int64)[st_b.orb]
            tot = verify.paf_sum([A, B, D]) + verify.paf_sum([B])
            if not (tot[0] == 4 * v and np.all(tot[1:] == T)): log('  WRONG autocorrelations'); continue
            Kb = None
            if a.border:
                Hm, Kb = PR.bordered_matrix(A, B, D)
                n4 = 4 * (v + 1)
                ok = Hm is not None and np.array_equal(Hm, Hm.T) and np.array_equal(Hm @ Hm.T, n4 * np.eye(n4, dtype=Hm.dtype))
                log(f"  hit: bordered symmetric Hadamard of order {n4} {ok}")
            elif a.plain:
                Hm = verify.goethals_seidel(*(verify.circulant(x) for x in (A, B, B, D)))
                ok = np.array_equal(Hm @ Hm.T, 4 * v * np.eye(4 * v, dtype=Hm.dtype))
                log(f"  hit: Hadamard (Goethals-Seidel, C = B) {ok}")
            else:
                Hm = PR.symmetric_matrix_v3(A, B, D)
                ok = np.array_equal(Hm, Hm.T) and np.array_equal(Hm @ Hm.T, 4 * v * np.eye(4 * v, dtype=Hm.dtype))
                log(f"  hit: symmetric Hadamard {ok}")
            if ok:
                rec = dict(n=v, gens=list(gens), symB=False, freeD=True, plain=bool(a.plain), rows=[int(A.sum()), int(B.sum()), int(D.sum())],
                           method=f'propus_chain d={d}', seqs=[''.join('+' if x > 0 else '-' for x in s) for s in (A, B, D)])
                if a.border: rec.update(border=True, order=4 * (v + 1), K=Kb.tolist())
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
    a = ap.parse_args()
    log = lambda s: print(s, flush=True)
    v = a.v
    ds = [a.d] if a.d else best_d(v, a.border)
    if not ds: log("v is prime: no compression"); return
    if a.sweep:
        jobs = []
        for S, gs in PR.subgroup_gens(v).items():
            Ob = len({min(t * h % v for h in S) for t in range(v)})
            if a.mino <= Ob <= a.maxo and (a.plain or (v - 1) not in S): jobs.append((Ob, gs))
        jobs.sort()
        log(f"v={v}: {len(jobs)} structures with {a.mino}..{a.maxo} orbits")
    else:
        jobs = [(0, [int(x) for x in a.gens.split(',')])]
    for _, gs in jobs:
        pats = [tuple(int(x) for x in a.rows.split(','))] if a.rows else [(2, 0, 0), (0, 0, 2)] if a.border else \
            [p for p in PR.patterns(v) if not a.plain or p[0] <= p[2]]
        st_ad, st_b = PR.structures(v, gs, False)
        if a.plain: st_ad = st_b
        Hs = {1}
        while True:
            n2 = Hs | {x * g % v for x in Hs for g in gs}
            if n2 == Hs: break
            Hs = n2
        kd = {p: len(residue_orbits(p, sorted({h % p for h in Hs}))[1]) for p in ds}
        small = [p for p in ds if kd[p] <= 10]
        d = max(small) if small else min(ds, key=lambda p: kd[p])    # the largest modulus whose compressions stay small
        for rows in pats:
            if min(st_ad.pool(rows[0]), st_b.pool(rows[1]), st_b.pool(rows[2])) == 0: continue
            if search(v, gs, rows, d, a, log):
                log(f"HIT v={v} gens={gs} rows={rows}"); return


if __name__ == '__main__':
    main()
