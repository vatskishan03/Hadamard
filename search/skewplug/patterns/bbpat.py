"""Big-big feasibility filter of defect-position patterns for the skew Paley plug-in shapes of ../sk.py (EXACT necessary condition).
Level-1 equation: sum_X PAF_{u_X}(k) = 0 (k != 0) with u = s + q e: |u| = 1 off the defects D_X, |u| in {q (type p), q-2 (type m)} on D_X,
constant on the H_K-orbits (block A: u(-x) = -u(x)).
Split PAF_X(k) = BB (both entries big) + (one big) + SS (both small):
  |sum_X BB_X(k)| <= N(k) := sum_X [ sum_{y in D_X} |u_y| ([y+k not in D_X] + [y-k not in D_X]) + #{x : x, x+k not in D_X} ].
A pattern = for each block the set of defect orbits with a type p/m.  A pattern is alive iff some orbit signs satisfy the bound at every lag class k != 0.
In the rigid regime (q-2)^2 > 2 E q + 4 |K| a lag with exactly one big-big pair is impossible (count test).
enum_patterns(sh) returns the alive patterns and the numbers tested / killed by the count test / killed by the signed test; rowfilter() then applies the
exact row-sum pattern test of ../sk_pat.py to the (allocation, m-split) class of each survivor.  A shape whose pattern list is empty (after both filters) is EXACTLY empty.
usage: python3 patterns/bbpat.py 'K' q 'hunits' [tlimit]     e.g.  python3 patterns/bbpat.py "[31]" 19 "[[5]]" 600     (lists the alive patterns, no join)"""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, '..')); sys.path.insert(0, HERE)
import numpy as np, itertools, math, time
from sk import SShape
from sk_pat import allocations


def identity_targets(sh):
    q, N = sh.q, sh.N
    out = []
    for ep in range(0, 4 * N // (q + 1) + 1):
        r = 4 * N - (q + 1) * ep
        if q > 3:
            if r % (q - 3): continue
            out.append((ep, r // (q - 3)))
        else:
            if r == 0: out.append((ep, None))
    return out


class PatCtx:
    def __init__(self, sh):
        self.sh = sh; K = sh.kCD; N = sh.N; self.N = N; self.q = sh.q
        self.pos = [np.array(p, dtype=np.int64) for p in sh.spos]       # orbits of H_K (same for A and free blocks)
        self.w = [len(p) for p in self.pos]
        self.nu = len(self.pos)
        self.zero = sh.zero_orb
        self.reps = [int(k) for k in sh.repsPM if k != 0]
        self.L = len(self.reps)
        self.SH = np.array([sh.shift[k] for k in self.reps])           # SH[l][j] = j + k_l
        # orbit incidence: C[i][j][l] = #{(a in O_i, b in O_j): b - a = k_l}
        ind = np.zeros((self.nu, N), dtype=np.int64)
        for i, p in enumerate(self.pos): ind[i, p] = 1
        self.ind = ind
        C = np.zeros((self.nu, self.nu, self.L), dtype=np.int64)
        for l in range(self.L):
            # b = a + k : sum_a ind_i[a] ind_j[a + k]
            C[:, :, l] = ind @ ind[:, self.SH[l]].T
        self.C = C
        neg = sh.negperm
        self.negorb = [int(sh.sorb[int(neg[int(p[0])])]) for p in self.pos]

    def free_cands(self, capP, capM, forced_ok=True):
        """all (tuple of (orbit, type)) with p-weight <= capP, m-weight <= capM (type 0 = p, 1 = m) for a free block"""
        res = []
        nu, w = self.nu, self.w
        def rec(i, cur, wp, wm):
            if i == nu: res.append((tuple(cur), wp, wm)); return
            rec(i + 1, cur, wp, wm)
            if wp + w[i] <= capP: cur.append((i, 0)); rec(i + 1, cur, wp + w[i], wm); cur.pop()
            if wm + w[i] <= capM: cur.append((i, 1)); rec(i + 1, cur, wp, wm + w[i]); cur.pop()
        rec(0, [], 0, 0)
        return res

    def A_cands(self, capP, capM):
        sh = self.sh; res = []; pairs = sh.apairs
        def rec(j, cur, wp, wm):
            if j == len(pairs): res.append((tuple(cur), wp, wm)); return
            rec(j + 1, cur, wp, wm)
            o, po = pairs[j]; wt = 2 * self.w[o]
            if wp + wt <= capP: cur.append((j, 0)); rec(j + 1, cur, wp + wt, wm); cur.pop()
            if wm + wt <= capM: cur.append((j, 1)); rec(j + 1, cur, wp, wm + wt); cur.pop()
        rec(0, [], 0, 0)
        return res

    def items_of(self, X, cand):
        """list of (orbit, type) 'items' (A: both orbits of a pair, the second with flipped sign)"""
        if X == 0:
            out = []
            for (j, ty) in cand:
                o, po = self.sh.apairs[j]; out.append((o, ty, +1)); out.append((po, ty, -1))
            return out
        return [(o, ty, +1) for (o, ty) in cand]

    def dset_vec(self, items):
        d = np.zeros(self.N, dtype=np.int64)
        for (o, ty, sg) in items: d[self.pos[o]] = 1
        return d

    def cvec(self, items):
        """lag count vector c_X(l) = |D_X cap (D_X - k_l)|"""
        d = self.dset_vec(items)
        return (d[None, :] * d[self.SH]).sum(1)


def check_pattern(ctx, blocks_items, want_signs=True, rigid=False):
    """blocks_items: list of 4 item lists.  Returns (count_ok, signed_ok, witness)"""
    q = ctx.q; N = ctx.N; L = ctx.L
    c = np.zeros(L, dtype=np.int64)
    ds = []
    for items in blocks_items:
        d = ctx.dset_vec(items); ds.append(d)
        c += (d[None, :] * d[ctx.SH]).sum(1)
    if rigid and (c == 1).any(): return False, False, None
    # bound N(k)
    Nk = np.zeros(L, dtype=np.int64)
    for X, items in enumerate(blocks_items):
        d = ds[X]; n = int(d.sum())
        mag = np.zeros(N, dtype=np.int64)
        for (o, ty, sg) in items: mag[ctx.pos[o]] = (q - 2) if ty else q
        dp = d[ctx.SH]                                  # d[j + k]
        # shift by -k: use inverse index
        dm = np.zeros_like(dp)
        for l in range(L):
            dm[l, ctx.SH[l]] = d                         # dm[l][j + k] = d[j]  -> dm[l][y] = d[y - k]
        bs = (mag[None, :] * ((1 - dp) + (1 - dm))).sum(1) * 1       # sum_y |u_y| ([y+k not in D] + [y-k not in D])  (y in D via mag)
        cc = (d[None, :] * dp).sum(1)
        ss = N - 2 * n + cc
        Nk += bs + ss
    if not want_signs: return True, True, None
    # signed test
    allitems = []; blk = []
    for X, items in enumerate(blocks_items):
        for it in items: allitems.append(it); blk.append(X)
    ni = len(allitems)
    if ni == 0: return True, True, ()
    mags = np.array([(q - 2) if ty else q for (o, ty, sg) in allitems], dtype=np.int64)
    # sign variables: one per block item, A pair items tied (second = -first); gauge: first free item of each block fixed +
    vars_ = []     # each variable controls a list of item indices with multipliers
    seen_pair = {}
    for idx, (o, ty, sg) in enumerate(allitems):
        X = blk[idx]
        if X == 0:
            if sg == +1: vars_.append([(idx, +1)])
            else: vars_[-1].append((idx, -1))
        else: vars_.append([(idx, +1)])
    nv = len(vars_)
    # gauge fixing: first variable of each block
    fixed = set(); firstvar = {}
    for vi, grp in enumerate(vars_):
        X = blk[grp[0][0]]
        if X not in firstvar: firstvar[X] = vi
    free = [vi for vi in range(nv) if vi not in firstvar.values()]
    ncomb = 1 << len(free)
    # build sign matrix (ncomb, ni)
    sgn = np.ones((ncomb, nv), dtype=np.int64)
    for b, vi in enumerate(free):
        sgn[:, vi] = 1 - 2 * ((np.arange(ncomb) >> b) & 1)
    S = np.zeros((ncomb, ni), dtype=np.int64)
    for vi, grp in enumerate(vars_):
        for (idx, mult) in grp: S[:, idx] = sgn[:, vi] * mult * mags[idx]
    orb = [it[0] for it in allitems]
    BB = np.zeros((ncomb, L), dtype=np.int64)
    for a in range(ni):
        for b in range(ni):
            if blk[a] != blk[b]: continue
            Cab = ctx.C[orb[a], orb[b]]
            if not Cab.any(): continue
            BB += (S[:, a] * S[:, b])[:, None] * Cab[None, :]
    ok = (np.abs(BB) <= Nk[None, :]).all(1)
    if ok.any(): return True, True, S[np.argmax(ok)].tolist()
    return True, False, None


def enum_patterns(sh, targets=None, maxleaves=10 ** 7, verbose=False, signs=True, tlimit=1e18):
    ctx = PatCtx(sh); q = sh.q; N = sh.N
    tg = identity_targets(sh) if targets is None else targets
    Etot = max([ep + (em or 0) for (ep, em) in tg] or [0])
    # the count test (a lag with exactly one big-big term) is valid iff the smallest term (q-2)^2 exceeds every possible bound N(k) <= 2 E q + 4 N
    rigid = (q - 2) ** 2 > 2 * Etot * q + 4 * N
    t0 = time.time()
    alive = []; ntest = 0; kill_count = 0; kill_sign = 0; timed_out = False
    for (ep, em) in tg:
        em_ = 0 if em is None else em
        # candidates per block with weights bounded by the targets
        capP, capM = ep, em_
        if sh.q == 3: capM = 0       # q = 3: e_m unconstrained by the identity (coefficient q - 3 = 0): type m allowed with any weight
        candsA = ctx.A_cands(capP, capM if sh.q > 3 else 10 ** 9)
        candsF = ctx.free_cands(capP, capM if sh.q > 3 else 10 ** 9)
        if verbose: print('target', (ep, em), 'A cands', len(candsA), 'free cands', len(candsF), flush=True)
        # index candidates by (wp, wm)
        byw = {}
        for i, (c, wp, wm) in enumerate(candsF): byw.setdefault((wp, wm), []).append(i)
        wkeys = sorted(byw)
        itemsF = [ctx.items_of(1, c) for (c, wp, wm) in candsF]
        cvF = [ctx.cvec(it) for it in itemsF]
        for (cA, wpA, wmA) in candsA:
            itA = ctx.items_of(0, cA); cvA = ctx.cvec(itA)
            # B <= C <= D by candidate index; total weights must hit (ep, em)
            for (wb, wc, wd) in itertools.product(wkeys, repeat=3):
                if wpA + wb[0] + wc[0] + wd[0] != ep: continue
                if em is not None and wmA + wb[1] + wc[1] + wd[1] != em: continue
                for ib in byw[wb]:
                    for ic in byw[wc]:
                        if ic < ib: continue
                        cAB = cvA + cvF[ib] + cvF[ic]
                        for id_ in byw[wd]:
                            if id_ < ic: continue
                            ntest += 1
                            if ntest > maxleaves or time.time() - t0 > tlimit: timed_out = True; break
                            c = cAB + cvF[id_]
                            if rigid and (c == 1).any(): kill_count += 1; continue
                            ok1, ok2, wit = check_pattern(ctx, [itA, itemsF[ib], itemsF[ic], itemsF[id_]], want_signs=signs, rigid=rigid)
                            if not ok1: kill_count += 1
                            elif not ok2: kill_sign += 1
                            else: alive.append(dict(A=cA, B=candsF[ib][0], C=candsF[ic][0], D=candsF[id_][0], tg=(ep, em), sign=wit))
                        if timed_out: break
                    if timed_out: break
                if timed_out: break
            if timed_out: break
        if timed_out: break
    return dict(alive=alive, ntest=ntest, kill_count=kill_count, kill_sign=kill_sign, timed_out=timed_out, secs=time.time() - t0, targets=tg)


def weights_ms(sh, A, B, C, D):
    """(al, ms) of a pattern: A pairs (j, type), free blocks ((orbit, type)...) -> position weights and m-type weights per block (A, B, C, D)"""
    pos = sh.spos; w = lambda o: len(pos[o])
    nA = sum(2 * w(sh.apairs[j][0]) for (j, ty) in A); mA = sum(2 * w(sh.apairs[j][0]) for (j, ty) in A if ty == 1)
    al = [nA]; ms = [mA]
    for blk in (B, C, D):
        al.append(sum(w(o) for (o, ty) in blk)); ms.append(sum(w(o) for (o, ty) in blk if ty == 1))
    return tuple(al), tuple(ms)


def rowfilter(sh, alive, cache=None):
    """exact row-sum pattern filter (sk_pat.row_patterns_skew) of the (allocation, m-split) class of every alive pattern; returns (alive', number removed)"""
    from sk_pat import row_patterns_skew
    cache = {} if cache is None else cache
    keep = []
    for a in alive:
        al, ms = weights_ms(sh, a['A'], a['B'], a['C'], a['D'])
        key = (al, ms)
        if key not in cache: cache[key] = bool(row_patterns_skew(sh, al, ms))
        if cache[key]: keep.append(a)
    return keep, len(alive) - len(keep)


if __name__ == '__main__':
    import json
    K = json.loads(sys.argv[1]); q = int(sys.argv[2]); hu = json.loads(sys.argv[3])
    sh = SShape(K, q, hu)
    r = enum_patterns(sh, verbose=True, tlimit=float(sys.argv[4]) if len(sys.argv) > 4 else 600)
    print('nu', sh.kCD.nu, 'targets', r['targets'], 'tested', r['ntest'], 'killed(count)', r['kill_count'], 'killed(signed)', r['kill_sign'], 'alive', len(r['alive']), 'timed_out', r['timed_out'], 'secs', round(r['secs'], 1))
    for a in r['alive'][:10]: print(a)
