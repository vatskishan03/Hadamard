"""Exact compression-chain search for skew Goethals-Seidel families on Z_v with a multiplier subgroup H.

This is the method METHOD.md describes for order 3964, specialised to skew families on cyclic groups.

**Compression.** For d | v, the compression of a block X to Z_d is x_j = sum_{t = j (mod d)} X(t). If the four
blocks are complementary, so are their compressions: sum_X PAF_x(k) = 4v [k = 0] on Z_d, in integers. Every
block is H-invariant, so its compression is constant on the H-orbits of Z_d. It is stored by its values on
orbit representatives (reduced coordinates).

**Search.**
1. **First level.** A dynamic program over the H-orbits of Z_v lists, for each block, every pair of
   compressions (to Z_d1 and to Z_d2) that some sign vector attains, with multiplicities. The Z_d1
   compressions of the four blocks are joined exactly.
2. **Second level.** For each solution of the first level, the Z_d2 compressions inside the four fibres are
   joined exactly.
3. **Final join.** For each consistent combination, every sign vector with that pair of compressions is listed
   (meet in the middle over the orbits). The four lists are joined on all autocorrelation classes (join4x).
   Hits are verified as skew-Hadamard matrices.

Nothing is sampled. With --exhaustive and no truncation, a run that finds nothing rules the structure and
pattern out; any truncation is reported.

  python3 chain.py v gens b,c,d [--d1 p --d2 q] [--exhaustive] [--maxfibre N] [--maxstates N]
"""
import argparse, json, math, os, struct, subprocess, sys, tempfile, time
from collections import defaultdict
import numpy as np
from pathlib import Path as _P; sys.path[:0] = [str(_P(__file__).resolve().parent), str(_P(__file__).resolve().parent.parent / 'verify')]
import workdir as W
import spence_lib as L
import verify


def h_orbits_mod(d, H):
    """H-orbits on Z_d (H acting by multiplication mod d): (orbit id per residue, representatives, sizes)."""
    oid = -np.ones(d, dtype=np.int64); reps = []
    for j in range(d):
        if oid[j] < 0:
            orb = {j * h % d for h in H}
            for x in orb: oid[x] = len(reps)
            reps.append(j)
    sizes = np.bincount(oid, minlength=len(reps))
    return oid, reps, sizes


class Chain:
    def __init__(self, v, gens, ds, skew=True):
        self.v = v; self.skew = skew
        G = L.AbGroup([v]); self.G = G
        self.st = st = L.Structure(G, [G.perm_from_matrix(np.array([[g]])) for g in gens])
        H = {1}
        while True:
            H2 = H | {h * g % v for h in H for g in gens}
            if H2 == H: break
            H = H2
        self.H = sorted(H)
        rep = np.array([np.flatnonzero(st.orb == o)[0] for o in range(st.O)])
        partner = st.orb[G.neg()[rep]]
        if skew:
            assert not np.any(partner[1:] == np.arange(1, st.O)), "-1 fixes an orbit"
        self.pairs = [(o, int(partner[o])) for o in range(1, st.O) if o < partner[o]] if skew else []
        self.ds = ds
        self.lev = []                       # per level: (d, oid, reps, sizes, contrib_red[O, h_d])
        for d in ds:
            oid, reps, sizes = h_orbits_mod(d, self.H)
            C = np.zeros((st.O, len(reps)), dtype=np.int64)
            t = np.arange(v)
            for r_i, r in enumerate(reps):
                mask = (t % d) == r
                C[:, r_i] = np.bincount(st.orb[mask], minlength=st.O)
            self.lev.append((d, oid, reps, sizes, C))
        self.Ctot = np.concatenate([l[4] for l in self.lev], axis=1)      # O x (sum h_d)
        self.splits = np.cumsum([0] + [len(l[2]) for l in self.lev])

    # ---- dynamic program over orbits: distinct joint compressions with multiplicities ----------------
    def states(self, role_rowsum, maxstates):
        """role_rowsum None = skew block A. Returns (keys [n x dim] int16, counts [n] float64).
        Rows are deduplicated through a byte view (much faster than np.unique(axis=0))."""
        Ct = self.Ctot.astype(np.int16)
        dim = Ct.shape[1]
        vt = np.dtype((np.void, 2 * dim))
        def dedup(K, cnt):
            key = np.ascontiguousarray(K).view(vt).ravel()
            u, idx, inv = np.unique(key, return_index=True, return_inverse=True)
            return K[idx], np.bincount(inv.ravel(), weights=cnt)
        if role_rowsum is None:
            K = Ct[0][None, :].copy(); cnt = np.ones(1)
            for o, q in self.pairs:
                it = Ct[o] - Ct[q]
                K, cnt = dedup(np.concatenate([K + it, K - it]), np.concatenate([cnt, cnt]))
                if len(K) > maxstates: return None
            return K, cnt
        K = (-Ct.sum(0))[None, :].astype(np.int16); cnt = np.ones(1)
        for o in range(self.st.O):
            K, cnt = dedup(np.concatenate([K, K + 2 * Ct[o]]), np.concatenate([cnt, cnt]))
            if len(K) > maxstates: return None
        w = self.lev[0][3]
        rs = K[:, :len(w)].astype(np.int64) @ w
        keep = rs == role_rowsum
        return K[keep], cnt[keep]

    def paf(self, level, red):
        """Compressed autocorrelations (shifts 1..d//2) of reduced vectors red [n x h_d]."""
        d, oid = self.lev[level][0], self.lev[level][1]
        X = red[:, oid]
        return np.stack([(X * np.roll(X, -k, axis=1)).sum(1) for k in range(1, d // 2 + 1)], axis=1)


class JoinTooLarge(Exception):
    pass


def join4(P, target=None, maxpairs=6e7):
    """All (i, j, k, l) with P0[i] + P1[j] + P2[k] + P3[l] == target (integer vectors). Pair-sum hash join.
    Raises JoinTooLarge when the stored pair sums (two smallest lists) would exceed maxpairs."""
    sizes = [len(p) for p in P]
    if min(sizes) == 0: return []
    perm = list(np.argsort(sizes))            # two smallest lists form the stored pair sums
    if float(sizes[perm[0]]) * sizes[perm[1]] > maxpairs: raise JoinTooLarge()
    if perm != [0, 1, 2, 3]:
        sub = join4([P[k] for k in perm], target, maxpairs)
        inv = np.argsort(perm)
        return [tuple(t[inv[k]] for k in range(4)) for t in sub]
    dim = P[0].shape[1]
    if target is None: target = np.zeros(dim, dtype=np.int64)
    rng = np.random.default_rng(12345)
    R = rng.integers(1, 2 ** 62, size=dim, dtype=np.int64).astype(np.uint64)
    hs = [(p.astype(np.int64).astype(np.uint64) * R).sum(1, dtype=np.uint64) for p in P]
    ht = (target.astype(np.int64).astype(np.uint64) * R).sum(dtype=np.uint64)
    left = (hs[0][:, None] + hs[1][None, :]).ravel()
    order = np.argsort(left, kind='stable'); ls = left[order]
    out = []
    for k0 in range(0, len(hs[2]), 2048):
        right = (ht - (hs[2][k0:k0 + 2048, None] + hs[3][None, :])).ravel()
        lo = np.searchsorted(ls, right, 'left'); hi = np.searchsorted(ls, right, 'right')
        for r in np.flatnonzero(hi > lo):
            for pos in range(lo[r], hi[r]):
                ij = order[pos]; i, j = divmod(int(ij), len(hs[1]))
                k, l = divmod(int(r), len(hs[3])); k += k0
                if np.array_equal(P[0][i] + P[1][j] + P[2][k] + P[3][l], target): out.append((i, j, k, l))
    return out


class Halves:
    """Meet-in-the-middle tables of one block: all half assignments with their joint compressions, sorted by hash."""
    def __init__(self, ch, role_rowsum):
        st, Ct = ch.st, ch.Ctot
        self.skew = role_rowsum is None
        if self.skew:
            self.items = list(ch.pairs); self.base = Ct[0].copy(); vecs = [Ct[o] - Ct[q] for o, q in self.items]
        else:
            self.items = list(range(st.O)); self.base = -Ct.sum(0); vecs = [2 * Ct[o] for o in self.items]
        m = len(self.items); self.h = h = m // 2; self.m = m; self.O = st.O
        self.R = np.random.default_rng(7).integers(1, 2 ** 62, size=Ct.shape[1], dtype=np.int64).astype(np.uint64)
        def half(lo, hi):
            n = hi - lo
            idx = np.arange(1 << n, dtype=np.int64)
            bits = ((idx[:, None] >> np.arange(n)) & 1).astype(np.int8)
            V = np.array(vecs[lo:hi], dtype=np.int64).reshape(n, -1)
            S = ((2 * bits.astype(np.int64) - 1) if self.skew else bits.astype(np.int64)) @ V
            return bits, S.astype(np.int32)
        self.b1, self.s1 = half(0, h); self.b2, self.s2 = half(h, m)
        self.k1 = (self.s1.astype(np.int64).astype(np.uint64) * self.R).sum(1, dtype=np.uint64)
        k2 = (self.s2.astype(np.int64).astype(np.uint64) * self.R).sum(1, dtype=np.uint64)
        self.order = np.argsort(k2); self.ks = k2[self.order]

    def fibre(self, key, cap):
        need = np.asarray(key, dtype=np.int64) - self.base
        kn = (need.astype(np.uint64) * self.R).sum(dtype=np.uint64)
        want = kn - self.k1
        lo = np.searchsorted(self.ks, want, 'left'); hi = np.searchsorted(self.ks, want, 'right')
        cnt = hi - lo
        tot = int(cnt.sum())
        if tot > cap: return None
        Y = np.zeros((0, self.O), dtype=np.int8)
        parts = []
        nz = np.flatnonzero(cnt)
        for c0 in range(0, len(nz), 1 << 16):                 # vectorised over all hash matches, in chunks
            ii = nz[c0:c0 + (1 << 16)]
            k = cnt[ii]
            I = np.repeat(ii, k)
            off = np.arange(int(k.sum())) - np.repeat(np.cumsum(k) - k, k)
            J = self.order[lo[I] + off]
            ok = np.all(self.s1[I].astype(np.int64) + self.s2[J] == need, axis=1)
            I, J = I[ok], J[ok]
            if len(I) == 0: continue
            B = np.concatenate([self.b1[I], self.b2[J]], axis=1).astype(np.int8)
            Yc = np.zeros((len(I), self.O), dtype=np.int8)
            if self.skew:
                Yc[:, 0] = 1
                for c, (o, q) in enumerate(self.items): Yc[:, o] = 2 * B[:, c] - 1; Yc[:, q] = -(2 * B[:, c] - 1)
            else:
                for c, o in enumerate(self.items): Yc[:, o] = 2 * B[:, c] - 1
            parts.append(Yc)
        return np.concatenate(parts) if parts else Y


JOIN_TIMEOUT = None          # seconds per final join (None: no limit); set by run(join_timeout=...)


class JoinTimeout(Exception):
    pass


def final_join(ch, lists, K=2):
    st = ch.st; pm = ch.v % 4
    Vs = []
    for Y in lists:
        Yf = Y.astype(np.float64)
        Vs.append(np.rint((st.class_values(Yf) - pm) / 4).astype(np.int16))
    # join4x builds the smaller of the pair products (0,1), (2,3): pair the smallest list with the largest
    order = np.argsort([len(x) for x in Vs]); perm = [order[0], order[3], order[1], order[2]]
    fd, path = tempfile.mkstemp(suffix='.bin', dir=W.tmpdir()); os.close(fd)
    try:
        with open(path, 'wb') as f:
            f.write(struct.pack('iiiii', *[len(Vs[k]) for k in perm], st.C)); f.write(np.full(st.C, -pm, dtype=np.int16).tobytes())
            for k in perm: f.write(np.ascontiguousarray(Vs[k]).tobytes())
        try:
            out = subprocess.run([W.bin('join4x'), path, str(K)], capture_output=True, text=True, check=True,
                                 timeout=JOIN_TIMEOUT).stdout
        except subprocess.TimeoutExpired:
            raise JoinTimeout()
    finally:
        os.remove(path)
    total = int(out.strip().split('TOTAL')[-1])
    sols = []
    for line in out.splitlines():
        if line.startswith('SOL'):
            ids = list(map(int, line.split()[1:]))
            q = [None] * 4
            for slot, k in enumerate(perm): q[k] = lists[k][ids[slot]]
            sols.append(q)
    return total, sols



def split_join(ch, lists, extra):
    """Final join through extra compression levels computed from the enumerated fibres (used when the joint DP
    over all levels does not fit). Each list is grouped by its compressions to the extra primes, the compressed
    autocorrelations are joined exactly, and join4x runs only on the consistent sub-lists. Exact."""
    groups = []
    for Y in lists:
        keys = np.concatenate([Y.astype(np.int64) @ e.lev[0][4] for e in extra], axis=1)
        u, inv = np.unique(keys, axis=0, return_inverse=True)
        groups.append((u, inv.ravel()))
    P = []
    for (u, _) in groups:
        parts, a = [], 0
        for e in extra:
            h = e.lev[0][4].shape[1]
            parts.append(e.paf(0, u[:, a:a + h])); a += h
        P.append(np.concatenate(parts, axis=1))
    total, sols = 0, []
    try:
        quads = join4(P)
    except JoinTooLarge:
        return final_join(ch, lists)
    for q in quads:
        sub = [lists[k][groups[k][1] == q[k]] for k in range(4)]
        t, s = final_join(ch, sub)
        total += t; sols += s
    return total, sols

def verify_sol(ch, Ys):
    seqs = [y.astype(np.int64)[ch.st.orb] for y in Ys]
    tot = verify.paf_sum(seqs)
    if not (tot[0] == 4 * ch.v and not tot[1:].any()): return None
    Hm = verify.goethals_seidel(*map(verify.circulant, seqs)).astype(np.float64)
    ok = np.array_equal(Hm @ Hm.T, 4 * ch.v * np.eye(4 * ch.v))
    if ch.skew: ok = ok and np.array_equal(Hm + Hm.T, 2 * np.eye(4 * ch.v))
    return seqs if ok else None


def run(v, gens, pat, d1=None, d2=None, exhaustive=False, maxfibre=3_000_000, maxstates=20_000_000, verbose=True, budget=3600,
        maxbranches=None, join_timeout=None, levels=None, split=False):
    global JOIN_TIMEOUT
    JOIN_TIMEOUT = join_timeout
    """pat = (b, c, d) for skew families (first block skew, row sum 1); pat = (a, b, c, d) for plain Hadamard GS."""
    t0 = time.time()
    skew = len(pat) == 3
    roles = ([None] if skew else []) + list(pat)
    # Screen each prime level with its own DP (levels whose state space is too large are skipped). An
    # empty join at any level rules the pattern out. The chain then uses the two feasible levels with the
    # fewest level-1 solutions, jointly if their joint DP fits, else the single best level.
    ps = [p for p in range(2, v) if v % p == 0 and all(p % r for r in range(2, int(p ** 0.5) + 1))]
    counts1 = {}
    for p in ps:
        c1 = Chain(v, gens, [p], skew=skew)
        Sd = [c1.states(r, maxstates) for r in roles]      # intermediate DP states can exceed the final count 100-fold
        if any(x is None for x in Sd): continue
        U = [np.unique(K_, axis=0).astype(np.int64) for K_, _ in Sd]
        try:
            n1 = 0 if any(len(u) == 0 for u in U) else len(join4([c1.paf(0, u) for u in U]))
        except JoinTooLarge:
            continue
        counts1[p] = n1
        if n1 == 0:
            print(f"  RESULT v={v} H={gens} rows=(1,{pat}): ruled out at d={p} (empty level-1 join) [{time.time() - t0:.0f}s]", flush=True)
            return dict(status='ruled out', found=[], branches=0, truncated=0, solutions=0, reason=f'd={p}')
    if not counts1:
        print(f"  state space too large at every level", flush=True); return dict(status='too-large')
    best = sorted(counts1, key=lambda p: counts1[p])
    S = None
    for levels_ in ([list(levels)] if levels else [[best[0], best[1]] if len(best) > 1 else [], [best[0]]]):
        if not levels_: continue
        ch = Chain(v, gens, levels_, skew=skew)
        S = [ch.states(r, maxstates) for r in roles]
        if all(x is not None for x in S): break
        S = None
    if S is None:
        print(f"  joint state space too large", flush=True); return dict(status='too-large')
    def project(level):
        a, b = ch.splits[level], ch.splits[level + 1]
        out = []
        for K_, c in S:
            u, inv = np.unique(K_[:, a:b], axis=0, return_inverse=True)
            out.append((u.astype(np.int64), inv.ravel(), np.bincount(inv.ravel(), weights=c)))
        return out
    L1 = 0; L2 = 1 if len(ch.ds) > 1 else None
    extra = [Chain(v, gens, [p], skew=skew) for p in ps if p not in ch.ds and p < v] if split else []   # levels the joint DP left out
    d1, d2 = ch.ds[0], (ch.ds[1] if L2 is not None else 0)
    proj = project(L1)
    P1 = [ch.paf(L1, u) for u, _, _ in proj]
    try:
        sols1 = join4(P1)
    except JoinTooLarge:
        print(f"  level-1 join too large", flush=True); return dict(status='too-large')
    if verbose:
        print(f"v={v} H={gens} |H|={len(ch.H)} O={ch.st.O} C={ch.st.C} rows=(1,{pat}) d1={d1} d2={d2}: "
              f"joint states {[len(s[0]) for s in S]}, level-1 vectors {[len(p[0]) for p in proj]}, level-1 solutions {len(sols1)} "
              f"[{time.time() - t0:.0f}s]", flush=True)
    # order branches by the number of sequences they contain (largest first)
    weight = [np.prod([proj[k][2][q[k]] for k in range(4)]) for q in sols1]
    branches = [sols1[i] for i in np.argsort(weight)[::-1]]
    truncated, found, combos_done, nsol, timed_out = 0, [], 0, 0, False
    halves = [None] * 4
    if max(ch.st.O, len(ch.pairs)) > 46:
        print('  too many orbits for the meet-in-the-middle expansion', flush=True); return dict(status='too-large')
    for bi, q in enumerate(branches):
        if time.time() - t0 > budget:
            truncated += 1; timed_out = True; break
        if maxbranches is not None and bi >= maxbranches:          # breadth-first passes: heaviest branches only
            truncated += 1; break
        if d2:
            fib = []
            for k in range(4):
                rows = np.flatnonzero(proj[k][1] == q[k])
                Kf, cf = S[k][0][rows], S[k][1][rows]
                fib.append((Kf, cf))
            a2, b2 = ch.splits[L2], ch.splits[L2 + 1]
            P2 = [ch.paf(L2, Kf[:, a2:b2].astype(np.int64)) for Kf, _ in fib]
            try:
                sols2 = join4(P2)
            except JoinTooLarge:
                truncated += 1; continue
            cands = [[fib[k][0][s[k]] for k in range(4)] for s in sols2]
            counts = [[fib[k][1][s[k]] for k in range(4)] for s in sols2]
            # largest fibres first: the expected yield of a final join grows like the product of the four fibre
            # sizes, its work like the larger pair product, so yield per unit of work favours big combinations
            order = sorted(range(len(cands)), key=lambda i: -float(np.prod([float(c) for c in counts[i]])))
            cands = [cands[i] for i in order]; counts = [counts[i] for i in order]
        else:
            cands = [[S[k][0][np.flatnonzero(proj[k][1] == q[k])[0]] for k in range(4)]]
            counts = [[proj[k][2][q[k]] for k in range(4)]]
        for keys, cnt in zip(cands, counts):
            if time.time() - t0 > budget:                 # the budget also bounds the work inside a branch
                truncated += 1; timed_out = True; break
            if max(cnt) > maxfibre:
                truncated += 1; continue
            lists = []
            for k in range(4):
                if halves[k] is None: halves[k] = Halves(ch, roles[k])
                Y = halves[k].fibre(keys[k], maxfibre)
                if Y is None: lists = None; break
                lists.append(Y)
            if lists is None: truncated += 1; continue
            combos_done += 1
            if any(len(Y) == 0 for Y in lists): continue
            try:
                total, sols = split_join(ch, lists, extra) if extra else final_join(ch, lists)
            except JoinTimeout:
                truncated += 1; continue
            nsol += total
            for Ys in sols:
                seqs = verify_sol(ch, Ys)
                if seqs is not None:
                    rec = dict(n=v, gens=list(gens), pattern=[int(x.sum()) for x in seqs], hadamard=True, skew=skew,
                               method=f'chain d1={d1} d2={d2}', seqs=[''.join('+' if x > 0 else '-' for x in s) for s in seqs])
                    with open(W.hits() if skew else W.hits('had_hits.jsonl'), 'a') as f: f.write(json.dumps(rec) + '\n')
                    found.append(rec)
            if found and not exhaustive: break
        if verbose and (bi < 5 or bi % 50 == 0 or found):
            print(f"  branch {bi + 1}/{len(branches)}: combos done {combos_done}, truncated {truncated}, verified {len(found)} "
                  f"[{time.time() - t0:.0f}s]", flush=True)
        if found and not exhaustive: break
    status = 'found' if found else ('ruled out' if truncated == 0 else
              (f'nothing found (time budget of {budget}s reached)' if timed_out else
               (f'nothing found in the first {maxbranches} of {len(branches)} branches' if maxbranches is not None and len(branches) > maxbranches
                else f'nothing found ({truncated} combinations truncated)')))
    print(f"  RESULT v={v} H={gens} rows=(1,{pat}): {status}; level-1 branches {len(branches)}, combos {combos_done}, solutions counted {nsol} "
          f"[{time.time() - t0:.0f}s]", flush=True)
    return dict(status=status, found=found, branches=len(branches), truncated=truncated, solutions=nsol)


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('v', type=int); ap.add_argument('gens'); ap.add_argument('pattern')
    ap.add_argument('--d1', type=int); ap.add_argument('--d2', type=int)
    ap.add_argument('--exhaustive', action='store_true')
    ap.add_argument('--maxfibre', type=float, default=3e6); ap.add_argument('--maxstates', type=float, default=2e7)
    ap.add_argument('--budget', type=float, default=3600); ap.add_argument('--maxbranches', type=int)
    ap.add_argument('--join-timeout', type=float)
    ap.add_argument('--split', action='store_true', help='group the final join by the compressions the joint DP left out')
    a = ap.parse_args()
    v = a.v; ps = [p for p in range(2, v) if v % p == 0 and all(p % q for q in range(2, int(p ** 0.5) + 1))]
    d1 = a.d1 or ps[0]; d2 = a.d2 if a.d2 is not None else (ps[1] if len(ps) > 1 else 0)
    run(v, [int(x) for x in a.gens.split(',')], tuple(int(x) for x in a.pattern.split(',')), d1, d2, a.exhaustive,
        int(a.maxfibre), int(a.maxstates), budget=a.budget, maxbranches=a.maxbranches, join_timeout=a.join_timeout,
        split=a.split)
