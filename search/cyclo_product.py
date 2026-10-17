"""Cyclotomic-product search for Goethals-Seidel families on Z_v, v = m q, skew or plain.

**Structure.** q is a prime not dividing m, e divides q - 1 and F = (q - 1)/e is odd, so that -1 is not an e-th
power. H = {h : h = 1 (mod m), h mod q an e-th power} has order F and -1 is not in H. Over Z_m x GF(q) an
H-invariant block is given, at each a in Z_m, by its value f(a) at (a, 0) and its values x_0(a) .. x_{e-1}(a) on
(a, C_j), where C_j = g^j <g^e> are the cyclotomic classes of index e (g a primitive root). For e = 2 this is the
Paley product of paley_config.py.

**Compression identity.** The compression of a block to Z_m is y(a) = f(a) + F s(a) with s(a) = sum_j x_j(a).
Summing the squares of the four compressions gives 4mq, which reads

    sum over the 4m positions of  F s^2 + 2 f s  =  4 m e.

A position with s != 0 is a spike. For q large compared with m only a few spikes fit; for v = 933 = 3 * 311 with
e = 10 (F = 31) exactly one spike fits, with s = +-2 and f = -sign(s), so y = +-61.

**Search.**
1. Stage 1 lists every quadruple of compressions (spikes included) whose autocorrelations sum to zero on Z_m,
   up to negation, translation and reversal of the blocks (block 0 stays fixed in the skew case).
2. For a chosen quadruple each block is enumerated position by position (class patterns with the prescribed f and
   s; in the skew case position 0 of block 0 is of skew type and the positions a, -a are paired).
3. Level B: the compression to GF(q) (sum over Z_m) gives, per block, the autocorrelations on the classes C_t,
   t < e/2. Level C: the autocorrelations at the shifts (b, C_c), 1 <= b <= (m-1)/2. Stage 1 + B + C is
   equivalent to the full Goethals-Seidel condition.
4. join3.c joins the four lists exactly: level-B group quadruples by hash, then the level-C vectors inside each
   quadruple. Each hit is rebuilt on Z_v and its autocorrelations are checked before it is recorded.

  python3 search/cyclo_product.py v --q Q --e E [--plain] [--kgen k] [--config I] [--K N] [--sub G] [--maxquads N]
      --plain      ordinary Goethals-Seidel quadruple (no skew block)
      --kgen k     also impose invariance under the multipliers <k> of Z_m (acting on the Z_m coordinate)
      --config I   which stage-1 class to search (default 0; --list prints them)
      --sub G      keep only the G largest level-B groups of every list (partial search)
      --bkey b,..  add the autocorrelations at the shifts (b, C_c) to the group key (finer groups)
  Hits go to $HADAMARD_WORK/skew_hits.jsonl (skew) or had_hits.jsonl (plain).
"""
import argparse, itertools, json, os, struct, subprocess, sys, tempfile, time
import numpy as np
from pathlib import Path as _P
sys.path[:0] = [str(_P(__file__).resolve().parent), str(_P(__file__).resolve().parent.parent / 'verify')]
import workdir as W
import verify


def primitive_root(p):
    t = p - 1; fs = [f for f in range(2, t + 1) if t % f == 0 and all(f % d for d in range(2, int(f ** 0.5) + 1))]
    return next(r for r in range(2, p) if all(pow(r, t // f, p) != 1 for f in fs))


class Structure:
    def __init__(self, v, q, e, skew=True, kgen=1):
        assert v % q == 0 and (v // q) % q, "q must divide v exactly once"
        self.v, self.q, self.e, self.m, self.skew = v, q, e, v // q, skew
        assert (q - 1) % e == 0, "need e | q - 1"
        assert not skew or ((q - 1) // e) % 2 == 1, "skew needs (q - 1)/e odd (-1 not an e-th power)"
        self.F = (q - 1) // e
        g = primitive_root(q); self.g = g
        ind = np.zeros(q, dtype=np.int64); x = 1
        for k in range(q - 1): ind[x] = k; x = x * g % q
        self.cls = np.full(q, -1, dtype=np.int64); self.cls[1:] = ind[1:] % e
        assert e % 2 == 0, "e must be even"
        self.half = int(self.cls[q - 1])                      # class of -1: e/2 if (q-1)/e is odd, else 0
        self.nb = e // 2 if self.half else e                  # level-B dimensions: classes up to negation
        ge = pow(g, e, q)
        self.h = next(t for t in range(1, v) if t % self.m == 1 % self.m and t % q == ge)
        self.reps = [0] + [int(np.flatnonzero(self.cls == c)[0]) for c in range(e)]
        # optional multiplier subgroup K_m = <kgen> of Z_m^*, acting on the Z_m coordinate
        m = self.m; self.kgen = kgen % m if m > 1 else 1
        K = {1}
        while m > 1:
            K2 = K | {k * self.kgen % m for k in K}
            if K2 == K: break
            K = K2
        self.K = sorted(K)
        self.korb = np.full(m, -1, dtype=np.int64); self.kreps = []
        for a in range(m):
            if self.korb[a] < 0:
                for k in self.K: self.korb[a * k % m] = len(self.kreps)
                self.kreps.append(a)
        if skew and m > 1:
            assert all(self.korb[(-a) % m] != self.korb[a] for a in range(1, m)), "-1 in K_m: no skew block"
        self.gens = [self.h] + ([next(t for t in range(1, v) if t % m == self.kgen and t % q == 1)] if len(self.K) > 1 else [])

    # ---- slice patterns: id = (fbit << e) | mask, f = 1 - 2 fbit, x_j = 1 - 2 bit_j ------------------------------
    def pattern(self, pid):
        e = self.e
        return 1 - 2 * ((pid >> e) & 1), [1 - 2 * ((pid >> j) & 1) for j in range(e)]

    def patterns(self, f, s):
        e = self.e; fb = 0 if f == 1 else 1
        return [(fb << e) | mk for mk in range(1 << e) if e - 2 * bin(mk).count('1') == s]

    def skew0(self):
        """Position 0 of the skew block: f = 1 and x_{j + e/2} = -x_j."""
        e, h = self.e, self.e // 2; out = []
        for mk in range(1 << h):
            xs = [1 - 2 * ((mk >> j) & 1) for j in range(h)]; xs = xs + [-t for t in xs]
            out.append(sum(1 << j for j in range(e) if xs[j] < 0))
        return out

    def partner(self, pid):
        """Pattern at -a forced by the pattern at a in the skew block: f -> -f, x_j -> -x_{j + e/2}."""
        f, xs = self.pattern(pid); e, h = self.e, self.e // 2
        ys = [-xs[(j + h) % e] for j in range(e)]
        return ((0 if -f == 1 else 1) << e) | sum(1 << j for j in range(e) if ys[j] < 0)

    def cc_table(self, used):
        """CC[i, j, c] = sum_y U_i(y + x_c) U_j(y) on GF(q), x_0 = 0 and x_{1+c} a representative of C_c."""
        q = self.q
        V = np.empty((len(used), q))
        for i, pid in enumerate(used):
            f, xs = self.pattern(pid); V[i, 0] = f; V[i, 1:] = np.asarray(xs)[self.cls[1:]]
        Fv = np.fft.rfft(V, axis=1); out = np.empty((len(used), len(used), len(self.reps)), dtype=np.int16)
        for i in range(len(used)):
            R = np.fft.irfft(Fv[i][None, :] * np.conj(Fv), n=q, axis=1)[:, self.reps]
            out[i] = np.rint(R).astype(np.int16)
        return out

    def sequence(self, slices):
        """+-1 sequence on Z_v from the pattern ids at the positions a = 0 .. m-1."""
        v, m, q = self.v, self.m, self.q; t = np.arange(v)
        P = np.array([[self.pattern(p)[0]] + self.pattern(p)[1] for p in slices])     # m x (e + 1)
        col = np.where(t % q == 0, 0, 1 + self.cls[t % q])
        return P[t % m, col]


# ---- stage 1: compressions on Z_m ---------------------------------------------------------------------------
def stage1_options(S, skewblock, maxopts=None):
    """Compressions y on Z_m of one block, constant on the K_m-orbits (skew block: f = 1, s = 0 at 0 and odd under
    negation), with their spike costs C (sum over the 4 blocks = 4 m e) and autocorrelations P at the shifts
    1 .. m/2. With maxopts, a RuntimeError is raised before generating more options than that."""
    m, e, F = S.m, S.e, S.F
    T = 4 * m * e
    cost = lambda s, f: F * s * s + 2 * f * s
    vals = [(1, 0), (-1, 0)] + [(f, s) for s in range(-e, e + 1, 2) if s for f in (1, -1) if cost(s, f) <= T]
    spikes = [(f, s) for f, s in vals if s]
    orbits = [[a for a in range(m) if S.korb[a] == o] for o in range(len(S.kreps))]
    def placements(units, budget):
        """(unit indices, (f, s) per spike) with total cost <= budget; units[i] = cost multiplier of unit i."""
        out = [((), ())]; tot = [1 << len(units)]
        def rec(start, used, cur, c):
            for i in range(start, len(units)):
                for f, s in spikes:
                    cc = c + units[i] * cost(s, f)
                    if cc > budget: continue
                    out.append((used + (i,), cur + ((f, s),)))
                    tot[0] += 1 << (len(units) - len(used) - 1)
                    if maxopts and tot[0] > maxopts: raise RuntimeError(f'stage 1: more than {maxopts:.0f} options per block')
                    rec(i + 1, used + (i,), cur + ((f, s),), cc)
        rec(0, (), (), 0)
        return out
    if skewblock:
        free, seen = [], {S.korb[0]}
        for o, orb in enumerate(orbits):
            if o in seen: continue
            neg = S.korb[(-orb[0]) % m]
            free.append((orb, orbits[neg])); seen |= {o, neg}
        units = [2 * len(orb) for orb, _ in free]
    else:
        free = [(orb, None) for orb in orbits]
        units = [len(orb) for orb, _ in free]
    Y, C = [], []
    for pos, fsl in placements(units, T):
        rest = [i for i in range(len(free)) if i not in pos]
        c0 = sum(units[i] * cost(s, f) for i, (f, s) in zip(pos, fsl))
        for signs in itertools.product((1, -1), repeat=len(rest)):
            y = [0] * m
            if skewblock: y[0] = 1
            assign = dict(zip(pos, fsl)); assign.update({i: (sg, 0) for i, sg in zip(rest, signs)})
            for i, (f, s) in assign.items():
                orb, neg = free[i]
                for a in orb: y[a] = f + F * s
                if neg is not None:
                    for a in neg: y[a] = -(f + F * s)
            Y.append(y); C.append(c0)
    Y = np.array(Y, dtype=np.int64).reshape(-1, m); C = np.array(C, dtype=np.int64)
    P = np.stack([(Y * np.roll(Y, -b, axis=1)).sum(1) for b in range(1, m // 2 + 1)], axis=1) if m > 1 else np.zeros((len(Y), 0), np.int64)
    return Y, C, P


def stage1(S, maxpairs=3e7, maxopts=3e6):
    """All quadruples of compressions (one class per orbit of negation, translation and reversal of the non-skew
    blocks), as tuples of 4 y-vectors."""
    m, e = S.m, S.e
    T = 4 * m * e
    options = lambda skewblock: stage1_options(S, skewblock)
    L = [options(S.skew and k == 0) for k in range(4)]
    if max(len(x[0]) for x in L) > maxopts: raise RuntimeError(f'stage 1: {max(len(x[0]) for x in L)} options per block')
    def canon_ids(Y):
        """Integer id of the class of each row under negation, translation and reversal."""
        best = None
        for sgn in (1, -1):
            for r in (1, -1):
                for t in (range(m) if len(S.K) == 1 else [0]):     # translations break K_m-invariance
                    Z = sgn * Y[:, (r * np.arange(m) + t) % m]
                    if best is None: best = Z.copy(); continue
                    lt = np.zeros(len(Y), bool); eq = np.ones(len(Y), bool)
                    for a in range(m):
                        lt |= eq & (Z[:, a] < best[:, a]); eq &= Z[:, a] == best[:, a]
                    best[lt] = Z[lt]
        u, inv = np.unique(best, axis=0, return_inverse=True)
        return inv.ravel()
    cid = [canon_ids(L[k][0]) for k in range(4)]
    if S.skew: cid[0] = np.arange(len(L[0][0]))                     # the skew block is not reduced
    R = np.random.default_rng(5).integers(1, 2 ** 61, size=max(1, m // 2))
    hk = [(L[k][2] * R[:L[k][2].shape[1]]).sum(1) for k in range(4)]
    buckets = [{int(c): np.flatnonzero(L[k][1] == c) for c in np.unique(L[k][1])} for k in range(4)]
    found = {}
    for c0, i0 in buckets[0].items():
        for c1, i1 in buckets[1].items():
            for c2, i2 in buckets[2].items():
                c3 = T - c0 - c1 - c2
                if c3 not in buckets[3]: continue
                i3 = buckets[3][c3]
                if len(i0) * len(i1) > maxpairs or len(i2) * len(i3) > maxpairs:
                    raise RuntimeError(f'stage 1: {len(i0) * len(i1)} x {len(i2) * len(i3)} pairs')
                A, B = np.meshgrid(i0, i1, indexing='ij'); A, B = A.ravel(), B.ravel()
                k01 = hk[0][A] + hk[1][B]; o = np.argsort(k01); ks = k01[o]
                Cg, Dg = np.meshgrid(i2, i3, indexing='ij'); Cg, Dg = Cg.ravel(), Dg.ravel()
                need = -(hk[2][Cg] + hk[3][Dg])
                lo = np.searchsorted(ks, need, 'left'); hi = np.searchsorted(ks, need, 'right')
                for r in np.flatnonzero(hi > lo):
                    for pos in range(lo[r], hi[r]):
                        a, b, c, d = A[o[pos]], B[o[pos]], Cg[r], Dg[r]
                        if np.any(L[0][2][a] + L[1][2][b] + L[2][2][c] + L[3][2][d]): continue
                        ids = (cid[0][a], cid[1][b], cid[2][c], cid[3][d])
                        key = (ids[0],) + tuple(sorted(ids[1:])) if S.skew else tuple(sorted(ids))
                        if key not in found:
                            found[key] = tuple(tuple(int(x) for x in L[k][0][j]) for k, j in enumerate((a, b, c, d)))
    return [found[k] for k in sorted(found)]


# ---- lists and levels -----------------------------------------------------------------------------------------
def fs_of(S, yv):
    s = int(round(yv / S.F)); f = yv - S.F * s
    assert f in (1, -1)
    return f, s


def members(S, y, skewblock):
    """Array [n x m] of pattern ids, position by position (constant on the K_m-orbits of positions)."""
    m = S.m; R = S.kreps; orb = S.korb
    if skewblock:
        pairs = []                                     # orbit representatives r with -r in a different orbit, r first
        seen = {0}
        for r in R:
            if orb[r] in seen: continue
            pairs.append(r); seen.add(orb[r]); seen.add(orb[(-r) % m])
        half = [S.patterns(*fs_of(S, y[r])) for r in pairs]
        rows = []
        for p0 in S.skew0():
            for combo in itertools.product(*half):
                row = [0] * m
                for a in range(m):
                    if orb[a] == orb[0]: row[a] = p0
                for r, p in zip(pairs, combo):
                    pp = S.partner(p)
                    for a in range(m):
                        if orb[a] == orb[r]: row[a] = p
                        elif orb[a] == orb[(-r) % m]: row[a] = pp
                rows.append(row)
        return np.array(rows, dtype=np.int64)
    cols = [np.array(S.patterns(*fs_of(S, y[r])), dtype=np.int16) for r in R]
    grids = np.meshgrid(*cols, indexing='ij')
    G = np.stack([g.ravel() for g in grids], axis=1)
    return G[:, orb]


def level_vectors(S, Pm, idx, CC, chunk=1 << 21):
    """Level-B and level-C vectors of the members Pm [n x m] (pattern ids -> rows of CC via idx). Level B is the
    GF(q)-compression (S.nb classes) plus, for each shift b in S.bkey, the autocorrelations at (b, C_c)."""
    m, e, h = S.m, S.e, S.nb
    n = len(Pm); B = np.empty((n, S.dkey), dtype=np.int16); C = np.empty((n, (m // 2) * e), dtype=np.int16)
    for c0 in range(0, n, chunk):
        X = idx[Pm[c0:c0 + chunk]]
        paf = [sum(CC[X[:, (a + b) % m], X[:, a], :].astype(np.int32) for a in range(m)) for b in range(m)]
        lb = sum(paf[b][:, 1:1 + h] for b in range(m))                 # sum over Z_m: classes C_0 .. C_{h-1}
        B[c0:c0 + chunk] = np.concatenate([lb] + [paf[b][:, 1:1 + e] for b in S.bkey], axis=1) if S.bkey else lb
        C[c0:c0 + chunk] = np.concatenate([paf[b][:, 1:1 + e] for b in range(1, m // 2 + 1)], axis=1)
    return B, C


def grouped(B, C):
    key = np.ascontiguousarray(B).view(np.dtype((np.void, 2 * B.shape[1]))).ravel()
    order = np.argsort(key, kind='stable'); ks = key[order]
    brk = np.flatnonzero(np.r_[True, ks[1:] != ks[:-1]])
    offs = np.r_[brk, len(order)].astype(np.int64)
    return B[order[brk]], offs, C[order], order.astype(np.int64)


def run(v, q, e, plain=False, config=0, K=5, sub=0, maxquads=-1, threads=None, listonly=False, kgen=1, bkey=()):
    t0 = time.time(); log = lambda s: print(s, flush=True)
    S = Structure(v, q, e, skew=not plain, kgen=kgen)
    S.bkey = [int(b) for b in bkey]; S.dkey = S.nb + S.e * len(S.bkey)
    quads = stage1(S)
    log(f"v={v} = {S.m}*{q}, e={e}, F={S.F}, H=<{', '.join(map(str, S.gens))}> of order {S.F * len(S.K)}, {'plain' if plain else 'skew'}: "
        f"{len(quads)} stage-1 classes [{time.time() - t0:.0f}s]")
    for i, qd in enumerate(quads):
        rs = [int(sum(y)) for y in qd]
        log(f"  config {i}: y = {qd}, row sums {rs}")
    if listonly or not quads: return []
    qd = quads[config]
    Pms = [members(S, qd[k], S.skew and k == 0).astype(np.int16) for k in range(4)]   # pattern ids < 2^(e+1)
    used = sorted(set(np.concatenate([np.unique(P) for P in Pms]).tolist()))
    idx = np.full(1 << (e + 1), -1, dtype=np.int64); idx[used] = np.arange(len(used))
    CC = S.cc_table(used)
    log(f"  config {config}: list sizes {[len(P) for P in Pms]}, {len(used)} slice patterns [{time.time() - t0:.0f}s]")
    parts, ngroups = [], []
    try:
        for k in range(4):                         # one block at a time, grouped and written to its own file
            B, C = level_vectors(S, Pms[k], idx, CC)
            Bg, offs, Cs, ids = grouped(B, C); del B, C
            fd, pk = tempfile.mkstemp(suffix=f'.b{k}', dir=W.tmpdir()); os.close(fd)
            with open(pk, 'wb') as f:
                f.write(struct.pack('iiq', len(Bg), 0, len(ids))); f.write(np.ascontiguousarray(Bg).tobytes())
                f.write(offs.tobytes()); f.write(np.ascontiguousarray(Cs).tobytes()); f.write(ids.tobytes())
            parts.append(pk); ngroups.append(len(Bg)); del Bg, offs, Cs, ids
        order = sorted(range(4), key=lambda k: ngroups[k])           # fewest level-B groups first (hashed pair)
        if sub:                                    # keep the `sub` largest level-B groups of every list (partial search)
            dc = (S.m // 2) * S.e
            for k0 in range(4):
                with open(parts[k0], 'rb') as f:
                    nG, _, N = struct.unpack('iiq', f.read(16))
                    Bg = np.frombuffer(f.read(2 * S.dkey * nG), np.int16).reshape(nG, S.dkey)
                    offs = np.frombuffer(f.read(8 * (nG + 1)), np.int64)
                    Cs = np.frombuffer(f.read(2 * dc * N), np.int16).reshape(N, dc)
                    ids = np.frombuffer(f.read(8 * N), np.int64)
                if nG <= sub: continue
                sz = np.diff(offs); pick = np.sort(np.argsort(-sz, kind='stable')[:sub])
                rows = np.concatenate([np.arange(offs[g], offs[g + 1]) for g in pick])
                with open(parts[k0], 'wb') as f:
                    f.write(struct.pack('iiq', len(pick), 0, len(rows))); f.write(np.ascontiguousarray(Bg[pick]).tobytes())
                    f.write(np.r_[0, np.cumsum(sz[pick])].astype(np.int64).tobytes())
                    f.write(np.ascontiguousarray(Cs[rows]).tobytes()); f.write(ids[rows].tobytes())
                log(f"  list {k0}: kept the {len(pick)} largest groups, {len(rows)} of {N} members")
                ngroups[k0] = len(pick)
            order = sorted(range(4), key=lambda k: ngroups[k])
        log(f"  level-B groups {ngroups}, join order {order} [{time.time() - t0:.0f}s]")
        fd, path = tempfile.mkstemp(suffix='.bin', dir=W.tmpdir()); os.close(fd)
        parts.append(path)
        with open(path, 'wb') as f:
            f.write(struct.pack('ii', S.dkey, (S.m // 2) * S.e))
            for k in order:
                with open(parts[k], 'rb') as g:
                    while True:
                        buf = g.read(1 << 26)
                        if not buf: break
                        f.write(buf)
        for k in range(4): os.remove(parts[k])
        parts = [path]
        env = dict(os.environ)
        if threads: env['OMP_NUM_THREADS'] = str(threads)
        p = subprocess.Popen([W.bin('join3'), path, str(K), str(maxquads)], stdout=subprocess.PIPE, text=True, env=env)
        found = []
        for line in p.stdout:
            if line.startswith('SOL'):
                ids = list(map(int, line.split()[1:]))
                mem = [None] * 4
                for slot, k in enumerate(order): mem[k] = ids[slot]
                seqs = [S.sequence(Pms[k][mem[k]].astype(np.int64)) for k in range(4)]
                tot = verify.paf_sum(seqs)
                ok = tot[0] == 4 * v and not tot[1:].any()
                if S.skew: ok = ok and seqs[0][0] == 1 and all(seqs[0][(-t) % v] == -seqs[0][t] for t in range(1, v))
                log(f"  hit: autocorrelations {'ok' if ok else 'WRONG'}, row sums {[int(x.sum()) for x in seqs]} [{time.time() - t0:.0f}s]")
                if ok:
                    rec = dict(n=v, gens=S.gens, pattern=[int(x.sum()) for x in seqs], hadamard=True, skew=S.skew,
                               method=f'cyclo_product q={q} e={e} kgen={kgen} config={config}',
                               seqs=[''.join('+' if z > 0 else '-' for z in x) for x in seqs])
                    with open(W.hits() if S.skew else W.hits('had_hits.jsonl'), 'a') as fh: fh.write(json.dumps(rec) + '\n')
                    found.append(rec)
            elif line.startswith('TOTAL'):
                log(f"  join3: {line.strip()} [{time.time() - t0:.0f}s]")
        p.wait()
    finally:
        for pth in parts:
            if os.path.exists(pth): os.remove(pth)
    return found


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('v', type=int); ap.add_argument('--q', type=int, required=True); ap.add_argument('--e', type=int, required=True)
    ap.add_argument('--plain', action='store_true'); ap.add_argument('--config', type=int, default=0)
    ap.add_argument('--K', type=int, default=5); ap.add_argument('--sub', type=int, default=0)
    ap.add_argument('--maxquads', type=int, default=-1); ap.add_argument('--threads', type=int)
    ap.add_argument('--list', action='store_true'); ap.add_argument('--kgen', type=int, default=1)
    ap.add_argument('--bkey', default='', help='comma-separated shifts b whose autocorrelations (b, C_c) join the group key')
    a = ap.parse_args()
    run(a.v, a.q, a.e, a.plain, a.config, a.K, a.sub, a.maxquads, a.threads, a.list, a.kgen, [b for b in a.bkey.split(',') if b])
