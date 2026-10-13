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

  python3 attack/cyclo_product.py v --q Q --e E [--plain] [--config I] [--K N] [--sub G] [--maxquads N]
      --plain      ordinary Goethals-Seidel quadruple (no skew block)
      --config I   which stage-1 class to search (default 0; --list prints them)
      --sub G      keep only the G largest level-B groups of list 0 (partial search, for plain structures)
  Hits go to $HADAMARD_WORK/skew_hits.jsonl (skew) or had_hits.jsonl (plain).
"""
import argparse, itertools, json, math, os, struct, subprocess, sys, tempfile, time
import numpy as np
from pathlib import Path as _P
sys.path[:0] = [str(_P(__file__).resolve().parent), str(_P(__file__).resolve().parent.parent / 'verify')]
import workdir as W
import verify


def primitive_root(p):
    t = p - 1; fs = [f for f in range(2, t + 1) if t % f == 0 and all(f % d for d in range(2, int(f ** 0.5) + 1))]
    return next(r for r in range(2, p) if all(pow(r, t // f, p) != 1 for f in fs))


class Structure:
    def __init__(self, v, q, e, skew=True):
        assert v % q == 0 and (v // q) % q, "q must divide v exactly once"
        self.v, self.q, self.e, self.m, self.skew = v, q, e, v // q, skew
        assert (q - 1) % e == 0 and ((q - 1) // e) % 2 == 1, "need e | q - 1 with (q - 1)/e odd"
        self.F = (q - 1) // e
        g = primitive_root(q); self.g = g
        ind = np.zeros(q, dtype=np.int64); x = 1
        for k in range(q - 1): ind[x] = k; x = x * g % q
        self.cls = np.full(q, -1, dtype=np.int64); self.cls[1:] = ind[1:] % e
        self.half = int(self.cls[q - 1])                      # class of -1, = e/2
        assert self.half == e // 2
        ge = pow(g, e, q)
        self.h = next(t for t in range(1, v) if t % self.m == 1 and t % q == ge)
        self.reps = [0] + [int(np.flatnonzero(self.cls == c)[0]) for c in range(e)]

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
def stage1(S, maxpairs=5e7):
    """All quadruples of compressions (one class per orbit of negation, translation and reversal of the non-skew
    blocks), as tuples of 4 y-vectors."""
    m, e, F = S.m, S.e, S.F
    T = 4 * m * e
    cost = lambda s, f: F * s * s + 2 * f * s
    vals = [(1, 0), (-1, 0)] + [(f, s) for s in range(-e, e + 1, 2) if s for f in (1, -1) if cost(s, f) <= T]
    def options(skewblock):
        Y, C = [], []
        if skewblock:
            for half in itertools.product(vals, repeat=m // 2):
                fs = [(1, 0)] + list(half) + [(-f, -s) for f, s in reversed(half)]
                c = sum(cost(s, f) for f, s in fs)
                if c <= T: Y.append([f + F * s for f, s in fs]); C.append(c)
        else:
            for fs in itertools.product(vals, repeat=m):
                c = sum(cost(s, f) for f, s in fs)
                if c <= T: Y.append([f + F * s for f, s in fs]); C.append(c)
        Y = np.array(Y, dtype=np.int64).reshape(-1, m); C = np.array(C, dtype=np.int64)
        P = np.stack([(Y * np.roll(Y, -b, axis=1)).sum(1) for b in range(1, m // 2 + 1)], axis=1) if m > 1 else np.zeros((len(Y), 0), np.int64)
        return Y, C, P
    L = [options(S.skew and k == 0) for k in range(4)]
    # canonical id of every option under negation, translation and reversal
    def canon_ids(Y):
        best = None
        for sgn in (1, -1):
            for r in (1, -1):
                for t in range(m):
                    Z = sgn * Y[:, (r * np.arange(m) + t) % m]
                    if best is None: best = Z.copy(); continue
                    lt = np.zeros(len(Y), bool); eq = np.ones(len(Y), bool)
                    for a in range(m):
                        lt |= eq & (Z[:, a] < best[:, a]); eq &= Z[:, a] == best[:, a]
                    best[lt] = Z[lt]
        return best
    R = np.random.default_rng(5).integers(1, 2 ** 61, size=1 + m // 2)
    def key(c, P): return c * R[0] + (P * R[1:]).sum(1)
    pairs = []
    for (i, j) in ((0, 1), (2, 3)):
        Ci, Cj = L[i][1], L[j][1]
        I, J = np.nonzero(Ci[:, None] + Cj[None, :] <= T)
        if len(I) > maxpairs: raise RuntimeError(f'stage 1: {len(I)} pairs')
        pairs.append((I, J, Ci[I] + Cj[J], L[i][2][I] + L[j][2][J]))
    I01, J01, c01, p01 = pairs[0]; I23, J23, c23, p23 = pairs[1]
    k01 = key(c01, p01); o = np.argsort(k01); ks = k01[o]
    need = key(T - c23, -p23)
    lo = np.searchsorted(ks, need, 'left'); hi = np.searchsorted(ks, need, 'right')
    sel = np.flatnonzero(hi > lo)
    quads = []
    for r in sel:
        for pos in range(lo[r], hi[r]):
            a = o[pos]
            quads.append((I01[a], J01[a], I23[r], J23[r]))
    if not quads: return []
    Q = np.array(quads, dtype=np.int64)
    canon = [canon_ids(L[k][0]) for k in range(4)]
    keys = {}
    for row in Q:
        cs = [tuple(canon[k][row[k]]) for k in range(4)]
        kk = (tuple(L[0][0][row[0]]),) + tuple(sorted(cs[1:])) if S.skew else tuple(sorted(cs))
        if kk not in keys: keys[kk] = tuple(tuple(int(x) for x in L[k][0][row[k]]) for k in range(4))
    return [keys[k] for k in sorted(keys)]


# ---- lists and levels -----------------------------------------------------------------------------------------
def fs_of(S, yv):
    s = int(round(yv / S.F)); f = yv - S.F * s
    assert f in (1, -1)
    return f, s


def members(S, y, skewblock):
    """Array [n x m] of pattern ids, position by position."""
    m = S.m
    if skewblock:
        half = [S.patterns(*fs_of(S, y[a])) for a in range(1, m // 2 + 1)]
        rows = []
        for p0 in S.skew0():
            for combo in itertools.product(*half):
                row = [p0] + list(combo) + [S.partner(p) for p in reversed(combo)]
                rows.append(row)
        return np.array(rows, dtype=np.int64)
    cols = [np.array(S.patterns(*fs_of(S, y[a])), dtype=np.int64) for a in range(m)]
    grids = np.meshgrid(*cols, indexing='ij')
    return np.stack([g.ravel() for g in grids], axis=1)


def level_vectors(S, Pm, idx, CC, chunk=1 << 21):
    """Level-B (e/2) and level-C ((m-1)/2 * e) vectors of the members Pm [n x m] (pattern ids -> rows of CC via idx)."""
    m, e, h = S.m, S.e, S.e // 2
    n = len(Pm); B = np.empty((n, h), dtype=np.int16); C = np.empty((n, (m // 2) * e), dtype=np.int16)
    for c0 in range(0, n, chunk):
        X = idx[Pm[c0:c0 + chunk]]
        paf = [sum(CC[X[:, (a + b) % m], X[:, a], :].astype(np.int32) for a in range(m)) for b in range(m)]
        lb = sum(paf[b][:, 1:1 + h] for b in range(m))                 # sum over Z_m: classes C_0 .. C_{h-1}
        B[c0:c0 + chunk] = lb
        C[c0:c0 + chunk] = np.concatenate([paf[b][:, 1:1 + e] for b in range(1, m // 2 + 1)], axis=1)
    return B, C


def grouped(B, C):
    key = np.ascontiguousarray(B).view(np.dtype((np.void, 2 * B.shape[1]))).ravel()
    order = np.argsort(key, kind='stable'); ks = key[order]
    brk = np.flatnonzero(np.r_[True, ks[1:] != ks[:-1]])
    offs = np.r_[brk, len(order)].astype(np.int64)
    return B[order[brk]], offs, C[order], order.astype(np.int64)


def run(v, q, e, plain=False, config=0, K=5, sub=0, maxquads=-1, threads=None, listonly=False):
    t0 = time.time(); log = lambda s: print(s, flush=True)
    S = Structure(v, q, e, skew=not plain)
    quads = stage1(S)
    log(f"v={v} = {S.m}*{q}, e={e}, F={S.F}, H=<{S.h}> of order {S.F}, {'plain' if plain else 'skew'}: "
        f"{len(quads)} stage-1 classes [{time.time() - t0:.0f}s]")
    for i, qd in enumerate(quads):
        rs = [int(sum(y)) for y in qd]
        log(f"  config {i}: y = {qd}, row sums {rs}")
    if listonly or not quads: return []
    qd = quads[config]
    Pms = [members(S, qd[k], S.skew and k == 0) for k in range(4)]
    used = sorted(set(np.concatenate([P.ravel() for P in Pms]).tolist()))
    idx = np.full(1 << (e + 1), -1, dtype=np.int64); idx[used] = np.arange(len(used))
    CC = S.cc_table(used)
    log(f"  config {config}: list sizes {[len(P) for P in Pms]}, {len(used)} slice patterns [{time.time() - t0:.0f}s]")
    G = []
    for k in range(4):
        B, C = level_vectors(S, Pms[k], idx, CC)
        G.append(grouped(B, C))
    order = sorted(range(4), key=lambda k: len(G[k][0]))            # fewest level-B groups first (hashed pair)
    if sub:
        Bg, offs, Cs, ids = G[order[0]]
        sz = np.diff(offs); pick = np.sort(np.argsort(-sz)[:sub])
        rows = np.concatenate([np.arange(offs[g], offs[g + 1]) for g in pick])
        G[order[0]] = (Bg[pick], np.r_[0, np.cumsum(sz[pick])].astype(np.int64), Cs[rows], ids[rows])
    log(f"  level-B groups {[len(G[k][0]) for k in range(4)]}, join order {order} [{time.time() - t0:.0f}s]")
    fd, path = tempfile.mkstemp(suffix='.bin', dir=W.tmpdir()); os.close(fd)
    try:
        with open(path, 'wb') as f:
            f.write(struct.pack('ii', S.e // 2, (S.m // 2) * S.e))
            for k in order:
                Bg, offs, Cs, ids = G[k]
                f.write(struct.pack('iiq', len(Bg), 0, len(ids))); f.write(np.ascontiguousarray(Bg).tobytes())
                f.write(offs.tobytes()); f.write(np.ascontiguousarray(Cs).tobytes()); f.write(ids.tobytes())
        del G
        env = dict(os.environ)
        if threads: env['OMP_NUM_THREADS'] = str(threads)
        p = subprocess.Popen([W.bin('join3'), path, str(K), str(maxquads)], stdout=subprocess.PIPE, text=True, env=env)
        found = []
        for line in p.stdout:
            if line.startswith('SOL'):
                ids = list(map(int, line.split()[1:]))
                mem = [None] * 4
                for slot, k in enumerate(order): mem[k] = ids[slot]
                seqs = [S.sequence(Pms[k][mem[k]]) for k in range(4)]
                tot = verify.paf_sum(seqs)
                ok = tot[0] == 4 * v and not tot[1:].any()
                if S.skew: ok = ok and seqs[0][0] == 1 and all(seqs[0][(-t) % v] == -seqs[0][t] for t in range(1, v))
                log(f"  hit: autocorrelations {'ok' if ok else 'WRONG'}, row sums {[int(x.sum()) for x in seqs]} [{time.time() - t0:.0f}s]")
                if ok:
                    rec = dict(n=v, gens=[S.h], pattern=[int(x.sum()) for x in seqs], hadamard=True, skew=S.skew,
                               method=f'cyclo_product q={q} e={e} config={config}',
                               seqs=[''.join('+' if z > 0 else '-' for z in x) for x in seqs])
                    with open(W.hits() if S.skew else W.hits('had_hits.jsonl'), 'a') as fh: fh.write(json.dumps(rec) + '\n')
                    found.append(rec)
            elif line.startswith('TOTAL'):
                log(f"  join3: {line.strip()} [{time.time() - t0:.0f}s]")
        p.wait()
    finally:
        os.remove(path)
    return found


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('v', type=int); ap.add_argument('--q', type=int, required=True); ap.add_argument('--e', type=int, required=True)
    ap.add_argument('--plain', action='store_true'); ap.add_argument('--config', type=int, default=0)
    ap.add_argument('--K', type=int, default=5); ap.add_argument('--sub', type=int, default=0)
    ap.add_argument('--maxquads', type=int, default=-1); ap.add_argument('--threads', type=int)
    ap.add_argument('--list', action='store_true')
    a = ap.parse_args()
    run(a.v, a.q, a.e, a.plain, a.config, a.K, a.sub, a.maxquads, a.threads, a.list)
