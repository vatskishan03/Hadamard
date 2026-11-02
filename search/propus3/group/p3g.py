"""p3g -- front end that lets the p3 engine (p3.c, unchanged) search propus quadruples (A, B, B, D) over a finite abelian
group G = prod Z_p^a (elementary p-parts with a Singer cycle acting, Z_p parts with units acting) with a cyclic multiplier
group H = <phi>, phi in Aut' = prod C_{p^a - 1}.  A is constant on <H,-1>-orbits (symmetric), B and D on H-orbits.

  PAF_A + 2 PAF_B + PAF_D = 4 v delta  on G  (optionally = TGT: bordered problems, TGT = class values),
  array: propus_group.propus_matrix (group-developed blocks, R[g,h] = [g+h = 0]).

Everything the engine reads is derived here from spence_lib.Structure (orbits, class tensors S) plus
  * theta: |X^(chi)|^2 = v + sum_c theta[c'][c] F_c(X), chi over representatives of the orbits of <H^T,-1> on the dual group,
  * mult_reps: the automorphisms phi_e (e in Aut'/H) as permutations of G (lossless canonical form of B),
  * neg: the permutation g -> -g (lossless canonical form of D under reversal).
A hit is certified by direct group PAF on G and by the matrix test H = H^T, H H^T = 4 v I (not by the engine's arithmetic).
"""
import sys, os, json, itertools, math, time
from math import gcd
import numpy as np
REPO = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..')     # search/ when this directory is search/propus3/group
sys.path[:0] = [REPO, os.path.join(REPO, '..', 'verify')]
import spence_lib as L
import propus_group as PG
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import p3run as P


def factorint(n):
    f, p = {}, 2
    while p * p <= n:
        while n % p == 0:
            f[p] = f.get(p, 0) + 1; n //= p
        p += 1
    if n > 1: f[n] = f.get(n, 0) + 1
    return f


def matpow(M, e, p):
    R = np.eye(M.shape[0], dtype=np.int64); B = M.copy() % p
    while e:
        if e & 1: R = (R @ B) % p
        B = (B @ B) % p; e >>= 1
    return R


def singer(p, a):
    N = p ** a - 1; rs = list(factorint(N)); I = np.eye(a, dtype=np.int64)
    for coeffs in itertools.product(range(p), repeat=a):
        if coeffs[0] == 0: continue
        M = np.zeros((a, a), dtype=np.int64)
        for i in range(a - 1): M[i + 1, i] = 1
        for i in range(a): M[i, a - 1] = coeffs[i]
        if not np.array_equal(matpow(M, N, p), I): continue
        if all(not np.array_equal(matpow(M, N // r, p), I) for r in rs): return M
    raise ValueError


def primroot(p):
    if p == 2: return 1
    for g in range(2, p):
        if all(pow(g, (p - 1) // r, p) != 1 for r in factorint(p - 1)): return g
    return 1


def build_group(v):
    fac = factorint(v); parts = []
    for p, a in sorted(fac.items()):
        parts.append((p, a, np.array([[primroot(p)]], dtype=np.int64) if a == 1 else singer(p, a)))
    mods = []
    for p, a, M in parts: mods += [p] * a
    return L.AbGroup(mods), parts


def aut_matrix(G, parts, e):
    k = len(G.mods); Mf = np.eye(k, dtype=np.int64); off = 0
    for (p, a, M), ei in zip(parts, e):
        Mf[off:off + a, off:off + a] = matpow(M, ei, p) if ei else np.eye(a, dtype=np.int64)
        off += a
    return Mf


class GBlock:
    def __init__(self, st, v):
        self.O, self.orb, self.sizes = st.O, st.orb, np.asarray(st.sizes, dtype=np.int64)
        self.reps = [int(np.flatnonzero(st.orb == i)[0]) for i in range(st.O)]
        S = np.rint(st.S).astype(np.int64)
        C = S.shape[0]
        self.const0 = np.array([int(np.trace(S[c])) for c in range(C)], dtype=np.int64)
        T = S + S.transpose(0, 2, 1)
        for c in range(C): np.fill_diagonal(T[c], 0)
        self.T = T
        self.v = v
        self._W = None
    seq = lambda self, Y: np.asarray(Y)[..., self.orb]
    def _dp(self):
        v, O = self.v, self.O
        W = np.zeros((O + 1, 2 * v + 1)); W[O, v] = 1.0
        for j in range(O - 1, -1, -1):
            s = int(self.sizes[j])
            W[j, s:] += W[j + 1, :2 * v + 1 - s]; W[j, :2 * v + 1 - s] += W[j + 1, s:]
        return W
    def pool(self, r):
        if self._W is None: self._W = self._dp()
        return 0 if abs(r) > self.v else int(round(self._W[0, r + self.v]))


class Ambient:
    """Ambient automorphism group  prod_i GammaL(1, p_i^{a_i})  (a_i = 1: units of Z_p, no Frobenius), elements = tuples
    (a_1..a_r, b_1..b_r): S_i^{a_i} F_i^{b_i}  with  F S F^-1 = S^p.  Matrices act on G = prod Z_p^a."""
    def __init__(self, G, parts):
        self.G, self.parts = G, parts
        self.n = [p ** a - 1 for p, a, M in parts]           # order of S_i
        self.k = [a for p, a, M in parts]                    # order of F_i
        self.p = [p for p, a, M in parts]
        self.r = len(parts)
        # Frobenius matrices
        self.Fm = []
        for (p, a, M) in parts:
            if a == 1:
                self.Fm.append(np.array([[1]], dtype=np.int64)); continue
            e0 = np.zeros(a, dtype=np.int64); e0[0] = 1
            cols = [(matpow(M, (i * p) % (p ** a - 1), p) @ e0) % p for i in range(a)]
            self.Fm.append(np.array(cols, dtype=np.int64).T)
        self.ident = tuple([0] * (2 * self.r))
    def mul(self, x, y):
        r = self.r; out = []
        for i in range(r):
            out.append((x[i] + pow(self.p[i], x[r + i], self.n[i]) * y[i]) % self.n[i])
        for i in range(r):
            out.append((x[r + i] + y[r + i]) % self.k[i])
        return tuple(out)
    def inv(self, x):
        # brute: x^(order-1)
        y = x; z = self.ident
        # compute order
        cur = x; order = 1
        while cur != self.ident:
            cur = self.mul(cur, x); order += 1
        res = self.ident
        for _ in range(order - 1): res = self.mul(res, x)
        return res
    def closure(self, gens, limit=1 << 16):
        S = {self.ident}; frontier = [self.ident]
        while frontier:
            nxt = []
            for x in frontier:
                for g in gens:
                    y = self.mul(x, g)
                    if y not in S:
                        S.add(y); nxt.append(y)
                        if len(S) > limit: return None
            frontier = nxt
        return S
    def matrix(self, x):
        G = self.G; k = len(G.mods); Mf = np.eye(k, dtype=np.int64); off = 0
        for i, (p, a, M) in enumerate(self.parts):
            blk = np.eye(a, dtype=np.int64)
            if x[i]: blk = matpow(M, x[i], p)
            if x[self.r + i]: blk = (blk @ matpow(self.Fm[i], x[self.r + i], p)) % p
            Mf[off:off + a, off:off + a] = blk; off += a
        return Mf
    def elements(self):
        return [tuple(t) for t in itertools.product(*[range(n) for n in self.n], *[range(k) for k in self.k])]
    def neg_elem(self):
        return tuple([(n // 2) if p != 2 else 0 for n, p in zip(self.n, self.p)] + [0] * self.r)
    def normalizer_cosets(self, Hset, hgens):
        """coset representatives (identity first) of N(H)/H in the ambient group."""
        els = self.elements()
        N = []
        for g in els:
            gi = self.inv(g)
            if all(self.mul(self.mul(g, h), gi) in Hset for h in hgens): N.append(g)
        seen, reps = set(), []
        for g in sorted(N, key=lambda t: (t != self.ident, t)):
            if g in seen: continue
            reps.append(g)
            for h in Hset: seen.add(self.mul(g, h))
        return reps, len(N)


class GStruct:
    def __init__(self, v, e=None, hgens=None, mult=None):
        """e: exponent vector of a cyclic H inside prod C_{p^a - 1}  OR  hgens: list of ambient elements (tuples) generating H;
        mult: list of ambient elements = representatives of N(H)/H (default: Aut'/H for the cyclic case)."""
        self.v = v
        G, parts = build_group(v); self.G, self.parts = G, parts
        amb = Ambient(G, parts); self.amb = amb
        ns = [p ** a - 1 for p, a, M in parts]; self.ns = ns
        if hgens is None:
            hgens = [tuple(list(e) + [0] * len(parts))]
        self.hgens = hgens
        self.e = list(e) if e is not None else None
        Hset = amb.closure(hgens)
        self.Hset = Hset; self.Hord = len(Hset)
        mats = [amb.matrix(h) for h in hgens]
        perms = [G.perm_from_matrix(M) for M in mats]
        self.Mf = mats[0]
        self.neg_perm = G.neg(); self.neg = self.neg_perm
        self.stb = L.Structure(G, perms); self.stad = L.Structure(G, perms + [self.neg_perm])
        assert np.array_equal(self.stb.cls, self.stad.cls)
        self.C = self.stad.C
        self.csize = self.stad.csize
        self.creps = self.stad.reps
        self.A = GBlock(self.stad, v); self.B = GBlock(self.stb, v)
        self.has_m1 = bool(self.stb.O == self.stad.O)
        assert not self.has_m1, "-1 in H: not handled"
        coords = G.coords
        duals = [G.index(coords @ M) for M in mats]
        cl, k = L.closure_partition(v, duals + [self.neg_perm])
        zero = cl[0]
        labels = sorted(set(cl.tolist()) - {int(zero)})
        assert len(labels) == self.C, 'dual class count differs'
        wreps = [int(np.flatnonzero(cl == lab)[0]) for lab in labels]
        mods = G.mods.astype(np.float64)
        theta = np.zeros((self.C, self.C))
        for cp, w in enumerate(wreps):
            ph = (coords * coords[w]) / mods
            co = np.cos(2 * np.pi * ph.sum(1))
            for c in range(self.C):
                theta[cp, c] = co[self.stad.cls == c].sum()
        self.theta = theta
        if mult is None:
            if self.e is not None and len(hgens) == 1:
                order = self.Hord
                seen, reps = set(), []
                for ee in itertools.product(*[range(n) for n in ns]):
                    key = min(tuple((ee[i] + kk * self.e[i]) % ns[i] for i in range(len(ns))) for kk in range(order))
                    if key in seen: continue
                    seen.add(key); reps.append(tuple(list(ee) + [0] * len(parts)))
                reps.sort(key=lambda t: (any(t), t))
                mult = reps
            else:
                mult, self.Nsize = amb.normalizer_cosets(Hset, hgens)
        self.mult_elems = mult
        self.mult_reps = [G.perm_from_matrix(amb.matrix(x)) for x in mult]
        self.fix = [g for g in range(v) if self.stb.sizes[self.stb.orb[g]] == 1]

    def perm_of(self, block, t):
        return np.array([block.orb[int(t[r])] for r in block.reps], dtype=np.int64)

    def gamma_order(self):
        f = len(self.fix)
        return len(self.mult_reps) * 2 * f * f

    def feasible_patterns(self):
        out = []
        v = self.v; lim = math.isqrt(4 * v)
        for ra in range(1, lim + 1, 2):
            for rb in range(1, lim + 1, 2):
                rest = 4 * v - ra * ra - 2 * rb * rb
                if rest < 1: continue
                rd = math.isqrt(rest)
                if rd * rd == rest and rd % 2:
                    pa, pb, pd = self.A.pool(ra), self.B.pool(rb), self.B.pool(rd)
                    if min(pa, pb, pd) > 0: out.append(dict(rows=(ra, rb, rd), pools=(pa, pb, pd)))
        return out

    # ---- exact checks -----------------------------------------------------------------------------------------
    def masks_to_seqs(self, ma, mb, md):
        ya = np.array([1 if (ma >> i) & 1 else -1 for i in range(self.A.O)], dtype=np.int64)
        yb = np.array([1 if (mb >> i) & 1 else -1 for i in range(self.B.O)], dtype=np.int64)
        yd = np.array([1 if (md >> i) & 1 else -1 for i in range(self.B.O)], dtype=np.int64)
        return self.A.seq(ya), self.B.seq(yb), self.B.seq(yd)

    def exact_check(self, A, B, D, tgt=None):
        G, v = self.G, self.v
        if not np.array_equal(A, A[self.neg_perm]): return False
        tot = PG.group_paf(G, A) + 2 * PG.group_paf(G, B) + PG.group_paf(G, D)
        want = np.zeros(v, dtype=np.int64) if tgt is None else np.asarray(tgt)
        want = want.copy(); want[0] = 4 * v
        return bool(np.array_equal(tot, want))

    def matrix_check(self, A, B, D):
        Hm = PG.propus_matrix(self.G, A, B, D)
        return bool(np.array_equal(Hm, Hm.T) and np.array_equal(Hm @ Hm.T, 4 * self.v * np.eye(4 * self.v, dtype=Hm.dtype)))

    # ---- Gamma orbit of a solution (orbit-mask space) -----------------------------------------------------------
    def gamma_images(self, ma, mb, md):
        A, B = self.A, self.B
        def pm(blk, perm, m):
            return sum(1 << i for i in range(blk.O) if (m >> int(perm[i])) & 1)
        gens = []
        for t in self.mult_reps[1:]:
            pa, pb = self.perm_of(A, t), self.perm_of(B, t)
            gens.append(lambda x, pa=pa, pb=pb: (pm(A, pa, x[0]), pm(B, pb, x[1]), pm(B, pb, x[2])))
        pneg = self.perm_of(B, self.neg_perm)
        gens.append(lambda x: (x[0], pm(B, pneg, x[1]), x[2]))
        gens.append(lambda x: (x[0], x[1], pm(B, pneg, x[2])))
        for f in self.fix:
            if f == 0: continue
            tr = self.G.translate(f)
            ptr = np.array([B.orb[int(tr[r])] for r in B.reps])
            gens.append(lambda x, ptr=ptr: (x[0], pm(B, ptr, x[1]), x[2]))
            gens.append(lambda x, ptr=ptr: (x[0], x[1], pm(B, ptr, x[2])))
        seen = {(ma, mb, md)}; frontier = [(ma, mb, md)]
        while frontier:
            nxt = []
            for x in frontier:
                for g in gens:
                    y = g(x)
                    if y not in seen: seen.add(y); nxt.append(y)
            frontier = nxt
        return seen

    def is_canonical(self, ma, mb, md):
        B = self.B
        def pm(perm, m): return sum(1 << i for i in range(B.O) if (m >> int(perm[i])) & 1)
        for t in self.mult_reps[1:]:
            if pm(self.perm_of(B, t), mb) < mb: return False
        return pm(self.perm_of(B, self.neg_perm), md) >= md


def write_problem(st, rows, path, **kw):
    P.write_problem(st, rows, path, **kw)


def run(st, rows, workdir, threads=1, ns=2, noreduce=False, nopsd=False, shard=None, ckpt=None, resume=False, tgt=None, count=False, extra=()):
    tag = '_'.join(map(str, st.e)) if st.e is not None else 'h' + str(abs(hash(tuple(st.hgens))) % 10**8)
    pf = os.path.join(workdir, f'prob_{st.v}_{tag}_{"_".join(map(str, rows))}_{os.getpid()}.txt')
    write_problem(st, rows, pf, canonB=not noreduce, canonD=not noreduce, psd=(not nopsd,) * 3, seed=st.v, tgt=tgt)
    return P.run_p3(pf, binary=os.path.join(HERE, 'p3'), threads=threads, ns=ns, shard=shard, ckpt=ckpt, resume=resume, count=count, extra=extra)


# ---------------------------------------------------------------------------------------------------------------------
# "dmult" reduction: the multiplier group Mult = N(H)/H (joint symmetry of A, B, D) is spent on D instead of B.
#   D minimal under the group generated by Mult and the reversal of D alone  (p3 -dperms, D list shrinks by ~2|Mult|),
#   B <= reversal of B alone.   Lossless: both reductions use independent symmetries of the solution set.
class _OP(list):
    pass


class _Shim:
    def __init__(self, st, orbit_perms):
        self.__dict__.update(st.__dict__)
        self._st = st
        self.mult_reps = [None] + [_OP(map(int, p)) for p in orbit_perms]
    def perm_of(self, block, t):
        if isinstance(t, _OP): return np.array(t, dtype=np.int64)
        return self._st.perm_of(block, t)


def dmult_perms(st):
    """(B-perms for p3 problem file [reversal only], extra D perms): orbit permutations (y'_i = y_{pi[i]})."""
    neg = st.neg_perm
    rev = st.perm_of(st.B, neg)
    dp = []
    for t in st.mult_reps[1:]:
        dp.append(st.perm_of(st.B, t))
    for t in st.mult_reps:
        dp.append(st.perm_of(st.B, neg[t]))          # reversal of D composed with the multiplier t (t = identity: the reversal itself)
    return [rev], dp


def write_problem_dmult(st, rows, path, dperm_path, psd=(True, True, True), seed=1, tgt=None):
    bp, dp = dmult_perms(st)
    shim = _Shim(st, bp)
    P.write_problem(shim, rows, path, canonB=True, canonD=True, psd=psd, seed=seed, tgt=tgt)
    with open(dperm_path, 'w') as f:
        f.write(f'{len(dp)}\n')
        for p in dp: f.write(' '.join(map(str, p.tolist())) + '\n')


def is_canonical_dmult(st, ma, mb, md):
    B = st.B
    def pm(perm, m): return sum(1 << i for i in range(B.O) if (m >> int(perm[i])) & 1)
    bp, dp = dmult_perms(st)
    if pm(bp[0], mb) < mb: return False
    if pm(bp[0], md) < md: return False
    for p in dp:
        if pm(p, md) < md: return False
    return True


def run_dmult(st, rows, workdir, threads=1, ns=2, nopsd=False, shard=None, ckpt=None, resume=False, extra=(), tag=''):
    pf = os.path.join(workdir, f'probdm_{st.v}_{tag}_{"_".join(map(str, rows))}_{os.getpid()}.txt')
    dpf = pf + '.dperms'
    write_problem_dmult(st, rows, pf, dpf, psd=(not nopsd,) * 3, seed=st.v)
    return P.run_p3(pf, binary=os.path.join(HERE, 'p3'), threads=threads, ns=ns, shard=shard, ckpt=ckpt, resume=resume, extra=list(extra) + ['-dperms', dpf])
