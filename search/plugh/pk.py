"""K-symmetric plug-in systems (orbit-constant K-level functions) for the zero-diagonal Szollosi array over G = K x (conference core of order q).

Block X = s (x) I + t (x) Q + e (x) J,  s,t,e : K -> Z.   Per K-orbit of the symmetry group Gamma = <(h, phi(h))>:
   s, e constant on orbits;   t(h k) = phi(h) t(k)  (phi(h) = -1 means h acts on F_q by a NON-square multiplier, Q -> -Q).
A, B (symmetric blocks) additionally carry the generator (-1, eps)  (eps = +1 for q = 1 mod 4, -1 for q = 3 mod 4: t(-k) = eps t(k)).
Orbit data per block, per orbit O (weight w = |O|):
   non-defect: e = 0, s = +-1 (x0 = s), t = tau*sigma_k  (tau = +-1)              -> 2*2 choices
   defect    : t = 0, e = +-1, type p (s = 0, x0 = e) or m (s = -2e, x0 = -e)     -> 2*2 choices
   A/B at k = 0 (zero diagonal x0(0) = 0): non-defect s = 0, t = tau (only if the orbit is t-consistent), defect type z (s = -e).
   Orbits whose stabiliser contains an element with phi = -1 have t = 0 forced: forced defect.
K-level system (PROVED, search/NOTES.md 81):  (I) sum (s s~ + q t t~) = N2 delta,  (J) sum (q e e~ + s e~ + e s~ - t t~) = 0,
(Q) sum (t s~ + eps s t~) = 0;  equivalent triple used here: (U) sum PAF_u = N2 delta - 4q b,  u = s + q e;  (I);  (Q).
N2 = 4|G| - 2 (unbordered), 4|G| + 2 (bordered, b = 1), |G| = |K| q.
"""
import numpy as np, math, itertools, sys
from math import gcd


class Lat:
    """incremental integer lattice in Z^d (python ints, row echelon with entries reduced modulo later pivots)"""
    def __init__(self, d):
        self.d = d; self.B = {}

    def copy(self):
        L = Lat(self.d); L.B = {c: r[:] for c, r in self.B.items()}; return L

    def _red(self, r, c0):
        for c in sorted(k for k in self.B if k > c0):
            if r[c]:
                b = self.B[c]; q = r[c] // b[c]
                if q: r = [x - q * y for x, y in zip(r, b)]
        return r

    def add(self, r):
        r = [int(x) for x in r]
        for c in range(self.d):
            if r[c] == 0: continue
            if c not in self.B:
                if r[c] < 0: r = [-x for x in r]
                self.B[c] = self._red(r, c); return
            b = self.B[c]
            while r[c] != 0:
                qt = b[c] // r[c]
                b = [x - qt * y for x, y in zip(b, r)]
                b, r = r, b
            if b[c] < 0: b = [-x for x in b]
            self.B[c] = self._red(b, c)
            r = self._red(r, c)

    def add_many(self, rows):
        for r in rows: self.add(r)

    def rows(self):
        return [self.B[c] for c in sorted(self.B)]

    def contains(self, v):
        v = [int(x) for x in v]
        for c in sorted(self.B):
            if v[c]:
                b = self.B[c]
                if v[c] % b[c]: return False
                q = v[c] // b[c]; v = [x - q * y for x, y in zip(v, b)]
        return not any(v)

    def logvol(self):
        """natural log of the covolume of the lattice in its span: log sqrt(det Gram) (exact big-int determinant via numpy slogdet on a float Gram of the reduced basis)"""
        rows = self.rows()
        if not rows: return 0.0
        Bm = np.array(rows, dtype=np.float64)
        sgn, ld = np.linalg.slogdet(Bm @ Bm.T)
        return 0.5 * ld


class KGrp:
    """abelian group K = Z_m1 x ... x Z_mr (elements indexed in C order); automorphisms by component-wise unit multiplication"""
    def __init__(self, mods):
        self.mods = [int(m) for m in mods]
        self.n = int(np.prod(self.mods))
        self.el = np.array(list(itertools.product(*[range(m) for m in self.mods])), dtype=np.int64).reshape(self.n, len(self.mods))
        w = np.ones(len(self.mods), dtype=np.int64)
        for i in range(len(self.mods) - 2, -1, -1): w[i] = w[i + 1] * self.mods[i + 1]
        self.w = w
        m = np.array(self.mods, dtype=np.int64)
        self.add = np.zeros((self.n, self.n), dtype=np.int64)
        for a in range(self.n):
            self.add[a] = (((self.el[a][None, :] + self.el) % m) * w).sum(1)
        self.neg = (((-self.el) % m) * w).sum(1)
        self.m = m

    def mulperm(self, units):
        """permutation k -> (u_1 k_1, ..., u_r k_r); units: int (all components) or tuple; or an r x r integer matrix (elementary abelian K only)"""
        if isinstance(units, np.ndarray) and units.ndim == 1: units = list(units)
        if isinstance(units, (list, tuple)) and len(units) and isinstance(units[0], (list, tuple, np.ndarray)):
            M = np.array(units, dtype=np.int64)
            assert len(set(self.mods)) == 1
            img = (self.el @ M.T) % self.m
            return (img * self.w).sum(1)
        if np.isscalar(units): units = [int(units)] * len(self.mods)
        u = np.array([int(x) for x in units], dtype=np.int64)
        return (((self.el * u[None, :]) % self.m) * self.w).sum(1)


def orbits_signed(N, gens, mulgroup=None):
    """gens: list of (perm array (or unit for cyclic Z_N), phi).  Returns (orb id array, sigma array (+-1), forced array per orbit, list of position arrays).
    orbit of k: {g k}; sigma_k = product of phi along a path from the representative; inconsistent (stabiliser element with phi = -1) -> forced."""
    gens = [((np.asarray(u) if not np.isscalar(u) else (u * np.arange(N)) % N), ph) for (u, ph) in gens]
    orb = -np.ones(N, dtype=np.int64); sig = np.zeros(N, dtype=np.int64)
    pos = []; forced = []
    for k0 in range(N):
        if orb[k0] >= 0: continue
        oid = len(pos); orb[k0] = oid; sig[k0] = 1
        stack = [k0]; mem = [k0]; bad = False
        while stack:
            k = stack.pop()
            for (u, ph) in gens:
                k2 = int(u[k]); s2 = sig[k] * ph
                if orb[k2] < 0:
                    orb[k2] = oid; sig[k2] = s2; stack.append(k2); mem.append(k2)
                elif sig[k2] != s2:
                    bad = True
        pos.append(np.array(sorted(mem), dtype=np.int64)); forced.append(bad)
    return orb, sig, np.array(forced), pos


class Kind:
    """block ensemble of one kind ('AB' or 'CD') on the orbits of its group"""
    def __init__(self, shape, kind):
        self.shape = shape; self.kind = kind
        N, q, eps = shape.N, shape.q, shape.eps
        gens = list(shape.hgens)
        if kind == 'AB': gens = gens + [(shape.negperm, eps)]
        orb, sig, forced, pos = orbits_signed(N, gens)
        self.orb, self.sig, self.forced, self.pos = orb, sig, forced, pos
        self.nu = len(pos)
        self.w = np.array([len(p) for p in pos])
        self.zero_unit = int(orb[0]) if (kind == 'AB' and not shape.plain) else -1
        # forced defects (t = 0 forced).  A/B zero unit: if t-consistent it may be non-defect with s = t... (s=0,e=0,t=tau)
        self.free = np.array([not f for f in forced])
        # factor per unit for counting (non-defect and defect variants have the same count)
        self.fac = np.array([4.0] * self.nu)
        if kind == 'AB':
            self.fac[self.zero_unit] = 2.0
        # forced weight and DP over free units
        self.fw = int(self.w[forced].sum())
        fi = [i for i in range(self.nu) if not forced[i]]
        self.fi = fi
        maxw = N
        dp = [[0] * (maxw + 1) for _ in range(len(fi) + 1)]
        dp[len(fi)][0] = 1
        for ii in range(len(fi) - 1, -1, -1):
            wi = int(self.w[fi[ii]])
            for s_ in range(maxw + 1):
                dp[ii][s_] = dp[ii + 1][s_] + (dp[ii + 1][s_ - wi] if s_ >= wi else 0)
        self.dp = dp
        self.logfac = float(np.log10(self.fac).sum())
        # A/B zero unit: if the zero orbit is forced it is the z-defect; non-forced: nondefect (s=0,t=tau) or defect z.
        # q=1 mod 4 AB: zero unit may be forced? only if hgens contain a phi=-1 element fixing 0 (always fixes 0) -> yes: any phi=-1 gen forces t(0)=0.

    def weights(self):
        return [n for n in range(self.fw, self.shape.N + 1) if self.dp[0][n - self.fw] > 0]

    def lognblocks(self, n):
        c = self.dp[0][n - self.fw]
        return self.logfac + math.log10(c) if c > 0 else None

    def sample(self, n, M, rng, dp=None, fn=None, fd=None):
        """M blocks with total defect weight n -> arrays s,t,e (M,N) int16.  dp/fn/fd: optional weighted DP (factor fn[i] for non-defect, fd[i] for defect of unit i)"""
        dpx = self.dp if dp is None else dp
        N, q = self.shape.N, self.shape.q
        s = np.zeros((M, N), np.int16); t = np.zeros((M, N), np.int16); e = np.zeros((M, N), np.int16)
        rem = np.full(M, n - self.fw)
        defect = np.zeros((M, self.nu), bool)
        for i in range(self.nu):
            if self.forced[i]: defect[:, i] = True
        fi = self.fi
        for ii, i in enumerate(fi):
            wi = int(self.w[i])
            nxt = np.array([dpx[ii + 1][r] for r in range(N + 1)], dtype=float)
            skip = nxt[np.clip(rem, 0, N)] * (1.0 if fn is None else fn[i])
            take_w = np.where(rem >= wi, nxt[np.clip(rem - wi, 0, N)], 0.0) * (1.0 if fd is None else fd[i])
            p = take_w / np.maximum(take_w + skip, 1e-300)
            take = rng.random(M) < p
            defect[:, i] = take; rem = np.where(take, rem - wi, rem)
        for i in range(self.nu):
            P = self.pos[i]
            bit_e = rng.integers(0, 2, M) * 2 - 1       # sign of e (defect) or s (non-defect)
            bit_t = rng.integers(0, 2, M) * 2 - 1       # tau (non-defect) or type bit (defect: +1 -> p, -1 -> m)
            d = defect[:, i]
            zu = (i == self.zero_unit)
            for_s = np.zeros(M, np.int16); for_e = np.zeros(M, np.int16)
            # non-defect: s = bit_e (or 0 at the A/B zero unit), e = 0, t = tau*sigma
            # defect: e = bit_e; p: s = 0; m: s = -2 e; (zero unit: z: s = -e)
            if zu:
                for_s = np.where(d, -bit_e, 0).astype(np.int16)
            else:
                for_s = np.where(d, np.where(bit_t == 1, 0, -2 * bit_e), bit_e).astype(np.int16)
            for_e = np.where(d, bit_e, 0).astype(np.int16)
            tt = np.where(d, 0, bit_t).astype(np.int16)
            s[:, P] = for_s[:, None]; e[:, P] = for_e[:, None]
            t[:, P] = tt[:, None] * self.sig[P][None, :].astype(np.int16)
        return s, t, e

    def enumerate_all(self, nmax):
        """all blocks (as s,t,e arrays) with defect weight <= nmax; for exact counting of small shapes"""
        N, q = self.shape.N, self.shape.q
        units = list(range(self.nu))
        outs = []
        # per-unit variants: list of (defect?, bit_e, bit_t)
        variants = []
        for i in units:
            v = []
            zu = (i == self.zero_unit)
            for d in ((True,) if self.forced[i] else (False, True)):
                for be in (1, -1):
                    for bt in ((1,) if (zu and not d) or (zu and d) else (1, -1)):
                        v.append((d, be, bt))
            variants.append(v)
        for combo in itertools.product(*variants):
            wt = sum(int(self.w[i]) for i, (d, _, _) in zip(units, combo) if d)
            if wt > nmax: continue
            s = np.zeros(N, np.int16); t = np.zeros(N, np.int16); e = np.zeros(N, np.int16)
            for i, (d, be, bt) in zip(units, combo):
                P = self.pos[i]; zu = (i == self.zero_unit)
                if zu:
                    s[P] = -be if d else 0; e[P] = be if d else 0
                    t[P] = 0 if d else bt * self.sig[P]
                else:
                    if d:
                        e[P] = be; s[P] = 0 if bt == 1 else -2 * be; t[P] = 0
                    else:
                        s[P] = be; e[P] = 0; t[P] = bt * self.sig[P]
            outs.append((s, t, e, wt))
        return outs


class Shape:
    """Kmods: int (cyclic Z_N) or list of moduli;  hgens: list of (units, phi) (units: int or tuple, component-wise);
    mode: 'zd' (zero-diagonal conference array, unbordered), 'bord' (bordered), 'plain' (plain +-1 Szollosi Hadamard array: no zero position, N2 = 4|G|)"""
    def __init__(self, Kmods, q, hgens=(), bordered=False, mode=None):
        if mode is None: mode = 'bord' if bordered else 'zd'
        bordered = (mode == 'bord')
        self.mode = mode; self.plain = (mode == 'plain')
        if np.isscalar(Kmods): Kmods = [int(Kmods)]
        self.Kg = KGrp(Kmods); N = self.Kg.n
        self.Kmods = list(Kmods)
        self.N, self.q = N, q
        self.eps = 1 if q % 4 == 1 else -1
        self.hgens_raw = [(u, int(ph)) for (u, ph) in hgens]
        self.hgens = [(self.Kg.mulperm(u), int(ph)) for (u, ph) in hgens]
        self.negperm = self.Kg.neg
        self.bordered = bordered
        self.m = N * q
        self.N2 = 4 * self.m + 2 if bordered else (4 * self.m if self.plain else 4 * self.m - 2)
        self.kAB = Kind(self, 'AB'); self.kCD = Kind(self, 'CD')
        add = self.Kg.add          # add[a,b] = a+b
        self.shift = [np.ascontiguousarray(add[:, k]) for k in range(N)]   # shift[k][j] = j + k
        gensPM = [(u, 1) for (u, ph) in self.hgens] + [(self.negperm, 1)]
        orbPM, _, _, posPM = orbits_signed(N, gensPM)
        self.repsPM = [int(p[0]) for p in posPM]
        orbQ, _, forcedQ, posQ = orbits_signed(N, list(self.hgens) + [(self.negperm, self.eps)])
        self.repsQ = [int(p[0]) for p, f in zip(posQ, forcedQ) if not f]
        self.R = len(self.repsPM) * 2 + len(self.repsQ)
        self.nU = len(self.repsPM)
        T = np.zeros(self.R)
        for i, k in enumerate(self.repsPM):
            T[i] = (self.N2 if k == 0 else 0) - (4 * q if bordered else 0)
            T[self.nU + i] = self.N2 if k == 0 else 0
        self.T = T

    def eqvec(self, s, t, e, which='all'):
        q = self.q; eps = self.eps
        s = s.astype(np.int32); t = t.astype(np.int32); e = e.astype(np.int32)
        u = s + q * e
        out = []
        cols_u = []; cols_i = []; cols_q = []
        for k in self.repsPM:
            sh = self.shift[k]
            cols_u.append((u * u[:, sh]).sum(1))
        if which == 'u': return np.stack(cols_u, 1)
        if which == 'u0':
            # level-1 system + the k = 0 coefficient of (I) expressed through s and the defect indicator alone:
            # sum_X sum_k (s^2 - q [e != 0]) = 0   (from (I)(0): sum s^2 + q (#nondefect) = 4|G|/... = N2, N2 = 4|K| q in plain mode)
            c0 = (s * s - q * (e != 0)).sum(1)
            return np.stack(cols_u + [c0], 1)
        for k in self.repsPM:
            sh = self.shift[k]
            cols_i.append((s * s[:, sh]).sum(1) + q * (t * t[:, sh]).sum(1))
        for k in self.repsQ:
            sh = self.shift[k]
            cols_q.append((t * s[:, sh]).sum(1) + eps * (s * t[:, sh]).sum(1))
        return np.stack(cols_u + cols_i + cols_q, 1)

    def full_eqs(self, s, t, e):
        """all-k K-level quantities (single block, arrays (N,)): U, I, J, Q as length-N arrays"""
        N, q, eps = self.N, self.q, self.eps
        s = s.astype(np.int64); t = t.astype(np.int64); e = e.astype(np.int64)
        u = s + q * e
        paf = lambda x: np.array([(x * x[self.shift[k]]).sum() for k in range(N)])
        ccf = lambda x, y: np.array([(x * y[self.shift[k]]).sum() for k in range(N)])   # sum_j x(j) y(j+k)
        U = paf(u); I = paf(s) + q * paf(t)
        J = q * paf(e) + ccf(s, e) + ccf(e, s) - paf(t)
        Q = ccf(t, s) + eps * ccf(s, t)
        return U, I, J, Q


def hgroup(Kg, gens):
    """closure of the automorphism group (component-wise units, tuples) generated by gens"""
    r = len(Kg.mods)
    gens = [tuple(g) if not np.isscalar(g) else (int(g),) * r for g in gens]
    one = tuple([1] * r)
    G = {one}; frontier = [one]
    while frontier:
        x = frontier.pop()
        for g in gens:
            y = tuple((a * b) % m for a, b, m in zip(x, g, Kg.mods))
            if y not in G: G.add(y); frontier.append(y)
    return sorted(G)


def unit_tuples(Kg):
    return [u for u in itertools.product(*[[a for a in range(1, m) if gcd(a, m) == 1] for m in Kg.mods])]


def unit_subgroups(Kg, maxgens=3):
    """all subgroups of the component-wise unit group of K (cyclic and non-cyclic), via closure of subsets of <= maxgens generators"""
    units = unit_tuples(Kg)
    subs = {}
    for r in range(0, maxgens + 1):
        for gens in itertools.combinations(units, r):
            G = tuple(hgroup(Kg, list(gens)))
            if G not in subs: subs[G] = gens
    return subs


def signed_subgroups(Kmods, maxgens=3):
    """all pairs (subgroup H of Aut_diag(K), homomorphism phi: H -> {+-1}) given as generator lists [(units, phi)]; dedupe by the closure"""
    Kg = Kmods if isinstance(Kmods, KGrp) else KGrp([Kmods] if np.isscalar(Kmods) else Kmods)
    out = {}
    mods = Kg.mods
    for G, gens in unit_subgroups(Kg, maxgens).items():
        for signs in itertools.product((1, -1), repeat=len(gens)):
            one = tuple([1] * len(mods))
            pairs = {(one, 1)}; frontier = [(one, 1)]
            gp = list(zip(gens, signs))
            while frontier:
                (x, sx) = frontier.pop()
                for (g, sg_) in gp:
                    y = (tuple((a * b) % m for a, b, m in zip(x, g, mods)), sx * sg_)
                    if y not in pairs: pairs.add(y); frontier.append(y)
            if len(pairs) != len(G): continue          # not a homomorphism
            key = frozenset(pairs)
            if key not in out: out[key] = (G, [(g, sg_) for g, sg_ in gp])
    return list(out.values())
