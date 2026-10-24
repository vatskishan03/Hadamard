"""K-level system of the REAL Goethals-Seidel-type Paley plug-ins (four FREE blocks) over G = K x F_r.

Block X = s I + t Q + e J  (s,t,e : K -> Z; Q the Paley matrix of F_r; s = x0 - e, t = (x1 - x2)/2, e = (x1 + x2)/2, x0 on K x {0}, x1 on K x S, x2 on K x N).
Target on G:  F(0,0) = 4|G|,  F(k,f) = phi(k) for (k,f) != 0  (phi depends on k only):
   mode 'gs'     : phi = 0                       (sum PAF = 4|G| delta,  H(4|G|))
   mode 'spence' : phi = -4                      (sum PAF = -4 off 0, |G| even, row sums (2,0,0,0), H(4(|G|+1)))
   mode 'bord'   : phi(k) = -Gamma(k), Gamma = 4 sum_{chi in T} chi  (T: list of characters of K, i.e. exponent vectors; Gamma depends on k mod L,
                   L = common kernel;  K' = L x F_r;  n = |G| + |T|, H(4n))
K-level equations (all k in K):   I(k) = sum_X (s s~ + r t t~)(k)            = (4|G| - phi(0)) [k = 0]
                                  J(k) = sum_X (r e e~ + s e~ + e s~ - t t~)(k) = phi(k)
                                  Q(k) = sum_X (t s~ + eps s t~)(k)           = 0           (eps = +1 if r = 1 mod 4 else -1)
   U := I + r J = sum_X PAF_u,  u = s + r e (compression of the block to K):  U(0) = 4|G| + (r-1) phi(0),  U(k) = r phi(k) (k != 0)
Defect identity (from J(0)):  (r+1) e_p + (r-3) e_m = 4|K| + phi(0)   [e_p: defect with s = 0, e_m: defect with s = -2e; weights = orbit sizes]
Row sums R_X = sum_k u(k): sum_X R_X^2 = sum_g F(g) = 4|G| + sum_{k != 0 or f != 0} phi(k).
"""
import numpy as np, itertools, math, sys, os
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from pk import KGrp, orbits_signed, Kind
from idx import GF, Abel, factor_prime_power


class GShape:
    def __init__(self, Kmods, r, hgens=(), mode='gs', T=None):
        if np.isscalar(Kmods): Kmods = [int(Kmods)]
        self.Kmods = [int(m) for m in Kmods]
        self.Kg = KGrp(self.Kmods); N = self.Kg.n
        assert factor_prime_power(r) is not None and r % 2 == 1
        self.N, self.q, self.r = N, r, r              # .q kept for pk.Kind
        self.eps = 1 if r % 4 == 1 else -1
        self.hgens_raw = [(u, int(ph)) for (u, ph) in hgens]
        self.hgens = [(self.Kg.mulperm(u), int(ph)) for (u, ph) in hgens]
        self.negperm = self.Kg.neg
        self.plain = True
        self.mode = mode; self.G = N * r
        # target phi on K
        if mode == 'gs': self.chars = []; phi = np.zeros(N, dtype=np.int64)
        else:
            if mode == 'spence': T = [tuple([0] * len(self.Kmods))]
            assert T, 'bordered mode needs the character set T'
            self.chars = [tuple(int(x) for x in c) for c in T]
            m = np.array(self.Kmods)
            ph = np.zeros(N, dtype=complex)
            for c in self.chars:
                ph += np.exp(2j * np.pi * ((self.Kg.el * np.array(c)[None, :]) / m[None, :]).sum(1))
            assert np.abs(ph.imag).max() < 1e-9, 'T is not closed under complex conjugation'
            Gam = 4 * ph.real
            assert np.abs(Gam - np.rint(Gam)).max() < 1e-9, 'T is not a union of Galois orbits (Gamma not integral)'
            phi = -np.rint(Gam).astype(np.int64)
        self.phi = phi; self.t_border = len(self.chars)
        assert (phi[self.negperm] == phi).all()
        for (u, p_) in self.hgens: assert (phi[u] == phi).all(), 'phi not invariant under H_K'
        self.T2 = 4 * self.G - int(phi[0])            # I(0) target = bound of every PSD channel
        self.kCD = Kind(self, 'CD')
        self.nu = self.kCD.nu
        add = self.Kg.add
        self.shift = [np.ascontiguousarray(add[:, k]) for k in range(N)]       # shift[k][j] = j + k
        self.shneg = [np.ascontiguousarray(add[:, int(self.negperm[k])]) for k in range(N)]   # j - k
        gensPM = [(u, 1) for (u, ph_) in self.hgens] + [(self.negperm, 1)]
        _, _, _, posPM = orbits_signed(N, gensPM)
        self.repsPM = [int(p[0]) for p in posPM]
        self.nU = len(self.repsPM)
        _, _, forcedQ, posQ = orbits_signed(N, list(self.hgens) + [(self.negperm, self.eps)])
        self.repsQ = [int(p[0]) for p, f in zip(posQ, forcedQ) if not f]
        Ut = np.zeros(self.nU, dtype=np.int64)
        for i, k in enumerate(self.repsPM):
            Ut[i] = (4 * self.G + (r - 1) * int(phi[0])) if k == 0 else r * int(phi[k])
        self.Utarget = Ut
        self.Rsq = 4 * self.G + int(phi.sum()) * r - int(phi[0]) * 1 * 1 + 0   # placeholder fixed below
        # sum_g F(g) = 4|G| + sum_{(k,f) != 0} phi(k) = 4|G| + r sum_k phi(k) - phi(0)
        self.Rsq = 4 * self.G + r * int(phi.sum()) - int(phi[0])

    # ------------------------------------------------------------------ equations
    def full_eqs(self, s, t, e):
        N, r, eps = self.N, self.r, self.eps
        s = np.asarray(s, dtype=np.int64); t = np.asarray(t, dtype=np.int64); e = np.asarray(e, dtype=np.int64)
        u = s + r * e
        paf = lambda x: np.array([(x * x[self.shift[k]]).sum() for k in range(N)])
        ccf = lambda x, y: np.array([(x * y[self.shift[k]]).sum() for k in range(N)])   # sum_j x(j) y(j+k)
        U = paf(u); I = paf(s) + r * paf(t)
        J = r * paf(e) + ccf(s, e) + ccf(e, s) - paf(t)
        Q = ccf(s, t) + eps * ccf(t, s)
        return U, I, J, Q

    def check_klevel(self, blocks):
        """exact check of (I), (J), (Q), (U) over ALL k for four blocks (s,t,e)"""
        tot = [0, 0, 0, 0]
        for (s, t, e) in blocks:
            for i, v in enumerate(self.full_eqs(s, t, e)): tot[i] = tot[i] + v
        U, I, J, Q = tot
        tI = np.zeros(self.N, dtype=np.int64); tI[0] = self.T2
        ok = bool((I == tI).all() and (J == self.phi).all() and (Q == 0).all() and (U == I + self.r * J).all())
        return ok

    # ------------------------------------------------------------------ G level
    def glevel_blocks(self, blocks):
        """K-level (s,t,e) -> +-1 arrays on G (index k*r + a, a in the GF(r) index = additive Z_p^e mixed radix)"""
        p, ee = factor_prime_power(self.r)
        F = GF(p, ee)
        sq = {int(F.mul[x, x]) for x in range(1, self.r)}
        cls = np.array([0 if a == 0 else (1 if a in sq else 2) for a in range(self.r)])
        out = []
        for (s, t, e) in blocks:
            s = np.asarray(s, dtype=np.int64); t = np.asarray(t, dtype=np.int64); e = np.asarray(e, dtype=np.int64)
            vals = np.stack([s + e, e + t, e - t], 1)            # (N, 3)
            out.append(vals[:, cls].reshape(-1))
        return out

    def gmods(self):
        p, ee = factor_prime_power(self.r)
        return list(self.Kmods) + [p] * ee

    def gtarget(self):
        """target F on G (index k*r + a): F(0) = 4|G|, F(k, a) = phi(k) otherwise"""
        F = np.repeat(self.phi, self.r).astype(np.int64)
        F[0] = 4 * self.G
        return F

    def feasible_totals(self):
        """total defect weights E allowed by the defect identity (with e_m integral, 0 <= e_m <= E)"""
        r, nK, ph0 = self.r, self.N, int(self.phi[0])
        ds = []
        for d in range(0, 4 * nK + 1):
            num = (r + 1) * d - 4 * nK - ph0
            if num % 4: continue
            em = num // 4
            if 0 <= em <= d: ds.append(d)
        return ds


def exact_paf_G(mods, seqs):
    """exact integer sum of PAFs on the abelian group Z_m1 x ... (mixed radix index), no floating point"""
    mods = list(mods); n = int(np.prod(mods))
    coords = np.array(np.unravel_index(np.arange(n), mods)).T
    m = np.array(mods); w = np.ones(len(mods), dtype=np.int64)
    for i in range(len(mods) - 2, -1, -1): w[i] = w[i + 1] * mods[i + 1]
    tot = np.zeros(n, dtype=np.int64)
    for g in range(n):
        sh = (((coords + coords[g][None, :]) % m) * w).sum(1)
        tot[g] = sum(int(np.dot(x.astype(np.int64), x[sh].astype(np.int64))) for x in seqs)
    return tot


def verify_g(sh, blocks):
    """re-verify a K-level solution at G level: four +-1 sequences, exact sum PAF == target"""
    seqs = sh.glevel_blocks(blocks)
    ok1 = all(set(np.unique(x).tolist()) <= {-1, 1} for x in seqs)
    paf = exact_paf_G(sh.gmods(), seqs)
    return bool(ok1 and (paf == sh.gtarget()).all()), seqs


def char_set_cyclic(Kmods, qorder_divs):
    """convenience: characters of K = Z_m trivial on qZ_m  (Q = Z_q), restricted to orders in `qorder_divs` (all divisors d of q): exponent vectors (m/q * a)"""
    assert len(Kmods) == 1
    m = Kmods[0]; out = []
    from math import gcd
    for q_, divs in qorder_divs:
        pass
    return out


def chars_of_quotient(Kmods, L_gens, keep=None):
    """all characters of K trivial on the subgroup generated by L_gens (list of K elements as tuples); keep(c) filters (e.g. nontrivial only)"""
    m = list(Kmods); out = []
    for c in itertools.product(*[range(x) for x in m]):
        if all(sum(ci * li / mi for ci, li, mi in zip(c, l, m)) % 1 < 1e-9 or abs(sum(ci * li / mi for ci, li, mi in zip(c, l, m)) % 1 - 1) < 1e-9 for l in L_gens):
            if keep is None or keep(c): out.append(tuple(c))
    return out


def enum_weight(sh, n, rok=None, cansign=True, cap=40_000_000, mw=-1, sok=None, eok=None):
    """all (s,e) configs (distinct) of defect weight exactly n over the orbits of H_K, exact (explicit list).  rok: set of allowed row sums R = sum u, or None;
    cansign: canonical block sign (first nonzero of (s,e) over the orbits is positive).  Returns (S, E) int8 (M, N) or None."""
    K = sh.kCD; nu, N = K.nu, sh.N
    w = [int(x) for x in K.w]; forced = [bool(x) for x in K.forced]
    free = [i for i in range(nu) if not forced[i]]
    base = sum(w[i] for i in range(nu) if forced[i])
    subsets = []
    def rec(j, rem, cur):
        if rem == 0: subsets.append(list(cur)); return
        if j == len(free) or rem < 0: return
        i = free[j]
        if w[i] <= rem:
            cur.append(i); rec(j + 1, rem - w[i], cur); cur.pop()
        rec(j + 1, rem, cur)
    if n - base < 0: return None
    rec(0, n - base, [])
    NOND = np.array([[1, 0], [-1, 0]], dtype=np.int8)
    DEF = np.array([[0, 1], [-2, 1], [0, -1], [2, -1]], dtype=np.int8)
    outS, outE = [], []; tot = 0
    for sub in subsets:
        dset = set(sub) | {i for i in range(nu) if forced[i]}
        opts = [DEF if i in dset else NOND for i in range(nu)]
        sizes = [len(o) for o in opts]
        M = int(np.prod(sizes, dtype=np.float64)); tot += M
        if tot > cap: raise RuntimeError('list too large')
        idx = np.indices(sizes).reshape(nu, -1).T
        S = np.zeros((M, N), np.int8); E = np.zeros((M, N), np.int8)
        for i in range(nu):
            sel = opts[i][idx[:, i]]; P = K.pos[i]
            S[:, P] = sel[:, 0:1]; E[:, P] = sel[:, 1:2]
        keep = np.ones(M, bool)
        if cansign:
            sg = np.zeros(M, np.int8)
            for i in range(nu):
                P0 = int(K.pos[i][0])
                v = np.where(S[:, P0] != 0, S[:, P0], E[:, P0]); sg = np.where(sg == 0, v, sg)
            keep &= sg > 0
        if mw >= 0:
            mtype = ((E != 0) & (S != 0)).sum(1)          # weight of m-type defects (s = -2e != 0)
            keep &= (mtype == mw)
        if sok is not None: keep &= np.isin(S.astype(np.int64).sum(1), np.array(sorted(sok), dtype=np.int64))
        if eok is not None: keep &= np.isin(E.astype(np.int64).sum(1), np.array(sorted(eok), dtype=np.int64))
        if rok is not None:
            R = (S.astype(np.int64) + sh.r * E.astype(np.int64)).sum(1)
            keep &= np.isin(R, np.array(sorted(rok), dtype=np.int64))
        outS.append(S[keep]); outE.append(E[keep])
    if not outS: return None
    S = np.concatenate(outS); E = np.concatenate(outE)
    return (S, E) if len(S) else None


def space_size(sh, n):
    """number of distinct (s,e) configs of defect weight exactly n over the orbits of H_K (no sign gauge): dp[0][n - fw] * (4^forced)"""
    from est import kind_dpu
    K = sh.kCD
    if not hasattr(sh, '_dpu'): sh._dpu = kind_dpu(K)
    dp, ff = sh._dpu
    if n - K.fw < 0 or n - K.fw >= len(dp[0]): return 0.0
    return dp[0][n - K.fw] * ff


def block_classes(sh, n, m=-1):
    """EXACT set of (sigma, eps) = (sum_k s, sum_k e) over all (s,e) configs of defect weight n (and m-type weight m if m >= 0), over the orbits of H_K"""
    K = sh.kCD; w = [int(x) for x in K.w]; forced = [bool(x) for x in K.forced]
    states = {(0, 0, 0, 0)}
    for i in range(K.nu):
        wi = w[i]; new = set()
        for (nu_, mu_, sg, ep) in states:
            if not forced[i]:
                for s in (1, -1): new.add((nu_, mu_, sg + s * wi, ep))
            if nu_ + wi <= n:
                for e in (1, -1):
                    if m < 0 or True: new.add((nu_ + wi, mu_, sg, ep + e * wi))        # p-type: s = 0
                    if m < 0 or mu_ + wi <= m: new.add((nu_ + wi, mu_ + wi if m >= 0 else 0, sg - 2 * e * wi, ep + e * wi))   # m-type: s = -2e
        states = new
    return {(sg, ep) for (nu_, mu_, sg, ep) in states if nu_ == n and (m < 0 or mu_ == m)}


def row_patterns_v1(sh, al, ms=(-1, -1, -1, -1), maxpairs=6e6, normalize=True):
    """feasible row-sum patterns ((sigma_i, eps_i) for the four blocks) of an allocation: exist tau_i (|tau| <= N - n_i, parity of N - n_i) with
       sum sigma^2 + r sum tau^2 = T2 (I(0)),  [eps_Q = +1: sum sigma tau = 0],  sum (sigma + r eps)^2 = Rsq (sum of U).  Sound (tau only restricted to parity/range).
    Returns a sorted list of tuples ((sg0,ep0),...,(sg3,ep3)) normalised by the block-negation gauge, or None if the search is too large."""
    r = sh.r; T2 = sh.T2; Rsq = sh.Rsq; eps = sh.eps
    arrs = []
    for i in range(4):
        cl = sorted(block_classes(sh, al[i], ms[i]))
        Wnd = sh.N - al[i]
        taus = np.arange(-Wnd, Wnd + 1, 2)
        if not cl: return []
        sg = np.array([c[0] for c in cl]); ep = np.array([c[1] for c in cl])
        SG = np.repeat(sg, len(taus)); EP = np.repeat(ep, len(taus)); TA = np.tile(taus, len(cl))
        arrs.append((SG.astype(np.int64), EP.astype(np.int64), TA.astype(np.int64)))
    def comb(a, b):
        n1, n2 = len(a[0]), len(b[0])
        if n1 * n2 > maxpairs: return None
        i1 = np.repeat(np.arange(n1), n2); i2 = np.tile(np.arange(n2), n1)
        return i1, i2
    c01 = comb(arrs[0], arrs[1]); c23 = comb(arrs[2], arrs[3])
    if c01 is None or c23 is None: return None
    def keys(a, b, c):
        i1, i2 = c
        sg = (a[0][i1], b[0][i2]); ep = (a[1][i1], b[1][i2]); ta = (a[2][i1], b[2][i2])
        A = sg[0] ** 2 + sg[1] ** 2 + r * (ta[0] ** 2 + ta[1] ** 2)
        B = sg[0] * ta[0] + sg[1] * ta[1] if eps == 1 else np.zeros(len(i1), np.int64)
        C = (sg[0] + r * ep[0]) ** 2 + (sg[1] + r * ep[1]) ** 2
        return A, B, C
    A1, B1, C1 = keys(arrs[0], arrs[1], c01); A2, B2, C2 = keys(arrs[2], arrs[3], c23)
    BIG = np.int64(1 << 22)
    k1 = (A1 * BIG + (B1 + (1 << 20))) * BIG + C1
    k2 = ((T2 - A2) * BIG + (-B2 + (1 << 20))) * BIG + (Rsq - C2)
    o2 = np.argsort(k2, kind='stable'); k2s = k2[o2]
    lo = np.searchsorted(k2s, k1, 'left'); hi = np.searchsorted(k2s, k1, 'right')
    sel = np.flatnonzero(hi > lo)
    out = set()
    for i in sel:
        for j in o2[lo[i]:hi[i]]:
            a0 = c01[0][i]; a1 = c01[1][i]; a2 = c23[0][j]; a3 = c23[1][j]
            tup = tuple((int(arrs[b][0][a]), int(arrs[b][1][a])) for b, a in enumerate((a0, a1, a2, a3)))
            if normalize:
                tup = tuple(((-s_, -e_) if (s_ < 0 or (s_ == 0 and e_ < 0)) else (s_, e_)) for (s_, e_) in tup)
            out.add(tup)
    # symmetry among blocks with equal (n, m): sort their pairs
    res = set()
    for tup in out:
        t = list(tup)
        groups = {}
        for b in range(4): groups.setdefault((al[b], ms[b]), []).append(b)
        for g, idx in groups.items():
            vals = sorted(t[b] for b in idx)
            for b, v in zip(idx, vals): t[b] = v
        res.add(tuple(t))
    return sorted(res)


def row_patterns2(sh, al, ms=(-1, -1, -1, -1), budget=3e6, normalize=True):
    """Same contract as row_patterns, but organised by the row sums R_i = sigma_i + r eps_i (sum R_i^2 = Rsq is tiny for spence / bordered): enumerate the R tuples first
    (meet in the middle), then the (sigma, eps) choices with these R_i, then the tau tuples (sum tau^2 = (T2 - sum sigma^2)/r, eps = +1: sum sigma tau = 0, tau parity/range).
    Returns the sorted pattern list, or None when the work budget is exceeded."""
    r = sh.r; T2 = sh.T2; Rsq = sh.Rsq; eps = sh.eps
    byR = []
    for i in range(4):
        d = {}
        for (sg, ep) in block_classes(sh, al[i], ms[i]): d.setdefault(sg + r * ep, []).append((sg, ep))
        if not d: return []
        byR.append(d)
    # R tuples by MITM on pairs
    left = {}
    for R0 in byR[0]:
        for R1 in byR[1]: left.setdefault(R0 * R0 + R1 * R1, []).append((R0, R1))
    tuples = []
    for R2 in byR[2]:
        for R3 in byR[3]:
            need = Rsq - R2 * R2 - R3 * R3
            for (R0, R1) in left.get(need, ()): tuples.append((R0, R1, R2, R3))
    work = 0; out = set()
    Wnd = [sh.N - al[i] for i in range(4)]
    for Rt in tuples:
        choices = [byR[i][Rt[i]] for i in range(4)]
        ncomb = 1
        for c in choices: ncomb *= len(c)
        work += ncomb
        if work > budget: return None
        for comb in itertools.product(*choices):
            sg = [c[0] for c in comb]
            q2 = T2 - sum(x * x for x in sg)
            if q2 < 0 or q2 % r: continue
            q = q2 // r
            # tau tuples: tau_i in {-W..W step 2} (parity of W), sum tau^2 = q, eps = +1: sum sigma_i tau_i = 0
            ok = False
            A = {}
            for t0 in range(-Wnd[0], Wnd[0] + 1, 2):
                for t1 in range(-Wnd[1], Wnd[1] + 1, 2):
                    s2 = t0 * t0 + t1 * t1
                    if s2 > q: continue
                    A.setdefault((s2, (sg[0] * t0 + sg[1] * t1) if eps == 1 else 0), True)
            for t2 in range(-Wnd[2], Wnd[2] + 1, 2):
                for t3 in range(-Wnd[3], Wnd[3] + 1, 2):
                    s2 = t2 * t2 + t3 * t3
                    if s2 > q: continue
                    if (q - s2, (-(sg[2] * t2 + sg[3] * t3)) if eps == 1 else 0) in A: ok = True; break
                if ok: break
            work += 400
            if ok: out.add(tuple(comb))
    res = set()
    for tup in out:
        t = list(tup)
        if normalize: t = [((-s_, -e_) if (s_ < 0 or (s_ == 0 and e_ < 0)) else (s_, e_)) for (s_, e_) in t]
        groups = {}
        for b in range(4): groups.setdefault((al[b], ms[b]), []).append(b)
        for g, idx in groups.items():
            vals = sorted(t[b] for b in idx)
            for b, v in zip(idx, vals): t[b] = v
        res.add(tuple(t))
    return sorted(res)


def row_patterns(sh, al, ms=(-1, -1, -1, -1), budget=3e7, normalize=True):
    """feasible exact row-sum patterns ((sigma_i, eps_i) x 4) of an allocation: row_patterns2 (organised by R tuples; fast), falling back to the older MITM row_patterns_v1; None if both exceed their budgets"""
    P = row_patterns2(sh, al, ms, budget=budget, normalize=normalize)
    if P is None: P = row_patterns_v1(sh, al, ms, normalize=normalize)
    return P


def chars_cyclic(m, D):
    """characters of K = Z_m of order d for each d in D (d | m): exponent vectors [(m/d) j], gcd(j, d) = 1.  With Q = Z_q = K / L the bgs.py structures take D a set of divisors of q (L = qZ_m)."""
    from math import gcd
    out = []
    for d in sorted(D):
        assert m % d == 0
        for j in range(1 if d > 1 else 0, d if d > 1 else 1):
            if d == 1 or gcd(j, d) == 1: out.append(((m // d) * j % m,))
    return out


# ---------------------------------------------------------------------------------------------- search/NOTES.md section 92: tau-exact row patterns
def block_classes3(sh, n, m=-1):
    """EXACT set of (sigma, eps, tau) over all (s,e,t) configs of defect weight n (and m-type weight m if m >= 0): tau = sum over NON-defect orbits of t * (sum of sig over the orbit)"""
    K = sh.kCD; w = [int(x) for x in K.w]; forced = [bool(x) for x in K.forced]
    Tw = [int(sum(int(K.sig[p]) for p in K.pos[o])) for o in range(K.nu)]
    states = {(0, 0, 0, 0, 0)}
    for i in range(K.nu):
        wi = w[i]; Ti = Tw[i]; new = set()
        for (nu_, mu_, sg, ep, ta) in states:
            if not forced[i]:
                for s in (1, -1):
                    for t in (1, -1): new.add((nu_, mu_, sg + s * wi, ep, ta + t * Ti))
            if nu_ + wi <= n:
                for e in (1, -1):
                    new.add((nu_ + wi, mu_, sg, ep + e * wi, ta))
                    if m < 0 or mu_ + wi <= m: new.add((nu_ + wi, mu_ + wi if m >= 0 else 0, sg - 2 * e * wi, ep + e * wi, ta))
        states = new
    return {(sg, ep, ta) for (nu_, mu_, sg, ep, ta) in states if nu_ == n and (m < 0 or mu_ == m)}


def row_patterns3(sh, al, ms=(-1, -1, -1, -1), maxpairs=6e6, normalize=True):
    """feasible exact row-sum patterns ((sigma_i, eps_i) x 4) with the EXACT tau sets (block_classes3).  Same contract as row_patterns_v1; None if too large."""
    r = sh.r; T2 = sh.T2; Rsq = sh.Rsq; eps = sh.eps
    arrs = []
    for i in range(4):
        cl = sorted(block_classes3(sh, al[i], ms[i]))
        if not cl: return []
        arrs.append((np.array([c[0] for c in cl], dtype=np.int64), np.array([c[1] for c in cl], dtype=np.int64), np.array([c[2] for c in cl], dtype=np.int64)))
    def comb(a, b):
        n1, n2 = len(a[0]), len(b[0])
        if n1 * n2 > maxpairs: return None
        return np.repeat(np.arange(n1), n2), np.tile(np.arange(n2), n1)
    c01 = comb(arrs[0], arrs[1]); c23 = comb(arrs[2], arrs[3])
    if c01 is None or c23 is None: return None
    def keys(a, b, c):
        i1, i2 = c
        sg = (a[0][i1], b[0][i2]); ep = (a[1][i1], b[1][i2]); ta = (a[2][i1], b[2][i2])
        A = sg[0] ** 2 + sg[1] ** 2 + r * (ta[0] ** 2 + ta[1] ** 2)
        B = sg[0] * ta[0] + sg[1] * ta[1] if eps == 1 else np.zeros(len(i1), np.int64)
        C = (sg[0] + r * ep[0]) ** 2 + (sg[1] + r * ep[1]) ** 2
        return A, B, C
    A1, B1, C1 = keys(arrs[0], arrs[1], c01); A2, B2, C2 = keys(arrs[2], arrs[3], c23)
    BIG = np.int64(1 << 22)
    k1 = (A1 * BIG + (B1 + (1 << 20))) * BIG + C1
    k2 = ((T2 - A2) * BIG + (-B2 + (1 << 20))) * BIG + (Rsq - C2)
    o2 = np.argsort(k2, kind='stable'); k2s = k2[o2]
    lo = np.searchsorted(k2s, k1, 'left'); hi = np.searchsorted(k2s, k1, 'right')
    out = set()
    for i in np.flatnonzero(hi > lo):
        for j in o2[lo[i]:hi[i]]:
            a0 = c01[0][i]; a1 = c01[1][i]; a2 = c23[0][j]; a3 = c23[1][j]
            tup = tuple((int(arrs[b][0][a]), int(arrs[b][1][a])) for b, a in enumerate((a0, a1, a2, a3)))
            if normalize: tup = tuple(((-s_, -e_) if (s_ < 0 or (s_ == 0 and e_ < 0)) else (s_, e_)) for (s_, e_) in tup)
            out.add(tup)
    res = set()
    for tup in out:
        t = list(tup); groups = {}
        for b in range(4): groups.setdefault((al[b], ms[b]), []).append(b)
        for g, idx in groups.items():
            vals = sorted(t[b] for b in idx)
            for b, v in zip(idx, vals): t[b] = v
        res.add(tuple(t))
    return sorted(res)


def row_patterns(sh, al, ms=(-1, -1, -1, -1), budget=3e7, normalize=True):
    """feasible exact row-sum patterns ((sigma_i, eps_i) x 4): the tau-exact version (row_patterns3) when it fits, else row_patterns2 / row_patterns_v1 (sound supersets)"""
    P = row_patterns3(sh, al, ms, normalize=normalize)
    if P is None: P = row_patterns2(sh, al, ms, budget=budget, normalize=normalize)
    if P is None: P = row_patterns_v1(sh, al, ms, normalize=normalize)
    return P


def enum_pattern(sh, n, m, sg, ep, cansign=True, cap=3_000_000, count_only=False):
    """EXACT list of all (s,e) configs with defect weight n, m-type weight m (m < 0: any), sum s = sg, sum e = ep (DP-pruned depth-first enumeration; no unfiltered blow-up).
    Returns (S, E) int8 arrays, or the count if count_only; raises RuntimeError if more than cap configs."""
    K = sh.kCD; nu, N = K.nu, sh.N
    w = [int(x) for x in K.w]; forced = [bool(x) for x in K.forced]
    opts = []        # per orbit: list of (dn, dm, ds, de, s, e)
    for i in range(nu):
        o = []
        if not forced[i]:
            for s in (1, -1): o.append((0, 0, s * w[i], 0, s, 0))
        for e in (1, -1):
            o.append((w[i], 0, 0, e * w[i], 0, e))
            o.append((w[i], w[i], -2 * e * w[i], e * w[i], -2 * e, e))
        opts.append(o)
    # backward counts: cnt[i][state] = number of completions using orbits i.. reaching the state (sums of dn, dm, ds, de)
    cnt = [None] * (nu + 1); cnt[nu] = {(0, 0, 0, 0): 1}
    for i in range(nu - 1, -1, -1):
        d = {}
        for (a, b, c, e_), v in cnt[i + 1].items():
            for (dn, dm, ds, de, _, _) in opts[i]:
                k = (a + dn, b + dm, c + ds, e_ + de)
                if k[0] > n: continue
                if m >= 0 and k[1] > m: continue
                d[k] = d.get(k, 0) + v
        cnt[i] = d
    target = (n, m if m >= 0 else None, sg, ep)
    def match(st):
        return st[0] == n and st[2] == sg and st[3] == ep and (m < 0 or st[1] == m)
    total = 0
    if m >= 0: total = cnt[0].get((n, m, sg, ep), 0)
    else: total = sum(v for k, v in cnt[0].items() if k[0] == n and k[2] == sg and k[3] == ep)
    if count_only: return total
    if total > cap: raise RuntimeError(f'pattern list too large ({total})')
    rows_s = []; rows_e = []
    # remaining feasibility for DFS: need cnt of suffix reaching (n - a, ...) ; use cnt[i+1] membership
    def dfs(i, a, b, c, e_, cs, ce):
        if i == nu:
            if a == n and c == sg and e_ == ep and (m < 0 or b == m): rows_s.append(list(cs)); rows_e.append(list(ce))
            return
        for (dn, dm, ds, de, s, e) in opts[i]:
            a2, b2, c2, e2 = a + dn, b + dm, c + ds, e_ + de
            if a2 > n or (m >= 0 and b2 > m): continue
            # suffix must be able to complete: need (n - a2, m - b2 or any, sg - c2, ep - e2) in cnt[i+1]
            need = (n - a2, (m - b2) if m >= 0 else None, sg - c2, ep - e2)
            if m >= 0:
                if need not in cnt[i + 1]: continue
            else:
                if not any(True for _ in [0]) : pass
                if not any(k[0] == need[0] and k[2] == need[2] and k[3] == need[3] for k in cnt[i + 1]): continue
            cs.append(s); ce.append(e); dfs(i + 1, a2, b2, c2, e2, cs, ce); cs.pop(); ce.pop()
    dfs(0, 0, 0, 0, 0, [], [])
    M = len(rows_s)
    S = np.zeros((M, N), np.int8); E = np.zeros((M, N), np.int8)
    if M:
        RS = np.array(rows_s, dtype=np.int8); RE = np.array(rows_e, dtype=np.int8)
        for i in range(nu):
            P = K.pos[i]; S[:, P] = RS[:, i:i + 1]; E[:, P] = RE[:, i:i + 1]
        if cansign:
            sgn = np.zeros(M, np.int8)
            for i in range(nu):
                P0 = int(K.pos[i][0]); v = np.where(S[:, P0] != 0, S[:, P0], E[:, P0]); sgn = np.where(sgn == 0, v, sgn)
            keep = sgn > 0; S = S[keep]; E = E[keep]
    return (S, E) if len(S) else None
