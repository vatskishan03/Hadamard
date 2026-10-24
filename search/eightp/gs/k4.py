"""k4.py -- K-level system for the REAL Goethals-Seidel-type Paley plug-in (four free blocks), modes a (unbordered), b (Spence), c (bordered with T).

G = K x F_r (r odd prime power), blocks X = s I + t Q + e J (Q Paley matrix of F_r), all four blocks free (kind 'CD' of pk.py).
Target (all modes):  sum_i PAF_i(g) = 4|G| delta_g - (1 - delta_g) Gamma(pi(g)),  pi: G -> Q = K/L,  Gamma = 4 sum_{chi in T} chi,
so on K-level  phi(k) = -Gamma(k mod L)  (phi(0) is the value on {0} x F_r^*).  Group ring: F = phi (x) J + (4|G| - phi(0)) delta, hence (PROVED below, checked by brute force)
   (I)  sum (s s~ + r t t~)                       = (4|G| - phi(0)) delta_K
   (J)  sum (r e e~ + s e~ + e s~ - t t~)         = phi(k)
   (Q)  sum (t s~ + eps s t~)                     = 0              eps = chi(-1) = +1 (r = 1 mod 4), -1 (r = 3 mod 4)
   (U) = (I) + r (J):  sum PAF_u = (4|G| - phi(0)) delta + r phi(k),  u = s + r e = x0 + (r-1)/2 (x1 + x2)... (u = s + r e)
Defect identity (J(0) with s = 0 (type p) or s = -2e (type m) at defects, non-defect e = 0, t^2 = 1):
   (r + 1) e_p + (r - 3) e_m = 4|K| + phi(0)        (e_p, e_m = number of defect positions k in K over the four blocks; weights = orbit sizes)
Row sums rho_i = sum_G X_i = sum_k u_i(k):  sum rho_i^2 = (4|G| - phi(0)) + r sum_k phi(k).
Modes:  a: phi = 0 (|G| = n);  b (Spence): phi = -4 everywhere (|G| = n - 1);  c: phi = -Gamma(k mod L), t = |T|, |G| = n - t.
"""
import numpy as np, math, itertools, sys, time
from math import gcd
from pk import *


def gamma_table(qmods, T):
    """Gamma on Q = prod Z_{q_i}, as dict {x tuple: int};  T = list of characters a (tuples); Gamma(x) = 4 sum_a cos(2 pi sum a_i x_i / q_i)"""
    out = {}
    for x in itertools.product(*[range(q) for q in qmods]):
        v = 0.0
        for a in T:
            v += math.cos(2 * math.pi * sum(ai * xi / qi for ai, xi, qi in zip(a, x, qmods)))
        v *= 4
        assert abs(v - round(v)) < 1e-8, ('Gamma not integer: T is not a union of Galois orbits', T)
        out[x] = int(round(v))
    return out


def is_galois_closed(qmods, T):
    if not qmods: return True
    E = 1
    for q in qmods: E = E * q // gcd(E, q)
    S = set(tuple(a) for a in T)
    for a in list(S):
        for u in range(1, E):
            if gcd(u, E) == 1:
                if tuple((u * ai) % qi for ai, qi in zip(a, qmods)) not in S: return False
    return True


def make_phi(Kmods, qspec):
    """qspec = (qmods, T): quotient K -> Q = prod Z_{q_i} by k -> (k_i mod q_i) for the first len(qmods) coordinates of Kmods (q_i | m_i); T character list.
    qspec None = mode a (phi = 0).  Returns phi array over K (C order), Gamma(0) = 4t"""
    Kg = KGrp(Kmods)
    if qspec is None: return np.zeros(Kg.n, dtype=np.int64), 0
    qmods, T = qspec
    assert all(Kmods[i] % qmods[i] == 0 for i in range(len(qmods)))
    gam = gamma_table(qmods, T)
    phi = np.zeros(Kg.n, dtype=np.int64)
    for k in range(Kg.n):
        x = tuple(int(Kg.el[k][i]) % qmods[i] for i in range(len(qmods)))
        phi[k] = -gam[x]
    return phi, len(T)


class Shape4(Shape):
    """pk.Shape (kind CD only used) with the real-GS target.  Attributes: phi, t (border size), r, N, m = |G|, n = |G| + t, D = identity target 4|K| + phi(0)"""
    def __init__(self, Kmods, r, hgens=(), qspec=None):
        super().__init__(Kmods, r, hgens, mode='plain')
        phi, t = make_phi(Kmods, qspec)
        self.phi = phi; self.tborder = t; self.qspec = qspec
        N = self.N; m = N * r
        # H-invariance of the target
        for (u, ph) in self.hgens:
            assert (phi[u] == phi).all(), 'target phi not invariant under the multiplier'
        self.N2 = 4 * m - int(phi[0])
        self.nord = m + t
        self.Did = 4 * N + int(phi[0])
        T = np.zeros(self.R)
        for i, k in enumerate(self.repsPM):
            T[i] = (self.N2 if k == 0 else 0) + r * int(phi[k])
            T[self.nU + i] = self.N2 if k == 0 else 0
        self.T = T
        self.rowsq = self.N2 + r * int(phi.sum())          # sum of squares of the four row sums

    def kind(self):
        return self.kCD


# ------------------------------------------------------------------ G-level brute force (independent of the s,t,e algebra)
from idx import GF, Abel, factor_prime_power
from gbrute import build as _gbuild, orbits as _gorbits, block_list as _gblock, paf_batch as _gpaf


def count_quads(V, w, T, seed=12345):
    """number of ordered 4-tuples of rows (with multiplicity weights w) of V (n x R int) summing to T; hash meet in the middle on unordered pairs"""
    V = np.asarray(V, dtype=np.int64); n = len(V)
    rng = np.random.default_rng(seed)
    coef = rng.integers(1, 2 ** 61, V.shape[1], dtype=np.int64) | 1
    h = (V * coef).sum(1); hT = int((np.asarray(T, dtype=np.int64) * coef).sum())
    w = np.asarray(w, dtype=np.float64)
    iu, ju = np.triu_indices(n)
    ph = h[iu] + h[ju]
    pw = w[iu] * w[ju] * np.where(iu == ju, 1.0, 2.0)
    order = np.argsort(ph, kind='stable'); ph = ph[order]; pw = pw[order]
    uniq, start = np.unique(ph, return_index=True)
    cs = np.concatenate([[0.0], np.cumsum(pw)])
    end = np.concatenate([start[1:], [len(ph)]])
    mass = cs[end] - cs[start]
    need = hT - uniq
    pos = np.searchsorted(uniq, need)
    pos = np.clip(pos, 0, len(uniq) - 1)
    ok = uniq[pos] == need
    return float((mass[ok] * mass[pos[ok]]).sum())


def dedupe(V, extra=None):
    key = V if extra is None else np.concatenate([V, extra[:, None]], 1)
    u, first, inv, cnt = np.unique(key, axis=0, return_index=True, return_inverse=True, return_counts=True)
    return first, cnt


def g_target(Kmods, r, phi):
    Kg = KGrp(Kmods); N = Kg.n
    tg = np.zeros(N * r, dtype=np.int64)
    for k in range(N):
        tg[k * r:(k + 1) * r] = phi[k]
    tg[0] = 4 * N * r
    return tg


def count_G4(Kmods, r, hgens, phi, maxpairs=4e7, verbose=False):
    """exact ordered count of quadruples of +-1 functions on G = K x F_r invariant under H_G with sum PAF = target (G-level brute force)"""
    F, G, maps, neg = _gbuild(Kmods, r, hgens)
    n = G.n
    lab, nor = _gorbits(n, maps)
    if nor > 17: raise RuntimeError(f'G-level list too large: {nor} orbits')
    X = _gblock(G, lab, nor, None)
    P = _gpaf(G, X)
    reps = [g for g in range(n) if g <= int(neg[g])]
    P = P[:, reps]
    tg = g_target(Kmods, r, phi)[reps]
    first, cnt = dedupe(P)
    Pu = P[first]
    nu = len(Pu)
    if nu * (nu + 1) / 2 > maxpairs: raise RuntimeError(f'too many pairs {nu}')
    if verbose: print(f'  G-level: |G|={n}, orbits {nor}, list {len(X)}, unique PAF vectors {nu}', flush=True)
    return count_quads(Pu, cnt, tg), (nor, len(X), nu)


# ------------------------------------------------------------------ K-level exact count
def count_K4(sh, maxpairs=4e7, verbose=False):
    """exact ordered count of K-level solutions (four blocks from the CD ensemble: (s,t,e) with all the equations (U), (I), (Q) at orbit representatives)"""
    r = sh.q; N = sh.N
    hi = 4 * N * 4
    K = sh.kCD
    L = K.enumerate_all(N)
    s = np.stack([b[0] for b in L]); t = np.stack([b[1] for b in L]); e = np.stack([b[2] for b in L])
    V = np.concatenate([sh.eqvec(s[i:i + 100000], t[i:i + 100000], e[i:i + 100000], 'all') for i in range(0, len(L), 100000)]).astype(np.int64)
    # energy pruning on k = 0 U column: each block's PAF_u(0) = sum u^2 >= 0 and the total is T[0]
    keep = V[:, 0] <= sh.T[0]
    V = V[keep]
    first, cnt = dedupe(V)
    Vu = V[first]
    nu = len(Vu)
    if verbose: print(f'  K-level: list {len(L)}, unique vectors {nu}', flush=True)
    if nu * (nu + 1) / 2 > maxpairs: raise RuntimeError(f'too many pairs {nu}')
    return count_quads(Vu, cnt, sh.T.astype(np.int64)), (len(L), nu)


def cyclic_unit_subgroups(Kmods, maxorb=None, minorb=None, twist=True, maxn=None):
    """all (cyclic subgroup <u> of the component-wise unit group of K, sign homomorphism) as generator lists [(u tuple, phi)], deduplicated by the closed signed set"""
    Kg = KGrp(Kmods); mods = Kg.mods
    units = unit_tuples(Kg)
    seen = {}
    one = tuple([1] * len(mods))
    for u in units:
        # closure of <u> with sign
        for sg in ((1, -1) if twist else (1,)):
            pairs = [(one, 1)]; x = one; s = 1
            ok = True
            while True:
                x = tuple((a * b) % m for a, b, m in zip(x, u, mods)); s = s * sg
                if x == one:
                    ok = (s == 1)
                    break
                pairs.append((x, s))
            if not ok: continue
            key = frozenset(pairs)
            if key not in seen: seen[key] = ([(u, sg)] if (u != one or sg != 1) else [], len(pairs))
    return list(seen.values())


def orbit_count(Kmods, hg):
    Kg = KGrp(Kmods)
    perms = [(Kg.mulperm(u), ph) for (u, ph) in hg]
    _, _, forced, pos = orbits_signed(Kg.n, perms)
    return len(pos), int(forced.sum())


def abelian_groups_pp(n):
    """all abelian groups of order n as sorted lists of prime-power moduli"""
    def factorint(n):
        f = {}; d = 2
        while d * d <= n:
            while n % d == 0: f[d] = f.get(d, 0) + 1; n //= d
            d += 1
        if n > 1: f[n] = f.get(n, 0) + 1
        return f
    f = factorint(n)
    def parts(a, mx=None):
        mx = a if mx is None else mx
        if a == 0: yield []; return
        for k in range(min(a, mx), 0, -1):
            for r in parts(a - k, k): yield [k] + r
    per = []
    for p, a in f.items():
        per.append([[p ** k for k in part] for part in parts(a)])
    out = []
    for combo in itertools.product(*per):
        out.append(sorted([m for c in combo for m in c]))
    return out if n > 1 else [[]]
