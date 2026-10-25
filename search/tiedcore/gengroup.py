"""Tied and untied exact joins for the zero-diagonal Szollosi array (plain or bordered) over a finite abelian group
G = product of factors F(q) (finite field, additive group) and Z(k) (cyclic ring), with multiplier groups acting
factorwise (field: multiplication by g^s; ring: by a unit u).  Plain: sum PAF = 0 off 0.  Bordered: sum PAF = -4 off 0
(G of even order; row sums +-1, +-1, 0, 0)."""
import numpy as np, itertools, sys, time
from math import gcd

class GF:
    def __init__(self, q):
        from sympy import factorint
        f = factorint(q); assert len(f) == 1
        self.p, self.k = list(f.items())[0]; self.q = q
        p = self.p
        if self.k == 1:
            self.elts = [(a,) for a in range(p)]
            self.add = lambda u, v: ((u[0] + v[0]) % p,)
            self.mul = lambda u, v: ((u[0] * v[0]) % p,)
            self.tr = lambda u: u[0] % p
        elif self.k == 2:
            # find irreducible t^2 - c1 t - c0
            for c1 in range(p):
                for c0 in range(1, p):
                    if all((t * t - c1 * t - c0) % p for t in range(p)):
                        break
                else:
                    continue
                break
            self.c1, self.c0 = c1, c0
            self.elts = [(a, b) for a in range(p) for b in range(p)]
            def mul(u, v):
                a, b = u; c, d = v
                bd = b * d
                return ((a * c + bd * c0) % p, (a * d + b * c + bd * c1) % p)
            self.mul = mul
            self.add = lambda u, v: ((u[0] + v[0]) % p, (u[1] + v[1]) % p)
            self.tr = lambda u: (2 * u[0] + c1 * u[1]) % p
        else:
            # generic GF(p^k): coefficient tuples of 1, t, ..., t^(k-1), t^k = sum c_i t^i (monic irreducible)
            k = self.k
            def polymulmod(u, v, cs, k=k, p=p):
                prod = [0] * (2 * k - 1)
                for i, a in enumerate(u):
                    if a:
                        for j, b in enumerate(v):
                            prod[i + j] = (prod[i + j] + a * b) % p
                for d in range(2 * k - 2, k - 1, -1):
                    c = prod[d]
                    if c:
                        prod[d] = 0
                        for i in range(k):
                            prod[d - k + i] = (prod[d - k + i] + c * cs[i]) % p
                return tuple(prod[:k])
            for cs in itertools.product(range(p), repeat=k):
                if cs[0] == 0: continue
                # irreducible iff t has multiplicative order dividing p^k - 1 and the field has no zero divisors:
                # test that x^(p^k) = x for x = t and that no element of a proper subfield polynomial divides
                ok = True
                elts = list(itertools.product(range(p), repeat=k))
                one_ = tuple([1] + [0] * (k - 1))
                for u in elts[1:]:
                    for v in elts[1:]:
                        if polymulmod(u, v, cs) == tuple([0] * k):
                            ok = False; break
                    if not ok: break
                if ok: break
            self.cs = cs
            self.elts = list(itertools.product(range(p), repeat=k))
            self.mul = lambda u, v: polymulmod(u, v, cs)
            self.add = lambda u, v: tuple((a + b) % p for a, b in zip(u, v))
            def tr(u, k=k, p=p):
                x = u; s_ = list(u)
                for _ in range(k - 1):
                    y = tuple([1] + [0] * (k - 1))
                    for _ in range(p): y = polymulmod(y, x, cs)
                    x = y
                    s_ = [(a + b) % p for a, b in zip(s_, x)]
                assert all(c == 0 for c in s_[1:])
                return s_[0] % p
            self.tr = tr
        self.zero = self.elts[0]
        one = tuple([1] + [0] * (self.k - 1))
        self.one = one
        # primitive element
        for g in self.elts[1:]:
            x = one; o = 0
            seen = True
            for o in range(1, q):
                x = self.mul(x, g)
                if x == one:
                    break
            if o == q - 1:
                self.g = g; break
        self.ind = {}
        x = one
        for k in range(q - 1):
            self.ind[x] = k; x = self.mul(x, self.g)
        self.pos = {e: i for i, e in enumerate(self.elts)}
        self.neg = {e: self.mul(e, self.pow(self.g, (q - 1) // 2)) if e != self.zero else e for e in self.elts}
    def pow(self, x, n):
        r = self.one
        for _ in range(n): r = self.mul(r, x)
        return r


class Factor:
    def __init__(self, kind, q):
        self.kind, self.q = kind, q
        if kind == 'F':
            F = GF(q); self.F = F
            self.size = q
            self.shape = (F.p,) * F.k
            self.order_mult = q - 1
            self.ind = np.array([F.ind.get(x, -1) for x in F.elts])
            self.inv = np.zeros(q - 1, dtype=np.int64)
            for x, k in F.ind.items(): self.inv[k] = F.pos[x]
            tr = np.array([[F.tr(F.mul(a, x)) for x in F.elts] for a in F.elts])
            self.chi = np.exp(2j * np.pi * tr / F.p)            # chi[a, x]
            self.neg = (q - 1) // 2 if q % 2 else 0               # exponent of -1 (char 2: -1 = 1)
        else:
            self.size = q; self.shape = (q,)
            self.units = [u for u in range(1, q) if gcd(u, q) == 1] if q > 1 else [0]
            self.order_mult = len(self.units)
            a = np.arange(q)
            self.chi = np.exp(2j * np.pi * np.outer(a, a) / q)
            self.neg = q - 1 if q > 2 else 1
    def act(self, m):
        """permutation of factor positions for multiplier m (exponent for F, unit for Z)."""
        if self.kind == 'F':
            P = np.arange(self.size)
            nz = self.ind >= 0
            P[nz] = self.inv[(self.ind[nz] + m) % (self.q - 1)]
            return P
        return (np.arange(self.size) * m) % self.size
    def mcompose(self, a, b):
        return (a + b) % (self.q - 1) if self.kind == 'F' else (a * b) % self.q
    def mident(self):
        return 0 if self.kind == 'F' else 1
    def mall(self):
        return list(range(self.q - 1)) if self.kind == 'F' else self.units
    def minv(self, a):
        if self.kind == 'F': return (-a) % (self.q - 1)
        return pow(a, -1, self.q)

class Group:
    def __init__(self, spec):
        """spec like 'F25xF29' or 'Z2xF7xF17'."""
        self.spec = spec
        self.f = []
        for part in spec.split('x'):
            self.f.append(Factor(part[0], int(part[1:])))
        self.sizes = [f.size for f in self.f]
        self.m = int(np.prod(self.sizes))
        # W[x, a] = prod chi_i[a_i, x_i]
        W = np.ones((1, 1), dtype=complex)
        for f in self.f:
            W = np.einsum('xa,yb->xyab', W, f.chi.T).reshape(W.shape[0] * f.size, W.shape[1] * f.size)
        self.W = W
        self.neg = self.perm(tuple(f.neg if f.kind == 'F' else (f.q - 1 if f.q > 1 else 0) for f in self.f), negate=True)
    def perm(self, mult, negate=False):
        """element permutation for a factorwise multiplier tuple."""
        P = np.zeros(1, dtype=np.int64)
        for f, mm in zip(self.f, mult):
            if negate and f.kind == 'Z':
                p = (-np.arange(f.size)) % f.size
            elif negate and f.kind == 'F':
                p = f.act(f.neg) if f.q % 2 else np.arange(f.size)
            else:
                p = f.act(mm)
            P = (P[:, None] * f.size + p[None, :]).reshape(-1)
        return P
    def mult_elems(self):
        return list(itertools.product(*[f.mall() for f in self.f]))
    def mcomp(self, a, b):
        return tuple(f.mcompose(x, y) for f, x, y in zip(self.f, a, b))
    def mident(self):
        return tuple(f.mident() for f in self.f)
    def closure(self, gens):
        S = {self.mident()}; frontier = list(S)
        while frontier:
            nxt = []
            for x in frontier:
                for g in gens:
                    y = self.mcomp(x, g)
                    if y not in S: S.add(y); nxt.append(y)
            frontier = nxt
        return frozenset(S)
    def subgroups(self, max2=True):
        cyc = {}
        for e in self.mult_elems():
            H = self.closure([e]); cyc.setdefault(H, [e])
        subs = dict(cyc)
        if max2:
            items = list(cyc.items())
            for (H1, g1), (H2, g2) in itertools.combinations(items, 2):
                if H1 <= H2 or H2 <= H1: continue
                H = self.closure(g1 + g2)
                subs.setdefault(H, g1 + g2)
        return subs
    def orbits(self, perms, zero_forced=False):
        m = self.m
        parent = np.arange(m)
        def find(x):
            r = x
            while parent[r] != r: r = parent[r]
            while parent[x] != r:
                nx = parent[x]; parent[x] = r; x = nx
            return r
        for P in perms:
            for x in range(m):
                a, b = find(x), find(int(P[x]))
                if a != b: parent[a] = b
        roots = {}; lab = np.zeros(m, dtype=np.int64)
        for x in range(m):
            if zero_forced and x == 0: lab[x] = -1; continue
            lab[x] = roots.setdefault(find(x), len(roots))
        O = np.zeros((len(roots), m), dtype=np.int64)
        for x in range(m):
            if lab[x] >= 0: O[lab[x], x] = 1
        return O, lab

class Engine:
    def __init__(self, G, base_perms):
        self.G = G
        Ob, _ = G.orbits(base_perms)
        self.dual_reps = [int(np.nonzero(Ob[r])[0][0]) for r in range(len(Ob))]
        Ok, labk = G.orbits(list(base_perms) + [G.neg])
        self.Ok, self.labk = Ok, labk
        self.key_reps = [int(np.nonzero(Ok[r])[0][0]) for r in range(len(Ok))]
        self.M = (Ob @ G.W.real[self.key_reps, :].T) / G.m
        self.Wd = G.W[:, self.dual_reps]
    def paf(self, X, O, chunk=1 << 15):
        A = O @ self.Wd
        out = np.zeros((len(X), len(self.key_reps)), dtype=np.int16)
        for s in range(0, len(X), chunk):
            F = X[s:s + chunk].astype(float) @ A
            P = (F.real ** 2 + F.imag ** 2) @ self.M
            R = np.rint(P)
            assert np.abs(P - R).max() < 1e-4, np.abs(P - R).max()
            out[s:s + chunk] = R
        return out
    def key_perm(self, P):
        return np.array([self.labk[int(P[x])] for x in self.key_reps])

def all_signs(b):
    idx = np.arange(2 ** b, dtype=np.int64)
    return (1 - 2 * ((idx[:, None] >> np.arange(b)) & 1)).astype(np.int8)

def signs_with_sum(O, X, rowsum):
    """filter sign vectors X (over orbits O) to those with total sum rowsum (None = no filter)."""
    if rowsum is None: return X
    sz = O.sum(1)
    return X[(X.astype(np.int64) @ sz) == rowsum]

def keyvec(n, seed=12345):
    return np.random.default_rng(seed).integers(1, 2 ** 40, size=n, dtype=np.int64)

def lkeys(P, R, chunk=1 << 18):
    out = np.zeros(len(P), dtype=np.int64)
    for s in range(0, len(P), chunk):
        out[s:s + chunk] = P[s:s + chunk].astype(np.int64) @ R
    return out

def match(k1, k2):
    o2 = np.argsort(k2, kind='stable'); k2s = k2[o2]
    lo = np.searchsorted(k2s, k1, 'left'); hi = np.searchsorted(k2s, k1, 'right')
    res = []
    for i in np.nonzero(hi > lo)[0]:
        for jj in range(lo[i], hi[i]):
            res.append((int(i), int(o2[jj])))
    return res

def full_check(G, seqs, bordered):
    tot = 0
    shape = tuple(s for f in G.f for s in f.shape)
    for f in seqs:
        a = np.asarray(f, dtype=float).reshape(shape)
        Fh = np.fft.fftn(a)
        tot = tot + np.real(np.fft.ifftn(np.abs(Fh) ** 2))
    tot = np.rint(tot).astype(np.int64).flatten()
    want = -4 if bordered else 0
    return int(tot[0]), bool(np.all(tot[1:] == want))

def tied_run(G, Hgens, bordered=False, coarseA=(), coarseB=(), maxC=24, maxAB=26, stop_at=1, verbose=False):
    m = G.m
    sigs = [G.perm(h) for h in Hgens]
    OC, _ = G.orbits(sigs)
    OA, _ = G.orbits(sigs + [G.neg] + [G.perm(c) for c in coarseA], zero_forced=True)
    OB, _ = G.orbits(sigs + [G.neg] + [G.perm(c) for c in coarseB], zero_forced=True)
    bA, bB, bC = len(OA), len(OB), len(OC)
    if bC > maxC or bA + bB > maxAB:
        return None, (bA, bB, bC)
    eng = Engine(G, sigs)
    rA = rC = None
    XA = all_signs(bA); XB = all_signs(bB); XC = all_signs(bC)
    if bordered:
        XA = XA[np.abs(XA.astype(np.int64) @ OA.sum(1)) == 1]
        XB = XB[np.abs(XB.astype(np.int64) @ OB.sum(1)) == 1]
        XC = XC[(XC.astype(np.int64) @ OC.sum(1)) == 0]
    if len(XA) == 0 or len(XB) == 0 or len(XC) == 0:
        return (0, None), (bA, bB, bC)
    PA = eng.paf(XA, OA); PB = eng.paf(XB, OB); PC = eng.paf(XC, OC)
    K = len(eng.key_reps); R = keyvec(K)
    kAB = (lkeys(PA, R)[:, None] + lkeys(PB, R)[None, :]).reshape(-1)
    z = eng.labk[0]
    target = np.full(K, -4 if bordered else 0, dtype=np.int64); target[z] = 4 * m - 2
    kT = int(target @ R)
    H = G.closure(Hgens)
    seen = set(); total = 0; first = None
    for n in G.mult_elems():
        cls = frozenset(G.mcomp(n, h) for h in H)
        if cls in seen: continue
        seen.add(cls)
        ninv = tuple(f.minv(x) for f, x in zip(G.f, n))
        kp = eng.key_perm(G.perm(ninv))
        Rd = np.zeros_like(R); np.add.at(Rd, kp, R)
        kS = lkeys(PC, R + Rd)
        for i, j in match(kAB, kT - kS):
            a, b = divmod(i, len(PB))
            if np.array_equal(PA[a].astype(np.int64) + PB[b] + PC[j] + PC[j][kp], target):
                total += 1
                if first is None:
                    A = XA[a] @ OA; B = XB[b] @ OB; Cc = XC[j] @ OC
                    P = G.perm(n); D = np.zeros(m, dtype=np.int64); D[P] = Cc
                    first = (n, full_check(G, [A, B, Cc, D], bordered), [A.tolist(), B.tolist(), Cc.tolist(), D.tolist()])
        if stop_at and total >= stop_at: break
    return (total, first), (bA, bB, bC)
