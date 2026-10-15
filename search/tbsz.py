#!/usr/bin/env python3
"""The character-bordered Szöllősi array: complex Hadamard matrices of order 2N from a group of order N - t
(complex-hadamard/TBORDERED.md).

G abelian of order m, phi: G -> Q onto an abelian group of order q with kernel K, T a set of t characters of Q that is
a union of Galois orbits, N = m + t, Gamma(a) = 4 sum_{chi in T} chi(a). Four +-1 sequences on G, A and B symmetric,
with

    sum_k PAF_k(g) = -Gamma(phi(g))   (g != 0),

and a border (Xb, Rh, Ch) of quaternary matrices of sizes 2t x 2t, 2t x 2q, 2q x 2t give K = [[Xb, R], [R', Kc]]: Kc is
Szöllősi's array over G, R[a, (i, g)] = Rh[a, (i, phi g)], R'[(i, g), b] = Ch[(i, phi g), b]. The border needs rows of
Rh and columns of Ch with Fourier support in T on each half, Rh Rh^* = Ch^* Ch = 2q I, and Xb = -Rh Mc^* Ch / 2q
quaternary, Mc Szöllősi's array over Q of the compressions (fibre sums over the cosets of K). t = 1, Q = 1 is the
bordered array (BORDERED.md); Q = Z_q with T all characters needs q | N; Q = Z_4 with T = {1, chi, chi-bar} needs
N = 3 (mod 4) and reaches primes N, where Z_N and Z_(N-1) have no usable multipliers.

Symmetry: H <= Aut(G) preserving K and T; C, D constant on H-orbits, A, B on <H, -1>-orbits; the condition holds per
<H, -1>-class. Search: per compression tuple on Q, every block is enumerated by lowe_enum with its fibre sums fixed and
its low-energy members kept, the lists are joined by kt4 with the per-class targets, and each hit is completed by the
border solver and checked exactly.

    python3 search/tbsz.py test
    python3 search/tbsz.py rank N --family z4cyc|v4n|v4t [--omin 12 --omax 60]
    python3 search/tbsz.py search N --family z4cyc --gen a,b --rows 2,0,2,2 [--tuple i] [--L 3e7] [--nsig 2e5]
    python3 search/tbsz.py certify --hits FILE [--k 0]            writes complex-hadamard/C<2N>.json
"""
import argparse, itertools, json, math, os, struct, subprocess, sys, time
from math import gcd, log10, pi
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import workdir as W
import spence_lib as L
import comp_est as CE
import lowe_sz as LZ

KINDS = ["L", "L", "A", "A"]                 # A, B symmetric; C, D arbitrary
AL = (1 + 1j) / 2
U4 = np.array([1, 1j, -1, -1j])


# ---------------------------------------------------------------------------------------------------- structure
class TStruct:
    def __init__(self, mods, mats, qmods, phi, T, name=""):
        """mods: G = Z_{m_1} x ...; mats: generators of H (integer matrices on coordinates); phi: coords -> Q coords;
        T: characters of Q as index tuples k, chi_k(a) = exp(2 pi i sum_j k_j a_j / q_j)"""
        self.G, self.Q, self.name = L.AbGroup(mods), L.AbGroup(qmods), name
        G, Q = self.G, self.Q
        self.mods, self.qmods = [int(x) for x in mods], [int(x) for x in qmods]
        self.mats = [np.asarray(M, np.int64) for M in mats]
        perms = [G.perm_from_matrix(M) for M in self.mats]
        self.stA = L.Structure(G, perms)
        self.stL = L.Structure(G, perms + [G.neg()])
        assert np.array_equal(self.stA.cls, self.stL.cls)
        self.C = self.stA.C
        self.qlab = Q.index(np.array([phi(c) for c in G.coords]))
        assert len(np.unique(self.qlab)) == Q.v, "phi is not onto"
        self.m, self.q = G.v, Q.v
        self.T = [tuple(int(x) for x in k) for k in T]
        self.t = len(self.T)
        self.N = self.m + self.t
        gam = np.array([4 * sum(self.chi(k, a) for k in self.T) for a in range(Q.v)])
        assert np.allclose(gam.imag, 0), "T is not closed under conjugation"
        self.Gam = np.rint(gam.real).astype(np.int64)
        Tg = -self.Gam[self.qlab]
        Tg[0] = 4 * self.m
        self.Tg = Tg                                         # the target function on G
        tv = np.zeros(self.C, np.int64)
        for c in range(self.C):
            vals = np.unique(Tg[self.stA.cls == c])
            if len(vals) != 1:
                raise ValueError("H does not preserve the target (K and T)")
            tv[c] = vals[0]
        self.tvec = tv
        self.one_in_T = tuple([0] * len(self.qmods)) in self.T
        self.rsum2 = 4 * self.t if self.one_in_T else 4 * self.N

    def chi(self, k, a):
        return np.exp(2j * np.pi * sum(kj * aj / qj for kj, aj, qj in zip(k, self.Q.coords[a], self.qmods)))

    def st(self, kind):
        return self.stL if kind == "L" else self.stA

    def patterns(self):
        """row sums (rA >= rB >= 0, rC >= rD >= 0), sum of squares rsum2, parity of m, nonempty pools"""
        par, R, out = self.m % 2, int(math.isqrt(self.rsum2)), []
        for ra in range(par, R + 1, 2):
            for rb in range(par, ra + 1, 2):
                for rc in range(par, R + 1, 2):
                    rest = self.rsum2 - ra * ra - rb * rb - rc * rc
                    rd = math.isqrt(max(rest, 0))
                    if rest < 0 or rd * rd != rest or rd > rc or rd % 2 != par:
                        continue
                    rows = (ra, rb, rc, rd)
                    if all(self.st(k).pool(r) > 0 for k, r in zip(KINDS, rows)):
                        out.append(rows)
        return out

    def gauss(self, rows, K=1500, seed=1):
        """plain Gaussian estimate (independent classes, lattice spacings, one class dropped for the row-sum identity)"""
        rng = np.random.default_rng(seed)
        lN, mus, vars_, lats = 0.0, [], [], []
        for kind, r in zip(KINDS, rows):
            st = self.st(kind)
            n = st.pool(r) / (2 if r == 0 else 1)
            if n < 1:
                return -math.inf
            lN += log10(n)
            V = st.class_values(st.sample(r, K, rng))
            mus.append(V.mean(0)); vars_.append(V.var(0))
            Vi = np.rint(V).astype(np.int64)
            lat = []
            for c in range(self.C):
                u = np.unique(Vi[:, c])
                g = int(np.gcd.reduce(np.diff(u))) if len(u) > 1 else 0
                lat.append((g, int(u[0] % g) if g else int(u[0])))
            lats.append(lat)
        lp, dropped = 0.0, False
        for c in range(self.C):
            mu = sum(x[c] for x in mus); var = sum(x[c] for x in vars_)
            g, res = 0, 0
            for lat in lats:
                g = gcd(g, lat[c][0]); res += lat[c][1]
            tg = int(self.tvec[c])
            if var < 1e-9:
                if abs(tg - mu) > 1e-6:
                    return -math.inf
                continue
            if g and (tg - res) % g:
                return -math.inf
            if not dropped:
                dropped = True; continue
            lp += -0.5 * log10(2 * pi * var) - (tg - mu) ** 2 / (2 * var) / math.log(10) + log10(max(g, 1))
        return lN + lp

    def describe(self):
        return dict(group=self.mods, quotient=self.qmods, T=[list(k) for k in self.T],
                    phi=[int(x) for x in self.qlab], automorphisms=[M.tolist() for M in self.mats])


class TComp:
    """compressions to Q: fibre classes (orbit images in Q), weights, the compressed target, feasible tuples and the
    conditioned lattice-Gaussian estimate (as search/gcomp.py, with per-class targets)"""

    def __init__(self, ts):
        self.ts, self.Q, self.qlab, self.d = ts, ts.Q, ts.qlab, ts.Q.v
        Q = self.Q
        self.add = np.array([[Q.index(Q.coords[a] + Q.coords[b]) for b in range(self.d)] for a in range(self.d)])
        self.fc, self.fid, self.w = {}, {}, {}
        for kind in ("A", "L"):
            st = ts.st(kind)
            cid, classes, fid, w = {}, [], [], []
            for o in range(st.O):
                img = frozenset(self.qlab[st.orb == o].tolist())
                if img not in cid:
                    cid[img] = len(classes); classes.append(sorted(img))
                fid.append(cid[img])
                assert st.sizes[o] % len(img) == 0
                w.append(int(st.sizes[o]) // len(img))
            if sorted(a for c in classes for a in c) != list(range(self.d)):
                raise ValueError("orbit images do not partition Q")
            self.fc[kind], self.fid[kind], self.w[kind] = classes, np.array(fid), np.array(w)
        self.tgt = np.array([ts.Tg[self.qlab == a].sum() for a in range(self.d)], np.int64)
        self.cache = {}

    def compressions(self, kind, r):
        key = ("c", kind, r)
        if key in self.cache:
            return self.cache[key]
        cc = []
        for f in range(len(self.fc[kind])):
            cnt = {0: 1}
            for o in np.flatnonzero(self.fid[kind] == f):
                wo, nxt = int(self.w[kind][o]), {}
                for s, c in cnt.items():
                    nxt[s + wo] = nxt.get(s + wo, 0) + c
                    nxt[s - wo] = nxt.get(s - wo, 0) + c
                cnt = nxt
            cc.append(cnt)
        sizes = [len(c) for c in self.fc[kind]]
        F = len(cc)
        R = [set() for _ in range(F + 1)]
        R[F] = {0}
        for f in range(F - 1, -1, -1):
            R[f] = {sizes[f] * v + s for v in cc[f] for s in R[f + 1]}
        out = {}
        def rec(f, rem, vals, n):
            if f == F:
                if rem == 0:
                    out[tuple(vals)] = n
                return
            for v, c in cc[f].items():
                rest = rem - sizes[f] * v
                if rest in R[f + 1]:
                    vals.append(v); rec(f + 1, rest, vals, n * c); vals.pop()
        rec(0, r, [], 1)
        self.cache[key] = out
        return out

    def expand(self, kind, vals):
        x = np.zeros(self.d, np.int64)
        for f, O in enumerate(self.fc[kind]):
            x[O] = vals[f]
        return x

    def paf(self, x):
        return np.array([int(np.dot(x, x[self.add[:, s]])) for s in range(self.d)])

    def groups(self, kind, r):
        key = ("g", kind, r)
        if key not in self.cache:
            G = {}
            for x in self.compressions(kind, r):
                G.setdefault(tuple(int(t) for t in self.paf(self.expand(kind, x))), []).append(x)
            self.cache[key] = G
        return self.cache[key]

    def feasible(self, rows):
        """[(members per block, multiplicity)]: groups of compressions (equal compressed PAF) solving the compressed
        equation; every combination of members is a feasible tuple"""
        Gs = [self.groups(k, r) for k, r in zip(KINDS, rows)]
        if any(not G for G in Gs):
            return []
        keys = [list(G) for G in Gs]
        left = {}
        for i, j in itertools.product(range(len(keys[0])), range(len(keys[1]))):
            left.setdefault(tuple(int(a + b) for a, b in zip(keys[0][i], keys[1][j])), []).append((i, j))
        out = []
        for k, l in itertools.product(range(len(keys[2])), range(len(keys[3]))):
            s = tuple(int(t - a - b) for t, a, b in zip(self.tgt, keys[2][k], keys[3][l]))
            for i, j in left.get(s, []):
                members = [Gs[0][keys[0][i]], Gs[1][keys[1][j]], Gs[2][keys[2][k]], Gs[3][keys[3][l]]]
                mult = 1
                for mm in members:
                    mult *= len(mm)
                out.append((members, mult))
        return out

    def sample(self, kind, vals, K, rng):
        """K uniform sign vectors (over the orbits) with the given fibre-class values"""
        st = self.ts.st(kind)
        Y = np.empty((K, st.O))
        for f, v in enumerate(vals):
            idx = np.flatnonzero(self.fid[kind] == f)
            ws = [int(self.w[kind][o]) for o in idx]
            n, W_ = len(ws), sum(ws)
            Tb = np.zeros((n + 1, 2 * W_ + 1)); Tb[n, W_] = 1.0
            for i in range(n - 1, -1, -1):
                w = ws[i]
                Tb[i, w:] += Tb[i + 1, :2 * W_ + 1 - w]
                Tb[i, :2 * W_ + 1 - w] += Tb[i + 1, w:]
            rem = np.full(K, v, np.int64)
            for i in range(n):
                w = ws[i]
                a, b = rem - w + W_, rem + w + W_
                cp = np.where((a >= 0) & (a <= 2 * W_), Tb[i + 1, np.clip(a, 0, 2 * W_)], 0.0)
                cm = np.where((b >= 0) & (b <= 2 * W_), Tb[i + 1, np.clip(b, 0, 2 * W_)], 0.0)
                plus = rng.random(K) * (cp + cm) < cp
                Y[:, idx[i]] = np.where(plus, 1.0, -1.0)
                rem = rem - np.where(plus, w, -w)
        return Y

    def stats(self, kind, x, r, K, rng):
        key = ("s", kind, x)
        if key not in self.cache:
            st = self.ts.st(kind)
            V = st.class_values(self.sample(kind, x, K, rng))
            Vi = np.rint(V).astype(np.int64)
            U = np.unique(Vi, axis=0)
            self.cache[key] = (self.compressions(kind, r)[x], V.mean(0), np.cov(V.T), Vi[0],
                               CE.echelon((U[1:] - U[0])[:120]))
        return self.cache[key]

    def estimate(self, tup, rows, K, rng):
        """log10 P(class sums = target | compressions), log10 of the number of block tuples with them"""
        C = self.ts.C
        mu, S, base, gens, ln = np.zeros(C), np.zeros((C, C)), np.zeros(C, np.int64), [], 0.0
        for kind, x, r in zip(KINDS, tup, rows):
            cnt, m_, c_, b_, B_ = self.stats(kind, x, r, K, rng)
            ln += log10(cnt); mu += m_; S += np.atleast_2d(c_); base += b_; gens += [list(rr) for rr in B_]
        B = CE.echelon(gens) if gens else []
        covol = math.exp(0.5 * math.log(CE.bareiss_det([[sum(p * q for p, q in zip(r1, r2)) for r2 in B]
                                                        for r1 in B]))) if B else 1.0
        t = self.ts.tvec
        if not CE.in_lattice(B, t - base):
            return -math.inf, ln
        w, Uv = np.linalg.eigh(S)
        keep = w > 1e-7 * max(w.max(), 1e-12)
        delta = t - mu
        proj = Uv[:, keep] @ (Uv[:, keep].T @ delta)
        if np.linalg.norm(delta - proj) > 1e-6 * (1 + np.linalg.norm(delta)):
            return -math.inf, ln
        qq = float(np.sum((Uv[:, keep].T @ delta) ** 2 / w[keep]))
        lp = (log10(covol) - 0.5 * int(keep.sum()) * log10(2 * pi) - 0.5 * float(np.sum(np.log10(w[keep])))
              - qq / 2 / math.log(10))
        return lp, ln


def tuples(tc, rows, K, rng, maxtup=400):
    """every feasible compression tuple with its conditioned estimate, best first"""
    out = []
    for members, mult in tc.feasible(rows):
        for tup in itertools.product(*members):
            lp, ln = tc.estimate(tup, rows, K, rng)
            if lp > -math.inf:
                out.append((lp + ln, tup, ln, lp))
            if len(out) >= maxtup:
                break
    out.sort(key=lambda x: -x[0])
    return out


def triage(ts, K=300, seed=1, maxtup=150):
    rng = np.random.default_rng(seed)
    tc = TComp(ts)
    res = []
    for rows in ts.patterns():
        g0 = ts.gauss(rows)
        fs = tc.feasible(rows)
        if not fs:
            res.append((rows, g0, None, None)); continue
        tot, bt = 0.0, None
        for members, mult in fs[:maxtup]:
            tup = tuple(mm[0] for mm in members)
            lp, ln = tc.estimate(tup, rows, K, rng)
            if lp == -math.inf:
                continue
            tot += mult * 10 ** (lp + ln)
            if bt is None or lp + ln > bt[0]:
                bt = (lp + ln, tup, [tc.compressions(k, r)[x] for k, r, x in zip(KINDS, rows, tup)])
        if len(fs) > maxtup and tot > 0:
            tot *= len(fs) / maxtup
        e = log10(tot) - sum(1 for r in rows if r == 0) * log10(2) if tot > 0 else -math.inf
        res.append((rows, g0, (e, len(fs)), bt))
    return res


# ---------------------------------------------------------------------------------------------------- families
def units(n):
    return [u for u in range(1, n) if gcd(u, n) == 1]


def blockdiag(*Ms):
    n = sum(M.shape[0] for M in Ms)
    out, i = np.zeros((n, n), np.int64), 0
    for M in Ms:
        k = M.shape[0]; out[i:i + k, i:i + k] = M; i += k
    return out


def make(N, family, gen):
    """family z4cyc: G = Z_(2^a) x Z_o (cyclic, m = N - 3 = 2^a o, a >= 2), Q = Z_4 (first coordinate mod 4),
    T = {1, chi, chi-bar}, H = <diag(gen)>"""
    if family == "z4cyc":
        m = N - 3
        o = m
        while o % 2 == 0:
            o //= 2
        two = m // o
        assert two >= 4, "needs 4 | N - 3"
        mods = [two, o] if o > 1 else [two]
        M = np.diag([int(g) for g in gen][:len(mods)]).astype(np.int64)
        return TStruct(mods, [M], [4], lambda c: (int(c[0]) % 4,), [(0,), (1,), (3,)],
                       f"Z{two}xZ{o} Q=Z4 H=<diag{tuple(gen)}>")
    if family in ("v4n", "v4t"):
        # G = Z_2^2 x Z_o (N - 3 = 4 o), Q = Z_2^2 (the first two coordinates), T = the three nontrivial characters
        # (v4n) or {1, chi_01, chi_10} (v4t); H = <blockdiag(GL2[gen[0]], gen[1])>
        m = N - 3
        assert m % 4 == 0 and (m // 4) % 2 == 1, "needs N - 3 = 4 o with o odd"
        o = m // 4
        M = blockdiag(np.array(GL2[int(gen[0])], np.int64), np.array([[int(gen[1])]], np.int64)) if o > 1 else \
            np.array(GL2[int(gen[0])], np.int64)
        T = [(0, 1), (1, 0), (1, 1)] if family == "v4n" else [(0, 0), (0, 1), (1, 0)]
        return TStruct([2, 2, o] if o > 1 else [2, 2], [M], [2, 2], lambda c: (int(c[0]), int(c[1])), T,
                       f"Z2^2xZ{o} Q=V4 T={'nontrivial' if family == 'v4n' else '1,c01,c10'} H=<{GL2[int(gen[0])]},{gen[1]}>")
    raise ValueError(family)


GL2 = [[[1, 0], [0, 1]], [[0, 1], [1, 0]], [[1, 1], [0, 1]], [[1, 0], [1, 1]], [[0, 1], [1, 1]], [[1, 1], [1, 0]]]


def family_gens(N, family):
    if family in ("v4n", "v4t"):
        o = (N - 3) // 4
        for i in range(len(GL2)):
            for u in (units(o) if o > 1 else [1]):
                yield (i, u)
        return
    if family == "z4cyc":
        m = N - 3
        o = m
        while o % 2 == 0:
            o //= 2
        two = m // o
        for a in units(two):
            for b in (units(o) if o > 1 else [1]):
                yield (a, b)


def perm_closure(perms, v, cap=100000):
    idp = tuple(range(v))
    seen, frontier = {idp}, [idp]
    gens = [tuple(int(x) for x in p) for p in perms]
    while frontier and len(seen) < cap:
        nf = []
        for p in frontier:
            for g in gens:
                q = tuple(g[i] for i in p)
                if q not in seen:
                    seen.add(q); nf.append(q)
        frontier = nf
    return seen


def rank(N, family, omin, omax, K=300):
    seen = set()
    for gen in family_gens(N, family):
        try:
            ts = make(N, family, gen)
        except ValueError:
            continue                                   # H does not preserve the target
        G = ts.G
        grp = perm_closure([G.perm_from_matrix(M) for M in ts.mats], G.v)
        if tuple(int(x) for x in G.neg()) in grp:
            continue                                   # -1 in H: A, B would be all there is
        key = frozenset(grp)
        if key in seen:
            continue
        seen.add(key)
        if not (omin <= ts.stA.O <= omax):
            continue
        head = f"N={N} {ts.name} |H|={len(grp)} O_A={ts.stA.O} O_L={ts.stL.O} C={ts.C}"
        for rows, g0, fe, bt in triage(ts, K):
            if fe is None:
                print(f"{head} rows {rows}: no feasible compression on Q (plain {g0:.2f})", flush=True)
            elif bt:
                print(f"{head} rows {rows}: {fe[1]} feasible groups, conditioned logE {fe[0]:.2f} (plain {g0:.2f}); "
                      f"best tuple {bt[1]} logE {bt[0]:.2f} pools {['%.2g' % p for p in bt[2]]}", flush=True)
            else:
                print(f"{head} rows {rows}: {fe[1]} feasible groups, no estimable tuple (plain {g0:.2f})", flush=True)


# ---------------------------------------------------------------------------------------------------- border
def devidx(G):
    """diff[g, h] = index of h - g, so that x[diff] is the group-developed matrix X[g, h] = x(h - g)"""
    v, r = G.v, len(G.mods)
    return G.index((G.coords[None, :, :] - G.coords[:, None, :]).reshape(-1, r)).reshape(v, v)


def szarray(G, A, B, C, D):
    """Szöllősi's array over the abelian group G (spence_lib.AbGroup) of complex sequences, as a complex matrix"""
    v = G.v
    diff = devidx(G)
    Rf = np.zeros((v, v)); Rf[np.arange(v), G.neg()] = 1
    X = AL * np.asarray(A)[diff] + AL.conjugate() * np.asarray(B)[diff]
    Y = AL * np.asarray(C)[diff] + AL.conjugate() * np.asarray(D)[diff] @ Rf
    return np.block([[X, Y], [Y.conj().T, -X.conj().T]])


def mu4(P, s=1.0):
    P = P / s
    return (np.abs(np.abs(P) - 1) < 1e-9) & ((np.abs(P.real) < 1e-9) | (np.abs(P.imag) < 1e-9))


def vt_vectors(ts):
    """mu_4-valued functions on Q with Fourier support in T"""
    Q = ts.Q
    others = [k for k in itertools.product(*[range(q) for q in ts.qmods]) if tuple(k) not in ts.T]
    ch = [np.array([ts.chi(k, a) for a in range(Q.v)]) for k in others]
    out = []
    for f in itertools.product(range(4), repeat=Q.v):
        x = U4[list(f)]
        if all(abs(np.vdot(c, x)) < 1e-9 for c in ch):
            out.append(x)
    return np.array(out)


def solve_border(ts, Mc, maxstart=400):
    """quaternary Xb (2t x 2t), Rh (2t x 2q), Ch (2q x 2t) with rows of Rh and columns of Ch in V_T (both halves),
    orthogonal, and Xb = -Rh Mc^* Ch / 2q quaternary; depth-first over Rh rows with the compatible columns kept"""
    q, n = ts.q, 2 * ts.t
    V1 = vt_vectors(ts)
    allv = np.array([np.concatenate([a, b]) for a in V1 for b in V1])
    vecs = allv[np.abs(allv[:, 0] - 1) < 1e-9]           # rows and columns may be normalised (unit multiples)
    ortho = np.abs(vecs @ vecs.conj().T) < 1e-9
    comp = mu4(vecs @ Mc.conj().T @ vecs.T, 2 * q)
    def clique(cands):
        def rec(ch, pool):
            if len(ch) == n:
                return ch
            if len(ch) + len(pool) < n:
                return None
            for i_, a in enumerate(pool):
                r = rec(ch + [a], [b for b in pool[i_ + 1:] if ortho[a, b]])
                if r:
                    return r
            return None
        return rec([], list(cands))
    us = np.flatnonzero(comp.sum(1) >= n)
    def rec(ch, pool):
        V = np.flatnonzero(comp[ch].all(0))
        if len(V) < n:
            return None
        if len(ch) == n:
            c = clique(V)
            return (ch, c) if c else None
        for i_, a in enumerate(pool):
            r = rec(ch + [a], [b for b in pool[i_ + 1:] if ortho[a, b]])
            if r:
                return r
        return None
    for u0 in us[:maxstart]:
        r = rec([u0], [u for u in us if u != u0 and ortho[u0, u]])
        if r:
            Rh, Ch = vecs[r[0]], vecs[r[1]].T
            return -(Rh @ Mc.conj().T @ Ch) / (2 * q), Rh, Ch
    return None


def to_exp(M):
    """quaternary complex matrix -> exponents e (entry i^e)"""
    M = np.asarray(M)
    return np.where(np.abs(M - 1) < 1e-9, 0, np.where(np.abs(M - 1j) < 1e-9, 1,
                    np.where(np.abs(M + 1) < 1e-9, 2, 3))).astype(int)


def compressions_of(ts, seqs):
    return [np.bincount(ts.qlab, weights=np.asarray(x, float), minlength=ts.q) for x in seqs]


def build(ts, seqs, border):
    """K as integer real and imaginary parts from the core sequences and the border (exponent matrices)"""
    G, m, q, t = ts.G, ts.m, ts.q, ts.t
    A, B, C, D = (np.asarray(x, np.int64) for x in seqs)
    diff = devidx(G)
    Rf = np.zeros((m, m), np.int64); Rf[np.arange(m), G.neg()] = 1
    XA, XB, XC, XD = A[diff], B[diff], C[diff], D[diff] @ Rf
    Xr, Xi = (XA + XB) // 2, (XA - XB) // 2
    Yr, Yi = (XC + XD) // 2, (XC - XD) // 2
    Kcr = np.block([[Xr, Yr], [Yr.T, -Xr.T]])
    Kci = np.block([[Xi, Yi], [-Yi.T, Xi.T]])
    re = np.array([1, 0, -1, 0]); im = np.array([0, 1, 0, -1])
    Xb, Rh, Ch = (np.asarray(e, np.int64) for e in border)
    cols = np.concatenate([ts.qlab, q + ts.qlab])                         # (i, g) -> i q + phi(g)
    Rr, Ri = re[Rh][:, cols], im[Rh][:, cols]
    Cr, Ci = re[Ch][cols, :], im[Ch][cols, :]
    Kr = np.block([[re[Xb], Rr], [Cr, Kcr]])
    Ki = np.block([[im[Xb], Ri], [Ci, Kci]])
    return Kr, Ki


def check(Kr, Ki):
    n = Kr.shape[0]
    if not np.all(np.abs(Kr) + np.abs(Ki) == 1):
        return False
    I = np.eye(n, dtype=np.int64)
    return bool(np.array_equal(Kr @ Kr.T + Ki @ Ki.T, n * I) and not (Ki @ Kr.T - Kr @ Ki.T).any())


def core_ok(ts, seqs):
    G = ts.G
    tot = np.zeros(G.v, np.int64)
    for x in seqs:
        x = np.asarray(x, np.int64)
        for s in range(G.v):
            tot[s] += int(np.dot(x, x[G.translate(s)]))
    return bool(np.array_equal(tot, ts.Tg))


def complete(ts, seqs):
    """border for a core solution; returns (K parts, border exponents) or None"""
    ys = compressions_of(ts, seqs)
    Mc = szarray(ts.Q, *ys)
    bs = solve_border(ts, Mc)
    if bs is None:
        return None
    border = [to_exp(M) for M in bs]
    Kr, Ki = build(ts, seqs, border)
    return (Kr, Ki), border


# ---------------------------------------------------------------------------------------------------- search
def enum_block(ts, tc, kind, r, vals, center, Wm, tau, maxkeep, tag, threads):
    st = ts.st(kind)
    d = W.WORK / "tmp"; d.mkdir(exist_ok=True)
    inp = d / f"tbsz_{ts.N}_{tag}_{os.getpid()}.in"; out = d / f"tbsz_{ts.N}_{tag}_{os.getpid()}.out"
    with open(inp, "wb") as f:
        f.write(struct.pack("4i", st.O, ts.C, int(r), int(min(maxkeep, 2 ** 31 - 1))))
        f.write(np.asarray(st.sizes, np.int32).tobytes()); f.write(np.rint(st.S).astype(np.int32).tobytes())
        f.write(np.asarray(center, np.float32).tobytes()); f.write(np.asarray(Wm, np.float32).tobytes())
        f.write(struct.pack("f", float(tau)))
        f.write(struct.pack("i", len(vals)))
        f.write(np.asarray(tc.fid[kind], np.int32).tobytes()); f.write(np.asarray(tc.w[kind], np.int32).tobytes())
        f.write(np.asarray(vals, np.int32).tobytes())
    p = subprocess.run([W.bin("lowe_enum"), str(inp), str(out), str(threads)], capture_output=True, text=True,
                       check=True)
    inp.unlink()
    with open(out, "rb") as f:
        n = struct.unpack("q", f.read(8))[0]
        masks = np.frombuffer(f.read(8 * n), np.uint64)
        vals_ = np.frombuffer(f.read(2 * n * ts.C), np.int16).reshape(n, ts.C)
        f.read(4 * n)
    out.unlink()
    return masks, vals_, p.stdout.strip()


def plan_tuple(ts, tc, tup, rows, Lmax, K, rng):
    """conditioned means/covariances, energy thresholds and retention per block (search/lowe_comp.py)"""
    C, t = ts.C, ts.tvec.astype(np.float64)
    pools = [tc.compressions(k, r)[x] for k, r, x in zip(KINDS, rows, tup)]
    mus, sigs = [], []
    for k, x in zip(KINDS, tup):
        V = ts.st(k).class_values(tc.sample(k, x, K, rng))
        mus.append(V.mean(0)); sigs.append(np.cov(V.T))
    Sig = sum(sigs)
    _, Sig_pinv = LZ.psd_pinv_sqrt(Sig)
    plan = []
    for k in range(4):
        Lk = pools[k] if pools[k] <= Lmax else Lmax
        Wk, Smp = LZ.psd_pinv_sqrt(Sig - sigs[k])
        zk = t - (sum(mus) - mus[k])
        V = ts.st(KINDS[k]).class_values(tc.sample(KINDS[k], tup[k], K, rng))
        Ep = np.einsum("ij,jk,ik->i", V - zk, Smp, V - zk)
        frac = Lk / pools[k]
        tau = float(np.quantile(Ep, min(1.0, frac))) if frac < 1 else 1e30
        ck = mus[k] + sigs[k] @ Sig_pinv @ (t - sum(mus))
        Ck = sigs[k] - sigs[k] @ Sig_pinv @ sigs[k]
        w, U = np.linalg.eigh((Ck + Ck.T) / 2); w = np.clip(w, 0, None)
        Xs = ck + (rng.standard_normal((K, C)) * np.sqrt(w)) @ U.T
        Es = np.einsum("ij,jk,ik->i", Xs - zk, Smp, Xs - zk)
        plan.append((Lk, tau, float(np.mean(Es <= tau)) if frac < 1 else 1.0, Wk, zk))
    return pools, plan


def search(ts, rows, which=(0,), Lmax=3e7, nsig=2e5, nS=5, K=4000, threads=4, log=None, maxsol=5, seed=1,
           estimate_only=False):
    t0 = time.time()
    rng = np.random.default_rng(seed)
    tc = TComp(ts)
    tl = tuples(tc, rows, 1500, rng)
    tot = sum(10 ** e for e, *_ in tl)
    print(f"N = {ts.N} {ts.name} |H-orbits| A {ts.stA.O} L {ts.stL.O} classes {ts.C} rows {rows}: {len(tl)} feasible "
          f"compression tuples with positive estimate, conditioned logE {log10(tot) if tot else -math.inf:.2f} "
          f"[{time.time() - t0:.0f}s]", flush=True)
    found = []
    for i in which:
        if i >= len(tl):
            break
        e, tup = tl[i][0], tl[i][1]
        pools, plan = plan_tuple(ts, tc, tup, rows, Lmax, K, rng)
        lret = sum(log10(max(p[2], 1e-12)) for p in plan)
        print(f"tuple {i}: {tup}: pools {[f'{p:.3g}' for p in pools]}, logE {e:.2f}, kept "
              f"{[f'{p[0]:.3g}' for p in plan]}, retention {[f'{p[2]:.3g}' for p in plan]}, logE in lists "
              f"{e + lret:.2f}", flush=True)
        if estimate_only:
            continue
        lists, done = [], {}
        for k in range(4):
            Lk, tau, ret, Wk, zk = plan[k]
            key = (KINDS[k], tup[k], Lk, tau, tuple(np.round(zk, 6)))
            if key in done:
                lists.append(done[key]); continue
            masks, vals, sline = enum_block(ts, tc, KINDS[k], rows[k], tup[k], zk, Wk, tau,
                                            Lk if Lk < pools[k] else 2 ** 31 - 1, "ABCD"[k], threads)
            print(f"  block {'ABCD'[k]}: {len(masks)} kept; {sline} [{time.time() - t0:.0f}s]", flush=True)
            lists.append((masks, vals)); done[key] = (masks, vals)
        if min(len(x[0]) for x in lists) == 0:
            continue
        order = [0, 2, 1, 3]
        keepdir = W.WORK / "lowe"; keepdir.mkdir(exist_ok=True)
        tag = f"tb{ts.N}_{'_'.join(map(str, ts.mods))}_{abs(hash(ts.name)) % 10 ** 8}_{'_'.join(map(str, rows))}_t{i}"
        inp = keepdir / f"kt4_{tag}.bin"
        np.savez(keepdir / f"masks_{tag}.npz", *[lists[k][0] for k in range(4)])
        with open(inp, "wb") as f:
            f.write(struct.pack("5i", *[len(lists[k][0]) for k in order], ts.C))
            f.write(ts.tvec.astype(np.int16).tobytes())
            for k in order:
                f.write(np.ascontiguousarray(lists[k][1], np.int16).tobytes())
        var = sum(np.asarray(lists[k][1][::max(1, len(lists[k][1]) // 200000)], np.float64).var(0) for k in order)
        Sx = [int(c) for c in np.argsort(-var)[:nS]]
        sf = Path(str(inp) + ".sig")
        n, rad = LZ.sigma_file(lists, Sx, ts.tvec, nsig, sf)
        print(f"  kt4: merge classes {Sx}, {n} split values by the solution law (radius {rad:.2f})", flush=True)
        with open(str(inp) + ".progress", "w") as errf:
            p = subprocess.run([W.bin("kt4"), str(inp), ",".join(map(str, Sx)), "@" + str(sf), str(maxsol)],
                               stdout=subprocess.PIPE, stderr=errf, text=True)
        tail = [l for l in p.stdout.splitlines() if l.startswith("TOTAL")]
        sols = [list(map(int, l.split()[1:5])) for l in p.stdout.splitlines() if l.startswith("SOL")]
        print(f"  kt4: {len(sols)} solutions; {tail[0] if tail else 'exit ' + str(p.returncode)} "
              f"[{time.time() - t0:.0f}s]", flush=True)
        for s in sols:
            idx = [0] * 4
            for pos, k in enumerate(order):
                idx[k] = s[pos]
            seqs = []
            for k in range(4):
                mk = int(lists[k][0][idx[k]])
                y = np.array([1 if (mk >> o) & 1 else -1 for o in range(ts.st(KINDS[k]).O)], np.int64)
                seqs.append(y[ts.st(KINDS[k]).orb])
            if not core_ok(ts, seqs):
                print("  core check FAILED", flush=True); continue
            res = complete(ts, seqs)
            if res is None:
                print("  core solution without a border", flush=True); continue
            (Kr, Ki), border = res
            ok = check(Kr, Ki)
            print(f"  {'verified' if ok else 'FAILED'}: complex Hadamard matrix of order {2 * ts.N} [{time.time() - t0:.0f}s]",
                  flush=True)
            if ok:
                rec = dict(N=ts.N, structure=ts.describe(), name=ts.name, rows=list(rows),
                           sequences=["".join("+" if v > 0 else "-" for v in x) for x in seqs],
                           border=dict(Xb=border[0].tolist(), Rh=border[1].tolist(), Ch=border[2].tolist()))
                found.append(rec)
                if log:
                    with open(log, "a") as f:
                        f.write(json.dumps(rec) + "\n")
        if found:
            break
    return found


# ---------------------------------------------------------------------------------------------------- certificates
def text(Kr, Ki):
    e = np.where(Kr == 1, 0, np.where(Ki == 1, 1, np.where(Kr == -1, 2, 3)))
    return "".join("".join(map(str, row)) + "\n" for row in e)


def certify(hits, out_dir=None, k=0):
    """write complex-hadamard/C<2N>.json from the k-th record of a hits file: the structure is rebuilt from the
    record, the border is solved again for the core sequences, and K is checked exactly"""
    import hashlib
    rec = [json.loads(l) for l in open(hits) if l.strip()][k]
    s = rec["structure"]
    G, Q = L.AbGroup(s["group"]), L.AbGroup(s["quotient"])
    lab = np.array(s["phi"], np.int64)
    ts = TStruct(s["group"], s["automorphisms"], s["quotient"],
                 lambda c: tuple(int(x) for x in Q.coords[lab[int(G.index(np.asarray(c)[None, :])[0])]]),
                 [tuple(t) for t in s["T"]], rec.get("name", ""))
    seqs = [np.array([1 if c == "+" else -1 for c in x], np.int64) for x in rec["sequences"]]
    assert core_ok(ts, seqs), "core condition fails"
    res = complete(ts, seqs)
    assert res is not None, "no border"
    (Kr, Ki), border = res
    assert check(Kr, Ki), "not a complex Hadamard matrix"
    txt = text(Kr, Ki)
    cert = {"order": 2 * ts.N, "construction": "szollosi-tbordered", "N": ts.N, "group": s["group"],
            "quotient": s["quotient"], "T": s["T"], "phi": s["phi"], "automorphisms": s["automorphisms"],
            "rows": rec["rows"], "sequences": rec["sequences"],
            "border": {key: [list(map(int, r)) for r in b] for key, b in zip(("Xb", "Rh", "Ch"), border)},
            "sha256": hashlib.sha256(txt.encode()).hexdigest()}
    out = Path(out_dir or Path(__file__).resolve().parent.parent / "complex-hadamard") / f"C{2 * ts.N}.json"
    out.write_text(json.dumps(cert, indent=1) + "\n")
    print(f"wrote {out}: order {2 * ts.N}, sha256 {cert['sha256']}")
    return out


# ---------------------------------------------------------------------------------------------------- tests
def test():
    """exhaustive cores on small groups (Q = Z_4, T = {1, chi, chi-bar}), each completed and checked"""
    import itertools as it
    for mods, N in (([8], 11), ([2, 4], 11), ([12], 15)):
        if len(mods) == 1:
            ts = TStruct(mods, [np.eye(1, dtype=np.int64)], [4], lambda c: (int(c[0]) % 4,), [(0,), (1,), (3,)])
        else:
            ts = TStruct(mods, [np.eye(2, dtype=np.int64)], [4], lambda c: (int(c[1]) % 4,), [(0,), (1,), (3,)])
        G, m = ts.G, ts.m
        syms = [x for x in it.product((1, -1), repeat=m) if all(x[g] == x[G.neg()[g]] for g in range(m))]
        arbs = list(it.product((1, -1), repeat=m))
        P = {}
        for x in arbs:
            xa = np.array(x)
            P[x] = tuple(int(np.dot(xa, xa[G.translate(s)])) for s in range(m))
        left = {}
        for a in syms:
            for b in syms:
                left.setdefault(tuple(u + v for u, v in zip(P[a], P[b])), []).append((a, b))
        tg = tuple(int(x) for x in ts.Tg)
        n = ok = 0
        for c in arbs:
            for d in arbs:
                for a, b in left.get(tuple(x - u - v for x, u, v in zip(tg, P[c], P[d])), []):
                    n += 1
                    if n % 2000 == 1 and ok < 3:
                        res = complete(ts, [a, b, c, d])
                        ok += res is not None and check(*res[0])
        print(f"G = {mods}, N = {N}: {n} core solutions; {ok} completed and verified (order {2 * N})", flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=["test", "rank", "search", "certify"])
    ap.add_argument("N", type=int, nargs="?")
    ap.add_argument("--family", default="z4cyc")
    ap.add_argument("--gen")
    ap.add_argument("--rows")
    ap.add_argument("--tuple", default="0")
    ap.add_argument("--omin", type=int, default=12); ap.add_argument("--omax", type=int, default=60)
    ap.add_argument("--L", type=float, default=3e7); ap.add_argument("--nsig", type=float, default=2e5)
    ap.add_argument("--threads", type=int, default=4); ap.add_argument("--K", type=int, default=4000)
    ap.add_argument("--log"); ap.add_argument("--estimate-only", action="store_true")
    ap.add_argument("--hits", help="certify: the hits file (JSON lines)"); ap.add_argument("--k", type=int, default=0)
    a = ap.parse_args()
    if a.mode == "test":
        test()
    elif a.mode == "certify":
        certify(a.hits, k=a.k)
    elif a.mode == "rank":
        rank(a.N, a.family, a.omin, a.omax)
    else:
        ts = make(a.N, a.family, [int(x) for x in a.gen.split(",")])
        rows = tuple(int(x) for x in a.rows.split(","))
        which = [int(x) for x in a.tuple.split(",")]
        search(ts, rows, which, a.L, a.nsig, 5, a.K, a.threads, a.log, estimate_only=a.estimate_only)


if __name__ == "__main__":
    main()
