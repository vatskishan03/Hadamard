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
    python3 search/tbsz.py rank N --family z4cyc|z3cyc|z6s|z6t|v4n|v4t|all3 [--omin 12 --omax 60]
    python3 search/tbsz.py search N --family z4cyc --gen a,b --rows 2,0,2,2 [--tuple i] [--L 3e7] [--nsig 2e5]
    python3 search/tbsz.py search N --family all3 --name "NAME" --rows ...       (NAME as printed by rank)
    python3 search/tbsz.py rankg N --group 2,2,3,13 --quot z6:0,2         (any G given by prime-power factors)
    python3 search/tbsz.py search N --family group --group 2,2,3,13 --quot z6:0,2 --mat "1,0,0,0;1,1,0,0;0,0,2,0;0,0,0,3" --rows 2,0,4,0
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
    """compressions to a quotient: fibre classes (orbit images), weights, the compressed target, feasible tuples and
    the conditioned lattice-Gaussian estimate (as search/gcomp.py, with per-class targets). The quotient is Q, or a
    finer cyclic quotient Z_d of a cyclic G (d a multiple of q: cond=d), which splits the pools further; then
    `toQ` maps its labels to Q. A per-coordinate quotient (cond = [d_1, ...]) may also be coarser than Q or cross it
    (toQ = None): the target sums over its fibres still condition the estimate."""

    def __init__(self, ts, cond=None, parent=None, pqlab=None):
        self.ts = ts
        if cond is None or (isinstance(cond, int) and cond == ts.q and len(ts.mods) == 1):
            self.Q, self.qlab, self.d = ts.Q, ts.qlab, ts.Q.v
            self.toQ = np.arange(ts.q)
        else:
            # the quotient prod_i Z_(d_i) of G = prod_i Z_(m_i), coordinate-wise reduction (d_i | m_i); an integer d
            # stands for the cyclic quotient Z_d of a cyclic G
            ds = [int(cond)] if not isinstance(cond, (list, tuple)) else [int(x) for x in cond]
            assert len(ds) == len(ts.mods) and all(mi % di == 0 for mi, di in zip(ts.mods, ds)), "cond: d_i | m_i"
            self.Q = L.AbGroup(ds)
            self.qlab = self.Q.index(ts.G.coords % np.array(ds))
            self.d = self.Q.v
            self.toQ = np.zeros(self.d, np.int64)
            for a in range(self.d):
                vals = set(ts.qlab[self.qlab == a].tolist())
                if len(vals) != 1:                     # a quotient that does not refine Q: estimates only, the
                    self.toQ = None                    # border is then solved for each solution
                    break
                self.toQ[a] = vals.pop()
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
        # parent: the compressions of the four blocks on Q (a feasible tuple there); the compressions on this finer
        # quotient are then only those that sum to them over the fibres of toQ
        # pqlab: the labels of a coarser quotient (default Q) on which `parent` lives
        self.parent = None
        if parent is not None:
            pq = ts.qlab if pqlab is None else np.asarray(pqlab)
            self.toP = np.zeros(self.d, np.int64)
            for a in range(self.d):
                vals = set(pq[self.qlab == a].tolist())
                assert len(vals) == 1, "the quotient does not refine the parent quotient"
                self.toP[a] = vals.pop()
            self.parent = [np.asarray(x, np.int64) for x in parent]

    def par(self, b):
        return None if self.parent is None else self.parent[b]

    def classcounts(self, kind):
        """per fibre class, {value: number of sign vectors on its orbits}"""
        key = ("cc", kind)
        if key not in self.cache:
            cc = []
            for f in range(len(self.fc[kind])):
                cnt = {0: 1}
                for o in np.flatnonzero(self.fid[kind] == f):
                    wo, nxt = int(self.w[kind][o]), {}
                    for s_, c in cnt.items():
                        nxt[s_ + wo] = nxt.get(s_ + wo, 0) + c
                        nxt[s_ - wo] = nxt.get(s_ - wo, 0) + c
                    cnt = nxt
                cc.append(cnt)
            self.cache[key] = cc
        return self.cache[key]

    def count(self, kind, x):
        """the number of sign vectors (over the orbits) with the fibre-class values x"""
        cc, out = self.classcounts(kind), 1.0
        for f, v in enumerate(x):
            out *= float(cc[f].get(int(v), 0))
        return out

    def compressions(self, kind, r, par=None):
        """every achievable compression (tuple of fibre-class values) with row sum r and the number of sign vectors
        giving it: per class a DP over its orbits, then the product over classes pruned by reachability (numpy).
        par: compressions on Q that the result must refine (sums over the fibres of toQ)"""
        key = ("c", kind, r, None if par is None else tuple(int(x) for x in par))
        if key in self.cache:
            return self.cache[key]
        cc0 = self.classcounts(kind)
        F = len(cc0)
        sizes0 = [len(c) for c in self.fc[kind]]
        order = list(range(F))
        if par is not None:
            par = np.asarray(par, np.int64)
            Cn0 = np.zeros((F, len(par)), np.int64)
            for f, O in enumerate(self.fc[kind]):
                for p_ in O:
                    Cn0[f, self.toP[p_]] += 1
            # classes ordered by the Q-points they meet, so that the points close early
            order = sorted(range(F), key=lambda f: (int(np.flatnonzero(Cn0[f])[0]), f))
            Cn = Cn0[order]
        cc = [cc0[f] for f in order]
        sizes = [sizes0[f] for f in order]
        R = [None] * (F + 1)
        R[F] = np.array([0], np.int64)
        for f in range(F - 1, -1, -1):
            v = np.array(sorted(cc[f]), np.int64)
            R[f] = np.unique((sizes[f] * v[:, None] + R[f + 1][None, :]).ravel())
        # energy bound: sum_k sum_a y_k(a)^2 = tgt(0), and every point of the quotient carries at least
        # (fibre size mod 2)^2 in each block, so one block has energy <= tgt(0) - 3 * (number of odd fibres)
        emax = int(self.tgt[0]) - 3 * int(np.sum(np.bincount(self.qlab, minlength=self.d) % 2))
        mins = [min(x * x for x in cc[f]) for f in range(F)]
        remE = [0] * (F + 1)
        for f in range(F - 1, -1, -1):
            remE[f] = remE[f + 1] + sizes[f] * mins[f]
        if par is not None:
            vmax = np.array([max(abs(x) for x in cc[f]) for f in range(F)], np.int64)
            # rmax[f]: the largest |contribution| of classes f.. to each Q-point; open[f]: points still met by them
            rmax = np.zeros((F + 1, len(par)), np.int64)
            for f in range(F - 1, -1, -1):
                rmax[f] = rmax[f + 1] + Cn[f] * vmax[f]
            openp = np.zeros((F + 1, len(par)), bool)
            for f in range(F - 1, -1, -1):
                openp[f] = openp[f + 1] | (Cn[f] > 0)
            S = np.zeros((1, len(par)), np.int64)
        X = np.zeros((1, 0), np.int64)
        part = np.zeros(1, np.int64)
        pe = np.zeros(1, np.int64)
        cnt = np.ones(1, np.float64)
        for f in range(F):
            v = np.array(sorted(cc[f]), np.int64)
            v = v[sizes[f] * v * v + remE[f + 1] <= emax] if len(v) else v
            c = np.array([float(cc[f][x]) for x in v])
            P = (part[:, None] + sizes[f] * v[None, :]).ravel()
            E = (pe[:, None] + sizes[f] * (v * v)[None, :]).ravel()
            keep = np.isin(r - P, R[f + 1]) & (E + remE[f + 1] <= emax)
            if par is not None:
                for a in np.flatnonzero(Cn[f]):
                    Sa = (S[:, a][:, None] + Cn[f, a] * v[None, :]).ravel()
                    if openp[f + 1, a]:
                        keep &= np.abs(par[a] - Sa) <= rmax[f + 1, a]
                    else:
                        keep &= Sa == par[a]
            i, j = np.divmod(np.flatnonzero(keep), len(v))
            X = np.concatenate([X[i], v[j][:, None]], 1)
            part, pe, cnt = P[keep], E[keep], cnt[i] * c[j]
            if par is not None:
                S = S[i] + Cn[f][None, :] * v[j][:, None]
        if par is not None:                            # back to the class order
            inv = np.argsort(np.array(order))
            X = X[:, inv]
        out = dict(zip(map(tuple, X.tolist()), cnt.tolist()))
        self.cache[key] = out
        return out

    def expand(self, kind, vals):
        x = np.zeros(self.d, np.int64)
        for f, O in enumerate(self.fc[kind]):
            x[O] = vals[f]
        return x

    def paf(self, x):
        return np.array([int(np.dot(x, x[self.add[:, s]])) for s in range(self.d)])

    def groups(self, kind, r, par=None):
        """compressions grouped by their compressed PAF (vectorised)"""
        key = ("g", kind, r, None if par is None else tuple(int(x) for x in par))
        if key not in self.cache:
            comps = list(self.compressions(kind, r, par))
            G = {}
            if comps:
                fidx = np.zeros(self.d, np.int64)
                for f, O in enumerate(self.fc[kind]):
                    fidx[O] = f
                X = np.array(comps, np.int64)
                for i0 in range(0, len(X), 1 << 20):
                    E = X[i0:i0 + (1 << 20)][:, fidx]
                    P = np.stack([(E * E[:, self.add[:, s]]).sum(1) for s in range(self.d)], 1)
                    U, inv = np.unique(P, axis=0, return_inverse=True)
                    inv = np.asarray(inv).ravel()
                    keys = [tuple(int(t) for t in u) for u in U]
                    for j, x in zip(inv, comps[i0:i0 + (1 << 20)]):
                        G.setdefault(keys[j], []).append(x)
            self.cache[key] = G
        return self.cache[key]

    def feasible(self, rows):
        """[(members per block, multiplicity)]: groups of compressions (equal compressed PAF) solving the compressed
        equation; every combination of members is a feasible tuple. The join is vectorised: the sums of A and B keys
        are encoded and sorted, the targets minus C and D keys are looked up"""
        Gs = [self.groups(k, r, self.par(b)) for b, (k, r) in enumerate(zip(KINDS, rows))]
        if any(not G for G in Gs):
            return []
        names = [list(G) for G in Gs]
        nAB, nCD = len(names[0]) * len(names[1]), len(names[2]) * len(names[3])
        if min(nAB, nCD) > 3e7:
            raise ValueError(f"too many compression groups to join ({len(names[0])} x {len(names[1])}, "
                             f"{len(names[2])} x {len(names[3])})")
        keys = [np.array(n_, np.int64).reshape(-1, self.d) for n_ in names]
        V = max(int(np.abs(k).max()) for k in keys) * 3 + int(np.abs(self.tgt).max()) + 1
        bits = int(V).bit_length() + 1
        if self.d * bits <= 62:
            Bb = 1 << bits
            def enc(M):
                M = M + Bb // 2
                out = np.zeros(M.shape[:-1], np.int64)
                for j in range(self.d):
                    out = out * Bb + M[..., j]
                return out
        else:                                          # rows as raw bytes: equality matching only
            def enc(M):
                M = np.ascontiguousarray(M, np.int64)
                return M.reshape(-1, self.d).view(np.dtype((np.void, 8 * self.d))).reshape(M.shape[:-1])
        # when one pair side is small, iterate over its pairs, vectorise over a third list and look up the fourth
        n4 = [len(k) for k in keys]
        plans = [(n4[p] * n4[q] * n4[r], (p, q, r, s_)) for p, q, r, s_ in
                 ((0, 1, 2, 3), (0, 1, 3, 2), (2, 3, 0, 1), (2, 3, 1, 0))]
        cost, (p, q, r, s_) = min(plans)
        if cost < nAB + nCD and n4[p] * n4[q] <= 200000:
            Ws = enc(keys[s_])
            ordW = np.argsort(Ws, kind="stable")
            Wsorted = Ws[ordW]
            out = []
            for i in range(n4[p]):
                for j in range(n4[q]):
                    E = enc(self.tgt[None, :] - keys[p][i][None, :] - keys[q][j][None, :] - keys[r])
                    lo = np.searchsorted(Wsorted, E, side="left"); hi = np.searchsorted(Wsorted, E, side="right")
                    for z in np.flatnonzero(hi > lo):
                        for pos in range(lo[z], hi[z]):
                            idx = [0] * 4
                            idx[p], idx[q], idx[r], idx[s_] = i, j, int(z), int(ordW[pos])
                            members = [Gs[b][names[b][idx[b]]] for b in range(4)]
                            mult = 1
                            for mm in members:
                                mult *= len(mm)
                            out.append((members, mult))
            return out
        # otherwise the smaller of the two pair sides is encoded and sorted, the other is streamed in chunks
        ab_small = nAB <= nCD
        if ab_small:
            L_ = enc(keys[0][:, None, :] + keys[1][None, :, :]).ravel()          # (i, j) -> i * n1 + j
            nsm = len(keys[1])
        else:
            L_ = enc(self.tgt[None, None, :] - keys[2][:, None, :] - keys[3][None, :, :]).ravel()
            nsm = len(keys[3])
        order = np.argsort(L_, kind="stable")
        Ls = L_[order]
        big = (2, 3) if ab_small else (0, 1)
        nb = len(keys[big[1]])
        out = []
        chunk = max(1, 4_000_000 // max(nb, 1))
        for k0 in range(0, len(keys[big[0]]), chunk):
            if ab_small:
                E = enc(self.tgt[None, None, :] - keys[2][k0:k0 + chunk, None, :] - keys[3][None, :, :]).ravel()
            else:
                E = enc(keys[0][k0:k0 + chunk, None, :] + keys[1][None, :, :]).ravel()
            lo = np.searchsorted(Ls, E, side="left"); hi = np.searchsorted(Ls, E, side="right")
            for f in np.flatnonzero(hi > lo):
                k, l = divmod(int(f), nb)
                for pos in range(lo[f], hi[f]):
                    i, j = divmod(int(order[pos]), nsm)
                    idx = (i, j, k0 + k, l) if ab_small else (k0 + k, l, i, j)
                    members = [Gs[b][names[b][idx[b]]] for b in range(4)]
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
            self.cache[key] = (self.count(kind, x), V.mean(0), np.cov(V.T), Vi[0], CE.diff_lattice(Vi, rng))
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
    fs = tc.feasible(rows)
    ntup = sum(mult for _, mult in fs)
    print(f"  on the quotient of order {tc.d}: {len(fs)} feasible groups, {ntup} feasible tuples"
          + (f" (estimating the first {maxtup})" if ntup > maxtup else ""), flush=True)
    for members, mult in fs:
        for tup in itertools.product(*members):
            lp, ln = tc.estimate(tup, rows, K, rng)
            if lp > -math.inf:
                out.append((lp + ln, tup, ln, lp))
            if len(out) >= maxtup:
                break
    out.sort(key=lambda x: -x[0])
    return out


def triage(ts, K=300, seed=1, maxtup=150, cond=None):
    rng = np.random.default_rng(seed)
    tc = TComp(ts, cond)
    res = []
    for rows in ts.patterns():
        fs = tc.feasible(rows)
        if not fs:
            res.append((rows, float("nan"), None, None)); continue          # exact: no feasible compression
        g0 = ts.gauss(rows)
        tot, bt = 0.0, None
        for members, mult in fs[:maxtup]:
            tup = tuple(mm[0] for mm in members)
            lp, ln = tc.estimate(tup, rows, K, rng)
            if lp == -math.inf:
                continue
            tot += mult * 10 ** (lp + ln)
            if bt is None or lp + ln > bt[0]:
                bt = (lp + ln, tup, [tc.count(k, x) for k, x in zip(KINDS, tup)])
        if len(fs) > maxtup and tot > 0:
            tot *= len(fs) / maxtup
        e = log10(tot) - sum(1 for r in rows if r == 0) * log10(2) if tot > 0 else -math.inf
        res.append((rows, g0, (e, len(fs)), bt))
    return res


def invariant_refinements(ts, maxorder=None):
    """coordinate-wise quotients prod Z_(d_i) (d_i | m_i) of G whose kernel is H-invariant and lies in ker(phi),
    of order between 2q and maxorder (default 12q), smallest first: the levels at which compression tuples on Q can
    be refined"""
    G, mods = ts.G, ts.mods
    maxorder = maxorder or 12 * ts.q
    zero = int(ts.qlab[G.index(np.zeros((1, len(mods)), np.int64))[0]])
    out = []
    for ds in itertools.product(*[[d for d in range(1, m + 1) if m % d == 0] for m in mods]):
        order = int(np.prod(ds))
        if not (2 * ts.q <= order <= maxorder):
            continue
        ok = True
        for i, di in enumerate(ds):
            g = np.zeros(len(mods), np.int64); g[i] = di % mods[i]
            if ts.qlab[G.index(g[None, :])[0]] != zero:
                ok = False; break
            for M in ts.mats:
                img = (M @ g) % np.array(mods)
                if any(img[j] % ds[j] for j in range(len(mods))):
                    ok = False; break
            if not ok:
                break
        if ok:
            out.append((order, list(ds)))
    return [d for _, d in sorted(out)]


def refined_triage(ts, cond, K=300, seed=1, ntop=3, maxfine=150, emin=-3.0):
    """per row pattern: the best `ntop` coarse tuples on Q (one per feasible group, swaps identified) refined on the
    quotient `cond`: (rows, coarse logE of those tuples, refined logE, number of them alive on cond, number refined)"""
    rng = np.random.default_rng(seed)
    tc0 = TComp(ts)
    res = []
    for rows in ts.patterns():
        fs = tc0.feasible(rows)
        if not fs:
            continue
        cands, seen = [], set()
        for members, mult in fs:
            tup = tuple(mm[0] for mm in members)
            k0 = {tup}
            if rows[2] == rows[3]:
                k0 |= {(x[0], x[1], x[3], x[2]) for x in list(k0)}
            if rows[0] == rows[1]:
                k0 |= {(x[1], x[0], x[2], x[3]) for x in list(k0)}
            if k0 & seen:
                continue
            seen |= k0
            lp, ln = tc0.estimate(tup, rows, K, rng)
            if lp > -math.inf:
                cands.append((lp + ln, tup))
        cands.sort(key=lambda x: -x[0])
        if not cands or cands[0][0] < emin:              # not worth refining
            continue
        ce, fe, alive = 0.0, 0.0, 0
        for e0, tup0 in cands[:ntop]:
            ce += 10 ** e0
            tcf = TComp(ts, cond, parent=[tc0.expand(k, x) for k, x in zip(KINDS, tup0)], pqlab=tc0.qlab)
            ffs = tcf.feasible(rows)
            tot = 0.0
            for members, mult in ffs[:maxfine]:
                lp, ln = tcf.estimate(tuple(mm[0] for mm in members), rows, K, rng)
                if lp > -math.inf:
                    tot += mult * 10 ** (lp + ln)
            if len(ffs) > maxfine and tot > 0:
                tot *= len(ffs) / maxfine
            fe += tot
            alive += tot > 0
        res.append((rows, log10(ce) if ce else -math.inf, log10(fe) if fe else -math.inf, alive,
                    min(ntop, len(cands))))
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
    T = {1, chi, chi-bar}, H = <diag(gen)>. (T = {chi, chi^2, chi^3} is the same family: multiplying the four
    sequences by the real character chi^2 o phi maps T to chi^2 T, and the border along with it.)"""
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
    if family in ("z6s", "z6t"):
        # G = Z_m cyclic (m = N - 5, 6 | m), Q = Z_6 (x mod 6), T = all characters but the sign (z6s) or but the
        # trivial one (z6t), t = 5: N = 5 (mod 6)
        m = N - 5
        assert m % 6 == 0, "needs N = 5 (mod 6)"
        T = [(k,) for k in range(6) if k != (3 if family == "z6s" else 0)]
        return TStruct([m], [np.array([[int(gen[0])]], np.int64)], [6], lambda c: (int(c[0]) % 6,), T,
                       f"Z{m} Q=Z6 T={'all-sgn' if family == 'z6s' else 'all-1'} H=<{gen[0]}>")
    if family == "z3cyc":
        # G = Z_m cyclic (m = N - 3, 3 | m), Q = Z_3 (x mod 3), T = all characters: target -12 on 3Z_m, 0 off it
        m = N - 3
        assert m % 3 == 0, "needs 3 | N"
        M = np.array([[int(gen[0])]], np.int64)
        return TStruct([m], [M], [3], lambda c: (int(c[0]) % 3,), [(0,), (1,), (2,)], f"Z{m} Q=Z3 H=<{gen[0]}>")
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


def all3(N):
    """every t = 3 shape for N = 3 (mod 4): (name, mods, matrix, qmods, phi, T), H generated by the matrix"""
    m = N - 3
    o = m
    while o % 2 == 0:
        o //= 2
    k = (m // o).bit_length() - 1                      # m = 2^k o
    U = units(o) if o > 1 else [1]
    z4 = [(0,), (1,), (3,)]
    v4 = {"nontrivial": [(0, 1), (1, 0), (1, 1)], "1,c01,c10": [(0, 0), (0, 1), (1, 0)]}
    def odd(M, u):
        return blockdiag(np.asarray(M, np.int64), np.array([[u]], np.int64)) if o > 1 else np.asarray(M, np.int64)
    oddm = [o] if o > 1 else []
    two = 2 ** k
    if k >= 2:                                         # cyclic 2-part, Q = Z_4
        for a in units(two):
            for u in U:
                yield (f"Z{two}xZ{o} Q=Z4 a={a} u={u}", [two] + oddm, odd([[a]], u), [4],
                       (lambda c: (int(c[0]) % 4,)), z4)
    if k >= 3:                                         # Z_2 x Z_(2^(k-1)): Q = Z_4 on the second factor, or Z_2^2
        h = two // 2
        for b in (0, 1):
            for d in units(h):
                for u in U:
                    yield (f"Z2xZ{h}xZ{o} Q=Z4 b={b} d={d} u={u}", [2, h] + oddm, odd([[1, b], [0, d]], u), [4],
                           (lambda c: (int(c[1]) % 4,)), z4)
        for b in (0, 1):
            for c in sorted({0, h // 2}):
                for d in units(h):
                    for u in U:
                        for tn, T in v4.items():
                            yield (f"Z2xZ{h}xZ{o} Q=V4 T={tn} b={b} c={c} d={d} u={u}", [2, h] + oddm,
                                   odd([[1, b], [c, d]], u), [2, 2], (lambda cc: (int(cc[0]), int(cc[1]) % 2)), T)
    if k == 2:                                         # Z_2^2 x Z_o, Q = Z_2^2
        for i, M2 in enumerate(GL2):
            for u in U:
                for tn, T in v4.items():
                    yield (f"Z2^2xZ{o} Q=V4 T={tn} M={M2} u={u}", [2, 2] + oddm, odd(M2, u), [2, 2],
                           (lambda c: (int(c[0]), int(c[1]))), T)
    if k == 3:                                         # Z_2^3 x Z_o, Q = Z_2^2 = Z_2^3 / <e_3>
        for M2 in GL2:
            for e, f in itertools.product((0, 1), repeat=2):
                M3 = [[M2[0][0], M2[0][1], 0], [M2[1][0], M2[1][1], 0], [e, f, 1]]
                for u in U:
                    for tn, T in v4.items():
                        yield (f"Z2^3xZ{o} Q=V4 T={tn} M={M3} u={u}", [2, 2, 2] + oddm, odd(M3, u), [2, 2],
                               (lambda c: (int(c[0]), int(c[1]))), T)


def rank_all3(N, omin, omax, K=300):
    seen = set()
    for name, mods, M, qmods, phi, T in all3(N):
        G = L.AbGroup(mods)
        try:
            perm = G.perm_from_matrix(M)
        except ValueError:
            continue
        grp = perm_closure([perm], G.v)
        if tuple(int(x) for x in G.neg()) in grp:
            continue
        key = (tuple(mods), tuple(qmods), tuple(map(tuple, T)), frozenset(grp))
        if key in seen:
            continue
        seen.add(key)
        try:
            ts = TStruct(mods, [M], qmods, phi, T, name)
        except (ValueError, AssertionError):
            continue
        if not (omin <= ts.stA.O <= omax):
            continue
        head = f"N={N} {name} |H|={len(grp)} O_A={ts.stA.O} O_L={ts.stL.O} C={ts.C}"
        try:
            res = triage(ts, K)
        except ValueError as ex:
            print(f"{head}: {ex}", flush=True); continue
        for rows, g0, fe, bt in res:
            if fe is None:
                print(f"{head} rows {rows}: no feasible compression on Q (plain {g0:.2f})", flush=True)
            elif bt:
                print(f"{head} rows {rows}: {fe[1]} feasible groups, conditioned logE {fe[0]:.2f} (plain {g0:.2f}); "
                      f"best tuple {bt[1]} logE {bt[0]:.2f} pools {['%.2g' % p for p in bt[2]]}", flush=True)
            else:
                print(f"{head} rows {rows}: {fe[1]} feasible groups, no estimable tuple (plain {g0:.2f})", flush=True)


def family_gens(N, family):
    if family == "z3cyc":
        for u in units(N - 3):
            yield (u,)
        return
    if family in ("z6s", "z6t"):
        for u in units(N - 5):
            yield (u,)
        return
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


# ---------------------------------------------------------------------------------------------------- any group
QUOTS = {"z6": ([6], [(0,), (1,), (2,), (4,), (5,)]),          # T = all characters but the sign (t = 5)
         "z4": ([4], [(0,), (1,), (3,)]),                       # T = {1, chi, chi-bar}
         "v4": ([2, 2], [(0, 0), (0, 1), (1, 0)]),              # T = {1, chi_01, chi_10}
         "z3": ([3], [(0,), (1,), (2,)]),                       # T = all characters
         "e33": ([3, 3], [(0, 0), (1, 0), (2, 0), (0, 1), (0, 2)]),  # T = {1} + two independent order-3 pairs (t = 5)
         "one": ([1], [(0,)]),                                  # t = 1: the bordered array (BORDERED.md)
         # unions of the character groups of two factors, V_T = {g(a_1) + h(a_2)}: t = q_1 + q_2 - 1
         "u44": ([4, 4], [(0, 0), (1, 0), (2, 0), (3, 0), (0, 1), (0, 2), (0, 3)]),   # t = 7
         "u24": ([2, 4], [(0, 0), (1, 0), (0, 1), (0, 2), (0, 3)])}                   # t = 5


def quotient(spec):
    """'z6:i,j' (Z2 from coordinate i, Z3 from coordinate j), 'z4:i', 'v4:i,j', 'z3:i', 'e33:i,j' (Z3 x Z3 from
    coordinates i and j) -> (qmods, phi, T)"""
    kind, idx = (spec.split(":") + [""])[:2]
    ix = [int(x) for x in idx.split(",") if x]
    qmods, T = QUOTS[kind]
    if kind == "one":
        phi = lambda c: (0,)
    elif kind == "z6":
        phi = lambda c: ((3 * (int(c[ix[0]]) % 2) + 4 * (int(c[ix[1]]) % 3)) % 6,)
    elif kind == "z4":
        phi = lambda c: (int(c[ix[0]]) % 4,)
    elif kind == "v4":
        phi = lambda c: (int(c[ix[0]]) % 2, int(c[ix[1]]) % 2)
    elif kind == "e33":
        phi = lambda c: (int(c[ix[0]]) % 3, int(c[ix[1]]) % 3)
    elif kind == "u44":
        phi = lambda c: (int(c[ix[0]]) % 4, int(c[ix[1]]) % 4)
    elif kind == "u24":
        phi = lambda c: (int(c[ix[0]]) % 2, int(c[ix[1]]) % 4)
    else:
        phi = lambda c: (int(c[ix[0]]) % 3,)
    return qmods, phi, T


def block_auts(mods):
    """every automorphism matrix of Z_(m_1) x ... x Z_(m_k), all m_i powers of one prime (M_ij in Z_(m_i), m_i | M_ij m_j)"""
    k = len(mods)
    ranges = [list(range(0, mods[i], mods[i] // gcd(mods[i], mods[j]) if mods[i] > mods[j] else 1))
              for i in range(k) for j in range(k)]
    G = L.AbGroup(mods)
    out = []
    for vals in itertools.product(*ranges):
        M = np.array(vals, np.int64).reshape(k, k)
        try:
            p = G.perm_from_matrix(M)
        except ValueError:
            continue
        if len(set(int(x) for x in p)) == G.v:
            out.append(M)
    return out


def group_auts(mods):
    """automorphisms of prod Z_(m_i) (each m_i a prime power) that act on each Sylow subgroup separately"""
    bydp = {}
    for i, m_ in enumerate(mods):
        p_ = min(d for d in range(2, m_ + 1) if m_ % d == 0)
        bydp.setdefault(p_, []).append(i)
    blocks = list(bydp.values())
    per = [block_auts([mods[i] for i in b]) for b in blocks]
    for combo in itertools.product(*per):
        M = np.zeros((len(mods), len(mods)), np.int64)
        for b, Mb in zip(blocks, combo):
            M[np.ix_(b, b)] = Mb
        yield M


def elementary_stabiliser(p, n, k, others):
    """automorphisms of Z_p^n (coordinates 0..n-1 of a block) that map the span of the first k coordinate forms to
    itself with a monomial action on those k forms: [[D, 0], [V, A]] with D monomial (entries +-1), V in Z_p^((n-k) x k),
    A in GL(n - k, p); combined with the automorphisms in `others` (list of (indices, matrices)) of the other blocks"""
    import itertools as it
    mons = []
    for perm in it.permutations(range(k)):
        for sg in it.product((1, p - 1), repeat=k):
            D = np.zeros((k, k), np.int64)
            for i, (j, e) in enumerate(zip(perm, sg)):
                D[i, j] = e
            mons.append(D)
    r = n - k
    GLr = [np.array(v, np.int64).reshape(r, r) for v in it.product(range(p), repeat=r * r)
           if r and round(np.linalg.det(np.array(v, float).reshape(r, r))) % p] if r else [np.zeros((0, 0), np.int64)]
    Vs = [np.array(v, np.int64).reshape(r, k) for v in it.product(range(p), repeat=r * k)]
    for D in mons:
        for V in Vs:
            for A in GLr:
                M = np.zeros((n, n), np.int64)
                M[:k, :k] = D; M[k:, :k] = V; M[k:, k:] = A
                yield M


def rank_group(N, mods, qspec, omin, omax, K=300, mats=None, cond=None, refine=0):
    """rank the structures H = <M> (M an automorphism preserving the target) on G = prod Z_(m_i), distinct up to an
    invariant of their orbits (sizes and images in Q, for H and <H, -1>); mats: candidate matrices (default: every
    automorphism acting on each Sylow subgroup)"""
    qmods, phi, T = quotient(qspec)
    G, Q = L.AbGroup(mods), L.AbGroup(qmods)
    assert G.v == N - len(T), "|G| must be N - t"
    qlab = Q.index(np.array([phi(c) for c in G.coords]))
    gam = np.array([round((4 * sum(np.exp(2j * pi * sum(kj * aj / qj for kj, aj, qj in zip(k, Q.coords[a], qmods)))
                                   for k in T)).real) for a in range(Q.v)], np.int64)
    tg = gam[qlab]
    neg = tuple(int(x) for x in G.neg())
    negp = np.array(neg)
    seen, seen_inv = set(), set()
    for M in (group_auts(mods) if mats is None else mats):
        p = np.asarray(G.perm_from_matrix(M))
        if not np.array_equal(tg[p], tg):
            continue
        # the structure depends only on the orbits of H = <p>: deduplicate by the orbit partition
        q_ = p.copy(); order = 1
        while not np.array_equal(q_, np.arange(G.v)) and order <= 5000:
            q_ = q_[p]; order += 1
        cyc = np.arange(G.v)
        x = np.arange(G.v)
        for _ in range(order):
            x = p[x]; cyc = np.minimum(cyc, x)
        key = cyc.tobytes()
        if key in seen:
            continue
        seen.add(key)
        pw = np.arange(G.v); hasneg = False
        for _ in range(order):
            pw = p[pw]
            if np.array_equal(pw, negp):
                hasneg = True; break
        if hasneg:
            continue
        grp = range(order)
        stA = L.Structure(G, [p])
        if not (omin <= stA.O <= omax):
            continue
        stL = L.Structure(G, [p, G.neg()])
        inv = tuple(tuple(sorted((int(st.sizes[o]), tuple(sorted(set(qlab[st.orb == o].tolist())))) for o in range(st.O)))
                    for st in (stA, stL))
        if inv in seen_inv:
            continue
        seen_inv.add(inv)
        name = f"G={'x'.join(map(str, mods))} Q={qspec} M={';'.join(','.join(map(str, r)) for r in M.tolist())}"
        ts = TStruct(mods, [M], qmods, phi, T, name)
        head = f"N={N} {name} |H|={len(grp)} O_A={ts.stA.O} O_L={ts.stL.O} C={ts.C}"
        try:
            res = triage(ts, K, cond=cond)
        except (ValueError, MemoryError) as ex:
            print(f"{head}: {ex}", flush=True)
            continue
        for rows, g0, fe, bt in res:
            if fe is None:
                print(f"{head} rows {rows}: no feasible compression on Q (plain {g0:.2f})", flush=True)
            elif bt:
                print(f"{head} rows {rows}: {fe[1]} feasible groups, conditioned logE {fe[0]:.2f} (plain {g0:.2f}); "
                      f"best tuple {bt[1]} logE {bt[0]:.2f} pools {['%.2g' % p for p in bt[2]]}", flush=True)
            else:
                print(f"{head} rows {rows}: {fe[1]} feasible groups, no estimable tuple (plain {g0:.2f})", flush=True)
        if refine:
            levels = invariant_refinements(ts)
            if not levels:
                print(f"{head}: no invariant refinement of Q", flush=True)
                continue
            try:
                for rows, ce, fe2, alive, nref in refined_triage(ts, levels[0], K, ntop=refine):
                    print(f"{head} rows {rows}: REFINED on {'x'.join(map(str, levels[0]))}: best {nref} coarse tuples "
                          f"10^{ce:.2f} -> 10^{fe2:.2f}, {alive} alive", flush=True)
            except (ValueError, MemoryError) as ex:
                print(f"{head}: refinement on {levels[0]} failed: {ex}", flush=True)


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


_VT_CACHE = {}


def vt_vectors(ts):
    """mu_4-valued functions on Q with Fourier support in T (meet in the middle over the two halves of Q: the
    constraints sum_a f(a) conj(chi(a)) = 0, chi not in T, are split into the sums over the halves and matched)"""
    key = (tuple(ts.qmods), tuple(ts.T))
    if key in _VT_CACHE:
        return _VT_CACHE[key]
    Q = ts.Q
    others = [k for k in itertools.product(*[range(q) for q in ts.qmods]) if tuple(k) not in ts.T]
    ch = np.array([[ts.chi(k, a) for a in range(Q.v)] for k in others]).reshape(len(others), Q.v)
    h = Q.v // 2
    XL = np.array(list(itertools.product(range(4), repeat=h)), np.int64).reshape(4 ** h, h)
    XR = np.array(list(itertools.product(range(4), repeat=Q.v - h)), np.int64).reshape(4 ** (Q.v - h), Q.v - h)
    SL = U4[XL] @ ch[:, :h].conj().T if len(others) else np.zeros((len(XL), 0))
    SR = U4[XR] @ ch[:, h:].conj().T if len(others) else np.zeros((len(XR), 0))
    def keys(S):
        Z = np.rint(np.concatenate([S.real, S.imag], 1) * 1e6).astype(np.int64)
        return [z.tobytes() for z in Z]
    table = {}
    for i, kk in enumerate(keys(SL)):
        table.setdefault(kk, []).append(i)
    out = []
    for j, kk in enumerate(keys(-SR)):
        for i in table.get(kk, ()):
            out.append(np.concatenate([U4[XL[i]], U4[XR[j]]]))
    out = np.array(out).reshape(-1, Q.v)
    out = out[np.lexsort(to_exp(out).T[::-1])] if len(out) else out       # the order of the old enumeration
    _VT_CACHE[key] = out
    return out


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


def _clique(V, n, tries=1, rng=None, maxnodes=200000):
    """n mutually orthogonal rows of V: depth-first, the candidate pool filtered by one vectorised inner product per
    chosen row (no Gram matrix); random start orders on later tries"""
    V = np.asarray(V)
    nodes = [0]
    def rec(ch, pool):
        if len(ch) == n:
            return ch
        if len(ch) + len(pool) < n:
            return None
        for k in range(len(pool)):
            nodes[0] += 1
            if nodes[0] > maxnodes:
                return None
            a = pool[k]
            rest = pool[k + 1:]
            if len(rest):
                rest = rest[np.abs(V[rest] @ V[a].conj()) < 1e-9]
            r = rec(ch + [int(a)], rest)
            if r:
                return r
        return None
    for t in range(tries):
        order = np.arange(len(V)) if (t == 0 or rng is None) else rng.permutation(len(V))
        nodes[0] = 0
        r = rec([], order)
        if r:
            return r
    return None


def solve_border2(ts, Mc, nbases=60, seed=1):
    """the border by its Xb rows: for an orthogonal quaternary basis Ch of W = V_T + V_T (columns), every border row
    is u(x) = x Ch^* Mc / 2t for a quaternary x (then u Mc^* Ch / 2q = x); keep the x whose u is quaternary and find
    2t orthogonal ones. The column bases are (v_i, c w_i), (v_i, -c w_i) for orthogonal quaternary bases (v_i), (w_i)
    of V_T and a unit c; up to nbases of them are tried."""
    q, t, n = ts.q, ts.t, 2 * ts.t
    rng = np.random.default_rng(seed)
    V1 = vt_vectors(ts)
    if len(V1) < t:
        return None
    bases = []
    for _ in range(4 * nbases):                       # orthogonal bases of V_T (t-cliques), random start orders
        c = _clique(V1[rng.permutation(len(V1))], t)
        if c is None:
            break
        bases.append(c)
        if len(bases) >= nbases:
            break
    if not bases:
        return None
    perm_bases = []
    for _ in range(4 * nbases):
        o = rng.permutation(len(V1))
        c = _clique(V1[o], t)
        if c is not None:
            perm_bases.append(V1[o][c])
        if len(perm_bases) >= nbases:
            break
    xs = np.array([np.concatenate([[1], U4[list(c)]]) for c in itertools.product(range(4), repeat=n - 1)])
    for b in range(nbases):
        B1 = perm_bases[b % len(perm_bases)]
        B2 = perm_bases[rng.integers(len(perm_bases))]
        c = U4[rng.integers(4)]
        cols = np.array([np.concatenate([B1[i], s_ * c * B2[i]]) for i in range(t) for s_ in (1, -1)]).T
        Uu = xs @ cols.conj().T @ Mc / (2 * t)              # rows u(x)
        good = np.flatnonzero(mu4(Uu).all(1))
        if len(good) < n:
            continue
        cl = _clique(xs[good], n, tries=3, rng=rng)
        if cl:
            return -xs[good][cl], Uu[good][cl], cols
    return None


def solve_border_exh(ts, Mc, maxn=80000, threads=3):
    """exhaustive border search (search/border_exh.c) over the normalised quaternary vectors of W: rows depth-first
    with their common compatible columns as bitsets, then an orthogonal 2t-clique of columns. Returns the border, or
    None when none exists (an exact negative for this Mc), or False when W has more than maxn candidates."""
    q, t = ts.q, ts.t
    V1 = vt_vectors(ts)
    if len(V1) < t:
        return None
    n = (len(V1) // 4) * len(V1)
    if n > maxn:
        return False
    allv = np.array([np.concatenate([a, b]) for a in V1 for b in V1])
    vecs = allv[np.abs(allv[:, 0] - 1) < 1e-9]
    P2 = 2 * (vecs @ Mc.conj().T)
    P = np.rint(P2.real).astype(np.int32), np.rint(P2.imag).astype(np.int32)
    assert np.allclose(P[0], P2.real) and np.allclose(P[1], P2.imag), "Mc is not half-integral"
    d = W.WORK / "tmp"; d.mkdir(exist_ok=True)
    path = d / f"border_exh_{ts.N}_{os.getpid()}.in"
    with open(path, "wb") as f:
        f.write(struct.pack("<5i", len(vecs), 2 * q, 2 * t, 4 * q, 0))
        f.write(to_exp(vecs).astype(np.uint8).tobytes())
        f.write(np.stack([P[0], P[1]], -1).astype("<i4").tobytes())
    try:
        env = dict(os.environ, OMP_NUM_THREADS=str(threads))
        out = subprocess.run([str(W.bin("border_exh")), str(path)], capture_output=True, text=True, env=env,
                             check=True).stdout
    finally:
        path.unlink(missing_ok=True)
    for line in out.splitlines():
        if line.startswith("SOLUTION"):
            tok = line.split()
            rows = [int(x) for x in tok[2:2 + 2 * t]]
            cols = [int(x) for x in tok[3 + 2 * t:3 + 4 * t]]
            Rh, Ch = vecs[rows], vecs[cols].T
            Xb = -(Rh @ Mc.conj().T @ Ch) / (2 * q)
            assert np.allclose(Rh @ Rh.conj().T, 2 * q * np.eye(2 * t)) and np.allclose(Ch.conj().T @ Ch, 2 * q * np.eye(2 * t))
            assert mu4(Xb).all()
            return Xb, Rh, Ch
    return None


def solve_border_restricted(ts, Mc, threads=3):
    """for W too large for the exhaustive search: the first border row is fixed to the all-ones vector of both halves
    (in W when 1 is in T), the columns are restricted to those compatible with it and the rows to those orthogonal to
    it, and that restriction is searched exhaustively (border_exh.c with separate row and column candidates).
    Returns the border or None (no border with this first row)."""
    q, t = ts.q, ts.t
    if tuple([0] * len(ts.qmods)) not in ts.T:
        return None
    V1 = vt_vectors(ts)
    n1 = len(V1)
    ones = np.flatnonzero(np.all(np.abs(V1 - 1) < 1e-9, axis=1))
    if len(ones) == 0:
        return None
    o = int(ones[0])
    # normalised candidates (first entry 1): pairs (a, b) of V1 with a[0] = 1
    A = np.flatnonzero(np.abs(V1[:, 0] - 1) < 1e-9)
    u1 = np.concatenate([V1[o], V1[o]])
    P1 = 2 * (u1 @ Mc.conj().T)
    # column candidates compatible with u1: P1 . (a, b) in 4q mu_4, computed per pair from the two halves
    ha = V1[A] @ P1[:q]; hb = V1 @ P1[q:]
    Z = ha[:, None] + hb[None, :]
    okc = (np.abs(np.abs(Z) - 4 * q) < 1e-6) & ((np.abs(Z.real) < 1e-6) | (np.abs(Z.imag) < 1e-6))
    ia, ib = np.nonzero(okc)
    cols = np.concatenate([V1[A][ia], V1[ib]], 1)
    # row candidates orthogonal to u1 (u1 itself first)
    ga = V1[A] @ V1[o].conj(); gb = V1 @ V1[o].conj()
    Y = ga[:, None] + gb[None, :]
    ra, rb = np.nonzero(np.abs(Y) < 1e-6)
    rows = np.concatenate([u1[None, :], np.concatenate([V1[A][ra], V1[rb]], 1)])
    P = 2 * (rows @ Mc.conj().T)
    Pi = np.rint(P.real).astype(np.int32), np.rint(P.imag).astype(np.int32)
    assert np.allclose(Pi[0], P.real) and np.allclose(Pi[1], P.imag), "Mc is not half-integral"
    d = W.WORK / "tmp"; d.mkdir(exist_ok=True)
    path = d / f"border_rc_{ts.N}_{os.getpid()}.in"
    with open(path, "wb") as f:
        f.write(struct.pack("<7i", -1, len(rows), len(cols), 2 * q, 2 * t, 4 * q, 1))
        f.write(to_exp(rows).astype(np.uint8).tobytes())
        f.write(np.stack([Pi[0], Pi[1]], -1).astype("<i4").tobytes())
        f.write(to_exp(cols).astype(np.uint8).tobytes())
        f.write(struct.pack("<i", 0))
    try:
        env = dict(os.environ, OMP_NUM_THREADS=str(threads))
        out = subprocess.run([str(W.bin("border_exh")), str(path)], capture_output=True, text=True, env=env,
                             check=True).stdout
    finally:
        path.unlink(missing_ok=True)
    for line in out.splitlines():
        if line.startswith("SOLUTION"):
            tok = line.split()
            ri = [int(x) for x in tok[2:2 + 2 * t]]
            ci = [int(x) for x in tok[3 + 2 * t:3 + 4 * t]]
            Rh, Ch = rows[ri], cols[ci].T
            Xb = -(Rh @ Mc.conj().T @ Ch) / (2 * q)
            assert np.allclose(Rh @ Rh.conj().T, 2 * q * np.eye(2 * t)) and np.allclose(Ch.conj().T @ Ch, 2 * q * np.eye(2 * t))
            assert mu4(Xb).all()
            return Xb, Rh, Ch
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
    bs = solve_border_exh(ts, Mc)
    if bs is False:                                    # too many candidates for the exhaustive search
        bs = solve_border_restricted(ts, Mc)
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


def as4(tc):
    """one TComp per block (a list of four when the blocks are listed on different quotients)"""
    return list(tc) if isinstance(tc, (list, tuple)) else [tc] * 4


def plan_tuple(ts, tc, tup, rows, Lmax, K, rng):
    """conditioned means/covariances, energy thresholds and retention per block (search/lowe_comp.py)"""
    C, t = ts.C, ts.tvec.astype(np.float64)
    tcs = as4(tc)
    pools = [tcs[b].count(k, x) for b, (k, x) in enumerate(zip(KINDS, tup))]
    mus, sigs = [], []
    for b, (k, x) in enumerate(zip(KINDS, tup)):
        V = ts.st(k).class_values(tcs[b].sample(k, x, K, rng))
        mus.append(V.mean(0)); sigs.append(np.cov(V.T))
    Sig = sum(sigs)
    _, Sig_pinv = LZ.psd_pinv_sqrt(Sig)
    plan = []
    for k in range(4):
        Lk = pools[k] if pools[k] <= Lmax else Lmax
        Wk, Smp = LZ.psd_pinv_sqrt(Sig - sigs[k])
        zk = t - (sum(mus) - mus[k])
        V = ts.st(KINDS[k]).class_values(tcs[k].sample(KINDS[k], tup[k], K, rng))
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


def tuple_border(ts, tc, tup, borders):
    """the border for the compressions on Q of a tuple (None: proved not to exist), cached by those compressions"""
    tcs = as4(tc)
    ys = [np.bincount(tcs[b].toQ, weights=tcs[b].expand(k, x).astype(float), minlength=ts.q)
          for b, (k, x) in enumerate(zip(KINDS, tup))]
    ykey = tuple(int(round(v)) for y in ys for v in y)
    if ykey not in borders:
        Mc = szarray(ts.Q, *ys)
        bs = solve_border_exh(ts, Mc)
        if bs is False:                                # too many candidates for the exhaustive search
            bs = solve_border_restricted(ts, Mc)
        borders[ykey] = bs
    return borders[ykey]


def search(ts, rows, which=(0,), Lmax=3e7, nsig=2e5, nS=5, K=4000, threads=4, log=None, maxsol=5, seed=1,
           estimate_only=False, cond=None, refine=0, maxtup=400, pcond=None, coarse_ab=False, by_rate=False):
    """refine > 0: the best `refine` tuples on the quotient pcond (default Q) are refined one by one on the finer
    quotient `cond` (their compressions there enumerated under the coarse ones), and the fine tuples of all of them
    are ranked together. coarse_ab: blocks A and B stay on the coarse quotient, so that one join covers every refinement
    of them that goes with given fine C and D (the fine tuples are grouped by their C and D)"""
    t0 = time.time()
    rng = np.random.default_rng(seed)
    borders = {}
    if refine and cond is not None:
        tc0 = TComp(ts, pcond)
        tl0 = tuples(tc0, rows, 1500, rng)
        tl, nref, seen0 = [], 0, set()
        for e0, tup0, *_ in tl0:
            if nref >= refine:
                break
            # swapping C and D (or A and B) when their row sums agree maps solutions to solutions
            k0 = {tup0}
            if rows[2] == rows[3]:
                k0 |= {(x[0], x[1], x[3], x[2]) for x in list(k0)}
            if rows[0] == rows[1]:
                k0 |= {(x[1], x[0], x[2], x[3]) for x in list(k0)}
            if k0 & seen0:
                continue
            seen0 |= k0
            if tc0.toQ is not None and ts.t > 1 and tuple_border(ts, tc0, tup0, borders) is None:
                print(f"  coarse tuple {tup0} (logE {e0:.2f}): no border", flush=True)
                continue
            nref += 1
            tcf = TComp(ts, cond, parent=[tc0.expand(k, x) for k, x in zip(KINDS, tup0)], pqlab=tc0.qlab)
            sub = tuples(tcf, rows, 1500, rng, maxtup)
            tot0 = sum(10 ** e for e, *_ in sub)
            print(f"  coarse tuple {tup0} (logE {e0:.2f}): {len(sub)} fine tuples with positive estimate, total "
                  f"logE {log10(tot0) if tot0 else -math.inf:.2f} [{time.time() - t0:.0f}s]", flush=True)
            tl += [(e, tup, ln, lp, tcf, tc0, tup0) for e, tup, ln, lp in sub]
        tl.sort(key=lambda x: -x[0])
        if coarse_ab:
            groups = {}
            for e, tup, ln, lp, tcf, tc0_, tup0 in tl:
                g = groups.setdefault((id(tcf), tup[2], tup[3]), [0.0, 0, tup, tcf, tc0_, tup0])
                g[0] += 10 ** e; g[1] += 1
            tl = sorted([(log10(sm), (tup0[0], tup0[1], tup[2], tup[3]), n, None, [tc0_, tc0_, tcf, tcf])
                         for sm, n, tup, tcf, tc0_, tup0 in groups.values()], key=lambda x: -x[0])
            print(f"  A and B on the coarse quotient: {len(tl)} groups of fine tuples by C and D, the best "
                  f"{', '.join(f'10^{x[0]:.2f} ({x[2]} tuples)' for x in tl[:3])}", flush=True)
    else:
        tc = TComp(ts, cond)
        tl = [(e, tup, ln, lp, tc, None, None) for e, tup, ln, lp in tuples(tc, rows, 1500, rng, maxtup)]
    if by_rate:
        # expected solutions per unit of join work: the work grows with the product of the C and D pools (plus a fixed
        # cost per tuple for the lists and the split values)
        def rate(x):
            tcs_ = as4(x[4])
            return x[0] - log10(tcs_[2].count(KINDS[2], x[1][2]) * tcs_[3].count(KINDS[3], x[1][3]) + 1e12)
        tl.sort(key=lambda x: -rate(x))
    tot = sum(10 ** e for e, *_ in tl)
    print(f"N = {ts.N} {ts.name} |H-orbits| A {ts.stA.O} L {ts.stL.O} classes {ts.C} rows {rows}: {len(tl)} feasible "
          f"compression tuples with positive estimate, conditioned logE {log10(tot) if tot else -math.inf:.2f} "
          f"[{time.time() - t0:.0f}s]", flush=True)
    found = []
    seen = set()
    done = {}
    for i in which:
        if i >= len(tl):
            break
        e, tup, tc = tl[i][0], tl[i][1], tl[i][4]
        tcs = as4(tc)
        # swapping C and D (equal row sums) or A and B (equal row sums) maps solutions to solutions
        keys = {tup}
        if rows[2] == rows[3]:
            keys |= {(t[0], t[1], t[3], t[2]) for t in list(keys)}
        if rows[0] == rows[1]:
            keys |= {(t[1], t[0], t[2], t[3]) for t in list(keys)}
        if keys & seen:
            print(f"tuple {i}: equivalent to an earlier tuple, skipped", flush=True)
            continue
        seen |= keys
        # the border depends only on the compressions: solve it first (when the quotient refines Q)
        if all(x.toQ is not None for x in tcs):
            bs = tuple_border(ts, tcs, tup, borders)
            if bs is None:
                print(f"tuple {i}: {tup}: no border found, skipped", flush=True)
                continue
        pools, plan = plan_tuple(ts, tc, tup, rows, Lmax, K, rng)
        lret = sum(log10(max(p[2], 1e-12)) for p in plan)
        print(f"tuple {i}: {tup}: pools {[f'{p:.3g}' for p in pools]}, logE {e:.2f}, kept "
              f"{[f'{p[0]:.3g}' for p in plan]}, retention {[f'{p[2]:.3g}' for p in plan]}, logE in lists "
              f"{e + lret:.2f}", flush=True)
        if estimate_only:
            continue
        lists = []
        if len(done) > 8:                              # keep the most recent lists only
            for key_ in list(done)[:-8]:
                del done[key_]
        for k in range(4):
            Lk, tau, ret, Wk, zk = plan[k]
            # a complete list depends only on the block's class values (and is shared by later tuples)
            key = ((id(tcs[k]), KINDS[k], tup[k], rows[k]) if Lk >= pools[k]
                   else (id(tcs[k]), KINDS[k], tup[k], Lk, tau, tuple(np.round(zk, 6))))
            if key in done:
                lists.append(done[key]); continue
            masks, vals, sline = enum_block(ts, tcs[k], KINDS[k], rows[k], tup[k], zk, Wk, tau,
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
        if n == 0:
            continue
        with open(str(inp) + ".progress", "w") as errf:
            p = subprocess.run([W.bin("kt4"), str(inp), ",".join(map(str, Sx)), "@" + str(sf), str(maxsol)],
                               stdout=subprocess.PIPE, stderr=errf, text=True)
        tail = [l for l in p.stdout.splitlines() if l.startswith("TOTAL")]
        sols = [list(map(int, l.split()[1:5])) for l in p.stdout.splitlines() if l.startswith("SOL")]
        if not sols:                                   # the join input and masks are large: keep them only for hits
            for f_ in (inp, sf, Path(str(inp) + ".progress"), keepdir / f"masks_{tag}.npz"):
                try:
                    f_.unlink()
                except OSError:
                    pass
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
    """exhaustive cores on small groups (Q = Z_4, T = {1, chi, chi-bar}; Q = Z_6, T = all characters but the sign),
    samples completed by the border search and checked"""
    import itertools as it
    z6 = [(0,), (1,), (2,), (4,), (5,)]
    for mods, N, qm in (([8], 11, 4), ([2, 4], 11, 4), ([12], 15, 4), ([6], 11, 6), ([12], 17, 6)):
        if qm == 6:
            ts = TStruct(mods, [np.eye(1, dtype=np.int64)], [6], lambda c: (int(c[0]) % 6,), z6)
        elif len(mods) == 1:
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
        cap, step = (3, 2000) if qm == 4 else (2, 97)
        for c in arbs:
            for d in arbs:
                for a, b in left.get(tuple(x - u - v for x, u, v in zip(tg, P[c], P[d])), []):
                    n += 1
                    if n % step == 1 and ok < cap:
                        res = complete(ts, [a, b, c, d])
                        ok += res is not None and check(*res[0])
        print(f"G = {mods}, Q = Z{qm}, N = {N}: {n} core solutions; {ok} completed and verified (order {2 * N})",
              flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=["test", "rank", "rankg", "search", "certify"])
    ap.add_argument("N", type=int, nargs="?")
    ap.add_argument("--family", default="z4cyc")
    ap.add_argument("--gen"); ap.add_argument("--name", help="search --family all3: the structure name from rank")
    ap.add_argument("--rows")
    ap.add_argument("--tuple", default="0")
    ap.add_argument("--omin", type=int, default=12); ap.add_argument("--omax", type=int, default=60)
    ap.add_argument("--L", type=float, default=3e7); ap.add_argument("--nsig", type=float, default=2e5)
    ap.add_argument("--threads", type=int, default=4); ap.add_argument("--K", type=int, default=4000)
    ap.add_argument("--log"); ap.add_argument("--estimate-only", action="store_true")
    ap.add_argument("--hits", help="certify: the hits file (JSON lines)"); ap.add_argument("--k", type=int, default=0)
    ap.add_argument("--refine", type=int, default=0, help="search with --cond: refine the best REFINE tuples on Q one "
                    "by one on the finer quotient")
    ap.add_argument("--by-rate", action="store_true", help="search: take the tuples in order of estimate per join "
                    "work (product of the C and D pools) instead of estimate")
    ap.add_argument("--coarse-ab", action="store_true", help="search --refine: list A and B on the coarse quotient "
                    "and join each fine (C, D) once with all of them")
    ap.add_argument("--pcond", help="search --refine: the coarse quotient (as --cond; default Q)")
    ap.add_argument("--maxtup", type=int, default=400, help="search: estimate at most this many tuples per quotient")
    ap.add_argument("--cond", help="search: condition on a finer quotient: d for Z_d of a cyclic G, or d1,d2,... "
                    "per coordinate (d_i | m_i), refining Q")
    ap.add_argument("--group", help="rankg, search --family group: G = Z_(m_1) x ..., prime-power m_i, e.g. 2,2,3,13")
    ap.add_argument("--quot", help="rankg, search --family group: z6:i,j | z4:i | v4:i,j | z3:i | e33:i,j | u44:i,j | "
                    "u24:i,j (coordinates of phi) | one (t = 1, the bordered array)")
    ap.add_argument("--mat", help="search --family group: the generator of H, rows separated by ';'")
    a = ap.parse_args()
    if a.mode == "test":
        test()
    elif a.mode == "certify":
        certify(a.hits, k=a.k)
    elif a.mode == "rankg":
        cond = None if a.cond is None else [int(x) for x in a.cond.split(",")]
        rank_group(a.N, [int(x) for x in a.group.split(",")], a.quot, a.omin, a.omax, cond=cond, refine=a.refine)
    elif a.mode == "rank":
        if a.family == "all3":
            rank_all3(a.N, a.omin, a.omax)
        else:
            rank(a.N, a.family, a.omin, a.omax)
    else:
        if a.family == "all3":                       # the structure is named as in the rank output
            spec = next(x for x in all3(a.N) if x[0] == a.name)
            ts = TStruct(spec[1], [spec[2]], spec[3], spec[4], spec[5], spec[0])
        elif a.family == "group":                    # as printed by rankg
            mods = [int(x) for x in a.group.split(",")]
            qmods, phi, T = quotient(a.quot)
            M = np.array([[int(v) for v in r.split(",")] for r in a.mat.split(";")], np.int64)
            ts = TStruct(mods, [M], qmods, phi, T, f"G={'x'.join(map(str, mods))} Q={a.quot} M={a.mat}")
        else:
            ts = make(a.N, a.family, [int(x) for x in a.gen.split(",")])
        rows = tuple(int(x) for x in a.rows.split(","))
        which = [int(x) for x in a.tuple.split(",")]
        cond = None if a.cond is None else ([int(x) for x in a.cond.split(",")] if "," in a.cond else int(a.cond))
        search(ts, rows, which, a.L, a.nsig, 5, a.K, a.threads, a.log, estimate_only=a.estimate_only, cond=cond,
               refine=a.refine, maxtup=a.maxtup,
               pcond=None if a.pcond is None else [int(x) for x in a.pcond.split(",")], coarse_ab=a.coarse_ab,
               by_rate=a.by_rate)


if __name__ == "__main__":
    main()
