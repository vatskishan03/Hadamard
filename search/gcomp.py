#!/usr/bin/env python3
"""Compression-conditioned estimates for the bordered Szöllősi array over non-cyclic groups (search/bsz_group.py).

search/comp_est.py conditions on the compressions to a cyclic quotient Z_d. Here G is any abelian group
Z_{m_1} x ... x Z_{m_r} with a symmetry group S of block-diagonal automorphisms (bsz_group.BStruct, e.g. the spin
structures on Z_2^2 x Z_m, N = 4m + 1), and the quotient Q is the product of a subset of the factors (Q = Z_2^2 by
default). The fibre sums x(a) = sum_{q(g) = a} X(g) satisfy sum_X PAF_{x_X} = t_Q on Q with

    t_Q(0) = 4|G| + t (|K| - 1),   t_Q(s) = t |K|  (s != 0),   t = -4,  K = ker(G -> Q),

every S-orbit lies over one S-orbit of Q (a fibre class) and meets each of its fibres in the same number of points,
so a block's compression is a vector of fibre-class values. As in comp_est: exact compressions with counts, the
tuples solving the compressed equation, and per tuple a lattice Gaussian of the class values of blocks sampled with
those compressions. For Q = Z_2^2 the compressed equation says sum_X x̂_X(χ)^2 = 4N for the three non-trivial
characters χ and sum_X r_X^2 = 4.

    python3 search/gcomp.py spin N [N ...] [--maxo 64] [--K 400]     every spin structure of N (G = Z_2^2 x Z_m)
"""
import argparse, itertools, math, sys
from math import log10
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import comp_est as CE
import spence_lib as L
import bsz_group as BG

KINDS = BG.KINDS                                    # L, L, A, A


class GComp:
    def __init__(self, bs, qcoords=(0, 1)):
        G = bs.G
        self.bs = bs
        qcoords = list(qcoords)
        self.Q = L.AbGroup([int(G.mods[i]) for i in qcoords])
        self.d = self.Q.v
        self.qlab = self.Q.index(G.coords[:, qcoords])
        self.K = G.v // self.d
        self.add = np.array([[self.Q.index(self.Q.coords[a] + self.Q.coords[b]) for b in range(self.d)]
                             for a in range(self.d)])
        self.fc, self.fid, self.w = {}, {}, {}
        for kind in ("A", "L"):
            st = bs.st(kind)
            cid, classes, fid, w = {}, [], [], []
            for o in range(st.O):
                img = frozenset(self.qlab[st.orb == o].tolist())
                if img not in cid:
                    cid[img] = len(classes); classes.append(sorted(img))
                fid.append(cid[img])
                assert st.sizes[o] % len(img) == 0
                w.append(int(st.sizes[o]) // len(img))
            seen = sorted(a for c in classes for a in c)
            assert seen == list(range(self.d)), "orbit images do not partition Q"
            self.fc[kind], self.fid[kind], self.w[kind] = classes, np.array(fid), np.array(w)
        self.cache = {}
        t = -4
        self.tgt = np.full(self.d, t * self.K, np.int64)
        self.tgt[0] = 4 * G.v + t * (self.K - 1)

    def compressions(self, kind, r):
        key = ("c", kind, r)
        if key in self.cache:
            return self.cache[key]
        cc = []
        for f in range(len(self.fc[kind])):
            cnt = {0: 1}
            for o in np.flatnonzero(self.fid[kind] == f):
                wo = int(self.w[kind][o]); nxt = {}
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

    def sample(self, kind, vals, K, rng):
        st = self.bs.st(kind)
        Y = np.empty((K, st.O))
        for f, v in enumerate(vals):
            idx = np.flatnonzero(self.fid[kind] == f)
            ws = [int(self.w[kind][o]) for o in idx]
            n = len(ws); W_ = sum(ws)
            T = np.zeros((n + 1, 2 * W_ + 1)); T[n, W_] = 1.0
            for i in range(n - 1, -1, -1):
                w = ws[i]
                T[i, w:] += T[i + 1, :2 * W_ + 1 - w]
                T[i, :2 * W_ + 1 - w] += T[i + 1, w:]
            rem = np.full(K, v, np.int64)
            for i in range(n):
                w = ws[i]
                a = rem - w + W_; b = rem + w + W_
                cp = np.where((a >= 0) & (a <= 2 * W_), T[i + 1, np.clip(a, 0, 2 * W_)], 0.0)
                cm = np.where((b >= 0) & (b <= 2 * W_), T[i + 1, np.clip(b, 0, 2 * W_)], 0.0)
                plus = rng.random(K) * (cp + cm) < cp
                Y[:, idx[i]] = np.where(plus, 1.0, -1.0)
                rem = rem - np.where(plus, w, -w)
        return Y

    def groups(self, kind, r):
        key = ("g", kind, r)
        if key not in self.cache:
            G = {}
            for x in self.compressions(kind, r):
                G.setdefault(tuple(int(t) for t in self.paf(self.expand(kind, x))), []).append(x)
            self.cache[key] = G
        return self.cache[key]

    def feasible(self, rows):
        """[(members per block, multiplicity)] solving the compressed equation (blocks A, B | C, D)"""
        Gs = [self.groups(k, r) for k, r in zip(KINDS, rows)]
        if any(not G for G in Gs):
            return []
        keys = [list(G) for G in Gs]
        left = {}
        for i, j in itertools.product(range(len(keys[0])), range(len(keys[1]))):
            s = tuple(int(a + b) for a, b in zip(keys[0][i], keys[1][j]))
            left.setdefault(s, []).append((i, j))
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

    def stats(self, kind, x, K, rng):
        key = ("s", kind, x)
        if key not in self.cache:
            st = self.bs.st(kind)
            Y = self.sample(kind, x, K, rng)
            V = st.class_values(Y)
            Vi = np.rint(V).astype(np.int64)
            self.cache[key] = (self.compressions(kind, sum(len(c) * v for c, v in zip(self.fc[kind], x)))[x],
                               V.mean(0), np.cov(V.T), Vi[0], CE.diff_lattice(Vi, rng))
        return self.cache[key]

    def estimate(self, tup, K, rng):
        C = self.bs.C
        mu = np.zeros(C); S = np.zeros((C, C)); base = np.zeros(C, np.int64); gens = []; ln = 0.0
        for kind, x in zip(KINDS, tup):
            cnt, m_, c_, b_, B_ = self.stats(kind, x, K, rng)
            ln += log10(cnt); mu += m_; S += np.atleast_2d(c_); base += b_; gens += [list(r) for r in B_]
        B = CE.echelon(gens) if gens else []
        covol = math.exp(0.5 * math.log(CE.bareiss_det([[sum(p * q for p, q in zip(r1, r2)) for r2 in B]
                                                        for r1 in B]))) if B else 1.0
        t = np.full(C, -4, np.int64)
        if not CE.in_lattice(B, t - base):
            return -math.inf, ln
        w, Uv = np.linalg.eigh(S)
        keep = w > 1e-7 * max(w.max(), 1e-12)
        delta = t - mu
        proj = Uv[:, keep] @ (Uv[:, keep].T @ delta)
        if np.linalg.norm(delta - proj) > 1e-6 * (1 + np.linalg.norm(delta)):
            return -math.inf, ln
        k_ = int(keep.sum())
        q = float(np.sum((Uv[:, keep].T @ delta) ** 2 / w[keep]))
        lp = log10(covol) - 0.5 * k_ * log10(2 * math.pi) - 0.5 * float(np.sum(np.log10(w[keep]))) - q / 2 / math.log(10)
        return lp, ln


def triage(bs, K=400, seed=1, maxpool=6e10, maxtup=200, qcoords=(0, 1)):
    rng = np.random.default_rng(seed)
    gc = GComp(bs, qcoords)
    out = []
    for name, rows in bs.patterns():
        g0 = bs.gauss(rows)
        fs = gc.feasible(rows)
        if not fs:
            out.append((rows, g0, None, None)); continue
        tot, bt = 0.0, None
        for members, mult in fs[:maxtup]:
            tup = tuple(mm[0] for mm in members)
            pools = [gc.compressions(k, r)[x] for k, r, x in zip(KINDS, rows, tup)]
            if max(pools) > maxpool:
                continue
            lp, ln = gc.estimate(tup, K, rng)
            if lp == -math.inf:
                continue
            tot += mult * 10 ** (lp + ln)
            if bt is None or lp + ln > bt[0]:
                bt = (lp + ln, tup, pools)
        if len(fs) > maxtup and tot > 0:
            tot *= len(fs) / maxtup
        # solutions up to negating the blocks of row sum 0 (the convention of bsz.gauss); the tuple estimates count
        # every sign choice, and negation maps a compression tuple to another one
        e = log10(tot) - sum(1 for r in rows if r == 0) * log10(2) if tot > 0 else -math.inf
        out.append((rows, g0, (e, len(fs)), bt))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=["spin"]); ap.add_argument("N", type=int, nargs="+")
    ap.add_argument("--maxo", type=int, default=64); ap.add_argument("--mino", type=int, default=4)
    ap.add_argument("--K", type=int, default=400); ap.add_argument("--maxpool", type=float, default=6e10)
    a = ap.parse_args()
    for N in a.N:
        if (N - 1) % 4:
            continue
        for mods, mats, tag, O in BG.spin_structures(N, a.maxo, a.mino):
            bs = BG.BStruct(N, mods, mats, tag)
            head = f"N={N} {tag} O_A={bs.stA.O} O_L={bs.stL.O} cls={bs.C}"
            try:
                res = triage(bs, a.K, maxpool=a.maxpool)
            except Exception as ex:
                print(f"{head}: {type(ex).__name__} {ex}", flush=True); continue
            for rows, g0, fe, bt in res:
                if fe is None:
                    print(f"{head} rows {rows}: dead on Z_2^2 (plain Gaussian {g0:.2f})", flush=True)
                elif bt:
                    print(f"{head} rows {rows}: {fe[1]} feasible groups, conditioned logE {fe[0]:.2f} (plain {g0:.2f}); "
                          f"best tuple {bt[1]} logE {bt[0]:.2f} pools {['%.2g' % p for p in bt[2]]}", flush=True)
                else:
                    print(f"{head} rows {rows}: {fe[1]} feasible groups, no estimable tuple (plain {g0:.2f})", flush=True)


if __name__ == "__main__":
    main()
