"""U/W split for the bordered cyclotomic product Z_6 x GF(r) with Q = Z_3 and T = all characters (NOTES section 18).

G = Z_6 x GF(r), K = {a = 0 mod 3} x GF(r), targets sum_i PAF_i(g) = -12 on K \\ 0 and 0 off K; the sequences are
invariant under the e-th powers C_0 of GF(r)^* (F = (r - 1)/e odd). Stage 1 (bgs.stage1) fixes the compressions
y_k(a) = f_k(a) + F s_k(a) on Z_6.

**The split.** Pair the positions a and a + 3 (alpha = a mod 3): (a0, a1) = (0, 3), (4, 1), (2, 5). With
U = (x(a0) + x(a1))/2 and W = (x(a0) - x(a1))/2, both 0/+-1 functions on Z_3 x GF(r) with complementary supports,

    sum_k PAF_{U_k}(beta, d) = -6 [beta = 0]   and   sum_k PAF_{W_k}(beta, d) = 0     for (beta, d) != 0,

and conversely any such pair (with the supports complementary and the stage-1 compressions) gives a solution. The
U-condition has 14 classes (12 independent): the support size, (0, C_t) for t < e/2, (1, 0) and (1, C_t) for
t < e. Each block is three slices (alpha = 0, 1, 2); a slice is the U-value at (alpha, 0) (fixed by stage 1) and
a pattern in {-1, 0, 1}^e on the classes, and its W-completions are the sign choices on the zero cells. The search
finds U-solutions with a four-list k-tree (search/kt4.c) and completes each by an exact W-join.

  python3 search/bgs_uw.py search R E [--L N] [--sigmas N] [--seed S] [--threads T]   (Z_6 x GF(R), class UW)
  python3 search/bgs_uw.py scan N [N ...] [--maxspikes K]                                (structures for the split)
"""
import argparse, itertools, json, math, os, struct, subprocess, sys, tempfile, time
import numpy as np
from pathlib import Path as _P
sys.path[:0] = [str(_P(__file__).resolve().parent), str(_P(__file__).resolve().parent.parent / 'verify')]
import workdir as W
import verify
import bgs

PAIRS = [(0, 3), (4, 1), (2, 5)]


class UW:
    def __init__(self, r, e):
        self.r, self.e = r, e
        self.S1 = bgs.Structure(6, r, e, 3, (1, 3))
        self.n = self.S1.n; self.F = (r - 1) // e
        quads = bgs.stage1(self.S1)
        assert quads, "stage 1 is empty"
        self.y = quads[0]
        g = bgs.primitive_root(r); ind = np.zeros(r, dtype=np.int64); x = 1
        for k in range(r - 1): ind[x] = k; x = x * g % r
        self.cls = np.full(r, -1, dtype=np.int64); self.cls[1:] = ind[1:] % e
        self.col = np.where(np.arange(r) == 0, 0, 1 + self.cls)          # cell of y in GF(r): 0 point, 1 + j class
        self.creps = [0] + [int(np.flatnonzero(self.cls == c)[0]) for c in range(e)]   # d = 0, x_{C_0} .. x_{C_{e-1}}
        # slices: per block and alpha, U point value, W point value, the distinct U-class patterns and their W lists
        self.slices = []
        for k in range(4):
            row = []
            for (a0, a1) in PAIRS:
                f0, s0 = self.fs(self.y[k][a0]); f1, s1 = self.fs(self.y[k][a1])
                pats = {}
                for m0 in range(1 << e):
                    x0 = [1 - 2 * ((m0 >> j) & 1) for j in range(e)]
                    if sum(x0) != s0: continue
                    for m1 in range(1 << e):
                        x1 = [1 - 2 * ((m1 >> j) & 1) for j in range(e)]
                        if sum(x1) != s1: continue
                        u = tuple((p + q) // 2 for p, q in zip(x0, x1)); w = tuple((p - q) // 2 for p, q in zip(x0, x1))
                        pats.setdefault(u, []).append(w)
                keys = sorted(pats)
                row.append(dict(up=(f0 + f1) // 2, wp=(f0 - f1) // 2, U=np.array(keys, dtype=np.int8),
                                Wl=[np.array(pats[u], dtype=np.int8) for u in keys]))
            self.slices.append(row)
        # global table of U-slice patterns (point value, classes) -> id, and correlation table
        allp = {}
        for k in range(4):
            for a in range(3):
                sl = self.slices[k][a]
                sl['ids'] = np.array([allp.setdefault((sl['up'],) + tuple(int(t) for t in u), len(allp)) for u in sl['U']], dtype=np.int32)
        self.upats = np.array(list(allp.keys()), dtype=np.int64)        # [P, 1 + e]
        self.CC = self.corr_table(self.upats)
        self.T = np.array([6 * r - 3] + [-6] * (e // 2) + [0] * (e + 1), dtype=np.int64)   # U-target
        self.TW = np.array([self.n] + [0] * (e // 2) + [0] * (e + 1), dtype=np.int64)       # W-target

    def fs(self, yv):
        s = int(round(yv / self.F)); f = yv - self.F * s
        assert f in (1, -1); return f, s

    def corr_table(self, pats):
        """CC[p, q, c] = sum_y P_p(y) P_q(y + d_c) over GF(r), d_c in (0, x_{C_0}, ..)."""
        V = pats[:, self.col].astype(float)                     # [P, r]
        Fv = np.fft.rfft(V, axis=1)
        out = np.empty((len(pats), len(pats), len(self.creps)), dtype=np.int16)
        for i in range(len(pats)):
            R = np.fft.irfft(np.conj(Fv[i])[None, :] * Fv, n=self.r, axis=1)[:, self.creps]   # sum_y P_i(y) P_q(y + d)
            out[i] = np.rint(R).astype(np.int16)
        return out

    def keys(self, ids):
        """ids [n, 3] global U-pattern ids of the slices alpha = 0, 1, 2 -> U-keys [n, 14] (int16)."""
        e = self.e; CC = self.CC
        intra = sum(CC[ids[:, a], ids[:, a]][:, :1 + e // 2].astype(np.int32) for a in range(3))
        cross = sum(CC[ids[:, a], ids[:, (a + 1) % 3]].astype(np.int32) for a in range(3))
        return np.concatenate([intra, cross], axis=1).astype(np.int16)

    def sample(self, k, L, rng):
        """L random U-vectors of block k: local pattern indices [L, 3] and global ids [L, 3]."""
        loc = np.stack([rng.integers(0, len(self.slices[k][a]['U']), L) for a in range(3)], axis=1)
        gid = np.stack([self.slices[k][a]['ids'][loc[:, a]] for a in range(3)], axis=1)
        return loc, gid

    # ---- full sequences ------------------------------------------------------------------------------------
    def block_array(self, k, uloc, wsel):
        """6 x (e + 1) array of x_k from local U-pattern indices (3) and W-pattern indices (3)."""
        P = np.zeros((6, self.e + 1), dtype=np.int64)
        for a, (a0, a1) in enumerate(PAIRS):
            sl = self.slices[k][a]
            Uv = np.r_[sl['up'], sl['U'][uloc[a]]].astype(np.int64)
            Wv = np.r_[sl['wp'], sl['Wl'][uloc[a]][wsel[a]]].astype(np.int64)
            P[a0] = Uv + Wv; P[a1] = Uv - Wv
        assert np.all(np.abs(P) == 1)
        return P

    def sequences(self, Ps):
        """+-1 sequences on Z_{6r} (CRT: t -> (t mod 6, t mod r)) from the 6 x (e+1) arrays."""
        v = 6 * self.r; t = np.arange(v)
        return [P[t % 6, self.col[t % self.r]] for P in Ps]

    def wjoin(self, ulocs):
        """All W-completions of the U-solution given by local indices ulocs[k][a]; returns list of wsel tuples."""
        e = self.e; Ls = []; Ks = []
        for k in range(4):
            lists = [self.slices[k][a]['Wl'][ulocs[k][a]] for a in range(3)]
            combos = list(itertools.product(*[range(len(l)) for l in lists]))
            pats = np.zeros((len(combos), 3, e + 1), dtype=np.int64)
            for i, cb in enumerate(combos):
                for a in range(3):
                    pats[i, a, 0] = self.slices[k][a]['wp']; pats[i, a, 1:] = lists[a][cb[a]]
            Ls.append(combos); Ks.append(self.wkeys(pats))
        A, B, C, D = Ks
        AB = (A[:, None, :] + B[None, :, :]).reshape(-1, A.shape[1])
        need = self.TW[None, :] - (C[:, None, :] + D[None, :, :]).reshape(-1, C.shape[1])
        index = {}
        for i, row in enumerate(map(bytes, AB)): index.setdefault(row, []).append(i)
        out = []
        for j, row in enumerate(map(bytes, need)):
            for i in index.get(row, ()):
                out.append((Ls[0][i // len(B)], Ls[1][i % len(B)], Ls[2][j // len(D)], Ls[3][j % len(D)]))
        return out

    def cell_matrices(self):
        """M[c][i, j] = #{y in cell i : y + d_c in cell j} (cells: 0 = {0}, 1 + j = C_j; d_c = 0, x_{C_0}, ..)."""
        if getattr(self, '_M', None) is None:
            r = self.r; cells = self.col
            M = np.zeros((len(self.creps), self.e + 1, self.e + 1), dtype=np.int64)
            for c, d in enumerate(self.creps):
                np.add.at(M[c], (cells, cells[(np.arange(r) + d) % r]), 1)
            self._M = M
        return self._M

    def wjoin_fast(self, ulocs, K=10 ** 6):
        """Same as wjoin, vectorised: block W-keys from per-slice intra and cross terms, then a hashed 4-list join."""
        e = self.e; M = self.cell_matrices(); h = 1 + e // 2
        blocks = []
        for k in range(4):
            Wm = []
            for a in range(3):
                sl = self.slices[k][a]; Wc = sl['Wl'][ulocs[k][a]].astype(np.int64)
                Wm.append(np.concatenate([np.full((len(Wc), 1), sl['wp'], dtype=np.int64), Wc], axis=1))   # [n_a, e+1]
            intra = [np.einsum('ni,cij,nj->nc', Wm[a], M[:h], Wm[a]) for a in range(3)]
            cross = [np.einsum('ni,cij,mj->nmc', Wm[a], M, Wm[(a + 1) % 3]) for a in range(3)]
            n0, n1, n2 = (len(x) for x in Wm)
            I = (intra[0][:, None, None, :] + intra[1][None, :, None, :] + intra[2][None, None, :, :])
            C = (cross[0][:, :, None, :] + cross[1][None, :, :, :] + cross[2].transpose(1, 0, 2)[:, None, :, :])
            Kb = np.concatenate([np.broadcast_to(I, (n0, n1, n2, h)), C], axis=3).reshape(-1, h + e + 1)
            idx = np.stack(np.unravel_index(np.arange(n0 * n1 * n2), (n0, n1, n2)), axis=1)
            blocks.append((Kb, idx))
        R = np.random.default_rng(7).integers(1, 1 << 40, size=h + e + 1)
        hA, hB, hC, hD = ((b[0] * R).sum(1) for b in blocks)
        hT = int((self.TW * R).sum())
        AB = (hA[:, None] + hB[None, :]).ravel(); o = np.argsort(AB, kind='stable'); ABs = AB[o]
        need = (hT - (hC[:, None] + hD[None, :])).ravel()
        lo = np.searchsorted(ABs, need, 'left'); hi = np.searchsorted(ABs, need, 'right')
        out = []
        nB, nD = len(hB), len(hD)
        for j in np.flatnonzero(hi > lo):
            for p in range(lo[j], hi[j]):
                i = o[p]; a, b = divmod(int(i), nB); c, d = divmod(int(j), nD)
                if np.array_equal(blocks[0][0][a] + blocks[1][0][b] + blocks[2][0][c] + blocks[3][0][d], self.TW):
                    out.append(tuple(tuple(int(t) for t in blocks[q][1][x]) for q, x in enumerate((a, b, c, d))))
                    if len(out) >= K: return out
        return out

    def wkeys(self, pats):
        """pats [n, 3, e+1] W-slices -> 14 class values (same layout as keys), int64."""
        e = self.e; r = self.r
        V = pats[:, :, self.col].astype(float)                 # [n, 3, r]
        Fv = np.fft.rfft(V, axis=2)
        def corr(a, b):
            return np.rint(np.fft.irfft(np.conj(Fv[:, a]) * Fv[:, b], n=r, axis=1)[:, self.creps]).astype(np.int64)
        intra = sum(corr(a, a)[:, :1 + e // 2] for a in range(3))
        cross = sum(corr(a, (a + 1) % 3) for a in range(3))
        return np.concatenate([intra, cross], axis=1)

    def check(self, Ps):
        """Bordered targets on Z_{6r}, then the border and H H^T = 4n I. Returns (ok, H)."""
        import spence_lib as L, bordered as BD
        r, n = self.r, self.n; v = 6 * r
        xs = self.sequences(Ps)
        tot = verify.paf_sum(xs)
        want = np.array([4 * v] + [-12 if (d % 6) % 3 == 0 else 0 for d in range(1, v)])
        if not np.array_equal(tot, want): return False, None
        G = L.AbGroup([3, 2 * r]); co = G.coords
        tt = np.array([next(s for s in range(int(c1), v, 2 * r) if s % 3 == c0) for c0, c1 in co])
        seqs = [xx[tt] for xx in xs]
        yc = [[int(s[co[:, 0] == a].sum()) for a in range(3)] for s in seqs]
        b = BD.border(yc)
        if b is None: return False, None
        H = BD.assemble(G, seqs, *b)
        return bool(np.array_equal(H @ H.T, 4 * n * np.eye(4 * n, dtype=np.int64))), H



class UWG:
    """General even m = 2m': G = Z_m x GF(r), K of index q in Z_m (x GF(r)), T given by the divisor set D, sequences
    invariant under the e-th powers of GF(r)* and under K_m = <kgen> on Z_m (bgs.Structure). Slices alpha in [0, m')
    pair the positions alpha and alpha + m'; U = (x(alpha) + x(alpha + m'))/2, W = (x(alpha) - x(alpha + m'))/2. With c
    the carry of alpha + beta past m',

        PAF_x(beta + j m', d) = 2 (PAF_U(beta, d) + (-1)^j PAF_W(beta, d)),
        PAF_U(beta, d) = sum_alpha Corr(U_alpha, U_{alpha+beta mod m'})(d),
        PAF_W(beta, d) = sum_alpha (-1)^c Corr(W_alpha, W_{alpha+beta mod m'})(d)   (negacyclic),

    so the bordered condition is sum_k PAF_U = (t(beta) + t(beta + m'))/4 and sum_k PAF_W = (t(beta) - t(beta + m'))/4,
    t the target of sum_k PAF_x. K_m maps pairs to pairs (k is odd); a slice taken to the second position of its image
    pair carries -W there. Free slices are the K_m-orbit representatives; an orbit that maps a slice to itself with a
    flip forces W = 0 on it. Keys: all (beta, d), beta in [0, m'), d in (0, x_{C_0}, ..) (redundant classes kept)."""

    def __init__(self, m, r, e, q, D, config=0, kgen=1):
        assert m % 2 == 0
        self.m, self.mh, self.r, self.e = m, m // 2, r, e
        self.S1 = bgs.Structure(m, r, e, q, D, kgen)
        self.n = self.S1.n; self.F = (r - 1) // e
        quads = bgs.stage1(self.S1)
        assert quads, "stage 1 is empty"
        self.quads = quads; self.y = quads[config]
        g = bgs.primitive_root(r); ind = np.zeros(r, dtype=np.int64); x = 1
        for k in range(r - 1): ind[x] = k; x = x * g % r
        self.cls = np.full(r, -1, dtype=np.int64); self.cls[1:] = ind[1:] % e
        self.col = np.where(np.arange(r) == 0, 0, 1 + self.cls)
        self.creps = [0] + [int(np.flatnonzero(self.cls == c)[0]) for c in range(e)]
        mh = self.mh
        # slice orbits under K_m: rep[alpha] = representative, sgn[alpha] = +-1 (W sign relative to the representative)
        self.rep = -np.ones(mh, dtype=np.int64); self.sgn = np.zeros(mh, dtype=np.int64); self.selfflip = set()
        Km = self.S1.K
        for a in range(mh):
            if self.rep[a] >= 0: continue
            for kk in Km:
                b = kk * a % m; a2, fl = b % mh, (-1 if b >= mh else 1)
                if self.rep[a2] < 0: self.rep[a2] = a; self.sgn[a2] = fl
                elif self.rep[a2] == a and self.sgn[a2] != fl: self.selfflip.add(a)
        self.reps = sorted(set(self.rep.tolist()))
        self.osize = {a: int((self.rep == a).sum()) for a in self.reps}
        self.slices = []
        for k in range(4):
            row = [None] * mh
            for a in self.reps:
                f0, s0 = self.fs(self.y[k][a]); f1, s1 = self.fs(self.y[k][a + mh])
                pats = {}
                for m0 in range(1 << e):
                    x0 = [1 - 2 * ((m0 >> j) & 1) for j in range(e)]
                    if sum(x0) != s0: continue
                    for m1 in range(1 << e):
                        x1 = [1 - 2 * ((m1 >> j) & 1) for j in range(e)]
                        if sum(x1) != s1: continue
                        u = tuple((p + q_) // 2 for p, q_ in zip(x0, x1)); w = tuple((p - q_) // 2 for p, q_ in zip(x0, x1))
                        if a in self.selfflip and (any(w) or f0 != f1): continue
                        pats.setdefault(u, []).append(w)
                keys = sorted(pats)
                row[a] = dict(up=(f0 + f1) // 2, wp=(f0 - f1) // 2, U=np.array(keys, dtype=np.int8).reshape(len(keys), e),
                              Wl=[np.array(pats[u], dtype=np.int8) for u in keys])
            self.slices.append(row)
        allp = {}
        for k in range(4):
            for a in self.reps:
                sl = self.slices[k][a]
                sl['ids'] = np.array([allp.setdefault((sl['up'],) + tuple(int(t) for t in u), len(allp)) for u in sl['U']], dtype=np.int32)
        self.upats = np.array(list(allp.keys()), dtype=np.int64)
        self.CC = UW.corr_table(self, self.upats)
        tg = lambda b, d: 4 * m * r if (b % m == 0 and d == 0) else -int(self.S1.Gamma[(b % m) % q])
        self.T = np.array([(tg(b, d) + tg(b + mh, d)) // 4 for b in range(mh) for d in range(1 + e)], dtype=np.int64)
        self.TW = np.array([(tg(b, d) - tg(b + mh, d)) // 4 for b in range(mh) for d in range(1 + e)], dtype=np.int64)
        assert all((tg(b, d) + tg(b + mh, d)) % 4 == 0 and (tg(b, d) - tg(b + mh, d)) % 4 == 0 for b in range(mh) for d in range(1 + e))
        npts = sum(abs(self.slices[k][self.rep[a]]['up']) for k in range(4) for a in range(mh))
        assert (self.T[0] - npts) % self.F == 0
        self.ncells = (self.T[0] - npts) // self.F          # nonzero U class cells over all blocks and slices

    fs = UW.fs
    cell_matrices = UW.cell_matrices

    def expand(self, gid_rep):
        """[n, #reps] global ids at the representatives -> [n, m'] ids at all slices."""
        col = {a: i for i, a in enumerate(self.reps)}
        return gid_rep[:, [col[int(self.rep[a])] for a in range(self.mh)]]

    def keys(self, ids):
        """ids [n, m'] global U-pattern ids of all slices -> U-keys [n, m' (1 + e)] int16, layout (beta, d)."""
        mh, E = self.mh, 1 + self.e; CC = self.CC
        out = np.empty((len(ids), mh * E), dtype=np.int16)
        for b in range(mh):
            acc = sum(CC[ids[:, a], ids[:, (a + b) % mh]].astype(np.int32) for a in range(mh))
            out[:, b * E:(b + 1) * E] = acc
        return out

    def sample(self, k, total, L, rng):
        """L U-vectors of block k with exactly `total` nonzero class cells (counted over all slices), uniform among
        them. Returns local pattern indices at the representatives, global ids at all slices, and the count."""
        e = self.e; R = self.reps
        nz = [(self.slices[k][a]['U'] != 0).sum(1) for a in R]
        cnt = [np.bincount(z, minlength=e + 1).astype(float) for z in nz]
        w8 = [self.osize[a] for a in R]
        dp = [np.zeros(total + 1) for _ in range(len(R) + 1)]; dp[len(R)][0] = 1.0
        for i in range(len(R) - 1, -1, -1):
            for tot in range(total + 1):
                dp[i][tot] = sum(cnt[i][j] * dp[i + 1][tot - j * w8[i]] for j in range(e + 1) if tot - j * w8[i] >= 0)
        N = dp[0][total]
        if N == 0: raise ValueError(f"block {k}: no U-vector with {total} nonzero class cells")
        byn = [[np.flatnonzero(nz[i] == j) for j in range(e + 1)] for i in range(len(R))]
        loc = np.empty((L, len(R)), dtype=np.int64); rem = np.full(L, total)
        for i in range(len(R)):
            w = np.array([[cnt[i][j] * dp[i + 1][tt - j * w8[i]] if tt - j * w8[i] >= 0 else 0.0 for j in range(e + 1)]
                          for tt in range(total + 1)])
            w = w / np.maximum(w.sum(1, keepdims=True), 1e-300)
            cw = np.cumsum(w, axis=1); u = rng.random(L)
            js = np.minimum((u[:, None] > cw[rem]).sum(1), e)
            for j in np.unique(js):
                sel = np.flatnonzero(js == j); pool = byn[i][j]
                loc[sel, i] = pool[rng.integers(0, len(pool), len(sel))]
            rem = rem - js * w8[i]
        assert np.all(rem == 0)
        gid = self.expand(np.stack([self.slices[k][a]['ids'][loc[:, i]] for i, a in enumerate(R)], axis=1))
        return loc, gid, N

    def slice_W(self, k, uloc):
        """Per representative: W options [n_i, e + 1] (point value first)."""
        out = []
        for i, a in enumerate(self.reps):
            sl = self.slices[k][a]; Wc = sl['Wl'][uloc[i]].astype(np.int64)
            out.append(np.concatenate([np.full((len(Wc), 1), sl['wp'], dtype=np.int64), Wc], axis=1))
        return out

    def wjoin_fast(self, ulocs, K=10 ** 6):
        e = self.e; M = self.cell_matrices(); mh = self.mh; E = 1 + e; R = self.reps
        col = {a: i for i, a in enumerate(R)}
        blocks = []
        for k in range(4):
            Wr = self.slice_W(k, ulocs[k]); shape = tuple(len(x) for x in Wr)
            Kb = np.zeros(shape + (mh * E,), dtype=np.int64)
            for b in range(mh):
                for a in range(mh):
                    a2 = (a + b) % mh; sg = (-1 if a + b >= mh else 1) * int(self.sgn[a]) * int(self.sgn[a2])
                    i1, i2 = col[int(self.rep[a])], col[int(self.rep[a2])]
                    if i1 == i2:
                        t = np.einsum('ni,cij,nj->nc', Wr[i1], M, Wr[i1]) * sg
                        sh = [1] * len(R); sh[i1] = shape[i1]
                        Kb[..., b * E:(b + 1) * E] += t.reshape(sh + [E])
                    else:
                        t = np.einsum('ni,cij,mj->nmc', Wr[i1], M, Wr[i2]) * sg
                        sh = [1] * len(R); sh[i1] = shape[i1]; sh[i2] = shape[i2]
                        if i1 < i2: Kb[..., b * E:(b + 1) * E] += t.reshape(sh + [E])
                        else: Kb[..., b * E:(b + 1) * E] += t.transpose(1, 0, 2).reshape(sh + [E])
            n = int(np.prod(shape))
            idx = np.stack(np.unravel_index(np.arange(n), shape), axis=1)
            blocks.append((Kb.reshape(n, mh * E), idx))
        Rr = np.random.default_rng(7).integers(1, 1 << 40, size=mh * E)
        hA, hB, hC, hD = ((b[0] * Rr).sum(1) for b in blocks)
        hT = int((self.TW * Rr).sum())
        AB = (hA[:, None] + hB[None, :]).ravel(); o = np.argsort(AB, kind='stable'); ABs = AB[o]
        need = (hT - (hC[:, None] + hD[None, :])).ravel()
        lo = np.searchsorted(ABs, need, 'left'); hi = np.searchsorted(ABs, need, 'right')
        out = []; nB, nD = len(hB), len(hD)
        for j in np.flatnonzero(hi > lo):
            for p in range(lo[j], hi[j]):
                i = o[p]; a, b = divmod(int(i), nB); c, d = divmod(int(j), nD)
                if np.array_equal(blocks[0][0][a] + blocks[1][0][b] + blocks[2][0][c] + blocks[3][0][d], self.TW):
                    out.append(tuple(tuple(int(t) for t in blocks[q_][1][x]) for q_, x in enumerate((a, b, c, d))))
                    if len(out) >= K: return out
        return out

    def block_array(self, k, uloc, wsel):
        mh = self.mh; P = np.zeros((self.m, self.e + 1), dtype=np.int64)
        col = {a: i for i, a in enumerate(self.reps)}
        for a in range(mh):
            i = col[int(self.rep[a])]; sl = self.slices[k][int(self.rep[a])]
            Uv = np.r_[sl['up'], sl['U'][uloc[i]]].astype(np.int64)
            Wv = int(self.sgn[a]) * np.r_[sl['wp'], sl['Wl'][uloc[i]][wsel[i]]].astype(np.int64)
            P[a] = Uv + Wv; P[a + mh] = Uv - Wv
        assert np.all(np.abs(P) == 1)
        return P

    def sequences(self, Ps):
        v = self.m * self.r; t = np.arange(v)
        return [P[t % self.m, self.col[t % self.r]] for P in Ps]

    def check_targets(self, Ps):
        xs = self.sequences(Ps); v = self.m * self.r
        tot = verify.paf_sum(xs)
        want = np.array([4 * v] + [-int(self.S1.Gamma[(d % self.m) % self.S1.q]) for d in range(1, v)])
        return np.array_equal(tot, want), xs

def nnz_profiles(S, k, total):
    """Per-slice nonzero-class-cell counts (n0, n1, n2) with sum `total` and their numbers of U-patterns."""
    cnt = [np.bincount((S.slices[k][a]['U'] != 0).sum(1), minlength=S.e + 1) for a in range(3)]
    out = []
    for n0 in range(S.e + 1):
        for n1 in range(S.e + 1):
            n2 = total - n0 - n1
            if 0 <= n2 <= S.e:
                c = int(cnt[0][n0]) * int(cnt[1][n1]) * int(cnt[2][n2])
                if c: out.append(((n0, n1, n2), c))
    return out


def sample_nnz(S, k, total, L, rng):
    """L U-vectors of block k with exactly `total` nonzero class cells, uniform among such vectors."""
    prof = nnz_profiles(S, k, total)
    w = np.array([c for _, c in prof], dtype=float); w /= w.sum()
    which = rng.choice(len(prof), size=L, p=w)
    loc = np.empty((L, 3), dtype=np.int64)
    nz = [(S.slices[k][a]['U'] != 0).sum(1) for a in range(3)]
    byn = [[np.flatnonzero(nz[a] == j) for j in range(S.e + 1)] for a in range(3)]
    for pi, (p, _) in enumerate(prof):
        sel = np.flatnonzero(which == pi)
        for a in range(3):
            pool = byn[a][p[a]]
            loc[sel, a] = pool[rng.integers(0, len(pool), len(sel))]
    gid = np.stack([S.slices[k][a]['ids'][loc[:, a]] for a in range(3)], axis=1)
    return loc, gid, float(sum(c for _, c in prof))


def gauss_estimate(S, totals, samples=200000, seed=1):
    """Full-covariance Gaussian estimate of the number of U-solutions whose blocks have the given numbers of nonzero
    class cells (the support coordinate is then fixed), and of the mean number of W-completions."""
    rng = np.random.default_rng(seed); Ks = []; logN = 0.0
    for k in range(4):
        loc, gid, N = sample_nnz(S, k, totals[k], samples, rng)
        Ks.append(S.keys(gid).astype(float)[:, 1:]); logN += math.log10(N)
    mu = sum(K.mean(0) for K in Ks); Sig = sum(np.cov(K.T) for K in Ks)
    T = S.T[1:].astype(float); diff = T - mu
    ev, U = np.linalg.eigh(Sig); keep = ev > 1e-6 * ev.max()
    z = U[:, keep].T @ diff
    dens = sum(-0.5 * math.log10(2 * math.pi * x) for x in ev[keep]) - 0.5 * float((z * z / ev[keep]).sum()) / math.log(10)
    off = [float(U[:, j] @ diff) for j in np.flatnonzero(~keep)]
    return logN + dens, int(keep.sum()), off


def run_kt4(Ks, T, S_idx, nsig, K=10 ** 6, threads=None, log=print, maxpairs=3 * 10 ** 8):
    fd, path = tempfile.mkstemp(suffix='.kt4', dir=W.tmpdir()); os.close(fd)
    try:
        with open(path, 'wb') as f:
            f.write(struct.pack('5i', *[len(k) for k in Ks], Ks[0].shape[1]))
            f.write(np.asarray(T, dtype=np.int16).tobytes())
            for k in Ks: f.write(np.ascontiguousarray(k, dtype=np.int16).tobytes())
        env = dict(os.environ)
        if threads: env['OMP_NUM_THREADS'] = str(threads)
        p = subprocess.Popen([W.bin('kt4'), path, ','.join(map(str, S_idx)), str(nsig), str(K), str(maxpairs)], stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE, text=True, env=env)
        sols = []
        for line in p.stdout:
            if line.startswith('SOL'): sols.append(tuple(map(int, line.split()[1:])))
            elif line.startswith('TOTAL'): log(f"  kt4: {line.strip()}")
        err = p.stderr.read(); p.wait()
        for line in err.splitlines()[-max(3, nsig + 2):]: log('  ' + line)
        return sols
    finally:
        os.remove(path)


def save_hit(S, Ps, H, tag):
    import hashlib
    xs = S.sequences(Ps)
    rec = dict(n=S.n, order=4 * S.n, construction='bordered-gs-uw', r=S.r, e=S.e, tag=tag,
               seqs=[''.join('+' if z > 0 else '-' for z in x) for x in xs],
               sha256=hashlib.sha256(('\n'.join(','.join(str(int(t)) for t in row) for row in H) + '\n').encode()).hexdigest())
    with open(W.hits('bgs_uw_hits.jsonl'), 'a') as fh: fh.write(json.dumps(rec) + '\n')
    return rec


_CTX = None


def _wjoin_chunk(chunk):
    S, locs = _CTX
    out = []
    for (a, b, c, d) in chunk:
        ul = [locs[0][a].astype(np.int64), locs[1][b].astype(np.int64), locs[2][c].astype(np.int64), locs[3][d].astype(np.int64)]
        for ws in S.wjoin_fast(ul, K=2):
            out.append((ul, ws))
    return len(chunk), out


def search(r, e, L, nsig, seed, totals=(12, 12, 12, 12), S_idx=(6, 7, 8, 9, 10), maxu=10 ** 6, threads=None, log=print,
           workers=4, stop_on_hit=True):
    t0 = time.time()
    S = UW(r, e)
    log(f"n = {S.n} = 6 * {r} + 3, e = {e}, F = {S.F}; stage-1 compressions {S.y}; U-patterns {len(S.upats)} [{time.time() - t0:.0f}s]")
    rng = np.random.default_rng(seed)
    locs, Ks = [], []
    for k in range(4):
        loc, gid, N = sample_nnz(S, k, totals[k], L, rng)
        Kk = np.concatenate([S.keys(gid[c:c + (1 << 20)]) for c in range(0, L, 1 << 20)])
        del gid
        locs.append(loc.astype(np.int16)); Ks.append(Kk); del loc
        log(f"  block {k}: {L} of {N:.3g} U-vectors with {totals[k]} nonzero class cells [{time.time() - t0:.0f}s]")
    sols = run_kt4(Ks, S.T, S_idx, nsig, maxu, threads, log)
    log(f"  {len(sols)} U-solutions [{time.time() - t0:.0f}s]")
    global _CTX
    _CTX = (S, locs)
    import multiprocessing as mp
    full = 0; hits = []; done = 0
    chunks = [sols[i:i + 500] for i in range(0, len(sols), 500)]
    with mp.get_context('fork').Pool(workers) as pool:
        for res in pool.imap_unordered(_wjoin_chunk, chunks):
            done += res[0]; full += len(res[1])
            for (ul, ws) in res[1]:
                if len(hits) >= 3: break
                Ps = [S.block_array(k, ul[k], ws[k]) for k in range(4)]
                ok, H = S.check(Ps)
                log(f"  full solution: H H^T = {4 * S.n} I: {ok} [{time.time() - t0:.0f}s]")
                if ok:
                    hits.append(save_hit(S, Ps, H, f'seed={seed} L={L} nsig={nsig}'))
                    cert = certificate(S, Ps, H)
                    d = W.WORK / 'bgs_uw_certs'; d.mkdir(exist_ok=True)
                    fn = d / f"{cert['order']}_{cert['sha256'][:12]}.json"
                    fn.write_text(json.dumps(cert, indent=1) + '\n')
                    log(f"  certificate written to {fn}")
            if done % 20000 < 500: log(f"  W-joins: {done} of {len(sols)} U-solutions, {full} completions [{time.time() - t0:.0f}s]")
            if hits and stop_on_hit: pool.terminate(); break
    log(f"  W-completions: {full} from {len(sols)} U-solutions [{time.time() - t0:.0f}s]")
    return hits



def certificate(S, Ps, H, nsols_note=''):
    """Certificate in the format of verify.build_bordered_gs: G = Z_3 x Z_{2r}, K generated by (0, 1)."""
    import hashlib, bordered as BD, spence_lib as L
    r = S.r; v = 6 * r
    xs = S.sequences(Ps)
    G = L.AbGroup([3, 2 * r]); co = G.coords
    tt = np.array([next(s for s in range(int(c1), v, 2 * r) if s % 3 == c0) for c0, c1 in co])
    seqs = [xx[tt] for xx in xs]
    yc = [[int(s[co[:, 0] == a].sum()) for a in range(3)] for s in seqs]
    X, Rh, Ch = BD.border(yc)
    Hc = BD.assemble(G, seqs, X, Rh, Ch)
    assert np.array_equal(Hc, H)
    sgn = lambda M: [''.join('+' if t > 0 else '-' for t in row) for row in M]
    g = bgs.primitive_root(r); h = pow(g, S.e, r)
    lam = next(l for l in range(1, 2 * r) if l % 2 == 1 and l % r == h)
    cert = dict(order=4 * S.n, construction='bordered-goethals-seidel', group=[3, 2 * r], subgroup_generators=[[0, 1]],
                symmetries=[[[1, 0], [0, lam]]], autocorrelation=dict(subgroup=-12, elsewhere=0),
                rowsums=[int(s.sum()) for s in seqs], sequences=sgn(seqs),
                border=dict(X=sgn(X), Rhat=sgn(Rh), Chat=sgn(Ch)))
    # verify.py develops A as a(h - g); bordered.assemble as a(g - h). Hash the verifier's matrix (checked below).
    Hv, _ = verify.build_bordered_gs(cert)
    assert np.array_equal(Hv @ Hv.T, 4 * S.n * np.eye(4 * S.n, dtype=np.int64)), "verifier's matrix is not Hadamard"
    sha = hashlib.sha256()
    for chunk in verify.csv_chunks(Hv): sha.update(chunk)
    cert['sha256'] = sha.hexdigest()
    return cert



def scan_even(ns, maxspikes=6, log=print):
    """Cyclic bordered cyclotomic products Z_m x GF(r) with m even (the U/W split applies), K_m trivial, at most
    `maxspikes` spikes: print every structure with a nonempty stage 1 (bgs.stage1 with reduced caps; structures with
    more than 22 position orbits or too many stage-1 options are reported as not covered)."""
    for n in ns:
        for (m, r, t, q, D) in bgs.candidate_structures(n, tmax=40, rmin=11, mmax=80):
            if m % 2: continue
            for e in range(2, r, 2):
                if (r - 1) % e or (r - 1) // e < 3: continue
                F = (r - 1) // e; budget = 4 * (m - t) * e; cmin = 4 * F - 4
                if budget < 0 or (budget and (budget < cmin or budget // cmin > maxspikes)): continue
                try:
                    S = bgs.Structure(m, r, e, q, D, 1)
                except AssertionError:
                    continue
                rec = dict(n=n, m=m, r=r, t=t, q=q, D=list(D), e=e)
                try:
                    quads = bgs.stage1(S, maxpairs=4e6, maxopts=4e5, maxspikes=maxspikes)
                except (RuntimeError, MemoryError) as ex:
                    log(json.dumps(dict(rec, uncovered=str(ex)[:60]))); continue
                if quads:
                    log(json.dumps(dict(rec, F=F, classes=len(quads), eqs=bgs.nequations(S))))

if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest='cmd', required=True)
    s = sub.add_parser('search'); s.add_argument('r', type=int); s.add_argument('e', type=int)
    s.add_argument('--L', type=int, default=10 ** 6); s.add_argument('--sigmas', type=int, default=1)
    s.add_argument('--seed', type=int, default=1); s.add_argument('--threads', type=int)
    s.add_argument('--totals', default=''); s.add_argument('--S', default='6,7,8,9,10'); s.add_argument('--maxu', type=int, default=10 ** 6)
    sc = sub.add_parser('scan'); sc.add_argument('n', type=int, nargs='+'); sc.add_argument('--maxspikes', type=int, default=6)
    a = ap.parse_args()
    log = lambda m: print(m, flush=True)
    if a.cmd == 'scan':
        scan_even(a.n, a.maxspikes, log)
    if a.cmd == 'search':
        tot = tuple(int(x) for x in a.totals.split(',')) if a.totals else (3 * a.e // 2,) * 4
        search(a.r, a.e, a.L, a.sigmas, a.seed, tot, tuple(int(x) for x in a.S.split(',')), a.maxu, a.threads, log)
