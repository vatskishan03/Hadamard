#!/usr/bin/env python3
"""Level 1 of the tower (the R level): exact compressions over E = G/Z11 = Z2^4 and their join, for one Q-tuple.

G -> E forgets the Z11 coordinate.  The compression r(e) = sum_z X(e, z) of a sequence is constant on the classes of E
(<H, -1>-orbits for A and B, H-orbits for C and D; six classes of sizes 1, 3, 3, 3, 3, 3 for each kind), a cell of the
sequence lies in a single class, so the number of sequences with a given r is a product over the classes.  `RLevel.vectors`
lists every r with the prescribed Q-compression, its multiplicity and its PAF classes on E; `join_r` joins the four blocks
exactly on the PAF equations pushed down to E (the target of the sum of PAFs summed over the Z11 fibres) and returns the
weighted number of sequence quadruples that survive; `rtuples.explicit_tuples` lists the surviving R-tuples themselves.
For the pattern (18, 0, 14, 14) there are 56,916 R-tuples (196 distinct weights).

    python3 rlevel.py i [18,0,14,14]        join_r of every normalised pattern (or of the given one): R-tuples and retention
"""
import sys, itertools, math, json, time
from collections import defaultdict
import numpy as np
from structure import Struct
from qlevel import tuples, rs

class RLevel:
    def __init__(self, S):
        self.S = S
        eco = S.coords[:, :-1]
        assert (S.M[-1, :-1] == 0).all() and (S.M[:-1, -1] == 0).all() and S.M[-1, -1] == 1
        self.eshape = tuple(int(x) for x in S.mods[:-1]); self.nE = int(np.prod(self.eshape))
        ecoords = np.array(list(itertools.product(*[range(x) for x in self.eshape])), dtype=np.int64)
        self.ecoords = ecoords
        eidx = lambda c: np.ravel_multi_index(tuple(c[..., i] % self.eshape[i] for i in range(len(self.eshape))), self.eshape)
        self.eidx = eidx
        self.piE = eidx(eco)                                           # G index -> E index
        ME = S.M[:-1, :-1]
        self.hE = eidx(ecoords @ ME.T); self.negE = eidx(-ecoords)
        self.add = eidx(ecoords[:, None, :] + ecoords[None, :, :])      # add[e,f]
        # classes on E for each kind
        self.cls = {}
        for kind in "AC":
            perms = [self.hE] + ([self.negE] if kind == "A" else [])
            self.cls[kind] = self._orbits(perms)
        self.paf_cls = self._orbits([self.hE, self.negE])                # classes of PAF(e'): <H,-1>-orbits on E
        self.npaf = int(self.paf_cls.max()) + 1
        self.paf_rep = np.array([np.flatnonzero(self.paf_cls == c)[0] for c in range(self.npaf)])
        # target pushed down
        tE = np.zeros(self.nE, dtype=np.int64); np.add.at(tE, self.piE, S.target)
        self.tE = tE[self.paf_rep]
        # Q map on E
        self.qE = (ecoords @ S.Lphi[:, :-1].T) % 2
        self.qE = self.qE[:, 0] * 2 + self.qE[:, 1]
        # per kind: classes, per-class distribution of r_c, contribution matrices
        self.info = {}
        for kind in "AC":
            orb = S.orbA if kind == "A" else S.orbC
            ncell = int(orb.max()) + 1
            cellclass = np.zeros(ncell, dtype=np.int64); cellsize = np.bincount(orb)
            for c in range(ncell):
                cl = set(self.cls[kind][self.piE[orb == c]].tolist()); assert len(cl) == 1
                cellclass[c] = cl.pop()
            ncl = int(self.cls[kind].max()) + 1
            clsize = np.bincount(self.cls[kind])
            dist = []
            for k in range(ncl):
                d = {0: 1}
                for c in np.flatnonzero(cellclass == k):
                    contrib = cellsize[c] // clsize[k]; assert cellsize[c] % clsize[k] == 0
                    nd = defaultdict(int)
                    for v, n in d.items():
                        nd[v + contrib] += n; nd[v - contrib] += n
                    d = nd
                dist.append(dict(d))
            reps = np.array([np.flatnonzero(self.cls[kind] == k)[0] for k in range(ncl)])
            # Q contribution of a class value: y(a) = sum_classes r_c * |{e in class c: q(e) = a}|
            qc = np.zeros((ncl, 4), dtype=np.int64)
            for k in range(ncl):
                for a in range(4): qc[k, a] = int(((self.cls[kind] == k) & (self.qE == a)).sum())
            self.info[kind] = dict(ncl=ncl, dist=dist, reps=reps, qc=qc)

    def _orbits(self, perms):
        n = self.nE; par = list(range(n))
        def find(x):
            while par[x] != x: par[x] = par[par[x]]; x = par[x]
            return x
        for pm in perms:
            for i in range(n):
                a, b = find(i), find(int(pm[i]))
                if a != b: par[a] = b
        lab = {}; out = np.zeros(n, dtype=np.int64)
        for i in range(n): out[i] = lab.setdefault(find(i), len(lab))
        return out

    def vectors(self, kind, y):
        """all class-value vectors r (n, ncl) with Q-compression y=(y0,y1) (y(a)=y1 for a != 0), multiplicities (float64), PAF class vectors (n, npaf)"""
        inf = self.info[kind]; ncl = inf["ncl"]
        vals = [np.array(sorted(d)) for d in inf["dist"]]; cnts = [np.array([inf["dist"][k][v] for v in vals[k]], dtype=np.float64) for k in range(ncl)]
        grids = np.meshgrid(*[np.arange(len(v)) for v in vals], indexing="ij")
        R = np.stack([vals[k][grids[k].ravel()] for k in range(ncl)], 1)
        W = np.ones(len(R))
        for k in range(ncl): W *= cnts[k][grids[k].ravel()]
        y0 = R @ inf["qc"][:, 0]; y1 = R @ inf["qc"][:, 1]; y2 = R @ inf["qc"][:, 2]; y3 = R @ inf["qc"][:, 3]
        keep = (y0 == y[0]) & (y1 == y[1]) & (y2 == y[1]) & (y3 == y[1])
        R, W = R[keep], W[keep]
        full = R[:, self.cls[kind]]                                      # (n, nE)
        F = np.fft.fftn(full.reshape((-1,) + self.eshape).astype(np.float64), axes=tuple(range(1, 1 + len(self.eshape))))
        P = np.fft.ifftn(np.abs(F) ** 2, axes=tuple(range(1, 1 + len(self.eshape)))).real.reshape(len(R), -1)
        P = np.rint(P[:, self.paf_rep]).astype(np.int64)
        return R, W, P

def group(P, W):
    """unique PAF vectors with total weight and count"""
    u, inv = np.unique(P, axis=0, return_inverse=True)
    inv = inv.ravel()
    return u, np.bincount(inv, weights=W), np.bincount(inv)

def join_r(RL, tup, verbose=True, want_tuples=False):
    kinds = ["A", "A", "C", "C"]
    L = [RL.vectors(k, y) for k, y in zip(kinds, tup)]
    G = [group(P, W) for (R, W, P) in L]
    if verbose: print("  vectors", [len(l[0]) for l in L], "distinct PAF", [len(g[0]) for g in G])
    UA, WA, NA = G[0]; UB, WB, NB = G[1]; UC, WC, NC = G[2]; UD, WD, NDd = G[3]
    rng = np.random.default_rng(7); key = rng.integers(1, 2 ** 20, size=UA.shape[1]).astype(np.int64)
    kC = UC @ key; kD = UD @ key
    tk = int(RL.tE @ key)
    # table of (C,D) pair sums: key -> (sum of weights, count of R-pairs)
    nC, nD = len(UC), len(UD)
    ks = (kC[:, None] + kD[None, :]).ravel(); ws = (WC[:, None] * WD[None, :]).ravel(); cs = (NC[:, None] * NDd[None, :]).ravel().astype(np.float64)
    o = np.argsort(ks); ks, ws, cs = ks[o], ws[o], cs[o]
    cw = np.concatenate([[0.0], np.cumsum(ws)]); cc = np.concatenate([[0.0], np.cumsum(cs)])
    totw = 0.0; totc = 0.0
    kA = UA @ key; kB = UB @ key
    for i in range(len(UA)):
        need = tk - (kA[i] + kB)
        lo = np.searchsorted(ks, need, "left"); hi = np.searchsorted(ks, need, "right")
        m = hi > lo
        if m.any():
            # the key is a random linear hash of the PAF vectors (a collision would add a spurious tuple; rtuples.explicit_tuples
            # does the same join with exact dictionary keys and returns the same count)
            totw += float((WA[i] * WB[m] * (cw[hi[m]] - cw[lo[m]])).sum()); totc += float((NA[i] * NB[m] * (cc[hi[m]] - cc[lo[m]])).sum())
    return totw, totc, [g[0].shape[0] for g in G], L

if __name__ == "__main__":
    nm = sys.argv[1]
    want = tuple(int(x) for x in sys.argv[2].split(",")) if len(sys.argv) > 2 else None
    S = Struct(nm); RL = RLevel(S)
    print("E shape", RL.eshape, "classes A", RL.info["A"]["ncl"], "C", RL.info["C"]["ncl"], "paf classes", RL.npaf, "target_E", RL.tE.tolist())
    sols, dA, dC, tg = tuples(S)
    seen = set()
    for t in sols:
        key = tuple(abs(rs(y)) for y in t)
        if key in seen: continue
        seen.add(key)
        if want is not None and key != want: continue
        t2 = tuple(y if rs(y) >= 0 else (-y[0], -y[1]) for y in t)
        if key[0] < key[1] or key[2] < key[3]: continue
        t0 = time.time()
        totw, totc, ng, L = join_r(RL, t2)
        pools = [ {"A": dA, "C": dC}[k][y] for k, y in zip("AACC", t2)]
        lp = sum(math.log10(p) for p in pools)
        print(nm, key, "R-tuples (weighted count of R-vector quadruples) %.6g" % totc, "sequence quadruples %.4g (log10 %.2f) of pools log10 %.2f -> retention log10 %.2f" % (totw, math.log10(totw) if totw > 0 else float("nan"), lp, math.log10(totw) - lp if totw > 0 else float("nan")), "%.0fs" % (time.time() - t0), flush=True)
