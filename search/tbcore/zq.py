"""Exact compression (quotient) machinery for zs.ZS structures."""
import itertools, math
from math import log10, pi, gcd
import numpy as np
from zs import *


class Quot:
    """quotient rho: G -> Q = prod Z_{qmods}, rho(g) = (L g) mod qmods (surjective homomorphism)."""
    def __init__(self, st, L, qmods):
        self.st = st
        self.L = np.asarray(L, dtype=np.int64)
        self.qmods = tuple(int(x) for x in qmods)
        self.q = int(np.prod(self.qmods))
        qm = np.array(self.qmods)
        lab = (st.coords @ self.L.T) % qm
        self.lab = np.ravel_multi_index(tuple(lab.T), self.qmods)
        assert len(np.unique(self.lab)) == self.q, "rho not onto"
        # homomorphism check on generators of G
        for j in range(st.d):
            pass
        self.kern = np.flatnonzero(self.lab == self.lab[0])
        self.ksz = len(self.kern)
        # H-invariance: the action of H on Q must be well defined
        self.qcoords = np.array(list(itertools.product(*[range(x) for x in self.qmods])))
        self.qneg = np.ravel_multi_index(tuple((-self.qcoords % qm).T), self.qmods)
        self.qperms = []
        for p in st.perms:
            qp = np.full(self.q, -1)
            for g in range(st.m):
                a = self.lab[g]; b = self.lab[p[g]]
                if qp[a] == -1: qp[a] = b
                elif qp[a] != b: raise ValueError("H does not preserve the kernel")
            self.qperms.append(qp)
        self.fA = st._orbits_q(self.qperms, self.q)
        self.fL = st._orbits_q(self.qperms + [self.qneg], self.q)
        # compressed class structure for PAF on Q: L-orbits on Q
        self.qrep = np.array([np.flatnonzero(self.fL == c)[0] for c in range(self.fL.max() + 1)])
        # target on Q: sum of Tg over fibres
        self.qT = np.zeros(self.q, dtype=np.int64)
        np.add.at(self.qT, self.lab, st.Tg)
        self.qTc = np.array([self.qT[self.qrep[c]] for c in range(len(self.qrep))])
        # per kind: for every block orbit o: fibre class f, weight w
        self.info = {}
        for kind in ("L", "A"):
            orb, sz, free = st.kind_orbits(kind)
            fc = self.fL if kind == "L" else self.fA
            nf = int(fc.max()) + 1
            fsize = np.bincount(fc)
            fid, w = [], []
            for o in free:
                pts = np.flatnonzero(orb == o)
                fs = set(fc[self.lab[pts]].tolist())
                assert len(fs) == 1
                f = fs.pop()
                assert len(pts) % fsize[f] == 0
                fid.append(f); w.append(len(pts) // fsize[f])
            self.info[kind] = (np.array(fid), np.array(w), nf, fsize)
        # fibre sum of the zero element A(0)=0 is already excluded

    def class_counts(self, kind):
        """per fibre class f: dict value -> number of sign assignments of the orbits over f"""
        fid, w, nf, fsize = self.info[kind]
        out = []
        for f in range(nf):
            cnt = {0: 1}
            for o in np.flatnonzero(fid == f):
                nxt = {}
                for t, c in cnt.items():
                    for s in (1, -1):
                        v = t + s * int(w[o]); nxt[v] = nxt.get(v, 0) + c
                cnt = nxt
            out.append(cnt)
        return out

    def comp_paf(self, y):
        """PAF classes on Q for compression vectors y (n x q, values per quotient element)"""
        F = np.fft.fftn(y.reshape((-1,) + self.qmods).astype(float), axes=tuple(range(1, len(self.qmods) + 1)))
        P = np.fft.ifftn(np.abs(F) ** 2, axes=tuple(range(1, len(self.qmods) + 1))).real
        P = P.reshape(P.shape[0], -1)
        return np.rint(P[:, self.qrep]).astype(np.int64)

    def vectors(self, kind, r=None, cap=5_000_000):
        """all compression vectors (per fibre class values) of the kind, with counts and the full-Q vector."""
        cc = self.class_counts(kind)
        fid, w, nf, fsize = self.info[kind]
        vals = [sorted(c.keys()) for c in cc]
        tot = 1
        for v in vals: tot *= len(v)
        if r is not None:
            # prune by row sum sum_f fsize[f]*y_f = r using meet-in-the-middle style enumeration (simple recursion)
            pass
        if tot > cap:
            return None, tot
        grids = np.array(list(itertools.product(*vals)), dtype=np.int64) if nf else np.zeros((1, 0), dtype=np.int64)
        if r is not None:
            rs = grids @ fsize
            grids = grids[rs == r]
        cnts = np.ones(len(grids))
        for f in range(nf):
            cnts *= np.array([cc[f][v] for v in grids[:, f]], dtype=float)
        return (grids, cnts), tot


def _orbits_q(self, perms, n):
    rows = np.concatenate([np.arange(n)] * len(perms)) if perms else np.zeros(0, int)
    cols = np.concatenate(perms) if perms else np.zeros(0, int)
    g = coo_matrix((np.ones(len(rows)), (rows, cols)), shape=(n, n))
    nc, lab = connected_components(g, directed=False)
    order = {}; out = np.zeros(n, dtype=np.int64)
    for i in range(n):
        l = lab[i]
        if l not in order: order[l] = len(order)
        out[i] = order[l]
    return out
ZS._orbits_q = _orbits_q


def feasible(st, qt, rows, cap=5_000_000, verbose=True):
    """exact join of the four compression lists: returns list of tuples (indices) and weights."""
    kinds = ["L", "L", "A", "A"]
    data = []
    for kind, r in zip(kinds, rows):
        res, tot = qt.vectors(kind, r, cap)
        if res is None:
            if verbose: print("  too many compressions", tot)
            return None
        data.append(res)
    fid = {k: qt.info[k] for k in ("L", "A")}
    # expand to full Q vectors
    full = []
    for kind, (g, c) in zip(kinds, data):
        fc = qt.fL if kind == "L" else qt.fA
        full.append((g[:, fc], c, g))
    P = [qt.comp_paf(f[0]) for f in full]
    n = [len(f[1]) for f in full]
    if verbose: print("  compressions per block:", n, "classes on Q:", len(qt.qrep))
    # join: (A+B) vs (T - C - D)
    PAB = (P[0][:, None, :] + P[1][None, :, :]).reshape(-1, P[0].shape[1]) if n[0] * n[1] <= 3_000_000 else None
    PCD = (P[2][:, None, :] + P[3][None, :, :]).reshape(-1, P[0].shape[1]) if n[2] * n[3] <= 3_000_000 else None
    if PAB is None or PCD is None:
        if verbose: print("  pair lists too large", n)
        return None
    tgt = qt.qTc
    # hash
    rng = np.random.default_rng(5)
    key = rng.integers(1, 2**31, size=P[0].shape[1]).astype(np.int64)
    kAB = PAB @ key
    kCD = (tgt[None, :] - PCD) @ key
    order = np.argsort(kAB)
    kAB_s = kAB[order]
    lo = np.searchsorted(kAB_s, kCD, side="left"); hi = np.searchsorted(kAB_s, kCD, side="right")
    sols = []
    wsum = 0.0
    for j in np.flatnonzero(hi > lo):
        for ii in range(lo[j], hi[j]):
            i = order[ii]
            ia, ib = divmod(i, n[1]); ic, idd = divmod(j, n[3])
            if np.array_equal(PAB[i] + PCD[j], tgt):
                sols.append((ia, ib, ic, idd))
                wsum += full[0][1][ia] * full[1][1][ib] * full[2][1][ic] * full[3][1][idd]
    return data, sols, wsum


def sample_cond(st, qt, kind, y, K, rng):
    """sample K orbit-sign vectors whose compression over qt equals the class vector y (len nf)"""
    fid, w, nf, fsize = qt.info[kind]
    orb, sz, free = st.kind_orbits(kind)
    signs = np.zeros((K, len(free)), dtype=np.int8)
    for f in range(nf):
        idxs = np.flatnonzero(fid == f)
        ws = [int(w[i]) for i in idxs]
        tot = sum(ws)
        # suffix DP
        fdp = [None] * (len(ws) + 1)
        cur = np.zeros(2 * tot + 1); cur[tot] = 1.0; fdp[len(ws)] = cur
        for i in range(len(ws) - 1, -1, -1):
            nxt = np.zeros(2 * tot + 1)
            nxt[ws[i]:] += fdp[i + 1][:2 * tot + 1 - ws[i]]
            nxt[:2 * tot + 1 - ws[i]] += fdp[i + 1][ws[i]:]
            fdp[i] = nxt
        r = int(y[f])
        if abs(r) > tot or fdp[0][r + tot] <= 0:
            return None
        rs = np.full(K, r, dtype=np.int64)
        for i, oi in enumerate(idxs):
            ok = (rs - ws[i] + tot >= 0) & (rs - ws[i] + tot <= 2 * tot)
            fp = np.where(ok, fdp[i + 1][np.clip(rs - ws[i] + tot, 0, 2 * tot)], 0.0)
            p = fp / fdp[i][rs + tot]
            u = rng.random(K) < p
            signs[:, oi] = np.where(u, 1, -1)
            rs = np.where(u, rs - ws[i], rs + ws[i])
    return signs


def cond_estimate(st, qt, tup_vecs, counts, K=3000, seed=1, nd=400):
    """lattice Gaussian estimate conditioned on the compression vectors tup_vecs (list of 4 class vectors)"""
    rng = np.random.default_rng(seed)
    kinds = ["L", "L", "A", "A"]
    mus, covs, base, diffs = [], [], None, []
    for kind, y in zip(kinds, tup_vecs):
        S = sample_cond(st, qt, kind, y, K, rng)
        if S is None:
            return None
        V = st.paf_classes(st.seq(kind, S))
        mus.append(V.mean(0)); covs.append(np.cov(V.T.astype(float), bias=True))
        base = V[0].copy() if base is None else base + V[0]
        idx = rng.integers(1, len(V), nd)
        diffs.extend(list(V[idx] - V[0]))
    mu = sum(mus); Sg = sum(covs)
    var = np.diag(Sg)
    tg = st.tvec.astype(float)
    fixed = [c for c in range(st.C) if var[c] <= 1e-9]
    for c in fixed:
        if abs(tg[c] - mu[c]) > 1e-6:
            return ("fixed-mismatch", c, tg[c], mu[c])
    free = [c for c in range(st.C) if var[c] > 1e-9]
    c0 = free[1:]
    # drop the class with the dependency: the dependency is Sum_c size_c PAF_c = const; dropping any one free class is ok
    D = [[int(dv[c]) for c in c0] for dv in diffs]
    B, piv = hnf_basis(D, len(c0))
    covol = 1
    for k, b in enumerate(B):
        covol *= b[piv[k]]
    t = [int(st.tvec[c] - base[c]) for c in c0]
    for k, b in enumerate(B):
        j = piv[k]
        if t[j] % b[j]:
            return ("coset",)
        q = t[j] // b[j]
        t = [x - q * y for x, y in zip(t, b)]
    if any(t):
        return ("coset",)
    pc = [c0[p] for p in piv]
    dvec = (tg - mu)[pc]
    Sk = Sg[np.ix_(pc, pc)]
    sign, logdet = np.linalg.slogdet(Sk)
    quad = float(dvec @ np.linalg.solve(Sk, dvec))
    dens = (-0.5 * len(pc) * math.log(2 * pi) - 0.5 * logdet - 0.5 * quad) / math.log(10)
    lN = sum(log10(c) for c in counts)
    return (lN + dens + log10(abs(covol)), lN, dens, log10(abs(covol)), len(pc), len(c0), quad)
