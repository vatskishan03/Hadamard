#!/usr/bin/env python3
"""The tower for G = Z2^4 x Z11: context, R-tuples, stage 1 (class-level join), tuple priority, tower.c specs, hit decoding.

Library module (drivers: produce.py, test_tower.py).  One R-tuple (exact, from rlevel/rtuples) fixes the row sums r_k of the
four blocks over E = Z2^4: the A and B fibres are then complete (10^5 to 10^6 members each) and C and D are huge (10^17).
Stage 1 (`stage1`) lists the A and B fibres completely, draws `ns` uniform members of the C and D fibres, groups all four
by the column-weight vector w (the Z11 level: six PAF columns) and joins the classes exactly with pjoin; every match is an
"(R,w)-tuple" (class of A, class of B, w of C, w of D).  `tuple_priority` ranks the (R,w)-tuples by the estimated size of
their fibre product.  One (R,w)-tuple is then resolved completely by tower.c (children at line 1 and line 2, fibre
enumeration, exact four-list join on all 35 PAF columns): `write_spec` writes its input, `run_tower` runs it and reads the
hits.  Every hit is a full solution of all 35 class equations; `verify` checks it again by the definition of the PAF over G.
"""
import os, math, struct, subprocess
import numpy as np
from structure import Struct
from rlevel import RLevel
from rtuples import explicit_tuples, RSampler
from qlevel import tuples, rs
from fibres import Fib, LET
from lines import line_matrix, line_targets
from pio import write_lists, run_pjoin, read_out

HERE = os.path.dirname(os.path.abspath(__file__))


def tower_tables(S, RL):
    """K[u][v][c'] = sum_e f_u(e) f_v(e + e'_c) for the 64 letters u, v (f_u = class signs on E) and the six class representatives e'_c;
    colmap[c' * 6 + s] = position of the PAF class of g = (e'_c, s) among the 35 columns (-1 for the class of 0)"""
    cl = RL.cls["C"]
    fu = LET[:, cl].astype(np.int64)                      # (64, 16): letter -> function on E (class signs)
    add = RL.add                                          # add[e,f] = e+f index
    reps = [int(np.flatnonzero(cl == c)[0]) for c in range(6)]
    K = np.zeros((64, 64, 6), dtype=np.int16)
    for ci, e1 in enumerate(reps):
        sh = add[np.arange(16), e1]                        # index of e + e'
        K[:, :, ci] = (fu[:, None, :] * fu[None, :, :][:, :, sh]).sum(2)
    colpos = {c: i for i, c in enumerate(S.cols)}
    cm = -np.ones(36, dtype=np.int32)
    G = S.coords
    for ci, e1 in enumerate(reps):
        for s in range(6):
            g = int(np.flatnonzero((RL.piE == e1) & (G[:, -1] == s))[0])
            cls = int(S.orbA[g])
            if cls in colpos: cm[ci * 6 + s] = colpos[cls]
    assert (cm >= 0).sum() == len(S.cols) == 35
    return K, cm


class Ctx:
    def __init__(self):
        self.S = Struct("i"); self.RL = RLevel(self.S); self.RS = RSampler(self.S, self.RL); self.F = Fib(self.S, self.RL)
        self.K, self.colmap = tower_tables(self.S, self.RL)
        sols, self.dA, self.dC, tg = tuples(self.S)
        self.sols = sols
        rows = line_matrix(self.S, self.RL); self.rowsM = list(rows.keys())
        self.row1 = np.array(self.rowsM[1]); self.row2 = np.array(self.rowsM[2])
        # line targets
        self.etas = [v[0] for v in rows.values()]
        self.tau = line_targets(self.S, self.RL, self.etas)
        self.tZ = np.array([int(self.S.target[self.S.coords[:, -1] == s].sum()) for s in range(6)])
        assert (self.tau[0] == self.tZ).all()
        self.symidx = np.array([min(z, 11 - z) for z in range(11)])

    def qtuple(self, pat):
        for t in self.sols:
            if tuple(abs(rs(y)) for y in t) == pat:
                return tuple(y if rs(y) >= 0 else (-y[0], -y[1]) for y in t)
        raise KeyError(pat)

    def rtuples(self, pat):
        """the explicit R-tuples of a row-sum pattern, in the order of the production run: returns (L, tl, wts, uniq)
        with uniq = indices into tl of the tuples with distinct log10 weight, heaviest first (the "R-tuple rank")"""
        tup = self.qtuple(pat); L, tl = explicit_tuples(self.RL, tup)
        wts = np.array([np.sum([math.log10(L[q][1][t[q]]) for q in range(4)]) for t in tl])
        order = np.argsort(-wts); seen = set(); uniq = []
        for i in order:
            key = round(float(wts[i]), 6)
            if key in seen: continue
            seen.add(key); uniq.append(int(i))
        return L, tl, wts, uniq


def cells_to_X_C(F, cells):
    n = len(cells); X = np.zeros((n, 6, 11), dtype=np.int8)
    for c in range(6):
        for z in range(11): X[:, c, z] = cells[:, F.cellC[c, z]]
    return X


class RTuple:
    """one exact R-tuple: complete A,B fibres and their w-classes (w = column weights over Z11)"""
    def __init__(self, ctx, pat, idx, L, tl):
        self.ctx = ctx; F = ctx.F
        self.L, self.tl, self.pat = L, tl, pat
        tt = tl[idx]; self.tt = tt; self.idx = idx
        self.r = [L[q][0][tt[q]] for q in range(4)]
        self.fib = [L[q][1][tt[q]] for q in range(4)]
        self.XA = F.enum_A(self.r[0]); self.XB = F.enum_A(self.r[1])
        self.wA = F.wvec_A(self.XA); self.wB = F.wvec_A(self.XB)
        self.uA, self.invA = np.unique(self.wA, axis=0, return_inverse=True); self.invA = self.invA.ravel()
        self.uB, self.invB = np.unique(self.wB, axis=0, return_inverse=True); self.invB = self.invB.ravel()
        self.cntA = np.bincount(self.invA); self.cntB = np.bincount(self.invB)


def pafz11(w11):
    return np.stack([(w11 * np.roll(w11, -s, axis=1)).sum(1) for s in range(6)], 1)


def stage1(ctx, R, ns, seed, workdir, nbits=4):
    """class-level join on the six Z11-level columns: A,B classes complete, C,D classes from ns uniform samples of the R-fibres;
    returns the array quads (a_class, b_class, c_index, d_index) and the class lists [wC, wD] (rows = w vectors)"""
    F = ctx.F; rng = np.random.default_rng(seed)
    cl = []
    for q in (2, 3):
        s = ctx.RS.sample("C", R.r[q], ns, rng)
        X = cells_to_X_C(F, s); w = F.wvec_C(X)
        u = np.unique(w, axis=0)
        cl.append(u)
    PA = pafz11(R.uA[:, ctx.symidx]); PB = pafz11(R.uB[:, ctx.symidx]); PC = pafz11(cl[0]); PD = pafz11(cl[1])
    fn = os.path.join(workdir, "s1_lists.bin"); write_lists(fn, [PA, PB, PC, PD], ctx.tZ)
    out = os.path.join(workdir, "s1_join.txt"); run_pjoin(fn, 6, nbits, 5, 0, 1, 6, out)
    hist, quads, stats = read_out(out, 6, 6)
    q = np.array([x[1:] for x in quads], dtype=np.int64).reshape(-1, 4)
    return q, cl


def acc_proxy(ctx, w, r):
    """log10 of a Gaussian-local-limit proxy for P(R = r | letters uniform on the weight class of each column)"""
    F = ctx.F
    mu = np.zeros(6); cov = np.zeros((6, 6)); lg = 0.0
    for wz in w:
        lets = LET[F.byw[int(wz)]].astype(np.float64)
        mu += lets.mean(0); cov += np.cov(lets.T, bias=True) if len(lets) > 1 else 0
        lg += math.log10(len(lets))
    d = np.asarray(r) - mu
    ev, U = np.linalg.eigh(cov + 1e-12 * np.eye(6)); keep = ev > 1e-6
    q = float(((U[:, keep].T @ d) ** 2 / ev[keep]).sum())
    return lg - 0.5 * float(np.log10(ev[keep]).sum()) - 0.5 * q / math.log(10)


def tuple_priority(ctx, R, cl, quads):
    """estimated log10 fibre product of each class-level tuple (C,D by the Gaussian proxy, A,B exact); only used for the order"""
    cacheC = {}; cacheD = {}
    pr = np.zeros(len(quads))
    for i, (a, b, c, d) in enumerate(quads):
        if c not in cacheC: cacheC[c] = acc_proxy(ctx, cl[0][c], R.r[2])
        if d not in cacheD: cacheD[d] = acc_proxy(ctx, cl[1][d], R.r[3])
        pr[i] = math.log10(R.cntA[a]) + math.log10(R.cntB[b]) + cacheC[c] + cacheD[d]
    return pr


def prepare_members(ctx, R):
    """per A and B member: its 35 PAF class values, its line-1 and line-2 labels (a1, a2); class -> member lists (ordA/stA, ordB/stB)"""
    S, F = ctx.S, ctx.F
    out = []
    for X in (R.XA, R.XB):
        paf = np.zeros((len(X), 35), dtype=np.int16)
        for s in range(0, len(X), 50000):
            paf[s:s + 50000] = S.paf_classes(S.expand(F.cells_A(X[s:s + 50000]), "A"))[:, S.cols]
        a1 = np.einsum("ncz,c->nz", X.astype(np.int64), ctx.row1).astype(np.int8)
        a2 = np.einsum("ncz,c->nz", X.astype(np.int64), ctx.row2).astype(np.int8)
        out.append((paf, a1, a2))
    R.mem = out
    for nm, inv, nc in (("A", R.invA, len(R.uA)), ("B", R.invB, len(R.uB))):
        o = np.argsort(inv, kind="stable"); cnt = np.bincount(inv, minlength=nc); st = np.concatenate([[0], np.cumsum(cnt)])
        setattr(R, "ord" + nm, o); setattr(R, "st" + nm, st)


def write_spec(ctx, R, fn, a, b, rC, wC, rD, wD, T=None, tau1=None, tau2=None):
    """input file of tower.c for the (R,w)-tuple (class a of A, class b of B, C with row sums rC and weights wC, D with rD, wD).
    T, tau1, tau2 default to the real targets (a planted test passes its own).  Returns the member index lists (into R.XA, R.XB)."""
    ia = R.ordA[R.stA[a]:R.stA[a + 1]]; ib = R.ordB[R.stB[b]:R.stB[b + 1]]
    with open(fn, "wb") as f:
        f.write(struct.pack("3i", 35, len(ia), len(ib)))
        for (paf, a1, a2), ix in ((R.mem[0], ia), (R.mem[1], ib)):
            rec = np.concatenate([paf[ix].view(np.uint8).reshape(len(ix), 70), a1[ix].view(np.uint8), a2[ix].view(np.uint8)], axis=1)
            f.write(np.ascontiguousarray(rec).tobytes())
        f.write(np.asarray(ctx.S.T if T is None else T, dtype=np.int32).tobytes())
        f.write(np.asarray(ctx.tau[1] if tau1 is None else tau1, dtype=np.int32).tobytes())
        f.write(np.asarray(ctx.tau[2] if tau2 is None else tau2, dtype=np.int32).tobytes())
        f.write(ctx.K.astype(np.int16).tobytes()); f.write(ctx.colmap.astype(np.int32).tobytes()); f.write(LET.astype(np.int8).tobytes())
        f.write(np.asarray(ctx.row1, dtype=np.int32).tobytes()); f.write(ctx.row2.astype(np.int32).tobytes())
        f.write(np.asarray(rC, dtype=np.int32).tobytes()); f.write(np.asarray(wC, dtype=np.int32).tobytes())
        f.write(np.asarray(rD, dtype=np.int32).tobytes()); f.write(np.asarray(wD, dtype=np.int32).tobytes())
    return ia, ib


def run_tower(spec, out, exe=None):
    """run tower.c on a spec; returns (hits, stats line, level-1 line); a hit is (index of A member, index of B member, 11 letters of C, 11 letters of D)"""
    exe = exe or os.path.join(HERE, "tower")
    subprocess.run([exe, spec, out], check=True)
    hits = []; stats = None; l1 = None
    for l in open(out):
        if l.startswith("HIT"):
            v = list(map(int, l.split()[1:])); hits.append((v[0], v[1], v[2:13], v[13:24]))
        elif l.startswith("# STATS"): stats = l.strip()
        elif l.startswith("# level1"): l1 = l.strip()
    return hits, stats, l1


def hit_to_sequences(ctx, h):
    """h = dict(A=(6,6) signs, B=(6,6), C=11 letters, D=11 letters) -> the four sequences over G (length 176)"""
    S, F = ctx.S, ctx.F
    XA = np.array(h["A"], dtype=np.int8)[None]; XB = np.array(h["B"], dtype=np.int8)[None]
    XC = LET[np.array(h["C"])].T[None].astype(np.int8); XD = LET[np.array(h["D"])].T[None].astype(np.int8)
    return [S.expand(F.cells_A(XA), "A")[0], S.expand(F.cells_A(XB), "A")[0], S.expand(F.cells_C(XC), "C")[0], S.expand(F.cells_C(XD), "C")[0]]


def verify(ctx, seqs):
    """exact: sum of the four PAFs over G (by the definition) equals the target at every g"""
    S = ctx.S
    tot = sum(S.paf_direct(x) for x in seqs)
    return bool((tot == S.target).all())
