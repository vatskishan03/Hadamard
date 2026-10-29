"""plain_spence3.py -- EXHAUSTIVE search of plain multiplier structures for Spence's bordered construction over the cyclic group Z_v, v = 2p-1 ODD.
Four +-1 sequences on Z_v constant on the orbits of a subgroup H <= Z_v^*, every row sum +1 after negating blocks (row sums are +-1 each, sum r^2 = 4),
sum_X PAF_X(g) = -4 for every g != 0  (=> Hadamard matrix of order 4(v+1) = 8p; VALIDATED construction in oddspence/).
Algorithm: orbits of H on Z_v; classes of nonzero orbits under +-H; per class a quadratic form PAF_c(x) = x^T M_c x on the orbit signs; list L = {x : row sum = +1};
exact 4-sum on L (hash join of pair sums by bucket).  Counts ordered quadruples exactly (up to the block-negation gauge: all row sums +1).
usage: plain_spence3.py v [Omax [Omin [maxL]]]   -> enumerates every subgroup H of Z_v^* with Omin <= #orbits <= Omax and reports L and the exact count.
A large pair bucket is split by the residue of an exact small-integer linear form, which is additive, and the pairs are generated in chunks.  (A residue of the 64-bit hash would not be additive: the int64 hash sums
wrap around, so the residues of A + B would differ from those of the target and matches would be lost.  An earlier version split by (64-bit hash mod nparts) with nparts not a power of two; the shapes it could have
affected were re-solved with this one, search/NOTES.md section 92.)
Validated against an independent dictionary join (test_plain.py) with and without forced splitting.
"""
import sys, math, itertools, time
import numpy as np
from sympy import factorint


def units(v): return [u for u in range(1, v) if math.gcd(u, v) == 1]


def closure(gens, v):
    S = {1}; fr = [1]
    while fr:
        x = fr.pop()
        for g in gens:
            y = x * g % v
            if y not in S: S.add(y); fr.append(y)
    return frozenset(S)


def sylow_parts(v):
    """the unit group of Z_v split into Sylow subgroups (as sets of elements of Z_v^*)"""
    U = units(v); n = len(U)
    f = factorint(n)
    parts = {}
    for p, e in f.items():
        pe = p ** e
        # elements of p-power order: u^(n/pe * ?)... the p-part = { u : u^pe = 1 }
        parts[p] = [u for u in U if pow(u, pe, v) == 1]
    return parts


def subgroups_of(elems, v):
    """all subgroups of the abelian p-group given by its element list (brute force over <=3 generators)"""
    subs = {frozenset([1])}
    el = [e for e in elems if e != 1]
    # cyclic subgroups, then joins
    cyc = {closure([e], v) for e in el}
    subs |= cyc
    cur = set(cyc)
    while True:
        new = set()
        for A in cur:
            for C in cyc:
                if not C <= A:
                    J = closure(list(A | C), v)
                    if J not in subs: new.add(J)
        if not new: break
        subs |= new; cur = new
    return subs


def all_subgroups(v):
    parts = sylow_parts(v)
    per = []
    for p, el in parts.items():
        per.append(list(subgroups_of(el, v)))
    out = {}
    for combo in itertools.product(*per):
        gens = []
        for S in combo: gens += list(S)
        H = closure(gens, v) if gens else frozenset([1])
        out[H] = True
    return list(out)


def orbits(v, H):
    lab = -np.ones(v, dtype=np.int64); reps = []
    Hl = sorted(H)
    for g in range(v):
        if lab[g] < 0:
            o = len(reps); reps.append(g)
            for h in Hl: lab[(g * h) % v] = o
    return lab, reps


def structure(v, H):
    lab, reps = orbits(v, H)
    O = len(reps); size = np.bincount(lab, minlength=O)
    # classes of nonzero orbits under negation
    neg = lab[(-np.array(reps)) % v]
    seen = set(); cls = []
    for a in range(1, O):
        if a in seen: continue
        seen.add(a); seen.add(int(neg[a])); cls.append(a)
    return lab, reps, O, size, cls


def forms(v, lab, reps, O, cls):
    """M_c[a, b] = #{h : lab[h] = a, lab[h + g_c] = b}  for the class representative g_c = reps[c]"""
    Ms = []
    idx = np.arange(v)
    for a in cls:
        g = reps[a]
        sh = (idx + g) % v
        M = np.zeros((O, O), dtype=np.int64)
        np.add.at(M, (lab[idx], lab[sh]), 1)
        Ms.append(M)
    return Ms


def enum_list(O, size, rowsum=1, chunk=1 << 18):
    """all x in {+-1}^O with sum size*x = rowsum, as an int8 array (nrows, O)"""
    out = []
    N = 1 << O
    bits = np.arange(O, dtype=np.int64)
    for s in range(0, N, chunk):
        ids = np.arange(s, min(N, s + chunk), dtype=np.int64)
        X = 1 - 2 * ((ids[:, None] >> bits[None, :]) & 1)
        rs = X @ size
        sel = rs == rowsum
        if sel.any(): out.append(X[sel].astype(np.int8))
    return np.concatenate(out) if out else np.zeros((0, O), np.int8)


def solve(V, w, target, maxpairs=6e6):
    """exact number of ordered quadruples of sequences (rows of V with multiplicities w) whose vectors sum to target.
    Unordered pairs are bucketed by the exact value s = V[i,0] + V[j,0] of the first class; a pair with s must meet a pair with target[0] - s; inside a bucket pair the 64-bit hash of the
    sum is joined by sorting, candidates verified exactly.  Memory O(pairs per bucket)."""
    L = len(V); C = V.shape[1]
    rng = np.random.default_rng(12345)
    coef = rng.integers(1, 2 ** 61, C, dtype=np.int64) | 1
    h = (V.astype(np.int64) * coef).sum(1)
    ht = int((np.asarray(target, dtype=np.int64) * coef).sum())
    w = np.asarray(w, dtype=np.int64); tg = np.asarray(target, dtype=np.int64)
    v0 = V[:, 0].astype(np.int64)
    order = np.argsort(v0, kind='stable'); vals, starts = np.unique(v0[order], return_index=True)
    ends = np.append(starts[1:], L)
    grp = {int(vals[k]): order[starts[k]:ends[k]] for k in range(len(vals))}
    CH = 1500000
    acoef = np.random.default_rng(777).integers(1, 10 ** 6, C, dtype=np.int64)
    rho_all = V.astype(np.int64) @ acoef          # exact (small integers); residues mod nparts are additive
    rt_all = int(tg @ acoef)
    def pairs_for(s, part=0, nparts=1):
        """all unordered pairs (i <= j) with v0[i] + v0[j] = s, restricted to hash residue `part` mod nparts; chunked generation (peak memory ~ CH entries); returns hash sums, i, j, weights"""
        outH = []; outI = []; outJ = []; outW = []
        def emit(I, J):
            ph = h[I] + h[J]
            if nparts > 1:
                m = ((rho_all[I] + rho_all[J]) % nparts) == part; I = I[m]; J = J[m]; ph = ph[m]
            if len(I) == 0: return
            outH.append(ph); outI.append(I.astype(np.int32)); outJ.append(J.astype(np.int32)); outW.append(w[I] * w[J] * np.where(I == J, 1, 2))
        for a, ga in grp.items():
            b = s - a
            if b < a or b not in grp: continue
            gb = grp[b]
            if a == b:
                n = len(ga); i0 = 0
                while i0 < n:
                    # rows i0..i1 of the upper triangle (j >= i): about (n - i) entries per row
                    i1 = i0; cnt = 0
                    while i1 < n and cnt + (n - i1) <= CH: cnt += n - i1; i1 += 1
                    if i1 == i0: i1 = i0 + 1
                    ii = np.concatenate([np.full(n - i, i, dtype=np.int64) for i in range(i0, i1)]); jj = np.concatenate([np.arange(i, n, dtype=np.int64) for i in range(i0, i1)])
                    emit(ga[ii], ga[jj]); i0 = i1
            else:
                step = max(1, CH // len(gb))
                for r0 in range(0, len(ga), step):
                    gr = ga[r0:r0 + step]
                    emit(np.repeat(gr, len(gb)), np.tile(gb, len(gr)))
        if not outH: return None
        return np.concatenate(outH), np.concatenate(outI), np.concatenate(outJ), np.concatenate(outW)
    def npairs(s):
        n = 0
        for a, ga in grp.items():
            b = s - a
            if b < a or b not in grp: continue
            n += len(ga) * (len(ga) + 1) // 2 if a == b else len(ga) * len(grp[b])
        return n
    total = 0; found = []
    svals = sorted(set(a + b for a in grp for b in grp if a <= b))
    t0v = int(tg[0])
    for s in svals:
        s2 = t0v - s
        if s2 < s: continue            # handle each unordered {s, s2} once (s2 >= s)
        if s2 not in set(svals): continue
        na, nb_ = npairs(s), npairs(s2)
        nparts = int(max(na, nb_) / maxpairs) + 1      # residue classes mod nparts of the exact integer a.V (additive; NOT of the wrapped 64-bit hash)
        for part in range(nparts):
            if nparts > 1:
                part2 = (rt_all - part) % nparts
            else: part2 = 0
            A = pairs_for(s, part, nparts)
            B = A if (s2 == s and part2 == part) else pairs_for(s2, part2, nparts)
            if A is None or B is None: continue
            Ah, Ai, Aj, Aw = A; Bh, Bi, Bj, Bw = B
            oB = np.argsort(Bh, kind='stable'); Bs = Bh[oB]
            need = ht - Ah
            lo = np.searchsorted(Bs, need, 'left'); hi = np.searchsorted(Bs, need, 'right')
            mult = 1 if s2 == s else 2      # ordered pairs of pair-buckets: (s,s2) and (s2,s) both
            for t in np.flatnonzero(hi > lo):
                for k in range(lo[t], hi[t]):
                    q = oB[k]
                    if np.array_equal(V[Ai[t]] + V[Aj[t]] + V[Bi[q]] + V[Bj[q]], tg):
                        c = int(Aw[t]) * int(Bw[q])
                        total += c if (s2 != s) else c
                        if s2 != s: total += c
                        if len(found) < 20: found.append((int(Ai[t]), int(Aj[t]), int(Bi[q]), int(Bj[q])))
    return total, found


def run(v, Omax, Omin=2, maxL=20000, verbose=True):
    subs = all_subgroups(v)
    shapes = []
    for H in subs:
        lab, reps, O, size, cls = structure(v, H)
        if Omin <= O <= Omax: shapes.append((O, len(H), H))
    shapes.sort(key=lambda x: (x[0], -x[1]))
    if verbose: print(f'v = {v}: {len(subs)} subgroups of Z_v^*, {len(shapes)} with {Omin} <= orbits <= {Omax}', flush=True)
    res = []
    for O, h, H in shapes:
        t0 = time.time()
        lab, reps, O, size, cls = structure(v, H)
        Ms = forms(v, lab, reps, O, cls)
        X = enum_list(O, size)
        rec = dict(v=v, O=O, h=h, gens=sorted(H)[:3], classes=len(cls), L=len(X))
        if len(X) == 0 or len(X) > maxL * 4:
            rec['status'] = 'noL' if len(X) == 0 else 'bigL'
            res.append(rec);
            if verbose: print(rec, flush=True)
            continue
        V = np.stack([np.einsum('ni,ij,nj->n', X.astype(np.int64), M, X.astype(np.int64)) for M in Ms], 1)
        target = np.full(len(cls), -4, dtype=np.int64)
        # quick necessary: PAF values mod 4 (PAF = v mod 4 when row sum..) skip
        Vu, inv = np.unique(V, axis=0, return_inverse=True)
        cnt = np.bincount(inv.ravel())
        tot, found = solve(Vu, cnt, target)
        rec.update(Lu=len(Vu), status='EMPTY' if tot == 0 else 'SOLUTION', count_rho_plus=int(tot), secs=round(time.time() - t0, 1))
        if tot:
            rec['found'] = [[X[np.flatnonzero(inv.ravel() == q)[0]].tolist() for q in quad] for quad in found[:1]]
            rec['orbit_sizes'] = size.tolist()
        res.append(rec)
        if verbose: print({k: val for k, val in rec.items() if k not in ('found',)}, flush=True)
        if tot: print('FOUND', rec, flush=True)
    return res


if __name__ == '__main__':
    v = int(sys.argv[1]); Omax = int(sys.argv[2]) if len(sys.argv) > 2 else 16; Omin = int(sys.argv[3]) if len(sys.argv) > 3 else 2
    maxL = int(float(sys.argv[4])) if len(sys.argv) > 4 else 20000
    run(v, Omax, Omin, maxL=maxL)
