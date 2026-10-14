#!/usr/bin/env python3
"""Low-energy k-tree for Szöllősi structures (cyclic on Z_N, or bordered on Z_(N-1)) whose blocks are too large to
join exactly.

Each block k (A, B symmetric; C, D arbitrary; all constant on the orbits of a multiplier group H) is enumerated in full
by the depth-first search of search/lowe_enum.c, and only its members of low energy are kept:

    E_k(x) = (x - z_k)^T Sigma_{-k}^+ (x - z_k),     z_k = t - sum_{j != k} mu_j,

where x is the member's vector of class values, mu_j, Sigma_j the mean and covariance of block j's class values over
its pool (from uniform samples), Sigma_{-k} = sum_{j != k} Sigma_j and t the target (0 cyclic, -4 bordered). Under the
Gaussian model, exp(-E_k(x)/2) is proportional to the density of the other three blocks' sum at t - x, so the kept
members are those that most likely complete to a solution. The kept lists are joined by search/kt4.c (the four-list
k-tree of the low-energy skew searches) with the pairing (A, C), (B, D); every solution is rebuilt and checked exactly.

The estimate printed before the join is E (Gaussian, whole structure) times the retention of each filtered block,
the probability under the model that a solution's member passes the filter (Monte Carlo from the conditional law
N(c_k, Sigma_k - Sigma_k Sigma^+ Sigma_k)).

Before any enumeration the compressions to Z_d (d a prime power dividing the group order, d <= 64) are joined exactly
(search/compress_check.py); an empty join rules the structure/pattern out (--no-precheck skips this). A solution's
compressions form a level-1 solution, so each block keeps only members whose compressions to Z_d (d <= 32) occur in
the projection of the level-1 join (fibre filters, applied in lowe_enum before the energy; --no-fibre skips them);
the statistics, the energy quantile and the retention are those of the allowed sub-pool.

    python3 search/lowe_sz.py N g [g2] --rows a,b,c,d [--bordered] [--L 2e7] [--nsig 2e5] [--S 5] [--estimate-only]
        [--energy full|units] [--time SECONDS] [--resume] [--no-precheck]
"""
import argparse, json, math, os, struct, subprocess, sys, time
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import workdir as W
from bnw import Struct
import ssxx
import bsz

KINDS = ["L", "L", "A", "A"]


def stats(st, kind, r, K, rng):
    Y = st.sample(kind, r, K, rng)
    V = st.values(kind, Y)
    return V.mean(0), np.cov(V.T), V


def psd_pinv_sqrt(M, rel=1e-9):
    w, U = np.linalg.eigh(M)
    keep = w > rel * w.max()
    inv = np.zeros_like(w); inv[keep] = 1 / w[keep]
    return (U * np.sqrt(inv)) @ U.T, (U * inv) @ U.T


def fibre_filters(st, rows, bordered, dmax=32, maxstates=2_000_000, seed=5):
    """allowed compressions of each block. For every prime power d | m (d <= dmax) whose level-1 join
    (search/compress_check.py) can be listed, a solution's four compressions to Z_d form a level-1 solution, so block
    k's compression lies in the projection of the join. Returns, per block, a list of (d, K, keys, mass): K[o] the
    64-bit key contribution of orbit o, keys the sorted keys of the allowed compressions, mass the fraction of the
    block's pool they hold; None when a join is empty (no solution)."""
    import cchain4 as CC, chain as CH
    kinds, T, _ = CC.MODES["bsz" if bordered else "ssxx"]
    v = st.m
    ps = [p for p in range(2, v) if v % p == 0 and all(p % q for q in range(2, int(p ** 0.5) + 1))]
    cand = sorted({p ** e for p in ps for e in range(1, 8) if v % (p ** e) == 0 and p ** e < v and p ** e <= dmax})
    rng = np.random.default_rng(seed)
    out = [[] for _ in range(4)]
    t = np.arange(v)
    for d in cand:
        c1 = CC.KChain(v, list(st.gens), [d], T, kinds)
        Sd = [c1.states(k, r, maxstates) for k, r in zip(kinds, rows)]
        if any(x is None for x in Sd):
            continue
        U = [K_.astype(np.int64) for K_, _ in Sd]
        try:
            J = CH.join4([c1.paf(0, u) for u in U], c1.target(0)) if all(len(u) for u in U) else []
        except CH.JoinTooLarge:
            continue
        if not J:
            return None
        reps = c1.lev[0][2]
        R = rng.integers(1, 2 ** 63, size=len(reps), dtype=np.int64).astype(np.uint64)
        for k in range(4):
            al = sorted({tt[k] for tt in J})
            mass = float(Sd[k][1][al].sum() / Sd[k][1].sum())
            if mass > 0.98:
                continue
            orb = st.orb(KINDS[k]); O = len(st.sizes(KINDS[k]))
            P = np.zeros((O, len(reps)), np.int64)
            for j, rr in enumerate(reps):
                P[:, j] = np.bincount(orb[(t % d) == rr], minlength=O)
            Ko = (P.astype(np.uint64) * R[None, :]).sum(1, dtype=np.uint64)
            keys = np.unique((U[k][al].astype(np.uint64) * R[None, :]).sum(1, dtype=np.uint64))
            out[k].append((d, Ko, keys, mass, P, R))
    return out


def fibre_ok(filt, Y):
    """members (rows of Y, signs on the orbits) whose compressions pass every filter of the block"""
    ok = np.ones(len(Y), bool)
    for d, Ko, keys, mass, P, R in filt:
        comp = np.rint(Y).astype(np.int64) @ P
        key = (comp.astype(np.uint64) * R[None, :]).sum(1, dtype=np.uint64)
        ok &= np.isin(key, keys)
    return ok


def fibre_sample(st, kind, r, K, rng, filt, maxdraw=40):
    """K uniform members of the block's allowed sub-pool (rejection), and the acceptance rate"""
    if not filt:
        return st.sample(kind, r, K, rng), 1.0
    got, n = [], 0
    for _ in range(maxdraw):
        Y = st.sample(kind, r, K, rng); n += len(Y)
        Y = Y[fibre_ok(filt, Y)]
        got.append(Y)
        if sum(len(y) for y in got) >= K:
            break
    Y = np.concatenate(got)[:K]
    acc = sum(len(y) for y in got) / max(n, 1)
    return Y, acc


def enum_block(st, kind, r, center, Wm, tau, maxkeep, tag, threads, filt=()):
    """run lowe_enum; returns (masks, values int16 [n, C], energies)"""
    d = W.WORK / "tmp"; d.mkdir(exist_ok=True)
    inp = d / f"lowe_{st.N}_{tag}_{os.getpid()}.in"; out = d / f"lowe_{st.N}_{tag}_{os.getpid()}.out"
    sz = st.sizes(kind); S = st.S[kind]
    with open(inp, "wb") as f:
        f.write(struct.pack("4i", len(sz), st.C, int(r), int(min(maxkeep, 2 ** 31 - 1))))
        f.write(np.asarray(sz, np.int32).tobytes()); f.write(np.asarray(S, np.int32).tobytes())
        f.write(np.asarray(center, np.float32).tobytes()); f.write(np.asarray(Wm, np.float32).tobytes())
        f.write(struct.pack("f", float(tau)))
        f.write(struct.pack("i", len(filt)))
        for d_, Ko, keys, *_ in filt:
            f.write(np.asarray(Ko, np.uint64).tobytes()); f.write(struct.pack("q", len(keys)))
            f.write(np.asarray(keys, np.uint64).tobytes())
    p = subprocess.run([W.bin("lowe_enum"), str(inp), str(out), str(threads)], capture_output=True, text=True, check=True)
    inp.unlink()
    with open(out, "rb") as f:
        n = struct.unpack("q", f.read(8))[0]
        masks = np.frombuffer(f.read(8 * n), np.uint64)
        vals = np.frombuffer(f.read(2 * n * st.C), np.int16).reshape(n, st.C)
        en = np.frombuffer(f.read(4 * n), np.float32)
    out.unlink()
    return masks, vals, en, p.stdout.strip()


def resume(N, gens, rows, bordered=False, nsig=2e6, nS=5, log=None, maxsol=5, energy="full", time_limit=None,
           kt4_threads=1):
    """rerun the join on the lists saved by an earlier run (WORK/lowe), with more split values"""
    st = Struct(N, gens, plain=not bordered)
    keepdir = W.WORK / "lowe"
    tagname = f"{N}{'b' if bordered else ''}_{'_'.join(map(str, gens))}_{'_'.join(map(str, rows))}" + \
        ("" if energy == "full" else f"_{energy}")
    inp = keepdir / f"kt4_{tagname}.bin"
    z = np.load(keepdir / f"masks_{tagname}.npz")
    masks = [z[f"arr_{k}"] for k in range(4)]
    with open(inp, "rb") as f:
        hdr = struct.unpack("5i", f.read(20)); C = hdr[4]
    vals, off = [], 20 + 2 * C
    for n in hdr[:4]:                               # memory-mapped, not copied
        vals.append(np.memmap(inp, np.int16, "r", off, (n, C))); off += 2 * C * n
    order = [0, 2, 1, 3]
    lists = [None] * 4
    for pos, k in enumerate(order):
        lists[k] = (masks[k], vals[pos])
    return join(st, N, gens, rows, bordered, lists, inp, nsig, nS, log, maxsol, time.time(), time_limit=time_limit,
                kt4_threads=kt4_threads)


def unit_classes(st):
    """classes whose shifts are prime to the group order (the others are compression classes)"""
    return [c for c, s in enumerate(st.reps) if math.gcd(int(s), st.m) == 1]


def run(N, gens, rows, bordered=False, L=2e7, nsig=2e5, nS=5, K=20000, threads=4, estimate_only=False, log=None,
        seed=1, maxsol=5, energy="full", time_limit=None, precheck=True, fibre=True, kt4_threads=1):
    """energy="units": the filter's energy uses only the unit classes. On exact solution sets of composite N the
    solution members sit at arithmetically fixed values on the compression classes, away from the Gaussian centre;
    without those classes they rank much higher at N = 119 and 267 (search/NOTES.md section 34). fibre: the blocks
    keep only members whose compressions occur in a level-1 solution (fibre_filters; section 37)."""
    t0 = time.time()
    if precheck:
        import compress_check
        r = compress_check.check("bsz" if bordered else "ssxx", N, list(gens), tuple(rows), 3_000_000, 64, stop=True)
        dead = [d for d, n in r.items() if n == 0]
        if dead:
            print(f"N = {N}{' bordered' if bordered else ''} H = <{','.join(map(str, gens))}> rows {tuple(rows)}: "
                  f"ruled out by the compression to Z_{dead[0]} (level-1 joins {r}) [{time.time() - t0:.0f}s]",
                  flush=True)
            return []
    st = Struct(N, gens, plain=not bordered)
    rng = np.random.default_rng(seed)
    C = st.C
    t = np.full(C, float(st.target))
    pools = [st.counts(k).get(r, 0) for k, r in zip(KINDS, rows)]
    if min(pools) == 0:
        print(f"N = {N} rows {rows}: empty block"); return []
    # fibre filters: each block keeps only members whose compressions occur in a level-1 solution; the statistics,
    # the energy quantile and the retention below are those of the allowed sub-pool
    filt = fibre_filters(st, rows, bordered) if fibre else [[] for _ in range(4)]
    if filt is None:
        print(f"N = {N} rows {rows}: a level-1 join is empty"); return []
    fid = [tuple((f[0], f[2].tobytes()) for f in filt[k]) for k in range(4)]
    mus, sigs, cache, samples, acc = [], [], {}, [], []
    for k, (kd, r) in enumerate(zip(KINDS, rows)):
        if (kd, r, fid[k]) not in cache:            # blocks with the same pool and filter share their statistics
            Y, a = fibre_sample(st, kd, r, K, rng, filt[k])
            V = st.values(kd, Y)
            cache[(kd, r, fid[k])] = (V.mean(0), np.cov(V.T), V, a)
        mu, sg, V, a = cache[(kd, r, fid[k])]
        mus.append(mu); sigs.append(sg); samples.append(V); acc.append(a)
    apools = [pools[k] * acc[k] for k in range(4)]
    Sig = sum(sigs)
    _, Sig_pinv = psd_pinv_sqrt(Sig)
    e = (ssxx.gauss(st, rows) if not bordered else bsz.gauss(st, rows))
    keep, plan = [], []
    for k in range(4):
        Lk = apools[k] if apools[k] <= L else L
        twin = [j for j in range(k) if (KINDS[j], rows[j], fid[j]) == (KINDS[k], rows[k], fid[k])
                and np.allclose(t - (sum(mus) - mus[j]), t - (sum(mus) - mus[k]))]
        if twin:                                    # same pool, same filter
            plan.append(plan[twin[0]]); continue
        Sm = Sig - sigs[k]
        U = unit_classes(st) if energy == "units" else []
        if U:
            Wsub, Spsub = psd_pinv_sqrt(Sm[np.ix_(U, U)])
            Wk = np.zeros((C, C)); Wk[np.ix_(U, U)] = Wsub
            Smp = np.zeros((C, C)); Smp[np.ix_(U, U)] = Spsub
        else:
            Wk, Smp = psd_pinv_sqrt(Sm)
        zk = t - (sum(mus) - mus[k])
        # energy quantile from prior samples (of the allowed sub-pool)
        V = samples[k]
        Ep = np.einsum("ij,jk,ik->i", V - zk, Smp, V - zk)
        frac = Lk / apools[k]
        tau = float(np.quantile(Ep, min(1.0, frac))) if frac < 1 else 1e30
        # retention: conditional law of the block's key in a solution
        ck = mus[k] + sigs[k] @ Sig_pinv @ (t - sum(mus))
        Ck = sigs[k] - sigs[k] @ Sig_pinv @ sigs[k]
        w, U = np.linalg.eigh((Ck + Ck.T) / 2); w = np.clip(w, 0, None)
        Xs = ck + (rng.standard_normal((K, C)) * np.sqrt(w)) @ U.T
        Es = np.einsum("ij,jk,ik->i", Xs - zk, Smp, Xs - zk)
        ret = float(np.mean(Es <= tau)) if frac < 1 else 1.0
        plan.append((Lk, tau, ret, Wk, zk))
    lret = sum(math.log10(max(p[2], 1e-12)) for p in plan)
    print(f"N = {N}{' bordered' if bordered else ''} H = <{','.join(map(str, gens))}> |H| = {len(st.H)} rows {rows}: "
          f"pools {[f'{p:.3g}' for p in pools]}, allowed {[f'{a:.3g}' for a in acc]} (fibres at "
          f"{sorted({f[0] for fk in filt for f in fk})}), logE {e:.2f}, kept {[f'{p[0]:.3g}' for p in plan]}, "
          f"retention {[f'{p[2]:.3g}' for p in plan]}, logE in lists {e + lret:.2f}", flush=True)
    if estimate_only:
        return e + lret
    lists, done = [], {}
    for k in range(4):
        Lk, tau, ret, Wk, zk = plan[k]
        key = (KINDS[k], rows[k], fid[k], Lk, tau, tuple(np.round(zk, 6)))
        if key in done:                             # same pool and same filter: the same kept list
            lists.append(done[key]); print(f"  block {'ABCD'[k]}: same list as before", flush=True); continue
        maxkeep = Lk if Lk < apools[k] else (int(min(L, 2 ** 31 - 1)) if filt[k] else 2 ** 31 - 1)
        masks, vals, en, sline = enum_block(st, KINDS[k], rows[k], zk, Wk, tau, maxkeep, "ABCD"[k], threads,
                                            filt[k])
        print(f"  block {'ABCD'[k]}: {len(masks)} kept; {sline} [{time.time() - t0:.0f}s]", flush=True)
        lists.append((masks, vals)); done[key] = (masks, vals)
    if min(len(x[0]) for x in lists) == 0:
        return []
    # kt4 input: order (A, C, B, D) so that the pairs are (A, C) and (B, D)
    order = [0, 2, 1, 3]
    keepdir = W.WORK / "lowe"; keepdir.mkdir(exist_ok=True)
    tagname = f"{N}{'b' if bordered else ''}_{'_'.join(map(str, gens))}_{'_'.join(map(str, rows))}" + \
        ("" if energy == "full" else f"_{energy}")
    inp = keepdir / f"kt4_{tagname}.bin"
    np.savez(keepdir / f"masks_{tagname}.npz", *[lists[k][0] for k in range(4)])
    with open(inp, "wb") as f:
        f.write(struct.pack("5i", *[len(lists[k][0]) for k in order], C))
        f.write(np.full(C, st.target, np.int16).tobytes())
        for k in order:
            f.write(np.ascontiguousarray(lists[k][1], np.int16).tobytes())
    # the join needs only the masks (to rebuild solutions) and subsamples of the class values (merge classes, split
    # law); dropping the full value arrays saves about 1 GB per kept list while kt4 runs
    lists = [(m, np.array(v[::max(1, len(v) // 200000)])) for m, v in lists]
    import gc; gc.collect()
    return join(st, N, gens, rows, bordered, lists, inp, nsig, nS, log, maxsol, t0, time_limit=time_limit,
                kt4_threads=kt4_threads)


def live_splits(lists, Sx, target, M=2 * 10 ** 7, seed=1):
    """split values that carry pairs on both sides: merge-class sums a + c of M random (A, C) pairs and t - b - d of M
    random (B, D) pairs, intersected, with their sample counts. Most lattice points near the centre carry no pair on
    one side (at N = 203 and 215, 99%); only these can hold a solution."""
    rng = np.random.default_rng(seed)
    base, off = 1 << 12, 2048
    def keys(L1, L2, sign, shift):
        us, cs = [], []
        for s0 in range(0, M, 2 * 10 ** 6):
            m = min(2 * 10 ** 6, M - s0)
            i = np.sort(rng.integers(0, len(L1), m)); j = np.sort(rng.integers(0, len(L2), m))
            rng.shuffle(j)
            P = shift + sign * (np.asarray(L1[i], np.int64)[:, Sx] + np.asarray(L2[j], np.int64)[:, Sx])
            k = np.zeros(m, np.int64)
            for c in range(len(Sx)):
                k = k * base + (P[:, c] + off)
            u, cnt = np.unique(k, return_counts=True)
            us.append(u); cs.append(cnt)
        u, inv = np.unique(np.concatenate(us), return_inverse=True)
        return u, np.bincount(inv, weights=np.concatenate(cs))
    tM = np.full(len(Sx), int(target), np.int64)
    uac, cac = keys(lists[0][1], lists[1][1], 1, 0)
    ubd, cbd = keys(lists[2][1], lists[3][1], -1, tM)
    common, ia, ib = np.intersect1d(uac, ubd, return_indices=True)
    P = np.zeros((len(common), len(Sx)), np.int64)
    k = common.copy()
    for c in range(len(Sx) - 1, -1, -1):
        P[:, c] = k % base - off; k //= base
    return P, cac[ia], cbd[ib]


def coset_basis(arrays, Sx, nvec=3000, seed=2):
    """reduced Hermite basis (rows, upper triangular in pivot order) of the lattice spanned by the differences of the
    values on the merge classes within each array (the lattice of the sums of one member of each array), and the sum of
    their first sampled members as base point"""
    from gauss_lat import hnf_basis
    rng = np.random.default_rng(seed)
    diffs, base = [], np.zeros(len(Sx), np.int64)
    for vals in arrays:
        X = np.asarray(vals[np.sort(rng.integers(0, len(vals), nvec))], np.int64)[:, Sx]
        diffs.append(X[1:] - X[0]); base += X[0]
    B, piv = hnf_basis(np.concatenate(diffs), len(Sx))
    B = [row if row[j] > 0 else [-x for x in row] for row, j in zip(B, piv)]
    for i in range(len(B)):                          # reduce each row modulo the pivots of the later rows
        for l in range(i + 1, len(B)):
            q = B[i][piv[l]] // B[l][piv[l]]
            if q:
                B[i] = [a - q * b for a, b in zip(B[i], B[l])]
    return B, piv, [int(x) for x in base]


def in_coset(x, B, piv, base):
    r = [int(a) - b for a, b in zip(x, base)]
    for row, j in zip(B, piv):
        if r[j] % row[j]:
            return False
        k = r[j] // row[j]
        if k:
            r = [a - k * b for a, b in zip(r, row)]
    return not any(r)


def coset_points(B, piv, base, m, half, cap):
    """the points of base + lattice in the box |x - m| <= half, at most cap (coordinates that are not pivots are
    fixed by the others when the lattice has lower rank, and checked at the end)"""
    X = np.array([base], np.int64)
    for row, j in zip(B, piv):
        g = row[j]
        lo = np.ceil((m[j] - half[j] - X[:, j]) / g).astype(np.int64)
        hi = np.floor((m[j] + half[j] - X[:, j]) / g).astype(np.int64)
        n = np.maximum(hi - lo + 1, 0)
        if n.sum() == 0 or n.sum() > cap:
            return np.zeros((0, len(m)), np.int64)
        rep = np.repeat(np.arange(len(X)), n)
        k = lo[rep] + (np.arange(int(n.sum())) - np.repeat(np.cumsum(n) - n, n))
        X = X[rep] + k[:, None] * np.asarray(row, np.int64)[None, :]
    return X[np.all(np.abs(X - np.asarray(m)) <= np.asarray(half), axis=1)]


def sigma_file(lists, Sx, target, nsig, path, seed=1):
    """split values for kt4, in order of the Gaussian law of a solution's split: with mu_k, Sigma_k the mean and
    covariance of each kept list (subsampled), a solution's (a + c) on the merge classes has mean
    m = mu_A + mu_C + (Sigma_A + Sigma_C) Sigma^+ (t - sum mu) and covariance V = (Sigma_A + Sigma_C)
    - (Sigma_A + Sigma_C) Sigma^+ (Sigma_A + Sigma_C). Only split values with pairs on both sides are used
    (live_splits); the nsig of them closest to m in the metric V^-1 are written."""
    mus, sigs = [], []
    for k in range(4):
        v = np.asarray(lists[k][1][::max(1, len(lists[k][1]) // 200000)], np.float64)
        mus.append(v.mean(0)); sigs.append(np.cov(v.T) if len(v) > 1 else np.zeros((v.shape[1], v.shape[1])))
    C = len(mus[0]); t = np.full(C, float(target))
    Sig = sum(sigs)
    w, U = np.linalg.eigh(Sig); inv = np.where(w > 1e-9 * w.max(), 1 / np.where(w > 0, w, 1), 0)
    Sp = (U * inv) @ U.T
    SAC = sigs[0] + sigs[1]                     # kt4 pairs its lists 0 and 1 (A and C)
    m = (mus[0] + mus[1] + SAC @ Sp @ (t - sum(mus)))[Sx]
    V = (SAC - SAC @ Sp @ SAC)[np.ix_(Sx, Sx)]
    Vi = np.linalg.inv(V + 1e-6 * np.eye(len(Sx)))
    P, nac, nbd = live_splits(lists, Sx, target, seed=seed)
    nlive = len(P)
    # Random pairs see only the bulk of each side. When the target lies off the centre on classes where the small
    # blocks carry the variance, a solution's split sits in the tails of both sides and is never sampled. So the
    # lattice points nearest the solution law are added too: the points of the (A + C) coset that also lie in the
    # (t - B - D) coset, in a box around m (kt4 skips a split value without pairs in about a millisecond).
    try:
        Ba, pa, ba = coset_basis([lists[0][1], lists[1][1]], Sx)
        Bb, pb, bsum = coset_basis([lists[2][1], lists[3][1]], Sx)
        bb = [int(target) - x for x in bsum]          # t - b - d = (t - b0 - d0) - (lattice); the lattice is symmetric
        sd = np.sqrt(np.clip(np.diag(V), 1, None))
        Q = np.zeros((0, len(Sx)), np.int64)
        for R in (1.0, 1.5, 2.0, 2.5, 3.0, 4.0):
            X = coset_points(Ba, pa, ba, m, R * sd, 4 * int(nsig) + 10 ** 5)
            if len(X) == 0 and R > 1.0:
                break
            Q = X
            if len(X) >= 2 * nsig:
                break
        if len(Q):
            Q = Q[[in_coset(x, Bb, pb, bb) for x in Q]]
        P = np.unique(np.concatenate([P, Q]), axis=0) if len(Q) else P
    except Exception as e:                          # (the lattice step is an addition; never lose the live values)
        print(f"  kt4: lattice split values skipped ({e})", flush=True)
    if len(P) == 0:
        np.savetxt(path, np.zeros((0, len(Sx)), np.int64), fmt="%d")
        return 0, 0.0, 0
    D = P - m
    d2 = np.einsum("ij,jk,ik->i", D, Vi, D)
    o = np.argsort(d2)[:int(nsig)]
    np.savetxt(path, P[o].astype(np.int64), fmt="%d")
    return len(o), float(np.sqrt(d2[o[-1]])), nlive


def join(st, N, gens, rows, bordered, lists, inp, nsig, nS, log, maxsol, t0, ordered=True, time_limit=None,
         kt4_threads=1):
    order = [0, 2, 1, 3]
    # merge classes: largest variance, from a subsample of every list (a full int64 copy would need gigabytes)
    var = sum(np.asarray(lists[k][1][::max(1, len(lists[k][1]) // 200000)], np.float64).var(0) for k in order)
    Sx = [int(c) for c in np.argsort(-var)[:nS]]
    sigarg = str(int(nsig))
    if ordered:
        sf = Path(str(inp) + ".sig")
        n, rad, nlive = sigma_file([lists[k] for k in order], Sx, st.target, nsig, sf)
        sigarg = "@" + str(sf)
        print(f"  kt4: merge classes {Sx}, {n} split values nearest the solution law ({nlive} sampled live, the rest "
              f"lattice points; Mahalanobis radius {rad:.2f})", flush=True)
        if n == 0:
            print(f"  kt4: no live split value [{time.time() - t0:.0f}s]", flush=True)
            return []
    else:
        print(f"  kt4: merge classes {Sx}, {int(nsig)} split values", flush=True)
    try:
        with open(str(inp) + ".progress", "w") as errf:            # one line per split value
            # one thread by default: each split value opens several parallel regions over ~1e5 pairs, and on a shared
            # machine their barriers cost more than the work (128 split values per second with one thread against 4
            # with four, N = 453, load 11 on 4 cores)
            env = dict(os.environ, OMP_NUM_THREADS=str(kt4_threads))
            p = subprocess.run([W.bin("kt4"), str(inp), ",".join(map(str, Sx)), sigarg, str(maxsol)],
                               stdout=subprocess.PIPE, stderr=errf, text=True, timeout=time_limit, env=env)
        out, rc = p.stdout, p.returncode
    except subprocess.TimeoutExpired as e:                          # the join's time budget is spent
        out = e.stdout.decode() if isinstance(e.stdout, bytes) else (e.stdout or "")
        rc = f"time limit {time_limit:.0f}s"
    tail = [l for l in out.splitlines() if l.startswith("TOTAL")]
    sols = [list(map(int, l.split()[1:5])) for l in out.splitlines() if l.startswith("SOL")]
    print(f"  kt4: {len(sols)} solutions; {tail[0] if tail else 'exit ' + str(rc)} [{time.time() - t0:.0f}s]",
          flush=True)
    found = []
    for s in sols:
        idx = [0] * 4
        for pos, k in enumerate(order):
            idx[k] = s[pos]
        seqs = []
        for k in range(4):
            m = int(lists[k][0][idx[k]])
            y = [1 if (m >> o) & 1 else -1 for o in range(len(st.sizes(KINDS[k])))]
            seqs.append(st.expand(KINDS[k], y))
        rec = {"N": N, "gens": list(gens), "subgroup_order": len(st.H), "rowsums": [int(x.sum()) for x in seqs],
               "method": "lowe_sz", "sequences": ["".join("+" if v > 0 else "-" for v in x) for x in seqs]}
        if bordered:
            Kr, Ki, b = bsz.build(seqs)
            rec["m"] = N - 1; rec["border"] = [[int(z.real), int(z.imag)] for z in b]
        else:
            ssxx.complex_hadamard(seqs); ssxx.goethals_seidel(seqs)
        print(f"  verified: complex Hadamard matrix of order {2 * N} (rows {rec['rowsums']})", flush=True)
        found.append(rec)
        if log:
            with open(log, "a") as f:
                f.write(json.dumps(rec) + "\n")
    return found


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("N", type=int)
    ap.add_argument("gens", type=int, nargs="+")
    ap.add_argument("--rows", required=True)
    ap.add_argument("--bordered", action="store_true")
    ap.add_argument("--L", type=float, default=2e7)
    ap.add_argument("--nsig", type=float, default=2e5)
    ap.add_argument("--S", type=int, default=5)
    ap.add_argument("--K", type=int, default=20000)
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--estimate-only", action="store_true")
    ap.add_argument("--resume", action="store_true", help="join the lists saved by an earlier run again")
    ap.add_argument("--log")
    ap.add_argument("--time", type=float, help="time budget of the kt4 join in seconds")
    ap.add_argument("--energy", choices=["full", "units"], default="full",
                    help="classes in the filter's energy: all, or only those prime to the group order")
    ap.add_argument("--no-precheck", action="store_true", help="skip the exact compression check")
    ap.add_argument("--no-fibre", action="store_true", help="do not restrict the blocks to allowed compressions")
    ap.add_argument("--kt4-threads", type=int, default=1, help="OpenMP threads of the kt4 join (default 1)")
    a = ap.parse_args()
    rows = tuple(int(x) for x in a.rows.split(","))
    if a.resume:
        resume(a.N, a.gens, rows, a.bordered, a.nsig, a.S, a.log, energy=a.energy, time_limit=a.time,
               kt4_threads=a.kt4_threads)
    else:
        run(a.N, a.gens, rows, a.bordered, a.L, a.nsig, a.S, a.K, a.threads, a.estimate_only, a.log, energy=a.energy,
            time_limit=a.time, precheck=not a.no_precheck, fibre=not a.no_fibre, kt4_threads=a.kt4_threads)


if __name__ == "__main__":
    main()
