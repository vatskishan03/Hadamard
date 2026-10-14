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
(search/compress_check.py); an empty join rules the structure/pattern out (--no-precheck skips this).

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


def enum_block(st, kind, r, center, Wm, tau, maxkeep, tag, threads):
    """run lowe_enum; returns (masks, values int16 [n, C], energies)"""
    d = W.WORK / "tmp"; d.mkdir(exist_ok=True)
    inp = d / f"lowe_{st.N}_{tag}_{os.getpid()}.in"; out = d / f"lowe_{st.N}_{tag}_{os.getpid()}.out"
    sz = st.sizes(kind); S = st.S[kind]
    with open(inp, "wb") as f:
        f.write(struct.pack("4i", len(sz), st.C, int(r), int(min(maxkeep, 2 ** 31 - 1))))
        f.write(np.asarray(sz, np.int32).tobytes()); f.write(np.asarray(S, np.int32).tobytes())
        f.write(np.asarray(center, np.float32).tobytes()); f.write(np.asarray(Wm, np.float32).tobytes())
        f.write(struct.pack("f", float(tau)))
    p = subprocess.run([W.bin("lowe_enum"), str(inp), str(out), str(threads)], capture_output=True, text=True, check=True)
    inp.unlink()
    with open(out, "rb") as f:
        n = struct.unpack("q", f.read(8))[0]
        masks = np.frombuffer(f.read(8 * n), np.uint64)
        vals = np.frombuffer(f.read(2 * n * st.C), np.int16).reshape(n, st.C)
        en = np.frombuffer(f.read(4 * n), np.float32)
    out.unlink()
    return masks, vals, en, p.stdout.strip()


def resume(N, gens, rows, bordered=False, nsig=2e6, nS=5, log=None, maxsol=5, energy="full", time_limit=None):
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
    return join(st, N, gens, rows, bordered, lists, inp, nsig, nS, log, maxsol, time.time(), time_limit=time_limit)


def unit_classes(st):
    """classes whose shifts are prime to the group order (the others are compression classes)"""
    return [c for c, s in enumerate(st.reps) if math.gcd(int(s), st.m) == 1]


def run(N, gens, rows, bordered=False, L=2e7, nsig=2e5, nS=5, K=20000, threads=4, estimate_only=False, log=None,
        seed=1, maxsol=5, energy="full", time_limit=None, precheck=True):
    """energy="units": the filter's energy uses only the unit classes. On exact solution sets of composite N the
    solution members sit at arithmetically fixed values on the compression classes, away from the Gaussian centre;
    without those classes they rank much higher at N = 119 and 267 (search/NOTES.md section 34)."""
    t0 = time.time()
    if precheck:
        import compress_check
        r = compress_check.check("bsz" if bordered else "ssxx", N, list(gens), tuple(rows), 200_000, 64, stop=True)
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
    mus, sigs, cache = [], [], {}
    for k, r in zip(KINDS, rows):
        if (k, r) not in cache:                     # blocks with the same pool share their statistics
            cache[(k, r)] = stats(st, k, r, K, rng)[:2]
        mu, sg = cache[(k, r)]
        mus.append(mu); sigs.append(sg)
    Sig = sum(sigs)
    _, Sig_pinv = psd_pinv_sqrt(Sig)
    e = (ssxx.gauss(st, rows) if not bordered else bsz.gauss(st, rows))
    keep, plan = [], []
    for k in range(4):
        Lk = pools[k] if pools[k] <= L else L
        twin = [j for j in range(k) if (KINDS[j], rows[j]) == (KINDS[k], rows[k])
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
        # energy quantile from prior samples
        Y = st.sample(KINDS[k], rows[k], K, rng); V = st.values(KINDS[k], Y)
        Ep = np.einsum("ij,jk,ik->i", V - zk, Smp, V - zk)
        frac = Lk / pools[k]
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
          f"pools {[f'{p:.3g}' for p in pools]}, logE {e:.2f}, kept {[f'{p[0]:.3g}' for p in plan]}, "
          f"retention {[f'{p[2]:.3g}' for p in plan]}, logE in lists {e + lret:.2f}", flush=True)
    if estimate_only:
        return e + lret
    lists, done = [], {}
    for k in range(4):
        Lk, tau, ret, Wk, zk = plan[k]
        key = (KINDS[k], rows[k], Lk, tau, tuple(np.round(zk, 6)))
        if key in done:                             # same pool and same filter: the same kept list
            lists.append(done[key]); print(f"  block {'ABCD'[k]}: same list as before", flush=True); continue
        masks, vals, en, sline = enum_block(st, KINDS[k], rows[k], zk, Wk, tau, Lk if Lk < pools[k] else 2 ** 31 - 1,
                                            "ABCD"[k], threads)
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
    return join(st, N, gens, rows, bordered, lists, inp, nsig, nS, log, maxsol, t0, time_limit=time_limit)


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
    if len(P) == 0:
        np.savetxt(path, np.zeros((0, len(Sx)), np.int64), fmt="%d")
        return 0, 0.0, 0
    D = P - m
    d2 = np.einsum("ij,jk,ik->i", D, Vi, D)
    o = np.argsort(d2)[:int(nsig)]
    np.savetxt(path, P[o].astype(np.int64), fmt="%d")
    return len(o), float(np.sqrt(d2[o[-1]])), len(P)


def join(st, N, gens, rows, bordered, lists, inp, nsig, nS, log, maxsol, t0, ordered=True, time_limit=None):
    order = [0, 2, 1, 3]
    # merge classes: largest variance, from a subsample of every list (a full int64 copy would need gigabytes)
    var = sum(np.asarray(lists[k][1][::max(1, len(lists[k][1]) // 200000)], np.float64).var(0) for k in order)
    Sx = [int(c) for c in np.argsort(-var)[:nS]]
    sigarg = str(int(nsig))
    if ordered:
        sf = Path(str(inp) + ".sig")
        n, rad, nlive = sigma_file([lists[k] for k in order], Sx, st.target, nsig, sf)
        sigarg = "@" + str(sf)
        print(f"  kt4: merge classes {Sx}, {n} of {nlive} live split values, nearest the solution law "
              f"(Mahalanobis radius {rad:.2f})", flush=True)
        if n == 0:
            print(f"  kt4: no live split value [{time.time() - t0:.0f}s]", flush=True)
            return []
    else:
        print(f"  kt4: merge classes {Sx}, {int(nsig)} split values", flush=True)
    try:
        with open(str(inp) + ".progress", "w") as errf:            # one line per split value
            p = subprocess.run([W.bin("kt4"), str(inp), ",".join(map(str, Sx)), sigarg, str(maxsol)],
                               stdout=subprocess.PIPE, stderr=errf, text=True, timeout=time_limit)
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
    a = ap.parse_args()
    rows = tuple(int(x) for x in a.rows.split(","))
    if a.resume:
        resume(a.N, a.gens, rows, a.bordered, a.nsig, a.S, a.log, energy=a.energy, time_limit=a.time)
    else:
        run(a.N, a.gens, rows, a.bordered, a.L, a.nsig, a.S, a.K, a.threads, a.estimate_only, a.log, energy=a.energy,
            time_limit=a.time, precheck=not a.no_precheck)


if __name__ == "__main__":
    main()
