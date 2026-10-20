#!/usr/bin/env python3
"""Low-energy k-tree (search/lowe_sz.py) on one compression class at a time (search/comp_est.py).

For a Szöllősi structure (N, H, row sums) and a divisor d of N, every solution has block compressions to Z_d that
satisfy the compressed equation exactly. Each feasible compression tuple is a separate search: the blocks are
listed with their compressions fixed (lowe_enum with fibre-class constraints: per fibre class, the signs of the
orbits over it have a fixed weighted sum), the energies and the retention use the statistics of blocks sampled
with these compressions, and the kept lists are joined by kt4 as in lowe_sz. Tuples are reduced by the symmetries
that act on compressions (translations of C or D by the H-fixed points, C(g) -> C(-g), D(g) -> D(-g), and the
multipliers of Z_N acting on all four blocks), and searched in order of their conditional estimate.

    python3 search/lowe_comp.py N g [g2 ...] --rows a,b,c,d --d 3 [--bordered] [--L 2e7] [--nsig 2e5] [--tuples k]
                                [--only i]

With --bordered the blocks live on Z_(N-1) and the target is -4 (search/bsz.py); hits are checked with bsz.build.
"""
import argparse, math, os, struct, subprocess, sys, time
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import workdir as W
import ssxx
import comp_est as CE
import lowe_sz as LZ
from bnw import Struct

KINDS = ssxx.KINDS


def canon_tuples(st, cp, tups):
    """orbit representatives of compression tuples under the symmetries acting on compressions"""
    m, d = st.m, cp.d
    fixed = sorted({g % d for g in range(m) if all((h * g - g) % m == 0 for h in st.H)})
    units = sorted({u % d for u in range(1, m) if math.gcd(u, m) == 1})
    def expand(k, v):
        return tuple(int(x) for x in cp.expand(KINDS[k], v))
    def images(t):
        xs = [expand(k, v) for k, v in enumerate(t)]
        for u in units:
            ys = [tuple(x[(u * a) % d] for a in range(d)) for x in xs]
            for f1 in fixed:
                for r1 in (1, -1):
                    for f2 in fixed:
                        for r2 in (1, -1):
                            c = tuple(ys[2][(r1 * a + f1) % d] for a in range(d))
                            dd = tuple(ys[3][(r2 * a + f2) % d] for a in range(d))
                            yield (ys[0], ys[1], c, dd)
    reps = {}
    for t in tups:
        key = min(images(t))
        reps.setdefault(key, []).append(t)
    out = []
    for key, members in reps.items():
        out.append((members[0], len(members)))
    return out


def enum_block_comp(st, cp, kind, r, vals, center, Wm, tau, maxkeep, tag, threads):
    d = W.WORK / "tmp"; d.mkdir(exist_ok=True)
    inp = d / f"lowec_{st.N}_{tag}_{os.getpid()}.in"; out = d / f"lowec_{st.N}_{tag}_{os.getpid()}.out"
    sz = st.sizes(kind); S = st.S[kind]
    with open(inp, "wb") as f:
        f.write(struct.pack("4i", len(sz), st.C, int(r), int(min(maxkeep, 2 ** 31 - 1))))
        f.write(np.asarray(sz, np.int32).tobytes()); f.write(np.asarray(S, np.int32).tobytes())
        f.write(np.asarray(center, np.float32).tobytes()); f.write(np.asarray(Wm, np.float32).tobytes())
        f.write(struct.pack("f", float(tau)))
        f.write(struct.pack("i", len(vals)))
        f.write(np.asarray(cp.fid[kind], np.int32).tobytes()); f.write(np.asarray(cp.w[kind], np.int32).tobytes())
        f.write(np.asarray(vals, np.int32).tobytes())
        f.write(struct.pack("i", 0))                                       # no fibre filters
    p = subprocess.run([W.bin("lowe_enum"), str(inp), str(out), str(threads)], capture_output=True, text=True,
                       check=True)
    inp.unlink()
    with open(out, "rb") as f:
        n = struct.unpack("q", f.read(8))[0]
        masks = np.frombuffer(f.read(8 * n), np.uint64)
        vals_ = np.frombuffer(f.read(2 * n * st.C), np.int16).reshape(n, st.C)
        en = np.frombuffer(f.read(4 * n), np.float32)
    out.unlink()
    return masks, vals_, en, p.stdout.strip()


def plan_tuple(st, cp, tup, L, K, rng):
    """conditioned statistics, energy thresholds and retention for one compression tuple"""
    C = st.C
    t = np.full(C, float(st.target))
    pools = [cp.compressions_cache[(KINDS[k], v)] for k, v in enumerate(tup)]
    mus, sigs, samples = [], [], []
    for k, v in enumerate(tup):
        Y = cp.sample(KINDS[k], v, K, rng)
        V = st.values(KINDS[k], Y)
        mus.append(V.mean(0)); sigs.append(np.cov(V.T)); samples.append(V)
    Sig = sum(sigs)
    _, Sig_pinv = LZ.psd_pinv_sqrt(Sig)
    plan = []
    for k in range(4):
        Lk = pools[k] if pools[k] <= L else L
        Sm = Sig - sigs[k]
        Wk, Smp = LZ.psd_pinv_sqrt(Sm)
        zk = t - (sum(mus) - mus[k])
        Y = cp.sample(KINDS[k], tup[k], K, rng); V = st.values(KINDS[k], Y)
        Ep = np.einsum("ij,jk,ik->i", V - zk, Smp, V - zk)
        frac = Lk / pools[k]
        tau = float(np.quantile(Ep, min(1.0, frac))) if frac < 1 else 1e30
        ck = mus[k] + sigs[k] @ Sig_pinv @ (t - sum(mus))
        Ck = sigs[k] - sigs[k] @ Sig_pinv @ sigs[k]
        w, U = np.linalg.eigh((Ck + Ck.T) / 2); w = np.clip(w, 0, None)
        Xs = ck + (rng.standard_normal((K, C)) * np.sqrt(w)) @ U.T
        Es = np.einsum("ij,jk,ik->i", Xs - zk, Smp, Xs - zk)
        ret = float(np.mean(Es <= tau)) if frac < 1 else 1.0
        plan.append((Lk, tau, ret, Wk, zk))
    return pools, plan


def run(N, gens, rows, d, L=2e7, nsig=2e5, nS=5, K=4000, threads=4, ntup=None, only=None, estimate_only=False,
        log=None, seed=1, maxsol=5, maxpool=6e10, bordered=False):
    t0 = time.time()
    st = Struct(N, gens, plain=not bordered)
    cp = CE.Comp(st, d)
    cp.compressions_cache = {}
    rng = np.random.default_rng(seed)
    tups, lists = CE.feasible_tuples(cp, N, rows)
    for kind, Lc in zip(KINDS, lists):
        for v, n in Lc.items():
            cp.compressions_cache[(kind, v)] = n
    reps = canon_tuples(st, cp, tups)
    scored = []
    for tup, mult in reps:
        lp, ln, k_, rank = CE.cond_estimate(st, cp, tup, 1500, rng)
        if lp == -math.inf:
            continue
        scored.append((lp + ln, tup, mult, ln, lp))
    scored.sort(key=lambda x: -x[0])
    tot = sum(mult * 10 ** e for e, _, mult, _, _ in scored)
    print(f"N = {N}{' bordered' if bordered else ''} H = <{','.join(map(str, gens))}> |H| = {len(st.H)} rows {rows}, "
          f"d = {d}: {len(tups)} feasible "
          f"compression tuples, {len(reps)} up to symmetry, {len(scored)} with positive estimate; conditioned logE "
          f"{math.log10(tot) if tot > 0 else -math.inf:.2f} [{time.time() - t0:.0f}s]", flush=True)
    found = []
    sel = list(enumerate(scored))
    if only is not None:
        sel = [sel[i] for i in only]
    elif ntup is not None:
        sel = sel[:ntup]
    for i, (e, tup, mult, ln, lp) in sel:
        pools, plan = plan_tuple(st, cp, tup, L, K, rng)
        lret = sum(math.log10(max(p[2], 1e-12)) for p in plan)
        print(f"tuple {i}: {tup} (x{mult}): pools {[f'{p:.3g}' for p in pools]}, logE {e:.2f}, kept "
              f"{[f'{p[0]:.3g}' for p in plan]}, retention {[f'{p[2]:.3g}' for p in plan]}, logE in lists "
              f"{e + lret:.2f}", flush=True)
        if estimate_only:
            continue
        if max(pools) > maxpool:
            print(f"  skipped: pool above {maxpool:.2g}", flush=True); continue
        lists_, done = [], {}
        for k in range(4):
            Lk, tau, ret, Wk, zk = plan[k]
            key = (KINDS[k], tup[k], Lk, tau, tuple(np.round(zk, 6)))
            if key in done:
                lists_.append(done[key]); print(f"  block {'ABCD'[k]}: same list as before", flush=True); continue
            masks, vals, en, sline = enum_block_comp(st, cp, KINDS[k], rows[k], tup[k], zk, Wk, tau,
                                                     Lk if Lk < pools[k] else 2 ** 31 - 1, "ABCD"[k], threads)
            print(f"  block {'ABCD'[k]}: {len(masks)} kept; {sline} [{time.time() - t0:.0f}s]", flush=True)
            lists_.append((masks, vals)); done[key] = (masks, vals)
        if min(len(x[0]) for x in lists_) == 0:
            continue
        order = [0, 2, 1, 3]
        keepdir = W.WORK / "lowe"; keepdir.mkdir(exist_ok=True)
        tagname = f"{N}{'b' if bordered else ''}c{d}_{'_'.join(map(str, gens))}_{'_'.join(map(str, rows))}_t{i}"
        inp = keepdir / f"kt4_{tagname}.bin"
        np.savez(keepdir / f"masks_{tagname}.npz", *[lists_[k][0] for k in range(4)])
        with open(inp, "wb") as f:
            f.write(struct.pack("5i", *[len(lists_[k][0]) for k in order], st.C))
            f.write(np.full(st.C, st.target, np.int16).tobytes())
            for k in order:
                f.write(np.ascontiguousarray(lists_[k][1], np.int16).tobytes())
        got = LZ.join(st, N, gens, rows, bordered, lists_, inp, nsig, nS, log, maxsol, t0)
        found += got
        if got:
            break
    return found


def run_propus(v, gens, rows, d, freeD=True, L=2e7, nsig=2e5, nS=5, K=4000, threads=4, only=None, ntup=None,
                estimate_only=False, seed=1, maxsol=5, maxpool=6e10):
    """symmetric Hadamard matrices of order 4v from propus triples (A symmetric, B twice, D symmetric or free) on one
    compression class at a time: kt4 joins (A + D) against (2B + 0), hits are checked as H = Hᵀ, H Hᵀ = 4v I and
    appended to $HADAMARD_WORK/sym_hits.jsonl (search/sym_certify.py writes the certificate)"""
    import propus, json
    t0 = time.time()
    scheme = CE.PROPUS_FREED if freeD else CE.PROPUS
    st = Struct(v, gens, plain=True)
    cp = CE.Comp(st, d)
    cp.compressions_cache = {}
    rng = np.random.default_rng(seed)
    tups, lists = CE.feasible_tuples(cp, v, rows, scheme)
    for (kind, w), Lc in zip(scheme, lists):
        for x, n in Lc.items():
            cp.compressions_cache[(kind, x)] = n
    scored = []
    for tup in tups:
        lp, ln, k_, rank = CE.cond_estimate(st, cp, tup, 1500, rng, scheme=scheme)
        if lp == -math.inf:
            continue
        pools = [cp.compressions_cache[(k, x)] for (k, w), x in zip(scheme, tup)]
        if max(pools) <= maxpool:
            scored.append((lp + ln, tup, pools))
    scored.sort(key=lambda x: -x[0])
    print(f"v = {v} H = <{','.join(map(str, gens))}> |H| = {len(st.H)} propus{' freeD' if freeD else ''} rows {rows}, "
          f"d = {d}: {len(tups)} feasible tuples, {len(scored)} searchable with positive estimate [{time.time() - t0:.0f}s]",
          flush=True)
    sel = list(enumerate(scored))
    if only is not None:
        sel = [sel[i] for i in only]
    elif ntup is not None:
        sel = sel[:ntup]
    t = np.zeros(st.C)
    for i, (e, tup, pools) in sel:
        # statistics per block, energies relative to the weighted sum of the others
        mus, sigs = [], []
        for (kind, w), x in zip(scheme, tup):
            V = CE.values(st, kind, cp.sample(kind, x, K, rng))
            mus.append(w * V.mean(0)); sigs.append(w * w * np.cov(V.T))
        plan = []
        for k, ((kind, w), x) in enumerate(zip(scheme, tup)):
            Lk = pools[k] if pools[k] <= L else L
            Sm = sum(sigs) - sigs[k]
            Wk, Smp = LZ.psd_pinv_sqrt(Sm)
            zk = (t - (sum(mus) - mus[k])) / w                  # centre in the block's own values
            V = CE.values(st, kind, cp.sample(kind, x, K, rng))
            Ep = np.einsum("ij,jk,ik->i", w * (V - zk), Smp, w * (V - zk))
            frac = Lk / pools[k]
            tau = float(np.quantile(Ep, min(1.0, frac))) if frac < 1 else 1e30
            plan.append((Lk, tau, Wk * w, zk))
        print(f"tuple {i}: {tup}: pools {[f'{p:.3g}' for p in pools]}, logE {e:.2f}, kept "
              f"{[f'{p[0]:.3g}' for p in plan]}", flush=True)
        if estimate_only:
            continue
        blocks = []
        for k, ((kind, w), x) in enumerate(zip(scheme, tup)):
            Lk, tau, Wk, zk = plan[k]
            masks, vals, en, sline = enum_block_comp(st, cp, kind, rows[k], x, zk, Wk, tau,
                                                     Lk if Lk < pools[k] else 2 ** 31 - 1, "ABD"[k], threads)
            print(f"  block {'ABD'[k]}: {len(masks)} kept; {sline} [{time.time() - t0:.0f}s]", flush=True)
            blocks.append((masks, vals))
        if min(len(b[0]) for b in blocks) == 0:
            continue
        (mA, vA), (mB, vB), (mD, vD) = blocks
        zero = (np.zeros(1, np.uint64), np.zeros((1, st.C), np.int16))
        lists_ = [(mA, vA), (mB, (2 * vB.astype(np.int32)).astype(np.int16)), (mD, vD), zero]
        order = [0, 2, 1, 3]
        keepdir = W.WORK / "lowe"; keepdir.mkdir(exist_ok=True)
        tagname = f"p{v}{'f' if freeD else ''}c{d}_{'_'.join(map(str, gens))}_{'_'.join(map(str, rows))}_t{i}"
        inp = keepdir / f"kt4_{tagname}.bin"
        with open(inp, "wb") as f:
            f.write(struct.pack("5i", *[len(lists_[k][0]) for k in order], st.C))
            f.write(np.zeros(st.C, np.int16).tobytes())
            for k in order:
                f.write(np.ascontiguousarray(lists_[k][1], np.int16).tobytes())
        var = sum(np.asarray(lists_[k][1][::max(1, len(lists_[k][1]) // 200000)], np.float64).var(0) for k in order)
        Sx = [int(c) for c in np.argsort(-var)[:nS]]
        sf = Path(str(inp) + ".sig")
        n_, rad, _ = LZ.sigma_file([lists_[k] for k in order], Sx, 0, nsig, sf)
        print(f"  kt4: merge classes {Sx}, {n_} split values (radius {rad:.2f})", flush=True)
        with open(str(inp) + ".progress", "w") as errf:
            p = subprocess.run([W.bin("kt4"), str(inp), ",".join(map(str, Sx)), "@" + str(sf), str(maxsol)],
                               stdout=subprocess.PIPE, stderr=errf, text=True)
        sols = [list(map(int, l.split()[1:5])) for l in p.stdout.splitlines() if l.startswith("SOL")]
        tail = [l for l in p.stdout.splitlines() if l.startswith("TOTAL")]
        print(f"  kt4: {len(sols)} solutions; {tail[0] if tail else 'exit ' + str(p.returncode)} [{time.time() - t0:.0f}s]",
              flush=True)
        found = []
        for s_ in sols:
            ia, idd, ib = s_[0], s_[1], s_[2]                    # kt4 order: A, D, 2B, 0
            seqs = []
            for (kind, w), (masks, _), ix in zip(scheme, blocks, (ia, ib, idd)):
                mm = int(masks[ix])
                y = [1 if (mm >> o) & 1 else -1 for o in range(len(st.sizes(kind)))]
                seqs.append(st.expand(kind, y))
            A, B, D = seqs
            Hm = propus.symmetric_matrix_v3(A, B, D) if freeD else propus.symmetric_matrix(A, B, D)
            ok = np.array_equal(Hm, Hm.T) and np.array_equal(Hm @ Hm.T, 4 * v * np.eye(4 * v, dtype=Hm.dtype))
            print(f"  hit: symmetric Hadamard of order {4 * v}: {ok}", flush=True)
            if ok:
                rec = dict(n=v, gens=list(gens), symB=False, freeD=freeD, plain=False,
                           rows=[int(A.sum()), int(B.sum()), int(D.sum())],
                           seqs=[''.join('+' if x > 0 else '-' for x in q_) for q_ in (A, B, D)], method="lowe_comp")
                with open(W.hits('sym_hits.jsonl'), 'a') as fh:
                    fh.write(json.dumps(rec) + '\n')
                found.append(rec)
        if found:
            return found
    return []


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("N", type=int); ap.add_argument("gens", type=int, nargs="+")
    ap.add_argument("--rows", required=True); ap.add_argument("--d", type=int, required=True)
    ap.add_argument("--L", type=float, default=2e7); ap.add_argument("--nsig", type=float, default=2e5)
    ap.add_argument("--S", type=int, default=5); ap.add_argument("--K", type=int, default=4000)
    ap.add_argument("--threads", type=int, default=4); ap.add_argument("--tuples", type=int)
    ap.add_argument("--only", type=int, nargs="+"); ap.add_argument("--estimate-only", action="store_true")
    ap.add_argument("--maxpool", type=float, default=6e10); ap.add_argument("--log")
    ap.add_argument("--bordered", action="store_true")
    ap.add_argument("--propus", action="store_true", help="symmetric Hadamard 4N: propus (A sym, B twice, D)")
    ap.add_argument("--symD", action="store_true", help="with --propus: D symmetric too")
    a = ap.parse_args()
    rows = tuple(int(x) for x in a.rows.split(","))
    if a.propus:
        run_propus(a.N, a.gens, rows, a.d, not a.symD, a.L, a.nsig, a.S, a.K, a.threads, a.only, a.tuples,
                   a.estimate_only, maxpool=a.maxpool)
        return
    run(a.N, a.gens, rows, a.d, a.L, a.nsig, a.S, a.K, a.threads, a.tuples, a.only, a.estimate_only, a.log,
        maxpool=a.maxpool, bordered=a.bordered)


if __name__ == "__main__":
    main()
