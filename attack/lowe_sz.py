#!/usr/bin/env python3
"""Low-energy k-tree for Szöllősi structures (cyclic on Z_N, or bordered on Z_(N-1)) whose blocks are too large to
join exactly.

Each block k (A, B symmetric; C, D arbitrary; all constant on the orbits of a multiplier group H) is enumerated in full
by the depth-first search of attack/lowe_enum.c, and only its members of low energy are kept:

    E_k(x) = (x - z_k)^T Sigma_{-k}^+ (x - z_k),     z_k = t - sum_{j != k} mu_j,

where x is the member's vector of class values, mu_j, Sigma_j the mean and covariance of block j's class values over
its pool (from uniform samples), Sigma_{-k} = sum_{j != k} Sigma_j and t the target (0 cyclic, -4 bordered). Under the
Gaussian model, exp(-E_k(x)/2) is proportional to the density of the other three blocks' sum at t - x, so the kept
members are those that most likely complete to a solution. The kept lists are joined by attack/kt4.c (the four-list
k-tree of the low-energy skew searches) with the pairing (A, C), (B, D); every solution is rebuilt and checked exactly.

The estimate printed before the join is E (Gaussian, whole structure) times the retention of each filtered block,
the probability under the model that a solution's member passes the filter (Monte Carlo from the conditional law
N(c_k, Sigma_k - Sigma_k Sigma^+ Sigma_k)).

    python3 attack/lowe_sz.py N g [g2] --rows a,b,c,d [--bordered] [--L 2e7] [--nsig 2e5] [--S 5] [--estimate-only]
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


def resume(N, gens, rows, bordered=False, nsig=2e6, nS=5, log=None, maxsol=5):
    """rerun the join on the lists saved by an earlier run (WORK/lowe), with more split values"""
    st = Struct(N, gens, plain=not bordered)
    keepdir = W.WORK / "lowe"
    tagname = f"{N}{'b' if bordered else ''}_{'_'.join(map(str, gens))}_{'_'.join(map(str, rows))}"
    inp = keepdir / f"kt4_{tagname}.bin"
    z = np.load(keepdir / f"masks_{tagname}.npz")
    masks = [z[f"arr_{k}"] for k in range(4)]
    with open(inp, "rb") as f:
        hdr = struct.unpack("5i", f.read(20)); C = hdr[4]
        f.read(2 * C)
        vals = [np.frombuffer(f.read(2 * C * n), np.int16).reshape(n, C) for n in hdr[:4]]
    order = [0, 2, 1, 3]
    lists = [None] * 4
    for pos, k in enumerate(order):
        lists[k] = (masks[k], vals[pos])
    return join(st, N, gens, rows, bordered, lists, inp, nsig, nS, log, maxsol, time.time())


def run(N, gens, rows, bordered=False, L=2e7, nsig=2e5, nS=5, K=20000, threads=4, estimate_only=False, log=None,
        seed=1, maxsol=5):
    t0 = time.time()
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
    tagname = f"{N}{'b' if bordered else ''}_{'_'.join(map(str, gens))}_{'_'.join(map(str, rows))}"
    inp = keepdir / f"kt4_{tagname}.bin"
    np.savez(keepdir / f"masks_{tagname}.npz", *[lists[k][0] for k in range(4)])
    with open(inp, "wb") as f:
        f.write(struct.pack("5i", *[len(lists[k][0]) for k in order], C))
        f.write(np.full(C, st.target, np.int16).tobytes())
        for k in order:
            f.write(np.ascontiguousarray(lists[k][1], np.int16).tobytes())
    return join(st, N, gens, rows, bordered, lists, inp, nsig, nS, log, maxsol, t0)


def join(st, N, gens, rows, bordered, lists, inp, nsig, nS, log, maxsol, t0):
    order = [0, 2, 1, 3]
    allv = np.concatenate([lists[k][1] for k in order]).astype(np.int64)
    Sx = [int(c) for c in np.argsort(-allv.var(0))[:nS]]
    del allv
    print(f"  kt4: merge classes {Sx}, {int(nsig)} split values", flush=True)
    p = subprocess.run([W.bin("kt4"), str(inp), ",".join(map(str, Sx)), str(int(nsig)), str(maxsol)],
                       capture_output=True, text=True)
    tail = [l for l in p.stdout.splitlines() if l.startswith("TOTAL")]
    sols = [list(map(int, l.split()[1:5])) for l in p.stdout.splitlines() if l.startswith("SOL")]
    print(f"  kt4: {len(sols)} solutions; {tail[0] if tail else p.stderr.strip()[-300:]} [{time.time() - t0:.0f}s]",
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
    a = ap.parse_args()
    rows = tuple(int(x) for x in a.rows.split(","))
    if a.resume:
        resume(a.N, a.gens, rows, a.bordered, a.nsig, a.S, a.log)
    else:
        run(a.N, a.gens, rows, a.bordered, a.L, a.nsig, a.S, a.K, a.threads, a.estimate_only, a.log)


if __name__ == "__main__":
    main()
