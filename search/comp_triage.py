#!/usr/bin/env python3
"""Compression-conditioned triage of Szöllősi structures for the open complex orders 2N (search/comp_est.py).

For every open N and every multiplier group H <= Z_N* with -1 not in H, every row-sum pattern is tested against
the compressions to every quotient Z_d (d | N) with a compression space of at most --maxspace vectors per block: a pattern with no feasible
compression tuple for some d has no solution (exact). For the others the conditioned estimate (summed over the
feasible tuples) is computed for each such d, and the smallest is reported. Prints the structures with conditioned logE >= --min and the pools of the best
tuple (the lists lowe_comp.py would enumerate).

    python3 search/comp_triage.py [N ...] [--nmax 700] [--min -1] [--omin 12] [--omax 72] [--K 800]
"""
import argparse, math, sys, time
from math import gcd, log10
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import ssxx, ssxx_sweep
import comp_est as CE
from bnw import Struct, subgroups

KINDS = CE.KINDS
HERE = Path(__file__).resolve().parent


def found_here():
    out = set()
    for p in (HERE.parent / "complex-hadamard").glob("C*.json"):
        try:
            out.add(int(p.stem[1:]) // 2)
        except ValueError:
            pass
    return out


def feasible_tuples_fast(cp, N, rows):
    """as comp_est.feasible_tuples: (A, B) sums sorted once, (C, D) sums looked up in chunks by binary search"""
    lists = [cp.compressions(k, r) for k, r in zip(KINDS, rows)]
    d = cp.d
    tgt = CE.comp_target(cp.st, d)
    keys = [list(L) for L in lists]
    if min(len(k) for k in keys) == 0:
        return [], lists
    P = [np.array([CE.cpaf(cp.expand(k, v)) for v in ks], np.int64).reshape(-1, d) for k, ks in zip(KINDS, keys)]
    B = 1 << 21
    def enc(M):                                   # M (..., d) -> int64 keys
        M = M + (B // 2)
        out = np.zeros(M.shape[:-1], np.int64)
        for j in range(d):
            out = out * B + M[..., j]
        return out
    ab = enc(P[0][:, None, :] + P[1][None, :, :]).ravel()
    order = np.argsort(ab, kind="stable"); abs_ = ab[order]
    nB = len(keys[1]); nD = len(keys[3])
    out = []
    chunk = max(1, 4_000_000 // max(nD, 1))
    for c0 in range(0, len(keys[2]), chunk):
        E = enc(tgt - P[2][c0:c0 + chunk, None, :] - P[3][None, :, :])          # (chunk, nD)
        flat = E.ravel()
        lo = np.searchsorted(abs_, flat, "left"); hi = np.searchsorted(abs_, flat, "right")
        for f in np.flatnonzero(hi > lo):
            ci, di = divmod(int(f), nD)
            for t in range(lo[f], hi[f]):
                ia, ib = divmod(int(order[t]), nB)
                out.append((keys[0][ia], keys[1][ib], keys[2][c0 + ci], keys[3][di]))
    return out, lists


def divisors(N):
    return [d for d in range(2, N) if N % d == 0]


def triage_N(N, args, rng):
    """per (H, pattern): exact feasibility under every quotient Z_d with at most args.maxfc fibre classes, and the
    smallest of the conditioned estimates over these d"""
    rows_out = []
    m = N - 1 if args.bordered else N
    for H, gens in subgroups(m).items():
        if (m - 1) in H or len(H) == 1:
            continue
        st = Struct(N, gens, plain=not args.bordered)
        OA = len(st.sizesA)
        if not (args.omin <= OA <= args.omax):
            continue
        cps = []
        for d in divisors(m):
            cp = CE.Comp(st, d)
            # size of the compression space of an arbitrary block: prod over fibre classes of (orbits over it + 1)
            space = 1
            for f in range(len(cp.fc["A"])):
                space *= int((cp.fid["A"] == f).sum()) + 1
            if space <= args.maxspace:
                cp.compressions_cache = {}
                cps.append(cp)
        if not cps:
            continue
        gt = m
        for h in H:
            gt = gcd(gt, h - 1)
        trivial_used = gt == 1 or gt == m or any(cp.d == gt for cp in cps)
        if args.bordered:
            import bsz
            pats = [r for _, r in bsz.patterns(st)]
        else:
            pats = ssxx.patterns(st)
        for pat in pats:
            per_d, dead = [], None
            for cp in cps:
                tups, lists = feasible_tuples_fast(cp, N, pat)
                if not tups:
                    dead = cp.d; break
                per_d.append((cp, tups, lists))
            if dead:
                continue
            ests = []
            for cp, tups, lists in per_d:
                for kind, L in zip(KINDS, lists):
                    for v, n in L.items():
                        cp.compressions_cache[(kind, v)] = n
                tot, best = 0.0, None
                # only tuples whose pools a low-energy run can enumerate count (the others cannot be searched)
                stups = [t for t in tups if max(cp.compressions_cache[(k, v)] for k, v in zip(KINDS, t)) <= args.maxpool]
                for tup in stups[:args.maxtup]:
                    lp, ln, k_, rank = CE.cond_estimate(st, cp, tup, args.K, rng)
                    if lp == -math.inf:
                        continue
                    e = lp + ln
                    tot += 10 ** e
                    if best is None or e > best[0]:
                        best = (e, tup, [cp.compressions_cache[(k, v)] for k, v in zip(KINDS, tup)])
                if len(stups) > args.maxtup and tot > 0:
                    tot *= len(stups) / args.maxtup
                ests.append((log10(tot) if tot > 0 else -math.inf, cp.d, len(tups), best))
            le, d, nt, best = min(ests, key=lambda x: x[0])
            if le == -math.inf or best is None:
                continue
            rows_out.append((le, len(H), gens, d, OA, tuple(int(x) for x in pat), nt, best,
                             [(x[1], round(x[0], 2)) for x in ests] + ([] if trivial_used else [("unused trivial", gt)])))
    return rows_out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("N", type=int, nargs="*")
    ap.add_argument("--nmax", type=int, default=700); ap.add_argument("--min", type=float, default=-1.0)
    ap.add_argument("--omin", type=int, default=12); ap.add_argument("--omax", type=int, default=72)
    ap.add_argument("--K", type=int, default=800); ap.add_argument("--top", type=int, default=4)
    ap.add_argument("--maxspace", type=float, default=3e5, help="largest compression space per block for a quotient")
    ap.add_argument("--maxtup", type=int, default=400)
    ap.add_argument("--bordered", action="store_true", help="the bordered array on Z_(N-1) (search/bsz.py)")
    ap.add_argument("--maxpool", type=float, default=6e10, help="only tuples whose pools are at most this count")
    args = ap.parse_args()
    rng = np.random.default_rng(1)
    Ns = args.N or [N for N in ssxx_sweep.candidates(args.nmax) if N not in found_here()]
    for N in Ns:
        t0 = time.time()
        try:
            res = triage_N(N, args, rng)
        except Exception as ex:
            print(N, "error", repr(ex), flush=True); continue
        res.sort(key=lambda r: -r[0])
        shown = 0
        for le, h, gens, d, OA, pat, nt, best, ests in res:
            if le < args.min or shown >= args.top:
                break
            e, tup, pools = best
            print(f"{N}{' bordered' if args.bordered else ''} cond logE {le:.2f} |H| {h} gens {gens} d {d} O_A {OA} rows {pat} tuples {nt} best {tup} "
                  f"logE {e:.2f} pools {['%.2g' % p for p in pools]} per-d {ests}", flush=True)
            shown += 1
        if not shown:
            print(f"{N} none ({len(res)} feasible structure-patterns) [{time.time() - t0:.0f}s]", flush=True)


if __name__ == "__main__":
    main()
