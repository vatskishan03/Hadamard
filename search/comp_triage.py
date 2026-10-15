#!/usr/bin/env python3
"""Compression-conditioned triage of Szöllősi structures for the open complex orders 2N (search/comp_est.py).

For every open N and every multiplier group H <= Z_N* with -1 not in H, every row-sum pattern is tested against
the compressions to every quotient Z_d (d | N) with at most --maxfc fibre classes: a pattern with no feasible
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
    """as comp_est.feasible_tuples, with the (C, D) side vectorised"""
    lists = [cp.compressions(k, r) for k, r in zip(KINDS, rows)]
    d = cp.d
    tgt = CE.comp_target(cp.st, d)
    keys = [list(L) for L in lists]
    P = [np.array([CE.cpaf(cp.expand(k, v)) for v in ks], np.int64).reshape(-1, d) for k, ks in zip(KINDS, keys)]
    if min(len(k) for k in keys) == 0:
        return [], lists
    B = 1 << 21
    def enc(M):
        M = M + (B // 2)
        out = np.zeros(len(M), np.int64)
        for j in range(d):
            out = out * B + M[:, j]
        return out
    ab = {}
    for i, a in enumerate(keys[0]):
        e = enc(P[0][i] + P[1])
        for j, x in enumerate(e):
            ab.setdefault(int(x), []).append((a, keys[1][j]))
    abk = np.array(sorted(ab), np.int64)
    out = []
    for i, c in enumerate(keys[2]):
        e = enc(tgt - P[2][i] - P[3])
        hit = np.flatnonzero(np.isin(e, abk))
        for j in hit:
            for a, b in ab[int(e[j])]:
                out.append((a, b, c, keys[3][j]))
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
            if len(cp.fc["A"]) <= args.maxfc:
                cp.compressions_cache = {}
                cps.append(cp)
        if not cps:
            continue
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
                for tup in tups[:args.maxtup]:
                    lp, ln, k_, rank = CE.cond_estimate(st, cp, tup, args.K, rng)
                    if lp == -math.inf:
                        continue
                    e = lp + ln
                    tot += 10 ** e
                    if best is None or e > best[0]:
                        best = (e, tup, [cp.compressions_cache[(k, v)] for k, v in zip(KINDS, tup)])
                if len(tups) > args.maxtup and tot > 0:
                    tot *= len(tups) / args.maxtup
                ests.append((log10(tot) if tot > 0 else -math.inf, cp.d, len(tups), best))
            le, d, nt, best = min(ests, key=lambda x: x[0])
            if le == -math.inf:
                continue
            rows_out.append((le, len(H), gens, d, OA, tuple(int(x) for x in pat), nt, best,
                             [(x[1], round(x[0], 2)) for x in ests]))
    return rows_out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("N", type=int, nargs="*")
    ap.add_argument("--nmax", type=int, default=700); ap.add_argument("--min", type=float, default=-1.0)
    ap.add_argument("--omin", type=int, default=12); ap.add_argument("--omax", type=int, default=72)
    ap.add_argument("--K", type=int, default=800); ap.add_argument("--top", type=int, default=4)
    ap.add_argument("--maxfc", type=int, default=4); ap.add_argument("--maxtup", type=int, default=400)
    ap.add_argument("--bordered", action="store_true", help="the bordered array on Z_(N-1) (search/bsz.py)")
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
