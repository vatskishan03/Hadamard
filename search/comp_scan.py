#!/usr/bin/env python3
"""Fast compression-conditioned scans of the product families (search/comp_est.py), one quotient per structure.

  bordered2q   complex orders 2N, bordered array on Z_(N-1) = Z_2 x Z_q (q prime), H <= Z_q* of odd order; quotient
               Z_2 (fibre sums = +-1 mod 2|H|)
  plain        complex orders 2N, N = d q (q prime, d in 3..15 coprime to q), H <= Z_q* of odd order acting trivially
               on Z_d (h = 1 mod d); quotient Z_d
  propus       symmetric Hadamard orders 4v of the Đoković–Seberry census, v = d q as above; the propus schemes
               PAF_A + 2 PAF_B + PAF_D (A symmetric, D symmetric or free); quotient Z_d

For every structure with --omin..--omax orbits and every row pattern: the exact feasible compression tuples (none:
the pattern has no solution), and the conditioned estimate summed over the tuples whose pools are at most --maxpool.
The best pattern is printed with its best tuple and pools. For propus the 2-adic densities of the weight-2 block are
not modelled, so its estimates rank structures but do not predict counts.

    python3 search/comp_scan.py bordered2q|plain|propus [N ...] [--nmax 1000] [--omin 15] [--omax 60]
"""
import argparse, math, sys
from pathlib import Path
import numpy as np
from sympy import isprime, primitive_root, divisors

sys.path.insert(0, str(Path(__file__).resolve().parent))
import comp_est as CE
import ssxx, ssxx_sweep, bsz, propus
from bnw import Struct

ROOT = Path(__file__).resolve().parent.parent


def found_complex():
    return {int(p.stem[1:]) // 2 for p in (ROOT / "complex-hadamard").glob("C*.json")}


def found_symmetric():
    return {int(p.name) // 4 for p in (ROOT / "symmetric").iterdir() if p.name.isdigit()}


def best_pattern(st, cp, N, pats, scheme, K, rng, maxpool, maxtup=150):
    best, nfeas = None, 0
    for rows in pats:
        tups, lists = CE.feasible_tuples(cp, N, rows, scheme)
        if not tups:
            continue
        nfeas += 1
        for (kind, w), L in zip(scheme, lists):
            for x, n in L.items():
                cp.compressions_cache[(kind, x)] = n
        tot, bt = 0.0, None
        for t in tups[:maxtup]:
            pools = [cp.compressions_cache[(k, c)] for (k, w), c in zip(scheme, t)]
            if max(pools) > maxpool:
                continue
            lp, ln, k_, r_ = CE.cond_estimate(st, cp, t, K, rng, scheme=scheme)
            if lp == -math.inf:
                continue
            tot += 10 ** (lp + ln)
            if bt is None or lp + ln > bt[0]:
                bt = (lp + ln, t, pools)
        if tot > 0 and (best is None or math.log10(tot) > best[0]):
            best = (math.log10(tot), tuple(int(r) for r in rows), bt)
    return best, nfeas


def report(head, best, nfeas):
    if best:
        e, rows, bt = best
        print(f"{head}: {nfeas} feasible patterns; logE {e:.2f} rows {rows}; best tuple {bt[1]} logE {bt[0]:.2f} pools "
              f"{['%.2g' % p for p in bt[2]]}", flush=True)
    else:
        print(f"{head}: {nfeas} feasible patterns, nothing searchable", flush=True)


def dq_structures(n, omin, omax, ds=(3, 5, 7, 9, 11, 13, 15)):
    """(d, q, |H|, generator) with n = d q, H <= Z_q* of odd order acting trivially on Z_d"""
    for d in ds:
        if n % d or not isprime(n // d) or math.gcd(d, n // d) != 1:
            continue
        q = n // d
        g = primitive_root(q)
        for h in divisors(q - 1):
            if h % 2 == 0 or h == 1:
                continue
            OA = d * (1 + (q - 1) // h)
            if not (omin <= OA <= omax):
                continue
            y = pow(g, (q - 1) // h, q)
            yield d, q, h, next(t for t in range(n) if t % d == 1 and t % q == y), OA


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=["bordered2q", "plain", "propus"])
    ap.add_argument("N", type=int, nargs="*"); ap.add_argument("--nmax", type=int, default=1000)
    ap.add_argument("--omin", type=int, default=15); ap.add_argument("--omax", type=int, default=60)
    ap.add_argument("--K", type=int, default=400); ap.add_argument("--maxpool", type=float, default=6e10)
    ap.add_argument("--census", default=None, help="file with the census values v (propus mode)")
    a = ap.parse_args()
    rng = np.random.default_rng(1)
    if a.mode in ("bordered2q", "plain"):
        Ns = a.N or [N for N in ssxx_sweep.candidates(a.nmax) if N not in found_complex()]
    else:
        Ns = a.N or [v for v in map(int, open(a.census).read().split()) if v < a.nmax and v not in found_symmetric()]
    for N in Ns:
        if a.mode == "bordered2q":
            m = N - 1
            if m % 2 or not isprime(m // 2):
                continue
            q = m // 2
            g = primitive_root(q)
            for h in divisors(q - 1):
                O = 2 * (1 + (q - 1) // h)
                if h % 2 == 0 or h == 1 or not (a.omin <= O <= a.omax):
                    continue
                y = pow(g, (q - 1) // h, q)
                x = y if y % 2 == 1 else y + q
                st = Struct(N, (x,), plain=False)
                cp = CE.Comp(st, 2); cp.compressions_cache = {}
                best, nfeas = best_pattern(st, cp, N, [r for _, r in bsz.patterns(st)], CE.SZ, a.K, rng, a.maxpool)
                report(f"N={N} bordered Z_2 x Z_{q} |H|={h} g={x} O={O}", best, nfeas)
        elif a.mode == "plain":
            for d, q, h, x, OA in dq_structures(N, a.omin, a.omax):
                st = Struct(N, (x,), plain=True)
                if (N - 1) in st.H:
                    continue
                cp = CE.Comp(st, d); cp.compressions_cache = {}
                best, nfeas = best_pattern(st, cp, N, ssxx.patterns(st), CE.SZ, a.K, rng, a.maxpool)
                report(f"N={N} = {d} x {q} |H|={h} g={x} O_A={OA}", best, nfeas)
        else:
            for d, q, h, x, OA in dq_structures(N, a.omin, a.omax):
                st = Struct(N, (x,), plain=True)
                if (N - 1) in st.H:
                    continue
                cp = CE.Comp(st, d); cp.compressions_cache = {}
                for name, scheme in (("symD", CE.PROPUS), ("freeD", CE.PROPUS_FREED)):
                    best, nfeas = best_pattern(st, cp, N, propus.patterns(N), scheme, a.K, rng, a.maxpool)
                    report(f"v={N} = {d} x {q} |H|={h} g={x} O_A={OA} {name}", best, nfeas)


if __name__ == "__main__":
    main()
