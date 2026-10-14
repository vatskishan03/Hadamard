#!/usr/bin/env python3
"""Sweep of complex (quaternary) Hadamard matrices of order 2N, N odd, through Szöllősi's (ssxx) construction
(attack/ssxx.py), over the N for which we know no earlier construction.

An order 2N counts as known when one of these gives it:
- Williamson-type matrices of order N (construction (3) of arXiv:2605.08661; Turyn): Table A.1 of Seberry-Yamada
  (1992) as parsed in attack/data/wt_known_1992.txt, good matrices of orders 35 and 127 (Đoković 1993), and the
  near-Williamson matrices of Kharaghani-Mohammadian-Tayfeh-Rezaie-Zaitsev (2026): every odd N <= 69 and
  73, 79, 91, 93, 95, 97, 103, 109, 127, 133, 181;
- a conference matrix of order 2N, C + iI: Paley (2N - 1 a prime power), Mathon (2N = q^2 (q + 2) + 1) and
  Seberry-Whiteman (2N = 5 * 9^(2t+1) + 1);
- a quaternary Legendre pair of length N - 1 (order 2(N-1) + 2): every even length <= 100 except 70, 76, 88, 92,
  94, 100 (Lebedev and Verma 2026), and the two families of Jedwab and Pender (2025): (q - 1)/2 for prime powers
  q = 1 (mod 4), and 2p for primes p with 2p - 1 a prime power;
- Szöllősi's order 94.
This is our reading of the literature, not a complete table (Table 11.2 of Seberry and Yamada's 2020 book, which
the 2026 papers cite, was not available to us).

    python3 attack/ssxx_sweep.py list [--max 1000]
    python3 attack/ssxx_sweep.py run N1 [N2 ...] [--maxo 40] [--maxpairs 2e11] [--log FILE]
"""
import argparse, contextlib, io, math, sys, time
from math import gcd
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import ssxx
from bnw import Struct, subgroups

DATA = Path(__file__).resolve().parent / "data" / "wt_known_1992.txt"
KMT = set(range(1, 70, 2)) | {73, 79, 91, 93, 95, 97, 103, 109, 127, 133, 181}
DJOKOVIC_1993 = {35, 127}
QLP_OPEN_TO_100 = {70, 76, 88, 92, 94, 100}


def prime_power(n):
    if n < 2:
        return False
    for p in range(2, math.isqrt(n) + 1):
        if n % p == 0:
            while n % p == 0:
                n //= p
            return n == 1
    return True


def is_prime(n):
    return n > 1 and all(n % p for p in range(2, math.isqrt(n) + 1))


def mathon_orders(nmax):
    """N with a Mathon conference matrix of order 2N = q^2 (q + 2) + 1 (q = 3 mod 4 and q + 2 prime powers), and
    Seberry-Whiteman's 5 * 9^(2t+1) + 1."""
    out = set()
    for q in range(3, int(round(nmax ** (1 / 3))) + 2):
        if q % 4 == 3 and prime_power(q) and prime_power(q + 2):
            n = q * q * (q + 2) + 1
            if n <= nmax:
                out.add(n // 2)
    t = 0
    while 5 * 9 ** (2 * t + 1) + 1 <= nmax:
        out.add((5 * 9 ** (2 * t + 1) + 1) // 2); t += 1
    return out


def qlp_known(l):
    if l <= 100:
        return l not in QLP_OPEN_TO_100
    if prime_power(2 * l + 1) and (2 * l + 1) % 4 == 1:
        return True
    return l % 2 == 0 and is_prime(l // 2) and prime_power(l - 1)


def known_reason(N):
    wt = {int(x) for x in DATA.read_text().split() if x.isdigit()}
    if N in wt:
        return "WT (Seberry-Yamada 1992, Table A.1)"
    if N in DJOKOVIC_1993:
        return "good matrices (Đoković 1993)"
    if N in KMT:
        return "near-Williamson (KMTZ 2026)"
    if prime_power(2 * N - 1):
        return "Paley conference matrix"
    if N in mathon_orders(2 * N):
        return "Mathon conference matrix"
    if qlp_known(N - 1):
        return "quaternary Legendre pair"
    if N == 47:
        return "Szöllősi 2026"
    return None


def candidates(nmax):
    return [N for N in range(3, nmax, 2) if known_reason(N) is None]


def plan(N, maxo=40, maxpairs=2e11):
    """Every (H, pattern) with -1 not in H, in order of the pair work of its exact join (cheapest first). The
    Gaussian model is not used: for composite N it underestimates product structures by orders of magnitude
    (154: 10^-1.7 predicted, 2352 found)."""
    out = []
    for H, gens in subgroups(N).items():
        if (N - 1) in H:
            continue
        st = Struct(N, gens, plain=True)
        if len(st.sizesA) > maxo:
            continue
        for pat in ssxx.patterns(st):
            sizes = [st.counts(k).get(r, 0) for k, r in zip(ssxx.KINDS, pat)]
            if min(sizes) == 0:
                continue
            cost = ssxx.pair_cost(sizes, ssxx.best_pairing(sizes))
            if cost <= maxpairs:
                out.append((cost, len(H), gens, pat))
    out.sort(key=lambda r: (r[0], -r[1]))
    return out


def run_exhaustive(N, maxo=40, maxpairs=2e11, log=None, budget_s=3600):
    t0 = time.time()
    jobs = plan(N, maxo, maxpairs)
    print(f"N = {N}: {len(jobs)} structure/pattern pairs within {maxpairs:.0g} pair work", flush=True)
    groups = {}
    for cost, h, gens, pat in jobs:                     # one Struct per group, patterns in cost order
        groups.setdefault(gens, []).append(pat)
    for gens, pats in groups.items():
        if time.time() - t0 > budget_s:
            print(f"N = {N}: budget exhausted", flush=True)
            break
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            found = ssxx.search(N, gens, K=2, log=log, pats=pats, maxpairs=maxpairs, stop=True,
                                deadline=t0 + budget_s)
        if found:
            print(buf.getvalue(), end="", flush=True)
            return found
    return []


def run(N, maxo=40, maxpairs=2e11, log=None, budget_s=3600):
    """Structures in order of the Gaussian estimate; stop at the first verified matrix."""
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        rows = ssxx.triage(N, maxo)
    t0 = time.time()
    tried = set()
    for e, h, gens, OA, OL, C, pat, pools in rows:
        if e < -2.5 or time.time() - t0 > budget_s:
            break
        if gens in tried:
            continue
        tried.add(gens)
        st = Struct(N, gens, plain=True)
        pats = [p for p in ssxx.patterns(st) if ssxx.gauss(st, p) > -2.5]
        pats.sort(key=lambda p: -ssxx.gauss(st, p))
        print(f"N = {N}: H = <{','.join(map(str, gens))}> (order {h}), {len(pats)} patterns", flush=True)
        found = ssxx.search(N, gens, K=2, log=log, pats=pats, maxpairs=maxpairs, stop=True)
        if found:
            return found
    return []


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("list"); p.add_argument("--max", type=int, default=1000)
    p = sub.add_parser("run"); p.add_argument("N", type=int, nargs="+"); p.add_argument("--maxo", type=int, default=40)
    p.add_argument("--maxpairs", type=float, default=2e11); p.add_argument("--log")
    p.add_argument("--budget", type=float, default=3600)
    p.add_argument("--exhaustive", action="store_true", help="cost-ordered exact joins, no Gaussian ranking")
    a = ap.parse_args()
    if a.cmd == "list":
        c = candidates(a.max)
        print(len(c), "odd N with no complex Hadamard matrix of order 2N known to us:")
        print(" ".join(map(str, c)))
    else:
        for N in a.N:
            r = known_reason(N)
            if r:
                print(f"N = {N}: known ({r}), skipped"); continue
            found = run_exhaustive(N, a.maxo, a.maxpairs, a.log, a.budget) if a.exhaustive else \
                run(N, a.maxo, a.maxpairs, a.log, a.budget)
            print(f"N = {N}: {'FOUND' if found else 'nothing'}", flush=True)


if __name__ == "__main__":
    main()
