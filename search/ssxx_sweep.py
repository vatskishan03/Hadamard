#!/usr/bin/env python3
"""Sweep of complex (quaternary) Hadamard matrices of order 2N, N odd, through Szöllősi's (ssxx) construction
(search/ssxx.py), over the N for which we know no earlier construction.

An order 2N counts as known when Table 11.2 of Seberry-Yamada (1992) gives it (N < 1000,
search/data/complex_known_1992.txt) or one of these gives it:
- Williamson-type matrices of order N (construction (3) of arXiv:2605.08661; Turyn): Table A.1 of Seberry-Yamada
  (1992) as parsed in search/data/wt_known_1992.txt, good matrices of orders 35 and 127 (Đoković 1993), and the
  near-Williamson matrices of Kharaghani-Mohammadian-Tayfeh-Rezaie-Zaitsev (2026): every odd N <= 69 and
  73, 79, 91, 93, 95, 97, 103, 109, 127, 133, 181;
- a conference matrix of order 2N, C + iI: Paley (2N - 1 a prime power), Mathon (2N = q^2 (q + 2) + 1) and
  Seberry-Whiteman (2N = 5 * 9^(2t+1) + 1);
- a quaternary Legendre pair of length N - 1 (order 2(N-1) + 2): every even length <= 100 except 70, 76, 88, 92,
  94, 100 (Lebedev and Verma 2026), and the two families of Jedwab and Pender (2025): (q - 1)/2 for prime powers
  q = 1 (mod 4), and 2p for primes p with 2p - 1 a prime power;
- Szöllősi's order 94;
- Ito's dicyclic Hadamard groups with Egan's construction (R. Egan, Des. Codes Cryptogr. 87 (2019) 67-74, Theorem 4.5,
  after Schmidt 1999, Corollary 3.6): BH(2v, 4) for v = 2^a 10^b 26^c m with 2m - 1 or 4m - 1 a prime power or a
  Williamson matrix of order m. For odd v this adds every v with 4v - 1 a prime power (Paley's matrix; see
  complex-hadamard/PALEY.md);
- Williamson-type matrices of order N from the corollaries of Miyamoto's lemma (Seberry-Yamada 1992, Corollary 8.8;
  the M-structure paper, Corollaries 16, 18 and 20) with inputs known today: N = q = 5 (mod 8) a prime power with
  Williamson-type matrices of order (q - 1)/4, q = 1 (mod 8) with a Hadamard matrix of order (q - 1)/2, and
  N = 2q + 1 with q = 1 (mod 4), q + 1 a symmetric conference order and 2q - 1 a prime power. Djokovic
  (arXiv:1008.2043, Props 3.2-3.4) states 509, 857 and 883; 269 follows from the near-Williamson matrices of order 67
  (2026). We could not reconstruct the construction behind Corollary 8.8(i) (NOTES.md section 19) and follow the
  literature in counting it;
- V. Zaitsev's page "Quaternary Hadamard matrices" (https://www.cs.uleth.ca/~zaitsevv/quaternary.html, last updated
  14 Oct 2026; with Kharaghani, Mohammadian and Tayfeh-Rezaie): explicit matrices of 122 orders, checked here
  (search/data/complex_zaitsev_2026.txt). It covers the 11 orders found here by search, concurrently;
- Seberry-Yamada (1992) Corollaries 11.15 (orders n^j (n + 1)) and 11.17 (an explicit list, which includes 105),
  and the conference matrices of their codes c2 ((h - 1)^2 + 1, h skew) and c5 ((m - 1)^s + 1, m a conference order).
This is our reading of the literature, not a complete table (Table 11.2 of Seberry and Yamada's 2020 book, which
the 2026 papers cite, was not available to us).

    python3 search/ssxx_sweep.py list [--max 1000]
    python3 search/ssxx_sweep.py run N1 [N2 ...] [--maxo 40] [--maxpairs 2e11] [--log FILE]
"""
import argparse, contextlib, io, math, sys, time
from math import gcd
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import ssxx
from bnw import Struct, subgroups

DATA = Path(__file__).resolve().parent / "data" / "wt_known_1992.txt"
CDATA = Path(__file__).resolve().parent / "data" / "complex_known_1992.txt"
ZDATA = Path(__file__).resolve().parent / "data" / "complex_zaitsev_2026.txt"     # orders 2N
KMT = set(range(1, 70, 2)) | {73, 79, 91, 93, 95, 97, 103, 109, 127, 133, 181}
DJOKOVIC_1993 = {35, 127}
QLP_OPEN_TO_100 = {70, 76, 88, 92, 94, 100}
# Seberry-Yamada (1992), Corollary 11.17: complex Hadamard matrices of orders 2c from Williamson-type matrices
SY1992_COR_11_17 = {33, 39, 53, 73, 81, 83, 89, 93, 101, 105, 109, 113, 125, 137, 149, 153, 173, 189, 193, 197,
                    233, 241, 243, 257, 277, 281, 293}


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
    for q in range(3, int(round(nmax ** (1 / 3))) + 4):
        if q % 4 == 3 and prime_power(q):
            for p in (q + 2, q - 2):                       # both forms appear (1992 Table 11.2 uses q - 2)
                if p > 1 and prime_power(p) and q * q * p + 1 <= nmax:
                    out.add((q * q * p + 1) // 2)
    t = 0
    while 5 * 9 ** (2 * t + 1) + 1 <= nmax:
        out.add((5 * 9 ** (2 * t + 1) + 1) // 2); t += 1
    return out


def cor_11_15(N):
    """2N = n^j (n + 1) with n + 1 the order of a symmetric conference matrix (Paley: n = 1 mod 4 a prime power)."""
    for n in range(5, 2 * N, 4):
        if not prime_power(n):
            continue
        if (n + 1) % 2 == 0 and (2 * N) % (n + 1) == 0:
            r = 2 * N // (n + 1)
            while r % n == 0:
                r //= n
            if r == 1:
                return True
    return False


def conference_c2_c5(n):
    """Symmetric conference matrices of order n from the 1992 survey's codes c2 and c5:
    c2: n = (h - 1)^2 + 1 with h a skew Hadamard order (every 4k with k < 89 is one);
    c5: n = (m - 1)^s + 1, s >= 2, m the order of a conference matrix (Paley q + 1, Mathon, 46, c2)."""
    for k in range(1, 89):
        if (4 * k - 1) ** 2 + 1 == n:
            return True
    base = {46}
    for q in range(5, 2000, 4):
        if prime_power(q):
            base.add(q + 1)
    for k in range(1, 89):
        base.add((4 * k - 1) ** 2 + 1)
    for m in base:
        x = m - 1
        y = x * x
        while y + 1 <= n:
            if y + 1 == n:
                return True
            y *= x
    return False


def qlp_known(l):
    if l <= 100:
        return l not in QLP_OPEN_TO_100
    if prime_power(2 * l + 1) and (2 * l + 1) % 4 == 1:
        return True
    return l % 2 == 0 and is_prime(l // 2) and prime_power(l - 1)


def load(path):
    return {int(l) for l in path.read_text().splitlines() if l.strip().isdigit()}


def known_reason(N):
    wt = load(DATA)
    cx = load(CDATA)
    if N in cx:
        return "complex Hadamard (Seberry-Yamada 1992, Table 11.2)"
    if N in wt:
        return "WT (Seberry-Yamada 1992, Table A.1)"
    if N in SY1992_COR_11_17:
        return "complex Hadamard (Seberry-Yamada 1992, Cor. 11.17)"
    if cor_11_15(N):
        return "complex Hadamard n^j (n + 1) (Seberry-Yamada 1992, Cor. 11.15)"
    if N in DJOKOVIC_1993:
        return "good matrices (Đoković 1993)"
    if N in KMT:
        return "near-Williamson (KMTZ 2026)"
    if prime_power(2 * N - 1):
        return "Paley conference matrix"
    if N in mathon_orders(2 * N):
        return "Mathon conference matrix"
    if conference_c2_c5(2 * N):
        return "conference matrix (Seberry-Yamada 1992, codes c2/c5)"
    if qlp_known(N - 1):
        return "quaternary Legendre pair"
    if N == 47:
        return "Szöllősi 2026"
    if prime_power(4 * N - 1):
        return "Ito's dicyclic Hadamard groups (Paley's matrix), Egan 2019"
    r = miyamoto_wt(N, wt)
    if r:
        return r
    if 2 * N in load(ZDATA):
        return "Zaitsev et al., quaternary Hadamard matrices page (14 Oct 2026), explicit matrix"
    return None


def miyamoto_wt(N, wt):
    """Williamson-type matrices of order N from the corollaries of Miyamoto's lemma with inputs known today"""
    W = set(wt) | KMT | DJOKOVIC_1993 | {33}
    pk = prime_power(N)
    if pk and N % 8 == 5 and (N - 1) // 4 in W:
        return f"Williamson-type, Miyamoto's corollary from Williamson-type of order {(N - 1) // 4}"
    if pk and N % 8 == 1:                          # a Hadamard matrix of order (N - 1)/2 = 0 (mod 4) exists for N < 1000
        return f"Williamson-type, Miyamoto's corollary from H({(N - 1) // 2})"
    q = (N - 1) // 2
    sym_conf = (q + 1) % 4 == 2 and (prime_power(q) or (q + 1) // 2 in mathon_orders(q + 1))
    if q % 4 == 1 and sym_conf and prime_power(2 * q - 1):
        return f"Williamson-type, Seberry-Yamada Corollary 18 with q = {q}"
    return None


def candidates(nmax):
    return [N for N in range(3, nmax, 2) if known_reason(N) is None]


def plan(N, maxo=40, maxpairs=2e11, maxenum=2 ** 31):
    """Every (H, pattern) with -1 not in H, in order of the pair work of its exact join (cheapest first). The
    Gaussian model is not used: for composite N it underestimates product structures by orders of magnitude
    (154: 10^-1.7 predicted, 2352 found)."""
    out = []
    for H, gens in subgroups(N).items():
        if (N - 1) in H:
            continue
        if sum(gcd(h - 1, N) for h in H) // len(H) > maxo:   # orbit count by Burnside: |Fix(h)| = gcd(h - 1, N)
            continue
        st = Struct(N, gens, plain=True)
        for pat in ssxx.patterns(st):
            sizes = [st.counts(k).get(r, 0) for k, r in zip(ssxx.KINDS, pat)]
            if min(sizes) == 0:
                continue
            cost = ssxx.pair_cost(sizes, ssxx.best_pairing(sizes))
            # the Gray-code enumeration of a block costs 2^O steps whatever its list size: cap it as well
            enum = sum(2.0 ** len(st.sizes(k)) for k in ssxx.KINDS)
            if cost <= maxpairs and enum <= maxenum:
                out.append((cost + enum, len(H), gens, pat))
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
        try:
            with contextlib.redirect_stdout(buf):
                found = ssxx.search(N, gens, K=2, log=log, pats=pats, maxpairs=maxpairs, stop=True,
                                    deadline=t0 + budget_s)
        except Exception as ex:
            print(f"  H = <{','.join(map(str, gens))}>: error {ex!r}", flush=True)
            continue
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
