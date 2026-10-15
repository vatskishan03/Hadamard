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
from bnw import Struct, subgroups

ROOT = Path(__file__).resolve().parent.parent


def found_complex():
    return {int(p.stem[1:]) // 2 for p in (ROOT / "complex-hadamard").glob("C*.json")}


def found_symmetric():
    return {int(p.name) // 4 for p in (ROOT / "symmetric").iterdir() if p.name.isdigit()}


def grouped_tuples(cp, N, rows, scheme):
    """feasible compression tuples, grouped by the compressed PAF vectors of the blocks: returns a list of
    (members per block, multiplicity) where every combination of members is feasible"""
    d = cp.d
    tgt = CE.comp_target(cp.st, d)
    groups = []
    for (kind, w), r in zip(scheme, rows):
        comps = cp.compressions(kind, r)
        for x, n in comps.items():
            cp.compressions_cache[(kind, x)] = n
        G = {}
        for x in comps:
            G.setdefault(tuple(w * CE.cpaf(cp.expand(kind, x))), []).append(x)
        groups.append(G)
    if any(not G for G in groups):
        return []
    keys = [np.array(list(G), np.int64).reshape(-1, d) for G in groups]
    names = [list(G) for G in groups]
    B = 1 << 22
    def enc(M):
        M = M + (B // 2)
        out = np.zeros(M.shape[:-1], np.int64)
        for j in range(d):
            out = out * B + M[..., j]
        return out
    h = len(scheme) // 2
    import itertools as _it
    # left half: all combinations of group keys of blocks 0..h-1
    left = {}
    for combo in _it.product(*[range(len(names[k])) for k in range(h)]):
        s_ = sum(keys[k][combo[k]] for k in range(h))
        left.setdefault(int(enc(s_[None, :])[0]), []).append(combo)
    lk = np.array(sorted(left), np.int64)
    out = []
    right_blocks = list(range(h, len(scheme)))
    if len(right_blocks) == 2:
        a_, b_ = right_blocks
        chunk = max(1, 2_000_000 // max(len(keys[b_]), 1))
        for i0 in range(0, len(keys[a_]), chunk):
            E = enc(tgt - keys[a_][i0:i0 + chunk, None, :] - keys[b_][None, :, :]).ravel()
            pos = np.searchsorted(lk, E)
            ok = (pos < len(lk)) & (lk[np.minimum(pos, len(lk) - 1)] == E)
            for f in np.flatnonzero(ok):
                i, j = divmod(int(f), len(keys[b_]))
                for combo in left[int(E[f])]:
                    out.append(tuple(combo) + (i0 + i, j))
    else:
        for i in range(len(keys[right_blocks[0]])):
            E = int(enc((tgt - keys[right_blocks[0]][i])[None, :])[0])
            for combo in left.get(E, []):
                out.append(tuple(combo) + (i,))
    res = []
    for combo in out:
        members = [groups[k][names[k][c]] for k, c in enumerate(combo)]
        mult = 1
        for mm in members:
            mult *= len(mm)
        res.append((members, mult))
    return res


def best_pattern(st, cp, N, pats, scheme, K, rng, maxpool, maxtup=150):
    """per pattern: the feasible compression groups; the estimate of a group is that of its first member tuple
    times the number of member tuples (members of a group differ by translations and reversals, which preserve
    the estimate); returns the best pattern"""
    best, nfeas = None, 0
    for rows in pats:
        gs = grouped_tuples(cp, N, rows, scheme)
        if not gs:
            continue
        nfeas += 1
        tot, bt = 0.0, None
        for members, mult in gs[:maxtup]:
            t = tuple(mm[0] for mm in members)
            pools = [cp.compressions_cache[(k, c)] for (k, w), c in zip(scheme, t)]
            if max(pools) > maxpool:
                continue
            lp, ln, k_, r_ = CE.cond_estimate(st, cp, t, K, rng, scheme=scheme)
            if lp == -math.inf:
                continue
            tot += mult * 10 ** (lp + ln)
            if bt is None or lp + ln > bt[0]:
                bt = (lp + ln, t, pools)
        if len(gs) > maxtup and tot > 0:
            tot *= len(gs) / maxtup
        if tot > 0 and (best is None or math.log10(tot) > best[0]):
            best = (math.log10(tot), tuple(int(r) for r in rows), bt)
    return best, nfeas


def best_pattern_old(st, cp, N, pats, scheme, K, rng, maxpool, maxtup=150):
    best, nfeas = None, 0
    for rows in pats:
        if scheme is CE.SZ:
            import comp_triage
            tups, lists = comp_triage.feasible_tuples_fast(cp, N, rows)
        else:
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


def space(cp):
    """compression space of an arbitrary block: prod over fibre classes of (orbits over the class + 1)"""
    out = 1
    for f in range(len(cp.fc["A"])):
        out *= int((cp.fid["A"] == f).sum()) + 1
    return out


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
    ap.add_argument("mode", choices=["bordered2q", "bordered", "plain", "propus"])
    ap.add_argument("N", type=int, nargs="*"); ap.add_argument("--nmax", type=int, default=1000)
    ap.add_argument("--omin", type=int, default=15); ap.add_argument("--omax", type=int, default=60)
    ap.add_argument("--K", type=int, default=400); ap.add_argument("--maxpool", type=float, default=6e10)
    ap.add_argument("--census", default=None, help="file with the census values v (propus mode)")
    ap.add_argument("--maxspace", type=float, default=5e6)
    a = ap.parse_args()
    rng = np.random.default_rng(1)
    if a.mode in ("bordered2q", "bordered", "plain"):
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
        elif a.mode == "bordered":
            m = N - 1
            for H, gens in subgroups(m).items():
                if len(H) == 1 or len(H) % 2 == 0 or (m - 1) in H:
                    continue
                st = Struct(N, gens, plain=False)
                O = len(st.sizesA)
                if not (a.omin <= O <= a.omax):
                    continue
                gt = m
                for h in H:
                    gt = math.gcd(gt, h - 1)
                # quotients on which H acts trivially, smallest first; every one is an exact filter
                ds = [d for d in range(2, gt + 1) if gt % d == 0]
                best_all, dead, used = None, None, []
                for d in ds:
                    cp = CE.Comp(st, d); cp.compressions_cache = {}
                    if space(cp) > a.maxspace:
                        continue
                    best, nfeas = best_pattern(st, cp, N, [r for _, r in bsz.patterns(st)], CE.SZ, a.K, rng,
                                               a.maxpool)
                    used.append(d)
                    if nfeas == 0:
                        dead = d; break
                    if best and (best_all is None or best[0] < best_all[0]):      # the most informative quotient
                        best_all = best
                head = f"N={N} bordered Z_{m} |H|={len(H)} gens={gens} O={O} quotients {used}"
                if dead:
                    print(f"{head}: no feasible pattern on Z_{dead}", flush=True)
                else:
                    report(head, best_all, -1)
        elif a.mode == "plain":
            for d, q, h, x, OA in dq_structures(N, a.omin, a.omax):
                st = Struct(N, (x,), plain=True)
                if (N - 1) in st.H:
                    continue
                cp = CE.Comp(st, d); cp.compressions_cache = {}
                if space(cp) > a.maxspace:
                    print(f"N={N} = {d} x {q} |H|={h}: compression space {space(cp):.2g} too large, skipped", flush=True)
                    continue
                best, nfeas = best_pattern(st, cp, N, ssxx.patterns(st), CE.SZ, a.K, rng, a.maxpool)
                report(f"N={N} = {d} x {q} |H|={h} g={x} O_A={OA}", best, nfeas)
        else:
            for d, q, h, x, OA in dq_structures(N, a.omin, a.omax):
                st = Struct(N, (x,), plain=True)
                if (N - 1) in st.H:
                    continue
                cp = CE.Comp(st, d); cp.compressions_cache = {}
                if space(cp) > a.maxspace:
                    print(f"v={N} = {d} x {q} |H|={h}: compression space {space(cp):.2g} too large, skipped", flush=True)
                    continue
                for name, scheme in (("symD", CE.PROPUS), ("freeD", CE.PROPUS_FREED)):
                    best, nfeas = best_pattern(st, cp, N, propus.patterns(N), scheme, a.K, rng, a.maxpool)
                    report(f"v={N} = {d} x {q} |H|={h} g={x} O_A={OA} {name}", best, nfeas)


if __name__ == "__main__":
    main()
