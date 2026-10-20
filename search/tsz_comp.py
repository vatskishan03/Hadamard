#!/usr/bin/env python3
"""Compression-conditioned triage of the twisted Szöllősi structures (search/tsz.py, complex-hadamard/TWISTED.md).

The twisted quadruples are (A, A∘μ, C, C∘μ) (full twist: A on the L-orbits, C on the H-orbits) and (A, B, C, C∘μ)
(half twist: A, B on the orbits of L' = <L, μ>). For a divisor d of m (m = N, or N - 1 for the bordered array) the
fibre sums x(a) = sum_{g = a (mod d)} X(g) satisfy sum_k PAF_{x_k} = the compressed target on Z_d (search/comp_est.py),
and the image X∘μ has the compression x∘μ̄ (μ̄ = μ mod d), whose PAF is s -> PAF_x(μ̄ s). So per independent block

    full:  A: P_xA(s) + P_xA(μ̄ s),  C: P_xC(s) + P_xC(μ̄ s)
    half:  A: P_xA,  B: P_xB,  C: P_xC(s) + P_xC(μ̄ s)

and a row pattern with no compression tuple solving the compressed equation has no solution. For the feasible tuples
the conditioned estimate is that of comp_est.cond_estimate: blocks sampled uniformly with their compressions, mean and
covariance of the twisted class values (tsz.TStruct.values), the lattice of the samples, the density at the target.
The Gaussian of tsz.py ignores both effects: it rates the bordered half twist of N = 203 (H of order 5) at 10^0.24,
while its Z_2 compressions have no solution.

    python3 search/tsz_comp.py N g [g2] --mu MU --kind half|full [--bordered]       triage one structure
    python3 search/tsz_comp.py scan [N ...] [--maxo 44] [--bordered]       every twisted structure of every open N
    python3 search/tsz_comp.py search N g --mu MU --kind half [--rows a,b,c] [--log FILE]
        exact search restricted to the feasible compressions: the P side is the explicit list of the A (and B) with
        compressions in a feasible group, and each C compression of the group is streamed by tsz.c with its fibre
        constraint (exhaustive for the pattern, since every solution has feasible compressions)
    python3 search/tsz_comp.py sweep [N ...] [--bordered] [--minE -1] [--maxstream 2e10] [--log FILE]
        triage, then search every pattern with conditioned estimate >= minE (or >= minE - 2 and a stream <= 1e9)
"""
import argparse, math, subprocess, sys
from math import gcd, log10
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import comp_est as CE
import tsz


def twin(p, mu, d):
    return p + p[(mu * np.arange(d)) % d]


class TwComp:
    """the compressions of the independent blocks of a twisted structure on Z_d"""

    def __init__(self, ts, d):
        assert ts.m % d == 0
        self.ts, self.d, self.mu = ts, d, ts.mu % d
        self.cpP = CE.Comp(ts.st2 if ts.kind == "half" else ts.st, d)
        self.cpC = CE.Comp(ts.st, d)
        self.cpP.compressions_cache, self.cpC.compressions_cache = {}, {}
        self.tgt = CE.comp_target(ts.st, d)
        # independent blocks: (comp, kind, side, twinned)
        if ts.kind == "full":
            self.blocks = [(self.cpP, "L", "P", True), (self.cpC, "A", "C", True)]
        else:
            self.blocks = [(self.cpP, "L", "P", False), (self.cpP, "L", "P", False), (self.cpC, "A", "C", True)]

    def space(self):
        """the largest compression space of a block: prod over fibre classes of the number of class values"""
        out = 1
        for cp, kind, _, _ in self.blocks:
            s = 1
            for cnt in cp.class_counts(kind):
                s *= len(cnt)
            out = max(out, s)
        return out

    def contrib(self, b, x):
        cp, kind, _, tw = self.blocks[b]
        p = CE.cpaf(cp.expand(kind, x))
        return twin(p, self.mu, self.d) if tw else p

    def groups(self, b, r):
        """compressions of block b with row sum r, grouped by their contribution to the compressed equation"""
        cp, kind, _, tw = self.blocks[b]
        cache = self.__dict__.setdefault("_groups", {})
        if (b, r) in cache:
            return cache[(b, r)]
        comps = cp.compressions(kind, r)
        for x, n in comps.items():
            cp.compressions_cache[(kind, x)] = n
        G = {}
        if comps:
            xs = list(comps)
            d = self.d
            fidx = np.zeros(d, np.int64)
            for f, O in enumerate(cp.fc[kind]):
                fidx[O] = f
            E = np.array(xs, np.int64)[:, fidx]                         # (n, d) compressions
            P = np.stack([(E * np.roll(E, -t, axis=1)).sum(1) for t in range(d)], 1)
            if tw:
                P = P + P[:, (self.mu * np.arange(d)) % d]
            U, inv = np.unique(P, axis=0, return_inverse=True)
            for i, x in zip(np.asarray(inv).ravel(), xs):
                G.setdefault(tuple(int(t) for t in U[i]), []).append(x)
        cache[(b, r)] = G
        return G

    def feasible(self, pat):
        """[(members per block, multiplicity)]: every combination of members solves the compressed equation"""
        Gs = [self.groups(b, r) for b, r in enumerate(pat)]
        if any(not G for G in Gs):
            return []
        keys = [list(G) for G in Gs]
        K = [np.array(k, np.int64).reshape(-1, self.d) for k in keys]
        last = {}
        for i, k in enumerate(keys[-1]):
            last.setdefault(k, []).append(i)
        out = []
        def emit(combo, j):
            members = [Gs[b][keys[b][i]] for b, i in enumerate(combo)] + [Gs[-1][keys[-1][j]]]
            mult = 1
            for mm in members:
                mult *= len(mm)
            out.append((members, mult))
        if len(K) == 2:
            R = self.tgt[None, :] - K[0]
            for i in range(len(R)):
                for j in last.get(tuple(int(t) for t in R[i]), []):
                    emit((i,), j)
        else:
            lastset = {k for k in last}
            for i in range(len(K[0])):
                R = self.tgt[None, :] - K[0][i][None, :] - K[1]          # (n_B, d)
                for j1 in range(len(R)):
                    k = tuple(int(t) for t in R[j1])
                    if k in lastset:
                        for j in last[k]:
                            emit((i, j1), j)
        return out

    def stats(self, b, x, K, rng):
        cp, kind, side, _ = self.blocks[b]
        cache = cp.__dict__.setdefault("tw_stats", {})
        key = (kind, side, x)
        if key not in cache:
            Y = cp.sample(kind, x, K, rng)
            V = self.ts.values(side, Y)
            Vi = np.rint(V).astype(np.int64)
            cache[key] = (cp.compressions_cache[(kind, x)], V.mean(0), np.cov(V.T), Vi[0], CE.diff_lattice(Vi, rng))
        return cache[key]

    def estimate(self, tup, K, rng):
        """(log10 P(class values = target | compressions), log10 of the number of block tuples)"""
        ts = self.ts
        CM = ts.CM
        mu = np.zeros(CM); S = np.zeros((CM, CM)); base = np.zeros(CM, np.int64); gens = []; ln = 0.0
        for b, x in enumerate(tup):
            cnt, m_, c_, b_, B_ = self.stats(b, x, K, rng)
            ln += log10(cnt); mu += m_; S += np.atleast_2d(c_); base += b_; gens += [list(r) for r in B_]
        B = CE.echelon(gens) if gens else []
        if B:
            G = [[sum(p * q for p, q in zip(r1, r2)) for r2 in B] for r1 in B]
            covol = math.exp(0.5 * math.log(CE.bareiss_det(G)))
        else:
            covol = 1.0
        t = np.full(CM, int(ts.target), np.int64)
        if not CE.in_lattice(B, t - base):
            return -math.inf, ln
        w, Uv = np.linalg.eigh(S)
        keep = w > 1e-7 * max(w.max(), 1e-12)
        delta = t - mu
        proj = Uv[:, keep] @ (Uv[:, keep].T @ delta)
        if np.linalg.norm(delta - proj) > 1e-6 * (1 + np.linalg.norm(delta)):
            return -math.inf, ln
        k_ = int(keep.sum())
        q = float(np.sum((Uv[:, keep].T @ delta) ** 2 / w[keep]))
        lp = log10(covol) - 0.5 * k_ * log10(2 * math.pi) - 0.5 * float(np.sum(np.log10(w[keep]))) - q / 2 / math.log(10)
        return lp, ln


def quotients(ts, maxspace=2e5, maxd=60):
    """divisors d of m on which H acts trivially (the fibre sums then lie on coarse progressions), smallest first"""
    m = ts.m
    g = m
    for h in ts.st.H:
        g = gcd(g, h - 1)
    out = []
    for d in range(2, min(g, maxd) + 1):
        if g % d == 0:
            tc = TwComp(ts, d)
            if tc.space() <= maxspace:
                out.append(tc)
    return out


def triage(ts, K=400, seed=1, maxpool=3e10, maxtup=200, quiet=False):
    """per pattern: the quotient that kills it, or the conditioned estimate on the most informative quotient"""
    rng = np.random.default_rng(seed)
    res = []
    tcs = quotients(ts)
    for pat in ts.patterns():
        g0 = ts.gauss(pat)
        dead, best = None, None
        for tc in tcs:
            fs = tc.feasible(pat)
            if not fs:
                dead = tc.d
                break
            tot, bt = 0.0, None
            for members, mult in fs[:maxtup]:
                tup = tuple(mm[0] for mm in members)
                pools = [tc.blocks[b][0].compressions_cache[(tc.blocks[b][1], x)] for b, x in enumerate(tup)]
                if max(pools) > maxpool:
                    continue
                lp, ln = tc.estimate(tup, K, rng)
                if lp == -math.inf:
                    continue
                tot += mult * 10 ** (lp + ln)
                if bt is None or lp + ln > bt[0]:
                    bt = (lp + ln, tup, pools)
            if len(fs) > maxtup and tot > 0:
                tot *= len(fs) / maxtup
            if ts.kind == "half" and pat[0] == pat[1]:
                tot /= 2                            # (A, B) and (B, A) are the same solution, as in tsz.gauss
            e = log10(tot) if tot > 0 else -math.inf
            if best is None or e < best[0]:
                best = (e, tc.d, len(fs), bt)
        res.append((pat, g0, dead, best))
        if not quiet:
            head = (f"N={ts.N}{'b' if ts.bordered else ''} {ts.kind} |H|={len(ts.st.H)} gens={ts.gens} mu={ts.mu} "
                    f"O_P={ts.OP} O_C={ts.OC} classes={ts.CM} rows {pat} sides {[f'{p:.3g}' for p in tsz.side_size(ts, pat)]}")
            if dead:
                print(f"{head}: dead on Z_{dead} (plain Gaussian {g0:.2f})", flush=True)
            elif best and best[3]:
                e, d, nf, bt = best
                print(f"{head}: Z_{d} {nf} feasible groups, conditioned logE {e:.2f} (plain {g0:.2f}); best tuple "
                      f"{bt[1]} logE {bt[0]:.2f} pools {['%.2g' % p for p in bt[2]]}", flush=True)
            else:
                print(f"{head}: no quotient (plain Gaussian {g0:.2f})", flush=True)
    return res


def comp_of(cp, kind, Y):
    """fibre-class values of sign vectors Y (rows) for the blocks of the given kind"""
    F = len(cp.fc[kind])
    M = np.zeros((Y.shape[1], F), np.int64)
    M[np.arange(Y.shape[1]), cp.fid[kind]] = cp.w[kind]
    return np.rint(Y @ M).astype(np.int64)


def search(ts, pat, threads=1, log=None, d=None, maxside=3e7, quiet=False):
    """exact search of one pattern restricted to its feasible compressions: the P side is the explicit list of the
    A (and B) whose compressions lie in a feasible group, and every C compression of the group is streamed by tsz.c
    with its fibre constraint. Returns the verified records."""
    import struct as _st, os, time, json, subprocess as sp
    import workdir as W
    tcs = [tc for tc in quotients(ts) if d is None or tc.d == d]
    best = None
    for tc in tcs:
        fs = tc.feasible(pat)
        if not fs:
            if not quiet:
                print(f"N={ts.N} rows {pat}: dead on Z_{tc.d}", flush=True)
            return []
        cost = 0
        cpC, kC = tc.blocks[-1][0], tc.blocks[-1][1]
        xs = {x for members, _ in fs for x in members[-1]}
        cost = sum(cpC.compressions_cache[(kC, x)] for x in xs)
        if best is None or cost < best[0]:
            best = (cost, tc, fs)
    if best is None:
        raise ValueError("no quotient")
    cost, tc, fs = best
    if not quiet:
        print(f"N={ts.N}{'b' if ts.bordered else ''} {ts.kind} gens={ts.gens} mu={ts.mu} rows {pat}: Z_{tc.d}, "
              f"{len(fs)} feasible groups, C stream {cost:.3g} (full pool {tsz.side_size(ts, pat)[1]:.3g})", flush=True)
    cpP, kP = tc.blocks[0][0], tc.blocks[0][1]
    nP = len(tc.blocks) - 1
    # translating C by an H-fixed point t (and D = C∘μ by μ^-1 t) maps solutions to solutions and shifts the
    # compression of C by t mod d: one C compression per translation class suffices (singleton fibre classes)
    cpC0, kC0 = tc.blocks[-1][0], tc.blocks[-1][1]
    shifts = [0]
    if all(len(c) == 1 for c in cpC0.fc[kC0]):
        fix = [t for t in range(ts.m) if all(((h - 1) * t) % ts.m == 0 for h in ts.st.H)]
        shifts = sorted({t % tc.d for t in fix})
    order = [c[0] for c in cpC0.fc[kC0]]
    def canon(x):
        if len(shifts) == 1:
            return tuple(x)
        full = np.zeros(tc.d, np.int64); full[order] = x
        return min(tuple(int(v) for v in np.roll(full, s)[order]) for s in shifts)
    # all P-side sign vectors with the needed row sums, their compressions and class values
    Pl = {}
    for r in set(pat[:nP]):
        masks, Y = tsz.all_signs(ts.sizesP, r)
        Pl[r] = (masks, comp_of(cpP, kP, Y), np.rint(ts.values("P", Y)).astype(np.int64))
    found = []
    t0 = time.time()
    for gi, (members, mult) in enumerate(fs):
        sel = []
        for b in range(nP):
            masks, comps, vals = Pl[pat[b]]
            ok = np.zeros(len(masks), bool)
            want = {tuple(x) for x in members[b]}
            for i, c in enumerate(map(tuple, comps)):
                ok[i] = c in want
            sel.append((masks[ok], vals[ok]))
        if nP == 1:
            codes, pvals = sel[0][0].astype(np.uint64), sel[0][1]
        else:
            (ma, va), (mb, vb) = sel
            if len(ma) * len(mb) > maxside:
                print(f"  group {gi}: {len(ma)} x {len(mb)} pairs, skipped", flush=True); continue
            I, J = np.repeat(np.arange(len(ma)), len(mb)), np.tile(np.arange(len(mb)), len(ma))
            codes = ma[I].astype(np.uint64) | (mb[J].astype(np.uint64) << np.uint64(32))
            pvals = va[I] + vb[J]
        if len(codes) == 0:
            continue
        cpC, kC = tc.blocks[-1][0], tc.blocks[-1][1]
        seen_c = set()
        for xc in members[-1]:
            if canon(xc) in seen_c:
                continue
            seen_c.add(canon(xc))
            inp = Path(W.tmpdir()) / f"tszc_{ts.N}_{ts.mu}_{os.getpid()}.in"
            with open(inp, "wb") as f:
                f.write(_st.pack("6i", 0, ts.OC, ts.CM, pat[0], pat[-1], ts.target))
                f.write(np.asarray(ts.sizesC, np.int32).tobytes())
                f.write(np.asarray(ts.TC, np.int32).tobytes())
                f.write(_st.pack("q", len(codes))); f.write(np.asarray(codes, np.uint64).tobytes())
                f.write(np.ascontiguousarray(pvals.astype(np.int32)).tobytes())
                F = len(cpC.fc[kC])
                f.write(_st.pack("i", F))
                f.write(np.asarray(cpC.fid[kC], np.int32).tobytes()); f.write(np.asarray(cpC.w[kC], np.int32).tobytes())
                f.write(np.asarray(xc, np.int32).tobytes())
            p = sp.run([W.bin("tsz"), str(inp), str(threads)], capture_output=True, text=True)
            inp.unlink()
            sols = [l.split() for l in p.stdout.splitlines() if l.startswith("SOL")]
            stats = [l for l in p.stdout.splitlines() if l.startswith("STATS")]
            if not quiet:
                print(f"  group {gi} xC {xc}: P side {len(codes)}, {stats[0] if stats else p.stderr.strip()[-200:]} "
                      f"[{time.time() - t0:.0f}s]", flush=True)
            for s_ in sols[:3]:
                cp_, mc = int(s_[1]), int(s_[2])
                yc = [1 if (mc >> o) & 1 else -1 for o in range(ts.OC)]
                if ts.kind == "full":
                    seqs = ts.full([1 if (cp_ >> o) & 1 else -1 for o in range(ts.OP)], yc)
                else:
                    ma_, mb_ = cp_ & 0xFFFFFFFF, cp_ >> 32
                    seqs = ts.full([1 if (ma_ >> o) & 1 else -1 for o in range(ts.OP)], yc,
                                   [1 if (mb_ >> o) & 1 else -1 for o in range(ts.OP)])
                rec = {"N": ts.N, "gens": list(ts.gens), "subgroup_order": len(ts.st.H), "twist": ts.kind,
                       "mu": ts.mu, "rowsums": [int(x.sum()) for x in seqs], "method": "tsz_comp",
                       "sequences": ["".join("+" if v > 0 else "-" for v in x) for x in seqs]}
                if ts.bordered:
                    import bsz
                    Kr, Ki, b = bsz.build(seqs)
                    rec["m"] = ts.N - 1
                    rec["border"] = [[int(z.real), int(z.imag)] for z in b]
                    print(f"  verified: complex Hadamard matrix of order {2 * ts.N}, bordered", flush=True)
                else:
                    import ssxx
                    ssxx.complex_hadamard(seqs); ssxx.goethals_seidel(seqs)
                    print(f"  verified: complex Hadamard matrix of order {2 * ts.N} and the Goethals-Seidel matrix "
                          f"of order {4 * ts.N}", flush=True)
                found.append(rec)
                if log:
                    with open(log, "a") as f:
                        f.write(json.dumps(rec) + "\n")
            if found:
                return found
    return found


def open_N():
    out = subprocess.run([sys.executable, str(Path(__file__).resolve().parent / "complex_status.py")],
                         capture_output=True, text=True).stdout
    line = [l for l in out.splitlines() if l.startswith("open (order 2N)")][0]
    return [int(x) // 2 for x in line.split(":")[1].split()]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("args", nargs="*")
    ap.add_argument("--mu", type=int); ap.add_argument("--kind", default="half")
    ap.add_argument("--bordered", action="store_true"); ap.add_argument("--rows")
    ap.add_argument("--maxo", type=int, default=44); ap.add_argument("--mino", type=int, default=1)
    ap.add_argument("--K", type=int, default=400); ap.add_argument("--maxpool", type=float, default=3e10)
    ap.add_argument("--threads", type=int, default=1); ap.add_argument("--log")
    ap.add_argument("--minE", type=float, default=-1.0); ap.add_argument("--maxstream", type=float, default=2e10)
    a = ap.parse_args()
    if a.args and a.args[0] == "scan":
        Ns = [int(x) for x in a.args[1:]] or open_N()
        for N in Ns:
            for ts in tsz.structures(N, a.maxo, a.mino, ("half",) if a.bordered else ("full", "half"), a.bordered):
                try:
                    triage(ts, a.K, maxpool=a.maxpool)
                except Exception as ex:
                    print(f"N={N} {ts.kind} gens={ts.gens} mu={ts.mu}: {type(ex).__name__} {ex}", flush=True)
        return
    if a.args and a.args[0] == "sweep":
        # triage every twisted structure of every N, and search a pattern exactly (restricted to its feasible
        # compressions) when its conditioned estimate is at least --minE, or at least --minE - 2 with a C stream of
        # at most 1e9
        Ns = [int(x) for x in a.args[1:]] or open_N()
        for N in Ns:
            for ts in tsz.structures(N, a.maxo, a.mino, ("half",) if a.bordered else ("full", "half"), a.bordered):
                try:
                    res = triage(ts, a.K, maxpool=a.maxpool)
                except Exception as ex:
                    print(f"N={N} {ts.kind} gens={ts.gens} mu={ts.mu}: {type(ex).__name__} {ex}", flush=True)
                    continue
                for pat, g0, dead, best in res:
                    if dead or not best or best[0] == -math.inf:
                        continue
                    e = best[0]
                    if e < a.minE - 2:
                        continue
                    try:
                        tcs = [tc for tc in quotients(ts)]
                        stream = min(sum(tc.blocks[-1][0].compressions_cache[(tc.blocks[-1][1], x)]
                                         for x in {x for members, _ in tc.feasible(pat) for x in members[-1]})
                                     for tc in tcs)
                        if e < a.minE and stream > 1e9:
                            continue
                        if stream > a.maxstream:
                            print(f"  N={N} rows {pat}: stream {stream:.3g} above --maxstream, not searched", flush=True)
                            continue
                        if search(ts, pat, a.threads, a.log):
                            print(f"N={N}: FOUND ({ts.kind} twist, H = <{ts.gens}>, mu = {ts.mu}, rows {pat})", flush=True)
                    except Exception as ex:
                        print(f"  N={N} rows {pat}: {type(ex).__name__} {ex}", flush=True)
        return
    if a.args and a.args[0] == "search":
        N, gens = int(a.args[1]), [int(g) for g in a.args[2:]]
        ts = tsz.TStruct(N, gens, a.mu, "half" if a.bordered else a.kind, a.bordered)
        pats = [tuple(int(x) for x in a.rows.split(","))] if a.rows else ts.patterns()
        for pat in pats:
            if search(ts, pat, a.threads, a.log):
                break
        return
    N, gens = int(a.args[0]), [int(g) for g in a.args[1:]]
    ts = tsz.TStruct(N, gens, a.mu, "half" if a.bordered else a.kind, a.bordered)
    triage(ts, a.K, maxpool=a.maxpool)


if __name__ == "__main__":
    main()
