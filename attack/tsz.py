#!/usr/bin/env python3
"""Twisted Szöllősi quadruples: complex Hadamard matrices of order 2N from two +-1 sequences A, C on Z_N.

Szöllősi's array (attack/ssxx.py; arXiv:2603.09572, Theorem 4) needs A, B, C, D on Z_N with A, B symmetric and
sum PAF = 0 off 0. Take a multiplier group H <= Z_N^* with -1 not in H, L = <H, -1>, and a unit mu with mu not in L
and mu^2 in L. Put

    B(t) = A(mu t),   D(t) = C(mu t),       A constant on the L-orbits, C constant on the H-orbits.

B is symmetric because A is, PAF_B(s) = PAF_A(mu s), PAF_D(s) = PAF_C(mu s), so with F = PAF_A + PAF_C

    sum PAF (s) = F(s) + F(mu s),

which is invariant under s -> mu s: the condition holds on the orbits of M = <L, mu> ("classes"), half as many as the
L-classes of the untwisted search, while the unknowns are two sequences instead of four. The row sums satisfy
r_A^2 + r_C^2 = 2N, so N must be a sum of two squares (N = 1 mod 4 for prime N).

The search is a two-list join: every A (L-orbits, row sum r_A) is hashed by its class values, every C (H-orbits, row
sum r_C) is enumerated by a depth-first search with incremental class values and looked up (attack/tsz.c).

    python3 attack/tsz.py triage N [N ...] [--maxo 40]          rank the twisted structures (Gaussian estimate)
    python3 attack/tsz.py search N g [g2] --mu MU [--rows rA,rC] [--log FILE]
    python3 attack/tsz.py test                                   small cases, every hit checked exactly
"""
import argparse, json, math, os, struct, subprocess, sys, time
from math import gcd, log10, pi
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import workdir as W
from bnw import Struct, subgroups, closure
import ssxx
import bsz


class TStruct:
    """kind "full": B = A o mu, D = C o mu (A on L-orbits, C on H-orbits), classes = orbits of <L, mu>.
    kind "half": A, B on the orbits of L' = <L, mu> (both symmetric and mu-invariant), D = C o mu, classes = L'-classes.
    Sides: "P" = the hashed side (A, or the pairs (A, B)), "C" = the streamed side."""

    def __init__(self, N, gens, mu, kind="full", bordered=False):
        st = Struct(N, gens, plain=not bordered)
        m = st.m                           # N (cyclic) or N - 1 (bordered: one border point, target -4)
        self.N, self.m, self.gens, self.mu, self.kind, self.bordered = N, m, tuple(gens), mu % m, kind, bordered
        self.target = st.target
        self.st = st
        L = st.L
        if self.mu in L or (self.mu * self.mu) % m not in L or gcd(self.mu, m) != 1:
            raise ValueError("mu must be a unit outside L = <H, -1> with mu^2 in L")
        if bordered and kind == "full":
            raise ValueError("the full twist forces r_A = r_B and r_C = r_D, impossible for row sums (2, 0, 0, 0)")
        lab = st.orbL                      # L-orbit label of every element (orbit 0 = {0}); class c = orbit c + 1
        SL, SH = st.S["L"], st.S["A"]
        if kind == "full":
            partner = [int(lab[(self.mu * s) % m]) - 1 for s in st.reps]
            seen, groups = set(), []
            for c in range(st.C):
                if c in seen:
                    continue
                p = partner[c]
                seen.add(c); seen.add(p)
                groups.append((c, p))
            self.groups = groups
            self.CM = len(groups)
            self.TP = np.stack([SL[c] + SL[p] if p != c else SL[c] for c, p in groups])
            self.TC = np.stack([SH[c] + SH[p] if p != c else SH[c] for c, p in groups])
            self.sizesP, self.sizesC = st.sizesL, st.sizesA
            self.orbP = st.orbL
        else:
            st2 = Struct(N, tuple(gens) + (self.mu,), plain=not bordered)
            self.st2 = st2
            self.CM = st2.C
            rows = []
            for s in st2.reps:
                c, p = int(lab[s]) - 1, int(lab[(self.mu * s) % m]) - 1
                rows.append(SH[c] + SH[p] if p != c else 2 * SH[c])
            self.TP = st2.S["L"]
            self.TC = np.stack(rows)
            self.sizesP, self.sizesC = st2.sizesL, st.sizesA
            self.orbP = st2.orbL
        self.OP, self.OC = len(self.sizesP), len(self.sizesC)
        self.OL, self.OA = self.OP, self.OC

    @staticmethod
    def _counts(sizes):
        cnt = {0: 1}
        for sz in sizes:
            nxt = {}
            for t, c in cnt.items():
                nxt[t + sz] = nxt.get(t + sz, 0) + c
                nxt[t - sz] = nxt.get(t - sz, 0) + c
            cnt = nxt
        return cnt

    def counts(self, side):
        return self._counts(self.sizesP if side == "P" else self.sizesC)

    def sample(self, side, r, K, rng):
        st = self.st2 if (side == "P" and self.kind == "half") else self.st
        kind = "L" if side == "P" else "A"
        return st.sample(kind, r, K, rng)

    def values(self, side, Y):
        T = (self.TP if side == "P" else self.TC).astype(np.float64)
        return np.stack([((Y @ T[m]) * Y).sum(1) for m in range(self.CM)], 1)

    def blocks(self, pat):
        """(side, row sum) of every independent block: full (A, C), half (A, B, C)"""
        return [("P", pat[0]), ("C", pat[1])] if self.kind == "full" else [("P", pat[0]), ("P", pat[1]), ("C", pat[2])]

    def patterns(self):
        N = self.N
        if self.bordered:                  # half twist: r_A^2 + r_B^2 + 2 r_C^2 = 4, so (2, 0, 0) up to A <-> B
            cP, cC = self.counts("P"), self.counts("C")
            return [(0, 2, 0)] if cP.get(0) and cP.get(2) and cC.get(0) else []
        okP = sorted(r for r in self.counts("P") if r > 0)
        okC = {r for r in self.counts("C") if r > 0}
        out = []
        if self.kind == "full":
            for a in okP:
                c2 = 2 * N - a * a
                if c2 <= 0:
                    break
                c = math.isqrt(c2)
                if c * c == c2 and c in okC:
                    out.append((int(a), int(c)))
        else:
            for c in sorted(okC):
                rest = 4 * N - 2 * c * c
                if rest <= 0:
                    break
                for a in okP:
                    b2 = rest - a * a
                    if b2 < a * a:
                        break
                    b = math.isqrt(b2)
                    if b * b == b2 and b in okP:
                        out.append((int(a), int(b), int(c)))
        return out

    def gauss(self, pat, K=3000, seed=1):
        """log10 of the expected number of solutions for row sums pat, Gaussian model on the classes (one class
        dropped for the row-sum identity, lattice of the class values from the samples)."""
        rng = np.random.default_rng(seed)
        lN, mus, vars_, lats = 0.0, [], [], []
        for side, r in self.blocks(pat):
            n = self.counts(side).get(r, 0)
            if n == 0:
                return -math.inf
            lN += log10(n)
            V = self.values(side, self.sample(side, r, K, rng))
            mus.append(V.mean(0)); vars_.append(V.var(0))
            Vi = np.rint(V).astype(np.int64)
            lat = []
            for m in range(self.CM):
                u = np.unique(Vi[:, m])
                g = int(np.gcd.reduce(np.diff(u))) if len(u) > 1 else 0
                lat.append((g, int(u[0] % g) if g else int(u[0])))
            lats.append(lat)
        if self.kind == "half" and pat[0] == pat[1]:
            lN -= log10(2)                 # (A, B) and (B, A) are the same solution
        lp, dropped = 0.0, False
        for m in range(self.CM):
            mu = sum(x[m] for x in mus); var = sum(x[m] for x in vars_)
            g, res = 0, 0
            for lat in lats:
                g = gcd(g, lat[m][0]); res += lat[m][1]
            if var < 1e-9:
                if abs(mu - self.target) > 1e-6:
                    return -math.inf
                continue
            if g and (self.target - res) % g:
                return -math.inf
            if not dropped:
                dropped = True; continue
            lp += -0.5 * log10(2 * pi * var) - (mu - self.target) ** 2 / (2 * var) / math.log(10) + log10(max(g, 1))
        return lN + lp

    def full(self, yp, yc, yb=None):
        """the four sequences on Z_m (m = N, or N - 1 when bordered)"""
        m, mu = self.m, self.mu
        t = np.arange(m)
        A = np.asarray(yp, np.int64)[self.orbP]
        C = np.asarray(yc, np.int64)[self.st.orbA]
        B = A[(mu * t) % m] if self.kind == "full" else np.asarray(yb, np.int64)[self.orbP]
        D = C[(mu * t) % m]
        return A, B, C, D


def twists(m, H, L):
    """representatives mu of the cosets of L in <L, mu> = M with [M : L] = 2, one per subgroup M (units mod m)"""
    units = [u for u in range(2, m) if gcd(u, m) == 1 and u not in L and (u * u) % m in L]
    out, seenM = [], set()
    for u in units:
        M = frozenset(L | {(u * x) % m for x in L})
        if M in seenM:
            continue
        seenM.add(M); out.append(u)
    return out


def sum_two_squares(n):
    for a in range(1, math.isqrt(n) + 1):
        b2 = n - a * a
        b = math.isqrt(b2)
        if b * b == b2:
            return True
    return False


def structures(N, maxo=40, mino=1, kinds=("full", "half"), bordered=False):
    """every twisted structure (TStruct) with mino <= O_C <= maxo"""
    out = []
    m = N - 1 if bordered else N
    for H, gens in subgroups(m).items():
        if (m - 1) in H:
            continue
        st = Struct(N, gens, plain=not bordered)
        OC = len(st.sizesA)
        if not (mino <= OC <= maxo):
            continue
        for mu in twists(m, H, st.L):
            for kind in kinds:
                if kind == "full" and (bordered or not sum_two_squares(2 * N)):
                    continue
                try:
                    out.append(TStruct(N, gens, mu, kind, bordered))
                except ValueError:
                    pass
    return out


def side_size(ts, pat):
    """(size of the hashed side, size of the streamed C pool)"""
    cP, cC = ts.counts("P"), ts.counts("C")
    if ts.kind == "full":
        return cP.get(pat[0], 0), cC.get(pat[1], 0)
    a, b = cP.get(pat[0], 0), cP.get(pat[1], 0)
    return (a * (a + 1) // 2 if pat[0] == pat[1] else a * b), cC.get(pat[2], 0)


def triage(N, maxo=40, mino=1, top=12, quiet=False, kinds=("full", "half"), bordered=False):
    rows = []
    for ts in structures(N, maxo, mino, kinds, bordered):
        for pat in ts.patterns():
            e = ts.gauss(pat)
            if e == -math.inf:
                continue
            rows.append((e, ts, pat, side_size(ts, pat)))
    rows.sort(key=lambda r: -r[0])
    if not quiet:
        print(f"N = {N}: {len(rows)} twisted structure/pattern pairs", flush=True)
        for e, ts, pat, sz in rows[:top]:
            print(f"  logE {e:6.2f} {ts.kind:4s}{'/b' if ts.bordered else '  '} |H| {len(ts.st.H):3d} gens {str(ts.gens):>10} mu {ts.mu:4d} O_P {ts.OP:3d} "
                  f"O_C {ts.OC:3d} classes {ts.CM:3d} rows {pat} sides {[f'{p:.3g}' for p in sz]}", flush=True)
    return rows


# ---------------------------------------------------------------------------------------------------- exact search
def all_signs(sizes, r):
    """every sign vector on orbits of the given sizes with row sum r, as (masks, Y)"""
    O = len(sizes)
    if O > 26:
        raise ValueError("too many orbits for explicit enumeration")
    masks = np.arange(1 << O, dtype=np.int64)
    Y = ((masks[:, None] >> np.arange(O)) & 1) * 2 - 1
    keep = Y @ np.asarray(sizes, np.int64) == r
    return masks[keep], Y[keep].astype(np.float64)


def write_input(ts, pat, path, maxpairs=3e7):
    """tsz.c input (see tsz.c): the full twist enumerates A in C; the half twist sends the pairs (A, B) explicitly"""
    with open(path, "wb") as f:
        if ts.kind == "full":
            f.write(struct.pack("6i", ts.OP, ts.OC, ts.CM, pat[0], pat[1], ts.target))
            f.write(np.asarray(ts.sizesP, np.int32).tobytes()); f.write(np.asarray(ts.sizesC, np.int32).tobytes())
            f.write(np.asarray(ts.TP, np.int32).tobytes()); f.write(np.asarray(ts.TC, np.int32).tobytes())
            return
        ma, Ya = all_signs(ts.sizesP, pat[0])
        mb, Yb = all_signs(ts.sizesP, pat[1])
        va = np.rint(ts.values("P", Ya)).astype(np.int64); vb = np.rint(ts.values("P", Yb)).astype(np.int64)
        if len(ma) * len(mb) > maxpairs:
            raise ValueError(f"{len(ma)} x {len(mb)} pairs")
        if pat[0] == pat[1]:
            I, J = np.triu_indices(len(ma))
        else:
            I, J = np.repeat(np.arange(len(ma)), len(mb)), np.tile(np.arange(len(mb)), len(ma))
        codes = (ma[I].astype(np.uint64)) | (mb[J].astype(np.uint64) << np.uint64(32))
        vals = (va[I] + vb[J]).astype(np.int32)
        f.write(struct.pack("6i", 0, ts.OC, ts.CM, pat[0], pat[2], ts.target))
        f.write(np.asarray(ts.sizesC, np.int32).tobytes())
        f.write(np.asarray(ts.TC, np.int32).tobytes())
        f.write(struct.pack("q", len(codes))); f.write(codes.tobytes()); f.write(np.ascontiguousarray(vals).tobytes())


def run_one(ts, pat, threads=4, log=None, frac=1.0, quiet=False):
    t0 = time.time()
    N = ts.N
    inp = W.WORK / "tmp" / f"tsz_{N}_{'_'.join(map(str, ts.gens))}_{ts.mu}_{ts.kind}_{'_'.join(map(str, pat))}_{os.getpid()}.in"
    inp.parent.mkdir(exist_ok=True)
    write_input(ts, pat, inp)
    p = subprocess.run([W.bin("tsz"), str(inp), str(threads), str(frac)], capture_output=True, text=True)
    inp.unlink()
    sols = [l.split() for l in p.stdout.splitlines() if l.startswith("SOL")]
    stats = [l for l in p.stdout.splitlines() if l.startswith("STATS")]
    if not quiet or sols:
        print(f"  N = {N} {ts.kind} H = <{','.join(map(str, ts.gens))}> mu = {ts.mu} rows {pat}: {len(sols)} solutions "
              f"{stats[0] if stats else p.stderr.strip()[-200:]} [{time.time() - t0:.0f}s]", flush=True)
    found = []
    for s in sols[:3]:
        cp, mc = int(s[1]), int(s[2])
        yc = [1 if (mc >> o) & 1 else -1 for o in range(ts.OC)]
        if ts.kind == "full":
            yp = [1 if (cp >> o) & 1 else -1 for o in range(ts.OP)]
            seqs = ts.full(yp, yc)
        else:
            ma, mb = cp & 0xFFFFFFFF, cp >> 32
            yp = [1 if (ma >> o) & 1 else -1 for o in range(ts.OP)]
            yb = [1 if (mb >> o) & 1 else -1 for o in range(ts.OP)]
            seqs = ts.full(yp, yc, yb)
        rec = {"N": N, "gens": list(ts.gens), "subgroup_order": len(ts.st.H), "twist": ts.kind, "mu": ts.mu,
               "rowsums": [int(x.sum()) for x in seqs], "method": "tsz",
               "sequences": ["".join("+" if v > 0 else "-" for v in x) for x in seqs]}
        if ts.bordered:
            Kr, Ki, b = bsz.build(seqs)               # raises unless K K* = 2N I
            rec["m"] = N - 1
            rec["border"] = [[int(z.real), int(z.imag)] for z in b]
            print(f"  verified: complex Hadamard matrix of order {2 * N} (K K* = {2 * N} I), bordered", flush=True)
        else:
            Kr, Ki = ssxx.complex_hadamard(seqs)      # raises unless K K* = 2N I
            ssxx.goethals_seidel(seqs)
            print(f"  verified: complex Hadamard matrix of order {2 * N} (K K* = {2 * N} I) and the real "
                  f"Goethals-Seidel matrix of order {4 * N}", flush=True)
        found.append(rec)
        if log:
            with open(log, "a") as f:
                f.write(json.dumps(rec) + "\n")
    return found, (stats[0] if stats else None)


def search(N, gens, mu, kind="full", pats=None, threads=4, log=None, frac=1.0, stop=True, bordered=False):
    ts = TStruct(N, gens, mu, kind, bordered)
    found = []
    for pat in pats or ts.patterns():
        f, _ = run_one(ts, pat, threads, log, frac)
        found += f
        if f and stop:
            break
    return found


def sweep(Ns, maxo=40, maxpool=3e10, maxside=3e7, budget=3600, threads=4, log=None, kinds=("full", "half"),
          bordered=False):
    """every twisted structure/pattern of every N within the size limits, cheapest first, until a solution"""
    for N in Ns:
        t0 = time.time()
        jobs = []
        for ts in structures(N, maxo, 1, kinds, bordered):
            for pat in ts.patterns():
                sp, sc = side_size(ts, pat)
                if sp == 0 or sc == 0 or sc > maxpool or sp > maxside:
                    continue
                jobs.append((sc + sp, ts, pat))
        jobs.sort(key=lambda j: j[0])
        print(f"N = {N}: {len(jobs)} twisted structure/pattern pairs", flush=True)
        hit, done = False, 0
        for cost, ts, pat in jobs:
            if time.time() - t0 > budget:
                print(f"N = {N}: budget exhausted after {done} of {len(jobs)}", flush=True)
                break
            try:
                f, st = run_one(ts, pat, threads, log, quiet=True)
            except Exception as ex:
                print(f"  N = {N} {ts.kind} <{ts.gens}> mu {ts.mu} {pat}: {type(ex).__name__} {ex}", flush=True)
                continue
            done += 1
            if f:
                print(f"N = {N}: FOUND ({ts.kind} twist, H = <{','.join(map(str, ts.gens))}>, mu = {ts.mu}, rows {pat})",
                      flush=True)
                hit = True
                break
        if not hit:
            print(f"N = {N}: nothing ({done} of {len(jobs)} done, {time.time() - t0:.0f}s)", flush=True)


def test():
    ok = 0
    for N in (5, 13, 15, 17, 21, 25, 29, 33, 35, 39, 41, 51, 55, 57, 61, 63, 65, 69, 75, 77):
        rows = triage(N, maxo=30, quiet=True)
        res = []
        for kind in ("full", "half"):
            rk = [r for r in rows if r[1].kind == kind]
            if not rk:
                res.append(f"{kind}: none")
                continue
            got = None
            for e, ts, pat, sz in rk[:4]:
                if sz[0] > 3e6 or sz[1] > 3e8:
                    continue
                f, stt = run_one(ts, pat, threads=2, quiet=True)
                if f:
                    got = (ts.gens, ts.mu, pat, round(e, 2), stt.split()[6] if stt else "?")
                    break
            res.append(f"{kind}: {got}")
            ok += got is not None
        print(f"N = {N}:", "; ".join(res), flush=True)
    print("test: structures with solutions:", ok)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["triage", "search", "sweep", "test"])
    ap.add_argument("args", nargs="*", type=int)
    ap.add_argument("--maxo", type=int, default=40)
    ap.add_argument("--mino", type=int, default=1)
    ap.add_argument("--mu", type=int)
    ap.add_argument("--kind", default="full")
    ap.add_argument("--kinds", default="full,half")
    ap.add_argument("--rows")
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--frac", type=float, default=1.0)
    ap.add_argument("--maxpool", type=float, default=3e10)
    ap.add_argument("--budget", type=float, default=3600)
    ap.add_argument("--log")
    ap.add_argument("--bordered", action="store_true", help="sequences on Z_(N-1), sum PAF = -4 (BORDERED.md)")
    a = ap.parse_args()
    kinds = tuple(a.kinds.split(","))
    if a.bordered:
        kinds = ("half",)
    if a.cmd == "triage":
        for N in a.args:
            triage(N, a.maxo, a.mino, kinds=kinds, bordered=a.bordered)
    elif a.cmd == "search":
        N, gens = a.args[0], a.args[1:]
        pats = [tuple(map(int, a.rows.split(",")))] if a.rows else None
        search(N, gens, a.mu, "half" if a.bordered else a.kind, pats, a.threads, a.log, a.frac, bordered=a.bordered)
    elif a.cmd == "sweep":
        sweep(a.args, a.maxo, a.maxpool, 3e7, a.budget, a.threads, a.log, kinds, a.bordered)
    else:
        test()


if __name__ == "__main__":
    main()
