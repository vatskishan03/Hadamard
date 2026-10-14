#!/usr/bin/env python3
"""Williamson-type matrices from twisted multiplier structures: the order-3 twist of Kharaghani, Mohammadian,
Tayfeh-Rezaie and Zaitsev (arXiv:2605.08661, section 4.4, type (ii)), on Z_N (near-Williamson matrices) or on
Z_(N-1) (the bordered construction of williamson-type/README.md).

Let u have multiplicative order 3 in Z_m^* and let K' <= Z_m^* be a subgroup with u not in <K', -1>. The first
block A is constant on the orbits of H_A = <u, K'> (and symmetric as well with --symA); B is symmetric and constant
on the orbits of K_B = <-1, K'>, and C = B o u, D = B o u^2 (C(g) = B(ug)). C and D are symmetric, and

    PAF_A(g) + PAF_B(g) + PAF_B(ug) + PAF_B(u^2 g)

is constant on the orbits of L = <u, -1, K'>. So the condition (= 0 on Z_N, = -4 on Z_(N-1), for g != 0) is a
two-list problem on the L-classes: the A-list against the B-list. Row sums: r_C = r_D = r_B, so the bordered case
needs (r_A, r_B) = (+-2, 0) and the plain case r_A^2 + 3 r_B^2 = 4N. The point of the bordered form: when N is prime
with 3 not dividing N - 1 (N = 71, 107, 167, ...), Z_N^* has no element of order 3, but Z_(N-1)^* may have one
(Z_70^*: u = 11).

The A-list is enumerated completely (attack/bnw_enum.c) and held in a hash table; the B-list is streamed in Gray
order and probed (attack/twist_join.c), so the search of a structure is exhaustive.

    python3 attack/wt_twist.py triage N [--plain] [--maxa 4e8]      structures, Gaussian estimates, list sizes
    python3 attack/wt_twist.py search N u [k ...] [--plain] [--symA] [--log FILE] [--passes P]
"""
import argparse, itertools, json, math, os, struct, subprocess, sys, time
from math import gcd, log10, pi
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import workdir as W
from bnw import closure, orbits, paf, mask_to_y, near_williamson_type, williamson_type, williamson_array


def order3(m):
    """Elements of order 3 in Z_m^*, one of each pair {u, u^2}."""
    out = []
    for u in range(2, m):
        if gcd(u, m) == 1 and pow(u, 3, m) == 1 and u % m != 1 and pow(u, 2, m) not in out:
            out.append(u)
    return out


def small_subgroups(m, maxorder=8):
    units = [x for x in range(1, m) if gcd(x, m) == 1]
    subs = {frozenset({1}): ()}
    for x in units:
        H = closure([x], m)
        if len(H) <= maxorder:
            subs.setdefault(H, (x,))
    return subs


class Twist:
    def __init__(self, N, u, kgens=(), plain=False, symA=False):
        self.N, self.u, self.kgens, self.plain, self.symA = N, u, tuple(kgens), plain, symA
        m = self.m = N if plain else N - 1
        self.target = 0 if plain else -4
        if gcd(u, m) != 1 or pow(u, 3, m) != 1 or u % m == 1:
            raise ValueError(f"{u} does not have order 3 mod {m}")
        self.HA = closure([u, *kgens] + ([m - 1] if symA else []), m)
        self.KB = closure([m - 1, *kgens], m)
        if u in self.KB:
            raise ValueError("u lies in <K', -1>: the twist is trivial")
        self.Lg = closure([u, m - 1, *kgens], m)
        self.orbA, self.sizesA = orbits(m, self.HA)
        self.orbB, self.sizesB = orbits(m, self.KB)
        orbL, sizesL = orbits(m, self.Lg)
        self.reps = [int(np.flatnonzero(orbL == k)[0]) for k in range(1, len(sizesL))]
        self.csize = sizesL[1:]
        self.C = len(self.reps)
        t = np.arange(m)
        SA = np.zeros((self.C, len(self.sizesA), len(self.sizesA)), np.int64)
        SB = np.zeros((self.C, len(self.sizesB), len(self.sizesB)), np.int64)
        for c, s in enumerate(self.reps):
            np.add.at(SA[c], (self.orbA, self.orbA[(t + s) % m]), 1)
            for e in range(3):
                np.add.at(SB[c], (self.orbB, self.orbB[(t + pow(u, e, m) * s) % m]), 1)
        self.S = {"A": SA, "B": SB}

    def sizes(self, kind):
        return self.sizesA if kind == "A" else self.sizesB

    def orb(self, kind):
        return self.orbA if kind == "A" else self.orbB

    def counts(self, kind):
        cnt = {0: 1}
        for s in self.sizes(kind):
            nxt = {}
            for t, c in cnt.items():
                nxt[t + s] = nxt.get(t + s, 0) + c
                nxt[t - s] = nxt.get(t - s, 0) + c
            cnt = nxt
        return cnt

    def patterns(self):
        if not self.plain:
            return [(2, 0)]
        out = []
        for b in range(1, math.isqrt(4 * self.N // 3) + 1, 2):
            a2 = 4 * self.N - 3 * b * b
            if a2 > 0 and math.isqrt(a2) ** 2 == a2 and math.isqrt(a2) % 2 == 1:
                out.append((math.isqrt(a2), b))
        return out

    def z2_ok(self, pat):
        """Bordered only: the quotient Z_m -> Z_2. u is odd, so C and D have B's even/odd sums, and
        (a_A - b_A)^2 + 3 (a_B - b_B)^2 = 4N is necessary."""
        if self.plain:
            return True
        opts = []
        for kind, r in zip("AB", pat):
            sz, orb = self.sizes(kind), self.orb(kind)
            par = [int(np.flatnonzero(orb == k)[0]) % 2 for k in range(len(sz))]
            ab = {(0, 0)}
            for s, p in zip(sz, par):
                ab = {(a + (s if p == 0 else 0) * e, b + (s if p == 1 else 0) * e) for a, b in ab for e in (1, -1)}
            opts.append(sorted({(a - b) ** 2 for a, b in ab if a + b == r}))
        return any(x + 3 * y == 4 * self.N for x in opts[0] for y in opts[1])

    def sample(self, kind, r, K, rng):
        sz = self.sizes(kind); O = len(sz)
        tab = [dict() for _ in range(O + 1)]
        tab[O][0] = 1
        for i in range(O - 1, -1, -1):
            d = {}
            for t, c in tab[i + 1].items():
                d[t + sz[i]] = d.get(t + sz[i], 0) + c
                d[t - sz[i]] = d.get(t - sz[i], 0) + c
            tab[i] = d
        Y = np.empty((K, O))
        for k in range(K):
            rem = r
            for i in range(O):
                cp = tab[i + 1].get(rem - sz[i], 0); cm = tab[i + 1].get(rem + sz[i], 0)
                if rng.random() * (cp + cm) < cp:
                    Y[k, i] = 1; rem -= sz[i]
                else:
                    Y[k, i] = -1; rem += sz[i]
        return Y

    def values(self, kind, Y):
        S = self.S[kind].astype(np.float64)
        return np.stack([((Y @ S[c]) * Y).sum(1) for c in range(self.C)], 1)

    def pools(self, pat):
        nA = self.counts("A").get(pat[0], 0)
        nB = self.counts("B").get(pat[1], 0)
        return nA, (nB // 2 if pat[1] == 0 else nB)

    def gauss(self, pat, K=3000, seed=1):
        """log10 of the Gaussian estimate: pool sizes times the class densities (with lattice spacings), one class
        dropped for the row-sum identity."""
        rng = np.random.default_rng(seed)
        nA, nB = self.pools(pat)
        if not nA or not nB:
            return -math.inf
        lN = log10(nA) + log10(nB)
        mus, vars_, lats = [], [], []
        for kind, r in zip("AB", pat):
            V = self.values(kind, self.sample(kind, r, K, rng))
            mus.append(V.mean(0)); vars_.append(V.var(0))
            Vi = np.rint(V).astype(np.int64)
            lat = []
            for c in range(self.C):
                u = np.unique(Vi[:, c])
                g = int(np.gcd.reduce(np.diff(u))) if len(u) > 1 else 0
                lat.append((g, int(u[0] % g) if g else int(u[0])))
            lats.append(lat)
        lp, dropped = 0.0, False
        for c in range(self.C):
            mu = mus[0][c] + mus[1][c]; var = vars_[0][c] + vars_[1][c]
            g = gcd(lats[0][c][0], lats[1][c][0]); res = lats[0][c][1] + lats[1][c][1]
            tg = self.target
            if var < 1e-9:
                if abs(tg - mu) > 1e-6:
                    return -math.inf
                continue
            if g and (tg - res) % g:
                return -math.inf
            if not dropped:
                dropped = True; continue
            lp += -0.5 * log10(2 * pi * var) - (tg - mu) ** 2 / (2 * var) / math.log(10) + log10(max(g, 1))
        return lN + lp

    def expand(self, kind, y):
        return np.asarray(y, np.int64)[self.orb(kind)]

    def quadruple(self, yA, yB):
        A = self.expand("A", yA); B = self.expand("B", yB)
        t = np.arange(self.m)
        return [A, B, B[(self.u * t) % self.m], B[(self.u * self.u * t) % self.m]]

    def tag(self):
        return f"{'p' if self.plain else 'b'}{self.N}_u{self.u}_k{'-'.join(map(str, self.kgens)) or '1'}" + \
            ("_sA" if self.symA else "")


def triage(N, plain=False, maxa=4e8, maxk=8):
    m = N if plain else N - 1
    rows = []
    for u in order3(m):
        for Kp, kg in small_subgroups(m, maxk).items():
            for symA in (False, True):
                try:
                    tw = Twist(N, u, kg, plain, symA)
                except ValueError:
                    continue
                for pat in tw.patterns():
                    if not tw.z2_ok(pat):
                        continue
                    e = tw.gauss(pat)
                    if e == -math.inf:
                        continue
                    nA, nB = tw.pools(pat)
                    rows.append((e, u, kg, symA, pat, len(tw.sizesA), len(tw.sizesB), tw.C, nA, nB))
    rows.sort(key=lambda r: -r[0])
    print(f"N = {N}, m = {m} ({'plain' if plain else 'bordered'}): {len(rows)} twisted structures")
    print(f"{'log10 E':>8} {'u':>4} {'K':>10} {'symA':>5} {'rows':>10} {'O_A':>4} {'O_B':>4} {'cls':>4} {'pool A':>10} {'pool B':>10}")
    for e, u, kg, symA, pat, OA, OB, C, nA, nB in rows:
        flag = "" if nA <= maxa else "  (A-list too large)"
        print(f"{e:8.2f} {u:4d} {str(kg):>10} {str(symA):>5} {str(pat):>10} {OA:4d} {OB:4d} {C:4d} {nA:10.3g} {nB:10.3g}{flag}")
    return rows


def enum_A(tw, r):
    d = W.WORK / "twist"; d.mkdir(exist_ok=True)
    base = d / f"{tw.tag()}_A_{r}"
    vals, masks = base.with_suffix(".v"), base.with_suffix(".m")
    if not vals.exists():
        inp = base.with_suffix(".in")
        S = tw.S["A"]; sz = tw.sizesA
        with open(inp, "wb") as f:
            f.write(struct.pack("4i", len(sz), tw.C, r, 0))
            f.write(np.asarray(sz, np.int32).tobytes()); f.write(np.asarray(S, np.int32).tobytes())
        subprocess.run([W.bin("bnw_enum"), str(inp), str(vals) + ".tmp", str(masks) + ".tmp"], check=True,
                       capture_output=True)
        os.replace(str(vals) + ".tmp", vals); os.replace(str(masks) + ".tmp", masks)
    return vals, masks


def search(N, u, kgens=(), plain=False, symA=False, log=None, passes=1, K=100):
    tw = Twist(N, u, kgens, plain, symA)
    print(f"N = {N}, m = {tw.m}, u = {u}, K' = <{','.join(map(str, kgens))}>, symA = {symA}: O_A = {len(tw.sizesA)}, "
          f"O_B = {len(tw.sizesB)}, classes = {tw.C}", flush=True)
    found = []
    for pat in tw.patterns():
        if not tw.z2_ok(pat):
            print(f"  rows {pat}: Z_2 test fails"); continue
        t0 = time.time()
        vals, masks = enum_A(tw, pat[0])
        nA = os.path.getsize(vals) // (2 * tw.C)
        print(f"  rows {pat}: A-list {nA} ({time.time() - t0:.0f}s), Gaussian estimate 10^{tw.gauss(pat):.2f}",
              flush=True)
        if nA == 0:
            continue
        d = W.WORK / "twist"
        inp = d / f"{tw.tag()}_B_{pat[1]}.in"
        S = tw.S["B"]; sz = tw.sizesB
        with open(inp, "wb") as f:
            f.write(struct.pack("4i", len(sz), tw.C, pat[1], 1 if pat[1] == 0 else 0))   # -B, -C, -D: same PAFs
            f.write(np.asarray(sz, np.int32).tobytes()); f.write(np.asarray(S, np.int32).tobytes())
            f.write(np.full(tw.C, tw.target, np.int16).tobytes())
        out = subprocess.run([W.bin("twist_join"), str(inp), str(vals), str(K), str(passes)], capture_output=True,
                             text=True, check=True)
        sys.stderr.write(out.stderr)
        sols = [tuple(map(int, ln.split()[1:3])) for ln in out.stdout.splitlines() if ln.startswith("SOL")]
        total = [int(ln.split()[1]) for ln in out.stdout.splitlines() if ln.startswith("TOTAL")]
        print(f"  rows {pat}: {total[0] if total else '?'} solutions ({time.time() - t0:.0f}s)", flush=True)
        M = np.fromfile(masks, np.uint64)
        for a, bmask in sols:
            yA = mask_to_y(M[a], len(tw.sizesA)); yB = mask_to_y(bmask, len(tw.sizesB))
            seqs = tw.quadruple(yA, yB)
            if plain:
                Ws = near_williamson_type(seqs)
            else:
                Ws = williamson_type(seqs)
            williamson_array(Ws)
            rec = {"N": N, "m": tw.m, "plain": plain, "twist": u, "kgens": list(kgens), "symA": symA,
                   "gens": [u, *kgens], "rowsums": [int(x.sum()) for x in seqs],
                   "sequences": ["".join("+" if v > 0 else "-" for v in x) for x in seqs]}
            found.append(rec)
            if log:
                with open(log, "a") as f:
                    f.write(json.dumps(rec) + "\n")
        if found:
            print(f"  {len(found)} verified (Williamson array H({4 * N}) checked); first rows "
                  f"{found[0]['rowsums']}", flush=True)
    return found


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("triage"); p.add_argument("N", type=int, nargs="+"); p.add_argument("--plain", action="store_true")
    p.add_argument("--maxa", type=float, default=4e8); p.add_argument("--maxk", type=int, default=8)
    p = sub.add_parser("search"); p.add_argument("N", type=int); p.add_argument("u", type=int)
    p.add_argument("k", type=int, nargs="*"); p.add_argument("--plain", action="store_true")
    p.add_argument("--symA", action="store_true"); p.add_argument("--log"); p.add_argument("--passes", type=int, default=1)
    p.add_argument("--K", type=int, default=100)
    a = ap.parse_args()
    if a.cmd == "triage":
        for N in a.N:
            triage(N, a.plain, a.maxa, a.maxk)
    else:
        search(a.N, a.u, a.k, a.plain, a.symA, a.log, a.passes, a.K)


if __name__ == "__main__":
    main()
