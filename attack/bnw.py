#!/usr/bin/env python3
"""Bordered near-Williamson (BNW) search: Williamson-type matrices of odd order N from sequences on Z_{N-1}.

Construction. Let m = N - 1 (even) and let x_1 (= A), x_2, x_3, x_4 (= B, C, D) be +-1 sequences on Z_m with B, C, D
symmetric (x(-g) = x(g)) and

    PAF_A(g) + PAF_B(g) + PAF_C(g) + PAF_D(g) = -4      for every g != 0.

Then sum r_k^2 = 4 for the row sums, so one row sum is +-2 (block j) and the others are 0. With M_k the group-developed
matrices (M_1 replaced by the back-developed M_1 R, which is symmetric) and s = r_j / 2,

    W_k = [[1, b_k e^T], [c_k e, M_k]],   c_k = 1 (k != j), c_j = -1,   b_k = s (k != j), b_j = -s,

are Williamson-type matrices of order N: W_k W_l^T = W_l W_k^T and sum W_k W_k^T = 4N I (the derivation is in
williamson-type/README.md; every result is checked numerically). With T-matrices of order t (Cooper-Wallis), this gives a
Hadamard matrix of order 4tN; t = 1, 3, 5, ... are elementary. The point of working on Z_{N-1} instead of Z_N: when N is
prime, Z_N^* may have no useful multiplier subgroup (179, 263), while Z_{N-1}^* does.

Symmetry: all four sequences constant on the orbits of a multiplier group H <= Z_m^*, B, C, D on the orbits of
L = <H, -1>. Then the condition holds on the L-orbits of Z_m \\ {0} ("classes").

    python3 attack/bnw.py triage N [--maxo 40]          rank the multiplier structures
    python3 attack/bnw.py search N g1 [g2] [--maxo 30]  exact search (complete pools, join4x), every row pattern
    python3 attack/bnw.py check FILE.json                rebuild and verify a solution
"""
import argparse, itertools, json, math, os, struct, subprocess, sys, time
from math import gcd, log10, pi
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import workdir as W


def closure(gens, m):
    S, todo = {1}, [1]
    while todo:
        x = todo.pop()
        for g in gens:
            y = x * g % m
            if y not in S:
                S.add(y); todo.append(y)
    return frozenset(S)


def subgroups(m):
    units = [u for u in range(1, m) if gcd(u, m) == 1]
    cyc = {}
    for u in units:
        cyc.setdefault(closure([u], m), (u,))
    subs = dict(cyc)
    items = list(cyc.items())
    for i in range(len(items)):
        for j in range(i + 1, len(items)):
            subs.setdefault(closure([items[i][1][0], items[j][1][0]], m), (items[i][1][0], items[j][1][0]))
    return subs


def orbits(m, H):
    lab = -np.ones(m, np.int64); sizes = []
    for t in range(m):
        if lab[t] < 0:
            O = sorted({t * h % m for h in H})
            lab[O] = len(sizes); sizes.append(len(O))
    return lab, np.array(sizes, np.int64)


class Struct:
    def __init__(self, N, gens, plain=False):
        """plain=False: bordered, sequences on Z_(N-1), target -4; plain=True: near-Williamson, sequences on Z_N,
        target 0 (Kharaghani, Mohammadian, Tayfeh-Rezaie, Zaitsev), Williamson-type (A R, B, C, D)."""
        self.N, self.plain, self.gens = N, plain, tuple(gens)
        self.m = N if plain else N - 1
        self.target = 0 if plain else -4
        m = self.m
        self.H = closure(gens, m)
        self.L = frozenset(self.H | {(-h) % m for h in self.H})
        self.orbA, self.sizesA = orbits(m, self.H)
        self.orbL, self.sizesL = orbits(m, self.L)
        cl = self.orbL.copy()
        self.reps = [int(np.flatnonzero(cl == k)[0]) for k in range(1, len(self.sizesL))]   # orbit 0 is {0}
        self.C = len(self.reps)
        self.csize = self.sizesL[1:]
        t = np.arange(m)
        self.S = {}
        for kind, orb, O in (("A", self.orbA, len(self.sizesA)), ("L", self.orbL, len(self.sizesL))):
            S = np.zeros((self.C, O, O), np.int64)
            for c, s in enumerate(self.reps):
                np.add.at(S[c], (orb, orb[(t + s) % m]), 1)
            self.S[kind] = S

    def sizes(self, kind):
        return self.sizesA if kind == "A" else self.sizesL

    def orb(self, kind):
        return self.orbA if kind == "A" else self.orbL

    def counts(self, kind):
        """number of sign vectors with each row sum"""
        cnt = {0: 1}
        for s in self.sizes(kind):
            nxt = {}
            for t, c in cnt.items():
                nxt[t + s] = nxt.get(t + s, 0) + c
                nxt[t - s] = nxt.get(t - s, 0) + c
            cnt = nxt
        return cnt

    def sample(self, kind, r, K, rng):
        """K uniform sign vectors with row sum r (backward DP over the orbits)."""
        sz = self.sizes(kind); O = len(sz)
        tab = [dict() for _ in range(O + 1)]
        tab[O][0] = 1
        for i in range(O - 1, -1, -1):
            d = {}
            for t, c in tab[i + 1].items():
                d[t + sz[i]] = d.get(t + sz[i], 0) + c
                d[t - sz[i]] = d.get(t - sz[i], 0) + c
            tab[i] = d
        if r not in tab[0]:
            return np.zeros((0, O))
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

    def expand(self, kind, y):
        return np.asarray(y, np.int64)[self.orb(kind)]


def paf(x):
    X = np.fft.fft(np.asarray(x, float))
    return np.rint(np.fft.ifft(X * np.conj(X)).real).astype(np.int64)


def patterns(st):
    """Row-sum assignments up to symmetry. Bordered: the block with row sum 2 is A or B. Plain: positive odd
    (a; b <= c <= d) with a^2 + b^2 + c^2 + d^2 = 4N (negating a sequence keeps its PAF)."""
    if not st.plain:
        return [("A", (2, 0, 0, 0)), ("B", (0, 2, 0, 0))]
    out, N = [], st.N
    odd = range(1, int(math.isqrt(4 * N)) + 1, 2)
    for a in odd:
        for b in odd:
            for c in odd:
                if not b <= c:
                    continue
                d2 = 4 * N - a * a - b * b - c * c
                if d2 >= c * c and math.isqrt(d2) ** 2 == d2:
                    out.append((f"{a}|{b},{c},{math.isqrt(d2)}", (a, b, c, math.isqrt(d2))))
    return out


def kinds():
    return ["A", "L", "L", "L"]


def z2_ok(st, rows):
    if st.plain:
        return True
    """Exact test on the quotient Z_m -> Z_2 (slice sums): each block's even/odd sums (a, b) with a + b = r must give
    sum (a - b)^2 = 4N. Enumerated per block by DP over orbits (an orbit lies in one coset of 2Z_m)."""
    m = st.m
    opts = []
    for kind, r in zip(kinds(), rows):
        sz, orb = st.sizes(kind), st.orb(kind)
        par = np.array([int(np.flatnonzero(orb == k)[0]) % 2 for k in range(len(sz))])
        ab = {(0, 0)}
        for s, p in zip(sz, par):
            ab = {(a + (s if p == 0 else 0) * e, b + (s if p == 1 else 0) * e) for a, b in ab for e in (1, -1)}
        opts.append(sorted({(a - b) ** 2 for a, b in ab if a + b == r}))
    tgt = 4 * st.N
    for combo in itertools.product(*opts):
        if sum(combo) == tgt:
            return True
    return False


def gauss(st, rows, K=3000, seed=1):
    """log10 of the Gaussian estimate of the number of solutions (pools x product of class densities, with
    lattice spacings, dropping the class fixed by the row-sum identity)."""
    rng = np.random.default_rng(seed)
    lN, mus, vars_, lats = 0.0, [], [], []
    for kind, r in zip(kinds(), rows):
        n = st.counts(kind).get(r, 0)
        if r == 0 and n:
            n //= 2                                   # global sign
        if n == 0:
            return -math.inf
        lN += log10(n)
        V = st.values(kind, st.sample(kind, r, K, rng))
        mus.append(V.mean(0)); vars_.append(V.var(0))
        Vi = np.rint(V).astype(np.int64)
        lat = []
        for c in range(st.C):
            u = np.unique(Vi[:, c])
            g = int(np.gcd.reduce(np.diff(u))) if len(u) > 1 else 0
            lat.append((g, int(u[0] % g) if g else int(u[0])))
        lats.append(lat)
    lp, dropped = 0.0, False
    for c in range(st.C):
        mu = sum(x[c] for x in mus); var = sum(x[c] for x in vars_)
        g, res = 0, 0
        for lat in lats:
            g = gcd(g, lat[c][0]); res += lat[c][1]
        tg = st.target
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


def triage(N, mino, maxo, maxo_l, plain=False):
    m = N if plain else N - 1
    rows = []
    for H, gens in subgroups(m).items():
        st = Struct(N, gens, plain)
        OA, OL = len(st.sizesA), len(st.sizesL)
        if not (mino <= OA <= maxo) or OL > maxo_l:
            continue
        for name, pat in patterns(st):
            if not z2_ok(st, pat):
                continue
            e = gauss(st, pat)
            if e == -math.inf:
                continue
            nA = st.counts("A").get(pat[0], 0); nL = st.counts("L").get(pat[1], 0); nL0 = st.counts("L").get(pat[2], 0)
            rows.append((e, len(H), gens, OA, OL, st.C, name, nA, nL, nL0))
    rows.sort(key=lambda r: -r[0])
    print(f"N = {N}, m = {m}: {len(rows)} admissible structure/pattern pairs")
    print(f"{'log10 E':>8} {'|H|':>4} {'gens':>12} {'O_A':>4} {'O_L':>4} {'cls':>4} {'2 in':>4} {'pool A':>10} {'pool B':>10} {'pool L0':>10}")
    for e, h, gens, OA, OL, C, name, nA, nL, nL0 in rows:
        print(f"{e:8.2f} {h:4d} {str(gens):>12} {OA:4d} {OL:4d} {C:4d} {name:>4} {nA:10.3g} {nL:10.3g} {nL0:10.3g}")
    return rows


def enum_block(st, kind, r, fix0, tag):
    """All sign vectors of a block with row sum r: (int16 class values, uint64 masks)."""
    sz = st.sizes(kind); S = st.S[kind]
    d = W.WORK / "bnw"; d.mkdir(exist_ok=True)
    base = d / f"{'p' if st.plain else 'b'}{st.N}_{'_'.join(map(str, st.gens))}_{tag}_{kind}_{r}_{fix0}"
    vals, masks = base.with_suffix(".v"), base.with_suffix(".m")
    if not vals.exists():
        inp = base.with_suffix(".in")
        with open(inp, "wb") as f:
            f.write(struct.pack("4i", len(sz), st.C, r, fix0))
            f.write(np.asarray(sz, np.int32).tobytes()); f.write(np.asarray(S, np.int32).tobytes())
        subprocess.run([W.bin("bnw_enum"), str(inp), str(vals) + ".tmp", str(masks) + ".tmp"], check=True,
                       capture_output=True)
        os.replace(str(vals) + ".tmp", vals); os.replace(str(masks) + ".tmp", masks)
    V = np.fromfile(vals, np.int16).reshape(-1, st.C)
    M = np.fromfile(masks, np.uint64)
    return V, M


def mask_to_y(mask, O):
    return np.array([-1 if (int(mask) >> i) & 1 else 1 for i in range(O)], np.int64)


def join4(lists, C, K=20, budget=None, target=-4):
    d = W.WORK / "bnw"; d.mkdir(exist_ok=True)
    path = d / f"join_{os.getpid()}.bin"
    with open(path, "wb") as f:
        f.write(struct.pack("5i", *[len(x) for x in lists], C))
        f.write(np.full(C, target, np.int16).tobytes())
        for x in lists:
            f.write(np.ascontiguousarray(x, np.int16).tobytes())
    cmd = [W.bin("join4x"), str(path), str(K)] + ([str(budget)] if budget else [])
    out = subprocess.run(cmd, capture_output=True, text=True, check=True).stdout
    path.unlink()
    sols = [tuple(map(int, ln.split()[1:5])) for ln in out.splitlines() if ln.startswith("SOL")]
    total = [int(ln.split()[1]) for ln in out.splitlines() if ln.startswith("TOTAL")]
    return sols, (total[0] if total else None)


def kt4(lists, C, nsig=4000, K=20, nmerge=5, target=-4):
    """Memory-bounded four-list k-tree (attack/kt4.c) over the split values of the nmerge classes with the largest
    spread, nearest the centre first. Exhaustive only when nsig covers every split value."""
    d = W.WORK / "bnw"; d.mkdir(exist_ok=True)
    path = d / f"kt4_{os.getpid()}.bin"
    with open(path, "wb") as f:
        f.write(struct.pack("5i", *[len(x) for x in lists], C))
        f.write(np.full(C, target, np.int16).tobytes())
        for x in lists:
            f.write(np.ascontiguousarray(x, np.int16).tobytes())
    spread = sum(np.asarray(x, np.float64).var(0) for x in lists)
    S = ",".join(map(str, np.argsort(-spread)[:nmerge]))
    out = subprocess.run([W.bin("kt4"), str(path), S, str(nsig), str(K)], capture_output=True, text=True,
                         check=True).stdout
    path.unlink()
    sols = [tuple(map(int, ln.split()[1:5])) for ln in out.splitlines() if ln.startswith("SOL")]
    total = [int(ln.split()[1]) for ln in out.splitlines() if ln.startswith("TOTAL")]
    return sols, (total[0] if total else None)


def check_bnw(seqs):
    """Verify the sequence conditions; return the index of the row-sum-2 block."""
    m = len(seqs[0])
    t = np.arange(m)
    for x in seqs[1:]:
        if not np.array_equal(x, x[(-t) % m]):
            raise ValueError("B, C, D must be symmetric")
    P = sum(paf(x) for x in seqs)
    if P[0] != 4 * m or np.any(P[1:] != -4):
        raise ValueError("sum of PAFs is not -4 off zero")
    r = [int(x.sum()) for x in seqs]
    j = [k for k in range(4) if abs(r[k]) == 2]
    if sorted(map(abs, r)) != [0, 0, 0, 2] or len(j) != 1:
        raise ValueError(f"row sums {r}")
    return j[0], r


def near_williamson_type(seqs):
    """(A R, B, C, D) from near-Williamson sequences on Z_N (B, C, D symmetric, PAF sum 0 off 0), checked."""
    N = len(seqs[0]); t = np.arange(N)
    for x in seqs[1:]:
        if not np.array_equal(x, x[(-t) % N]):
            raise ValueError("B, C, D must be symmetric")
    P = sum(paf(x) for x in seqs)
    if P[0] != 4 * N or np.any(P[1:] != 0):
        raise ValueError("sum of PAFs is not 0 off zero")
    Ms = [np.asarray(x, np.int64)[(t[None, :] - t[:, None]) % N] for x in seqs]
    R = np.zeros((N, N), np.int64); R[t, (-t) % N] = 1
    Ms[0] = Ms[0] @ R
    I = np.eye(N, dtype=np.int64)
    if not all(np.array_equal(Ms[k] @ Ms[l].T, Ms[l] @ Ms[k].T) for k in range(4) for l in range(k + 1, 4)) or \
            not np.array_equal(sum(x @ x.T for x in Ms), 4 * N * I):
        raise ValueError("not Williamson-type")
    return Ms


def williamson_type(seqs):
    """Bordered Williamson-type matrices W_1..W_4 of order m + 1 (checked)."""
    j, r = check_bnw(seqs)
    m = len(seqs[0]); N = m + 1
    t = np.arange(m)
    Ms = [np.asarray(x, np.int64)[(t[None, :] - t[:, None]) % m] for x in seqs]
    R = np.zeros((m, m), np.int64); R[t, (-t) % m] = 1
    Ms[0] = Ms[0] @ R
    s = r[j] // 2
    Ws = []
    for k in range(4):
        b = -s if k == j else s
        c = -1 if k == j else 1
        Wk = np.empty((N, N), np.int64)
        Wk[0, 0] = 1; Wk[0, 1:] = b; Wk[1:, 0] = c; Wk[1:, 1:] = Ms[k]
        Ws.append(Wk)
    I = np.eye(N, dtype=np.int64)
    if not all(np.array_equal(Ws[k] @ Ws[l].T, Ws[l] @ Ws[k].T) for k in range(4) for l in range(k + 1, 4)):
        raise ValueError("not pairwise amicable")
    if not np.array_equal(sum(x @ x.T for x in Ws), 4 * N * I):
        raise ValueError("sum W W^T != 4N I")
    return Ws


def williamson_array(Ws):
    A, B, C, D = Ws
    H = np.block([[A, B, C, D], [-B, A, -D, C], [-C, D, A, -B], [-D, -C, B, A]])
    n = H.shape[0]
    if not np.array_equal(H @ H.T, n * np.eye(n, dtype=np.int64)):
        raise ValueError("Williamson array is not Hadamard")
    return H


def search(N, gens, pats=None, K=20, budget=None, log=None, sub=None, seed=1, nsig=None, plain=False):
    st = Struct(N, gens, plain)
    OA, OL = len(st.sizesA), len(st.sizesL)
    print(f"N = {N}, H = <{','.join(map(str, gens))}> of order {len(st.H)}: O_A = {OA}, O_L = {OL}, classes = {st.C}",
          flush=True)
    found = []
    for name, pat in patterns(st):
        if pats and name not in pats:
            continue
        if not z2_ok(st, pat):
            print(f"  pattern {name}: Z_2 slice test fails"); continue
        t0 = time.time()
        lists, masks = [], []
        for k, (kind, r) in enumerate(zip(kinds(), pat)):
            V, M = enum_block(st, kind, r, 1 if r == 0 else 0, f"{k}" if not plain else "p")
            if sub and len(V) > sub:                          # uniform sample of the pool
                keep = np.sort(np.random.default_rng(seed + k).choice(len(V), int(sub), replace=False))
                V, M = V[keep], M[keep]
            lists.append(V); masks.append(M)
        sizes = [len(x) for x in lists]
        print(f"  pattern {name}: lists {sizes} ({time.time() - t0:.1f}s)", flush=True)
        if min(sizes) == 0:
            continue
        sols, total = kt4(lists, st.C, nsig, K, target=st.target) if nsig else \
            join4(lists, st.C, K, budget, target=st.target)
        print(f"  pattern {name}: {total} solutions ({time.time() - t0:.1f}s)", flush=True)
        for sol in sols:
            seqs = []
            for k, (kind, idx) in enumerate(zip(kinds(), sol)):
                y = mask_to_y(masks[k][idx], len(st.sizes(kind)))
                seqs.append(st.expand(kind, y))
            Ws = near_williamson_type(seqs) if plain else williamson_type(seqs)
            rec = {"N": N, "m": st.m, "plain": plain, "gens": list(gens), "subgroup_order": len(st.H), "pattern": name,
                   "rowsums": [int(x.sum()) for x in seqs],
                   "sequences": ["".join("+" if v > 0 else "-" for v in x) for x in seqs]}
            found.append(rec)
            if log:
                with open(log, "a") as f:
                    f.write(json.dumps(rec) + "\n")
        if found:
            print(f"  verified {len(found)} Williamson-type quadruples of order {N}", flush=True)
            break
    return found


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    a = sub.add_parser("triage"); a.add_argument("N", type=int, nargs="+")
    a.add_argument("--mino", type=int, default=1); a.add_argument("--maxo", type=int, default=60)
    a.add_argument("--maxol", type=int, default=60); a.add_argument("--plain", action="store_true")
    b = sub.add_parser("search"); b.add_argument("N", type=int); b.add_argument("gens", type=int, nargs="+")
    b.add_argument("--pat", nargs="*"); b.add_argument("--K", type=int, default=20)
    b.add_argument("--budget", type=float); b.add_argument("--log")
    b.add_argument("--sub", type=float, help="sample each list down to this size"); b.add_argument("--seed", type=int, default=1)
    b.add_argument("--nsig", type=int, help="use the k-tree (kt4) over this many split values instead of join4x")
    b.add_argument("--plain", action="store_true", help="near-Williamson on Z_N instead of bordered on Z_(N-1)")
    c = sub.add_parser("check"); c.add_argument("file")
    args = ap.parse_args()
    if args.cmd == "triage":
        for N in args.N:
            triage(N, args.mino, args.maxo, args.maxol, args.plain)
    elif args.cmd == "search":
        search(args.N, args.gens, args.pat, args.K, int(args.budget) if args.budget else None, args.log, args.sub, args.seed, args.nsig, args.plain)
    else:
        rec = json.loads(Path(args.file).read_text())
        seqs = [np.array([1 if ch == "+" else -1 for ch in s]) for s in rec["sequences"]]
        Ws = near_williamson_type(seqs) if rec.get("plain") else williamson_type(seqs)
        H = williamson_array(Ws)
        print(f"Williamson-type matrices of order {len(Ws[0])} verified; Williamson array H({H.shape[0]}) Hadamard")


if __name__ == "__main__":
    main()
