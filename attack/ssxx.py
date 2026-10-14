#!/usr/bin/env python3
"""Complex Hadamard matrices of order 2N from four circulant +-1 matrices with two symmetric blocks
(Szöllősi, "A complex Hadamard matrix of order 94", arXiv:2603.09572, Theorem 4).

If A, B, C, D are circulant +-1 matrices of order N with A A^T + B B^T + C C^T + D D^T = 4N I, and A and B are
symmetric, then with R the back-identity

    K = [[ (A+B)/2 + i(A-B)/2,          (C+DR)/2 + i(C-DR)/2 ],
         [ (C^T+DR)/2 - i(C^T-DR)/2,   -(A+B)/2 + i(A-B)/2  ]]

has entries in {1, -1, i, -i} and K K* = 2N I. (Symmetry type (ssxx); the near-Williamson matrices of
Kharaghani-Mohammadian-Tayfeh-Rezaie-Zaitsev are the case (sssx) after reordering, Williamson matrices (ssss).)

Search: all four blocks constant on the orbits of a multiplier group H <= Z_N^*; A and B are then constant on the
orbits of <H, -1>. The condition sum PAF = 0 off 0 holds on the <H, -1>-classes, the four blocks are complete lists
(attack/bnw_enum.c) and the join is exact (attack/join4x.c). Every hit is turned into K and checked exactly
(Gaussian integers), and the real Goethals-Seidel matrix of order 4N is checked as well.

    python3 attack/ssxx.py triage N [N ...] [--maxo 40]
    python3 attack/ssxx.py search N g [g2] [--log FILE] [--K 20]
"""
import argparse, itertools, json, math, os, struct, subprocess, sys, time
from math import gcd, log10, pi
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import workdir as W
from bnw import Struct, subgroups, paf, mask_to_y, join4

KINDS = ["L", "L", "A", "A"]          # A, B symmetric; C, D arbitrary H-invariant


def patterns(st):
    """Row sums (a, b; c, d) with a^2 + b^2 + c^2 + d^2 = 4N, a <= b (A, B symmetric), c <= d, all positive (a block
    and its negative have the same PAF), restricted to the values the orbit structure allows."""
    N = st.N
    ok = {k: {r for r in st.counts(k) if r > 0} for k in ("A", "L")}
    out = []
    for a in sorted(ok["L"]):
        for b in sorted(ok["L"]):
            if b < a:
                continue
            for c in sorted(ok["A"]):
                d2 = 4 * N - a * a - b * b - c * c
                if d2 < c * c:
                    continue
                d = math.isqrt(d2)
                if d * d == d2 and d in ok["A"]:
                    out.append((a, b, c, d))
    return out


def gauss(st, rows, K=3000, seed=1):
    rng = np.random.default_rng(seed)
    lN, mus, vars_, lats = 0.0, [], [], []
    for kind, r in zip(KINDS, rows):
        n = st.counts(kind).get(r, 0)
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
        if var < 1e-9:
            if abs(mu) > 1e-6:
                return -math.inf
            continue
        if g and (-res) % g:
            return -math.inf
        if not dropped:
            dropped = True; continue
        lp += -0.5 * log10(2 * pi * var) - mu ** 2 / (2 * var) / math.log(10) + log10(max(g, 1))
    return lN + lp


def triage(N, maxo=40, mino=1):
    rows = []
    for H, gens in subgroups(N).items():
        if (N - 1) in H:
            continue                                     # -1 in H: C, D symmetric too (Williamson)
        st = Struct(N, gens, plain=True)
        OA, OL = len(st.sizesA), len(st.sizesL)
        if not (mino <= OA <= maxo):
            continue
        for pat in patterns(st):
            e = gauss(st, pat)
            if e == -math.inf:
                continue
            pools = [st.counts(k).get(r, 0) for k, r in zip(KINDS, pat)]
            rows.append((e, len(H), gens, OA, OL, st.C, pat, pools))
    rows.sort(key=lambda r: -r[0])
    print(f"N = {N}: {len(rows)} (ssxx) structure/pattern pairs")
    print(f"{'log10 E':>8} {'|H|':>4} {'gens':>10} {'O_A':>4} {'O_L':>4} {'cls':>4} {'rows':>18} pools")
    for e, h, gens, OA, OL, C, pat, pools in rows[:25]:
        print(f"{e:8.2f} {h:4d} {str(gens):>10} {OA:4d} {OL:4d} {C:4d} {str(pat):>18} {[f'{p:.3g}' for p in pools]}")
    return rows


def enum(st, kind, r):
    d = W.WORK / "ssxx"; d.mkdir(exist_ok=True)
    base = d / f"s{st.N}_{'_'.join(map(str, st.gens))}_{kind}_{r}"
    vals, masks = base.with_suffix(".v"), base.with_suffix(".m")
    if not vals.exists():
        inp = base.with_suffix(".in")
        sz = st.sizes(kind); S = st.S[kind]
        with open(inp, "wb") as f:
            f.write(struct.pack("4i", len(sz), st.C, r, 0))
            f.write(np.asarray(sz, np.int32).tobytes()); f.write(np.asarray(S, np.int32).tobytes())
        subprocess.run([W.bin("bnw_enum"), str(inp), str(vals) + ".tmp", str(masks) + ".tmp"], check=True,
                       capture_output=True)
        os.replace(str(vals) + ".tmp", vals); os.replace(str(masks) + ".tmp", masks)
    return np.fromfile(vals, np.int16).reshape(-1, st.C), np.fromfile(masks, np.uint64)


def circ(x):
    n = len(x); t = np.arange(n)
    return np.asarray(x, np.int64)[(t[None, :] - t[:, None]) % n]


def complex_hadamard(seqs):
    """Szöllősi's K from (A, B, C, D), A and B symmetric; returns K as (real, imag) int64 arrays, checked."""
    A, B, C, D = (circ(x) for x in seqs)
    N = A.shape[0]; t = np.arange(N)
    if not (np.array_equal(A, A.T) and np.array_equal(B, B.T)):
        raise ValueError("A and B must be symmetric")
    if not np.array_equal(A @ A.T + B @ B.T + C @ C.T + D @ D.T, 4 * N * np.eye(N, dtype=np.int64)):
        raise ValueError("sum X X^T != 4N I")
    R = np.zeros((N, N), np.int64); R[t, (N - 1 - t)] = 1
    DR = D @ R
    Xr, Xi = (A + B) // 2, (A - B) // 2
    Yr, Yi = (C + DR) // 2, (C - DR) // 2
    Zr, Zi = (C.T + DR) // 2, -(C.T - DR) // 2
    Kr = np.block([[Xr, Yr], [Zr, -Xr]])
    Ki = np.block([[Xi, Yi], [Zi, Xi]])
    n = 2 * N
    if not np.all((np.abs(Kr) + np.abs(Ki)) == 1):
        raise ValueError("entries not in {1, -1, i, -i}")
    # K K^* = (Kr + i Ki)(Kr^T - i Ki^T) = Kr Kr^T + Ki Ki^T + i (Ki Kr^T - Kr Ki^T)
    G_re = Kr @ Kr.T + Ki @ Ki.T
    G_im = Ki @ Kr.T - Kr @ Ki.T
    if not (np.array_equal(G_re, n * np.eye(n, dtype=np.int64)) and not G_im.any()):
        raise ValueError("K K* != 2N I")
    return Kr, Ki


def goethals_seidel(seqs):
    A, B, C, D = (circ(x) for x in seqs)
    N = A.shape[0]; t = np.arange(N)
    R = np.zeros((N, N), np.int64); R[t, (N - 1 - t)] = 1
    H = np.block([[A, B @ R, C @ R, D @ R], [-B @ R, A, D.T @ R, -C.T @ R],
                  [-C @ R, -D.T @ R, A, B.T @ R], [-D @ R, C.T @ R, -B.T @ R, A]])
    if not np.array_equal(H @ H.T, 4 * N * np.eye(4 * N, dtype=np.int64)):
        raise ValueError("Goethals-Seidel matrix is not Hadamard")
    return H


def quaternary_string(Kr, Ki):
    sym = {(1, 0): "0", (0, 1): "1", (-1, 0): "2", (0, -1): "3"}      # exponent of i
    return ["".join(sym[(int(a), int(b))] for a, b in zip(Kr[i], Ki[i])) for i in range(Kr.shape[0])]


def pair_cost(sizes, order):
    a, b, c, d = (sizes[k] for k in order)
    return a * b + c * d


def best_pairing(sizes):
    """join4x pairs lists (1, 2) against (3, 4): choose the order with the least pair work."""
    return min(([0, 1, 2, 3], [0, 2, 1, 3], [0, 3, 1, 2]), key=lambda o: pair_cost(sizes, o))


def search(N, gens, K=20, log=None, pats=None, maxpairs=2e11, stop=False):
    st = Struct(N, gens, plain=True)
    if (N - 1) in st.H:
        raise ValueError("-1 in H")
    print(f"N = {N}, H = <{','.join(map(str, gens))}> of order {len(st.H)}: O_A = {len(st.sizesA)}, "
          f"O_L = {len(st.sizesL)}, classes = {st.C}", flush=True)
    found = []
    for pat in patterns(st):
        if pats and pat not in pats:
            continue
        t0 = time.time()
        lists, masks = [], []
        for kind, r in zip(KINDS, pat):
            V, M = enum(st, kind, r)
            lists.append(V); masks.append(M)
        sizes = [len(x) for x in lists]
        if min(sizes) == 0:
            continue
        order = best_pairing(sizes)
        if pair_cost(sizes, order) > maxpairs:
            print(f"  rows {pat}: lists {sizes}, skipped (pair work {pair_cost(sizes, order):.2g})", flush=True)
            continue
        sols, total = join4([lists[k] for k in order], st.C, K, None, target=0)
        print(f"  rows {pat}: lists {sizes}, {total} solutions, estimate 10^{gauss(st, pat):.2f} "
              f"({time.time() - t0:.1f}s)", flush=True)
        for sol in sols:
            idx = [0] * 4
            for pos, k in enumerate(order):
                idx[k] = sol[pos]
            seqs = [st.expand(kind, mask_to_y(masks[k][i], len(st.sizes(kind))))
                    for k, (kind, i) in enumerate(zip(KINDS, idx))]
            Kr, Ki = complex_hadamard(seqs)
            goethals_seidel(seqs)
            rec = {"N": N, "gens": list(gens), "subgroup_order": len(st.H), "rowsums": [int(x.sum()) for x in seqs],
                   "sequences": ["".join("+" if v > 0 else "-" for v in x) for x in seqs]}
            found.append(rec)
            if log:
                with open(log, "a") as f:
                    f.write(json.dumps(rec) + "\n")
        if found:
            print(f"  verified: complex Hadamard matrix of order {2 * N} (K K* = {2 * N} I) and the real "
                  f"Goethals-Seidel matrix of order {4 * N}", flush=True)
            if stop:
                break
    return found


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("triage"); p.add_argument("N", type=int, nargs="+"); p.add_argument("--maxo", type=int, default=40)
    p = sub.add_parser("search"); p.add_argument("N", type=int); p.add_argument("g", type=int, nargs="+")
    p.add_argument("--K", type=int, default=20); p.add_argument("--log")
    a = ap.parse_args()
    if a.cmd == "triage":
        for N in a.N:
            triage(N, a.maxo)
    else:
        search(a.N, a.g, a.K, a.log)


if __name__ == "__main__":
    main()
