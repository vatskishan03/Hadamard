#!/usr/bin/env python3
"""Bordered Szöllősi array: complex (quaternary) Hadamard matrices of order 2N from four +-1 sequences on Z_(N-1).

Let m = N - 1 and let A, B, C, D be +-1 sequences on Z_m with A and B symmetric (x(-g) = x(g)) and

    PAF_A(g) + PAF_B(g) + PAF_C(g) + PAF_D(g) = -4      for every g != 0.

Summing over g gives r_A^2 + r_B^2 + r_C^2 + r_D^2 = 4, so one row sum is +-2 and the others are 0. Write
alpha = (1 + i)/2, X' = alpha A + conj(alpha) B and Y' = alpha C + conj(alpha) D R (group-developed matrices, R the
matrix of g -> -g); these are Szöllősi's blocks (arXiv:2603.09572, Theorem 4) on Z_m. Border them,

    X = [[x, s e^T], [t e, X']],   Y = [[y, u e^T], [v e, Y']],      K = [[X, Y], [Y^*, -X^*]],

with x, s, t, y, u, v in {1, i, -1, -i}. Since X' is symmetric, X' Y' = Y' X', X'X'^* + Y'Y'^* = 2N I - 2J and e is an
eigenvector of X', Y' and their adjoints (eigenvalues rho_X = alpha r_A + conj(alpha) r_B, rho_Y = alpha r_C +
conj(alpha) r_D), K K^* = 2N I reduces to five scalar equations:

    s v = u t,   x u + s rho_Y = y s + u rho_X,   t y + v rho_X = v x + t rho_Y,
    x conj(t) + s conj(rho_X) + y conj(v) + u conj(rho_Y) = 0,   conj(x) s + conj(t) rho_X + conj(y) u + conj(v) rho_Y = 0.

`border()` solves them by enumeration for each row-sum pattern; `build()` makes K and checks K K^* = 2N I exactly.

Search (as search/bnw.py, kinds L, L, A, A): all four sequences constant on the orbits of H <= Z_m^*, A and B on
the orbits of <H, -1>; the sum of PAFs is -4 on the <H, -1>-classes of Z_m minus 0; complete lists (bnw_enum.c) and
an exact join (join4x.c).

    python3 search/bsz.py test                         small exhaustive checks (N = 7, 11)
    python3 search/bsz.py triage N [N ...] [--maxo 44]
    python3 search/bsz.py search N g [g2] [--log FILE] [--maxpairs 2e10]
"""
import argparse, itertools, json, math, sys, time
from math import gcd, log10, pi
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from bnw import Struct, subgroups, mask_to_y, join4, enum_block

KINDS = ["L", "L", "A", "A"]
UNITS = [1, 1j, -1, -1j]
ALPHA = (1 + 1j) / 2


_BORDER = {}


def border(rows):
    """All (x, s, t, y, u, v) in {1, i, -1, -i}^6 solving the border equations for row sums rows = (rA, rB, rC, rD)."""
    rows = tuple(int(r) for r in rows)
    if rows not in _BORDER:
        _BORDER[rows] = _border(rows)
    return _BORDER[rows]


def _border(rows):
    rA, rB, rC, rD = rows
    rX = ALPHA * rA + ALPHA.conjugate() * rB
    rY = ALPHA * rC + ALPHA.conjugate() * rD
    out = []
    for x, s, t, y, u, v in itertools.product(UNITS, repeat=6):
        eqs = (s * v - u * t,
               x * u + s * rY - y * s - u * rX,
               t * y + v * rX - v * x - t * rY,
               x * t.conjugate() + s * rX.conjugate() + y * v.conjugate() + u * rY.conjugate(),
               x.conjugate() * s + t.conjugate() * rX + y.conjugate() * u + v.conjugate() * rY)
        if all(abs(e) < 1e-9 for e in eqs):
            out.append((x, s, t, y, u, v))
    return out


def dev(x):
    """group-developed matrix over Z_m: M[g, h] = x(h - g)"""
    m = len(x); t = np.arange(m)
    return np.asarray(x, np.int64)[(t[None, :] - t[:, None]) % m]


def build(seqs, b=None):
    """K (as real and imaginary int64 arrays) from (A, B, C, D) on Z_m; checks every condition exactly."""
    A, B, C, D = (np.asarray(x, np.int64) for x in seqs)
    m = len(A); N = m + 1; t = np.arange(m)
    if not (np.array_equal(A, A[(-t) % m]) and np.array_equal(B, B[(-t) % m])):
        raise ValueError("A and B must be symmetric")
    rows = tuple(int(x.sum()) for x in (A, B, C, D))
    if b is None:
        sols = border(rows)
        if not sols:
            raise ValueError(f"no border for row sums {rows}")
        b = sols[0]
    x, s, t_, y, u, v = b
    MA, MB, MC, MD = dev(A), dev(B), dev(C), dev(D)
    R = np.zeros((m, m), np.int64); R[t, (-t) % m] = 1
    MDR = MD @ R
    Xr, Xi = (MA + MB) // 2, (MA - MB) // 2                    # X' = alpha A + conj(alpha) B
    Yr, Yi = (MC + MDR) // 2, (MC - MDR) // 2                  # Y' = alpha C + conj(alpha) D R

    def bordered(c00, top, left, Mr, Mi):
        Kr = np.empty((N, N), np.int64); Ki = np.empty((N, N), np.int64)
        Kr[0, 0], Ki[0, 0] = int(c00.real), int(c00.imag)
        Kr[0, 1:], Ki[0, 1:] = int(top.real), int(top.imag)
        Kr[1:, 0], Ki[1:, 0] = int(left.real), int(left.imag)
        Kr[1:, 1:], Ki[1:, 1:] = Mr, Mi
        return Kr, Ki

    XR, XI = bordered(x, s, t_, Xr, Xi)
    YR, YI = bordered(y, u, v, Yr, Yi)
    # K = [[X, Y], [Y^*, -X^*]]
    Kr = np.block([[XR, YR], [YR.T, -XR.T]])
    Ki = np.block([[XI, YI], [-YI.T, XI.T]])
    n = 2 * N
    if not np.all(np.abs(Kr) + np.abs(Ki) == 1):
        raise ValueError("entries not in {1, -1, i, -i}")
    I = np.eye(n, dtype=np.int64)
    if not (np.array_equal(Kr @ Kr.T + Ki @ Ki.T, n * I) and not (Ki @ Kr.T - Kr @ Ki.T).any()):
        raise ValueError("K K* != 2N I")
    return Kr, Ki, b


def paf(x):
    X = np.fft.fft(np.asarray(x, float))
    return np.rint(np.fft.ifft(X * np.conj(X)).real).astype(np.int64)


def test():
    """Exhaustive search over Z_6 and Z_10 (N = 7, 11): every solution must give a complex Hadamard matrix."""
    for m in (6, 10):
        t = np.arange(m)
        allv = [np.array(v) for v in itertools.product((1, -1), repeat=m)]
        sym = [v for v in allv if np.array_equal(v, v[(-t) % m])]
        P = {id(v): paf(v) for v in allv}
        tgt = np.full(m, -4); tgt[0] = 4 * m
        from collections import defaultdict
        cd = defaultdict(list)
        for c in allv:
            for d in allv:
                cd[tuple(P[id(c)] + P[id(d)])].append((c, d))
        n = 0; pats = defaultdict(int)
        for a in sym:
            for b in sym:
                need = tuple(tgt - P[id(a)] - P[id(b)])
                for c, d in cd.get(need, []):
                    rows = tuple(int(x.sum()) for x in (a, b, c, d))
                    build([a, b, c, d])
                    n += 1; pats[rows] += 1
        print(f"m = {m} (N = {m + 1}): {n} quadruples, all give K K* = {2 * (m + 1)} I; row-sum patterns {dict(pats)}")
    for rows in [(2, 0, 0, 0), (0, 2, 0, 0), (0, 0, 2, 0), (0, 0, 0, 2), (-2, 0, 0, 0), (0, 0, -2, 0)]:
        print(rows, len(border(rows)), "border solutions, e.g.", border(rows)[:1])


# ---------------------------------------------------------------- search over multiplier structures

def patterns(st):
    """the row-sum-2 block (sign fixed by negating that block)"""
    out = []
    for j, name in enumerate("ABCD"):
        rows = [0, 0, 0, 0]; rows[j] = 2
        if border(tuple(rows)) and st.counts(KINDS[j]).get(2, 0) and st.counts("L").get(0, 0) is not None:
            if all(st.counts(KINDS[k]).get(rows[k], 0) for k in range(4)):
                out.append((name, tuple(rows)))
    return out


def z2_ok(st, rows):
    """compressions to Z_2 (m even): slice sums (a, b), a + b = r, need sum (a - b)^2 = 4N"""
    opts = []
    for kind, r in zip(KINDS, rows):
        sz, orb = st.sizes(kind), st.orb(kind)
        par = np.array([int(np.flatnonzero(orb == k)[0]) % 2 for k in range(len(sz))])
        ab = {(0, 0)}
        for s_, p_ in zip(sz, par):
            ab = {(a + (s_ if p_ == 0 else 0) * e, b + (s_ if p_ == 1 else 0) * e) for a, b in ab for e in (1, -1)}
        opts.append(sorted({(a - b) ** 2 for a, b in ab if a + b == r}))
    return any(sum(c) == 4 * st.N for c in itertools.product(*opts))


def gauss(st, rows, K=3000, seed=1):
    rng = np.random.default_rng(seed)
    lN, mus, vars_, lats = 0.0, [], [], []
    for kind, r in zip(KINDS, rows):
        n = st.counts(kind).get(r, 0)
        if r == 0 and n:
            n //= 2
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
        tg = -4
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


def triage(N, maxo=44, mino=1, show=25):
    rows = []
    for H, gens in subgroups(N - 1).items():
        st = Struct(N, gens, plain=False)
        OA, OL = len(st.sizesA), len(st.sizesL)
        if not (mino <= OA <= maxo):
            continue
        for name, pat in patterns(st):
            if not z2_ok(st, pat):
                continue
            e = gauss(st, pat)
            if e == -math.inf:
                continue
            pools = [st.counts(k).get(r, 0) // (2 if r == 0 else 1) for k, r in zip(KINDS, pat)]
            rows.append((e, len(H), gens, OA, OL, st.C, name, pools))
    rows.sort(key=lambda r: -r[0])
    print(f"N = {N}, m = {N - 1}: {len(rows)} bordered (ssxx) structure/pattern pairs")
    for e, h, gens, OA, OL, C, name, pools in rows[:show]:
        print(f"{e:8.2f} {h:4d} {str(gens):>12} O_A {OA:3d} O_L {OL:3d} cls {C:3d} 2 in {name}  pools "
              f"{[f'{p:.3g}' for p in pools]}")
    return rows


def search(N, gens, K=20, log=None, maxpairs=2e10, pats=None, stop=True):
    st = Struct(N, gens, plain=False)
    print(f"N = {N}, m = {N - 1}, H = <{','.join(map(str, gens))}> of order {len(st.H)}: O_A = {len(st.sizesA)}, "
          f"O_L = {len(st.sizesL)}, classes = {st.C}", flush=True)
    found = []
    for name, pat in patterns(st):
        if pats and name not in pats:
            continue
        if not z2_ok(st, pat):
            print(f"  2 in {name}: excluded by the Z_2 compression", flush=True); continue
        t0 = time.time()
        lists, masks = [], []
        for k, (kind, r) in enumerate(zip(KINDS, pat)):
            V, M = enum_block(st, kind, r, 1 if r == 0 else 0, f"sz{k}")
            lists.append(V); masks.append(M)
        sizes = [len(x) for x in lists]
        if min(sizes) == 0:
            continue
        order = min(([0, 1, 2, 3], [0, 2, 1, 3], [0, 3, 1, 2]),
                    key=lambda o: sizes[o[0]] * sizes[o[1]] + sizes[o[2]] * sizes[o[3]])
        work = sizes[order[0]] * sizes[order[1]] + sizes[order[2]] * sizes[order[3]]
        if work > maxpairs:
            print(f"  2 in {name}: lists {sizes}, skipped (pair work {work:.2g})", flush=True); continue
        sols, total = join4([lists[k] for k in order], st.C, K, None, target=-4)
        print(f"  2 in {name}: lists {sizes}, {total} solutions ({time.time() - t0:.1f}s)", flush=True)
        for sol in sols:
            idx = [0] * 4
            for pos, k in enumerate(order):
                idx[k] = sol[pos]
            seqs = [st.expand(kind, mask_to_y(masks[k][i], len(st.sizes(kind)))) for k, (kind, i) in
                    enumerate(zip(KINDS, idx))]
            Kr, Ki, b = build(seqs)
            rec = {"N": N, "m": N - 1, "gens": list(gens), "subgroup_order": len(st.H),
                   "rowsums": [int(x.sum()) for x in seqs], "border": [[int(z.real), int(z.imag)] for z in b],
                   "sequences": ["".join("+" if v > 0 else "-" for v in x) for x in seqs]}
            found.append(rec)
            if log:
                with open(log, "a") as f:
                    f.write(json.dumps(rec) + "\n")
        if found:
            print(f"  verified: complex Hadamard matrix of order {2 * N} (K K* = {2 * N} I)", flush=True)
            if stop:
                break
    return found


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("test")
    p = sub.add_parser("triage"); p.add_argument("N", type=int, nargs="+"); p.add_argument("--maxo", type=int, default=44)
    p.add_argument("--show", type=int, default=25)
    p = sub.add_parser("search"); p.add_argument("N", type=int); p.add_argument("g", type=int, nargs="+")
    p.add_argument("--K", type=int, default=20); p.add_argument("--log"); p.add_argument("--maxpairs", type=float, default=2e10)
    a = ap.parse_args()
    if a.cmd == "test":
        test()
    elif a.cmd == "triage":
        for N in a.N:
            triage(N, a.maxo, show=a.show)
    else:
        search(a.N, a.g, a.K, a.log, a.maxpairs)


if __name__ == "__main__":
    main()
