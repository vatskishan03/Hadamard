#!/usr/bin/env python3
"""Bordered Szöllősi array (search/bsz.py) over an arbitrary abelian group G of order N - 1, with symmetry from
any group of automorphisms, in particular the spin structures of G = Z_2^2 x Z_m (N = 4m + 1): the order-3
automorphism alpha(x1, x2) = (x2, x1 + x2) of Z_2^2 coupled with a unit of Z_m.

Four +-1 sequences A, B, C, D on G, A and B symmetric (x(-g) = x(g)), sum of PAFs -4 off 0, row sums one +-2 and
three 0, give the complex Hadamard matrix K = [[X, Y], [Y^*, -X^*]] of order 2N with
X = [[x, s e^T], [t e, alpha A + conj(alpha) B]], Y = [[y, u e^T], [v e, alpha C + conj(alpha) D R]] (group-developed
matrices over G, R the matrix of g -> -g, border from bsz.border). Every hit is rebuilt and checked exactly.

Symmetry: S = <automorphisms>; C, D constant on S-orbits, A, B on <S, -1>-orbits; the condition holds on the
<S, -1>-classes. Complete lists (bnw_enum.c), exact join (join4x.c).

    python3 search/bsz_group.py spin N [--maxo 50] [--maxpairs 2e10] [--log FILE] [--rank]
    python3 search/bsz_group.py test
"""
import argparse, json, math, os, struct, subprocess, sys, time
from math import gcd, log10, pi
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import workdir as W
import spence_lib as L
from bnw import join4, mask_to_y
import bsz

KINDS = bsz.KINDS          # L, L, A, A


class BStruct:
    def __init__(self, N, mods, mats, tag):
        self.N, self.mods, self.mats, self.tag = N, list(mods), [np.asarray(M, np.int64) for M in mats], tag
        self.G = L.AbGroup(self.mods)
        assert self.G.v == N - 1
        perms = [self.G.perm_from_matrix(M) for M in self.mats]
        self.stA = L.Structure(self.G, perms)
        self.stL = L.Structure(self.G, perms + [self.G.neg()])
        assert np.array_equal(self.stA.cls, self.stL.cls)
        self.C = self.stA.C

    def st(self, kind):
        return self.stL if kind == "L" else self.stA

    def count(self, kind, r):
        n = int(round(self.st(kind).pool(r)))
        return n

    def patterns(self):
        out = []
        for j, name in enumerate("ABCD"):
            rows = [0, 0, 0, 0]; rows[j] = 2
            if bsz.border(tuple(rows)) and all(self.count(KINDS[k], rows[k]) for k in range(4)):
                out.append((name, tuple(rows)))
        return out

    def gauss(self, rows, K=2000, seed=1):
        rng = np.random.default_rng(seed)
        lN, mus, vars_, lats = 0.0, [], [], []
        for kind, r in zip(KINDS, rows):
            n = self.count(kind, r)
            if r == 0 and n:
                n //= 2
            if n == 0:
                return -math.inf
            lN += log10(n)
            st = self.st(kind)
            V = st.class_values(st.sample(r, K, rng))
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

    def enum(self, kind, r):
        st = self.st(kind)
        d = W.WORK / "bszg"; d.mkdir(exist_ok=True)
        base = d / f"b{self.N}_{self.tag}_{kind}_{r}"
        vals, masks = base.with_suffix(".v"), base.with_suffix(".m")
        if not vals.exists():
            inp = base.with_suffix(".in")
            S = np.rint(st.S).astype(np.int32)
            with open(inp, "wb") as f:
                f.write(struct.pack("4i", st.O, self.C, r, 1 if r == 0 else 0))
                f.write(np.asarray(st.sizes, np.int32).tobytes()); f.write(S.tobytes())
            subprocess.run([W.bin("bnw_enum"), str(inp), str(vals) + ".tmp", str(masks) + ".tmp"], check=True,
                           capture_output=True)
            os.replace(str(vals) + ".tmp", vals); os.replace(str(masks) + ".tmp", masks)
        return np.fromfile(vals, np.int16).reshape(-1, self.C), np.fromfile(masks, np.uint64)


def group_dev(x, G):
    diff = G.index(G.coords[None, :, :] - G.coords[:, None, :])          # diff[g, h] = h - g
    return np.asarray(x, np.int64)[diff]


def build(seqs, G, b=None):
    """complex Hadamard matrix of order 2N from (A, B, C, D) on G (checked exactly)"""
    A, B, C, D = (np.asarray(x, np.int64) for x in seqs)
    m = G.v; N = m + 1; neg = G.neg()
    if not (np.array_equal(A, A[neg]) and np.array_equal(B, B[neg])):
        raise ValueError("A and B must be symmetric")
    rows = tuple(int(x.sum()) for x in (A, B, C, D))
    if b is None:
        sols = bsz.border(rows)
        if not sols:
            raise ValueError(f"no border for row sums {rows}")
        b = sols[0]
    x, s, t_, y, u, v = b
    MA, MB, MC, MD = (group_dev(z, G) for z in (A, B, C, D))
    R = np.zeros((m, m), np.int64); R[np.arange(m), neg] = 1
    MDR = MD @ R
    Xr, Xi = (MA + MB) // 2, (MA - MB) // 2
    Yr, Yi = (MC + MDR) // 2, (MC - MDR) // 2

    def bordered(c00, top, left, Mr, Mi):
        Kr = np.empty((N, N), np.int64); Ki = np.empty((N, N), np.int64)
        Kr[0, 0], Ki[0, 0] = int(c00.real), int(c00.imag)
        Kr[0, 1:], Ki[0, 1:] = int(top.real), int(top.imag)
        Kr[1:, 0], Ki[1:, 0] = int(left.real), int(left.imag)
        Kr[1:, 1:], Ki[1:, 1:] = Mr, Mi
        return Kr, Ki

    XR, XI = bordered(x, s, t_, Xr, Xi)
    YR, YI = bordered(y, u, v, Yr, Yi)
    Kr = np.block([[XR, YR], [YR.T, -XR.T]])
    Ki = np.block([[XI, YI], [-YI.T, XI.T]])
    n = 2 * N
    if not np.all(np.abs(Kr) + np.abs(Ki) == 1):
        raise ValueError("entries not in {1, -1, i, -i}")
    I = np.eye(n, dtype=np.int64)
    if not (np.array_equal(Kr @ Kr.T + Ki @ Ki.T, n * I) and not (Ki @ Kr.T - Kr @ Ki.T).any()):
        raise ValueError("K K* != 2N I")
    return Kr, Ki, b


def search(bs, K=20, log=None, maxpairs=2e10, stop=True, quiet=False):
    found = []
    for name, pat in bs.patterns():
        t0 = time.time()
        pools = [bs.count(k, r) // (2 if r == 0 else 1) for k, r in zip(KINDS, pat)]
        po = min(([0, 1, 2, 3], [0, 2, 1, 3], [0, 3, 1, 2]),
                 key=lambda o: pools[o[0]] * pools[o[1]] + pools[o[2]] * pools[o[3]])
        pwork = pools[po[0]] * pools[po[1]] + pools[po[2]] * pools[po[3]]
        if pwork > maxpairs or max(pools) > 3e8:
            print(f"  2 in {name}: pools {[f'{x:.3g}' for x in pools]}, skipped (pair work {pwork:.2g})", flush=True)
            continue
        lists, masks = [], []
        for kind, r in zip(KINDS, pat):
            V, M = bs.enum(kind, r)
            lists.append(V); masks.append(M)
        sizes = [len(x) for x in lists]
        if min(sizes) == 0:
            continue
        order = min(([0, 1, 2, 3], [0, 2, 1, 3], [0, 3, 1, 2]),
                    key=lambda o: sizes[o[0]] * sizes[o[1]] + sizes[o[2]] * sizes[o[3]])
        work = sizes[order[0]] * sizes[order[1]] + sizes[order[2]] * sizes[order[3]]
        if work > maxpairs:
            print(f"  2 in {name}: lists {sizes}, skipped (pair work {work:.2g})", flush=True); continue
        sols, total = join4([lists[k] for k in order], bs.C, K, None, target=-4)
        if not quiet or total:
            print(f"  2 in {name}: lists {sizes}, {total} solutions ({time.time() - t0:.1f}s)", flush=True)
        for sol in sols:
            idx = [0] * 4
            for pos, k in enumerate(order):
                idx[k] = sol[pos]
            seqs = [mask_to_y(masks[k][i], bs.st(kind).O)[bs.st(kind).orb] for k, (kind, i) in
                    enumerate(zip(KINDS, idx))]
            Kr, Ki, b = build(seqs, bs.G)
            rec = {"N": N_of(bs), "group": bs.mods, "automorphisms": [M.tolist() for M in bs.mats], "tag": bs.tag,
                   "rowsums": [int(x.sum()) for x in seqs], "border": [[int(z.real), int(z.imag)] for z in b],
                   "sequences": ["".join("+" if v > 0 else "-" for v in x) for x in seqs]}
            found.append(rec)
            if log:
                with open(log, "a") as f:
                    f.write(json.dumps(rec) + "\n")
        if found:
            print(f"  verified: complex Hadamard matrix of order {2 * bs.N} (K K* = {2 * bs.N} I)", flush=True)
            if stop:
                break
    return found


def N_of(bs):
    return bs.N


# ---------------------------------------------------------------- spin structures on Z_2^2 x Z_m

ALPHA2 = np.array([[0, 1], [1, 1]], np.int64)        # (x1, x2) -> (x2, x1 + x2), order 3


def units(m):
    return [u for u in range(1, m) if gcd(u, m) == 1] if m > 1 else [1]


def spin_structures(N, maxo=50, mino=4):
    """S = <(alpha^a, u1), (1, u2)> on G = Z_2^2 x Z_m, m = (N - 1)/4; one representative per orbit partition."""
    assert (N - 1) % 4 == 0
    m = (N - 1) // 4
    mods = [2, 2, m] if m > 1 else [2, 2]
    G = L.AbGroup(mods)

    def mat(a, u):
        M = np.eye(len(mods), dtype=np.int64)
        M[:2, :2] = np.linalg.matrix_power(ALPHA2, a) % 2
        if m > 1:
            M[2, 2] = u
        return M

    seen, out = set(), []
    U = units(m)
    for a in (1, 0):
        for u1 in U:
            for u2 in [1] + U:
                if a == 0 and u1 != 1:
                    continue
                mats = [mat(a, u1)] + ([mat(0, u2)] if u2 != 1 else [])
                perms = [G.perm_from_matrix(M) for M in mats]
                orb, O = L.closure_partition(G.v, perms)
                key = tuple(orb.tolist())
                if key in seen or not (mino <= O <= maxo):
                    seen.add(key); continue
                seen.add(key)
                out.append((mods, mats, f"s{a}-{u1}-{u2}", O))
    return out


def run_spin(N, maxo=50, maxpairs=2e10, log=None, rank=False, K=20, minE=-6.0):
    rows = []
    for mods, mats, tag, O in spin_structures(N, maxo):
        bs = BStruct(N, mods, mats, tag)
        for name, pat in bs.patterns():
            e = bs.gauss(pat)
            if e == -math.inf:
                continue
            pools = [bs.count(k, r) // (2 if r == 0 else 1) for k, r in zip(KINDS, pat)]
            rows.append((e, tag, bs.stA.O, bs.stL.O, bs.C, name, pools, mods, mats))
    rows.sort(key=lambda r: -r[0])
    print(f"N = {N}: {len(rows)} spin structure/pattern pairs (G = Z_2^2 x Z_{(N - 1) // 4})", flush=True)
    for e, tag, OA, OL, C, name, pools, mods, mats in rows[:30]:
        print(f"{e:8.2f} {tag:>14} O_A {OA:3d} O_L {OL:3d} cls {C:3d} 2 in {name} pools {[f'{p:.3g}' for p in pools]}",
              flush=True)
    if rank:
        return rows
    done = set()
    for e, tag, OA, OL, C, name, pools, mods, mats in rows:
        if e < minE or tag in done:
            continue
        done.add(tag)
        bs = BStruct(N, mods, mats, tag)
        print(f"N = {N}, {tag}: O_A = {OA}, O_L = {OL}, classes = {C}, best estimate 10^{e:.2f}", flush=True)
        if search(bs, K, log, maxpairs):
            return True
    return False


def test():
    """N = 13, 17, 21, 25, 29 (G = Z_2^2 x Z_3 .. Z_7): every spin structure with a solution must verify."""
    for N in (13, 17, 21, 25, 29, 33, 37, 41, 45):
        ok = run_spin(N, maxo=40, maxpairs=1e9, minE=-99)
        print(f"N = {N}: {'found' if ok else 'none'}\n", flush=True)


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("test")
    p = sub.add_parser("spin"); p.add_argument("N", type=int, nargs="+"); p.add_argument("--maxo", type=int, default=50)
    p.add_argument("--maxpairs", type=float, default=2e10); p.add_argument("--log"); p.add_argument("--rank", action="store_true")
    p.add_argument("--minE", type=float, default=-6.0)
    a = ap.parse_args()
    if a.cmd == "test":
        test()
    else:
        for N in a.N:
            run_spin(N, a.maxo, a.maxpairs, a.log, a.rank, minE=a.minE)


if __name__ == "__main__":
    main()
