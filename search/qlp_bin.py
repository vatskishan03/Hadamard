#!/usr/bin/env python3
"""Quaternary Legendre pairs (A, B) of even length l with A binary and B reversible (search/NOTES.md; the shape of
Verma's length-64 pair, arXiv:2609.25573, where A is binary).

A qLP is a pair A, B in {1, -1, i, -i}^l with PAF_A(k) + PAF_B(k) = -2 for k != 0 (Kotsireas-Winterhof). Write
B = ((r + s) + i(r - s))/2 with r, s in {+-1}^l (Gray map, Lebedev arXiv:2609.04589 (5)). If r and s are symmetric
(x_j = x_-j) the cross term of PAF_B vanishes (Lebedev's Lemma 3) and PAF_B = (PAF_r + PAF_s)/2. With A binary the
condition becomes

    2 PAF_A(k) + PAF_r(k) + PAF_s(k) = -4        (k != 0),

and summing over k gives 2 r_A^2 + r_r^2 + r_s^2 = 4: r_A = 0 and (r_r, r_s) = (2, 0) up to sign and order, which is
the normalization sum A = 0, sum B = 1 + i. A is not symmetric in general, so this is outside the reversible family
that Lebedev searched exhaustively at lengths 64, 70 and 76.

Search: A constant on the orbits of a multiplier group H <= Z_l^*, r and s on the orbits of <H, -1>; the condition
holds on the <H, -1>-classes. The (r, s) pairs are tabulated (all of them) and the A-list is probed against the
table (exact, by a numpy sort join).

    python3 search/qlp_bin.py triage l [--maxa 3e8]
    python3 search/qlp_bin.py search l g [g2] [--log FILE]
"""
import argparse, json, math, os, struct, subprocess, sys, time
from math import gcd, log10, pi
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import workdir as W
from bnw import Struct, subgroups, mask_to_y


def gauss(st, K=3000, seed=1):
    """log10 Gaussian estimate for rows (r_A, r_r, r_s) = (0, 2, 0); A weighted by 2."""
    rng = np.random.default_rng(seed)
    pools = [st.counts("A").get(0, 0) // 2, st.counts("L").get(2, 0), st.counts("L").get(0, 0) // 2]
    if min(pools) == 0:
        return -math.inf, pools
    lN = sum(log10(p) for p in pools)
    mus, vars_, lats = [], [], []
    for (kind, r), w in zip((("A", 0), ("L", 2), ("L", 0)), (2, 1, 1)):
        V = w * st.values(kind, st.sample(kind, r, K, rng))
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
            if abs(-4 - mu) > 1e-6:
                return -math.inf, pools
            continue
        if g and (-4 - res) % g:
            return -math.inf, pools
        if not dropped:
            dropped = True; continue
        lp += -0.5 * log10(2 * pi * var) - (-4 - mu) ** 2 / (2 * var) / math.log(10) + log10(max(g, 1))
    return lN + lp, pools


def triage(l, maxa=3e8):
    rows = []
    for H, gens in subgroups(l).items():
        if (l - 1) in H:
            continue
        st = Struct(l, gens, plain=True)
        e, pools = gauss(st)
        if e == -math.inf:
            continue
        rows.append((e, len(H), gens, len(st.sizesA), len(st.sizesL), st.C, pools))
    rows.sort(key=lambda r: -r[0])
    print(f"l = {l}: {len(rows)} structures (A binary on H-orbits, r, s symmetric)")
    for e, h, gens, OA, OL, C, pools in rows[:20]:
        flag = "" if pools[0] <= maxa and pools[1] * pools[2] <= 3e8 else "  (too large)"
        print(f"{e:8.2f} |H|={h:3d} {str(gens):>10} O_A={OA:3d} O_L={OL:3d} cls={C:3d} pools={pools}{flag}")
    return rows


def enum(st, kind, r, fix0):
    d = W.WORK / "qlp"; d.mkdir(exist_ok=True)
    base = d / f"q{st.N}_{'_'.join(map(str, st.gens))}_{kind}_{r}_{fix0}"
    vals, masks = base.with_suffix(".v"), base.with_suffix(".m")
    if not vals.exists():
        inp = base.with_suffix(".in")
        sz = st.sizes(kind); S = st.S[kind]
        with open(inp, "wb") as f:
            f.write(struct.pack("4i", len(sz), st.C, r, fix0))
            f.write(np.asarray(sz, np.int32).tobytes()); f.write(np.asarray(S, np.int32).tobytes())
        subprocess.run([W.bin("bnw_enum"), str(inp), str(vals) + ".tmp", str(masks) + ".tmp"], check=True,
                       capture_output=True)
        os.replace(str(vals) + ".tmp", vals); os.replace(str(masks) + ".tmp", masks)
    return np.fromfile(vals, np.int16).reshape(-1, st.C), np.fromfile(masks, np.uint64)


def paf_c(x):
    X = np.fft.fft(np.asarray(x, complex))
    return np.fft.ifft(X * np.conj(X))


def check_qlp(A, B):
    l = len(A)
    P = paf_c(A) + paf_c(B)
    ok = np.allclose(P[1:], -2) and abs(P[0] - 2 * l) < 1e-6
    # exact: integer arithmetic over Z[i]
    a = np.array([[int(round(z.real)), int(round(z.imag))] for z in A]); b = np.array([[int(round(z.real)), int(round(z.imag))] for z in B])
    for k in range(1, l):
        re = im = 0
        for X in (a, b):
            Y = np.roll(X, -k, axis=0)          # X_{j+k}
            re += int((X[:, 0] * Y[:, 0] + X[:, 1] * Y[:, 1]).sum())       # X_j conj(X_{j+k})
            im += int((X[:, 1] * Y[:, 0] - X[:, 0] * Y[:, 1]).sum())
        if (re, im) != (-2, 0):
            return False
    return ok


def search(l, gens, log=None, K=5):
    st = Struct(l, gens, plain=True)
    print(f"l = {l}, H = <{','.join(map(str, gens))}> of order {len(st.H)}: O_A = {len(st.sizesA)}, "
          f"O_L = {len(st.sizesL)}, classes = {st.C}; estimate 10^{gauss(st)[0]:.2f}", flush=True)
    t0 = time.time()
    VA, MA = enum(st, "A", 0, 1)
    VR, MR = enum(st, "L", 2, 0)
    VS, MS = enum(st, "L", 0, 1)
    print(f"  lists A {len(VA)}, r {len(VR)}, s {len(VS)} ({time.time() - t0:.0f}s)", flush=True)
    # table of r + s values
    RS = (VR[:, None, :].astype(np.int32) + VS[None, :, :]).reshape(-1, st.C)
    key = lambda V: np.ascontiguousarray(V.astype(np.int16)).view(np.dtype((np.void, 2 * st.C))).ravel()
    krs = key(RS); order = np.argsort(krs); krs = krs[order]
    found = []
    step = 2_000_000
    for s0 in range(0, len(VA), step):
        tgt = (-4 - 2 * VA[s0:s0 + step].astype(np.int32))
        kt = key(tgt)
        pos = np.searchsorted(krs, kt)
        pos[pos >= len(krs)] = 0
        hit = np.flatnonzero(krs[pos] == kt)
        for h in hit:
            a = s0 + h; rsi = order[pos[h]]; ri, si = divmod(int(rsi), len(VS))
            yA = mask_to_y(MA[a], len(st.sizesA)); yr = mask_to_y(MR[ri], len(st.sizesL)); ys = mask_to_y(MS[si], len(st.sizesL))
            A = st.expand("A", yA).astype(complex); r = st.expand("L", yr); s = st.expand("L", ys)
            B = ((r + s) + 1j * (r - s)) / 2
            if check_qlp(A, B):
                rec = {"l": l, "gens": list(gens), "A": "".join("0" if v > 0 else "2" for v in A.real),
                       "B": "".join({1: "0", 1j: "1", -1: "2", -1j: "3"}[complex(z)] for z in B)}
                found.append(rec)
                print(f"  qLP of length {l} found and verified exactly", flush=True)
                if log:
                    with open(log, "a") as f:
                        f.write(json.dumps(rec) + "\n")
                if len(found) >= K:
                    return found
    print(f"  {len(found)} pairs ({time.time() - t0:.0f}s)", flush=True)
    return found


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("triage"); p.add_argument("l", type=int, nargs="+"); p.add_argument("--maxa", type=float, default=3e8)
    p = sub.add_parser("search"); p.add_argument("l", type=int); p.add_argument("g", type=int, nargs="+")
    p.add_argument("--log")
    a = ap.parse_args()
    if a.cmd == "triage":
        for l in a.l:
            triage(l, a.maxa)
    else:
        search(a.l, a.g, a.log)


if __name__ == "__main__":
    main()
