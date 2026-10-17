#!/usr/bin/env python3
"""RRR / difference-map search (Elser) for four +-1 sequences with a flat summed power spectrum.

Problems (all on Z_m, four sequences x_1..x_4):
  gs      : sum_k |X_k(w)|^2 = 4m for w != 0                       (Goethals-Seidel quadruple, H(4m))
  sz      : as gs, x_1, x_2 symmetric                                (Szöllősi, complex Hadamard 2m)
  spence  : sum_k |X_k(w)|^2 = 4(m+1) for w != 0, row sums (2,0,0,0) (Spence, H(4(m+1)), m even)
  bsz     : as spence, x_1, x_2 symmetric                             (bordered Szöllősi, complex Hadamard 2(m+1))
Constraint A: each x_k in {+-1}^m with prescribed row sum r_k (symmetric blocks: symmetric), projection = keep the
r-dependent number of largest entries (+1) after symmetrising. Constraint B: per frequency w != 0 the vector
(X_1(w), .., X_4(w)) in C^4 lies on the sphere of radius sqrt(target) (projection = rescaling); at w = 0 the row sums
are left to A. Iteration (RRR, beta): x <- x + beta (P_A(2 P_B(x) - x) - P_B(x)). A solution is detected when
P_A(x) satisfies the exact integer condition.
"""
import argparse, math, sys, time
import numpy as np


def proj_A(x, rows, symm):
    """nearest +-1 vectors with the given row sums (symmetric blocks: on orbit pairs {t, -t})"""
    m = x.shape[1]
    out = np.empty_like(x)
    for k in range(4):
        v = x[k]
        r = rows[k]
        if symm[k]:
            w = (v + v[(-np.arange(m)) % m]) / 2
            # orbits: {0} and pairs {t, m - t} (and {m/2} if m even)
            reps = [0] + [t for t in range(1, m) if t < m - t] + ([m // 2] if m % 2 == 0 else [])
            sizes = np.array([1] + [2] * (len([t for t in range(1, m) if t < m - t])) + ([1] if m % 2 == 0 else []))
            vals = w[reps]
            # choose signs s_o in {+-1} with sum s_o * size_o = r maximising sum s_o * size_o * vals_o:
            # greedy by value per orbit with a small DP over sizes (sizes are 1 or 2)
            n = len(reps)
            # number of + orbits of each size: sum(size * s) = r  ->  2 * sum_{+} size - total = r
            total = sizes.sum(); need = (r + total) // 2           # sum of sizes of the + orbits
            o1 = [i for i in range(n) if sizes[i] == 1]; o2 = [i for i in range(n) if sizes[i] == 2]
            o1.sort(key=lambda i: -vals[i]); o2.sort(key=lambda i: -vals[i])
            best, arg = -1e300, None
            for a in range(len(o1) + 1):                         # a size-1 orbits positive
                b2 = need - a
                if b2 < 0 or b2 % 2 or b2 // 2 > len(o2):
                    continue
                b = b2 // 2
                sc = vals[o1[:a]].sum() - vals[o1[a:]].sum() + 2 * (vals[o2[:b]].sum() - vals[o2[b:]].sum())
                if sc > best:
                    best, arg = sc, (a, b)
            s = -np.ones(n)
            if arg is None:
                s = np.sign(vals + 1e-12)
            else:
                a, b = arg
                s[o1[:a]] = 1; s[o2[:b]] = 1
            y = np.empty(m)
            for i, t in enumerate(reps):
                y[t] = s[i]; y[(-t) % m] = s[i]
            out[k] = y
        else:
            npos = (m + r) // 2
            idx = np.argsort(-v)
            y = -np.ones(m); y[idx[:npos]] = 1
            out[k] = y
    return out


def proj_B(x, target):
    X = np.fft.fft(x, axis=1)
    nrm = np.sqrt((np.abs(X[:, 1:]) ** 2).sum(0))
    scale = np.where(nrm > 1e-12, math.sqrt(target) / np.maximum(nrm, 1e-12), 1.0)
    X[:, 1:] *= scale
    return np.fft.ifft(X, axis=1).real


def check(y, target_off):
    s = sum(np.rint(np.fft.ifft(np.abs(np.fft.fft(v)) ** 2).real) for v in y)
    return bool(np.all(s[1:] == target_off)), s


def run(m, kind, rows, beta=0.5, iters=10 ** 6, seed=1, report=10000, timeout=None):
    symm = [kind in ("sz", "bsz")] * 2 + [False, False]
    target = 4 * m if kind in ("gs", "sz") else 4 * (m + 1)
    target_off = 0 if kind in ("gs", "sz") else -4
    rng = np.random.default_rng(seed)
    x = rng.standard_normal((4, m))
    t0 = time.time()
    for it in range(iters):
        pb = proj_B(x, target)
        pa = proj_A(2 * pb - x, rows, symm)
        x = x + beta * (pa - pb)
        if it % 10 == 0:
            y = proj_A(x, rows, symm)
            ok, s = check(y, target_off)
            if ok:
                return it, y, time.time() - t0
        if report and it % report == 0:
            y = proj_A(x, rows, symm)
            ok, s = check(y, target_off)
            err = float(((s[1:] - target_off) ** 2).sum())
            print(f"  it {it}: error {err:.0f} [{time.time() - t0:.0f}s]", flush=True)
        if timeout and time.time() - t0 > timeout:
            break
    return None, None, time.time() - t0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("m", type=int)
    ap.add_argument("kind", choices=["gs", "sz", "spence", "bsz"])
    ap.add_argument("--rows", required=True)
    ap.add_argument("--beta", type=float, default=0.5)
    ap.add_argument("--iters", type=int, default=10 ** 6)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--report", type=int, default=0)
    ap.add_argument("--timeout", type=float)
    a = ap.parse_args()
    rows = [int(v) for v in a.rows.split(",")]
    it, y, sec = run(a.m, a.kind, rows, a.beta, a.iters, a.seed, a.report, a.timeout)
    if it is None:
        print(f"m = {a.m} {a.kind} rows {rows}: no solution ({sec:.0f}s)")
    else:
        print(f"m = {a.m} {a.kind} rows {rows}: SOLUTION at iteration {it} ({sec:.1f}s)")
        for v in y:
            print("".join("+" if t > 0 else "-" for t in v))


if __name__ == "__main__":
    main()
