"""Skew Goethals-Seidel search on Paley-product structures for any spike configuration
(generalises paley_config.py, which handles N_w = 1 and 2).

v = m q with q = 3 (mod 4) prime, q not dividing m, H = {h = 1 (mod m), h a square (mod q)}. Over Z_m x Z_q
every H-invariant block is X = u (x) e + w (x) G + z (x) chi, and at each a in Z_m either
    z-type: w(a) = 0, u(a) = +-1, z(a) = +-1,   or   spike: z(a) = 0, w(a) = +-1, u(a) in {0, -2 w(a)}.
With y = u + q w (the compression to Z_m) the Goethals-Seidel condition is equivalent to
    (Y)  sum_k y_k y_k* = 4v                                   (off zero: sum_k PAF(y_k)(s) = 0)
    (Z1) sum_k z_k z_k* = sum_k (u_k w_k* + w_k u_k* + q w_k w_k*)
    (Z2) sum_k (z_k u_k* - u_k z_k*) = 0
in Z[Z_m]; (Y) is (I) + q (II) of the e/G/chi decomposition, (Z1) is (II) and (Z2) is (III). Block 0 is of skew
type: u_0(0) = 1, u_0 and w_0 antisymmetric, z_0 symmetric. The number of spikes N_w and of spikes with
u = -2w, N_-, satisfy (q + 1) N_w - 4 N_- = 4m.

A configuration fixes the spikes: (block, position, sign, u). Stage Y joins the y-lists of the four blocks
exactly (join4x) on the (m-1)/2 classes; stage Z takes each y-quadruple, computes the target of (Z1) and joins
the z-lists on PAF and cross classes. Lists larger than --lmax are uniform samples, so a run without hits
proves nothing; every hit is lifted to Z_v and verified as a skew-Hadamard matrix.

The configurations tried have sum_k w_k w_k* = N_w (no q^2 terms off zero): singletons at 0 in blocks 1-3,
pairs {0, d} in two blocks with opposite sign products, and an antipodal pair {a, -a} in block 0 with a pair
{0, 2a} of equal signs in another block.

  python3 search/pcg.py v [--nw N] [--lmax L] [--maxu U] [--seed S] [--stop]
      -> hits appended to $HADAMARD_WORK/skew_hits.jsonl
"""
import argparse, itertools, json, sys, time
import numpy as np
from pathlib import Path as _P
sys.path[:0] = [str(_P(__file__).resolve().parent)]
import workdir as W
from paley_config import pafs, accs, join, lift, check
from paley_product import structures
SKEWU = False
LOWE = 0


def divisor_reps(m):
    return [d for d in range(1, m) if m % d == 0]


def configs(m, q, nw):
    """Spike lists [(block, a, sign)] with sum_k w_k w_k* = nw at 0 and nothing off zero."""
    out = []
    if nw <= 3:
        out.append([(k, 0, 1) for k in range(1, nw + 1)])
    for d in divisor_reps(m):
        pp = [(1, 0, 1), (1, d, 1), (2, 0, 1), (2, d, -1)]                 # pair-pair, opposite products
        if nw == 4: out.append(pp)
        if nw == 5: out.append(pp + [(3, 0, 1)])
        if 2 * d % m:
            ap = [(0, d, 1), (0, (-d) % m, -1), (1, 0, 1), (1, 2 * d % m, 1)]   # antipodal pair + equal pair
            if nw == 4: out.append(ap)
            if nw == 5: out.append(ap + [(2, 0, 1)])
            if nw == 6: out.append(ap + [(2, 0, 1), (3, 0, 1)])
    return out


def with_minus(cfg, nm):
    """All placements of nm spikes with u = -2w (antipodal pairs in block 0 count twice and move together)."""
    units, used = [], set()
    for i, (k, a, s) in enumerate(cfg):
        if i in used: continue
        if k == 0:
            j = next(jj for jj, (kk, aa, ss) in enumerate(cfg) if kk == 0 and jj != i and jj not in used)
            units.append((i, j)); used |= {i, j}
        else:
            units.append((i,)); used.add(i)
    return [set(i for u in comb for i in u) for r in range(len(units) + 1)
            for comb in itertools.combinations(units, r) if sum(len(u) for u in comb) == nm]


def rows_pm(n, lmax, rng):
    if n <= 22 and 2 ** n <= lmax:
        return (1 - 2 * ((np.arange(2 ** n)[:, None] >> np.arange(n)) & 1)).astype(np.int8)
    S = rng.choice(np.array([-1, 1], np.int8), size=(int(lmax * 1.05) + 16, n))
    return np.unique(S, axis=0)[:lmax]


def ylist(m, q, k, spikes, minus, lmax, rng):
    """u- and w-arrays of block k for the configuration (u free +-1 off the spikes)."""
    h = (m - 1) // 2
    w = np.zeros(m, np.int64); u0 = np.zeros(m, np.int64); fixed = np.zeros(m, bool)
    for i, (kk, a, s) in enumerate(spikes):
        if kk != k: continue
        w[a] = s; u0[a] = -2 * s if i in minus else 0; fixed[a] = True
    if k == 0:
        free = [a for a in range(1, h + 1) if not fixed[a]]
        S = rows_pm(len(free), lmax, rng)
        U = np.zeros((len(S), m), np.int64); U[:, 0] = 1
        U[:, fixed] = u0[fixed]
        for j, a in enumerate(free): U[:, a] = S[:, j]; U[:, m - a] = -S[:, j]
    elif SKEWU and fixed.sum() == 1 and fixed[0]:
        S = rows_pm(h, lmax, rng)                       # u skew-type about the spike: no q-layer terms
        U = np.zeros((len(S), m), np.int64); U[:, 0] = u0[0]
        for j in range(h): U[:, j + 1] = S[:, j]; U[:, m - 1 - j] = -S[:, j]
    else:
        free = [a for a in range(m) if not fixed[a]]
        S = rows_pm(len(free), lmax, rng)
        U = np.zeros((len(S), m), np.int64); U[:, fixed] = u0[fixed]; U[:, free] = S
    return U, w


def lowe(Z, free, steps):
    """Greedy descent of sum_s PAF(s)^2 from the rows of Z, flipping one free position per row and step."""
    m = Z.shape[1]; h = (m - 1) // 2; sh = np.arange(1, h + 1)
    Z = Z.astype(np.int32)
    for _ in range(steps):
        P = np.stack([(Z * np.roll(Z, -s, axis=1)).sum(1) for s in sh], 1)
        E = (P * P).sum(1); best = np.zeros(len(Z), np.int64); bi = np.full(len(Z), -1)
        for i in free:
            d = -2 * Z[:, i:i + 1] * (Z[:, (i + sh) % m] + Z[:, (i - sh) % m])
            g = E - ((P + d) ** 2).sum(1)
            up = g > best; best[up] = g[up]; bi[up] = i
        r = np.flatnonzero(bi >= 0)
        if not len(r): break
        Z[r, bi[r]] *= -1
    return np.unique(Z.astype(np.int64), axis=0)


def zlist(m, k, w, lmax, rng):
    h = (m - 1) // 2
    spike = w != 0
    if k == 0:
        free = [a for a in range(1, h + 1) if not spike[a]]
        S = rows_pm(len(free) + 1, lmax, rng)
        Z = np.zeros((len(S), m), np.int64); Z[:, 0] = S[:, 0]
        for j, a in enumerate(free): Z[:, a] = S[:, j + 1]; Z[:, m - a] = S[:, j + 1]
    else:
        free = [a for a in range(m) if not spike[a]]
        S = rows_pm(len(free), lmax, rng)
        Z = np.zeros((len(S), m), np.int64); Z[:, free] = S
        if LOWE and 2 ** len(free) > lmax: Z = lowe(Z, free, LOWE)
    return Z


def corr(x, y, h):
    """(x y*)(s) = sum_b x(b + s) y(b), s = 1..h."""
    return np.array([(np.roll(x, -s) * y).sum() for s in range(1, h + 1)], np.int64)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('v', type=int)
    ap.add_argument('--nw', type=int, default=0, help='only this N_w (default: every admissible N_w <= 6)')
    ap.add_argument('--lmax', type=int, default=1 << 17)
    ap.add_argument('--lzmax', type=int, default=0, help='z-list cap (default: --lmax)')
    ap.add_argument('--maxu', type=int, default=30, help='y-quadruples tried per configuration')
    ap.add_argument('--seed', type=int, default=1)
    ap.add_argument('--stop', action='store_true')
    ap.add_argument('--skewu', action='store_true', help='u skew-type about a lone spike at 0')
    ap.add_argument('--lowe', type=int, default=0, help='greedy low-PAF steps for sampled z-lists')
    a = ap.parse_args()
    v = a.v; rng = np.random.default_rng(a.seed); lz = a.lzmax or a.lmax
    global SKEWU, LOWE
    SKEWU, LOWE = a.skewu, a.lowe
    log = lambda s: print(s, flush=True)
    found = 0
    for m, q, hgen in structures(v, 10 ** 6):
        h = (m - 1) // 2
        for nw in range(1, 7):
            nm4 = (q + 1) * nw - 4 * m
            if nm4 < 0 or nm4 % 4 or nm4 // 4 > nw or (a.nw and nw != a.nw): continue
            nm = nm4 // 4
            for cfg in configs(m, q, nw):
                for minus in with_minus(cfg, nm):
                    t0 = time.time()
                    Y = [ylist(m, q, k, cfg, minus, a.lmax, rng) for k in range(4)]
                    P = [pafs(U + q * w[None, :], h) for U, w in Y]
                    ty, sy = join(P, h, K=a.maxu)
                    desc = ' '.join(f"{'ABCD'[k]}{x}{'+' if s > 0 else '-'}{'*' if i in minus else ''}" for i, (k, x, s) in enumerate(cfg))
                    log(f"v={v} m={m} q={q} N_w={nw} N_-={nm} [{desc}]: y-lists {[len(U) for U, _ in Y]}, y-quadruples {ty} [{time.time() - t0:.0f}s]")
                    for qi, sol in enumerate(sy):
                        U = [Y[k][0][sol[k]] for k in range(4)]; Wv = [Y[k][1] for k in range(4)]
                        R = sum(corr(U[k], Wv[k], h) + corr(Wv[k], U[k], h) + q * corr(Wv[k], Wv[k], h) for k in range(4))
                        Zl = [zlist(m, k, Wv[k], lz, rng) for k in range(4)]
                        V = [np.concatenate([pafs(Z, h), accs(Z, U[k][None, :], h)], axis=1) for k, Z in enumerate(Zl)]
                        V[0] = V[0].copy(); V[0][:, :h] -= R[None, :]
                        tz, sz = join(V, 2 * h, K=4)
                        hits = 0
                        for s4 in sz:
                            Zs = [Zl[k][s4[k]] for k in range(4)]
                            seqs = lift(v, m, q, U, Wv, Zs)
                            if check(v, seqs):
                                hits += 1; found += 1
                                rec = dict(n=v, gens=[hgen], pattern=[int(x.sum()) for x in seqs], hadamard=True, skew=True,
                                           method=f'pcg N_w={nw} [{desc}]', seqs=[''.join('+' if x > 0 else '-' for x in sq) for sq in seqs])
                                with open(W.hits(), 'a') as f: f.write(json.dumps(rec) + '\n')
                        log(f"  y-quadruple {qi}: z-lists {[len(Z) for Z in Zl]}, z-solutions {tz}, verified {hits} [{time.time() - t0:.0f}s]")
                        if hits and a.stop:
                            log(f"v={v}: HIT"); return
    log(f"v={v}: done, {found} verified hits")


if __name__ == '__main__':
    main()
