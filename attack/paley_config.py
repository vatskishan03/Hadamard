"""Exact search for skew Goethals-Seidel families on Paley-product structures through configurations over Z_m.

Structure (see paley_product.py): v = m q, q = 3 (mod 4) prime, H = {h = 1 (mod m), h a square (mod q)}. Every block is
    X = u (x) 1 + w (x) G + z (x) chi        (G = all-ones and chi = Legendre symbol on Z_q, u, w, z on Z_m)
where at each a in Z_m either w(a) = 0, u(a), z(a) = +-1 ("z-type"), or z(a) = 0, w(a) = +-1, u(a) in {0, -2 w(a)}.
With N_w = number of w-type positions in all four blocks and N_- = those with u = -2w, the identity coefficient of
the G-equation gives (q + 1) N_w - 4 N_- = 4m.

For N_w = 1 (q = 4m - 1 or 4m + 3) and N_w = 2 with the two w-positions in different blocks (q = 2m + 1, N_- = 1),
the size bound |sum_k PAF_u(s)| <= 4m < 2q forces the decoupled system (all over Z_m, s != 0):
    sum_k PAF_{u_k}(s) = 0,   sum_k PAF_{z_k}(s) = 0,   sum_k [z_k * u_k~ - u_k * z_k~](s) = 0,
with u_A skew-type, z_A symmetric (block A has no w-type position), the w-positions translated to a = 0, and
    N_w = 1: u_j skew-type away from 0, u_j(0) = 0 (q = 4m - 1) or -2 (q = 4m + 3);
    N_w = 2: u_1(0) = 0, u_2(0) = -2 and u_1(s) + u_1(-s) = -(u_2(s) + u_2(-s)).
N_w = 1 also forces r_B^2 + r_C^2 = 4m - 2 for the row sums of the two free blocks, so 2m - 1 must be a sum of two
squares (false for m = 11, 17: no solutions for v = 473, 517, 1139, 1207).

Stage 1 lists every u-quadruple (exact join); stage 2 joins, for each u-quadruple, the z-lists on the 2(m-1)/2
classes (PAF_z and the cross terms) with join4x. Every hit is rebuilt over Z_v and verified as a skew-Hadamard
matrix. With --all, every u-quadruple is processed (an exhaustive search of the configuration space).

  python3 attack/paley_config.py v [--maxu N] [--all]       -> hits in $HADAMARD_WORK/skew_hits.jsonl
"""
import argparse, itertools, json, math, os, struct, subprocess, sys, time
import numpy as np
from pathlib import Path as _P
sys.path[:0] = [str(_P(__file__).resolve().parent), str(_P(__file__).resolve().parent.parent / 'verify')]
import workdir as W
import verify
from paley_product import structures


def pafs(X, h):
    """PAF(s), s = 1..h, of the rows of X (n x m)."""
    return np.stack([(X * np.roll(X, -s, axis=1)).sum(1) for s in range(1, h + 1)], axis=1)


def accs(Z, U, h):
    """sum_a z(a+s) u(a) - u(a+s) z(a), s = 1..h (U broadcast against Z)."""
    return np.stack([(np.roll(Z, -s, axis=1) * U - np.roll(U, -s, axis=1) * Z).sum(1) for s in range(1, h + 1)], axis=1)


def allpm(n):
    return (1 - 2 * ((np.arange(2 ** n)[:, None] >> np.arange(n)) & 1)).astype(np.int8)


def skew_type(m, x0=1):
    h = (m - 1) // 2
    S = allpm(h); X = np.zeros((len(S), m), np.int8); X[:, 0] = x0
    for a in range(1, h + 1): X[:, a] = S[:, a - 1]; X[:, m - a] = -S[:, a - 1]
    return X


def symmetric(m, zero_at_0=False):
    h = (m - 1) // 2
    S = allpm(h + (0 if zero_at_0 else 1)); X = np.zeros((len(S), m), np.int8)
    if zero_at_0:
        for a in range(1, h + 1): X[:, a] = S[:, a - 1]; X[:, m - a] = S[:, a - 1]
    else:
        X[:, 0] = S[:, 0]
        for a in range(1, h + 1): X[:, a] = S[:, a]; X[:, m - a] = S[:, a]
    return X


def free(m, zero_at_0=False, at0=None):
    if zero_at_0 or at0 is not None:
        S = allpm(m - 1); X = np.zeros((len(S), m), np.int8); X[:, 1:] = S
        X[:, 0] = 0 if at0 is None else at0
        return X
    return allpm(m)


def canonical_translates(X):
    """One representative (lexicographically smallest translate) per translation class of the rows of X."""
    m = X.shape[1]
    best = X.copy()
    for t in range(1, m):
        R = np.roll(X, -t, axis=1)
        key_b = best.astype(np.int64) @ (3 ** np.arange(m - 1, -1, -1, dtype=np.int64))
        key_r = R.astype(np.int64) @ (3 ** np.arange(m - 1, -1, -1, dtype=np.int64))
        best = np.where((key_r < key_b)[:, None], R, best)
    return np.unique(best, axis=0)


def join(lists, C, K, budget=8000000):
    """Exact 4-list join (join4x) of int16 class vectors to target 0; returns (total, [(i,j,k,l), ...])."""
    sz = [len(x) for x in lists]
    if min(sz) == 0: return 0, []
    order = np.argsort(sz); perm = [order[0], order[3], order[1], order[2]]
    path = os.path.join(W.tmpdir(), f'pc_{os.getpid()}.bin')
    with open(path, 'wb') as f:
        f.write(struct.pack('iiiii', *[sz[k] for k in perm], C)); f.write(np.zeros(C, np.int16).tobytes())
        for k in perm: f.write(np.ascontiguousarray(lists[k].astype(np.int16)).tobytes())
    r = subprocess.run([W.bin('join4x'), path, str(K), str(budget)], capture_output=True, text=True, check=True)
    os.remove(path)
    total = int(r.stdout.split('TOTAL')[-1])
    sols = []
    for line in r.stdout.splitlines():
        if line.startswith('SOL'):
            ids = list(map(int, line.split()[1:])); q = [None] * 4
            for slot, k in enumerate(perm): q[k] = ids[slot]
            sols.append(tuple(q))
    return total, sols


def lift(v, m, q, U, Wv, Z):
    """Sequences over Z_v from (u, w, z) over Z_m: X(a, b) = u(a)[b=0] + w(a) + z(a) chi(b)."""
    chi = np.zeros(q, np.int64)
    for x in range(1, q): chi[x * x % q] = 1
    chi[(chi == 0) & (np.arange(q) != 0)] = -1
    t = np.arange(v); a, b = t % m, t % q
    return [(u[a] * (b == 0) + w[a] + z[a] * chi[b]).astype(np.int64) for u, w, z in zip(U, Wv, Z)]


def check(v, seqs):
    tot = verify.paf_sum(seqs)
    if not (tot[0] == 4 * v and not tot[1:].any()): return False
    Hm = verify.goethals_seidel(*map(verify.circulant, seqs)).astype(np.float64)
    return np.array_equal(Hm @ Hm.T, 4 * v * np.eye(4 * v)) and np.array_equal(Hm + Hm.T, 2 * np.eye(4 * v))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('v', type=int); ap.add_argument('--maxu', type=int, default=200)
    ap.add_argument('--all', action='store_true'); ap.add_argument('--nsol', type=int, default=2)
    a = ap.parse_args()
    v = a.v
    log = lambda s: print(s, flush=True)
    cases = []
    for m, q, hgen in structures(v, 10 ** 6):
        for Nw in (1, 2):
            Nm = ((q + 1) * Nw - 4 * m)
            if Nm % 4 == 0 and 0 <= Nm // 4 <= Nw: cases.append((m, q, hgen, Nw, Nm // 4))
    if not cases: sys.exit(f"{v}: no Paley-product structure with N_w = 1 or 2 ((q+1)N_w - 4N_- = 4m has no solution)")
    for m, q, hgen, Nw, Nm in cases:
        h = (m - 1) // 2
        t0 = time.time()
        UA = skew_type(m)
        if Nw == 1:
            if not any(x * x + y * y == 4 * m - 2 for x in range(1, 2 * m, 2) for y in range(1, 2 * m, 2)):
                log(f"v={v} m={m} q={q} N_w=1: 4m - 2 = {4 * m - 2} is not a sum of two odd squares, so no solution"); continue
            U1 = skew_type(m, x0=(0 if Nm == 0 else -2))
            U2 = free(m); U3 = canonical_translates(free(m))
            P = [pafs(UA, h), pafs(U1, h), pafs(U2, h), pafs(U3, h)]
            lists_u = [UA, U1, U2, U3]
        else:
            # joint (u1, u2): u1(0) = 0, u2(0) = -2, symmetric parts opposite
            # per pair {a, -a}: u1 = (s, s) forces u2 = (-s, -s); u1 = (s, -s) leaves u2 = (t, -t) free: 6 options
            opts = []
            for p1, p2 in itertools.product((1, -1), repeat=2):
                for c2 in ([(-p1, -p1)] if p1 == p2 else [(1, -1), (-1, 1)]): opts.append(((p1, p2), c2))
            if h > 8: sys.exit("N_w = 2 stage 1 is materialised in memory; m <= 17 only")
            idx = np.array(list(itertools.product(range(6), repeat=h)), dtype=np.int64)
            U1 = np.zeros((len(idx), m), np.int8); U2 = np.zeros((len(idx), m), np.int8); U2[:, 0] = -2
            for p in range(h):
                o = np.array(opts, dtype=np.int8)[idx[:, p]]          # (n, 2, 2)
                U1[:, p + 1] = o[:, 0, 0]; U1[:, m - 1 - p] = o[:, 0, 1]
                U2[:, p + 1] = o[:, 1, 0]; U2[:, m - 1 - p] = o[:, 1, 1]
            U12 = np.concatenate([U1, U2], axis=1)
            U3 = canonical_translates(free(m))
            P = [pafs(UA, h), pafs(U1, h) + pafs(U2, h), pafs(U3, h), np.zeros((1, h), np.int64)]
            lists_u = [UA, U12, U3, np.zeros((1, m), np.int8)]
        total_u, sols_u = join(P, h, K=10 ** 7)
        log(f"v={v} m={m} q={q} N_w={Nw} N_-={Nm}: stage 1 lists {[len(x) for x in lists_u]}, u-quadruples {total_u} [{time.time()-t0:.0f}s]")
        if not total_u: continue
        ZA = symmetric(m); Z0 = free(m, zero_at_0=True); Zf = free(m)
        found = 0; tried = 0
        for qi, sol in enumerate(sols_u):
            if not a.all and tried >= a.maxu: break
            tried += 1
            if Nw == 1:
                u = [lists_u[k][sol[k]].astype(np.int64) for k in range(4)]
                zl = [ZA, Z0, Zf, Zf]
            else:
                u12 = lists_u[1][sol[1]]
                u = [lists_u[0][sol[0]].astype(np.int64), u12[:m].astype(np.int64), u12[m:].astype(np.int64),
                     lists_u[2][sol[2]].astype(np.int64)]
                zl = [ZA, Z0, Z0, Zf]
            V = [np.concatenate([pafs(Z.astype(np.int64), h), accs(Z.astype(np.int64), u[k][None, :], h)], axis=1)
                 for k, Z in enumerate(zl)]
            tz, sz = join(V, 2 * h, K=a.nsol)
            if tz:
                Wv = [np.zeros(m, np.int64) for _ in range(4)]
                Wv[1][0] = 1
                if Nw == 2: Wv[2][0] = 1
                for s in sz:
                    Z = [zl[k][s[k]].astype(np.int64) for k in range(4)]
                    seqs = lift(v, m, q, u, Wv, Z)
                    ok = check(v, seqs)
                    log(f"  u-quadruple {qi}: {tz} z-solutions; lifted to Z_{v}: skew-Hadamard {ok}")
                    if ok:
                        found += 1
                        rec = dict(n=v, gens=[hgen], pattern=[int(x.sum()) for x in seqs], hadamard=True, skew=True,
                                   method=f'paley_config N_w={Nw}', seqs=[''.join('+' if x > 0 else '-' for x in sq) for sq in seqs])
                        with open(W.hits(), 'a') as f: f.write(json.dumps(rec) + '\n')
                if found and not a.all: break
            if qi < 5 or qi % 50 == 0:
                log(f"  u-quadruple {qi}: z-solutions {tz} [{time.time()-t0:.0f}s]")
        log(f"v={v} N_w={Nw}: {tried} of {total_u} u-quadruples tried, {found} verified hits [{time.time()-t0:.0f}s]")


if __name__ == '__main__':
    main()
