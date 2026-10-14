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

  python3 search/paley_config.py v [--maxu N] [--all]       -> hits in $HADAMARD_WORK/skew_hits.jsonl
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


def classvec(Z, u, h, chunk=1 << 20):
    """[PAF_z(s), cross(z, u)(s)] for s = 1..h, as int16, computed in chunks (lists of 2^23 rows stay small)."""
    out = np.empty((len(Z), 2 * h), np.int16)
    for c0 in range(0, len(Z), chunk):
        X = Z[c0:c0 + chunk].astype(np.int64)
        out[c0:c0 + chunk, :h] = pafs(X, h)
        out[c0:c0 + chunk, h:] = accs(X, u[None, :], h)
    return out


def allpm(n):
    """All 2^n sign vectors of length n (row i has -1 where bit j of i is set), built column by column."""
    idx = np.arange(2 ** n, dtype=np.uint32)
    out = np.empty((2 ** n, n), np.int8)
    for j in range(n): out[:, j] = 1 - 2 * ((idx >> np.uint32(j)) & np.uint32(1)).astype(np.int8)
    return out


class Orbits:
    """Orbits of a multiplier group K <= Z_m^* (K = {1} by default) on Z_m, and generators of K-invariant sequences.
    A skew-type sequence needs -a outside K a for every a != 0; pairs lists the orbit pairs {O, -O}."""
    def __init__(self, m, K=(1,)):
        self.m = m; K = sorted(set(K))
        oid = -np.ones(m, np.int64); reps = []
        for a in range(m):
            if oid[a] < 0:
                for k in K: oid[a * k % m] = len(reps)
                reps.append(a)
        self.oid, self.reps, self.n = oid, reps, len(reps)
        neg = [int(oid[(-r) % m]) for r in reps]
        self.selfneg = [i for i in range(1, self.n) if neg[i] == i]
        self.pairs = [(i, neg[i]) for i in range(1, self.n) if i < neg[i]]
        self.trivial = len(K) == 1

    def expand(self, V):                          # values on orbits (n_rows x n) -> sequences (n_rows x m)
        return V[:, self.oid]

    def skew_type(self, x0=1):
        if self.selfneg: return np.zeros((0, self.m), np.int8)
        S = allpm(len(self.pairs)); V = np.zeros((len(S), self.n), np.int8); V[:, 0] = x0
        for j, (i, k) in enumerate(self.pairs): V[:, i] = S[:, j]; V[:, k] = -S[:, j]
        return self.expand(V)

    def symmetric(self, zero_at_0=False):
        free_ = [i for i, _ in self.pairs] + self.selfneg
        S = allpm(len(free_) + (0 if zero_at_0 else 1)); V = np.zeros((len(S), self.n), np.int8)
        off = 0
        if not zero_at_0: V[:, 0] = S[:, 0]; off = 1
        for j, i in enumerate(free_): V[:, i] = S[:, j + off]
        for i, k in self.pairs: V[:, k] = V[:, i]
        return self.expand(V)

    def free(self, zero_at_0=False, at0=None):
        S = allpm(self.n - 1); V = np.zeros((len(S), self.n), np.int16); V[:, 1:] = S
        V[:, 0] = 0 if zero_at_0 else (at0 if at0 is not None else 1)
        X = self.expand(V)
        if zero_at_0 or at0 is not None: return X
        return np.concatenate([X, -X])            # both values at 0

    def fixed(self, base, pos_vals, skew=False, sym=False, zero_pos=()):
        """K-invariant sequences of the given kind ('free', or skew-type with skew=True, or symmetric with sym=True)
        with prescribed values at singleton-orbit positions pos_vals {a: value} and zeros at zero_pos."""
        fixed_orb = {}
        for a_, val in list(pos_vals.items()) + [(z_, 0) for z_ in zero_pos]:
            i = int(self.oid[a_ % self.m])
            if self.reps[i] != a_ % self.m and np.count_nonzero(self.oid == i) > 1: return np.zeros((0, self.m), np.int16)
            fixed_orb[i] = val
        if skew:
            free_ = []
            for i, k in self.pairs:
                if i in fixed_orb or k in fixed_orb:
                    vi = fixed_orb.get(i, None); vk = fixed_orb.get(k, None)
                    if vi is not None and vk is not None and vi != -vk: return np.zeros((0, self.m), np.int16)
                    if vi is None: fixed_orb[i] = -vk
                    if vk is None: fixed_orb[k] = -fixed_orb[i]
                else: free_.append((i, k))
            S = allpm(len(free_)); V = np.zeros((len(S), self.n), np.int16); V[:, 0] = fixed_orb.get(0, 1)
            for i, val in fixed_orb.items(): V[:, i] = val
            for j, (i, k) in enumerate(free_): V[:, i] = S[:, j]; V[:, k] = -S[:, j]
            return self.expand(V)
        if sym:
            groups = [(0,)] + [(i, k) for i, k in self.pairs] + [(i,) for i in self.selfneg]
            free_ = []
            for g_ in groups:
                vals = {fixed_orb[i] for i in g_ if i in fixed_orb}
                if len(vals) > 1: return np.zeros((0, self.m), np.int16)
                if vals:
                    for i in g_: fixed_orb[i] = vals.pop() if len(vals) else fixed_orb[g_[0]] if g_[0] in fixed_orb else fixed_orb[i]
                else: free_.append(g_)
            S = allpm(len(free_)); V = np.zeros((len(S), self.n), np.int16)
            for i, val in fixed_orb.items(): V[:, i] = val
            for j, g_ in enumerate(free_):
                for i in g_: V[:, i] = S[:, j]
            return self.expand(V)
        free_ = [i for i in range(self.n) if i not in fixed_orb]
        S = allpm(len(free_)); V = np.zeros((len(S), self.n), np.int16)
        for i, val in fixed_orb.items(): V[:, i] = val
        for j, i in enumerate(free_): V[:, i] = S[:, j]
        return self.expand(V)

    def joint12(self):
        """(u1, u2) with u1(0) = 0, u2(0) = -2 and u1(s) + u1(-s) = -(u2(s) + u2(-s)): 6 options per orbit pair,
        2 per self-negative orbit (where u2 = -u1)."""
        opts = []
        for p1, p2 in itertools.product((1, -1), repeat=2):
            for c2 in ([(-p1, -p1)] if p1 == p2 else [(1, -1), (-1, 1)]): opts.append(((p1, p2), c2))
        O6 = np.array(opts, dtype=np.int8)
        npair, nself = len(self.pairs), len(self.selfneg)
        N = 6 ** npair * 2 ** nself
        idx = np.arange(N, dtype=np.int64)
        V1 = np.zeros((N, self.n), np.int8); V2 = np.zeros((N, self.n), np.int8); V2[:, 0] = -2
        r = idx.copy()
        for (i, k) in self.pairs:
            o = O6[r % 6]; r //= 6
            V1[:, i] = o[:, 0, 0]; V1[:, k] = o[:, 0, 1]; V2[:, i] = o[:, 1, 0]; V2[:, k] = o[:, 1, 1]
        for i in self.selfneg:
            b = 1 - 2 * (r % 2); r //= 2
            V1[:, i] = b; V2[:, i] = -b
        return self.expand(V1), self.expand(V2)


def skew_type(m, x0=1):
    return Orbits(m).skew_type(x0)


def symmetric(m, zero_at_0=False):
    return Orbits(m).symmetric(zero_at_0)


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


def join(lists, C, K, budget=8000000, target=None):
    """Exact 4-list join (join4x) of int16 class vectors to target (default 0); returns (total, [(i,j,k,l), ...])."""
    sz = [len(x) for x in lists]
    if min(sz) == 0: return 0, []
    order = np.argsort(sz); perm = [order[0], order[3], order[1], order[2]]
    path = os.path.join(W.tmpdir(), f'pc_{os.getpid()}.bin')
    T = np.zeros(C, np.int16) if target is None else np.asarray(target, dtype=np.int16)
    with open(path, 'wb') as f:
        f.write(struct.pack('iiiii', *[sz[k] for k in perm], C)); f.write(T.tobytes())
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


def gf_model(p, k):
    """GF(p^k) = F_p[x]/(f) with f irreducible: (multiplication function on coordinate vectors, quadratic character
    on all p^k coordinate vectors in mixed-radix order, the multiplication matrix of a generator of the squares)."""
    import itertools as it
    def mulpoly(a, b, f):
        r = [0] * (2 * k - 1)
        for i, x in enumerate(a):
            for j, y in enumerate(b): r[i + j] = (r[i + j] + x * y) % p
        for d in range(2 * k - 2, k - 1, -1):            # reduce by monic f (f[k] = 1)
            c = r[d]
            if c:
                for i in range(k + 1): r[d - k + i] = (r[d - k + i] - c * f[i]) % p
        return tuple(r[:k])
    for f in it.product(range(p), repeat=k):
        f = list(f) + [1]
        if f[0] == 0: continue
        if all(sum(f[i] * pow(x, i, p) for i in range(k + 1)) % p for x in range(p)) or k > 3:
            elems = list(it.product(range(p), repeat=k))
            nz = [e for e in elems if any(e)]
            # find a primitive element
            for g in nz:
                seen, x = set(), (1,) + (0,) * (k - 1)
                for _ in range(p ** k - 1):
                    seen.add(x); x = mulpoly(x, g, f)
                if len(seen) == p ** k - 1 and x == (1,) + (0,) * (k - 1): break
            else:
                continue
            sq = set(); x = (1,) + (0,) * (k - 1); g2 = mulpoly(g, g, f)
            for _ in range((p ** k - 1) // 2): sq.add(x); x = mulpoly(x, g2, f)
            # coordinates (b_0..b_{k-1}) in mixed-radix order with b_0 slowest
            chi = np.array([0 if not any(e) else (1 if e in sq else -1) for e in elems], np.int64)
            M = np.array([mulpoly(tuple(int(i == j) for i in range(k)), g2, f) for j in range(k)], np.int64).T
            return chi, M
    raise ValueError('no irreducible polynomial found')


def lift_field(m, p, k, U, Wv, Z):
    """Sequences on Z_m x GF(p^k) (mixed radix, Z_m slowest): X(a, b) = u(a)[b = 0] + w(a) + z(a) chi(b)."""
    chi, M = gf_model(p, k)
    Q = p ** k
    a = np.repeat(np.arange(m), Q); b = np.tile(np.arange(Q), m)
    return [(u[a] * (b == 0) + w[a] + z[a] * chi[b]).astype(np.int64) for u, w, z in zip(U, Wv, Z)], M


def check_group(mods, seqs):
    cert = dict(group=list(mods), sequences=[''.join('+' if x > 0 else '-' for x in sq) for sq in seqs],
                rowsums=[int(x.sum()) for x in seqs])
    Hm, _ = verify.build_goethals_seidel_group(cert)
    n = Hm.shape[0]; Hm = Hm.astype(np.float64)
    return np.array_equal(Hm @ Hm.T, n * np.eye(n)) and np.array_equal(Hm + Hm.T, 2 * np.eye(n))


def check(v, seqs):
    tot = verify.paf_sum(seqs)
    if not (tot[0] == 4 * v and not tot[1:].any()): return False
    Hm = verify.goethals_seidel(*map(verify.circulant, seqs)).astype(np.float64)
    return np.array_equal(Hm @ Hm.T, 4 * v * np.eye(4 * v)) and np.array_equal(Hm + Hm.T, 2 * np.eye(4 * v))

def general_search(v, m, q, hgen, Nw, Nm, a, log, field=None):
    """One w-position in each of the blocks 1..N_w (translated to 0, w = +1 after negating the block), none in A.
    For N_w <= 3 this is forced (two w-positions in one block, or a pair in A, leave a q^2 term in the
    autocorrelation of the compressions that only another such pair could cancel). No decoupling is assumed:
      stage 1: y_k = u_k + q w_k (the Z_m-compressions) with sum_k PAF_{y_k} = 0;
      stage 2: z-quadruples with sum_k PAF_{z_k}(s) = R(s) = sum_{k <= N_w} (u_k(s) + u_k(-s)) and zero cross terms."""
    h = (m - 1) // 2; t0 = time.time()
    orb = getattr(a, 'orbits', None) or Orbits(m)
    YA = orb.skew_type()
    lists, roles = [YA], ['A']
    for k in range(1, 4):
        if k <= Nw:
            big = q - 2 if k <= Nm else q              # the first N_- w-blocks have u(0) = -2
            lists.append(orb.free(at0=big)); roles.append(f'w{"-" if k <= Nm else "+"}')
        else:
            X = orb.free()
            lists.append(canonical_translates(X) if (k == 3 and orb.trivial) else X); roles.append('free')
    P = [pafs(Y.astype(np.int64), h) for Y in lists]
    total_u, sols = join(P, h, K=(10 ** 6 if a.all else a.maxu))
    log(f"v={v} m={m} q={q} N_w={Nw} N_-={Nm} (general, roles {roles}): stage 1 lists {[len(x) for x in lists]}, "
        f"y-quadruples {total_u} [{time.time()-t0:.0f}s]")
    if not total_u: return 0
    ZA = orb.symmetric(); Z0 = orb.free(zero_at_0=True); Zf = orb.free()
    zl = [ZA] + [Z0 if k <= Nw else Zf for k in range(1, 4)]
    found = 0
    for qi, sol in enumerate(sols):
        y = [lists[k][sol[k]].astype(np.int64) for k in range(4)]
        Wv = [np.zeros(m, np.int64) for _ in range(4)]
        for k in range(1, Nw + 1): Wv[k][0] = 1
        u = [y[k] - q * Wv[k] for k in range(4)]
        R = sum(np.array([u[k][s_] + u[k][(-s_) % m] for s_ in range(1, h + 1)]) for k in range(1, Nw + 1))
        V = [classvec(Z, u[k], h) for k, Z in enumerate(zl)]
        target = np.concatenate([R, np.zeros(h, np.int64)]) if Nw else np.zeros(2 * h, np.int64)
        tz, sz = join(V, 2 * h, K=a.nsol, target=target)
        for s_ in sz:
            Z = [zl[k][s_[k]].astype(np.int64) for k in range(4)]
            if field:
                p_, k_ = field
                seqs, Mg = lift_field(m, p_, k_, u, Wv, Z)
                ok = check_group([m] + [p_] * k_, seqs)
            else:
                seqs = lift(v, m, q, u, Wv, Z)
                ok = check(v, seqs)
            log(f"  y-quadruple {qi}: {tz} z-solutions; lifted to {'Z_%d x GF(%d^%d)' % (m, *field) if field else 'Z_%d' % v}: skew-Hadamard {ok}")
            if ok:
                found += 1
                rec = dict(n=v, gens=(hgen if isinstance(hgen, list) else [hgen]), pattern=[int(x.sum()) for x in seqs], hadamard=True, skew=True,
                           method=f'paley_config general N_w={Nw}', seqs=[''.join('+' if x > 0 else '-' for x in sq) for sq in seqs])
                if field:
                    rec.update(group=[m] + [field[0]] * field[1], automorphisms=[
                        np.block([[np.eye(1, dtype=np.int64), np.zeros((1, field[1]), np.int64)],
                                  [np.zeros((field[1], 1), np.int64), Mg]]).tolist()])
                with open(W.hits(), 'a') as f: f.write(json.dumps(rec) + '\n')
        if qi < 5 or qi % 50 == 0: log(f"  y-quadruple {qi}: z-solutions {tz} [{time.time()-t0:.0f}s]")
        if found and not a.all: break
    log(f"v={v} N_w={Nw} (general): {min(len(sols), qi + 1)} of {total_u} y-quadruples tried, {found} verified hits [{time.time()-t0:.0f}s]")
    return found

def layouts(m, q, Nw, Nm, orb):
    """w-layouts with paired w-positions for N_w = 4, 5: each layout is a list of 4 dicts {position: (w, u)}
    (block 0 = A). The q^2-sized terms of PAF_y must cancel, so pairs come as {0, d} in two blocks, or as {a, -a} in A
    against {0, 2a} in another block; any remaining w-position sits alone at 0 in a free block. Multipliers normalise
    d (or a) to the smallest K-fixed position."""
    fixedpos = [a_ for a_ in range(1, m) if np.count_nonzero(orb.oid == orb.oid[a_]) == 1]
    if not fixedpos: return []
    d = fixedpos[0]
    out = []
    singles = Nw - 4
    if singles not in (0, 1): return []
    def utypes(n):
        return [t for t in itertools.product((0, 1), repeat=n)]
    # type (ii): blocks 1, 2 with pairs {0, d}; optional single in block 3
    for s1 in (1, -1):
        for s2 in (1, -1):
            npos = 4 + singles
            for ut in utypes(npos):
                if sum(ut) != Nm: continue
                L = [dict(), {0: (1, ut[0]), d: (s1, ut[1])}, {0: (1, ut[2]), d: (s2, ut[3])}, dict()]
                if singles: L[3] = {0: (1, ut[4])}
                out.append(('pairs', L))
    # type (i): A pair {a, -a} (w(a) = +1, w(-a) = -1, equal u-types), block 1 pair {0, 2a}
    a_ = d; two = (2 * a_) % m
    if two in fixedpos or np.count_nonzero(orb.oid == orb.oid[two]) == 1:
        for s1 in (1, -1):
            for ut in utypes(3 + singles):
                nm = 2 * ut[0] + sum(ut[1:])
                if nm != Nm: continue
                L = [{a_: (1, ut[0]), (-a_) % m: (-1, ut[0])}, {0: (1, ut[1]), two: (s1, ut[2])}, dict(), dict()]
                if singles: L[2] = {0: (1, ut[3])}
                out.append(('Apair', L))
    return out


def layout_search(v, m, q, hgen, Nw, Nm, a, log):
    h = (m - 1) // 2; t0 = time.time()
    orb = getattr(a, 'orbits', None) or Orbits(m)
    Ls = layouts(m, q, Nw, Nm, orb)
    log(f"v={v} m={m} q={q} N_w={Nw} N_-={Nm}: {len(Ls)} layouts")
    total_found = 0
    for name, L in Ls:
        lists, Wv, U0 = [], [], []
        for k in range(4):
            spec = L[k]
            yv = {p_: (q - 2) * wv if ut else q * wv for p_, (wv, ut) in spec.items()}
            if k == 0:
                Y = orb.fixed(None, yv, skew=True) if spec else orb.skew_type()
            else:
                Y = orb.fixed(None, yv) if spec else orb.free()
            lists.append(Y)
            wk = np.zeros(m, np.int64)
            for p_, (wv, ut) in spec.items(): wk[p_] = wv
            Wv.append(wk)
        if any(len(Y) == 0 for Y in lists): continue
        P = [pafs(Y.astype(np.int64), h) for Y in lists]
        tot, sols = join(P, h, K=(10 ** 6 if a.all else a.maxu))
        log(f"  layout {name} {[sorted((p_, wv, ut) for p_, (wv, ut) in sp.items()) for sp in L]}: lists {[len(x) for x in lists]}, y-quadruples {tot}")
        for qi, sol in enumerate(sols):
            y = [lists[k][sol[k]].astype(np.int64) for k in range(4)]
            u = [y[k] - q * Wv[k] for k in range(4)]
            R = np.zeros(h, np.int64)
            for k in range(4):
                for s_ in range(1, h + 1):
                    R[s_ - 1] += sum(Wv[k][b] * (u[k][(b + s_) % m] + u[k][(b - s_) % m]) for b in np.flatnonzero(Wv[k]))
                    R[s_ - 1] += q * int(np.dot(Wv[k], np.roll(Wv[k], -s_)))
            zl = []
            for k in range(4):
                zero = list(np.flatnonzero(Wv[k]))
                if k == 0: Z = orb.fixed(None, {}, sym=True, zero_pos=zero) if zero else orb.symmetric()
                else: Z = orb.fixed(None, {}, zero_pos=zero) if zero else orb.free()
                zl.append(Z)
            V = [classvec(Z, u[k], h) for k, Z in enumerate(zl)]
            tz, sz = join(V, 2 * h, K=a.nsol, target=np.concatenate([R, np.zeros(h, np.int64)]))
            for s_ in sz:
                Z = [zl[k][s_[k]].astype(np.int64) for k in range(4)]
                seqs = lift(v, m, q, u, Wv, Z)
                ok = check(v, seqs)
                log(f"    y-quadruple {qi}: {tz} z-solutions; lifted: skew-Hadamard {ok}")
                if ok:
                    total_found += 1
                    rec = dict(n=v, gens=(hgen if isinstance(hgen, list) else [hgen]), pattern=[int(x.sum()) for x in seqs],
                               hadamard=True, skew=True, method=f'paley_config layout N_w={Nw}',
                               seqs=[''.join('+' if x > 0 else '-' for x in sq) for sq in seqs])
                    with open(W.hits(), 'a') as f: f.write(json.dumps(rec) + '\n')
            if qi < 3: log(f"    y-quadruple {qi}: z-solutions {tz} [{time.time()-t0:.0f}s]")
            if total_found and not a.all: break
        if total_found and not a.all: break
    log(f"v={v} N_w={Nw} (layouts): {total_found} verified hits [{time.time()-t0:.0f}s]")
    return total_found


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('v', type=int); ap.add_argument('--maxu', type=int, default=200)
    ap.add_argument('--all', action='store_true'); ap.add_argument('--nsol', type=int, default=2)
    ap.add_argument('--general', action='store_true', help='use the general one-w-per-block search (N_w <= 3)')
    ap.add_argument('--field', help='p,k: q = p^k a prime power = 3 (mod 4); group Z_m x GF(q) with m = v / q')
    ap.add_argument('--kgen', help='generators (mod m) of an extra multiplier group K on Z_m; u, z constant on its orbits')
    ap.add_argument('--korder', type=int, default=1,
                    help='also require u, z constant on the orbits of the subgroup of Z_m^* of this odd order (m prime)')
    a = ap.parse_args()
    v = a.v
    log = lambda s: print(s, flush=True)
    if a.field:
        p_, k_ = map(int, a.field.split(',')); q = p_ ** k_; m = v // q
        assert v == m * q and q % 4 == 3 and m % p_
        for Nw in (1, 2, 3):
            Nm = (q + 1) * Nw - 4 * m
            if Nm % 4 == 0 and 0 <= Nm // 4 <= Nw:
                general_search(v, m, q, 0, Nw, Nm // 4, a, log, field=(p_, k_))
        return
    cases = []
    for m, q, hgen in structures(v, 10 ** 6):
        for Nw in (1, 2, 3, 4, 5):
            Nm = ((q + 1) * Nw - 4 * m)
            if Nm % 4 == 0 and 0 <= Nm // 4 <= Nw: cases.append((m, q, hgen, Nw, Nm // 4))
    if not cases: sys.exit(f"{v}: no Paley-product structure with N_w = 1 or 2 ((q+1)N_w - 4N_- = 4m has no solution)")
    for m, q, hgen, Nw, Nm in cases:
        h = (m - 1) // 2
        K = [1]
        if a.korder > 1:
            if (m - 1) % a.korder or any(m % p == 0 for p in range(2, int(m ** 0.5) + 1)):
                log(f"v={v} m={m}: --korder {a.korder} needs m prime and korder | m - 1"); continue
            g = next(g for g in range(2, m) if all(pow(g, (m - 1) // f, m) != 1 for f in range(2, m) if (m - 1) % f == 0
                                                     and all(f % r for r in range(2, f))))
            k0 = pow(g, (m - 1) // a.korder, m); K = sorted({pow(k0, e, m) for e in range(a.korder)})
        if a.kgen:
            gens_ = [int(x) % m for x in a.kgen.split(',')]; K = {1}
            while True:
                K2 = K | {k * g % m for k in K for g in gens_}
                if K2 == K: break
                K = K2
            K = sorted(K)
            if len(K) % 2 == 0: log(f"v={v} m={m}: K={K} has even order"); continue
        orb = Orbits(m, K)
        a.orbits = orb
        if orb.selfneg: log(f"v={v} m={m} K={K}: -1 fixes an orbit, no skew-type block"); continue
        # the K-symmetric version of the generators (identical to the plain ones when K = {1})
        gen_skew = (lambda x0=1: orb.skew_type(x0)) if not orb.trivial else (lambda x0=1: skew_type(m, x0))
        gen_sym = (lambda zero_at_0=False: orb.symmetric(zero_at_0)) if not orb.trivial else (lambda zero_at_0=False: symmetric(m, zero_at_0))
        gen_free = (lambda zero_at_0=False: orb.free(zero_at_0)) if not orb.trivial else (lambda zero_at_0=False: free(m, zero_at_0))
        canon = (lambda X: X) if not orb.trivial else canonical_translates
        hgen_all = [hgen] + [next(x for x in range(1, v) if x % q == 1 and x % m == k) for k in K if k != 1][:1]
        if Nw >= 4:
            layout_search(v, m, q, hgen_all, Nw, Nm, a, log); continue
        if Nw == 3 or a.general:
            general_search(v, m, q, hgen_all, Nw, Nm, a, log); continue
        t0 = time.time()
        UA = gen_skew()
        if Nw == 1:
            if not any(x * x + y * y == 4 * m - 2 for x in range(1, 2 * m, 2) for y in range(1, 2 * m, 2)):
                log(f"v={v} m={m} q={q} N_w=1: 4m - 2 = {4 * m - 2} is not a sum of two odd squares, so no solution"); continue
            U1 = gen_skew(x0=(0 if Nm == 0 else -2))
            U2 = gen_free(); U3 = canon(gen_free())
            P = [pafs(UA, h), pafs(U1, h), pafs(U2, h), pafs(U3, h)]
            lists_u = [UA, U1, U2, U3]
            total_u, sols_u = join(P, h, K=(10 ** 6 if a.all else a.maxu))
            sols_u = [tuple(lists_u[k][sol[k]] for k in range(4)) for sol in sols_u]
        else:
            # joint (u1, u2): u1(0) = 0, u2(0) = -2, symmetric parts opposite
            # per pair {a, -a}: u1 = (s, s) forces u2 = (-s, -s); u1 = (s, -s) leaves u2 = (t, -t) free: 6 options
            opts = []
            for p1, p2 in itertools.product((1, -1), repeat=2):
                for c2 in ([(-p1, -p1)] if p1 == p2 else [(1, -1), (-1, 1)]): opts.append(((p1, p2), c2))
            O6 = np.array(opts, dtype=np.int8)                                  # (6, 2, 2)
            U3 = canon(gen_free())
            PA, P3 = pafs(UA, h), pafs(U3, h)
            total_u, sols_u, lists_u = 0, [], None
            N12 = 6 ** h if orb.trivial else 1; step = 1 << 21
            for c0 in range(0, N12, step):                                     # the joint (u1, u2) list in chunks
                if orb.trivial:
                    idx = np.arange(c0, min(N12, c0 + step), dtype=np.int64)
                    U1 = np.zeros((len(idx), m), np.int8); U2 = np.zeros((len(idx), m), np.int8); U2[:, 0] = -2
                    r = idx.copy()
                    for p in range(h):
                        o = O6[r % 6]; r //= 6
                        U1[:, p + 1] = o[:, 0, 0]; U1[:, m - 1 - p] = o[:, 0, 1]
                        U2[:, p + 1] = o[:, 1, 0]; U2[:, m - 1 - p] = o[:, 1, 1]
                else:
                    U1, U2 = orb.joint12()
                U12 = np.concatenate([U1, U2], axis=1)
                t_c, s_c = join([PA, pafs(U1, h) + pafs(U2, h), P3, np.zeros((1, h), np.int64)], h,
                                K=(10 ** 6 if a.all else a.maxu))
                total_u += t_c
                for sol in s_c:
                    sols_u.append((UA[sol[0]], U12[sol[1]], U3[sol[2]]))
                if not a.all and len(sols_u) >= a.maxu: break
            lists_u = None
        log(f"v={v} m={m} q={q} K={K if len(K) < 12 else len(K)} N_w={Nw} N_-={Nm}: stage 1 u-quadruples {total_u}"
            f"{'' if a.all or Nw == 1 else ' (stopped at --maxu)'} [{time.time()-t0:.0f}s]")
        if not total_u: continue
        ZA = gen_sym(); Z0 = gen_free(zero_at_0=True); Zf = gen_free()
        found = 0; tried = 0
        for qi, sol in enumerate(sols_u):
            if not a.all and tried >= a.maxu: break
            tried += 1
            if Nw == 1:
                u = [x.astype(np.int64) for x in sol]
                zl = [ZA, Z0, Zf, Zf]
            else:
                uA, u12, u3 = sol
                u = [uA.astype(np.int64), u12[:m].astype(np.int64), u12[m:].astype(np.int64), u3.astype(np.int64)]
                zl = [ZA, Z0, Z0, Zf]
            V = [classvec(Z, u[k], h) for k, Z in enumerate(zl)]
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
                        rec = dict(n=v, gens=hgen_all, pattern=[int(x.sum()) for x in seqs], hadamard=True, skew=True,
                                   method=f'paley_config N_w={Nw}', seqs=[''.join('+' if x > 0 else '-' for x in sq) for sq in seqs])
                        with open(W.hits(), 'a') as f: f.write(json.dumps(rec) + '\n')
                if found and not a.all: break
            if qi < 5 or qi % 50 == 0:
                log(f"  u-quadruple {qi}: z-solutions {tz} [{time.time()-t0:.0f}s]")
        log(f"v={v} N_w={Nw}: {tried} of {total_u} u-quadruples tried, {found} verified hits [{time.time()-t0:.0f}s]")


if __name__ == '__main__':
    main()
