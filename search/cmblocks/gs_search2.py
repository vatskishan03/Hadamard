"""Faster exhaustive search of the linear GS-type family (see gs_search.py).  Rows 2 are reduced modulo the symmetry
group fixing row 1 = (A, B, C, D): permutations of (columns, variables) together, phase scalings of a column with the
inverse scaling of its variable, the unit scaling of row 2, and global complex conjugation (phases e -> -e);
rows 3 and 4 are then found for every canonical row 2."""
import sys, itertools, time
import numpy as np
rng = np.random.default_rng(12345)
mode = sys.argv[1] if len(sys.argv) > 1 else "complex"
u = rng.normal(size=4) + 1j * rng.normal(size=4)
v = rng.normal(size=4) + 1j * rng.normal(size=4)
if mode == "real": v = u.conj()
S = np.array([[0, 1], [1, 0]], complex)
OPTS, MATS = [], []
for var in range(4):
    for tau in range(4):
        x, xp = {0: (u[var], v[var]), 1: (v[var], u[var]), 2: (v[var].conjugate(), u[var].conjugate()), 3: (u[var].conjugate(), v[var].conjugate())}[tau]
        for rho in range(2):
            for e in range(4):
                OPTS.append((var, tau, rho, e)); MATS.append((1j ** e) * np.diag([x, xp]) @ (S if rho else np.eye(2)))
MATS = np.array(MATS); NOPT = len(OPTS)
IDX = {o: i for i, o in enumerate(OPTS)}
VAR = np.array([o[0] for o in OPTS])
W = rng.integers(1, 2 ** 62, size=64, dtype=np.int64) | 1
def hashes(arr):            # arr (..., K) complex -> int64 hash of rounded values
    r = np.concatenate([arr.real, arr.imag], axis=-1)
    q = np.rint(r * 1e5).astype(np.int64)
    return (q * W[:q.shape[-1]]).sum(axis=-1)
def tabs_for(prev_rows):
    tabs = []
    for j in range(4):
        cols = [MATS[r[j]] @ np.conj(np.swapaxes(MATS, -1, -2)) for r in prev_rows]    # each (NOPT,2,2)
        tabs.append(np.concatenate([c.reshape(NOPT, 4) for c in cols], axis=1))         # NOPT x (4 nprev)
    return tabs
vv = (VAR[:, None] != VAR[None, :])
def rows_orthogonal_to(prev_rows):
    t0, t1, t2, t3 = tabs_for(prev_rows)
    L = t0[:, None, :] + t1[None, :, :]                  # 128 x 128 x K
    Rr = -(t2[:, None, :] + t3[None, :, :])
    hL, hR = hashes(L).ravel(), hashes(Rr).ravel()
    okL, okR = vv.ravel(), vv.ravel()
    iL = np.flatnonzero(okL); iR = np.flatnonzero(okR)
    order = np.argsort(hL[iL]); sl = hL[iL][order]; il = iL[order]
    lo = np.searchsorted(sl, hR[iR], 'left'); hi = np.searchsorted(sl, hR[iR], 'right')
    cnt = hi - lo
    out = []
    for k in np.flatnonzero(cnt):
        c_, d_ = divmod(int(iR[k]), NOPT)
        for m_ in range(lo[k], hi[k]):
            a_, b_ = divmod(int(il[m_]), NOPT)
            if len({VAR[a_], VAR[b_], VAR[c_], VAR[d_]}) == 4:
                out.append((a_, b_, c_, d_))
    return out
def opt(var, tau=0, rho=0, e=0): return IDX[(var, tau, rho, e)]
row1 = tuple(opt(j) for j in range(4))
t0 = time.time()
R2 = rows_orthogonal_to([row1])
print(mode, "rows orthogonal to row 1:", len(R2), f"({time.time()-t0:.1f}s)", flush=True)
# ---- symmetry reduction of row 2 ----
def relabel(o, perm):                  # variable permutation
    var, tau, rho, e = OPTS[o]; return IDX[(perm[var], tau, rho, e)]
def gen_colperm(p):                    # columns j -> p[j], variables j -> p[j]
    def f(row):
        new = [None] * 4
        for j in range(4): new[p[j]] = relabel(row[j], p)
        return tuple(new)
    return f
def gen_phase(j, d):                   # column j scaled by i^d, variable j rescaled by i^-d
    def f(row):
        new = list(row)
        for jj in range(4):
            var, tau, rho, e = OPTS[row[jj]]
            e2 = e
            if jj == j: e2 += d
            if var == j: e2 += (-d if tau in (0, 1) else d)
            new[jj] = IDX[(var, tau, rho, e2 % 4)]
        return tuple(new)
    return f
def rowscale(d):
    def f(row):
        return tuple(IDX[(OPTS[o][0], OPTS[o][1], OPTS[o][2], (OPTS[o][3] + d) % 4)] for o in row)
    return f
def conjall(row):
    return tuple(IDX[(OPTS[o][0], OPTS[o][1], OPTS[o][2], (-OPTS[o][3]) % 4)] for o in row)
gens = [gen_colperm([1, 0, 2, 3]), gen_colperm([1, 2, 3, 0])] + [gen_phase(j, 1) for j in range(4)] + [rowscale(1), conjall]
# NB: conjall is a symmetry only if the instance (u, v) is closed under it; the formal array statement is instance independent
# (random instance only tests identities), so the symmetry holds on the formal level.  Verify closure of R2 as a safeguard.
NOREDUCE = len(sys.argv) > 2 and sys.argv[2] == "noreduce"
R2s = set(R2)
for g_ in ([] if NOREDUCE else gens):
    for r in R2:
        assert g_(r) in R2s, "R2 not closed under a symmetry generator"
parent = {r: r for r in R2}
def find(x):
    while parent[x] != x:
        parent[x] = parent[parent[x]]; x = parent[x]
    return x
for r in R2:
    for g_ in ([] if NOREDUCE else gens):
        a, b = find(r), find(g_(r))
        if a != b: parent[a] = b
reps = sorted({find(r) for r in R2})
print("canonical row-2 classes:", len(reps), flush=True)
total = 0; found = []
for n2, r2 in enumerate(reps):
    R3 = rows_orthogonal_to([row1, r2])
    # rows 3 and 4: row 3 orthogonal to rows 1, 2; row 4 orthogonal to 1, 2, 3
    for r3 in R3:
        R4 = rows_orthogonal_to([row1, r2, r3])
        for r4 in R4:
            if all(len({VAR[row[j]] for row in (row1, r2, r3, r4)}) == 4 for j in range(4)):
                total += 1
                if len(found) < 4: found.append((row1, r2, r3, r4))
    print(f"  class {n2+1}/{len(reps)}: |R3| = {len(R3)}, arrays so far {total}  ({time.time()-t0:.0f}s)", flush=True)
print(mode, "arrays with orthogonal rows and Latin columns (row 2 canonical):", total)
for f in found[:3]:
    for row in f: print("   ", [OPTS[o] for o in row])
    print()
