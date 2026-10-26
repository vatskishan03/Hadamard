"""Independent exact check of a solution record (four sequences on G = Z2^3 x Z29) of the n = 955 core.

usage: python3 check_core.py SOLUTION.json [...]
Uses nothing of zs.py / zq.py: integer arithmetic on the group table only. Checks (all must hold):
  * A, B, C, D have length 232 (index ((e0*2+e1)*2+e2)*29+z), entries +-1, except A(0) = B(0) = 0;
  * A and B symmetric (x(-g) = x(g));
  * invariance under h = (tau, x16): (e, z) -> (tau e, 16 z), tau = [[0,0,1],[1,0,1],[0,1,0]] over GF(2) (order 7), 16 of order 7 mod 29;
  * row sums (21, 21, 6, 6) and coset sums over Q = Z2^3 (y(0), y(a)) = (0, 3), (0, 3), (-1, 1), (-1, 1);
  * the periodic autocorrelations PAF_A + PAF_B + PAF_C + PAF_D = 926 at 0, -28 on K0 minus {0} (K0 = 0 x Z29), +4 off K0, by the
    group table (no FFT); equivalently (Fourier, checked too) sum_k |X_k^(psi)|^2 = 954 = 4N - 2 for the characters psi of G that
    are not lifts of a nontrivial character of Q (the trivial one included: sum of the squared row sums) and 26 = 4t - 2 for the
    7 lifts (N = 232 + 7 = 239, t = 7; the zeros of A and B at 0 lower both constants by 2).
"""
import json, sys, itertools
import numpy as np

mods = (2, 2, 2, 29); m = 232
coords = np.array(list(itertools.product(*[range(x) for x in mods])), dtype=np.int64)
mo = np.array(mods)
def idx(c): return np.ravel_multi_index(tuple((c % mo).T), mods)
add = idx(coords[:, None, :] + coords[None, :, :])          # add[g, h] = index of g + h
neg = idx(-coords)
tau = np.array([[0,0,1],[1,0,1],[0,1,0]])
Hm = np.zeros((4, 4), dtype=np.int64); Hm[:3, :3] = tau; Hm[3, 3] = 16
hperm = idx(coords @ Hm.T)
assert len(set(hperm.tolist())) == m
# order of h
q = np.arange(m); k = 0
while True:
    q = hperm[q]; k += 1
    if (q == np.arange(m)).all(): break
assert k == 7, k

def check(path):
    sol = json.load(open(path))
    X = [np.array(sol[c], dtype=np.int64) for c in "ABCD"]
    ok = True
    def rep(name, cond):
        nonlocal ok
        print("  %-62s %s" % (name, "OK" if cond else "FAIL")); ok &= bool(cond)
    rep("lengths 232, entries +-1 (A, B zero exactly at 0)",
        all(len(x) == m for x in X) and all((np.abs(x[1:]) == 1).all() for x in X[:2]) and X[0][0] == 0 and X[1][0] == 0
        and all((np.abs(x) == 1).all() for x in X[2:]))
    rep("A, B symmetric", all(np.array_equal(x, x[neg]) for x in X[:2]))
    rep("invariant under h = (tau, x16), order 7", all(np.array_equal(x, x[hperm]) for x in X))
    rep("row sums (21, 21, 6, 6)", [int(x.sum()) for x in X] == [21, 21, 6, 6])
    e = coords[:, 0] * 4 + coords[:, 1] * 2 + coords[:, 2]
    ys = [[int(x[e == a].sum()) for a in range(8)] for x in X]
    want = [[0] + [3] * 7, [0] + [3] * 7, [-1] + [1] * 7, [-1] + [1] * 7]
    rep("coset sums over Q = Z2^3: (0,3), (0,3), (-1,1), (-1,1)", ys == want)
    # PAF by the group table: PAF(h) = sum_g x(g) x(g + h)
    paf = np.zeros(m, dtype=np.int64)
    for x in X:
        paf += np.array([int((x * x[add[:, h]]).sum()) for h in range(m)])
    inK = (coords[:, :3] == 0).all(1)
    tgt = np.where(inK, -28, 4); tgt[0] = 926
    rep("sum PAF = 926 at 0, -28 on K0 \\ {0}, +4 off K0", np.array_equal(paf, tgt))
    # Fourier form
    P = np.zeros(mods)
    for x in X: P += np.abs(np.fft.fftn(x.reshape(mods).astype(float))) ** 2
    P = P.reshape(-1)
    lift = np.zeros(m, dtype=bool)                          # characters psi_(a, b), a in Z2^3 (a != 0) with b = 0
    for i, c in enumerate(coords):
        lift[i] = (c[3] == 0) and c[:3].any()
    rep("sum |X^(psi)|^2 = 26 on the 7 lifts, 954 elsewhere", np.allclose(P[lift], 26) and np.allclose(P[~lift], 954) and lift.sum() == 7)
    return ok

if not sys.argv[1:]:
    print(__doc__); sys.exit(2)
allok = True
for p in sys.argv[1:]:
    print(p)
    r = check(p); allok &= r
    print("  -> ALL CHECKS PASSED" if r else "  -> FAILED")
sys.exit(0 if allok else 1)
