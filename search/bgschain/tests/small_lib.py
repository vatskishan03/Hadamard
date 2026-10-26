"""Small planted instances + independent brute force for validating the engine."""
import sys, os, itertools, json, copy
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
import hlib

AL = [[0, 1], [1, 1]]; I2 = [[1, 0], [0, 1]]
def mat(A, a, b, extra=None):
    M = np.zeros((4, 4), dtype=int); M[:2, :2] = A; M[2, 2] = a; M[3, 3] = b; return M.tolist()


def planted(mods, nq, gens, seed, same_y=(), ntries=1):
    """random H-invariant sequences; explicit target = sum of their PAF; returns cfg (without levels)"""
    cfg0 = dict(mods=mods, nq=nq, gens=gens, ys=[[0] * 4] * 4, t_in=0, t_out=0)
    P0 = hlib.Problem(dict(cfg0, tgt=None) if False else {**cfg0, 'tgt': [0] * int(np.prod(mods))}) if False else None
    G = hlib.Group(mods); perms = [G.perm_from_matrix(A) for A in gens]
    orb, O = hlib.closure_partition(G.v, perms)
    rng = np.random.default_rng(seed)
    cos = np.ravel_multi_index(tuple(G.coords[:, :nq].T), tuple(mods[:nq]))
    k = int(np.prod(mods[:nq]))
    seqs = []
    for i in range(4):
        if i in same_y and i > 0 and False: pass
        b = rng.integers(0, 2, O); x = (2 * b - 1)[orb]; seqs.append(x)
    for a, b in same_y:        # force equal coset sums by drawing until equal
        for _ in range(100000):
            x = (2 * rng.integers(0, 2, O) - 1)[orb]
            if np.array_equal(np.bincount(cos, weights=x, minlength=k), np.bincount(cos, weights=seqs[a], minlength=k)): seqs[b] = x; break
    ys = [np.bincount(cos, weights=x, minlength=k).round().astype(int).tolist() for x in seqs]
    v = G.v; tgt = np.zeros(v, dtype=np.int64)
    for s in range(v):
        sh = G.translate(s)
        tgt[s] = sum(int(x @ x[sh]) for x in seqs)
    return dict(mods=mods, nq=nq, gens=gens, ys=ys, tgt=tgt.tolist(), name='planted%d' % seed), seqs


def brute(cfg, maxn=6000):
    """all ordered quadruples of orbit-sign masks with the given coset sums whose PAF sum equals tgt (independent code path)"""
    P = hlib.Problem({**cfg, 'levels': []})
    lists = []
    for i in range(4):
        ms = np.arange(1 << P.O, dtype=np.uint64)
        Y = 2 * ((ms[:, None] >> np.arange(P.O, dtype=np.uint64)) & np.uint64(1)).astype(np.int64) - 1
        X = Y[:, P.orb]
        cs = np.zeros((len(ms), P.k), dtype=np.int64)
        for c in range(P.k): cs[:, c] = X[:, P.cos == c].sum(1)
        keep = np.all(cs == np.array(P.ys[i])[None, :], axis=1)
        lists.append(ms[keep])
    n = [len(l) for l in lists]
    assert max(n) <= maxn, n
    # PAF class values by direct summation over G (independent of the S tensors)
    Vs = []
    for l in lists:
        Y = 2 * ((l[:, None] >> np.arange(P.O, dtype=np.uint64)) & np.uint64(1)).astype(np.int64) - 1
        X = Y[:, P.orb]
        V = np.zeros((len(l), P.C), dtype=np.int64)
        for c, d in enumerate(P.reps):
            sh = P.G.translate(d)
            V[:, c] = (X * X[:, sh]).sum(1)
        Vs.append(V)
    rng = np.random.default_rng(5); key = rng.integers(1, 2 ** 40, size=P.C)
    K = [V @ key for V in Vs]
    tk = int(P.tgtc @ key)
    L = np.add.outer(K[0], K[1]).ravel(); R = np.add.outer(K[2], K[3]).ravel()
    order = np.argsort(R, kind='stable'); Rs = R[order]
    need = tk - L
    lo = np.searchsorted(Rs, need, 'left'); hi = np.searchsorted(Rs, need, 'right')
    sols = []
    for h in np.flatnonzero(hi > lo):
        a, b = divmod(int(h), n[1])
        for p in range(lo[h], hi[h]):
            c, d = divmod(int(order[p]), n[3])
            if np.array_equal(Vs[0][a] + Vs[1][b] + Vs[2][c] + Vs[3][d], P.tgtc):
                sols.append((int(lists[0][a]), int(lists[1][b]), int(lists[2][c]), int(lists[3][d])))
    return sols, n


def planted_sym(mods, nq, gens, nu, seed):
    """planted instance with blocks (0,1) and (2,3) exchanged by the normalizer element nu (x_1 = x_0 o nu, x_3 = x_2 o nu): the
    target is nu-invariant, so the solution set is invariant under nu"""
    G = hlib.Group(mods); perms = [G.perm_from_matrix(A) for A in gens]
    orb, O = hlib.closure_partition(G.v, perms)
    nup = G.perm_from_matrix(nu)
    rng = np.random.default_rng(seed)
    cos = np.ravel_multi_index(tuple(G.coords[:, :nq].T), tuple(mods[:nq])); k = int(np.prod(mods[:nq]))
    x0 = (2 * rng.integers(0, 2, O) - 1)[orb]; x2 = (2 * rng.integers(0, 2, O) - 1)[orb]
    seqs = [x0, x0[nup], x2, x2[nup]]
    ys = [np.bincount(cos, weights=x, minlength=k).round().astype(int).tolist() for x in seqs]
    tgt = np.zeros(G.v, dtype=np.int64)
    for s in range(G.v):
        sh = G.translate(s); tgt[s] = sum(int(x @ x[sh]) for x in seqs)
    return dict(mods=mods, nq=nq, gens=gens, ys=ys, tgt=tgt.tolist(), name='plantedsym%d' % seed), seqs
