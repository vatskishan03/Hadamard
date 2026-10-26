"""Helpers for the NON-SPLIT / general-subgroup validation: independent (pure python / tuple) group arithmetic, planted instances for
arbitrary S_gens, enumeration of all coset-sum patterns consistent on Q, and a brute force that shares no code with the engine's
level / quotient machinery (only the orbit partition is taken from hlib.closure_partition, which the legacy tests validated)."""
import sys, os, itertools, math
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
import hlib

# Teichmueller-type generators: block of GL(m, Z4) of order 2^m - 1 (companion matrices of basic primitive polynomials)
T2 = [[0, 3], [1, 3]]                              # order 3 on Z4^2  (x^2 = 3x + 3)
T3 = [[2, 2, 1], [1, 0, 0], [2, 3, 3]]             # order 7 on Z4^3  (the 839 structure's block)


def blockdiag(Tm, mods_rest_mult):
    """matrix on Z4^m x (cyclic factors): Tm on the first m coordinates, the given multipliers on the others"""
    m = len(Tm); r = m + len(mods_rest_mult)
    M = np.zeros((r, r), dtype=int); M[:m, :m] = Tm
    for i, u in enumerate(mods_rest_mult): M[m + i, m + i] = u
    return M.tolist()


# ---------------------------------------------------------------- independent group arithmetic (tuples)
def elems(mods):
    return list(itertools.product(*[range(m) for m in mods]))


def add(a, b, mods):
    return tuple((x + y) % m for x, y, m in zip(a, b, mods))


def indep_cosets(mods, gens):
    """labels of the cosets of <gens> in order of the smallest element (mixed radix); pure python"""
    E = elems(mods); idx = {e: i for i, e in enumerate(E)}
    K = {tuple([0] * len(mods))}; fr = list(K)
    while fr:
        x = fr.pop()
        for g in gens:
            y = add(x, tuple(g), mods)
            if y not in K: K.add(y); fr.append(y)
    lab = [-1] * len(E); q = 0
    for i, e in enumerate(E):
        if lab[i] < 0:
            for k in K: lab[idx[add(e, k, mods)]] = q
            q += 1
    return np.array(lab), q, len(K)


def planted(mods, S_gens, gens, seed, same_y=()):
    """random H-invariant sequences; explicit target = sum of their PAF (exact, computed by direct summation)"""
    G = hlib.Group(mods); perms = [G.perm_from_matrix(A) for A in gens]
    orb, O = hlib.closure_partition(G.v, perms)
    rng = np.random.default_rng(seed)
    cos, k, _ = indep_cosets(mods, S_gens)
    seqs = [(2 * rng.integers(0, 2, O) - 1)[orb] for _ in range(4)]
    for a, b in same_y:
        for _ in range(100000):
            x = (2 * rng.integers(0, 2, O) - 1)[orb]
            if np.array_equal(np.bincount(cos, weights=x, minlength=k), np.bincount(cos, weights=seqs[a], minlength=k)): seqs[b] = x; break
    ys = [np.bincount(cos, weights=x, minlength=k).round().astype(int).tolist() for x in seqs]
    tgt = np.zeros(G.v, dtype=np.int64)
    for s in range(G.v):
        sh = G.translate(s); tgt[s] = sum(int(x @ x[sh]) for x in seqs)
    return dict(mods=mods, S_gens=S_gens, gens=gens, ys=ys, tgt=tgt.tolist(), name='planted%d' % seed), seqs


def real_cfg(mods, S_gens, gens, ys, t_in, t_out):
    return dict(mods=mods, S_gens=S_gens, gens=gens, ys=ys, t_in=t_in, t_out=t_out, name='real')


def all_patterns(mods, S_gens, gens, t_in, t_out, sorted_only=True):
    """every quadruple of coset-sum vectors (constant on the H-orbits of Q) with sum_i PAF_{y_i} = target on Q (the level-0 equations).
    Returned as lists of 4 lists (one value per coset).  Uses hlib's Q level only to evaluate the Q-level PAF / target; the
    Q-level correctness is validated separately against the independent quotient computation."""
    cfg0 = dict(mods=mods, S_gens=S_gens, gens=gens, ys=[[0] * 1] * 4, t_in=t_in, t_out=t_out)
    G = hlib.Group(mods); cos, k, ssz = indep_cosets(mods, S_gens)
    cfg0['ys'] = [[0] * k] * 4
    P = hlib.Problem(cfg0)
    Q = P.Qlevel; cs = P.v // k
    vals = list(range(-cs, cs + 1, 2))
    feats = {}
    for assign in itertools.product(vals, repeat=Q.O):
        yo = np.array(assign, dtype=np.int64)
        feats[assign] = Q.paf(yo[None, :])[0].astype(np.int64)
    keys = list(feats)
    pair = {}
    for a in keys:
        for b in keys:
            pair.setdefault(tuple(feats[a] + feats[b]), []).append((a, b))
    out = []
    T = Q.T
    for key, lst in pair.items():
        need = tuple(T - np.array(key))
        if need in pair:
            for (a, b) in lst:
                for (c, d) in pair[need]:
                    t = (a, b, c, d)
                    if sorted_only and list(t) != sorted(t): continue
                    out.append([[int(t_i[Q.orb[g]]) for g in range(Q.GK.v)] for t_i in t])
    return out


def brute_all(cfg, maxn=8000, maxmasks=1 << 24):
    """all ordered quadruples of orbit-sign masks (bit o = +1) with the given coset sums whose PAF sum equals tgt.  Independent of the
    engine's levels: coset sums by direct summation over the independently computed cosets, PAF by direct translation sums."""
    mods = cfg['mods']; G = hlib.Group(mods)
    perms = [G.perm_from_matrix(A) for A in cfg['gens']]
    orb, O = hlib.closure_partition(G.v, perms)
    cos, k, ssz = indep_cosets(mods, cfg['S_gens'])
    v = G.v
    if 'tgt' in cfg: tgt = np.array(cfg['tgt'], dtype=np.int64)
    else:
        tgt = np.where(cos == 0, cfg['t_in'], cfg['t_out']).astype(np.int64); tgt[0] = 4 * v
    assert (1 << O) <= maxmasks, O
    # coset-count matrix cnt[o, c]
    cnt = np.zeros((O, k), dtype=np.int64); np.add.at(cnt, (orb, cos), 1)
    lists = []
    ysets = {tuple(y) for y in cfg['ys']}
    cache = {}
    for i in range(4):
        key = tuple(cfg['ys'][i])
        if key in cache: lists.append(cache[key]); continue
        keep = []
        B = 1 << 18
        for s in range(0, 1 << O, B):
            ms = np.arange(s, min(s + B, 1 << O), dtype=np.uint64)
            Y = 2 * ((ms[:, None] >> np.arange(O, dtype=np.uint64)) & np.uint64(1)).astype(np.int64) - 1
            cs = Y @ cnt
            ok = np.all(cs == np.array(key)[None, :], axis=1)
            keep.append(ms[ok])
        lists.append(np.concatenate(keep)); cache[key] = lists[-1]
    n = [len(l) for l in lists]
    assert max(n) <= maxn, n
    # class representatives of <H,-1> acting on G (independent: orbits of the permutation group generated by perms and negation)
    neg = np.array([[(-c) % m for c, m in zip(e, mods)] for e in G.coords.tolist()])
    negp = G.index(neg)
    cl, _ = hlib.closure_partition(v, perms + [negp])
    reps = [int(np.flatnonzero(cl == c)[0]) for c in sorted(set(cl.tolist()) - {cl[0]})]
    Vs = []
    for l in lists:
        Y = 2 * ((l[:, None] >> np.arange(O, dtype=np.uint64)) & np.uint64(1)).astype(np.int64) - 1
        X = Y[:, orb]
        V = np.zeros((len(l), len(reps)), dtype=np.int64)
        for c, d in enumerate(reps):
            # translate by d, independently: element index of (e + d)
            sh = np.array([i for i in range(v)]); sh = G.index(G.coords + G.coords[d][None, :])
            V[:, c] = (X * X[:, sh]).sum(1)
        Vs.append(V)
    # also verify the whole PAF identity on hits over ALL d (not only the representatives) below
    tgtc = np.array([tgt[r] for r in reps])
    rng = np.random.default_rng(5); key = rng.integers(1, 2 ** 40, size=len(reps))
    K = [V @ key for V in Vs]; tk = int(tgtc @ key)
    R = np.add.outer(K[2], K[3]).ravel()
    order = np.argsort(R, kind='stable'); Rs = R[order]; del R
    order = order.astype(np.int32)
    sols = []
    CH = max(1, int(4e6 // max(n[1], 1)))
    for a0 in range(0, n[0], CH):
        Lc = np.add.outer(K[0][a0:a0 + CH], K[1]).ravel()
        need = tk - Lc
        lo = np.searchsorted(Rs, need, 'left'); hi = np.searchsorted(Rs, need, 'right')
        for h in np.flatnonzero(hi > lo):
            a, b = divmod(int(h), n[1]); a += a0
            for p in range(lo[h], hi[h]):
                c, d = divmod(int(order[p]), n[3])
                if np.array_equal(Vs[0][a] + Vs[1][b] + Vs[2][c] + Vs[3][d], tgtc):
                    sols.append((int(lists[0][a]), int(lists[1][b]), int(lists[2][c]), int(lists[3][d])))
    # exact re-verification over all of G (not only the class reps)
    for ms in sols[:200]:
        tot = np.zeros(v, dtype=np.int64)
        for m in ms:
            x = (2 * ((m >> np.arange(O)) & 1) - 1)[orb]
            for s in range(v): tot[s] += int(x @ x[G.index(G.coords + G.coords[s][None, :])])
        assert np.array_equal(tot, tgt), "brute-force class reduction is inconsistent"
    return sols, n


def planted_sym(mods, S_gens, gens, nu, seed):
    """planted instance with blocks (0,1) and (2,3) exchanged by the automorphism nu (x1 = x0 o nu, x3 = x2 o nu)"""
    G = hlib.Group(mods); perms = [G.perm_from_matrix(A) for A in gens]
    orb, O = hlib.closure_partition(G.v, perms)
    nup = G.perm_from_matrix(nu)
    rng = np.random.default_rng(seed)
    cos, k, _ = indep_cosets(mods, S_gens)
    x0 = (2 * rng.integers(0, 2, O) - 1)[orb]; x2 = (2 * rng.integers(0, 2, O) - 1)[orb]
    seqs = [x0, x0[nup], x2, x2[nup]]
    ys = [np.bincount(cos, weights=x, minlength=k).round().astype(int).tolist() for x in seqs]
    tgt = np.zeros(G.v, dtype=np.int64)
    for s in range(G.v):
        sh = G.translate(s); tgt[s] = sum(int(x @ x[sh]) for x in seqs)
    return dict(mods=mods, S_gens=S_gens, gens=gens, ys=ys, tgt=tgt.tolist(), name='plantedsym%d' % seed), seqs
