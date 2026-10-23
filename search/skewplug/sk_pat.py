"""exact row-sum patterns and allocations of the skew plug-in (all filters EXACT / sound supersets)"""
import numpy as np, itertools, math
from collections import defaultdict
from gsk import block_classes


def allocations(sh, nmax_block=None):
    """all (nA, nB, nC, nD, mA..) allocations passing the defect identity: (r+1) e_p + (r-3) e_m = 4|K| (H_K trivial: weights = positions); B<=C<=D sorted.
    returns list of (al, ms) with ms = (mA, mB, mC, mD)  (m-type counts per block, A even)"""
    N, r = sh.N, sh.r
    out = []
    Wf = sh.kCD.weights(); WA = sh.A_weights()
    for nA in WA:
        for alF in itertools.combinations_with_replacement(Wf, 3):
            E = nA + sum(alF)
            num = (r + 1) * E - 4 * N
            if num % 4: continue
            em = num // 4
            if not (0 <= em <= E): continue
            out.append(((nA,) + alF, em))
    return out


def m_splits(al, em, wA=2):
    """all (mA, mB, mC, mD) with 0 <= m_X <= n_X, mA even (A pairs have equal types), sum = em, canonical for equal (n_X) among B, C, D (sorted)"""
    res = set()
    nA, nB, nC, nD = al
    for mA in range(0, nA + 1, 2):
        for mB in range(0, nB + 1):
            for mC in range(0, nC + 1):
                mD = em - mA - mB - mC
                if 0 <= mD <= nD:
                    t = [(nB, mB), (nC, mC), (nD, mD)]
                    t.sort()
                    res.add((mA, tuple(x[1] for x in t)) if False else (mA,) + tuple(m for (_, m) in sorted([(nB, mB), (nC, mC), (nD, mD)], key=lambda z: (z[0], z[1]))))
    return sorted(res)


def row_patterns_skew(sh, al, ms, budget=5e6):
    """exact feasible (sigma_X, eps_X) patterns of B, C, D (sign-normalised) of the allocation al = (nA, nB, nC, nD) with m-split ms:
    sum_X (sigma^2 + r tau^2) = 4|G|, sum_X (sigma + r eps)^2 = 4|G|, A: sigma = 1, eps = 0, tau_A odd |tau_A| <= N - nA; tau_X parity/range (N - n_X).
    Returns sorted list of ((sB,eB),(sC,eC),(sD,eD)) with the block-sign gauge normalised (sigma > 0 or sigma == 0 and eps > 0)."""
    r, N = sh.r, sh.N; T2 = sh.T2
    cls = [None] + [sorted(block_classes(sh, al[i], ms[i])) for i in (1, 2, 3)]
    if any(len(c) == 0 for c in cls[1:]): return []
    Wn = [N - al[i] for i in range(4)]
    # tau sets (parity of Wn, range Wn): sums of squares reachable by the A block with given tau: all tau in range with parity
    taus = [list(range(-w, w + 1, 2)) for w in Wn]
    # Rsq condition
    out = set()
    Rsq = sh.Rsq - 1       # R_A = 1 (sigma_A = 1, eps_A = 0)
    # tau-quadruple squares multiset
    tsq = defaultdict(set)
    for ta in taus[0]:
        for tb in taus[1]:
            tsq[ta * ta + tb * tb].add(None)
    sq2 = defaultdict(int)
    def pairs_sq(ts1, ts2):
        d = set()
        for a in ts1:
            for b in ts2: d.add(a * a + b * b)
        return d
    S01 = pairs_sq(taus[0], taus[1]); S23 = pairs_sq(taus[2], taus[3])
    def tau_ok(rem):
        # need 7 * (x + y) = rem with x in S01, y in S23
        if rem < 0 or rem % r: return False
        v = rem // r
        return any((v - x) in S23 for x in S01)
    for cB in cls[1]:
        for cC in cls[2]:
            for cD in cls[3]:
                RB = cB[0] + r * cB[1]; RC = cC[0] + r * cC[1]; RD = cD[0] + r * cD[1]
                if RB * RB + RC * RC + RD * RD != Rsq: continue
                s2 = 1 + cB[0] ** 2 + cC[0] ** 2 + cD[0] ** 2
                if not tau_ok(T2 - s2): continue
                tup = []
                for c in (cB, cC, cD):
                    tup.append(c if (c[0] > 0 or (c[0] == 0 and c[1] > 0)) else (-c[0], -c[1]))
                out.add(tuple(tup))
    # symmetry among equal (n, m) of B, C, D: sort
    res = set()
    for tup in out:
        t = list(tup)
        groups = defaultdict(list)
        for b in range(3): groups[(al[b + 1], ms[b + 1])].append(b)
        for g, idx in groups.items():
            vals = sorted(t[b] for b in idx)
            for b, v in zip(idx, vals): t[b] = v
        res.add(tuple(t))
    return sorted(res)
