r"""Exact fibre test for spin structures over a composite Z_v (spin_skew.py).

For a shift s in the subgroup dZ_v (d | v, 1 < d < v) the periodic autocorrelation splits over the d cosets
a + dZ_v ("fibres"):  P_X(s) = sum_a P_{X_a}(s), X_a(y) = X(a + d y).  Hence on the <K,-1>-class of such an s

    P_A(s) = sum_a P_{A_a}(s),        T_B(s) = sum_j P_B(mu^j s) = sum_a sum_j P_{B_a}(mu^j s),

both sums of independent fibre terms.  H permutes the fibres (the fibre g a contributes the term of fibre a at the
pulled-back shift g^-1 s) and the stabiliser of a fibre acts on it, so the attainable values follow from a DP over fibre orbits, exactly; the row sum of
B (b) and the skew type of A are carried along.  The spin condition needs P_A(s) + T_B(s) = 0, so if the attainable
sets of P_A(s) and -T_B(s) are disjoint for one such class, the structure has no solution.

  python3 search/spin_fibre.py v hgens mu      (hgens comma-separated)
"""
import itertools
import math
import sys
from math import gcd

import numpy as np


def closure(gens, v):
    S = {1}
    while True:
        T = S | {s * g % v for s in S for g in gens}
        if T == S:
            return sorted(S)
        S = T


def fibre_values(v, d, G, shifts, skew):
    """attainable (sum of fibre PAF terms at `shifts`, row sum) over G-invariant +-1 sequences on Z_v
    (skew: x(0) = 1, x(-t) = -x(t)); returns dict rowsum -> set of values"""
    # orbits of G on Z_v
    orb = -np.ones(v, dtype=np.int64)
    obs = []
    for x in range(v):
        if orb[x] < 0:
            o = sorted({x * g % v for g in G})
            for y in o:
                orb[y] = len(obs)
            obs.append(o)
    # fibres = cosets of dZ_v, indexed by a in 0..d-1; G permutes them: a -> g a mod d
    fib_orbits, seen = [], set()
    for a in range(d):
        if a in seen:
            continue
        o = sorted({a * g % d for g in G})
        seen |= set(o)
        fib_orbits.append(o)
    total = {(0, 0)}
    done = set()
    for fo in fib_orbits:
        a = fo[0]
        if a in done:
            continue
        done |= set(fo)
        negf = sorted({(-a) * g % d for g in G})
        paired = skew and negf != fo           # A on the fibre orbit -F is determined by A on F (A(-x) = -A(x))
        if paired:
            done |= set(negf)
        pts = [a + d * y for y in range(v // d)]
        vars_ = sorted({int(orb[x]) for x in pts})        # G-orbits meeting this fibre
        tied = {}
        free = []
        if skew and not paired:
            neg = {o: int(orb[(-obs[o][0]) % v]) for o in vars_}
            for o in vars_:
                if obs[o][0] == 0 or o in tied:
                    continue
                free.append(o)
                tied[o] = (o, 1)
                tied[neg[o]] = (o, -1)
        else:
            free = list(vars_)
        # fibre a' = g a contributes P_{X_a}(g^-1 s) at a shift s (X(g x) = X(x)), so collect the pulled-back shifts
        trans = {}
        for g in G:
            trans.setdefault(a * g % d, g)
        pulled = [(sh * pow(trans[a2], -1, v)) % v for a2 in fo for sh in shifts]
        out = set()
        mult = 2 if paired else 1
        for signs in itertools.product((1, -1), repeat=len(free)):
            sv = dict(zip(free, signs))
            x = np.zeros(v // d, dtype=np.int64)
            for y in range(v // d):
                o = int(orb[pts[y]])
                if obs[o][0] == 0:
                    x[y] = 1
                elif skew and not paired:
                    r, sg = tied[o]
                    x[y] = sg * sv[r]
                else:
                    x[y] = sv[o]
            val = 0
            for sh in pulled:
                k = (sh // d) % (v // d)
                val += int(np.dot(x, np.roll(x, -k)))
            rs = 0 if paired else len(fo) * int(x.sum())
            out.add((mult * val, rs))
        total = {(t0 + a0, t1 + a1) for (t0, t1) in total for (a0, a1) in out}
    res = {}
    for val, rs in total:
        res.setdefault(rs, set()).add(val)
    return res


def main():
    v = int(sys.argv[1])
    hg = [int(x) for x in sys.argv[2].split(',')]
    mu = int(sys.argv[3])
    H = closure(hg, v)
    K = closure(hg + [mu], v)
    assert v - 1 not in K
    b = math.isqrt((4 * v - 1) // 3)
    assert 1 + 3 * b * b == 4 * v
    # <K,-1>-class representatives that lie in proper subgroups
    seen, reps = set(), []
    for x in range(1, v):
        if x in seen:
            continue
        cl = {x * k % v for k in K} | {(-x * k) % v for k in K}
        seen |= cl
        if gcd(x, v) > 1:
            reps.append(x)
    dead = False
    for s in reps:
        d = gcd(s, v)
        PA = fibre_values(v, d, K, [s], True)
        TB = fibre_values(v, d, H, [s, s * mu % v, s * mu * mu % v], False)
        pa = set().union(*PA.values())
        tb = TB.get(b, set()) | TB.get(-b, set())
        ok = sorted(x for x in pa if -x in tb)
        print(f'class of {s} (gcd {d}): P_A in {sorted(pa)[:12]}{"..." if len(pa) > 12 else ""}; '
              f'T_B (row sum +-{b}) from {min(tb) if tb else None}; matching P_A values: {ok}')
        dead |= not ok
    print('RESULT', v, hg, mu, 'infeasible' if dead else 'passes')


if __name__ == '__main__':
    main()
