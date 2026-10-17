r"""Two-block skew-Hadamard matrices of order 2(w + 1) with a general second block (multiplier-invariant search).

Let G be an abelian group of odd order w, a a +-1 function of skew type (a(0) = 1, a(-x) = -a(x)) and b ANY +-1
function with row sum 1, such that PAF_a(s) + PAF_b(s) = -2 for every s != 0.  With A, B group-developed
(A + A^T = 2I, AA^T + BB^T = (2w + 2)I - 2J, row sums 1),

        [  1    1    e    e  ]
    H = [ -1    1    e   -e  ]
        [ -e^T -e^T  A    B  ]
        [ -e^T  e^T -B^T  A^T]

is a skew-Hadamard matrix of order 2(w + 1): H + H^T = 2I needs only A + A^T = 2I (the off-diagonal blocks pair B
with -B^T), and every row product reduces to the PAF condition and the row sums.  search/szekeres.py (NOTES §15)
required b symmetric (the classical Szekeres difference sets); here b is only H-invariant, which doubles the free
signs of b, so its 767 structures without a symmetric b and the others are open again.

Search: a ranges over the H-invariant skew functions (one sign per pair of opposite H-orbits), b over the H-invariant
functions with row sum 1; lists of PAF vectors at the <H,-1>-orbit representatives (gdfs.c), joined by join4x with
target (-2, ..., -2). Every hit is rebuilt as a matrix and checked.

With --sym: symmetric Hadamard matrices of order 2(w + 1) from a Legendre pair with A symmetric (NOTES §29): a and
b with row sums -1, a symmetric, PAF_a + PAF_b = -2, in [[X, E], [E^T, M]] with X = [[1, 1], [1, -1]],
E = [[e, e], [e, -e]], M = [[A, B], [B^T, -A]] (symmetric since A is).

Usage: python3 search/szek2.py GROUP "HGENS" [--sym] [--norm GENS] [--maxlist N] [--k K]
"""
import argparse
import os
import struct
import subprocess
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gsgroup import Group, Orbits, read_list, group_matrix  # noqa: E402

WORK = os.path.join(os.environ.get('HADAMARD_WORK', '/tmp/hadamard-work'), 'tmp')


def build(G, a, b):
    w = G.t
    A, B = group_matrix(G, a), group_matrix(G, b)
    e = np.ones((1, w), dtype=np.int64)
    one = np.ones((1, 1), dtype=np.int64)
    top = np.hstack([one, one, e, e])
    sec = np.hstack([-one, one, e, -e])
    r3 = np.hstack([-e.T, -e.T, A, B])
    r4 = np.hstack([-e.T, e.T, -B.T, A.T])
    return np.vstack([top, sec, r3, r4])


def build_sym(G, a, b):
    w = G.t
    A, B = group_matrix(G, a), group_matrix(G, b)
    e = np.ones((1, w), dtype=np.int64)
    X = np.array([[1, 1], [1, -1]], dtype=np.int64)
    E = np.vstack([np.hstack([e, e]), np.hstack([e, -e])])
    M = np.block([[A, B], [B.T, -A]])
    return np.block([[X, E], [E.T, M]])


def check_sym(H):
    n = H.shape[0]
    return bool(np.array_equal(H, H.T) and np.array_equal(H @ H.T, n * np.eye(n, dtype=np.int64)))


def check(H):
    n = H.shape[0]
    I = np.eye(n, dtype=np.int64)
    return bool(np.array_equal(H + H.T, 2 * I) and np.array_equal(H @ H.T, n * I))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('group')
    ap.add_argument('hgens')
    ap.add_argument('--norm', default='')
    ap.add_argument('--maxlist', type=float, default=3e7)
    ap.add_argument('--k', type=int, default=2)
    ap.add_argument('--gdfs', default='/tmp/hadamard-work/bin/gdfs')
    ap.add_argument('--join4x', default='/tmp/hadamard-work/bin/join4x')
    ap.add_argument('--dry', action='store_true')
    ap.add_argument('--sym', action='store_true')
    a = ap.parse_args()
    os.makedirs(WORK, exist_ok=True)
    G = Group(a.group)
    w = G.t
    ident = '/'.join(('1' if k == 1 else ';'.join(','.join('1' if i == j else '0' for j in range(k)) for i in range(k)))
                     for n, k in G.comps)
    hg = [x for x in a.hgens.split('|') if x and x != '1'] or [ident]
    H = G.closure([G.aut_perm(x) for x in hg])
    O = Orbits(G, H)
    varA = O.variables('sym' if a.sym else 'skew')
    varB = O.variables('plain')
    sizesB = [int(sum(O.size[o] for o, c in v)) for v in varB]
    # row sum 1 count (upper bound on the b list)
    dp = {0: 1}
    for z in sizesB:
        nd = {}
        for s, c in dp.items():
            nd[s + z] = nd.get(s + z, 0) + c
            nd[s - z] = nd.get(s - z, 0) + c
        dp = nd
    rsb = -1 if a.sym else 1
    nb = dp.get(rsb, 0)
    if a.sym:
        sizesA = [int(sum(O.size[o] for o, c in v)) for v in varA]
        dpa = {0: 1}
        for z in sizesA:
            nd = {}
            for s_, c in dpa.items():
                nd[s_ + z] = nd.get(s_ + z, 0) + c
                nd[s_ - z] = nd.get(s_ - z, 0) + c
            dpa = nd
        na = dpa.get(-1, 0)
    else:
        na = 2 ** len(varA)
    print(f'G={a.group} w={w} |H|={len(H)} H-orbits {O.no} classes {O.npr} a-pairs {len(varA)} b-orbits {len(varB)} '
          f'|A|={na} |B|={nb}', file=sys.stderr)
    if a.dry or nb == 0:
        print(f'RESULT {a.group} {a.hgens} classes={O.npr} nA={na} nB={nb} skipped={"dry" if a.dry else "noB"}')
        return
    if na > a.maxlist or nb > a.maxlist:
        print(f'RESULT {a.group} {a.hgens} classes={O.npr} nA={na} nB={nb} skipped=too_large')
        return
    tag = f'szek2_{a.group.replace("^", "p")}_{os.getpid()}'
    lists = []
    for nm, var, target in (('A', varA, -1 if a.sym else 0), ('B', varB, rsb)):
        tp, lp = os.path.join(WORK, f'{tag}_{nm}.tab'), os.path.join(WORK, f'{tag}_{nm}.lst')
        O.write_tables(tp, var, 0, target, 1)
        r = subprocess.run([a.gdfs, tp, lp], capture_output=True, text=True)
        if r.returncode:
            sys.exit(r.stderr)
        n, D, vec, masks = read_list(lp)
        os.remove(tp); os.remove(lp)
        lists.append((n, vec[:, :O.npr].copy(), masks, var))     # PAF part only (no cross vector)
    D = O.npr
    fn = os.path.join(WORK, f'{tag}_join.bin')
    with open(fn, 'wb') as f:
        f.write(struct.pack('5i', lists[0][0], lists[1][0], 1, 1, D))
        f.write(np.full(D, -2, dtype=np.int16).tobytes())
        f.write(np.ascontiguousarray(lists[0][1]).tobytes())
        f.write(np.ascontiguousarray(lists[1][1]).tobytes())
        f.write(b'\0' * (2 * D))
        f.write(b'\0' * (2 * D))
    j = subprocess.run([a.join4x, fn, str(a.k)], capture_output=True, text=True)
    os.remove(fn)
    sols = [list(map(int, ln.split()[1:5])) for ln in j.stdout.splitlines() if ln.startswith('SOL')]
    tot = [ln for ln in j.stdout.splitlines() if ln.startswith('TOTAL')]
    nsol = 0
    for s in sols:
        seqs = []
        for (n, vec, masks, var), idx, base in zip(lists, s[:2], (1, 0)):
            m = int(masks[idx])
            vals = [-1 if m >> jv & 1 else 1 for jv in range(len(var))]
            x = O.seq_from(var, vals)
            if base and not a.sym:
                x[0] = 1
            seqs.append(x)
        Hm = build_sym(G, seqs[0], seqs[1]) if a.sym else build(G, seqs[0], seqs[1])
        ok = check_sym(Hm) if a.sym else check(Hm)
        nsol += ok
        print(f'SOL {a.group} {a.hgens} checked={ok}')
        print('a', ' '.join(map(str, seqs[0].tolist())))
        print('b', ' '.join(map(str, seqs[1].tolist())))
        sys.stdout.flush()
    print(f'RESULT {a.group} {a.hgens} classes={O.npr} nA={lists[0][0]} nB={lists[1][0]} {tot[0] if tot else "?"} '
          f'checked={nsol}')


if __name__ == '__main__':
    main()
