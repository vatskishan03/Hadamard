r"""Skew Goethals-Seidel quadruples (A, A o mu, A o mu^2, D): three skew blocks in one multiplier orbit.

Let G be abelian of odd order v, H a group of automorphisms of odd order, mu an automorphism with mu^3 in H, and
K = <H, mu> with -1 not in K.  If A is of skew type and H-invariant and D is K-invariant, then the blocks
A, A o mu, A o mu^2, D have periodic autocorrelations summing to

    T_A(s) + P_D(s),    T_A(s) = P_A(s) + P_A(mu s) + P_A(mu^2 s),

and both terms are K-invariant (and even in s).  So the quadruple is complementary iff T_A = -P_D on the
<K, -1>-classes: an exact join of two lists, the 2^(pairs) skew blocks A and the K-invariant blocks D with row sum d.
The row sums are (1, 1, 1, d), so d^2 = 4v - 3, i.e. v = k^2 + k + 1 and d = 2k + 1.  The first block is of skew type,
so the Goethals-Seidel array is a skew-Hadamard matrix of order 4v.  (Spin structures, attack/spin_skew.py, permute
the three non-skew blocks instead and need 1 + 3b^2 = 4v.)

Usage: python3 search/tri_skew.py GROUP "HGENS" "MU" [--maxlist N] [--k K] [--dry]
       python3 search/tri_skew.py --sweep v [--maxlist N]     (cyclic Z_v: every odd H and every admissible mu)
Automorphisms are given as in gsgroup.py ('|' between generators); for cyclic groups a unit u.
"""
import argparse
import itertools
import json
import math
import os
import struct
import subprocess
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'verify'))
from gsgroup import Group, Orbits, read_list, group_matrix  # noqa: E402
import workdir as W  # noqa: E402

WORK = str(W.WORK / 'tmp')


def count_rowsum(sizes, target):
    dp = {0: 1}
    for z in sizes:
        nd = {}
        for s, c in dp.items():
            nd[s + z] = nd.get(s + z, 0) + c
            nd[s - z] = nd.get(s - z, 0) + c
        dp = nd
    return dp.get(target, 0)


def build_gs(G, seqs):
    A, B, C, D = (group_matrix(G, x) for x in seqs)
    R = np.zeros((G.t, G.t), dtype=np.int64)
    R[np.arange(G.t), G.neg] = 1
    return np.block([[A, B @ R, C @ R, D @ R],
                     [-B @ R, A, D.T @ R, -C.T @ R],
                     [-C @ R, -D.T @ R, A, B.T @ R],
                     [-D @ R, C.T @ R, -B.T @ R, A]])


def run(spec, hgens, mugen, maxlist=5e7, k=2, dry=False, quiet=False):
    G = Group(spec)
    v = G.t
    d = math.isqrt(4 * v - 3)
    if d * d != 4 * v - 3:
        return dict(status='no-rowsum', v=v)
    hg = [x for x in hgens.split('|') if x and x != '1']
    Hp = G.closure([G.aut_perm(x) for x in hg]) if hg else G.closure([np.arange(G.t)])
    mu = G.aut_perm(mugen)
    Kp = G.closure([G.aut_perm(x) for x in hg] + [mu])
    if len(Kp) != 3 * len(Hp):
        return dict(status='mu-order', v=v)
    OH, OK = Orbits(G, Hp), Orbits(G, Kp)      # raises if -1 in K (an orbit closed under negation)
    varA, varD = OH.variables('skew'), OK.variables('plain')
    sizesD = [int(sum(OK.size[o] for o, c in var)) for var in varD]
    nA, nD = 2 ** len(varA), count_rowsum(sizesD, d)
    info = dict(v=v, spec=spec, H=hgens, mu=mugen, h=len(Hp), K=len(Kp), OH=OH.no, OK=OK.no, classes=OK.npr,
                pairs=len(varA), nA=nA, nD=nD, d=d)
    if dry or nD == 0 or nA > maxlist or nD > maxlist or len(varA) > 63 or len(varD) > 63:
        info['status'] = 'dry' if dry else ('noD' if nD == 0 else 'too_large')
        return info
    # map: element -> index of its <H,-1>-class among OH.reps
    hcls = {}
    for i, r in enumerate(OH.reps):
        o = OH.orb[r]
        hcls[o] = i
        hcls[int(OH.negorb[o])] = i
    mupow = [np.arange(G.t), mu, mu[mu]]
    idx = np.array([[hcls[OH.orb[mp[s]]] for mp in mupow] for s in OK.reps], dtype=np.int64)   # (npr_K, 3)
    os.makedirs(WORK, exist_ok=True)
    tag = f'tri_{spec.replace("^", "p")}_{os.getpid()}'
    lists = []
    for nm, O, var, target in (('A', OH, varA, 0), ('D', OK, varD, d)):
        tp, lp = os.path.join(WORK, f'{tag}_{nm}.tab'), os.path.join(WORK, f'{tag}_{nm}.lst')
        O.write_tables(tp, var, 0, target, 1)
        r = subprocess.run([W.bin('gdfs'), tp, lp], capture_output=True, text=True)
        if r.returncode:
            sys.exit(r.stderr)
        n, Dm, vec, masks = read_list(lp)
        os.remove(tp)
        os.remove(lp)
        if nm == 'A':
            P = vec[:, :OH.npr].astype(np.int32)
            T = P[:, idx[:, 0]] + P[:, idx[:, 1]] + P[:, idx[:, 2]]
            lists.append((n, T.astype(np.int16), masks, var, O))
        else:
            lists.append((n, vec[:, :OK.npr].copy(), masks, var, O))
    C = OK.npr
    fn = os.path.join(WORK, f'{tag}_join.bin')
    with open(fn, 'wb') as f:
        f.write(struct.pack('5i', lists[0][0], lists[1][0], 1, 1, C))
        f.write(np.zeros(C, dtype=np.int16).tobytes())
        f.write(np.ascontiguousarray(lists[0][1]).tobytes())
        f.write(np.ascontiguousarray(lists[1][1]).tobytes())
        f.write(b'\0' * (2 * C))
        f.write(b'\0' * (2 * C))
    j = subprocess.run([W.bin('join4x'), fn, str(k)], capture_output=True, text=True)
    os.remove(fn)
    sols = [list(map(int, ln.split()[1:5])) for ln in j.stdout.splitlines() if ln.startswith('SOL')]
    info['joined'] = len(sols)
    info['checked'] = 0
    hits = []
    for s in sols:
        seqs = []
        for (n, vec, masks, var, O), ix, skew in zip(lists, s[:2], (True, False)):
            m = int(masks[ix])
            vals = [-1 if m >> jv & 1 else 1 for jv in range(len(var))]
            x = O.seq_from(var, vals)
            if skew:
                x[0] = 1
            seqs.append(x)
        A, D = seqs
        quad = [A, A[mu], A[mu[mu]], D]          # (A o mu)(g) = A(mu g)
        Hm = build_gs(G, quad)
        n4 = 4 * v
        ok = bool(np.array_equal(Hm + Hm.T, 2 * np.eye(n4, dtype=np.int64))
                  and np.array_equal(Hm @ Hm.T, n4 * np.eye(n4, dtype=np.int64)))
        info['checked'] += ok
        if ok:
            hits.append([''.join('+' if z > 0 else '-' for z in x) for x in quad])
            if not quiet:
                print('HIT', json.dumps(dict(v=v, spec=spec, H=hgens, mu=mugen, seqs=hits[-1])), flush=True)
    info['status'] = 'done'
    info['hits'] = hits
    return info


def sweep(v, maxlist, k):
    units = [u for u in range(1, v) if math.gcd(u, v) == 1]

    def closure(gens):
        S = {1}
        fr = [1]
        while fr:
            x = fr.pop()
            for g in gens:
                y = x * g % v
                if y not in S:
                    S.add(y)
                    fr.append(y)
        return frozenset(S)

    def order(u):
        x, e = u, 1
        while x != 1:
            x, e = x * u % v, e + 1
        return e
    odd = [u for u in units if order(u) % 2 == 1]
    cyc = {}
    for u in odd:
        cyc.setdefault(closure([u]), u)
    subs = {c: [g] for c, g in cyc.items()}
    grow = True
    while grow:
        grow = False
        for (c1, g1), (c2, g2) in itertools.product(list(subs.items()), list(cyc.items())):
            if c2 <= c1:
                continue
            c = closure(g1 + [g2])
            if c not in subs:
                subs[c] = g1 + [g2]
                grow = True
    for Hs, gens in sorted(subs.items(), key=lambda t: -len(t[0])):
        seen = set()
        for mu in units:
            if mu in Hs or pow(mu, 3, v) not in Hs:
                continue
            Ks = closure(list(gens) + [mu])
            if v - 1 in Ks or Ks in seen:
                continue
            seen.add(Ks)
            for m in (mu, pow(mu, 2, v)):              # mu and mu^-1 (mod H) give mirror structures; try both
                try:
                    info = run(str(v), '|'.join(map(str, gens)) or '1', str(m), maxlist, k, quiet=False)
                except ValueError:
                    continue
                print(json.dumps({k_: info[k_] for k_ in info if k_ != 'hits'}), flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('group', nargs='?')
    ap.add_argument('hgens', nargs='?', default='1')
    ap.add_argument('mu', nargs='?')
    ap.add_argument('--maxlist', type=float, default=5e7)
    ap.add_argument('--k', type=int, default=2)
    ap.add_argument('--dry', action='store_true')
    ap.add_argument('--sweep', type=int, default=0)
    a = ap.parse_args()
    if a.sweep:
        sweep(a.sweep, a.maxlist, a.k)
        return
    info = run(a.group, a.hgens, a.mu, a.maxlist, a.k, a.dry)
    print(json.dumps({k_: info[k_] for k_ in info if k_ != 'hits'}))


if __name__ == '__main__':
    main()
