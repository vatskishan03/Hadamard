r"""Symmetric orthogonal designs OD(4t; mu, 4t - mu) in Balonin form over Z_t, invariant under a multiplier group H of
odd order (see search/symod_engine.py for the conditions and the plug-in).

1. x templates: `symtpl.c` lists the H-invariant ternary x parts of A (symmetric), B and D with their PAF vectors and
   weights; join4x finds PAF_A + 2 PAF_B + PAF_D = 0 off the peak with w_A + 2 w_B + w_D = mu (at most --bk of them).
   They are reduced modulo the multipliers of Z_t, the sign of each block and the reversal of B and of D.
2. Row-sum screen: rho with rho_A^2 + 2 rho_B^2 + rho_D^2 = 4t - mu and r_A rho_A + 2 r_B rho_B + r_D rho_D = 0, each
   rho_X a signed sum of the sizes of the free variables of block X.
3. y lists by `odmc.c` (block 5 = symmetric for A, block 2 for B and D), cached per (block type, x part, row sum),
   with their (PAF, X) vectors at the pair-orbit representatives; join4x on (A, 2B, D, {0}) with target 0.
4. Each solution is decoded with `odmc --decode` and checked: P and N symmetric with complementary supports,
   P^2 = mu I, N^2 = (4t - mu) I, PN + NP = 0.

Usage: python3 search/symod_join.py t "h1,h2" mu [--xsub d] [--bk N] [--k K] [--first] [--maxtpl N] [--dry]
"""
import argparse
import os
import struct
import subprocess
import sys

import numpy as np

WORK = os.path.join(os.environ.get('HADAMARD_WORK', '/tmp/hadamard-work'), 'tmp')


def closure(gens, t):
    H = {1}
    fr = [1]
    while fr:
        x = fr.pop()
        for g in gens:
            y = x * g % t
            if y not in H:
                H.add(y)
                fr.append(y)
    return sorted(H)


def balonin(t, codes, var):
    """the coefficient matrix of one variable (1 = x, 2 = y) in the Balonin array"""
    def circ(row):
        idx = (np.arange(t)[None, :] - np.arange(t)[:, None]) % t
        return np.asarray(row)[idx]
    sel = [np.where(np.abs(np.array(c)) == var, np.sign(np.array(c)), 0) for c in codes]
    A, B, D = (circ(s) for s in sel)
    R = np.fliplr(np.eye(t, dtype=np.int64))
    BR, RB, DR, RD = B @ R, R @ B, D @ R, R @ D
    return np.block([[-A, BR, BR, DR], [BR, RD, A, -RB], [BR, A, -RD, RB], [DR, -RB, RB, A]])


def check_design(t, codes, mu):
    P, N = balonin(t, codes, 1), balonin(t, codes, 2)
    n = 4 * t
    ok = (np.array_equal(P, P.T) and np.array_equal(N, N.T) and np.all(np.abs(P) + np.abs(N) == 1)
          and np.array_equal(P @ P, mu * np.eye(n, dtype=np.int64))
          and np.array_equal(N @ N, (n - mu) * np.eye(n, dtype=np.int64))
          and not np.any(P @ N + N @ P))
    return ok


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('t', type=int)
    ap.add_argument('gens')
    ap.add_argument('mu', type=int)
    ap.add_argument('--xsub', type=int, default=0)
    ap.add_argument('--bk', type=int, default=200000)
    ap.add_argument('--k', type=int, default=2)
    ap.add_argument('--first', action='store_true')
    ap.add_argument('--maxtpl', type=int, default=10 ** 7)
    ap.add_argument('--dry', action='store_true')
    ap.add_argument('--symtpl', default='/tmp/hadamard-work/bin/symtpl')
    ap.add_argument('--odmc', default='/tmp/hadamard-work/bin/odmc')
    ap.add_argument('--join4x', default='/tmp/hadamard-work/bin/join4x')
    a = ap.parse_args()
    t, mu = a.t, a.mu
    os.makedirs(WORK, exist_ok=True)
    gens = [int(x) for x in a.gens.split(',')]
    H = closure(gens, t)
    orbit = {}
    for g in range(t):
        orbit[g] = sorted({h * g % t for h in H})
    tag = f'sym_{t}_{mu}_{a.gens.replace(",", "_")}_{os.getpid()}'
    fb, fm = os.path.join(WORK, tag + '_tpl.bin'), os.path.join(WORK, tag + '_masks.bin')
    r = subprocess.run([a.symtpl, str(t), a.gens, str(mu), str(a.xsub), fb, fm], capture_output=True, text=True)
    print(r.stderr.strip(), file=sys.stderr)
    if r.returncode:
        sys.exit('symtpl failed')
    lists = []
    with open(fm, 'rb') as f:
        for b in range(3):
            n, = struct.unpack('i', f.read(4))
            masks = np.frombuffer(f.read(16 * n), dtype=np.uint64).reshape(n, 2)
            nv, = struct.unpack('i', f.read(4))
            var = [struct.unpack('2i', f.read(8)) for _ in range(nv)]
            lists.append((masks, var))
    os.remove(fm)
    j = subprocess.run([a.join4x, fb, str(a.bk), '8000000', '1'], capture_output=True, text=True)
    os.remove(fb)
    sols = [list(map(int, ln.split()[1:5])) for ln in j.stdout.splitlines() if ln.startswith('SOL')]
    tot = [ln for ln in j.stdout.splitlines() if ln.startswith('TOTAL')]
    print(f'x triples from join4x: {len(sols)} ({tot[0] if tot else "?"})', file=sys.stderr)

    # H-orbits and the transformations, all acting as signed permutations of the H-orbits
    oid = [-1] * t
    oreps = []
    for g in range(t):
        if oid[g] < 0:
            for x in orbit[g]:
                oid[x] = len(oreps)
            oreps.append(g)
    no = len(oreps)

    def operm(f):                                   # the permutation of H-orbits induced by g -> f(g)
        return np.array([oid[f(r0)] for r0 in oreps], dtype=np.int64)

    units = [u for u in range(1, t) if np.gcd(u, t) == 1]
    ucos, useen = [], set()                         # multipliers modulo H (H fixes every H-orbit)
    for u in units:
        k = tuple(operm(lambda g: g * u % t))
        if k not in useen:
            useen.add(k)
            ucos.append(np.array(k))
    # translations by elements fixed by H keep H-invariance (all of Z_t when H is trivial)
    fixed = [c for c in range(t) if all(h * c % t == c for h in H)]
    bd_perms = []
    for c in fixed:
        for rev in (1, t - 1):
            bd_perms.append(operm(lambda g: (rev * g + c) % t))
    pow3 = 3 ** np.arange(no - 1, -1, -1, dtype=np.int64) if no <= 39 else None
    if pow3 is None:
        sys.exit('too many H-orbits for the template code')

    def code(V):                                    # (N, no) in {-1, 0, 1} -> lexicographic int64 code
        return (V + 1).astype(np.int64) @ pow3

    def to_orbits(b, idx):                           # (N,) candidate indices -> (N, no) orbit vectors
        masks, var = lists[b]
        V = np.zeros((len(idx), no), dtype=np.int64)
        pl = masks[idx, 0].astype(np.uint64)
        mi = masks[idx, 1].astype(np.uint64)
        for jv, reps in enumerate(var):
            bp = ((pl >> np.uint64(jv)) & np.uint64(1)).astype(np.int64)
            bm = ((mi >> np.uint64(jv)) & np.uint64(1)).astype(np.int64)
            for r0 in reps:
                if r0 >= 0:
                    V[:, oid[r0]] = bp - bm
        return V

    def canon_block(V, sym):
        """per-block canonical code: min over sign (and, for B and D, translations by fixed points and reversal)"""
        perms = [np.arange(no)] if sym else bd_perms
        best = None
        for pm in perms:
            W = np.zeros_like(V)
            W[:, pm] = V
            for sg in (1, -1):
                c = code(sg * W)
                best = c if best is None else np.minimum(best, c)
        return best

    S = np.array(sols, dtype=np.int64).reshape(-1, 4)
    XA, XB, XD = to_orbits(0, S[:, 0]), to_orbits(1, S[:, 1]), to_orbits(2, S[:, 2])
    keys = None
    for pm in ucos:
        img = []
        for V, sym in ((XA, True), (XB, False), (XD, False)):
            W = np.zeros_like(V)
            W[:, pm] = V
            img.append(canon_block(W, sym))
        K = np.stack(img, axis=1)
        if keys is None:
            keys = K
        else:
            better = (K[:, 0] < keys[:, 0]) | ((K[:, 0] == keys[:, 0]) & ((K[:, 1] < keys[:, 1]) |
                     ((K[:, 1] == keys[:, 1]) & (K[:, 2] < keys[:, 2]))))
            keys[better] = K[better]
    uniq = np.unique(keys, axis=0) if len(keys) else np.zeros((0, 3), dtype=np.int64)

    def decode_code(c):                              # code -> sequence over Z_t
        digs = []
        for _ in range(no):
            digs.append(int(c % 3) - 1)
            c //= 3
        digs = digs[::-1]
        seq = [0] * t
        for i, v in enumerate(digs):
            if v:
                for g in orbit[oreps[i]]:
                    seq[g] = v
        return tuple(seq)

    templates = [tuple(decode_code(int(c)) for c in row) for row in uniq]
    print(f'{len(templates)} x templates up to equivalence', file=sys.stderr)

    # free variable sizes per block (H-invariant y parts on the complement of the x part)
    def free_sizes(seq, sym):
        done, sizes = set(), []
        for g in range(t):
            if seq[g] or g in done:
                continue
            o = set(orbit[g])
            if sym:
                o |= set(orbit[(-g) % t])
            done |= o
            sizes.append(len(o))
        return sizes

    def signed_sums(sizes):
        s = {0}
        for z in sizes:
            s = {x + z for x in s} | {x - z for x in s}
        return s

    work = []
    for tpl in templates:
        r = [sum(tpl[0]), sum(tpl[1]), sum(tpl[2])]
        att = [signed_sums(free_sizes(tpl[0], True)), signed_sums(free_sizes(tpl[1], False)),
               signed_sums(free_sizes(tpl[2], False))]
        rhos = []
        for ra in att[0]:
            for rb in att[1]:
                rem = 4 * t - mu - ra * ra - 2 * rb * rb
                if rem < 0:
                    continue
                rd = int(round(rem ** 0.5))
                for rdd in {rd, -rd}:
                    if rdd * rdd == rem and rdd in att[2] and r[0] * ra + 2 * r[1] * rb + r[2] * rdd == 0:
                        rhos.append((ra, rb, rdd))
        # equivalent row-sum vectors: negating every y part, and negating a whole block whose x part is zero
        flips = [i for i in range(3) if not any(tpl[i])]
        def orbit_rho(rv):
            out = set()
            for gsg in (1, -1):
                for mask in range(1 << len(flips)):
                    w = [gsg * x for x in rv]
                    for bit, i in enumerate(flips):
                        if mask >> bit & 1:
                            w[i] = -w[i]
                    out.add(tuple(w))
            return min(out)
        rhos = sorted(set(orbit_rho(rv) for rv in rhos))
        if rhos:
            work.append((tpl, rhos))
    print(f'{len(work)} templates pass the row-sum screen', file=sys.stderr)
    if a.dry:
        for tpl, rhos in work[:30]:
            print('TPL', [[(g, v) for g, v in enumerate(s) if v] for s in tpl], rhos[:4], len(rhos))
        return

    cache = {}

    def entries(seq, blk):
        e = [f'{blk}:{g}:{v}' for g, v in enumerate(seq) if v]
        return ','.join(e) or f'{blk}:0:0'

    def get_list(seq, blk, rs):
        key = (blk, seq, rs)
        if key in cache:
            return cache[key]
        fn = os.path.join(WORK, f'{tag}_{len(cache)}.lst')
        subprocess.run([a.odmc, str(t), a.gens, str(blk), entries(seq, blk), str(rs), fn], check=True,
                       capture_output=True)
        with open(fn, 'rb') as f:
            n, D = struct.unpack('2i', f.read(8))
            vec = np.frombuffer(f.read(2 * n * D), dtype=np.int16).reshape(n, D)
            masks = np.frombuffer(f.read(8 * n), dtype=np.uint64)
        os.remove(fn)
        cache[key] = (n, D, vec, masks)
        return cache[key]

    total = 0
    for ti, (tpl, rhos) in enumerate(work[:a.maxtpl]):
        nsol = 0
        for rho in rhos:
            L = [get_list(tpl[0], 5, rho[0]), get_list(tpl[1], 2, rho[1]), get_list(tpl[2], 3, rho[2])]
            if any(x[0] == 0 for x in L):
                continue
            D = L[0][1]
            fn = os.path.join(WORK, f'{tag}_join.bin')
            with open(fn, 'wb') as f:
                f.write(struct.pack('5i', L[0][0], L[1][0], L[2][0], 1, D))
                f.write(b'\0' * (2 * D))
                f.write(L[0][2].tobytes())
                f.write((2 * L[1][2].astype(np.int32)).astype(np.int16).tobytes())
                f.write(L[2][2].tobytes())
                f.write(b'\0' * (2 * D))
            jj = subprocess.run([a.join4x, fn, str(a.k)], capture_output=True, text=True)
            os.remove(fn)
            ss = [list(map(int, ln.split()[1:5])) for ln in jj.stdout.splitlines() if ln.startswith('SOL')]
            tt = [ln for ln in jj.stdout.splitlines() if ln.startswith('TOTAL')]
            print(f'template {ti} rho {rho}: lists {[x[0] for x in L]} {tt[0] if tt else "?"}', file=sys.stderr)
            for s in ss:
                codes = []
                for (seq, blk), lst, idx in zip(((tpl[0], 5), (tpl[1], 2), (tpl[2], 3)), L, s[:3]):
                    out = subprocess.run([a.odmc, str(t), a.gens, str(blk), entries(seq, blk), '0', '--decode',
                                          str(int(lst[3][idx]))], capture_output=True, text=True).stdout.split()
                    # odmc codes: +-2 = x part, +-3 = y part -> +-1 = x, +-2 = y
                    codes.append([int(c) // 2 if abs(int(c)) == 2 else 2 * (int(c) // 3) for c in out])
                ok = check_design(t, codes, mu)
                print(f'SOL t={t} mu={mu} template {ti} rho {rho} checked={ok}')
                for name, c in zip('ABD', codes):
                    print(name, ' '.join(map(str, c)))
                sys.stdout.flush()
                nsol += 1
        total += nsol
        if nsol and a.first:
            break
    print(f'total solutions {total}; {len(cache)} distinct y lists', file=sys.stderr)


if __name__ == '__main__':
    main()
