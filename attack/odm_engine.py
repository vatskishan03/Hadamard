"""Exhaustive search for skew-type OD(4t; 1, m, 4t - m - 1) in Goethals-Seidel form over Z_t that are invariant under a
multiplier group H (Seberry's Theorem 7 then gives skew-Hadamard matrices of order m n (n - 1) for every amicable
order n with 4t = m n; see attack/NOTES.md §31).

1. Orbits of H on Z_t; H must not map any s != 0 to -s (block 1 is skew).
2. All H-invariant b parts: block 1 a skew ternary combination of pair orbits, blocks 2-4 ternary combinations of
   orbits, total weight m, PAF sum m delta (meet in the middle on the PAF vectors at the pair-orbit representatives),
   reduced modulo permutations of blocks 2-4, negations and the multipliers of Z_t.
3. Row-sum screen: some (0, r2, r3, r4) with r2^2 + r3^2 + r4^2 = 4t - m - 1, sum_i rowsum(u_i) r_i = 0, attainable
   by H-invariant c parts on the free positions of each block.
4. `odmx` (C) on each surviving template: exact search for the c parts.

Usage: python3 attack/odm_engine.py t "h1,h2" m [--odmx PATH] [--maxtpl N] [--maxsol K]
"""
import argparse
import itertools
import subprocess
import sys
from collections import defaultdict


def closure(gens, t):
    H = {1}
    frontier = [1]
    while frontier:
        x = frontier.pop()
        for g in gens:
            y = x * g % t
            if y not in H:
                H.add(y)
                frontier.append(y)
    return sorted(H)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('t', type=int)
    ap.add_argument('gens')
    ap.add_argument('m', type=int)
    ap.add_argument('--odmx', default='/tmp/hadamard-work/bin/odmx')
    ap.add_argument('--maxtpl', type=int, default=10 ** 6)
    ap.add_argument('--maxsol', type=int, default=3)
    ap.add_argument('--dry', action='store_true')
    a = ap.parse_args()
    t, m = a.t, a.m
    H = closure([int(x) for x in a.gens.split(',')], t)
    orb = [-1] * t
    orbits = []
    for g in range(t):
        if orb[g] < 0:
            o = sorted({h * g % t for h in H})
            for x in o:
                orb[x] = len(orbits)
            orbits.append(o)
    for g in range(1, t):
        if orb[g] == orb[(-g) % t]:
            sys.exit(f'H maps {g} to -{g}: not admissible')
    no = len(orbits)
    neg = [orb[(-orbits[i][0]) % t] for i in range(no)]
    pairs = [i for i in range(1, no) if i < neg[i]]          # pair orbits (i, neg[i]) with i < neg[i]
    reps = [orbits[i][0] for i in pairs]                       # PAF representatives (one per +- pair of orbits)
    print(f't={t} |H|={len(H)} orbits {no} (sizes {sorted(len(o) for o in orbits)}), pair orbits {len(pairs)}, m={m}',
          file=sys.stderr)

    def seq_from(vals):                                      # vals: per orbit value
        s = [0] * t
        for i, v in enumerate(vals):
            if v:
                for x in orbits[i]:
                    s[x] = v
        return s

    # T[o1][o2][k] = #{x in orbit o1 : x + reps[k] in orbit o2}
    T = [[[0] * len(reps) for _ in range(no)] for _ in range(no)]
    for o1 in range(no):
        for x in orbits[o1]:
            for k, r in enumerate(reps):
                T[o1][orb[(x + r) % t]][k] += 1

    def paf_of(nz):                                          # nz: list of (orbit, value)
        p = [0] * len(reps)
        for o1, v1 in nz:
            for o2, v2 in nz:
                row = T[o1][o2]
                if v1 * v2 == 1:
                    for k in range(len(reps)):
                        p[k] += row[k]
                else:
                    for k in range(len(reps)):
                        p[k] -= row[k]
        return tuple(p)

    sizes = [len(o) for o in orbits]

    def gen(items, budget, skew):
        """weight-bounded ternary assignments over the given orbit items; skew: items are pair orbits"""
        out = []
        def rec(j, w, nz):
            if j == len(items):
                out.append((w, nz))
                return
            rec(j + 1, w, nz)
            i = items[j]
            z = sizes[i] * (2 if skew else 1)
            if w + z <= budget:
                for v in (1, -1):
                    nz2 = nz + ([(i, v), (neg[i], -v)] if skew else [(i, v)])
                    rec(j + 1, w + z, nz2)
        rec(0, 0, [])
        return out

    def to_full(nz):
        full = [0] * no
        for i, v in nz:
            full[i] = v
        return tuple(full)

    def block1_cands():
        return [(w, paf_of(nz), to_full(nz)) for w, nz in gen(pairs, m, True)]

    def blockk_cands():
        return [(w, paf_of(nz), to_full(nz)) for w, nz in gen(list(range(no)), m, False)]

    B1 = block1_cands()
    Bk = blockk_cands()
    print(f'b candidates: block 1 {len(B1)}, blocks 2-4 {len(Bk)}', file=sys.stderr)
    # meet in the middle: (1, 2) vs (3, 4)
    left = defaultdict(list)
    for w1, p1, v1 in B1:
        for w2, p2, v2 in Bk:
            if w1 + w2 <= m:
                left[(w1 + w2, tuple(x + y for x, y in zip(p1, p2)))].append((v1, v2))
    units = [u for u in range(1, t) if all(u * k % t for k in range(1, t))]
    seen = set()
    templates = []
    for w3, p3, v3 in Bk:
        for w4, p4, v4 in Bk:
            w = w3 + w4
            if w > m:
                continue
            key = (m - w, tuple(-(x + y) for x, y in zip(p3, p4)))
            for v1, v2 in left.get(key, ()):
                tpl = (v1, v2, v3, v4)
                # canonical form: multipliers (act on orbits), block 2-4 permutations, negations
                best = None
                for u in units:
                    perm = [orb[orbits[i][0] * u % t] for i in range(no)]
                    imgs = []
                    for blk in tpl:
                        nb = [0] * no
                        for i, v in enumerate(blk):
                            nb[perm[i]] = v
                        imgs.append(tuple(nb))
                    b1 = min(imgs[0], tuple(-x for x in imgs[0]))
                    rest = sorted(min(b, tuple(-x for x in b)) for b in imgs[1:])
                    k = (b1,) + tuple(rest)
                    if best is None or k < best:
                        best = k
                if best in seen:
                    continue
                seen.add(best)
                templates.append(best)
    print(f'{len(templates)} b templates up to equivalence', file=sys.stderr)
    cw = 4 * t - m - 1
    ok_tpl = []
    for tpl in templates:
        seqs = [seq_from(b) for b in tpl]
        if seqs[0][0] != 0:
            continue
        ru = [sum(s) for s in seqs]
        # attainable row sums of H-invariant c parts per block
        att = []
        for i, s in enumerate(seqs):
            free = [o for o in range(no) if not s[orbits[o][0]] and not (i == 0 and o == orb[0])]
            if i == 0:
                att.append({0})
                continue
            sizes = [len(orbits[o]) for o in free]
            sums = {0}
            for z in sizes:
                sums = {x + z for x in sums} | {x - z for x in sums}
            att.append(sums)
        feas = []
        lim = int(cw ** 0.5) + 1
        for r2 in range(-lim, lim + 1):
            if r2 not in att[1]:
                continue
            for r3 in range(-lim, lim + 1):
                if r3 not in att[2]:
                    continue
                rem = cw - r2 * r2 - r3 * r3
                if rem < 0:
                    continue
                r4 = int(round(rem ** 0.5))
                for r4s in {r4, -r4}:
                    if r4s * r4s == rem and r4s in att[3] and ru[1] * r2 + ru[2] * r3 + ru[3] * r4s == 0:
                        feas.append((r2, r3, r4s))
        if feas:
            ok_tpl.append((tpl, feas))
    print(f'{len(ok_tpl)} templates pass the row-sum screen', file=sys.stderr)
    found = 0
    for k, (tpl, feas) in enumerate(ok_tpl[:a.maxtpl]):
        parts = []
        for i, b in enumerate(tpl):
            s = seq_from(b)
            for x in range(t):
                if s[x]:
                    parts.append(f'{i + 1}:{x}:{s[x]}')
        tstr = ','.join(parts)
        if a.dry:
            print(tstr, feas)
            continue
        r = subprocess.run([a.odmx, str(t), ','.join(map(str, H[1:2] or [1])) if False else a.gens, tstr,
                            str(a.maxsol)], capture_output=True, text=True)
        nsol = r.stdout.count('SOL')
        last = [ln for ln in r.stderr.splitlines() if ln.startswith('rows')]
        print(f'template {k}: {tstr}  -> {nsol} solutions  ({len(last)} row vectors)', file=sys.stderr)
        if nsol:
            found += nsol
            print(f'# template {tstr}')
            print(r.stdout, end='')
            sys.stdout.flush()
    print(f'total solutions {found}', file=sys.stderr)


if __name__ == '__main__':
    main()
