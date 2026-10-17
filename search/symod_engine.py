r"""Exhaustive search for symmetric orthogonal designs OD(4t; mu, 4t - mu) in Balonin (propus) form over Z_t that are
invariant under a multiplier group H.

Balonin array (A symmetric circulant, B and D any circulants, R the back identity):

    [ -A    BR    BR    DR ]
    [ BR    RD    A    -RB ]
    [ BR    A    -RD    RB ]
    [ DR   -RB    RB    A  ]

is symmetric, and orthogonal whenever AA^T + 2BB^T + DD^T is scalar. With every block split into an x part and a
y part (disjoint supports covering Z_t), A_x and A_y symmetric, xP + yN is a symmetric OD of type (mu, 4t - mu) iff
    PAF(A_x) + 2 PAF(B_x) + PAF(D_x) = mu delta,
    PAF(A_y) + 2 PAF(B_y) + PAF(D_y) = (4t - mu) delta,
    X_A(s) + 2 X_B(s) + X_D(s) = 0 for s != 0,   X(s) = sum_g x(g) [y(g + s) + y(g - s)].
With a symmetric Hadamard matrix of order h, 4t = mu h, the plug-in P (x) J_{h-1} + N (x) D_h (D_h the core) is a
symmetric Hadamard matrix of order mu h (h - 1) (search/NOTES.md §31, "weighted variant").

H may be any subgroup of Z_t^*; A is invariant under <H, -1>, B and D under H. All conditions are invariant under
<H, -1>, so they are imposed at one representative of each <H, -1>-orbit on Z_t \ {0}.

1. x templates: ternary orbit combinations with w_A + 2 w_B + w_D = mu and PAF sum mu delta (meet in the middle on
   the PAF vectors), reduced modulo multipliers of Z_t, the sign of each block and the reversal of B and of D.
2. Row-sum screen: rho with Q(rho) = 4t - mu and B(r, rho) = 0 (Q = diag(1, 2, 1)), each rho_X a signed sum of the
   free orbit sizes of block X.
3. For each template and row-sum vector: the lists of y parts with their (PAF, X) vectors, joined exactly
   (vec_A + 2 vec_B + vec_D = 0).

Usage: python3 search/symod_engine.py t "h1,h2" mu [--xsub d] [--maxsol K] [--dry] [--maxlist N]
"""
import argparse
import sys
from collections import defaultdict

import numpy as np


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


def orbits_of(G, t):
    orb = [-1] * t
    obs = []
    for g in range(t):
        if orb[g] < 0:
            o = sorted({h * g % t for h in G})
            for x in o:
                orb[x] = len(obs)
            obs.append(o)
    return orb, obs


class Structure:
    def __init__(self, t, gens):
        self.t = t
        self.H = closure(gens, t)
        self.Hpm = closure(gens + [t - 1], t)
        self.orbH, self.obH = orbits_of(self.H, t)
        self.orbP, self.obP = orbits_of(self.Hpm, t)
        self.reps = [o[0] for o in self.obP if o[0] != 0]
        self.K = len(self.reps)
        # tables over H-orbits: Tp[o1, o2, k] = #{g in o1 : g + s_k in o2}, Tm with g - s_k
        n = len(self.obH)
        Tp = np.zeros((n, n, self.K), dtype=np.int64)
        Tm = np.zeros((n, n, self.K), dtype=np.int64)
        for o1, o in enumerate(self.obH):
            for g in o:
                for k, s in enumerate(self.reps):
                    Tp[o1, self.orbH[(g + s) % t], k] += 1
                    Tm[o1, self.orbH[(g - s) % t], k] += 1
        self.Tp, self.Tm = Tp, Tm
        # <H,-1>-orbits as unions of H-orbits
        self.pm_members = [sorted({self.orbH[g] for g in o}) for o in self.obP]
        self.sizeH = np.array([len(o) for o in self.obH], dtype=np.int64)
        self.sizeP = np.array([len(o) for o in self.obP], dtype=np.int64)

    def lift(self, vals, block):
        """orbit values (A: over <H,-1>-orbits; B, D: over H-orbits) -> values over H-orbits"""
        if block == 'A':
            out = np.zeros(len(self.obH), dtype=np.int64)
            for i, v in enumerate(vals):
                for j in self.pm_members[i]:
                    out[j] = v
            return out
        return np.asarray(vals, dtype=np.int64)

    def seq(self, vals, block):
        hv = self.lift(vals, block)
        s = np.zeros(self.t, dtype=np.int64)
        for j, o in enumerate(self.obH):
            s[o] = hv[j]
        return s


def paf_vecs(S, F):
    """F: (N, nH) values over H-orbits -> (N, K) PAF at the representatives"""
    return np.einsum('na,abk,nb->nk', F, S.Tp, F, optimize=True)


def ternary_cands(S, block, maxw, items):
    """all ternary assignments on the given orbit items with weight <= maxw: (vals over orbits of the block, weight,
    PAF vector)"""
    nob = len(S.obP) if block == 'A' else len(S.obH)
    sizes = S.sizeP if block == 'A' else S.sizeH
    out_v, out_w = [], []
    def rec(j, w, cur):
        if j == len(items):
            out_v.append(tuple(cur)); out_w.append(w)
            return
        rec(j + 1, w, cur)
        i = items[j]
        if w + sizes[i] <= maxw:
            for v in (1, -1):
                cur[i] = v
                rec(j + 1, w + int(sizes[i]), cur)
            cur[i] = 0
    rec(0, 0, [0] * nob)
    V = np.array(out_v, dtype=np.int64).reshape(len(out_v), nob)
    F = np.stack([S.lift(v, block) for v in V]) if len(V) else np.zeros((0, len(S.obH)), dtype=np.int64)
    P = paf_vecs(S, F) if len(F) else np.zeros((0, S.K), dtype=np.int64)
    return V, np.array(out_w, dtype=np.int64), P


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('t', type=int)
    ap.add_argument('gens')
    ap.add_argument('mu', type=int)
    ap.add_argument('--xsub', type=int, default=0, help='restrict the x part to the subgroup of this order')
    ap.add_argument('--maxsol', type=int, default=4)
    ap.add_argument('--maxlist', type=int, default=30000000)
    ap.add_argument('--dry', action='store_true')
    ap.add_argument('--first', action='store_true', help='stop after the first template with solutions')
    a = ap.parse_args()
    t, mu = a.t, a.mu
    gens = [int(x) for x in a.gens.split(',')] if a.gens not in ('1', '') else [1]
    S = Structure(t, gens)
    print(f't={t} |H|={len(S.H)} |H+-|={len(S.Hpm)} H-orbits {len(S.obH)} <H,-1>-orbits {len(S.obP)} '
          f'reps {S.K} mu={mu}', file=sys.stderr)
    if a.xsub:
        step = t // a.xsub
        okH = [all(g % step == 0 for g in o) for o in S.obH]
        okP = [all(g % step == 0 for g in o) for o in S.obP]
    else:
        okH = [True] * len(S.obH)
        okP = [True] * len(S.obP)
    VA, WA, PA = ternary_cands(S, 'A', mu, [i for i in range(len(S.obP)) if okP[i]])
    VB, WB, PB = ternary_cands(S, 'B', mu // 2, [i for i in range(len(S.obH)) if okH[i]])
    VD, WD, PD = ternary_cands(S, 'D', mu, [i for i in range(len(S.obH)) if okH[i]])
    print(f'x candidates: A {len(VA)}, B {len(VB)}, D {len(VD)}', file=sys.stderr)
    dD = defaultdict(list)
    for i in range(len(VD)):
        dD[(int(WD[i]),) + tuple(int(x) for x in PD[i])].append(i)
    # templates
    units = [u for u in range(1, t) if np.gcd(u, t) == 1]
    def img(vals, block, u):
        ob, orb = (S.obP, S.orbP) if block == 'A' else (S.obH, S.orbH)
        out = [0] * len(ob)
        for i, v in enumerate(vals):
            if v:
                out[orb[ob[i][0] * u % t]] = int(v)
        return tuple(out)
    def canon_block(vals, block):
        c = [vals, tuple(-x for x in vals)]
        if block != 'A':
            r = img(vals, block, t - 1)
            c += [r, tuple(-x for x in r)]
        return min(c)
    seen = set()
    templates = []
    for ia in range(len(VA)):
        for ib in range(len(VB)):
            w = int(WA[ia]) + 2 * int(WB[ib])
            if w > mu:
                continue
            need = (mu - w,) + tuple(int(x) for x in (-PA[ia] - 2 * PB[ib]))
            for idd in dD.get(need, ()):
                tpl = (tuple(VA[ia]), tuple(VB[ib]), tuple(VD[idd]))
                best = None
                for u in units:
                    k = (canon_block(img(tpl[0], 'A', u), 'A'), canon_block(img(tpl[1], 'B', u), 'B'),
                         canon_block(img(tpl[2], 'D', u), 'D'))
                    if best is None or k < best:
                        best = k
                if best not in seen:
                    seen.add(best)
                    templates.append(best)
    print(f'{len(templates)} x templates up to equivalence', file=sys.stderr)
    # row-sum screen
    def signed_sums(sizes):
        s = {0}
        for z in sizes:
            s = {x + z for x in s} | {x - z for x in s}
        return s
    def rowsum(vals, block):
        sizes = S.sizeP if block == 'A' else S.sizeH
        return int(sum(int(v) * int(sizes[i]) for i, v in enumerate(vals)))
    def free(vals):
        return [i for i, v in enumerate(vals) if v == 0]
    work = []
    for tpl in templates:
        r = [rowsum(tpl[0], 'A'), rowsum(tpl[1], 'B'), rowsum(tpl[2], 'D')]
        att = [signed_sums([int(S.sizeP[i]) for i in free(tpl[0])]),
               signed_sums([int(S.sizeH[i]) for i in free(tpl[1])]),
               signed_sums([int(S.sizeH[i]) for i in free(tpl[2])])]
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
        if rhos:
            work.append((tpl, r, sorted(set(rhos))))
    print(f'{len(work)} templates pass the row-sum screen', file=sys.stderr)
    if a.dry:
        for tpl in templates[:40]:
            print('ALL', tpl, 'r', [rowsum(tpl[0], 'A'), rowsum(tpl[1], 'B'), rowsum(tpl[2], 'D')])
    if a.dry:
        for tpl, r, rhos in work[:50]:
            print('TPL', tpl, 'r', r, 'rhos', rhos[:6], len(rhos))
        return
    total = 0
    for ti, (tpl, r, rhos) in enumerate(work):
        nsol = solve_template(S, tpl, rhos, a, ti)
        total += nsol
        if nsol and a.first:
            break
    print(f'total solutions {total}', file=sys.stderr)


def ylist(S, xvals, block, rho, maxlist):
    """all +-1 y parts on the free orbits of a block with row sum rho: (Y over orbits, vec = (PAF_y, X))"""
    sizes = S.sizeP if block == 'A' else S.sizeH
    fr = [i for i, v in enumerate(xvals) if v == 0]
    f = len(fr)
    if 2 ** f > maxlist * 64:
        raise RuntimeError(f'block {block}: {f} free orbits, too many')
    xH = S.lift(xvals, block)
    # linear X part: L[j, k] = sum_{o1} x(o1) (Tp + Tm)[o1, j, k] over H-orbits j
    L = np.einsum('a,ajk->jk', xH, S.Tp + S.Tm)
    outY, outV = [], []
    chunk = 1 << min(f, 20)
    fsz = np.array([int(sizes[i]) for i in fr], dtype=np.int64)
    for start in range(0, 1 << f, chunk):
        idx = np.arange(start, min(1 << f, start + chunk), dtype=np.int64)
        bits = ((idx[:, None] >> np.arange(f)[None, :]) & 1).astype(np.int64)
        Y = 1 - 2 * bits
        keep = (Y @ fsz) == rho
        Y = Y[keep]
        if not len(Y):
            continue
        full = np.zeros((len(Y), len(sizes)), dtype=np.int64)
        full[:, fr] = Y
        FH = np.stack([S.lift(v, block) for v in full]) if block == 'A' else full
        P = paf_vecs(S, FH)
        X = FH @ L
        outY.append(full); outV.append(np.concatenate([P, X], axis=1))
    if not outY:
        return np.zeros((0, len(sizes)), dtype=np.int64), np.zeros((0, 2 * S.K), dtype=np.int64)
    return np.concatenate(outY), np.concatenate(outV)


def solve_template(S, tpl, rhos, a, ti):
    nsol = 0
    for rho in rhos:
        LA = ylist(S, tpl[0], 'A', rho[0], a.maxlist)
        LB = ylist(S, tpl[1], 'B', rho[1], a.maxlist)
        LD = ylist(S, tpl[2], 'D', rho[2], a.maxlist)
        if not (len(LA[0]) and len(LB[0]) and len(LD[0])):
            continue
        # join: vA + 2 vB + vD = 0; hash D, loop over the smaller of A x B
        dD = defaultdict(list)
        for i, v in enumerate(LD[1]):
            dD[v.tobytes()].append(i)
        VA, VB = LA[1], 2 * LB[1]
        for ia in range(len(VA)):
            need = -(VA[ia][None, :] + VB)
            for ib in range(len(VB)):
                hit = dD.get(need[ib].tobytes())
                if hit:
                    for idd in hit:
                        nsol += 1
                        if nsol <= a.maxsol:
                            report(S, tpl, (LA[0][ia], LB[0][ib], LD[0][idd]), ti, rho)
        print(f'template {ti} rho {rho}: lists {len(LA[0])} {len(LB[0])} {len(LD[0])}, solutions so far {nsol}',
              file=sys.stderr)
    return nsol


def report(S, tpl, ys, ti, rho):
    seqs = []
    for blk, xv, yv in zip('ABD', tpl, ys):
        xs = S.seq(xv, blk)
        ysq = S.seq(yv, blk)
        seqs.append((xs, ysq))
    # codes: +-1 = +-x, +-2 = +-y
    codes = [[int(xs[g]) if xs[g] else 2 * int(ysq[g]) for g in range(S.t)] for xs, ysq in seqs]
    print('SOL template', ti, 'rho', rho)
    for name, c in zip('ABD', codes):
        print(name, ' '.join(str(x) for x in c))
    sys.stdout.flush()


if __name__ == '__main__':
    main()
