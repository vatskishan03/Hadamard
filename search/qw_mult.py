#!/usr/bin/env python3
"""Quasi-Williamson matrices of odd order w (Balonin-Djokovic, arXiv:1508.00640, (27)-(28)): circulant A, B, C, D with
sum PAF = 4w delta and a*b - b*a + c*d - d*c = 0 (cross-correlations), constant on the orbits of a multiplier group H
of odd order (so -1 is not in H and the matrices need not be symmetric). Pair join: (a, b) keyed by
(PAF_a + PAF_b, antisymmetric part of a*b) against (c, d). Every odd-order multiplier group with at most MAXO orbits.
Quasi-Williamson matrices of order w give an NG pair of length 2w (search/miyamoto_blocks.py ng_williamson).

    python3 search/qw_mult.py w1,w2,... [MAXO]     # e.g. 67,93,107,127,235 15
"""
import itertools, sys, math, time
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from collections import defaultdict

def closure(gens, w):
    S = {1}; todo = [1]
    while todo:
        x = todo.pop()
        for g in gens:
            y = x * g % w
            if y not in S: S.add(y); todo.append(y)
    return sorted(S)

def odd_subgroups(w):
    units = [x for x in range(1, w) if math.gcd(x, w) == 1]
    seen = {}
    for g in units:
        H = tuple(closure([g], w))
        if len(H) % 2 == 1 and H not in seen: seen[H] = g
    for g1, g2 in itertools.combinations(list(seen.values()), 2):
        H = tuple(closure([g1, g2], w))
        if len(H) % 2 == 1 and H not in seen: seen[H] = (g1, g2)
    return seen

def search(w, H, maxpairs=3e7, log=print):
    seen, orbits = set(), []
    for x in range(w):
        if x not in seen:
            o = sorted({x * h % w for h in H}); seen |= set(o); orbits.append(o)
    O = len(orbits)
    reps = [o[0] for o in orbits]
    ind = np.zeros((O, w), np.int64)
    for i, o in enumerate(orbits): ind[i, o] = 1
    S = np.array(list(itertools.product((1, -1), repeat=O)), np.int64) @ ind
    F = np.fft.fft(S.astype(np.float64), axis=1)
    paf = np.rint(np.fft.ifft(np.abs(F) ** 2, axis=1).real).astype(np.int64)[:, reps]
    rs = S.sum(1)
    neg = [orbits.index(next(o for o in orbits if (-r) % w in o)) for r in reps]
    quads = sorted({tuple(sorted(abs(x) for x in c)) for c in itertools.product(range(1, math.isqrt(4 * w) + 1, 2), repeat=4)
                    if sum(x * x for x in c) == 4 * w})
    target = np.zeros(O, np.int64); target[0] = 4 * w
    by_r = defaultdict(list)
    for i, r in enumerate(rs): by_r[abs(int(r))].append(i)
    def cross_anti(I, J):
        # a*b (k) = sum_i a_i b_{i+k} = IFFT(conj(FA) FB); antisym part on orbit reps: f(k) - f(-k)
        X = np.rint(np.fft.ifft(np.conj(F[I]) * F[J], axis=1).real).astype(np.int64)
        Xr = X[:, reps]
        return Xr - Xr[:, neg]
    for quad in quads:
        for perm in set(itertools.permutations(quad)):
            ra, rb, rc, rd = perm
            if (ra, rb) > (rc, rd): continue
            LA, LB, LC, LD = (np.array(by_r[r], np.int64) for r in perm)
            if len(LA) * len(LB) > maxpairs or len(LC) * len(LD) > maxpairs:
                log(f"  rows {perm}: pairs too many ({len(LA)*len(LB):.2e}); skipped"); continue
            IA, IB = np.repeat(LA, len(LB)), np.tile(LB, len(LA))
            IC, ID = np.repeat(LC, len(LD)), np.tile(LD, len(LC))
            if len(IA) == 0 or len(IC) == 0: continue
            KL = np.concatenate([paf[IA] + paf[IB], cross_anti(IA, IB)], axis=1)
            KR = np.concatenate([target - paf[IC] - paf[ID], -cross_anti(IC, ID)], axis=1)
            d = {}
            for idx, row in enumerate(KL): d.setdefault(row.tobytes(), idx)
            for idx, row in enumerate(KR):
                k = d.get(row.tobytes())
                if k is not None:
                    A, B, C, D = S[IA[k]], S[IB[k]], S[IC[idx]], S[ID[idx]]
                    return (A, B, C, D), perm
    return None, None

def check(W):
    A, B, C, D = W; w = len(A)
    circ = lambda x: np.array([np.roll(x, i) for i in range(w)])
    a, b, c, d = map(circ, W)
    M = np.block([[a, b, c, d], [-b, a, -d, c], [-c.T, d.T, a.T, -b.T], [-d.T, -c.T, b.T, a.T]])
    return np.array_equal(M @ M.T, 4 * w * np.eye(4 * w, dtype=np.int64))

if __name__ == '__main__':
    maxO = int(sys.argv[2]) if len(sys.argv) > 2 else 12
    for w in [int(x) for x in sys.argv[1].split(',')]:
        subs = odd_subgroups(w)
        cand = []
        for H, g in subs.items():
            seen, O = set(), 0
            for x in range(w):
                if x not in seen: seen |= {x * h % w for h in H}; O += 1
            cand.append((O, len(H), g, H))
        cand.sort()
        print(f"=== w = {w}: odd subgroups {[(O, h) for O, h, g, H in cand]}", flush=True)
        for O, h, g, H in cand:
            if O > maxO or O < 3: continue
            t = time.time()
            W, perm = search(w, H)
            print(f"  |H| = {h} (gen {g}), {O} orbits: {'FOUND ' + str(perm) if W is not None else 'none'} [{time.time() - t:.1f}s]", flush=True)
            if W is not None:
                print('  check:', check(W), flush=True)
                import workdir as Wd
                with open(Wd.hits(f'qw_mult_{w}.txt'), 'w') as f:
                    for x in W: f.write(''.join('+' if v > 0 else '-' for v in x) + '\n')
                break
