# independent corroboration of the large-|K| negative (n = 1981 shape: gs, K = Z283, r = 7, H = <4>): pure-numpy level 2 (full I and Q equations over all t signs, 4-list join by dictionary, no row-sum or PSD filter)
# on sampled level-1 instances versus gl2.c (they must agree; both are expected to give 0 completions there)
import sys, time, itertools, random; sys.path.insert(0, '.')
import numpy as np
from gsk import *; from gl1 import *; from gl2io import run_gl2
sh = GShape([283], 7, [([4], 1)], 'gs'); l1 = L1(sh); K = sh.kCD
ds = set(sh.feasible_totals()); W = K.weights()
lists = {n: enum_weight(sh, n) for n in W}; lists = {n: v for n, v in lists.items() if v is not None}
allocs = [al for al in itertools.combinations_with_replacement(sorted(lists), 4) if sum(al) in ds]
random.seed(5); random.shuffle(allocs)
insts = []; t0 = time.time()
for al in allocs:
    order = sorted(range(4), key=lambda i: len(lists[al[i]][0]))
    specs = [BlockSpec(al[i], S=lists[al[i]][0], E=lists[al[i]][1]) for i in order]
    if max(float(specs[0].L) * specs[2].L, float(specs[1].L) * specs[3].L) > 5e8: continue
    cnt, S, E, st = l1.run(specs, maxout=40)
    for i in range(min(cnt, 40)): insts.append([(S[i, b], E[i, b]) for b in range(4)])
    if len(insts) >= 600 or time.time() - t0 > 100: break
print('instances', len(insts), flush=True)
def py_level2(blocks):
    """count t quadruples (t = tau sigma on free orbits) with I(k) = 0 (k != 0), Q(k) = 0 for the reps, plus I(0) = T2 (sum of squares identity)"""
    r = sh.r; eps = sh.eps; N = sh.N
    feats = []; sizes = []
    for (s, e) in blocks:
        s64 = s.astype(np.int64)
        free = [i for i in range(K.nu) if e[int(K.pos[i][0])] == 0]
        M = 1 << len(free); ids = np.arange(M)
        t = np.zeros((M, N), np.int64)
        for j, i in enumerate(free):
            P = K.pos[i]; bit = 1 - 2 * ((ids >> j) & 1); t[:, P] = bit[:, None] * K.sig[P][None, :]
        cols = []
        for k in sh.repsPM:
            if k == 0: continue
            sk = sh.shift[k]
            cols.append((s64 * s64[sk]).sum() + r * (t * t[:, sk]).sum(1))
        for k in sh.repsQ:
            sk = sh.shift[k]
            cols.append((s64[None, :] * t[:, sk]).sum(1) + eps * (t * s64[sk][None, :]).sum(1))
        # I(0): sum s^2 + r sum t^2
        cols.append((s64 * s64).sum() + r * (t * t).sum(1))
        feats.append(np.stack(cols, 1))
    tgt = np.zeros(feats[0].shape[1], np.int64); tgt[-1] = sh.T2
    d = {}
    for i, f1 in enumerate(feats[0]):
        for j, f2 in enumerate(feats[1]): d.setdefault((f1 + f2).tobytes(), 0); d[(f1 + f2).tobytes()] += 1
    tot = 0
    for i, f3 in enumerate(feats[2]):
        for j, f4 in enumerate(feats[3]): tot += d.get((tgt - f3 - f4).tobytes(), 0)
    return tot
res, sols = run_gl2(sh, insts, maxsol=3)
bad = 0; nz = 0; t0 = time.time(); ncmp = 0
for i, bl in enumerate(insts[:200]):
    pc = py_level2(bl); cc = res[i]['compl']; ncmp += 1
    if pc != cc: bad += 1; print('MISMATCH', i, pc, cc, res[i]['raw'])
    nz += (pc > 0 or cc > 0)
print('compared', ncmp, 'mismatches', bad, 'nonzero', nz, 'gl2 status', {st: sum(1 for d in res.values() if d['status'] == st) for st in set(d['status'] for d in res.values())}, f'{time.time()-t0:.0f}s')
