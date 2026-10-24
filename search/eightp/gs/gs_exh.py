"""Exhaustive two-level solver for the real GS-type plug-ins (four free blocks): level 1 (gl1.c, explicit lists, exact join) + level 2 (gl2.c).
Counts ORDERED K-level solutions (blocks in a fixed order, canonical block signs, t-gauge included) and returns E_G = 16 * (that) = the number of ordered
quadruples of +-1 sequences on G invariant under the plug-in group (comparable with gbf.count_G).
"""
import numpy as np, itertools, math, time, sys, json
from collections import Counter
from gsk import *
from gl1 import *
from gl2io import *


def perm_mult(al):
    m = math.factorial(4)
    for v in Counter(al).values(): m //= math.factorial(v)
    return m


def l1_chunks(l1, specs, maxout, per_class=30000.0):
    """level-1 join of explicit lists, recursively splitting block 0 until each call returns <= maxout instances; yields (S, E)"""
    cnt, S, E, st = l1.run(specs, per_class=per_class, maxout=maxout)
    if cnt <= maxout:
        yield cnt, S, E
        return
    L0 = specs[0].L
    if L0 <= 1: raise RuntimeError('single outer element yields too many instances')
    h = L0 // 2
    for lo, hi in ((0, h), (h, L0)):
        sp = [BlockSpec(specs[0].nw, S=specs[0].S[lo:hi], E=specs[0].E[lo:hi], cansign=specs[0].cansign)] + specs[1:]
        yield from l1_chunks(l1, sp, maxout, per_class)


def exhaust(sh, max_inst=3_000_000, tmax=600.0, want_sols=10, chunk=2500, verbose=False, maxfib=2000000, tl2=60.0):
    t0 = time.time(); l1 = L1(sh)
    ds = set(sh.feasible_totals()); W = sh.kCD.weights()
    lists = {}
    for n in W:
        r = enum_weight(sh, n)
        if r is not None: lists[n] = r
    allocs = [al for al in itertools.combinations_with_replacement(sorted(lists), 4) if sum(al) in ds]
    tot = dict(allocs=len(allocs), n1=0, compl=0, nonzero_inst=0, skipped=[], flags=[])
    sols = []
    for al in allocs:
        if time.time() - t0 > tmax: tot['skipped'].append((al, 'timeout')); continue
        order = sorted(range(4), key=lambda i: len(lists[al[i]][0]))
        specs = [BlockSpec(al[i], S=lists[al[i]][0], E=lists[al[i]][1]) for i in order]
        pairs = max(float(specs[0].L) * specs[2].L, float(specs[1].L) * specs[3].L)
        if pairs > 3e9: tot['skipped'].append((al, f'join too big {pairs:.2g}')); continue
        cnt0, _, _, _ = l1.run(specs, maxout=1)
        if cnt0 > max_inst: tot['skipped'].append((al, f'{cnt0} level-1 instances')); continue
        mult = perm_mult(al); ncomp = 0; nn = 0; flagged = 0
        for cnt, S, E in l1_chunks(l1, specs, 200000):
            for c0 in range(0, len(S), chunk):
                insts = [[(S[i, b], E[i, b]) for b in range(4)] for i in range(c0, min(len(S), c0 + chunk))]
                res, ss = run_gl2(sh, insts, maxsol=want_sols, tlimit=tl2, maxfib=maxfib)
                for i, d in res.items():
                    ncomp += d['compl']; nn += (d['compl'] > 0)
                    if d['status'] in ('TOOBIG', 'TIMEOUT'): flagged += 1
                if len(sols) < want_sols:
                    for (i, ts) in ss[:want_sols - len(sols)]:
                        sols.append([(insts[i][X][0].astype(np.int64), ts[X].astype(np.int64), insts[i][X][1].astype(np.int64)) for X in range(4)])
        tot['n1'] += cnt0 * mult; tot['compl'] += ncomp * mult; tot['nonzero_inst'] += nn * mult
        if flagged: tot['flags'].append((al, flagged))
        if verbose: print(f'  alloc {al} mult {mult}: level-1 {cnt0}, completions {ncomp}  [{time.time()-t0:.0f}s]', flush=True)
    tot['E_G'] = 16 * tot['compl']; tot['N1_G'] = 16 * tot['n1']
    tot['exhaustive'] = (not tot['skipped']) and (not tot['flags'])
    tot['secs'] = time.time() - t0
    for s in sols: assert sh.check_klevel(s)
    return tot, sols


if __name__ == '__main__':
    Kmods = json.loads(sys.argv[1]); r = int(sys.argv[2]); mode = sys.argv[3]
    hg = [(tuple(u) if isinstance(u, list) else u, ph) for u, ph in json.loads(sys.argv[4])] if len(sys.argv) > 4 else []
    T = json.loads(sys.argv[5]) if len(sys.argv) > 5 else None
    sh = GShape(Kmods, r, hg, mode, T)
    tot, sols = exhaust(sh, verbose=True)
    print('RESULT', {k: v for k, v in tot.items()})
