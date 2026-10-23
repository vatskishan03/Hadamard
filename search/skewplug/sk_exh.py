"""Exhaustive two-level pipeline for the skew plug-in (explicit lists, exact join): returns the exact ordered number of K-level solutions (= G-level quadruples) of a small shape.
E = 8 * sum over (nA ; ordered (nB,nC,nD)) of  level1(cansign B,C,D) completions (x2 t-gauge inside sl2)."""
import numpy as np, itertools, math, time, sys, json
from collections import Counter
from sk import *
from gl1 import L1, BlockSpec
from sl2io import run_sl2, header_text


def perm3(al):
    c = Counter(al); m = 6
    for v in c.values(): m //= math.factorial(v)
    return m


def l1_all(l1, specs, maxout=200000, per_class=30000.0):
    cnt, S, E, st = l1.run(specs, per_class=per_class, maxout=maxout)
    if cnt <= maxout:
        yield cnt, S, E; return
    L0 = specs[0].L
    if L0 <= 1: raise RuntimeError('single outer element yields too many instances')
    h = L0 // 2
    for lo, hi in ((0, h), (h, L0)):
        sp = [BlockSpec(specs[0].nw, S=specs[0].S[lo:hi], E=specs[0].E[lo:hi], cansign=False)] + specs[1:]
        yield from l1_all(l1, sp, maxout, per_class)


def exhaust(sh, tmax=900.0, want_sols=8, chunk=2000, maxfib=2000000, tl2=60.0, verbose=False, max_inst=4_000_000, pairs_lo=0.0, pairs_hi=3e9):
    t0 = time.time(); l1 = L1(sh)
    ds = set(sh.feasible_totals()); W = sh.kCD.weights()
    listsF = {}; toolargeF = set(); toolargeA = set()
    for n in W:
        try:
            r = enum_weight(sh, n, cansign=True, cap=8_000_000)
        except RuntimeError:
            toolargeF.add(n); continue
        if r is not None: listsF[n] = r
    listsA = {}
    for n in sh.A_weights():
        try:
            r = sh.enum_A(n, cap=8_000_000)
        except RuntimeError:
            toolargeA.add(n); continue
        if r is not None: listsA[n] = r
    hdr = header_text(sh)
    tot = dict(allocs=0, n1=0, compl=0, skipped=[], flags=[])
    sols = []
    for nA in sorted(set(listsA) | toolargeA):
        for alF in itertools.combinations_with_replacement(sorted(set(listsF) | toolargeF), 3):
            if nA + sum(alF) not in ds: continue
            tot['allocs'] += 1
            if nA in toolargeA or any(x in toolargeF for x in alF): tot['skipped'].append(((nA,) + alF, 'list too large')); continue
            if time.time() - t0 > tmax: tot['skipped'].append(((nA,) + alF, 'timeout')); continue
            al = (nA,) + alF
            lists = [listsA[nA]] + [listsF[n] for n in alF]
            order = sorted(range(4), key=lambda i: len(lists[i][0]))
            specs = [BlockSpec(al[i], S=lists[i][0], E=lists[i][1], cansign=False) for i in order]
            pairs = max(float(specs[0].L) * specs[2].L, float(specs[1].L) * specs[3].L)
            if pairs > pairs_hi: tot['skipped'].append((al, f'join too big {pairs:.2g}')); continue
            if pairs <= pairs_lo: continue          # (already done in an earlier run with pairs_hi = pairs_lo)
            cnt0, _, _, _ = l1.run(specs, maxout=1)
            if cnt0 > max_inst: tot['skipped'].append((al, f'{cnt0} level-1 instances')); continue
            inv = [order.index(b) for b in range(4)]
            mult = perm3(alF); ncomp = 0; flagged = 0
            for cnt, S, E in l1_all(l1, specs):
                for c0 in range(0, len(S), chunk):
                    insts = [[(S[i, inv[b]], E[i, inv[b]]) for b in range(4)] for i in range(c0, min(len(S), c0 + chunk))]
                    res, ss = run_sl2(sh, insts, maxsol=want_sols, tlimit=tl2, maxfib=maxfib, hdr=hdr)
                    for i, d in res.items():
                        ncomp += d['compl']
                        if d['status'] in ('TOOBIG', 'TIMEOUT'): flagged += 1
                    if len(sols) < want_sols:
                        for (i, ts) in ss[:want_sols - len(sols)]:
                            sols.append([(insts[i][X][0].astype(np.int64), ts[X].astype(np.int64), insts[i][X][1].astype(np.int64)) for X in range(4)])
            tot['n1'] += cnt0 * mult; tot['compl'] += ncomp * mult
            if flagged: tot['flags'].append((al, flagged))
            if verbose: print(f'  alloc {al} mult {mult}: level-1 {cnt0}, completions {ncomp} [{time.time()-t0:.0f}s]', flush=True)
    tot['E_G'] = 8 * tot['compl']; tot['N1_G'] = 8 * tot['n1']
    tot['exhaustive'] = (not tot['skipped']) and (not tot['flags'])
    tot['secs'] = time.time() - t0
    for s in sols: assert sh.check_solution(s), 'K-level check failed'
    return tot, sols
