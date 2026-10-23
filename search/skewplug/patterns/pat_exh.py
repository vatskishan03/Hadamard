"""Exhaustive level 1 + level 2 restricted to ONE defect-position pattern (A: pairs with types; B, C, D: orbits with types): explicit lists of all sign
configurations, exact hash join (../gl1.c), exact completion of the t signs (../gl2.c, 'sl2').  solve_pattern() is EXACT for the pattern: it returns the number
of level-1 instances (canonical block signs) and the number of level-2 completions, and the solutions found.
(Library module; used by run_pat.py, val_pat.py, val_n1.py, val_sym.py, posctl_1100.py.)"""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, '..')); sys.path.insert(0, HERE)
import numpy as np, itertools, time, math
from sk import SShape
from gl1 import L1, BlockSpec
from sl2io import run_sl2, header_text
from sk_exh import l1_all


def free_list(sh, items, cansign=True):
    """items: ((orbit, type),...) of the block; all sign configurations over the H_K-orbits; returns S, E int8 (M, N)"""
    K = sh.kCD; nu = K.nu; N = sh.N
    dt = {o: ty for (o, ty) in items}
    opts = []
    for o in range(nu):
        if o in dt:
            opts.append([(0, 1), (0, -1)] if dt[o] == 0 else [(-2, 1), (2, -1)])
        else: opts.append([(1, 0), (-1, 0)])
    M = 1 << nu
    ids = np.arange(M, dtype=np.int64)
    S = np.zeros((M, N), np.int8); E = np.zeros((M, N), np.int8)
    for o in range(nu):
        bit = ((ids >> o) & 1)
        a = np.array(opts[o], dtype=np.int8)[bit]
        P = K.pos[o]; S[:, P] = a[:, 0:1]; E[:, P] = a[:, 1:2]
    if cansign:
        o0 = 0
        P0 = int(K.pos[o0][0]); v = np.where(E[:, P0] != 0, E[:, P0], S[:, P0])
        keep = v > 0
        S, E = S[keep], E[keep]
    return S, E


def A_list(sh, items):
    """items: ((pair index, type),...); all skew-type configs with exactly these defect pairs"""
    N = sh.N; pairs = sh.apairs; pos = sh.spos; npair = len(pairs)
    dt = {j: ty for (j, ty) in items}
    opts = []
    for j in range(npair):
        if j in dt: opts.append([(0, 1), (0, -1)] if dt[j] == 0 else [(-2, 1), (2, -1)])
        else: opts.append([(1, 0), (-1, 0)])
    M = 1 << npair; ids = np.arange(M, dtype=np.int64)
    S = np.zeros((M, N), np.int8); E = np.zeros((M, N), np.int8); S[:, 0] = 1
    for j, (o, po) in enumerate(pairs):
        a = np.array(opts[j], dtype=np.int8)[(ids >> j) & 1]
        S[:, pos[o]] = a[:, 0:1]; E[:, pos[o]] = a[:, 1:2]
        S[:, pos[po]] = -a[:, 0:1]; E[:, pos[po]] = -a[:, 1:2]
    return S, E


def weight_of(sh, items, A=False):
    if A: return sum(2 * len(sh.spos[sh.apairs[j][0]]) for (j, ty) in items)
    return sum(len(sh.spos[o]) for (o, ty) in items)


def solve_pattern(sh, l1, hdr, A, B, C, D, want_sols=2, chunk=1000, tl2=120.0, maxfib=2000000, max_inst=3_000_000, per_class=30000.0, verbose=False):
    t0 = time.time()
    lists = [A_list(sh, A)] + [free_list(sh, b) for b in (B, C, D)]
    wts = [weight_of(sh, A, True)] + [weight_of(sh, b) for b in (B, C, D)]
    order = sorted(range(4), key=lambda i: len(lists[i][0]))
    specs = [BlockSpec(wts[i], S=lists[i][0], E=lists[i][1], cansign=False) for i in order]
    cnt0, _, _, _ = l1.run(specs, maxout=1, per_class=per_class)
    rec = dict(n1=int(cnt0), lens=[len(x[0]) for x in lists], compl=0, flagged=0, found=0, status='ok')
    if cnt0 > max_inst: rec['status'] = f'too many level-1 instances ({cnt0})'; return rec, []
    inv = [order.index(b) for b in range(4)]
    sols = []
    for cnt, S, E in l1_all(l1, specs, per_class=per_class):
        for c0 in range(0, len(S), chunk):
            insts = [[(S[i, inv[b]], E[i, inv[b]]) for b in range(4)] for i in range(c0, min(len(S), c0 + chunk))]
            res, ss = run_sl2(sh, insts, maxsol=want_sols, tlimit=tl2, maxfib=maxfib, hdr=hdr)
            for i, d in res.items():
                rec['compl'] += d['compl']
                if d['status'] in ('TOOBIG', 'TIMEOUT'): rec['flagged'] += 1
            for (i, ts) in ss[:max(0, want_sols - len(sols))]:
                sols.append([(insts[i][X][0].astype(np.int64), ts[X].astype(np.int64), insts[i][X][1].astype(np.int64)) for X in range(4)])
    rec['found'] = len(sols); rec['secs'] = round(time.time() - t0, 1)
    if rec['flagged']: rec['status'] = 'FLAGGED(level 2 timeout/toobig)'
    return rec, sols
