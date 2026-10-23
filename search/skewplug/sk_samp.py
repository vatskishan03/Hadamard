"""sampler / probe for the skew plug-in at an allocation: level 1 (gl1.c: exact A list x sampled B, C, D lists, class-restricted join) + level 2 (sl2)."""
import numpy as np, time, math, json, sys, os, random
from sk import SShape
from sk_pat import allocations, m_splits, row_patterns_skew
from gl1 import L1, BlockSpec
from sl2io import run_sl2, header_text


def level1(sh, l1, al, ms, pat, LB, LC, LD, seed, A_list=None, maxout=200000, per_class=30000.0, hseed=1, cache=None):
    """pat = ((sB,eB),(sC,eC),(sD,eD)) or None.  returns (count, S, E, stats); blocks in order A,B,C,D"""
    nA, nB, nC, nD = al
    if A_list is None:
        A_list = sh.enum_A(nA, ms[0])
    if A_list is None: return 0, None, None, None
    SA, EA = A_list
    specs = [BlockSpec(nA, S=SA, E=EA, cansign=False)]
    for i, (L, nn) in enumerate(((LB, nB), (LC, nC), (LD, nD))):
        fs = fe = None
        if pat is not None: fs, fe = pat[i]
        specs.append(BlockSpec(nn, L=L, seed=seed * 7919 + 101 * (i + 1), mw=ms[i + 1], cansign=False, fs=fs, fe=fe))
    return l1.run(specs, per_class=per_class, maxout=maxout, hseed=hseed)


def level2(sh, S, E, hdr, chunk=500, stop_first=True, tl2=30.0, maxfib=1500000, want=5):
    res_all = {}; sols = []; ncomp = 0; nok = 0; tl = 0.0; stat = {}
    for c0 in range(0, len(S), chunk):
        insts = [[(S[i, b], E[i, b]) for b in range(4)] for i in range(c0, min(len(S), c0 + chunk))]
        res, ss = run_sl2(sh, insts, maxsol=want, tlimit=tl2, maxfib=maxfib, stop_first=stop_first, hdr=hdr)
        for i, d in res.items():
            stat[d['status']] = stat.get(d['status'], 0) + 1
            ncomp += d['compl']; nok += d['compl'] > 0; tl += d.get('t', 0.0)
        for (i, ts) in ss[:want - len(sols)]:
            sols.append([(insts[i][X][0].astype(np.int64), ts[X].astype(np.int64), insts[i][X][1].astype(np.int64)) for X in range(4)])
    return dict(ncomp=ncomp, nok=nok, t2=tl, stat=stat), sols
