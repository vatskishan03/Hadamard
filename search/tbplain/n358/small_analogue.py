#!/usr/bin/env python3
"""Positive end-to-end control on the small analogue G = Z2^4 x Z3 (N = 51, order 102): the same code path, exhaustively.

structure "i3" (same H and phi as QH(358), Z11 replaced by Z3, 11 PAF columns instead of 35) -> Q-level tuples (qlevel) ->
all sequences of each kind with the Q-compression of the tuple (pools of up to 2^18 members, listed completely) -> exact
four-list join on all 11 columns with pjoin -> every full match checked by the definition of the PAF over G -> border and
assembly (../assemble.py, structure "102") -> a quaternary Hadamard matrix of order 102 with K K* = 102 I exactly.
The patterns are taken in order until three cores are found; the pair of pools of a pattern must be small enough for the join.

    python3 small_analogue.py --out DIR            writes DIR/complex-hadamard/C102.json
Needs ./pjoin.  About 15 seconds.
"""
import sys, os, json, time, argparse, tempfile
import numpy as np
from structure import Struct
from qlevel import tuples, rs, cell_q
from pio import write_lists, run_pjoin, read_out
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from assemble import assemble


def all_pool(S, kind, y):
    """every sign vector on the cells of a kind with Q-compression y = (y0, y1), by enumeration of all 2^ncells patterns"""
    n = int((S.orbA if kind == "A" else S.orbC).max()) + 1
    w0, w1 = cell_q(S, kind)
    out = []
    N = 1 << n
    for s0 in range(0, N, 1 << 18):
        idx = np.arange(s0, min(N, s0 + (1 << 18)), dtype=np.int64)
        bits = ((idx[:, None] >> np.arange(n)[None, :]) & 1) * 2 - 1
        ok = (bits @ w0 == y[0]) & (bits @ w1 == y[1])
        if ok.any(): out.append(bits[ok].astype(np.int8))
    return np.concatenate(out) if out else np.zeros((0, n), np.int8)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    S = Struct("i3")
    sols, dA, dC, tg = tuples(S)
    seen = {}
    for t in sols:
        t2 = tuple(y if rs(y) >= 0 else (-y[0], -y[1]) for y in t)
        seen.setdefault(tuple(abs(rs(y)) for y in t), t2)
    print("i3 m", S.m, "N", S.N, "cells A/C", S.nA, S.nC, "columns", len(S.cols), "ordered patterns", len(seen), flush=True)
    tmp = tempfile.mkdtemp(prefix="small_"); fn, out = os.path.join(tmp, "small.bin"), os.path.join(tmp, "small.out")
    found = []; ncores = 0
    for key in sorted(seen):
        tup = seen[key]
        pools = [all_pool(S, k, y) for k, y in zip("AACC", tup)]
        sizes = [len(p) for p in pools]
        print(" pattern", key, "pool sizes", sizes, flush=True)
        if min(sizes) == 0 or sizes[2] * sizes[0] > 6e10: continue
        V = [S.paf_classes(S.expand(p, k))[:, S.cols] for p, k in zip(pools, "AACC")]
        write_lists(fn, V, S.T)
        run_pjoin(fn, len(S.cols), 6, 3, 0, 1, len(S.cols), out)
        h, q, st = read_out(out, len(S.cols), len(S.cols))
        print("   full matches:", h.get(len(S.cols), 0), st, flush=True)
        ncores += h.get(len(S.cols), 0)
        for (l, i, j, k, l4) in q:
            seqs = [S.expand(pools[0][i:i + 1], "A")[0], S.expand(pools[1][j:j + 1], "A")[0], S.expand(pools[2][k:k + 1], "C")[0], S.expand(pools[3][l4:l4 + 1], "C")[0]]
            tot = sum(S.paf_direct(x) for x in seqs)
            assert (tot == S.target).all(), "join solution fails the exact PAF check"
            found.append(["".join("+" if v == 1 else "-" for v in x) for x in seqs])
        if len(found) >= 3: break
    for f in (fn, out): os.path.exists(f) and os.remove(f)
    os.rmdir(tmp)
    print("cores found: %d (all verified by the definition of the PAF over G)" % ncores)
    assert found, "no core found"
    cert = assemble("102", found[0], a.out)
    print("assembled: order %d, K K* = %d I exact; sha256 %s" % (cert["order"], cert["order"], cert["sha256"]))
    print("SMALL ANALOGUE PASSED")
