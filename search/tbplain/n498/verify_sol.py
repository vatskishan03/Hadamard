#!/usr/bin/env python3
"""Reconstruct and exactly verify the solutions printed by band (SOL ida idb idc idd ga gb gc gd): the four sequences over G.

A SOL line gives the object ids of A, B, C, D (see plevel.py) and the Gray-code index of the sign vector of each M part
(sign i = -1 iff bit i of gray(index) is set, gray(m) = m ^ (m >> 1), variables in the order of mcomp.variables).  The sequences are
rebuilt from the slices and the signs and checked by the group convolution over G = Z2 x Z3 x Z41 (structure.check_solution_full:
A, B symmetric and sum_k PAF_k = target at all 246 elements), independently of the search code.  Writes the verified ones as
[{"ids", "gray", "seqs"}] to the --out file (the input of ../assemble.py, one record per SOL line, in the order of the file).

    python3 verify_sol.py OUT/b16.sol [more .sol files] --out verified_solutions.json
"""
import sys, json, numpy as np
from structure import *
import plevel as PL, mcomp as MC
NV = PL.NV
def obj_slices(k, oid):
    if k in (0, 1):
        cen = 1 if k == 0 else 0
        s = oid // NV; y = oid % NV
        u0, u1, u2 = PL.SYMIDX[cen][s], y, int(PL.RF[y])
        typ = 'S'
    else:
        c2 = oid % NV; c1 = (oid // NV) % NV; c0 = oid // (NV * NV)
        u0, u1, u2 = NV + c0, c1, c2
        typ = 'C'
    return typ, MC.pslices(typ, u0, u1, u2)
def reconstruct(ids, gidx):
    seqs = []; rows = []
    for k in range(4):
        typ, p = obj_slices(k, ids[k])
        vs = MC.variables(typ, p)
        g = gidx[k] ^ (gidx[k] >> 1)
        s = np.array([-1 if (g >> i) & 1 else 1 for i in range(len(vs))])
        m = MC.mslices(vs, s)
        Rw = MC.to_words(p, m)
        X = words_to_full(Rw)
        seqs.append(X); rows.append(Rw)
    return seqs, rows
def verify(seqs):
    return is_symmetric(seqs[0]) and is_symmetric(seqs[1]) and check_solution_full(seqs)
if __name__ == "__main__":
    args = sys.argv[1:]
    if "--out" not in args: sys.exit(__doc__)
    outfn = args[args.index("--out") + 1]; args = [x for x in args if x not in ("--out", outfn)]
    out = []
    for fn in args:
        for line in open(fn):
            if not line.startswith("SOL"): continue
            v = list(map(int, line.split()[1:]))
            seqs, rows = reconstruct(v[:4], v[4:])
            ok = verify(seqs)
            print("SOL", v, "exact check:", ok, flush=True)
            if ok:
                out.append(dict(ids=v[:4], gray=v[4:], seqs=["".join('+' if x == 1 else '-' for x in X) for X in seqs]))
    json.dump(out, open(outfn, "w"))
    print("verified", len(out))
