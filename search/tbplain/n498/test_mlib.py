#!/usr/bin/env python3
"""mlib.c (the M part of band.c) against numpy: sign lists, partial-agreement counts and the exact four-list join with a planted solution.

For 4 random draws of four objects (two symmetric A/B type, two C/D type, at most 12 free signs each):
 (1) the lists made by gen_list (Gray order) hold the same multiset of 13-column vectors as the numpy lists of mcomp.mlist;
 (2) count_partial (hash join of the first j columns, j = 2, 4, 7, against a planted target) equals a numpy count by sorting;
 (3) join4 on all 13 columns, with the target the sum of one planted member of each list, returns the planted quadruple.
Needs mlib.so (gcc -O2 -march=native -shared -fPIC -o mlib.so mlib.c -lm).

    python3 test_mlib.py
"""
import ctypes, os
import numpy as np
import plevel as PL, mcomp as MC

lib = ctypes.CDLL(os.path.join(os.path.dirname(os.path.abspath(__file__)), "mlib.so"))
lib.minit()
LANES = 16
class Block(ctypes.Structure):
    _fields_ = [("n", ctypes.c_int), ("nv", (ctypes.c_int * 2) * 32), ("np", ctypes.c_int * 32), ("diag", ctypes.c_int16 * LANES), ("K", ((ctypes.c_int8 * LANES) * 32) * 32)]
P16 = np.ctypeslib.ndpointer(dtype=np.int16, flags="C_CONTIGUOUS")
PI32 = np.ctypeslib.ndpointer(dtype=np.int32, flags="C_CONTIGUOUS")
lib.build_block.argtypes = [ctypes.c_int, ctypes.POINTER(ctypes.c_int), ctypes.POINTER(Block)]; lib.build_block.restype = ctypes.c_int
lib.gen_list.argtypes = [ctypes.POINTER(Block), ctypes.c_void_p]
lib.join4.argtypes = [P16, ctypes.c_long, P16, ctypes.c_long, P16, ctypes.c_long, P16, ctypes.c_long, P16, PI32, ctypes.c_long, ctypes.POINTER(ctypes.c_long)]; lib.join4.restype = ctypes.c_long
lib.count_partial.argtypes = [P16, ctypes.c_long, P16, ctypes.c_long, P16, ctypes.c_long, P16, ctypes.c_long, PI32, ctypes.c_int, P16]; lib.count_partial.restype = ctypes.c_long

def build(typ, p):
    """typ 0 = C/D type, 1 = symmetric A/B type; p: 3x9 int slices (P part)"""
    B = Block(); arr = (ctypes.c_int * 27)(*[int(x) for x in np.asarray(p).ravel()])
    assert lib.build_block(typ, arr, ctypes.byref(B)) >= 0
    return B
def gen(B):
    out = np.zeros(((1 << B.n), LANES), dtype=np.int16); lib.gen_list(ctypes.byref(B), out.ctypes.data); return out
def join(Vs, T):
    Vs = [np.ascontiguousarray(v, dtype=np.int16) for v in Vs]; out = np.zeros(4000, dtype=np.int32); nm = ctypes.c_long(0)
    ns = lib.join4(Vs[0], len(Vs[0]), Vs[1], len(Vs[1]), Vs[2], len(Vs[2]), Vs[3], len(Vs[3]), np.ascontiguousarray(T, dtype=np.int16), out, 1000, ctypes.byref(nm))
    return out[:4 * min(ns, 1000)].reshape(-1, 4), ns
def partial(Vs, cols, T):
    Vs = [np.ascontiguousarray(v, dtype=np.int16) for v in Vs]
    return lib.count_partial(Vs[0], len(Vs[0]), Vs[1], len(Vs[1]), Vs[2], len(Vs[2]), Vs[3], len(Vs[3]), np.ascontiguousarray(cols, dtype=np.int32), len(cols), np.ascontiguousarray(T, dtype=np.int16))

rng = np.random.default_rng(11)
NV = PL.NV
def rand_obj(typ):
    if typ == 0: return (0, PL.NV + rng.integers(PL.NV), rng.integers(PL.NV), rng.integers(PL.NV))
    s0 = PL.SYMIDX[1][rng.integers(19)]; y = rng.integers(PL.NV)
    return (1, s0, y, PL.RF[y])
for trial in range(4):
    blocks = []
    for typ in (1, 1, 0, 0):
        while True:
            t, u0, u1, u2 = rand_obj(typ)
            p = MC.pslices('C' if typ == 0 else 'S', u0, u1, u2)
            B = build(typ, p)
            if B.n <= 12: break
        blocks.append((typ, p, B))
    Vc = []
    for typ, p, B in blocks:
        Vcl = gen(B)[:, :13]
        vs, Q, S, Vp = MC.mlist('C' if typ == 0 else 'S', p)
        assert B.n == len(vs)
        a = np.sort(Vcl.astype(np.int64) @ np.arange(1, 14) * 1000003 + (Vcl.astype(np.int64) ** 2).sum(1))
        b = np.sort(Vp @ np.arange(1, 14) * 1000003 + (Vp ** 2).sum(1))
        assert (a == b).all(), "list mismatch"
        Vc.append(gen(B))
    Vp = [MC.mlist('C' if typ == 0 else 'S', p)[3] for typ, p, B in blocks]
    pl = [int(rng.integers(len(v))) for v in Vc]
    Tn = sum(Vc[k][pl[k]].astype(np.int64) for k in range(4)).astype(np.int16)
    for j in (2, 4, 7):
        cols = np.arange(j)
        c_c = partial(Vc, cols, Tn)
        wr = np.random.default_rng(5).integers(1, 2 ** 40, size=j).astype(np.int64) * 2 + 1
        hA, hB, hC, hD = [(v[:, :j].astype(np.int64) @ wr) for v in Vp]
        kac = (hA[:, None] + hC[None, :]).ravel(); hT = int(Tn[:j].astype(np.int64) @ wr)
        kbd = (hT - (hB[:, None] + hD[None, :])).ravel()
        ks = np.sort(kac)
        c_n = int((np.searchsorted(ks, kbd, 'right') - np.searchsorted(ks, kbd, 'left')).sum())
        print("n", [b[2].n for b in blocks], "j", j, "C count", c_c, "numpy", c_n)
        assert c_c == c_n
    sols, ns = join(Vc, Tn)
    found = any(tuple(x) == tuple(pl) for x in sols)
    print("full join solutions", ns, "planted", pl, "found planted", found)
    assert found
print("mlib test OK")
