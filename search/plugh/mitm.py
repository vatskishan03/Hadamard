import ctypes, numpy as np, os
_lib = ctypes.CDLL(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'libcmitm.so'))
_lib.mitm4.restype = ctypes.c_long
_V = ctypes.c_void_p; _I = ctypes.c_int; _L = ctypes.c_long
_lib.mitm4.argtypes = [_I, _I, _I, _I, _I, _V, _V, _V, _V, _V, _I, _V, _L, ctypes.c_uint64, _I]
def mitm4(VA, VB, VC, VD, target, logR=0, maxout=50_000_000, seed=1, verbose=0):
    """all (a,b,c,d) with VA[a]+VB[b]+VC[c]+VD[d]==target (int16 rows). returns (n_found_total, array (n_stored,4) uint32)"""
    Vs = [np.ascontiguousarray(v, dtype=np.int16) for v in (VA, VB, VC, VD)]
    d = Vs[0].shape[1]
    assert all(v.shape[1] == d for v in Vs) and d <= 64
    tgt = np.ascontiguousarray(target, dtype=np.int16)
    out = np.zeros((min(maxout, 1 << 20), 4), dtype=np.uint32) if maxout > (1 << 20) else np.zeros((maxout, 4), dtype=np.uint32)
    cap = len(out)
    while True:
        n = _lib.mitm4(len(Vs[0]), len(Vs[1]), len(Vs[2]), len(Vs[3]), d,
                       *[v.ctypes.data_as(ctypes.c_void_p) for v in Vs], tgt.ctypes.data_as(ctypes.c_void_p),
                       logR, out.ctypes.data_as(ctypes.c_void_p), cap, seed, verbose)
        if n < 0: raise MemoryError('mitm4 table allocation failed')
        if n <= cap: return n, out[:n].copy()
        if cap >= maxout: return n, out[:cap].copy()
        cap = min(maxout, max(cap * 4, n + 10)); out = np.zeros((cap, 4), dtype=np.uint32)
