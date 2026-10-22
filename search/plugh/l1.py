import ctypes, numpy as np, os, math
HERE = os.path.dirname(os.path.abspath(__file__))
_lib = ctypes.CDLL(os.path.join(HERE, 'libl1.so'))
_V = ctypes.c_void_p; _I = ctypes.c_int; _L = ctypes.c_long
_lib.l1_run.restype = _L
_lib.l1_run.argtypes = [_I, _I, _I, _V, _V, _I, _V, _V, _V, _V, _V, _I, _V, _I, _I, _I, _L, _L, ctypes.c_uint64, _I, _V, _V, _I, _V, _V, _I, _L, _V, _V, _V, _V, _V]
_lib.l1_gen.restype = None
_lib.l1_gen.argtypes = [_I, _I, _I, _V, _V, _V, _V, _V, _I, _V, _I, _I, ctypes.c_uint64, ctypes.c_uint64, _V, _V]
from est import kind_dpu


class KindC:
    """C-side description of the C/D kind of a plain Shape"""
    def __init__(self, sh):
        K = sh.kCD; self.K = K; self.sh = sh
        dp, ff = kind_dpu(K)
        self.dpw = np.ascontiguousarray(np.array(dp, dtype=np.float64))     # (nfree+1) x (N+1)
        self.nfree = len(K.fi); self.fw = int(K.fw)
        self.freeidx = np.ascontiguousarray(np.array(K.fi, dtype=np.int32)) if self.nfree else np.zeros(1, np.int32)
        self.w = np.ascontiguousarray(K.w.astype(np.int32)); self.forced = np.ascontiguousarray(K.forced.astype(np.int32))
        self.poff = np.zeros(K.nu + 1, np.int32); self.plist = []
        for i, P in enumerate(K.pos): self.poff[i + 1] = self.poff[i] + len(P); self.plist += [int(x) for x in P]
        self.plist = np.ascontiguousarray(np.array(self.plist, dtype=np.int32))
        self.reps = np.ascontiguousarray(np.array(sh.repsPM, dtype=np.int32))
        self.target = np.ascontiguousarray(sh.T[:sh.nU].astype(np.int16))
        self.n = sh.N; self.q = sh.q; self.nu = K.nu

    def gen(self, nw, seed, j):
        s = np.zeros(self.n, np.int8); e = np.zeros(self.n, np.int8)
        _lib.l1_gen(self.n, self.q, self.nu, self.w.ctypes.data, self.forced.ctypes.data, self.poff.ctypes.data, self.plist.ctypes.data, self.dpw.ctypes.data,
                    self.nfree, self.freeidx.ctypes.data, self.fw, nw, seed, j, s.ctypes.data, e.ctypes.data)
        return s, e

    def run(self, nC, nD, LC, seed, SA, EA, SB, EB, per_class=30000.0, maxout=1500000, LD=None):
        if LD is None: LD = LC
        SA = np.ascontiguousarray(SA, np.int8); EA = np.ascontiguousarray(EA, np.int8); SB = np.ascontiguousarray(SB, np.int8); EB = np.ascontiguousarray(EB, np.int8)
        pairs = max(len(SA) * LC, len(SB) * LD)
        lgM = max(1, int(math.ceil(math.log2(max(2.0, pairs / per_class)))))
        lgM = min(lgM, 28)
        n = self.n
        oSC = np.zeros((maxout, n), np.int8); oEC = np.zeros((maxout, n), np.int8); oSD = np.zeros((maxout, n), np.int8); oED = np.zeros((maxout, n), np.int8)
        oAB = np.zeros((maxout, 2), np.uint32)
        cnt = _lib.l1_run(n, self.q, len(self.reps), self.reps.ctypes.data, self.target.ctypes.data, self.nu, self.w.ctypes.data, self.forced.ctypes.data,
                          self.poff.ctypes.data, self.plist.ctypes.data, self.dpw.ctypes.data, self.nfree, self.freeidx.ctypes.data, self.fw, nC, nD, LC, LD, seed,
                          len(SA), SA.ctypes.data, EA.ctypes.data, len(SB), SB.ctypes.data, EB.ctypes.data, lgM, maxout,
                          oSC.ctypes.data, oEC.ctypes.data, oSD.ctypes.data, oED.ctypes.data, oAB.ctypes.data)
        m = min(cnt, maxout)
        return int(cnt), oSC[:m], oEC[:m], oSD[:m], oED[:m], oAB[:m]
