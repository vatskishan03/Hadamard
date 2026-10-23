"""python wrapper of gl1.c (level 1 of the real GS-type plug-ins, four free blocks)"""
import ctypes, numpy as np, os, math
from est import kind_dpu
HERE = os.path.dirname(os.path.abspath(__file__))
_lib = ctypes.CDLL(os.path.join(HERE, 'libgl1.so'))
_V = ctypes.c_void_p; _I = ctypes.c_int; _L = ctypes.c_long; _U = ctypes.c_uint64
_lib.gl1_run.restype = _L
_lib.gl1_run.argtypes = [_I, _I, _I, _V, _V, _I, _V, _V, _V, _V, _V, _I, _V, _I, _V, _V, _V, _V, _V, _V, _V, _I, _V, _V, _V, _I, _I, _V, _V, _I, _L, _V, _V, _V, _U, _L]
_lib.gl1_gen.restype = None
_lib.gl1_gen.argtypes = [_I, _I, _I, _V, _V, _V, _V, _V, _I, _V, _I, _I, _I, _U, _U, _V, _I, _I, _V, _V, _I, _I, _I, _I, _V, _V]


class BlockSpec:
    """one block's list: generated (nw, L, seed) or explicit (S, E arrays (m, N) int8); rok: set of allowed row sums or None; cansign: canonical block sign"""
    def __init__(self, nw, L=None, seed=0, S=None, E=None, rok=None, cansign=True, mw=-1, sok=None, eok=None, fs=None, fe=None):
        self.fs = fs; self.fe = fe; self.sok = sok; self.eok = eok; self.mw = mw; self.nw = nw; self.seed = seed; self.S = S; self.E = E; self.rok = rok; self.cansign = cansign
        self.L = len(S) if S is not None else int(L)


class L1:
    def __init__(self, sh):
        self.sh = sh; K = sh.kCD; self.K = K
        dp, ff = kind_dpu(K)
        self.dpw = np.ascontiguousarray(np.array(dp, dtype=np.float64))
        self.nfree = len(K.fi); self.fw = int(K.fw)
        self.freeidx = np.ascontiguousarray(np.array(K.fi, dtype=np.int32)) if self.nfree else np.zeros(1, np.int32)
        self.w = np.ascontiguousarray(K.w.astype(np.int32)); self.forced = np.ascontiguousarray(K.forced.astype(np.int32))
        self.poff = np.zeros(K.nu + 1, np.int32); plist = []
        for i, P in enumerate(K.pos): self.poff[i + 1] = self.poff[i] + len(P); plist += [int(x) for x in P]
        self.plist = np.ascontiguousarray(np.array(plist, dtype=np.int32))
        self.n = sh.N; self.r = sh.r; self.nu = K.nu
        self.shift = np.ascontiguousarray(np.stack([sh.shift[k] for k in sh.repsPM]).astype(np.int32))
        self.target = np.ascontiguousarray(sh.Utarget.astype(np.int32))
        self.RM = sh.N * sh.r + 2
        assert sh.N <= 640
        assert len(sh.repsPM) <= 96, 'more than 96 U features (H_K-orbits of +-k): not supported by gl1.c'

    def rok_table(self, rok):
        if rok is None: return None
        t = np.zeros(2 * self.RM + 1, np.uint8)
        for R in rok:
            if abs(R) <= self.RM: t[R + self.RM] = 1
        return t

    def sok_table(self, ok, M):
        if ok is None: return None
        t = np.zeros(2 * M + 1, np.uint8)
        for R in ok:
            if abs(R) <= M: t[R + M] = 1
        return t

    def gen(self, nw, seed, j, rok=None, cansign=True, mw=-1, sok=None, eok=None, fs=None, fe=None):
        s = np.zeros(self.n, np.int8); e = np.zeros(self.n, np.int8)
        t = self.rok_table(rok); ts = self.sok_table(sok, 2 * self.n + 2); te = self.sok_table(eok, self.n + 2)
        _lib.gl1_gen(self.n, self.r, self.nu, self.w.ctypes.data, self.forced.ctypes.data, self.poff.ctypes.data, self.plist.ctypes.data, self.dpw.ctypes.data,
                     self.nfree, self.freeidx.ctypes.data, self.fw, nw, mw, seed, j, (t.ctypes.data if t is not None else None), self.RM, int(cansign), (ts.ctypes.data if ts is not None else None), (te.ctypes.data if te is not None else None), 2 * self.n + 2, self.n + 2, (-100000 if fs is None else fs), (-100000 if fe is None else fe), s.ctypes.data, e.ctypes.data)
        return s, e

    def run(self, specs, per_class=30000.0, maxout=200000, hseed=1, stop_after=0):
        """specs: 4 BlockSpec (0,1 outer = small lists; 2,3 inner = large).  returns (count, S (m,4,n), E (m,4,n), stats)"""
        n = self.n
        keep = []
        fsv = np.array([-100000 if s.fs is None else s.fs for s in specs], np.int32); fev = np.array([-100000 if s.fe is None else s.fe for s in specs], np.int32)
        nw = np.array([s.nw for s in specs], np.int32); mwv = np.array([s.mw for s in specs], np.int32); L = np.array([s.L for s in specs], np.int64); seeds = np.array([s.seed for s in specs], np.uint64)
        Sp = (ctypes.c_void_p * 4)(); Ep = (ctypes.c_void_p * 4)(); Rp = (ctypes.c_void_p * 4)(); SOp = (ctypes.c_void_p * 4)(); EOp = (ctypes.c_void_p * 4)()
        cans = np.array([int(s.cansign) for s in specs], np.int32)
        for b, s in enumerate(specs):
            if s.S is not None:
                S_ = np.ascontiguousarray(s.S, np.int8); E_ = np.ascontiguousarray(s.E, np.int8); keep += [S_, E_]
                Sp[b] = S_.ctypes.data; Ep[b] = E_.ctypes.data
            t = self.rok_table(s.rok)
            if t is not None: keep.append(t); Rp[b] = t.ctypes.data
            ts = self.sok_table(s.sok, 2 * self.n + 2); te = self.sok_table(s.eok, self.n + 2)
            if ts is not None: keep.append(ts); SOp[b] = ts.ctypes.data
            if te is not None: keep.append(te); EOp[b] = te.ctypes.data
        pairs = max(float(L[0]) * L[2], float(L[1]) * L[3])
        lgM = max(1, int(math.ceil(math.log2(max(2.0, pairs / per_class))))); lgM = min(lgM, 30)
        oS = np.zeros((maxout, 4, n), np.int8); oE = np.zeros((maxout, 4, n), np.int8); stats = np.zeros(8)
        cnt = _lib.gl1_run(n, self.r, len(self.sh.repsPM), self.shift.ctypes.data, self.target.ctypes.data, self.nu, self.w.ctypes.data, self.forced.ctypes.data,
                           self.poff.ctypes.data, self.plist.ctypes.data, self.dpw.ctypes.data, self.nfree, self.freeidx.ctypes.data, self.fw,
                           nw.ctypes.data, mwv.ctypes.data, L.ctypes.data, seeds.ctypes.data, ctypes.addressof(Sp), ctypes.addressof(Ep), ctypes.addressof(Rp), self.RM, cans.ctypes.data, ctypes.addressof(SOp), ctypes.addressof(EOp), 2 * self.n + 2, self.n + 2, fsv.ctypes.data, fev.ctypes.data,
                           lgM, maxout, oS.ctypes.data, oE.ctypes.data, stats.ctypes.data, hseed, stop_after)
        m = min(cnt, maxout)
        return int(cnt), oS[:m], oE[:m], stats
