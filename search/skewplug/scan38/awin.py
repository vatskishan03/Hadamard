"""window-aware allocation generator: all allocations (nA, nB <= nC <= nD) with total defect weight E in [Emax - erel, Emax] passing the defect identity (same as sk_pat.allocations restricted to the window, without enumerating the rest)"""
import numpy as np
def alloc_window(sh, erel=5, want_all=False, erange=None):
    N, r = sh.N, sh.r
    Wf = sorted(set(int(x) for x in sh.kCD.weights())); WA = sorted(set(int(x) for x in sh.A_weights()))
    if not Wf or not WA: return [], None, None
    size = WA[-1] + 3 * Wf[-1] + 1
    def arr(ws):
        a = np.zeros(size, dtype=bool); a[ws] = True; return a
    def conv(a, b):
        out = np.zeros(size, dtype=bool)
        for i in np.nonzero(a)[0]:
            out[i:] |= b[:size - i]
        return out
    aF = arr(Wf); a2 = conv(aF, aF); a3 = conv(a2, aF); aT = conv(arr(WA), a3)
    valid = []
    for E in np.nonzero(aT)[0]:
        E = int(E); num = (r + 1) * E - 4 * N
        if num % 4: continue
        em = num // 4
        if 0 <= em <= E: valid.append(E)
    if not valid: return [], None, None
    Emax = max(valid); Emin = min(valid)
    lo = Emin if want_all else Emax - erel
    hi = Emax
    if erange is not None: lo, hi = erange
    Wf_arr = np.array(Wf, dtype=np.int64)
    inW = np.zeros(size + 1, dtype=bool); inW[Wf_arr] = True
    out = []
    for E in valid:
        if E < lo or E > hi: continue
        em = ((r + 1) * E - 4 * N) // 4
        for nA in WA:
            rem = E - nA
            if rem < 3 * Wf[0]: continue
            for i, w1 in enumerate(Wf):
                if 3 * w1 > rem: break
                w2 = Wf_arr[i:]
                w3 = rem - w1 - w2
                ok = (w3 >= w2) & (w3 >= 0)
                if not ok.any(): continue
                idx = np.nonzero(ok)[0]
                good = inW[w3[idx]]
                for j in idx[good]:
                    out.append(((nA, w1, int(w2[j]), int(w3[j])), em))
    return out, Emin, Emax


class AllocSampler:
    """lazy random allocations in an energy range [elo, ehi] (no enumeration): E uniform among valid totals, nA uniform among A weights, (w1 <= w2 <= w3) uniform among the solutions for the drawn (E, nA, w1)"""
    def __init__(self, sh, elo, ehi):
        N, r = sh.N, sh.r
        self.N, self.r = N, r
        self.Wf = sorted(set(int(x) for x in sh.kCD.weights())); self.WA = sorted(set(int(x) for x in sh.A_weights()))
        size = self.WA[-1] + 3 * self.Wf[-1] + 2
        self.Wf_arr = np.array(self.Wf, dtype=np.int64)
        self.inW = np.zeros(size + 1, dtype=bool); self.inW[self.Wf_arr] = True
        self.valid = []
        for E in range(max(elo, 0), min(ehi, size - 1) + 1):
            num = (r + 1) * E - 4 * N
            if num % 4: continue
            em = num // 4
            if 0 <= em <= E and E >= self.WA[0] + 3 * self.Wf[0]: self.valid.append(E)
    def em_of(self, E): return ((self.r + 1) * E - 4 * self.N) // 4
    def draw(self, rng):
        for _ in range(500):
            if not self.valid: return None
            E = rng.choice(self.valid); nA = rng.choice(self.WA); rem = E - nA
            if rem < 3 * self.Wf[0]: continue
            w1 = rng.choice(self.Wf)
            if 3 * w1 > rem: continue
            w2 = self.Wf_arr[self.Wf_arr >= w1]
            w3 = rem - w1 - w2
            ok = (w3 >= w2) & (w3 >= 0) & (w3 < len(self.inW))
            idx = np.nonzero(ok)[0]
            if len(idx) == 0: continue
            idx = idx[self.inW[w3[idx]]]
            if len(idx) == 0: continue
            j = idx[rng.randrange(len(idx))]
            return ((nA, w1, int(w2[j]), int(w3[j])), self.em_of(E))
        return None


def valid_totals(sh):
    """(valid E list, Emin, Emax) from the reachable sums of the orbit weights and the defect identity (cheap: boolean convolutions)"""
    N, r = sh.N, sh.r
    Wf = sorted(set(int(x) for x in sh.kCD.weights())); WA = sorted(set(int(x) for x in sh.A_weights()))
    if not Wf or not WA: return [], None, None
    size = WA[-1] + 3 * Wf[-1] + 1
    def arr(ws):
        a = np.zeros(size, dtype=bool); a[ws] = True; return a
    def conv(a, b):
        out = np.zeros(size, dtype=bool)
        for i in np.nonzero(a)[0]:
            out[i:] |= b[:size - i]
        return out
    aF = arr(Wf); aT = conv(arr(WA), conv(conv(aF, aF), aF))
    valid = []
    for E in np.nonzero(aT)[0]:
        E = int(E); num = (r + 1) * E - 4 * N
        if num % 4: continue
        em = num // 4
        if 0 <= em <= E: valid.append(E)
    return (valid, valid[0], valid[-1]) if valid else ([], None, None)

class AllocSamplerE(AllocSampler):
    def __init__(self, sh, Elist):
        N, r = sh.N, sh.r
        self.N, self.r = N, r
        self.Wf = sorted(set(int(x) for x in sh.kCD.weights())); self.WA = sorted(set(int(x) for x in sh.A_weights()))
        size = self.WA[-1] + 3 * self.Wf[-1] + 2
        self.Wf_arr = np.array(self.Wf, dtype=np.int64)
        self.inW = np.zeros(size + 1, dtype=bool); self.inW[self.Wf_arr] = True
        self.valid = list(Elist)
