"""Exact quotient-level data for the 573 spin structure.  KEY FACT (proved below, checked by assertion): every K-orbit (H-orbit) of G lies over a
single K-orbit (H-orbit) R of Q=G/N, and meets every fibre over R in the same number c_O of points.  Hence the compressed vectors
a = (sum of A over fibre)  are  a_R = sum_{O over R} x_O c_O  : the set of attainable compressions is a PRODUCT over Q-orbits R of independent
value sets S_R (with exact multiplicities = number of sign patterns), cut only by the row sum  sum_R |R| a_R = rowsum.
So no DP over partial sums is needed (the old 400k cap was an artefact)."""
import numpy as np, math
from collections import defaultdict
from struct573 import build
from ss_core import closure_partition

def subgroup_from_flags(S, e2, e11, e13):
    G = S.G; c = G.coords
    ok = np.ones(G.v, bool)
    # N = (Z_2^2 if e2) x (Z_11 if e11) x (Z_13 if e13): elements (a, m) with  a=0 unless e2; m = 0 mod 11 unless e11; m = 0 mod 13 unless e13
    if not e2: ok &= (c[:, 0] == 0) & (c[:, 1] == 0)
    if not e11: ok &= (c[:, 2] % 11 == 0)
    if not e13: ok &= (c[:, 2] % 13 == 0)
    return np.flatnonzero(ok)

class Level:
    def __init__(self, S, Nidx, name=''):
        G = S.G; v = G.v
        self.S, self.name, self.N = S, name, np.asarray(Nidx)
        self.sizeN = len(self.N)
        cos = -np.ones(v, dtype=np.int64); reps = []
        for g in range(v):
            if cos[g] < 0:
                cos[G.index(G.coords[g] + G.coords[self.N])] = len(reps); reps.append(g)
        self.cos, self.reps = cos, np.array(reps); nQ = self.nQ = len(reps)
        self.Qadd = cos[G.index(G.coords[self.reps][:, None, :] + G.coords[self.reps][None, :, :])]
        act = lambda p: cos[p[self.reps]]
        for p in S.Hperms + [S.mu]:
            assert np.array_equal(cos[p], act(p)[cos]), 'N not invariant'
        self.pH = [act(p) for p in S.Hperms]; self.pmu = act(S.mu); self.pmu2 = self.pmu[self.pmu]
        self.negQ = cos[S.neg[self.reps]]
        self.orbH, self.RH = closure_partition(nQ, self.pH)
        self.orbK, self.RK = closure_partition(nQ, self.pH + [self.pmu])
        cl, _ = closure_partition(nQ, self.pH + [self.pmu, self.negQ])
        z = cl[0]; labs = sorted(set(cl.tolist()) - {z})
        self.qreps = [int(np.flatnonzero(cl == l)[0]) for l in labs]; self.D = len(self.qreps)
        self.qclass_size = [int((cl == l).sum()) for l in labs]
        self.wA = np.bincount(self.orbK, minlength=self.RK); self.wB = np.bincount(self.orbH, minlength=self.RH)
        # contributions
        self.valsA = self._value_sets(S.orbK, S.kK, self.orbK, self.RK)
        self.valsB = self._value_sets(S.orbH, S.kH, self.orbH, self.RH)
        # quadratic forms for each class
        self.MA = np.array([self._pm(self.orbK, self.RK, q) for q in self.qreps])
        self.MB = np.array([self._pm(self.orbH, self.RH, q) + self._pm(self.orbH, self.RH, int(self.pmu[q])) + self._pm(self.orbH, self.RH, int(self.pmu2[q])) for q in self.qreps])
        self.target = -4 * self.sizeN
        self.norm = 4 * S.G.v // 1 + 4 - 4 * self.sizeN            # 4n - 4|N|, n = v+1

    def _pm(self, orb, R, q):
        M = np.zeros((R, R), dtype=np.int64)
        np.add.at(M, (orb, orb[self.Qadd[:, q]]), 1)
        return M

    def _value_sets(self, orbG, kG, orbQ, R):
        S = self.S
        cnt = np.zeros((kG, self.nQ), dtype=np.int64)
        np.add.at(cnt, (orbG, self.cos), 1)
        contribs = [[] for _ in range(R)]
        for o in range(kG):
            nz = np.flatnonzero(cnt[o])
            rs = set(int(orbQ[q]) for q in nz)
            assert len(rs) == 1, 'orbit over several Q-orbits'
            r = rs.pop()
            cs = set(int(cnt[o, q]) for q in np.flatnonzero(orbQ == r))
            assert len(cs) == 1, 'non-constant fibre count on Q-orbit'
            contribs[r].append(cs.pop())
        res = []
        for r in range(R):
            d = {0: 1}
            for c in contribs[r]:
                nd = defaultdict(int)
                for val, w in d.items():
                    nd[val + c] += w; nd[val - c] += w
                d = nd
            res.append(dict(sorted(d.items())))
        self._contribs = contribs if orbG is S.orbK else getattr(self, '_contribs', None)
        return res

    def info(self):
        szA = [len(x) for x in self.valsA]; szB = [len(x) for x in self.valsB]
        return dict(name=self.name, nQ=self.nQ, sizeN=self.sizeN, RK=self.RK, RH=self.RH, D=self.D,
                    listA=float(np.prod([float(x) for x in szA])), listB=float(np.prod([float(x) for x in szB])),
                    qclass=self.qclass_size)

LEVELS = {  # name: (e2, e11, e13) flags of N
    'Z2^2': (0, 1, 1), 'Z11': (1, 0, 1), 'Z13': (1, 1, 0), 'Z143': (1, 0, 0), 'Z2^2xZ11': (0, 0, 1), 'Z2^2xZ13': (0, 1, 0)}
def make(S, name):
    return Level(S, subgroup_from_flags(S, *LEVELS[name]), name)

if __name__ == '__main__':
    S = build(3)
    for name in LEVELS:
        L = make(S, name)
        i = L.info()
        print(name, {k: (round(math.log10(i[k]), 2) if k.startswith('list') else i[k]) for k in i if k != 'qclass'}, 'qclass', i['qclass'])
        print('   A value-set sizes', [len(x) for x in L.valsA], ' B', [len(x) for x in L.valsB][:60])
