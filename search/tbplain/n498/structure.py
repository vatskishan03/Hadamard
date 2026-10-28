#!/usr/bin/env python3
"""The structure behind QH(498): G = Z2 x Z3 x Z41 (order 246), N = 249, t = 3, Q = Z3, T = all three characters of Q.

Tables for the search, written from complex-hadamard/TBORDERED.md.  H = <(1, 1, 10)> of order 5 (10 has order 5 modulo 41),
A and B symmetric, A, B, C, D invariant under H.  Conditions: sum_k PAF_k(g) = -Gamma(phi(g)), Gamma = 4 sum_{chi in T} chi =
12 on ker(phi) and 0 elsewhere, i.e. -12 for g != 0 with b = 0, 0 for b != 0, 4 * 246 = 984 at g = 0.
Conventions: element g = (a, b, z) has index ((a * 3) + b) * 41 + z; row u = (a, b) has index r = 3a + b (0..5); a sequence is six
rows of 9-bit words, bit o of the word w set iff the sequence is -1 on the Z41-orbit o (under z -> 10 z), orbit 0 = {0}.
`classvec` gives the 28 PAF class values of a sequence by table look-up (CORR), `paf_full` the PAF by the definition and
`check_solution_full` the exact condition over G.  Running the file prints the orbit tables and self-tests the tables.

    python3 structure.py
"""
import itertools, math
import numpy as np

P = 41
HG = 10
assert pow(HG, 5, P) == 1 and pow(HG, 1, P) != 1

# ---- Z41 orbits under multiplication by 10
orb_of = [-1] * P
orbits = []
for z in range(P):
    if orb_of[z] < 0:
        o = []
        w = z
        while True:
            o.append(w)
            orb_of[w] = len(orbits)
            w = (w * HG) % P
            if w == z:
                break
        orbits.append(sorted(o))
NO = len(orbits)  # 9
assert NO == 9 and orbits[0] == [0]
OSIZE = np.array([len(o) for o in orbits])
OREP = np.array([o[0] for o in orbits])
NEGO = np.array([orb_of[(-o[0]) % P] for o in orbits])  # orbit of -z
assert (NEGO[NEGO] == np.arange(NO)).all() and NEGO[0] == 0 and all(NEGO[o] != o for o in range(1, NO))
orb_of = np.array(orb_of)

# ---- G elements, rows
def gidx(a, b, z):
    return ((a % 2) * 3 + (b % 3)) * P + (z % P)

GA = np.array([g // (3 * P) for g in range(246)])
GB = np.array([(g // P) % 3 for g in range(246)])
GZ = np.array([g % P for g in range(246)])
M = 246
NEG = np.array([gidx(-GA[g], -GB[g], -GZ[g]) for g in range(M)])
ADD = np.zeros((M, M), dtype=np.int16)
for g in range(M):
    ADD[g] = [gidx(GA[g] + GA[h], GB[g] + GB[h], GZ[g] + GZ[h]) for h in range(M)]

def row_index(a, b):
    return 3 * (a % 2) + (b % 3)

ROWA = np.array([r // 3 for r in range(6)])
ROWB = np.array([r % 3 for r in range(6)])
RADD = np.array([[row_index(ROWA[u] + ROWA[s], ROWB[u] + ROWB[s]) for u in range(6)] for s in range(6)])  # RADD[s][u] = u+s
RNEG = np.array([row_index(-ROWA[u], -ROWB[u]) for u in range(6)])

# ---- H-orbits and (H,-1)-classes on G; PAF classes
def h_act(g):
    return gidx(GA[g], GB[g], GZ[g] * HG)

orbH = -np.ones(M, dtype=int)
for g in range(M):
    if orbH[g] < 0:
        k = orbH.max() + 1
        w = g
        while orbH[w] < 0:
            orbH[w] = k
            w = h_act(w)
NH = orbH.max() + 1
assert NH == 54
# (s, O') -> index 9*r + O' is the H-orbit of g = (row r, z in O'); check consistency
for g in range(M):
    r = row_index(GA[g], GB[g])
    assert orbH[g] == orbH[gidx(GA[g], GB[g], OREP[orb_of[GZ[g]]])]
def horb_index(r, o):  # canonical id 9r+o
    return 9 * r + o
# PAF classes: (r, o) ~ (RNEG[r], NEGO[o])
cls_of = {}
classes = []
for r in range(6):
    for o in range(9):
        if (r, o) in cls_of:
            continue
        k = len(classes)
        classes.append([(r, o), (int(RNEG[r]), int(NEGO[o]))])
        cls_of[(r, o)] = k
        cls_of[(int(RNEG[r]), int(NEGO[o]))] = k
NCL = len(classes)  # 28
assert NCL == 28
CLS54 = np.array([cls_of[(r, o)] for r in range(6) for o in range(9)])  # class of H-orbit 9r+o
CLSIZE = np.zeros(NCL, dtype=int)
for g in range(M):
    CLSIZE[CLS54[9 * row_index(GA[g], GB[g]) + orb_of[GZ[g]]]] += 1
ZERO_CLASS = int(CLS54[0])  # row 0, orbit 0 = the zero element

# target of sum_k PAF_k
def target_of_class(c):
    (r, o) = classes[c][0]
    b = ROWB[r]
    if c == ZERO_CLASS:
        return 4 * M
    return -12 if b == 0 else 0

TARGET = np.array([target_of_class(c) for c in range(NCL)])
COLS = np.array([c for c in range(NCL) if c != ZERO_CLASS])  # 27 free columns
REPCLS = np.array([9 * classes[c][0][0] + classes[c][0][1] for c in range(NCL)])  # a representative 54-index per class

# ---- corr table: CORR[x*512+y, o'] = sum_z X(z) Y(z+g_{o'}) for words x,y (sign vectors constant on Z41-orbits)
SGN = np.array([[1 - 2 * ((w >> o) & 1) for o in range(9)] for w in range(512)], dtype=np.int64)  # 512 x 9, +1/-1
NN = np.zeros((9, 9, 9), dtype=np.int64)  # NN[o', o, o2] = #{z in o : z + g_{o'} in o2}
for op in range(9):
    g = int(OREP[op])
    for z in range(P):
        NN[op, orb_of[z], orb_of[(z + g) % P]] += 1
CORR = np.zeros((512, 512, 9), dtype=np.int16)
for op in range(9):
    CORR[:, :, op] = SGN @ NN[op] @ SGN.T
CORR = CORR.reshape(512 * 512, 9)

ROWSUM = SGN @ OSIZE  # row sums of the 512 words (odd)
REFL = np.array([sum(((w >> o) & 1) << NEGO[o] for o in range(9)) for w in range(512)])  # reflect z -> -z
SYMW = np.array([w for w in range(512) if REFL[w] == w])  # 32 symmetric words (5 free orbit bits)
assert len(SYMW) == 32


def paf54(R):
    """R: (n, 6) int array of words -> (n, 54) PAF values on the H-orbits, index 9*s + o' (s = row index of the shift)"""
    R = np.asarray(R)
    n = R.shape[0]
    out = np.zeros((n, 54), dtype=np.int32)
    for s in range(6):
        for u in range(6):
            idx = R[:, u].astype(np.int64) * 512 + R[:, RADD[s][u]]
            out[:, 9 * s:9 * s + 9] += CORR[idx]
    return out


def classvec(R):
    """(n, 28) PAF values on the 28 classes"""
    p = paf54(R)
    return p[:, REPCLS]


def words_to_full(Rw):
    """one sequence: 6 words -> +-1 vector over G (length 246)"""
    X = np.zeros(M, dtype=np.int64)
    for g in range(M):
        r = row_index(GA[g], GB[g])
        X[g] = SGN[Rw[r], orb_of[GZ[g]]]
    return X


def paf_full(X):
    """brute-force periodic autocorrelation of a +-1 vector over G (independent of the tables above)"""
    X = np.asarray(X, dtype=np.int64)
    return np.array([int((X * X[ADD[g]]).sum()) for g in range(M)])


def is_symmetric(X):
    return bool((X == X[NEG]).all())


def check_solution_full(Xs):
    """Xs: list of 4 full +-1 vectors; returns True iff sum PAF = -Gamma off 0 (exact, group convolution)"""
    tot = sum(paf_full(X) for X in Xs)
    tgt = np.array([4 * M if g == 0 else (-12 if GB[g] == 0 else 0) for g in range(M)])
    return bool((tot == tgt).all())


if __name__ == "__main__":
    print("orbits of Z41 under x10:", orbits)
    print("sizes", OSIZE.tolist(), "neg", NEGO.tolist())
    print("classes", NCL, "sizes", CLSIZE.tolist())
    print("target", TARGET.tolist())
    # self-test of the tables
    rng = np.random.default_rng(1)
    Rw = rng.integers(0, 512, size=(5, 6))
    cv = classvec(Rw)
    for i in range(5):
        X = words_to_full(Rw[i])
        pf = paf_full(X)
        for c in range(NCL):
            g = None
            r, o = classes[c][0]
            g = gidx(ROWA[r], ROWB[r], OREP[o])
            assert pf[g] == cv[i, c], (i, c, pf[g], cv[i, c])
        # PAF constant on classes
        for g in range(M):
            c = CLS54[9 * row_index(GA[g], GB[g]) + orb_of[GZ[g]]]
            assert pf[g] == cv[i, c]
    print("table self-test OK")
