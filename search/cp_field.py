"""Cyclotomic products with a prime-power fibre: Z_m x GF(q), q = p^k (k >= 1), for cp_lowe.py and cp_hunt.py.

FieldStructure has the interface of cyclo_product.Structure (classes of index e of GF(q)*, slice patterns, the
correlation table CC, K_m = <kgen> on Z_m, twists (a, y) -> (kgen a, g^t y)), but the fibre is the additive group
Z_p^k of GF(q). The group is G = Z_m x Z_p^k, listed in mixed-radix order: the element (a, c_0, .., c_(k-1)) has index
a q + ravel(c), where c are the coordinates of y in the basis 1, w, .., w^(k-1) of GF(q) over GF(p) (w a root of a
primitive polynomial). Multiplication by w^t is a k x k matrix over Z_p, so every multiplier of H is an automorphism
of G given by an integer matrix on the coordinates (the form verify.build_goethals_seidel_group checks).

With k = 1 this is the prime field in coordinates, which lets the group code be checked against cyclo_product.
"""
import itertools, math
import numpy as np
from pathlib import Path as _P
import sys
sys.path[:0] = [str(_P(__file__).resolve().parent)]
import cyclo_product as CP


def mpow(M, e, p):
    R = np.eye(len(M), dtype=np.int64); B = M % p
    while e:
        if e & 1: R = R @ B % p
        B = B @ B % p; e >>= 1
    return R


def primitive_companion(p, k):
    """Companion matrix of a primitive polynomial x^k + c_(k-1) x^(k-1) + .. + c_0 over GF(p): multiplication by w."""
    q = p ** k
    primes = [f for f in range(2, q) if (q - 1) % f == 0 and all(f % r for r in range(2, int(f ** 0.5) + 1))]
    I = np.eye(k, dtype=np.int64)
    for coeffs in itertools.product(range(p), repeat=k):
        if coeffs[0] == 0: continue
        M = np.zeros((k, k), dtype=np.int64)
        if k > 1: M[1:, :-1] = np.eye(k - 1, dtype=np.int64)
        M[:, -1] = [(-c) % p for c in coeffs]
        if k == 1 and (M[0, 0] in (0, 1) and q > 2): continue
        if all(not np.array_equal(mpow(M, (q - 1) // f, p), I) for f in primes) and np.array_equal(mpow(M, q - 1, p), I):
            return M
    raise ValueError(f'no primitive polynomial of degree {k} over GF({p})')


class FieldStructure:
    def __init__(self, v, p, k, e, skew=True, kgen=1):
        q = p ** k
        assert v % q == 0 and (v // q) % p, "q = p^k must divide v exactly"
        assert (q - 1) % e == 0 and e % 2 == 0
        assert not skew or ((q - 1) // e) % 2 == 1, "skew needs (q - 1)/e odd"
        self.v, self.p, self.k, self.q, self.e, self.m, self.skew = v, p, k, q, e, v // q, skew
        self.F = (q - 1) // e
        self.Wm = primitive_companion(p, k)                       # multiplication by w on coordinates
        # coordinates of the elements, and w^i for i = 0 .. q-2
        self.coords = np.array(np.unravel_index(np.arange(q), (p,) * k)).T.astype(np.int64)    # [q, k]
        self.ravel = lambda C: np.ravel_multi_index(tuple((np.asarray(C) % p).T), (p,) * k)
        x = np.zeros(k, dtype=np.int64); x[0] = 1                 # the element 1 = (1, 0, .., 0)
        self.exp = np.zeros(q - 1, dtype=np.int64); self.log = np.full(q, -1, dtype=np.int64)
        for i in range(q - 1):
            idx = int(self.ravel(x[None, :])[0]); self.exp[i] = idx; self.log[idx] = i
            x = self.Wm @ x % p
        assert (self.log[1:] >= 0).all() if False else (np.sort(self.exp) == np.arange(1, q)).all() or True
        assert len(set(self.exp.tolist())) == q - 1 and self.log[self.ravel(np.zeros((1, k), np.int64))[0]] == -1
        self.cls = np.full(q, -1, dtype=np.int64)
        nz = self.log >= 0
        self.cls[nz] = self.log[nz] % e
        self.zero = int(self.ravel(np.zeros((1, k), np.int64))[0])
        minus1 = int(self.ravel((-self.coords[self.exp[0]])[None, :])[0])
        self.half = int(self.cls[minus1])
        self.nb = e // 2 if self.half else e
        self.reps = [self.zero] + [int(self.exp[c]) for c in range(e)]     # w^c represents C_c
        self.g = None                                             # twists act through Wm (see mul_matrix)
        m = self.m; self.kgen = kgen % m if m > 1 else 1
        K = {1}
        while m > 1:
            K2 = K | {kk * self.kgen % m for kk in K}
            if K2 == K: break
            K = K2
        self.K = sorted(K)
        self.korb = np.full(m, -1, dtype=np.int64); self.kreps = []
        for a in range(m):
            if self.korb[a] < 0:
                for kk in self.K: self.korb[a * kk % m] = len(self.kreps)
                self.kreps.append(a)
        if skew and m > 1:
            assert all(self.korb[(-a) % m] != self.korb[a] for a in range(1, m)), "-1 in K_m: no skew block"
        # addition table on indices (q <= 2500 here)
        C = self.coords
        self.add = self.ravel((C[:, None, :] + C[None, :, :]).reshape(-1, k)).reshape(q, q)
        self.neg_f = self.ravel(-C)
        self.bkey = []; self.dkey = self.nb

    # ---- patterns (identical to cyclo_product.Structure) --------------------------------------------------------
    pattern = CP.Structure.pattern
    patterns = CP.Structure.patterns
    skew0 = CP.Structure.skew0
    partner = CP.Structure.partner

    def mul_matrix(self, t):
        """Coordinates matrix of the multiplication by w^t."""
        return mpow(self.Wm, t % (self.q - 1), self.p)

    def cc_table(self, used):
        """CC[i, j, c] = sum_y U_i(y + x_c) U_j(y) over GF(q), x_0 = 0 and x_(1+c) = w^c."""
        q, p, k = self.q, self.p, self.k
        V = np.empty((len(used), q))
        for i, pid in enumerate(used):
            f, xs = self.pattern(pid)
            V[i, self.zero] = f
            nz = self.cls >= 0
            V[i, nz] = np.asarray(xs)[self.cls[nz]]
        Fv = np.fft.fftn(V.reshape((len(used),) + (p,) * k), axes=tuple(range(1, k + 1)))
        out = np.empty((len(used), len(used), len(self.reps)), dtype=np.int16)
        for i in range(len(used)):
            R = np.fft.ifftn(Fv[i][None] * np.conj(Fv), axes=tuple(range(1, k + 1))).real.reshape(len(used), q)
            out[i] = np.rint(R[:, self.reps]).astype(np.int16)
        return out

    def sequence(self, slices):
        """+-1 sequence on G (mixed-radix order, Z_m slowest) from the pattern ids at the positions a = 0 .. m-1."""
        P = np.array([[self.pattern(pp)[0]] + self.pattern(pp)[1] for pp in slices])     # m x (e + 1)
        col = np.where(self.cls < 0, 0, 1 + self.cls)                                    # per fibre element
        return P[:, col].reshape(-1)

    # ---- the group G = Z_m x Z_p^k ------------------------------------------------------------------------------
    @property
    def mods(self):
        return ([self.m] if self.m > 1 else []) + [self.p] * self.k

    def paf_sum(self, seqs):
        """Sum of the periodic autocorrelations on G, flattened (index 0 is the identity)."""
        shape = tuple(self.mods)
        tot = 0
        for x in seqs:
            X = np.fft.fftn(np.asarray(x, float).reshape(shape))
            tot = tot + np.rint(np.fft.ifftn(np.abs(X) ** 2).real).astype(np.int64)
        return tot.reshape(-1)

    def neg(self):
        """Index permutation of t -> -t on G."""
        a = np.repeat(np.arange(self.m), self.q); y = np.tile(np.arange(self.q), self.m)
        return ((-a) % self.m) * self.q + self.neg_f[y]

    def automorphisms(self, twist=0):
        """Integer matrices (on the coordinates of G) of the generators of H: w^e on the fibre, and (kgen, w^twist)."""
        k = self.k; out = []
        def block(a_mult, M):
            if self.m == 1: return M.tolist()
            B = np.zeros((k + 1, k + 1), dtype=np.int64); B[0, 0] = a_mult; B[1:, 1:] = M
            return B.tolist()
        out.append(block(1, self.mul_matrix(self.e)))
        if len(self.K) > 1: out.append(block(self.kgen, self.mul_matrix(twist)))
        return out

    def perm_of(self, M):
        """Index permutation of G for the coordinate matrix M."""
        mods = np.array(self.mods, dtype=np.int64)
        coords = np.array(np.unravel_index(np.arange(self.v), tuple(int(x) for x in mods))).T.astype(np.int64)
        return np.ravel_multi_index(tuple(((coords @ np.array(M, dtype=np.int64).T) % mods).T), tuple(int(x) for x in mods))


def fibre_structures(v, skew, maxpat=2e5):
    """(p, k, e, F) with q = p^k exactly dividing v, k >= 2, e | q - 1 even, F = (q-1)/e >= 3 (odd when skew) and at most
    maxpat class patterns per position."""
    out = []
    n, p = v, 2
    fac = {}
    while p * p <= n:
        while n % p == 0: fac[p] = fac.get(p, 0) + 1; n //= p
        p += 1
    if n > 1: fac[n] = fac.get(n, 0) + 1
    for p, k in fac.items():
        if k < 2: continue
        q = p ** k
        for e in range(2, q, 2):
            if (q - 1) % e: continue
            F = (q - 1) // e
            if F < 3 or (skew and F % 2 == 0) or math.comb(e, e // 2) > maxpat: continue
            out.append((p, k, e, F))
    return out
