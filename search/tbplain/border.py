#!/usr/bin/env python3
"""Exhaustive Q-level border of the character-bordered Szöllősi array (complex-hadamard/TBORDERED.md) for a small quotient Q.

Given the four compressions y_A, ..., y_D of a core (functions on Q: fibre sums of the four sequences), find quaternary
matrices Rh (2t x 2q) and Ch (2q x 2t), as exponent matrices (entry i^e), such that the rows of Rh and the columns of Ch lie
in V_T (each half of a row is a function on Q whose Fourier transform vanishes off T), Rh Rh* = Ch* Ch = 2q I, and
Xb = -Rh Mc* Ch / 2q is quaternary, Mc = [[X, Y], [Y*, -X*]], X = a X_A + conj(a) X_B, Y = a X_C + conj(a) X_D R, a = (1 + i)/2,
X_Z[a, b] = y_Z(b - a), R = [a -> -a].  The search is exhaustive over the rows and columns normalised to a first entry 1
(a unit scalar of a row or column changes nothing), in lexicographic order of the exponents, depth first, bitset pruned; the
border returned is the first one in that order.  It depends on the compressions only.  Q = Z2^2 with T the three nontrivial
characters (q = 4, t = 3, 324 candidates, the border of QH(358)) and Q = Z3 with T all characters (q = 3, t = 3, 1024
candidates, QH(498)) are used; any Q given by its moduli and any list T of characters works.  Library module for assemble.py;
running the file finds and checks the borders of the two recorded compression tuples.

    python3 border.py
"""
import itertools, time
import numpy as np

AL = (1 + 1j) / 2


class Quot:
    """Q = Z_m1 x ... x Z_mr (mixed radix, last factor fastest), T = list of characters given by their coordinate tuples"""
    def __init__(self, mods, T):
        self.mods = [int(m) for m in mods]; self.q = int(np.prod(self.mods)); self.T = [tuple(t) for t in T]
        self.el = np.array(list(itertools.product(*[range(m) for m in self.mods])), dtype=np.int64)       # element index -> coordinates
        idx = lambda c: np.ravel_multi_index(tuple(c[..., i] % self.mods[i] for i in range(len(self.mods))), tuple(self.mods))
        self.sub = idx(self.el[None, :, :] - self.el[:, None, :])          # sub[a, b] = index of b - a
        self.neg = idx(-self.el)
        self.R = np.zeros((self.q, self.q), dtype=np.int64); self.R[np.arange(self.q), self.neg] = 1
        self.chars = [np.exp(2j * np.pi * (self.el @ np.array([t[i] / self.mods[i] for i in range(len(self.mods))]))) for t in itertools.product(*[range(m) for m in self.mods])]
        self.charkeys = list(itertools.product(*[range(m) for m in self.mods]))
        self.t = len(self.T)

    def Mc(self, ys):
        """the 2q x 2q array Mc of the four compressions (each a vector of length q)"""
        XA, XB, XC, XD = (np.asarray(y, dtype=np.int64)[self.sub] for y in ys)
        X = AL * XA + np.conj(AL) * XB
        Y = AL * XC + np.conj(AL) * (XD @ self.R)
        return np.block([[X, Y], [Y.conj().T, -X.conj().T]])

    def half_vectors(self):
        """all f: Q -> {1, i, -1, -i} (exponent tuples) with f^(chi) = 0 for every character chi not in T"""
        out = []
        for e in itertools.product(range(4), repeat=self.q):
            z = 1j ** np.array(e)
            if all(abs((z * np.conj(ch)).sum()) < 1e-9 for ch, key in zip(self.chars, self.charkeys) if key not in self.T):
                out.append(e)
        return out

    def candidates(self):
        """rows of Rh (and, transposed, columns of Ch): exponent vectors of length 2q, both halves in V_T, first entry 1; lexicographic"""
        V = self.half_vectors()
        return np.array([e1 + e2 for e1 in V for e2 in V if e1[0] == 0], dtype=np.int64)


def isunit(z):
    return (np.abs(np.abs(z) - 1) < 1e-9) & ((np.abs(z.real) < 1e-9) | (np.abs(z.imag) < 1e-9))


def border_search(Qt, ys):
    """first border in the search order, as (Rh, Ch) exponent matrices, or None if the compression tuple has no border"""
    q, t = Qt.q, Qt.t
    Mc = Qt.Mc(ys)
    cand = Qt.candidates()
    Z = 1j ** cand
    W = Z @ Mc.conj().T                                 # r Mc*
    Vals = -(W @ Z.T) / (2 * q)                         # Vals[r, c] = -(r Mc* c) / 2q, c a column vector with the entries of Z[c]
    OK = isunit(Vals)
    ORT = np.abs(Z @ Z.conj().T) < 1e-9                  # r r'* = 0 (rows) and likewise for columns
    n = len(Z)
    okbits = [int(''.join('1' if x else '0' for x in OK[r][::-1]), 2) for r in range(n)]
    ortbits = [int(''.join('1' if x else '0' for x in ORT[r][::-1]), 2) for r in range(n)]
    k = 2 * t
    def find_cols(mask):
        def rec(chosen, cand_mask):
            if len(chosen) == k: return list(chosen)
            if bin(cand_mask).count('1') < k - len(chosen): return None
            m = cand_mask
            while m:
                l = m & -m; c = l.bit_length() - 1; m ^= l
                r = rec(chosen + [c], m & ortbits[c])
                if r: return r
            return None
        return rec([], mask)
    sols = []
    def dfs(chosen, cand_mask, colmask):
        if len(chosen) == k:
            cols = find_cols(colmask)
            if cols: sols.append((list(chosen), cols)); return True
            return False
        m = cand_mask
        while m:
            l = m & -m; r = l.bit_length() - 1; m ^= l
            cm = colmask & okbits[r]
            if bin(cm).count('1') < k: continue
            if dfs(chosen + [r], m & ortbits[r], cm): return True
        return False
    full = (1 << n) - 1
    dfs([], full, full)
    if not sols: return None
    rows, cols = sols[0]
    return cand[rows], cand[cols].T


def finish(Qt, ys, Rh, Ch):
    """Xb = -Rh Mc* Ch / 2q as an exponent matrix (checks that every entry is a unit), and the exact orthogonality of the border"""
    q = Qt.q
    Mc = Qt.Mc(ys)
    Rm = 1j ** Rh; Cm = 1j ** Ch
    assert np.abs(Rm @ Rm.conj().T - 2 * q * np.eye(len(Rm))).max() < 1e-9 and np.abs(Cm.conj().T @ Cm - 2 * q * np.eye(len(Cm[0]))).max() < 1e-9
    Xb = -(Rm @ Mc.conj().T @ Cm) / (2 * q)
    assert isunit(Xb).all()
    return np.rint(np.angle(Xb) / (np.pi / 2)).astype(int) % 4


if __name__ == "__main__":
    import json, sys
    # the compressions of the two recorded hits: (y(0), y(a != 0)) per block for QH(358) (spike form, Q = Z2^2), the fibre sums for QH(498) (Z3)
    Q4 = Quot([2, 2], [(0, 1), (1, 0), (1, 1)]); Q3 = Quot([3], [(0,), (1,), (2,)])
    cases = [("358", Q4, [[6, 4, 4, 4], [0, 0, 0, 0], [2, 4, 4, 4], [2, 4, 4, 4]]),
             ("498", Q3, [[2, 0, 0], [0, 0, 0], [2, 0, 0], [2, 0, 0]])]
    for name, Qt, ys in cases:
        t0 = time.time()
        r = border_search(Qt, ys)
        assert r is not None
        Xb = finish(Qt, ys, *r)
        print("QH(%s): border found, Rh %s Ch %s Xb %s, %.1fs" % (name, r[0].shape, r[1].shape, Xb.shape, time.time() - t0))
