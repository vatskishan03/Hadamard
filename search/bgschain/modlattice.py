"""Necessary condition mod q = p^k:  sum_b V_b[i_b] = T  with V_b[i] in V_b[0] + <differences>  =>  T - sum_b V_b[0]  in the Z/q-span of all differences of all blocks.
(Howell-style echelon form over Z/p^k; vectorised bulk reduction of all difference rows; exact for the modulus used.)"""
import numpy as np

def val(x, p, k):
    x = int(x) % (p ** k)
    if x == 0: return k
    v = 0
    while x % p == 0: x //= p; v += 1
    return v

class ModSpan:
    def __init__(self, p, k, n):
        self.p, self.k, self.q, self.n = p, k, p ** k, n
        self.basis = {}                      # col -> (row np.int64 array, valuation of the pivot); pivot entry = p^v exactly

    def _normalise(self, x, col):
        """scale so that x[col] = p^v"""
        p, q = self.p, self.q
        v = val(x[col], p, self.k); u = (int(x[col]) // p ** v) % q
        inv = pow(u, -1, q)
        return (x * inv) % q, v

    def insert(self, x):
        q, p = self.q, self.p
        stack = [np.asarray(x, dtype=np.int64) % q]
        while stack:
            x = stack.pop()
            for col in range(self.n):
                if x[col] % q == 0: continue
                vx = val(x[col], p, self.k)
                if col not in self.basis:
                    x2, v = self._normalise(x, col); self.basis[col] = (x2, v)
                    if v > 0: stack.append((x2 * p ** (self.k - v)) % q)
                    break
                b, vb = self.basis[col]
                if vx >= vb:
                    x = (x - (int(x[col]) // p ** vb) * b) % q
                else:
                    x2, v = self._normalise(x, col)
                    self.basis[col] = (x2, v)
                    if v > 0: stack.append((x2 * p ** (self.k - v)) % q)
                    # old pivot row is re-inserted after eliminating column col
                    nb = (b - (int(b[col]) // p ** v) * x2) % q
                    stack.append(nb)
                    break
            # fallthrough: x reduced to zero

    def reduce_rows(self, X):
        """reduce all rows of X (N x n) against the basis; returns the rows that are not in the span (after reduction)"""
        q, p = self.q, self.p
        X = X % q
        for col in sorted(self.basis):
            b, vb = self.basis[col]
            c = X[:, col]
            # rows whose entry at col has valuation < vb cannot be reduced here: leave them (they are not in the span)
            ok = (c % (p ** vb)) == 0
            qq = np.where(ok, c // p ** vb, 0)
            X = (X - np.outer(qq, b)) % q
        return X

    def contains(self, t):
        q, p = self.q, self.p
        t = np.asarray(t, dtype=np.int64) % q
        for col in sorted(self.basis):
            b, vb = self.basis[col]
            if t[col] % (p ** vb) != 0: return False
            t = (t - (int(t[col]) // p ** vb) * b) % q
        for col in range(self.n):
            if col not in self.basis and t[col] % q != 0 and False: pass
        return not np.any(t % q)

def span_of_differences(Vs, p, k, rounds=200):
    n = Vs[0].shape[1]; S = ModSpan(p, k, n); q = p ** k
    Ds = [((v - v[0]) % q).astype(np.int64) for v in Vs]
    rng = np.random.default_rng(1)
    for D in Ds:
        for i in rng.choice(len(D), size=min(40, len(D)), replace=False): S.insert(D[i])
    for D in Ds:
        for _ in range(rounds):
            R = S.reduce_rows(D)
            bad = np.flatnonzero(np.any(R % q != 0, axis=1))
            if len(bad) == 0: break
            for i in bad[:20]: S.insert(R[i])
    return S

def check(Vs, T, moduli=((2, 6), (3, 3), (5, 2), (7, 2), (13, 2))):
    base = sum(v[0] for v in Vs).astype(np.int64); t = np.asarray(T, dtype=np.int64) - base
    out = {}
    for p, k in moduli:
        S = span_of_differences(Vs, p, k)
        out[p ** k] = S.contains(t)
    return out
