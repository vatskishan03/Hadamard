"""Index-r plug-in classes: G = K x F_q, H = {1} x S_r, S_r = subgroup of index r in F_q^* (acts trivially on K).
Block X <-> K-functions x_b, b = 0..r:  x_0 = X on K x {0},  x_{1+j} = X on K x C_j,  C_j = g^j S_r (g primitive).
Group-ring picture: X = sum_b x_b (x) B_b,  B_0 = delta_0,  B_{1+j} = 1_{C_j}.
PAF convention (as plug.Abel.paf): PAF_X(g) = sum_h X(h) X(h+g).  Then
  PAF_X(k, y) = sum_{b,c} CC_{bc}(k) N_{bc}(y),  CC_{bc}(k) = sum_{k'} x_b(k') x_c(k'+k),  N_{bc}(y) = sum_z B_b(z) B_c(z+y)
and N_{bc}(y) = M[a, b, c] for y in class a (a = 0: y = 0; a = 1+j: y in C_j).  K-level system:
  R_a(k) := sum_X sum_{bc} M[a,b,c] CC^X_{bc}(k) = T_a(k).
"""
import numpy as np, itertools, math
from math import gcd


def factor_prime_power(q):
    for p in range(2, q + 1):
        if q % p == 0:
            e = 0; t = q
            while t % p == 0: t //= p; e += 1
            if t != 1: return None
            return p, e
    return None


class GF:
    """GF(p^e); element index = sum d[t] p^(e-1-t), d[t] = coefficient of x^t  (so the additive group is Z_p^e in C order)."""
    def __init__(self, p, e):
        self.p, self.e, self.q = p, e, p ** e
        q = self.q
        if e == 1:
            self.mul = (np.arange(q)[:, None] * np.arange(q)[None, :]) % p
        else:
            red = None
            for coeffs in itertools.product(range(p), repeat=e):   # x^e = sum coeffs[t] x^t  <=> poly x^e - sum
                if coeffs[0] == 0: continue
                # irreducible iff x^e - sum has no factor: test via brute force that x has multiplicative order dividing q-1 and
                # that the field is a field (every nonzero element invertible)  -> build table and test
                tab = self._table(coeffs)
                if tab is None: continue
                red = coeffs; self.mul = tab; break
            assert red is not None
        self.add = None
        self._prim()

    def dig(self, i):
        e, p = self.e, self.p
        return [(i // p ** (e - 1 - t)) % p for t in range(e)]

    def idx(self, d):
        e, p = self.e, self.p
        return sum(int(d[t]) * p ** (e - 1 - t) for t in range(e))

    def _table(self, coeffs):
        p, e, q = self.p, self.e, self.q
        digs = [self.dig(i) for i in range(q)]
        tab = np.zeros((q, q), dtype=np.int64)
        for a in range(q):
            da = digs[a]
            for b in range(a, q):
                db = digs[b]
                prod = [0] * (2 * e - 1)
                for i in range(e):
                    for j in range(e):
                        prod[i + j] += da[i] * db[j]
                for t in range(2 * e - 2, e - 1, -1):
                    c = prod[t]; prod[t] = 0
                    for s in range(e): prod[t - e + s] += c * coeffs[s]
                v = self.idx([x % p for x in prod[:e]])
                tab[a, b] = v; tab[b, a] = v
        # field test: no zero divisors
        for a in range(1, q):
            if len(set(tab[a, 1:].tolist())) != q - 1: return None
        return tab

    def _prim(self):
        q = self.q
        one = self.idx([1] + [0] * (self.e - 1))
        self.one = one
        for g in range(2, q) if q > 2 else [1]:
            x = one; order = 0
            seen = set()
            while True:
                x = int(self.mul[x, g]); order += 1
                if x == one: break
                if order > q: break
            if order == q - 1:
                self.g = g; return
        if q == 2: self.g = 1; return
        raise RuntimeError('no primitive element')

    def neg(self, y):
        return self.idx([(-d) % self.p for d in self.dig(y)])


class Abel:
    def __init__(self, mods):
        self.mods = list(mods); self.n = int(np.prod(self.mods)) if self.mods else 1
        self.el = np.array(list(itertools.product(*[range(m) for m in self.mods])), dtype=np.int64).reshape(self.n, len(self.mods))
        w = np.ones(len(self.mods), dtype=np.int64)
        for i in range(len(self.mods) - 2, -1, -1): w[i] = w[i + 1] * self.mods[i + 1]
        self.w = w
        m = np.array(self.mods, dtype=np.int64)
        self.add = np.zeros((self.n, self.n), dtype=np.int64)
        for a in range(self.n):
            self.add[a] = (((self.el[a][None, :] + self.el) % m) * w).sum(1)
        self.neg = (((-self.el) % m) * w).sum(1)

    def paf(self, x):
        a = np.asarray(x, dtype=float).reshape(self.mods)
        f = np.fft.fftn(a)
        return np.rint(np.fft.ifftn(f * np.conj(f)).real).astype(np.int64).reshape(-1)

    def ccf(self, x, y):
        """CCF_{x,y}(g) = sum_h x(h) y(h+g)"""
        a = np.asarray(x, dtype=float).reshape(self.mods); b = np.asarray(y, dtype=float).reshape(self.mods)
        return np.rint(np.fft.ifftn(np.conj(np.fft.fftn(a)) * np.fft.fftn(b)).real).astype(np.int64).reshape(-1)


class Plug:
    def __init__(self, Kmods, q, r, bordered=False):
        pe = factor_prime_power(q); assert pe is not None
        self.p, self.e = pe; self.q = q; self.r = r
        assert (q - 1) % r == 0
        self.L = (q - 1) // r                      # |S_r| = |H|
        self.bordered = bordered
        self.Kmods = list(Kmods); self.K = Abel(self.Kmods); self.kap = self.K.n
        self.F = GF(self.p, self.e)
        F = self.F
        # classes
        self.cls = np.zeros(q, dtype=np.int64)       # 0 at y=0, 1+j on C_j
        x = F.one; self.rep = [0] * (r + 1); self.ind = np.full(q, -1, dtype=np.int64)
        for i in range(q - 1):
            self.cls[x] = 1 + (i % r); self.ind[x] = i
            if i < r: self.rep[1 + i] = x
            x = int(F.mul[x, F.g])
        self.rep[0] = 0
        # shift of -1:  -C_j = C_{j+delta}
        m1 = F.neg(F.one)
        self.delta = int((self.ind[m1]) % r)
        self.m1_in_S = (self.delta == 0)
        self.G = Abel(self.Kmods + [self.p] * self.e)
        self.m = self.kap * q
        self._tensor()
        self.negk = self.K.neg
        # K-reps modulo +-: used for equations
        self.reps = [k for k in range(self.kap) if k <= self.negk[k]]

    def _tensor(self):
        q, r, F = self.q, self.r, self.F
        B = np.zeros((r + 1, q), dtype=np.int64)
        for y in range(q): B[self.cls[y], y] = 1
        # N_{bc}(y) = sum_z B_b(z) B_c(z+y)
        add = np.zeros((q, q), dtype=np.int64)
        for a in range(q):
            da = F.dig(a)
            for b in range(q):
                add[a, b] = F.idx([(x + y) % F.p for x, y in zip(da, F.dig(b))])
        N = np.zeros((r + 1, r + 1, q), dtype=np.int64)
        for b in range(r + 1):
            for c in range(r + 1):
                for y in range(q):
                    N[b, c, y] = int((B[b] * B[c][add[:, y]]).sum())
        M = np.zeros((r + 1, r + 1, r + 1), dtype=np.int64)
        for a in range(r + 1):
            M[a] = N[:, :, self.rep[a]]
        # constancy on classes
        for y in range(q):
            assert (N[:, :, y] == M[self.cls[y]]).all()
        self.M = M; self.N = N; self.B = B

    # ---- block handling
    def expand(self, x):
        """x: (r+1, kap) -> full block on G, index k*q + y"""
        x = np.asarray(x, dtype=np.int64)
        out = x[self.cls, :].T   # (kap, q)
        return np.ascontiguousarray(out).reshape(-1)

    def target(self):
        """T_a(k) array (r+1, kap): unbordered (4m-2) delta; bordered -4 off 0, 4n-2 at 0 (n = |G|)"""
        T = np.zeros((self.r + 1, self.kap), dtype=np.int64)
        if not self.bordered:
            T[0, 0] = 4 * self.m - 2
        else:
            T[:, :] = -4; T[0, 0] = 4 * self.m - 2
        return T

    def CC(self, x):
        """CC[b,c,k] = sum_{k'} x_b(k') x_c(k'+k)"""
        r1 = self.r + 1
        out = np.zeros((r1, r1, self.kap), dtype=np.int64)
        for b in range(r1):
            for c in range(r1):
                out[b, c] = self.K.ccf(x[b], x[c])      # ccf(x,y)(g) = sum_h x(h) y(h+g)
        return out

    def Rlevel(self, blocks):
        """sum over blocks of sum_bc M[a,b,c] CC_bc(k)  -> (r+1, kap)"""
        tot = np.zeros((self.r + 1, self.kap), dtype=np.int64)
        for x in blocks:
            cc = self.CC(np.asarray(x))
            tot += np.einsum('abc,bck->ak', self.M, cc)
        return tot

    def residual(self, blocks):
        return self.Rlevel(blocks) - self.target()

    def Glevel(self, blocks):
        tot = np.zeros(self.G.n, dtype=np.int64)
        for x in blocks: tot += self.G.paf(self.expand(x))
        return tot

    def check_expansion(self, blocks):
        """K-level R_a(k) must equal the G-level PAF sum at (k, y), y in class a"""
        Rl = self.Rlevel(blocks); Gl = self.Glevel(blocks).reshape(self.kap, self.q)
        pred = Rl[self.cls, :].T
        return bool((pred == Gl).all())

    # ---- A/B symmetry
    def sym_ok(self, x):
        x = np.asarray(x)
        r = self.r
        for b in range(r + 1):
            bn = b if b == 0 else 1 + ((b - 1 + self.delta) % r)
            if not (x[bn, self.negk] == x[b]).all(): return False
        return True

    def sym_free_positions(self):
        """list of orbits of (b, k) under the A/B symmetry (b,k) ~ (b', -k): returns list of lists of (b,k)"""
        r = self.r; seen = set(); orbs = []
        for b in range(r + 1):
            for k in range(self.kap):
                if (b, k) in seen: continue
                bn = b if b == 0 else 1 + ((b - 1 + self.delta) % r)
                o = {(b, k), (bn, int(self.negk[k]))}
                seen |= o; orbs.append(sorted(o))
        return orbs

    def rowsums(self, x):
        return np.asarray(x).sum(1)

    # ---- balance data
    def sigma(self, x):
        return np.asarray(x)[1:].sum(0)

    def u_trivial(self, x):
        x = np.asarray(x)
        return x[0] + self.L * x[1:].sum(0)

    def balance_sum(self, blocks):
        """sum_{X,k} sigma (L sigma + 2 x0)   (must equal 4 r |K| unbordered, 4 r (|K| - 1) bordered... see derivation)"""
        s = 0
        for x in blocks:
            x = np.asarray(x); sg = x[1:].sum(0)
            s += int((sg * (self.L * sg + 2 * x[0])).sum())
        return s


def random_blocks(P, rng, sym=True):
    """four random +-1 blocks with AB symmetric/zero at (0,0)"""
    blocks = []
    for X in range(4):
        x = rng.choice([-1, 1], size=(P.r + 1, P.kap)).astype(np.int64)
        if X < 2 and sym:
            for orb in P.sym_free_positions():
                v = x[orb[0][0], orb[0][1]]
                for (b, k) in orb: x[b, k] = v
            x[0, 0] = 0
        blocks.append(x)
    return blocks


if __name__ == '__main__':
    rng = np.random.default_rng(5)
    for (Km, q, r) in [([3], 7, 3), ([5], 13, 3), ([3], 13, 4), ([7], 17, 4), ([3, 3], 7, 6), ([3], 25, 4), ([5], 25, 3), ([3], 9, 4), ([5], 19, 6), ([3], 13, 6), ([4], 5, 2), ([3], 7, 2)]:
        P = Plug(Km, q, r, bordered=(np.prod(Km) % 2 == 0))
        bl = random_blocks(P, rng)
        ok1 = P.check_expansion(bl)
        ok2 = all(P.sym_ok(b) for b in bl[:2])
        print(Km, 'q=', q, 'r=', r, 'L=', P.L, 'delta=', P.delta, 'expansion ok:', ok1, 'sym ok:', ok2, flush=True)
