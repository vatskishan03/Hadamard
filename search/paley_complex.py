#!/usr/bin/env python3
"""Quaternary (complex) Hadamard matrices of order (q + 1)/2 from Paley's skew Hadamard matrix of order q + 1,
q = p^k = 3 (mod 4) a prime power (complex-hadamard/PALEY.md).

Index PG(1, q) by the finite points x = 0..q-1 (representative (x, 1); for k > 1, x encodes the polynomial
x_0 + x_1 t + ... + x_(k-1) t^(k-1) by its base-p digits) and infinity = q (representative (1, 0)). With chi the
quadratic character of GF(q) (chi(0) = 0), S[u, v] = chi(det(u, v)) is Paley's skew conference matrix and
H = I + S his skew Hadamard matrix of order q + 1. An element g of SL(2, q) maps the representative of u to
lambda_u(g) times the representative of g(u), and P_g: e_u -> chi(lambda_u) e_(g(u)) satisfies P_g S P_g^T = S, so P_g
commutes with H. Since chi(-1) = -1, P_(-I) = -I, so for g = [[0, 1], [-1, 0]] (u -> -1/u, of order 4) P_g^2 = -I. In
the basis e_f, P_g e_f (f the smaller point of each pair {u, -1/u}, and of {0, infinity}) H = [[A, B], [-B, A]], and
K = (A + i B)/(1 + i) is a quaternary Hadamard matrix of order (q + 1)/2: K[k, l] = (H[f_k, f_l] + i s_l H[f_k,
g(f_l)])/(1 + i) with P_g e_(f_l) = s_l e_(g(f_l)).

    python3 search/paley_complex.py build 811             # check order 406 and print its SHA-256
    python3 search/paley_complex.py certify 811 [...]     # also write complex-hadamard/C406.json
    python3 search/paley_complex.py list                  # the orders 2N < 2000 (N odd) with 4N - 1 a prime power
"""
import hashlib, json, sys
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parent.parent


def prime_power(n):
    """(p, k) with n = p^k, or None"""
    for p in range(2, int(n ** 0.5) + 2):
        if n % p == 0:
            k = 0
            while n % p == 0:
                n //= p; k += 1
            return (p, k) if n == 1 else None
    return (n, 1) if n > 1 else None


class GF:
    """GF(p^k) on the integers 0..q-1 (base-p digits = polynomial coefficients modulo the monic irreducible f, the
    lexicographically first one with a primitive root t), with log and exp tables"""

    def __init__(self, p, k):
        self.p, self.k, self.q = p, k, p ** k
        q = self.q
        self.digits = np.array([[(x // p ** i) % p for i in range(k)] for x in range(q)], np.int64)
        self.pw = p ** np.arange(k)
        if k == 1:
            self.poly = [0, 1]
            g = next(g for g in range(2, p) if all(pow(g, (p - 1) // r, p) != 1 for r in self._primes(p - 1)))
            seq = [1]
            for _ in range(q - 2):
                seq.append(seq[-1] * g % p)
        else:
            for c in self._polys():
                seq = self._powers_of_t(c)
                if seq is not None:
                    self.poly = list(c) + [1]
                    break
        self.exp = np.array(seq + seq, np.int64)             # exp[i] = t^i, doubled to avoid a modulo
        self.log = np.full(q, -1, np.int64)
        self.log[np.array(seq)] = np.arange(q - 1)
        assert (self.log[1:] >= 0).all(), "not a primitive element"

    @staticmethod
    def _primes(n):
        out, d = [], 2
        while d * d <= n:
            if n % d == 0:
                out.append(d)
                while n % d == 0:
                    n //= d
            d += 1
        return out + ([n] if n > 1 else [])

    def _polys(self):
        p, k = self.p, self.k
        for code in range(p ** k):
            yield tuple((code // p ** i) % p for i in range(k))     # c_0, ..., c_(k-1) of t^k + ... + c_0

    def _powers_of_t(self, c):
        """powers 1, t, t^2, ... of t modulo t^k + c_(k-1) t^(k-1) + ... + c_0, if t has order q - 1"""
        p, k, q = self.p, self.k, self.q
        v = [1] + [0] * (k - 1)
        seq, seen = [], set()
        for _ in range(q - 1):
            code = sum(d * p ** i for i, d in enumerate(v))
            if code in seen:
                return None
            seen.add(code); seq.append(code)
            top = v[-1]                                          # multiply by t
            v = [0] + v[:-1]
            v = [(v[i] - top * c[i]) % p for i in range(k)]
        return seq if len(seen) == q - 1 and v == [1] + [0] * (k - 1) else None

    def add(self, x, y):
        return ((self.digits[x] + self.digits[y]) % self.p) @ self.pw

    def neg(self, x):
        return ((-self.digits[x]) % self.p) @ self.pw

    def sub(self, x, y):
        return self.add(x, self.neg(y))

    def mul(self, x, y):
        x, y = np.asarray(x), np.asarray(y)
        out = self.exp[self.log[np.where(x == 0, 1, x)] + self.log[np.where(y == 0, 1, y)]]
        return np.where((x == 0) | (y == 0), 0, out)

    def inv(self, x):
        return self.exp[(self.q - 1 - self.log[x]) % (self.q - 1)]

    def chi(self, x):
        x = np.asarray(x)
        return np.where(x == 0, 0, np.where(self.log[np.where(x == 0, 1, x)] % 2 == 0, 1, -1))


def paley_complex(q):
    """the quaternary Hadamard matrix of order (q + 1)/2 as an exponent matrix (entry i^e), and H, F"""
    pk = prime_power(q)
    if pk is None or q % 4 != 3:
        raise ValueError("q must be a prime power = 3 (mod 4)")
    F = GF(*pk)
    n = q + 1
    # representatives: (x, 1) for x < q, (1, 0) for infinity
    a = np.array(list(range(q)) + [1], np.int64)
    b = np.array([1] * q + [0], np.int64)
    # S[u, v] = chi(a_u b_v - b_u a_v)
    det = F.sub(F.mul(a[:, None], b[None, :]), F.mul(b[:, None], a[None, :]))
    H = F.chi(det).astype(np.int64) + np.eye(n, dtype=np.int64)
    # g = [[0, 1], [-1, 0]]: (x, y) -> (y, -x); the image point and lambda with g(rep u) = lambda rep(g u)
    img, lam = np.zeros(n, np.int64), np.zeros(n, np.int64)
    for u in range(n):
        w1, w2 = int(b[u]), int(F.neg(a[u]))
        if w2 != 0:
            img[u] = F.mul(w1, F.inv(w2)); lam[u] = w2
        else:
            img[u] = q; lam[u] = w1
    s = F.chi(lam)
    first = [u for u in range(n) if u < img[u]]
    assert len(first) == n // 2 and all(img[img[u]] == u for u in range(n))
    f = np.array(first)
    A = H[np.ix_(f, f)]
    B = s[f][None, :] * H[np.ix_(f, img[f])]
    # (A + iB)/(1 + i): (1, 1) -> 1, (-1, 1) -> i, (-1, -1) -> -1, (1, -1) -> -i
    E = np.select([(A == 1) & (B == 1), (A == -1) & (B == 1), (A == -1) & (B == -1)], [0, 1, 2], 3)
    return E, H, F


def check(E):
    """exact check: E is the exponent matrix of a quaternary matrix K with K K* = m I"""
    m = E.shape[0]
    re = np.select([E == 0, E == 2], [1.0, -1.0], 0.0)
    im = np.select([E == 1, E == 3], [1.0, -1.0], 0.0)
    # K K* = (re re^T + im im^T) + i (im re^T - re im^T); float64 sums of at most m terms of size 1 are exact
    G1 = re @ re.T + im @ im.T
    G2 = im @ re.T - re @ im.T
    return np.array_equal(G1, m * np.eye(m)) and not G2.any()


def canonical(E):
    return "\n".join("".join(str(int(e)) for e in row) for row in E) + "\n"


def build(q, write=False):
    E, H, F = paley_complex(q)
    m = E.shape[0]
    ok = check(E)
    sha = hashlib.sha256(canonical(E).encode()).hexdigest()
    print(f"q = {q} = {F.p}^{F.k}: quaternary Hadamard matrix of order {m}: {'verified' if ok else 'FAILED'}; "
          f"sha256 {sha}", flush=True)
    if not ok:
        raise SystemExit(1)
    if write:
        cert = {"order": m, "N": m // 2, "construction": "paley", "q": q, "p": F.p, "k": F.k, "sha256": sha}
        if F.k > 1:
            cert["poly"] = F.poly
        path = ROOT / "complex-hadamard" / f"C{m}.json"
        if path.exists():
            raise SystemExit(f"{path} exists; not overwritten")
        path.write_text(json.dumps(cert, indent=1) + "\n")
        print(f"  wrote {path}", flush=True)
    return sha


def main():
    if len(sys.argv) < 2:
        print(__doc__); return
    cmd = sys.argv[1]
    if cmd in ("build", "certify"):
        for x in sys.argv[2:]:
            build(int(x), write=(cmd == "certify"))
    elif cmd == "list":
        for N in range(1, 1000, 2):
            q = 4 * N - 1
            if prime_power(q):
                print(2 * N, q, prime_power(q))
    else:
        print(__doc__)


if __name__ == "__main__":
    main()
