"""Write a Yamada (Kiyasu) certificate: a Hadamard matrix of order 4(q + 2) from one of the skew-Hadamard
certificates in skew/ (doubled k times), with q = 2s - 3 a prime or a prime power, s the skew order. For q = p^k
(k > 1) the certificate also gives p and a primitive polynomial of degree k over GF(p), and r and xi are encoded as
integers whose base-p digits are polynomial coefficients (see verify.build_yamada).

  python3 attack/make_yamada.py SKEW_ORDER K [OUTDIR]

The certificate nests the skew certificate (and k skew-doubling steps). verify/verify.py rebuilds everything:
skew matrix -> doublings -> Yamada's construction -> H H^T = nI -> SHA-256. The matrix itself is not stored."""
import hashlib, json, math, sys
from pathlib import Path
ROOT = Path(__file__).resolve().parent.parent
sys.path[:0] = [str(ROOT / "verify")]
import numpy as np
import verify


def is_prime(n):
    return n > 1 and all(n % d for d in range(2, math.isqrt(n) + 1))


def primitive(q, r):
    mul = lambda x, y: ((x[0] * y[0] + r * x[1] * y[1]) % q, (x[0] * y[1] + x[1] * y[0]) % q)
    def power(x, e):
        y = (1, 0)
        while e:
            if e & 1: y = mul(y, x)
            x, e = mul(x, x), e >> 1
        return y
    N = q * q - 1; ps = [p for p in range(2, math.isqrt(N) + 2) if N % p == 0 and is_prime(p)]
    m = N
    for p in ps:
        while m % p == 0: m //= p
    if m > 1: ps.append(m)
    for b in range(1, q):
        for a in range(q):
            if all(power((a, b), N // p) != (1, 0) for p in ps): return [a, b]


def prime_power(q):
    for p in range(2, math.isqrt(q) + 2):
        if q % p == 0:
            k, x = 0, q
            while x % p == 0: x //= p; k += 1
            return (p, k) if x == 1 else None
    return (q, 1)


def primitive_poly(p, k):
    """Coefficients c_0 .. c_{k-1} of a primitive x^k + c_{k-1} x^{k-1} + ... + c_0 over GF(p) (x generates GF(p^k)*)."""
    import itertools
    q = p ** k
    for coeffs in itertools.product(range(p), repeat=k):
        if coeffs[0] == 0: continue
        e, seen = 1, 0
        for _ in range(q - 1):
            d = [0] + [(e // p ** j) % p for j in range(k)]
            top = d[k]; e = sum(((d[j] - top * coeffs[j]) % p) * p ** j for j in range(k))
            seen += 1
            if e == 1: break
        if e == 1 and seen == q - 1: return list(coeffs)


def field_params(q, p, poly):
    """A non-square r of GF(q) (the element x) and a primitive xi of GF(q^2) = GF(q)[w]/(w^2 - r)."""
    k = len(poly)
    digits = np.array([[(e // p ** j) % p for j in range(k)] for e in range(q)], dtype=np.int64)
    wts = p ** np.arange(k)
    enc = lambda d: int((np.asarray(d) % p) @ wts)
    exp = np.empty(q - 1, dtype=np.int64); log = -np.ones(q, dtype=np.int64); e = 1
    for i in range(q - 1):
        exp[i], log[e] = e, i
        d = np.zeros(k + 1, dtype=np.int64); d[1:] = digits[e]; top = d[k]; e = enc(d[:k] - top * np.asarray(poly))
    add = lambda u, v: enc(digits[u] + digits[v])
    fm = lambda u, v: 0 if u == 0 or v == 0 else int(exp[(log[u] + log[v]) % (q - 1)])
    r = int(exp[1])
    mul = lambda x, y: (add(fm(x[0], y[0]), fm(r, fm(x[1], y[1]))), add(fm(x[0], y[1]), fm(x[1], y[0])))
    def power(x, e):
        y = (1, 0)
        while e:
            if e & 1: y = mul(y, x)
            x, e = mul(x, x), e >> 1
        return y
    N = q * q - 1; ps = [f for f in range(2, math.isqrt(N) + 2) if N % f == 0 and is_prime(f)]
    m = N
    for f in ps:
        while m % f == 0: m //= f
    if m > 1: ps.append(m)
    for b in range(1, q):
        for a in range(q):
            if all(power((a, b), N // f) != (1, 0) for f in ps): return r, [a, b]


def main():
    s0, k = int(sys.argv[1]), int(sys.argv[2])
    base = json.loads((ROOT / "skew" / str(s0) / "certificate.json").read_text())
    Q = base
    for _ in range(k):
        Q = {"order": 2 * Q["order"], "construction": "skew-doubling", "S": Q}
        H, _ = verify.BUILDERS["skew-doubling"](Q)
        sha = hashlib.sha256()
        for ch in verify.csv_chunks(H): sha.update(ch)
        Q["sha256"] = sha.hexdigest(); Q["skew"] = True
        del H
    s = Q["order"]; q = 2 * s - 3
    pk = prime_power(q)
    if not (pk and q % 8 == 5): sys.exit(f"q = {q} is not a prime power = 5 (mod 8)")
    order = 4 * (q + 2)
    if pk[1] == 1:
        r = next(x for x in range(2, q) if pow(x, (q - 1) // 2, q) == q - 1)
        xi = primitive(q, r)
        cert = {"order": order, "construction": "yamada", "q": q, "r": r, "xi": xi, "Q": Q}
    else:
        poly = primitive_poly(*pk)
        r, xi = field_params(q, pk[0], poly)
        cert = {"order": order, "construction": "yamada", "q": q, "p": pk[0], "poly": poly, "r": r, "xi": xi, "Q": Q}
    H, _ = verify.BUILDERS["yamada"](cert)
    assert H.shape[0] == order and verify.is_pm1(H) and verify.gram_is_scalar(H, order), "construction failed"
    sha = hashlib.sha256()
    for ch in verify.csv_chunks(H): sha.update(ch)
    cert["sha256"] = sha.hexdigest()
    out = Path(sys.argv[3]) if len(sys.argv) > 3 else ROOT / str(order)
    out.mkdir(parents=True, exist_ok=True)
    (out / "certificate.json").write_text(json.dumps(cert, indent=1) + "\n")
    print(f"order {order} = 4*{q + 2}: q = {q}, r = {r}, xi = {xi}, skew order {s} = 2^{k} * {s0}, sha256 {cert['sha256']}")


if __name__ == "__main__":
    main()
