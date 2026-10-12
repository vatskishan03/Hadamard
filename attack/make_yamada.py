"""Write a Yamada (Kiyasu) certificate: a Hadamard matrix of order 4(q + 2) from one of the skew-Hadamard
certificates in skew/ (doubled k times), with q = 2s - 3 prime, s the skew order.

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
    if not (is_prime(q) and q % 8 == 5): sys.exit(f"q = {q} is not a prime = 5 (mod 8)")
    r = next(x for x in range(2, q) if pow(x, (q - 1) // 2, q) == q - 1)
    xi = primitive(q, r)
    order = 4 * (q + 2)
    cert = {"order": order, "construction": "yamada", "q": q, "r": r, "xi": xi, "Q": Q}
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
