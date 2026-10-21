#!/usr/bin/env python3
"""Symmetric Hadamard matrices of order 4q, q = 1 (mod 8) a prime power, from a symmetric Hadamard matrix of order
(q - 1)/2: Seberry and Yamada's Corollary 29(ii) Williamson-type matrices (W_3 = W_4) in Seberry and Balonin's propus
array (verify.py construction `symmetric-propus-miyamoto`; NOTES.md (top level) section 29).

    python3 search/sym_propus_miyamoto.py targets          census-undecided v = q with an input we can build
    python3 search/sym_propus_miyamoto.py certify v [...]  writes symmetric/<4v>/certificate.json

Inputs of order n = (q - 1)/2: Sylvester's [[1, 1], [1, -1]]^k times one symmetric Hadamard matrix of Paley's (the
symmetric member of the amicable pair of order p + 1, p = 3 (mod 4), or Paley II of order 2(p + 1), p = 1 (mod 4)),
a certificate of this repository (symmetric/876 for q = 1753), or a published propus quadruple (orders 116, 236 and 604
for q = 233, 1889 and 2417).
"""
import hashlib, json, math, sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "verify"))
sys.path.insert(0, str(ROOT / "search"))
import verify  # noqa: E402
from make_yamada import primitive_poly  # noqa: E402

REPO_INPUTS = {876: "symmetric/876/certificate.json"}


def _sym(half, v):
    return sorted({x % v for x in half} | {-x % v for x in half})


# Published symmetric Hadamard matrices of orders 4 * 29, 4 * 59 and 4 * 151 as propus quadruples (A, B, B, D) over
# Z_v, A symmetric; the blocks list the positions of -1 (orbit representatives under the multipliers for v = 151).
SPORADIC = {
    116: {"source": "Di Matteo, Djokovic, Kotsireas, Spec. Matrices 3 (2015), arXiv:1503.04226v2, section 4, the first of "
                    "the fifteen SDSs for v = 29 (written [A, D, B, C = B])",
          "v": 29, "multipliers": [1], "subgroup_order": 1,
          "blocks": [[0, 4, 5, 6, 7, 9, 13, 16, 20, 22, 23, 24, 25], [0, 1, 2, 5, 6, 8, 11, 12, 14, 16, 22],
                     [0, 1, 2, 3, 5, 8, 10, 13, 14, 15, 18, 22, 25]]},
    236: {"source": "Djokovic, arXiv:2101.05429v2, section 4, the first PDF for v = 59, (59; 28, 29, 29, 22; 49) (s**): "
                    "X_1 = X_1' + (-X_1'), X_2 = X_3, X_4",
          "v": 59, "multipliers": [1], "subgroup_order": 1,
          "blocks": [_sym([1, 3, 5, 8, 10, 11, 13, 15, 16, 20, 21, 22, 26, 29], 59),
                     [0, 1, 3, 7, 8, 9, 10, 12, 13, 15, 16, 19, 21, 22, 25, 27, 29, 34, 35, 36, 37, 38, 39, 40, 44, 51,
                      54, 55, 58],
                     [0, 3, 4, 5, 7, 8, 14, 15, 16, 17, 19, 22, 24, 27, 28, 33, 39, 49, 53, 54, 55, 56]]},
    604: {"source": "Balonin, Djokovic, Inf. Control Syst. 2018 no. 4, arXiv:1803.08787v1, Table 1, "
                    "(151; 71, 71, 71, 66; 128), H = {1, 8, 19, 59, 64}",
          "v": 151, "multipliers": [8], "subgroup_order": 5,
          "blocks": [[0, 2, 5, 6, 7, 11, 15, 17, 23, 27, 30, 34, 37, 51, 68],
                     [0, 1, 2, 3, 4, 14, 17, 23, 27, 28, 34, 47, 51, 68, 87],
                     [0, 1, 2, 3, 4, 5, 7, 10, 29, 34, 46, 47, 51, 68]]},
}


def is_prime(n):
    return n > 1 and all(n % d for d in range(2, math.isqrt(n) + 1))


def prime_power(n):
    for p in range(2, math.isqrt(n) + 2):
        if n % p == 0:
            k, x = 0, n
            while x % p == 0:
                x //= p
                k += 1
            return (p, k) if x == 1 else None
    return (n, 1) if n > 1 else None


def undecided():
    p = ROOT / "search" / "data" / "symmetric_census_undecided.txt"
    return [int(x) for x in p.read_text().splitlines() if x.strip() and not x.startswith("#")]


def base_spec(m):
    """A symmetric Hadamard matrix of order m from Paley's constructions over a prime field, or None."""
    if m in REPO_INPUTS:
        cert = json.loads((ROOT / REPO_INPUTS[m]).read_text())
        return {"kind": "certificate", "path": REPO_INPUTS[m], "sha256": cert["sha256"]}
    if is_prime(m - 1) and (m - 1) % 4 == 3:
        return {"kind": "paley-amicable", "q": m - 1}
    if m % 2 == 0 and is_prime(m // 2 - 1) and (m // 2 - 1) % 4 == 1:
        return {"kind": "paley-ii", "q": m // 2 - 1}
    if m in SPORADIC:
        s = SPORADIC[m]
        return {"kind": "propus", **s, "rowsums": _rowsums(s)}
    return None


def _rowsums(s):
    v, gens = s["v"], s["multipliers"]
    H = {1}
    while True:
        grown = H | {h * g % v for h in H for g in gens}
        if grown == H:
            break
        H = grown
    return [v - 2 * len({r * h % v for r in reps for h in H}) for reps in s["blocks"]]


def input_spec(n):
    """Sylvester^a times one base matrix (fewest Sylvester factors first)."""
    a = 0
    while n % 2 == 0 or a == 0:
        m = n >> a
        b = base_spec(m) if m > 1 else None
        if b is not None:
            return b if a == 0 else {"kind": "kron", "factors": [{"kind": "sylvester", "k": a}, b]}
        if m % 2:
            break
        a += 1
    return None


def q_spec(q):
    p, k = prime_power(q)
    if k == 1:
        g = next(x for x in range(2, q) if all(pow(x, (q - 1) // r, q) != 1
                                               for r in range(2, q) if (q - 1) % r == 0 and is_prime(r)))
        return {"q": q, "g": g}
    return {"q": {"p": p, "poly": primitive_poly(p, k)}}


def targets():
    out = []
    for v in undecided():
        pk = prime_power(v)
        if pk and v % 8 == 1:
            s = input_spec((v - 1) // 2)
            out.append((v, s))
    return out


def certify(v):
    s = input_spec((v - 1) // 2)
    if s is None:
        raise SystemExit(f"no input for v = {v}")
    cert = {"order": 4 * v, "construction": "symmetric-propus-miyamoto", "symmetric": True, **q_spec(v), "input": s}
    H, _ = verify.BUILDERS[cert["construction"]](cert)
    n = H.shape[0]
    assert n == 4 * v and verify.is_pm1(H) and verify.gram_is_scalar(H, n) and (H == H.T).all(), "not symmetric Hadamard"
    sha = hashlib.sha256()
    for chunk in verify.csv_chunks(H):
        sha.update(chunk)
    cert["sha256"] = sha.hexdigest()
    d = ROOT / "symmetric" / str(4 * v)
    if (d / "certificate.json").exists():
        raise SystemExit(f"{d}/certificate.json exists; not overwritten")
    d.mkdir(exist_ok=True)
    (d / "certificate.json").write_text(json.dumps(cert, indent=1) + "\n")
    print(f"symmetric {4 * v} = 4 * {v}: H = H^T, H H^T = {n} I verified; input {json.dumps(s)}; sha256 {cert['sha256']}",
          flush=True)


if __name__ == "__main__":
    if sys.argv[1] == "targets":
        for v, s in targets():
            print(v, (v - 1) // 2, json.dumps(s) if s else "no input")
    elif sys.argv[1] == "certify":
        for v in map(int, sys.argv[2:]):
            certify(v)
