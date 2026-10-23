#!/usr/bin/env python3
"""Symmetric Hadamard matrices of order 4(q + 2), q = 5 (mod 8) a prime power, from ANY skew-Hadamard matrix of order
(q + 3)/2: Yamada's (Kiyasu's) matrix in quaternary form, made symmetric (verify.py construction
`symmetric-yamada-kiyasu`; NOTES.md section 31).

    python3 search/sym_yk.py targets          census-undecided v = q + 2 and the skew input we use
    python3 search/sym_yk.py certify v [...]  writes symmetric/<4v>/certificate.json

Skew inputs, all checked by verify.py inside the certificate: Paley's I + C of order p + 1 (p = 3 mod 4), the
Goethals-Seidel quadruples of data/yk_skew_inputs.json (published data as stored in SageMath, or found by
search/skew_exhaustive.py), the two-block matrices of Szekeres's (q = 5 mod 8 prime) and Whiteman's (q = p^2,
p = 5 mod 8) complementary difference sets, certificates of skew/, and doublings [[S, S], [-S^T, S^T]].
"""
import hashlib, json, sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "verify"))
sys.path.insert(0, str(ROOT / "search"))
import verify  # noqa: E402
import make_yamada as MY  # noqa: E402
import sym_propus_miyamoto as SPM  # noqa: E402

# v: skew input of order (v + 1)/2, as (kind, parameter, doublings)
RECIPES = {295: ("gs", 37, 0), 319: ("paley", 19, 3), 391: ("gs", 49, 0), 463: ("gs", 29, 1), 655: ("paley", 163, 1),
           679: ("cds", 169, 0), 703: ("paley", 43, 3), 823: ("gs", 103, 0), 943: ("gs", 59, 1),
           1063: ("gs", 133, 0), 1183: ("gs", 37, 2), 1231: ("paley", 307, 1), 1303: ("gs", 163, 0),
           1495: ("cds", 373, 0), 1599: ("paley", 199, 2), 1639: ("file", "skew/820/certificate.json", 0),
           1671: ("file", "skew/836/certificate.json", 0), 1903: ("file", "skew/476/certificate.json", 1),
           1951: ("paley", 487, 1), 1975: ("gs", 247, 0),
           1999: ("paley", 499, 1), 2071: ("file", "skew/1036/certificate.json", 0), 2143: ("gs", 67, 2),
           2223: ("cds", 277, 1), 2239: ("paley", 139, 3), 2271: ("paley", 283, 2),
           2295: ("file", "skew/1148/certificate.json", 0), 2335: ("gs", 73, 2)}


def _digest(H):
    sha = hashlib.sha256()
    for chunk in verify.csv_chunks(H):
        sha.update(chunk)
    return sha.hexdigest()


def _seal(cert):
    cert["sha256"] = _digest(verify.BUILDERS[cert["construction"]](cert)[0])
    verify.build_checked(cert)
    return cert


def gs_cert(w):
    """Goethals-Seidel skew certificate of order 4w from data/yk_skew_inputs.json (first row of skew type)."""
    d = json.loads((ROOT / "search" / "data" / "yk_skew_inputs.json").read_text())[str(w)]
    seqs = [[1 if c == "+" else -1 for c in s] for s in d["seqs"]]
    k = next(i for i, x in enumerate(seqs) if x[0] == 1 and all(x[-t % w] == -x[t] for t in range(1, w)))
    seqs = [seqs[k]] + seqs[:k] + seqs[k + 1:]
    return _seal({"order": 4 * w, "construction": "goethals-seidel", "skew": True, "group": f"Z_{w}", "v": w,
                  "multiplier": 1, "subgroup_order": 1, "rowsums": [sum(x) for x in seqs],
                  "blocks": [[t for t in range(w) if x[t] < 0] for x in seqs], "source": d["source"]})


def cds_cert(q):
    """Two-block skew certificate of order 2(q + 1) from complementary difference sets in GF(q): Szekeres (q = p = 5
    mod 8, A = C_0 + C_1, B = C_0 + C_3 over the fourth powers) or Whiteman (q = p^2, p = 5 mod 8, A = C_0 + ... + C_3,
    B = C_0 + C_1 + C_6 + C_7 over the eighth powers); C_i = rho^i C_0."""
    p, k = MY.prime_power(q)
    assert p % 8 == 5 and k in (1, 2)
    poly = MY.primitive_poly(p, k)                            # x^k = -(c_0 + ... + c_{k-1} x^{k-1}), x primitive
    def xmul(c):
        c = [0] + list(c)
        top = c.pop()
        return tuple((c[j] - top * poly[j]) % p for j in range(k))
    pw, c = [], tuple([1] + [0] * (k - 1))
    for _ in range(q - 1):
        pw.append(c)
        c = xmul(c)
    e = 4 if k == 1 else 8
    cls = lambda i: {pw[j] for j in range(i, q - 1, e)}
    A = cls(0) | cls(1) | ((cls(2) | cls(3)) if k == 2 else set())
    B = cls(0) | (cls(3) if k == 1 else cls(1) | cls(6) | cls(7))
    elts = [tuple(int(x) for x in t) for t in __import__("itertools").product(range(p), repeat=k)]  # mixed radix
    return _seal({"order": 2 * (q + 1), "construction": "two-block", "kind": "skew", "skew": True,
                  "group": f"{p}^{k}" if k > 1 else str(p), "a": [-1 if t in A else 1 for t in elts],
                  "b": [-1 if t in B else 1 for t in elts],
                  "source": "Szekeres's complementary difference sets (1969)" if k == 1 else
                            "Whiteman's complementary difference sets (1971)"})


def skew_cert(recipe):
    kind, x, doublings = recipe
    if kind == "paley":
        Q = _seal({"order": x + 1, "construction": "paley-skew", "skew": True, "q": x})
    elif kind == "gs":
        Q = gs_cert(x)
    elif kind == "cds":
        Q = cds_cert(x)
    else:
        Q = json.loads((ROOT / x).read_text())
    for _ in range(doublings):
        Q = _seal({"order": 2 * Q["order"], "construction": "skew-doubling", "skew": True, "S": Q})
    return Q


def targets():
    done = {int(p.name) // 4 for p in (ROOT / "symmetric").iterdir() if p.name.isdigit()}
    out = []
    for v in SPM.undecided():
        pk = SPM.prime_power(v - 2)
        if pk and (v - 2) % 8 == 5 and v not in done:
            out.append((v, RECIPES.get(v)))
    return out


def certify(v):
    q = v - 2
    p, k = MY.prime_power(q)
    assert q % 8 == 5 and v in RECIPES, f"no recipe for v = {v}"
    Q = skew_cert(RECIPES[v])
    assert Q["order"] == (q + 3) // 2
    if k == 1:
        r = next(x for x in range(2, q) if pow(x, (q - 1) // 2, q) == q - 1)
        cert = {"order": 4 * v, "construction": "symmetric-yamada-kiyasu", "symmetric": True, "q": q, "r": r,
                "xi": MY.primitive(q, r)}
    else:
        poly = MY.primitive_poly(p, k)
        r, xi = MY.field_params(q, p, poly)
        cert = {"order": 4 * v, "construction": "symmetric-yamada-kiyasu", "symmetric": True, "q": q, "p": p,
                "poly": poly, "r": r, "xi": xi}
    cert["Q"] = Q
    H, _ = verify.BUILDERS[cert["construction"]](cert)
    m = H.shape[0]
    assert m == 4 * v and verify.is_pm1(H) and verify.gram_is_scalar(H, m) and (H == H.T).all(), "not symmetric Hadamard"
    cert["sha256"] = _digest(H)
    d = ROOT / "symmetric" / str(4 * v)
    if (d / "certificate.json").exists():
        raise SystemExit(f"{d} exists")
    d.mkdir(exist_ok=True)
    (d / "certificate.json").write_text(json.dumps(cert, indent=1) + "\n")
    print(f"symmetric {4 * v} = 4 * {v}: H = H^T, H H^T = {m} I verified; skew input {RECIPES[v]}; "
          f"sha256 {cert['sha256']}", flush=True)


if __name__ == "__main__":
    if sys.argv[1:2] == ["targets"]:
        for v, r in targets():
            print(v, 4 * v, (v + 1) // 2, r if r else "no skew-Hadamard matrix of this order in hand")
    elif sys.argv[1:2] == ["certify"]:
        for v in sys.argv[2:]:
            certify(int(v))
    else:
        raise SystemExit(__doc__)
