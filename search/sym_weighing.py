#!/usr/bin/env python3
"""Symmetric Hadamard matrices of order 4q, q = 1 (mod 8) a prime power, from ANY Hadamard matrix K of order
n = (q - 1)/2 (verify.py construction `symmetric-propus-weighing`; NOTES.md section 30): the U-part of Seberry and
Yamada's Corollary 29(ii) in the propus array, with the V-part replaced by a symmetric weighing matrix of order 4n and
weight q built from K.

    python3 search/sym_weighing.py targets          census-undecided v = q not reached with a symmetric input
    python3 search/sym_weighing.py certify v [...]  writes symmetric/<4v>/certificate.json

K is given by a complete certificate: Cooper-Wallis (T-matrices of t-matrices/T<t>.txt, or of order 5 below, and
Williamson matrices of order w = (r + 1)/2 from Turyn's construction over GF(r), r = 2w - 1 a prime power = 1 (mod 4),
through Yamada's circulants search/symwt.yamada_pair), Miyamoto's construction over such a K, or a certificate file
(data/h716_alpoge.json: the matrix of order 716 of Alpöge, Voinov and Reynolds-Haertle, 2026, as a bordered
Goethals-Seidel array).
"""
import hashlib, json, sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "verify"))
sys.path.insert(0, str(ROOT / "search"))
import verify  # noqa: E402
import sym_propus_miyamoto as SPM  # noqa: E402

T5 = ["++000", "00+-0", "0000+", "00000"]        # T-matrices of order 5 (classical; found by exhaustive search)

# n = (q - 1)/2: ("cw", t, w) | ("miyamoto", p, recipe of order p - 1) | ("file", path)
RECIPES = {404: ("cw", 101, 1), 428: ("cw", 107, 1), 476: ("cw", 17, 7), 596: ("miyamoto", 149, ("cw", 37, 1)),
           716: ("file", "search/data/h716_alpoge.json"), 956: ("file", "search/data/skew956_dk.json"), 996: ("file", "skew/996/certificate.json"),
           1076: ("miyamoto", 269, ("cw", 67, 1)), 1196: ("cw", 23, 13), 1220: ("cw", 5, 61)}


def _signs(row):
    return "".join("+" if x > 0 else "-" for x in row)


def _digest(H):
    sha = hashlib.sha256()
    for chunk in verify.csv_chunks(H):
        sha.update(chunk)
    return sha.hexdigest()


def williamson(w):
    """first rows of Williamson matrices I + B, I - B, A, A of order w (Turyn; A, B = Yamada's circulants over
    GF(2w - 1), B of zero diagonal, A^2 + B^2 = (2w - 1) I)"""
    if w == 1:
        return ["+"] * 4
    import symwt
    al, be = symwt.yamada_pair(2 * w - 1)
    e = [1] + [0] * (w - 1)
    return [_signs([x + y for x, y in zip(e, be)]), _signs([x - y for x, y in zip(e, be)]), _signs(al), _signs(al)]


def k_cert(recipe):
    kind = recipe[0]
    if kind == "file":
        return json.loads((ROOT / recipe[1]).read_text())
    if kind == "cw":
        t, w = recipe[1], recipe[2]
        T = T5 if t == 5 else [ln.strip() for ln in (ROOT / "t-matrices" / f"T{t}.txt").read_text().splitlines()
                              if ln.strip()]
        cert = {"order": 4 * t * w, "construction": "cooper-wallis", "t": t, "w": w, "t_matrices": T,
                "williamson": williamson(w)}
    else:
        p, inner = recipe[1], recipe[2]
        cert = {"order": 4 * p, "construction": "miyamoto", "q": p, "conference_matrix": "paley", "K": k_cert(inner)}
    cert["sha256"] = _digest(verify.BUILDERS[cert["construction"]](cert)[0])
    verify.build_checked(cert)
    return cert


def targets():
    """census-undecided v = q = 1 (mod 8), a prime power, without a symmetric input of order (q - 1)/2"""
    return [(v, RECIPES.get((v - 1) // 2)) for v, s in SPM.targets() if s is None]


def certify(v):
    n = (v - 1) // 2
    if n not in RECIPES:
        raise SystemExit(f"no Hadamard matrix of order {n} in RECIPES")
    cert = {"order": 4 * v, "construction": "symmetric-propus-weighing", "symmetric": True, **SPM.q_spec(v),
            "K": k_cert(RECIPES[n])}
    H, _ = verify.BUILDERS[cert["construction"]](cert)
    m = H.shape[0]
    assert m == 4 * v and verify.is_pm1(H) and verify.gram_is_scalar(H, m) and (H == H.T).all(), "not symmetric Hadamard"
    cert["sha256"] = _digest(H)
    d = ROOT / "symmetric" / str(4 * v)
    if (d / "certificate.json").exists():
        raise SystemExit(f"{d} exists")
    d.mkdir(exist_ok=True)
    (d / "certificate.json").write_text(json.dumps(cert, indent=1) + "\n")
    print(f"symmetric {4 * v} = 4 * {v}: H = H^T, H H^T = {m} I verified; K {RECIPES[n]}; sha256 {cert['sha256']}",
          flush=True)


if __name__ == "__main__":
    if sys.argv[1:2] == ["targets"]:
        for v, r in targets():
            print(v, 4 * v, (v - 1) // 2, r if r else "no Hadamard matrix of this order known to us")
    elif sys.argv[1:2] == ["certify"]:
        for v in sys.argv[2:]:
            certify(int(v))
    else:
        raise SystemExit(__doc__)
