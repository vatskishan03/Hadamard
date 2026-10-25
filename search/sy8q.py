#!/usr/bin/env python3
"""Certificates for Seberry-Yamada's H(8q) (verify.py construction `seberry-yamada-8q`): a Hadamard matrix of order 8q
for a prime q = 3 (mod 4) from symmetric Williamson-type matrices of order n = (q - 1)/2 (Seberry-Yamada 1992,
Lemma 8.12, Corollaries 8.14 and 8.18; search/NOTES.md section 95).

    python3 sy8q.py wt q RECIPE         certificate with the Williamson-type matrices packed in it; RECIPE is a
                                        search/symwt.py recipe, e.g. python3 sy8q.py wt 1087 "cor27i(541)"
    python3 sy8q.py c29 q INPUT         certificate whose Williamson-type matrices verify.py rebuilds by Corollary 29(ii)
                                        from the symmetric Hadamard matrix of order (n - 1)/2 in INPUT (a certificate
                                        path relative to the repository), e.g. python3 sy8q.py c29 1619
                                        symmetric/404/certificate.json

The certificate is printed (JSON) after verify.py has built the matrix and checked H H^T = 8q I; the primitive roots
are the least ones (sympy). Run from this directory."""
import base64
import hashlib
import json
import sys
from pathlib import Path

import numpy as np
from sympy import primitive_root

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "verify"))
sys.path.insert(0, str(ROOT / "search"))
import verify as V  # noqa: E402


def pack_symmetric(X):
    """inverse of verify.unpack_symmetric: the upper triangle, row by row, 1 for +1, np.packbits, base64"""
    X = np.asarray(X)
    return base64.b64encode(np.packbits((X[np.triu_indices(X.shape[0])] > 0).astype(np.uint8)).tobytes()).decode()


def finish(cert):
    H, _ = V.BUILDERS[cert["construction"]](cert)
    if not (V.is_pm1(H) and V.gram_is_scalar(H, H.shape[0])):
        raise SystemExit("not a Hadamard matrix")
    sha = hashlib.sha256()
    for chunk in V.csv_chunks(H):
        sha.update(chunk)
    cert["sha256"] = sha.hexdigest()
    return cert


def main():
    mode, q = sys.argv[1], int(sys.argv[2])
    cert = {"order": 8 * q, "construction": "seberry-yamada-8q", "q": q, "g": int(primitive_root(q))}
    if mode == "wt":
        import symwt
        recipe = sys.argv[3]
        cert["wt_source"] = f'search/symwt.py: build("{recipe}")'
        cert["wt"] = [pack_symmetric(w) for w in symwt.build(recipe)]
    elif mode == "c29":
        inp = sys.argv[3]
        n = (q - 1) // 2
        cert["wt_corollary29"] = {"q": n, "g": int(primitive_root(n)), "input": inp,
                                  "input_sha256": json.loads((ROOT / inp).read_text())["sha256"]}
    else:
        raise SystemExit(__doc__)
    print(json.dumps(finish(cert), indent=1))


if __name__ == "__main__":
    main()
