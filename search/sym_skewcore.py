#!/usr/bin/env python3
"""Symmetric Hadamard matrices of order 4(q + 2) from a skew-propus quadruple of length l = (q + 3)/8 (NOTES.md §29).

A skew-propus quadruple (A3, A1, A1, A0) of odd length l (A3 symmetric, A0 of skew type, PAF_A3 + 2 PAF_A1 + PAF_A0 =
4l delta) gives, in the Balonin array, a symmetric orthogonal design xP + yN of type (1, 4l - 1), so S = I + NP is a
skew-Hadamard matrix of order 4l with the involutive anti-automorphism P (PSP = S^T), which fixes two points. Put n = 4l -
1 and normalise S at one fixed point: its core B (order n, skew, B J = 0) satisfies R B R = B^T for the reversal R once
the other points are labelled so that P acts as x -> -x on Z_n. Seberry-Yamada's Lemma 26 with U = (I, 0, BR, BR) and
V = (M, N, R, R), M and N Yamada's symmetric circulants of order n = (q + 1)/2 with M^2 + N^2 = qI (q = 2n - 1 a prime
= 1 mod 4), gives symmetric Williamson-type matrices of order q + 2 with W_3 = W_4, and the propus array a symmetric
Hadamard matrix of order 4(q + 2). This is Seberry-Yamada's Corollary 27(i) with a skew core in place of Paley's.

    python3 search/sym_skewcore.py check A3 A1 A0     build and check (+/- strings of length l)
    python3 search/sym_skewcore.py certify A3 A1 A0   also write symmetric/<4(q+2)>/certificate.json
"""
import hashlib, json, sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path[:0] = [str(ROOT / "verify"), str(ROOT / "search")]
import verify  # noqa: E402
from make_yamada import primitive  # noqa: E402


def seq(s):
    return np.array([1 if c == "+" else -1 for c in s], dtype=np.int64)


def certify(a3, a1, a0, write):
    l = len(a3)
    q = 8 * l - 3
    r = next(x for x in range(2, q) if pow(x, (q - 1) // 2, q) == q - 1)
    cert = {"order": 4 * (q + 2), "construction": "symmetric-propus-skewcore", "symmetric": True,
            "skew_propus": [a3, a1, a0], "q": q, "r": r, "xi": list(primitive(q, r))}
    H, _ = verify.BUILDERS[cert["construction"]](cert)
    m = H.shape[0]
    assert verify.is_pm1(H) and verify.gram_is_scalar(H, m) and (H == H.T).all(), "not symmetric Hadamard"
    sha = hashlib.sha256()
    for chunk in verify.csv_chunks(H):
        sha.update(chunk)
    cert["sha256"] = sha.hexdigest()
    print(f"symmetric {m} = 4 * {q + 2} (q = {q}, l = {l}): H = H^T, H H^T = {m} I verified; sha256 {cert['sha256']}",
          flush=True)
    if write:
        d = ROOT / "symmetric" / str(m)
        d.mkdir(exist_ok=True)
        (d / "certificate.json").write_text(json.dumps(cert, indent=1) + "\n")


if __name__ == "__main__":
    certify(*sys.argv[2:5], write=sys.argv[1] == "certify")
