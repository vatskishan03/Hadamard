#!/usr/bin/env python3
"""Symmetric Hadamard matrices of order 4(l + 1) from a quaternary Legendre pair (A, B) of even length l with A real
and B symmetric (B(-g) = B(g)): decoding B(g) = ((1 + i) p(g) + (1 - i) q(g))/2 with p, q = +-1 gives a bordered
propus quadruple (p, A, A, q) over Z_l, since PAF_A + PAF_B = -2 off 0 becomes PAF_p + 2 PAF_A + PAF_q = -4, and p, q
are symmetric because B is (verify.py construction `propus-bordered`, search/NOTES.md section 27).

    python3 search/sym_qlp.py check            decode the stored pairs and check the identities
    python3 search/sym_qlp.py certify 65       writes symmetric/260/certificate.json and sequences.txt

PAIRS holds published pairs as digit strings (digit e stands for i^e).
"""
import hashlib, json, sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "verify"))
import verify  # noqa: E402

PAIRS = {64: {"source": "H. Verma, A quaternary Legendre pair of length 64, arXiv:2609.25573 (22 September 2026), Table 1",
              "A": "02022220 02022222 02200020 00202220 20220002 02002000 00220022 00022022",
              "B": "02332221 10130032 10300312 13202110 30112023 12130030 12300310 11222332"}}


def decode(l):
    d = PAIRS[l]
    A = np.array([1j ** int(c) for c in d["A"].replace(" ", "")])
    B = np.array([1j ** int(c) for c in d["B"].replace(" ", "")])
    assert len(A) == len(B) == l
    paf = lambda x, s: complex(np.sum(x * np.conj(np.roll(x, -s))))
    assert all(abs(paf(A, s) + paf(B, s) + 2) < 1e-9 for s in range(1, l)), "not a quaternary Legendre pair"
    assert np.allclose(A.imag, 0), "A is not real"
    t = np.arange(l)
    assert np.allclose(B, B[-t % l]), "B is not symmetric"
    a = np.rint(A.real).astype(np.int64)
    p = np.rint((B.real + B.imag)).astype(np.int64)            # B = ((1 + i) p + (1 - i) q)/2
    q = np.rint((B.real - B.imag)).astype(np.int64)
    assert set(np.unique(p)) <= {-1, 1} and set(np.unique(q)) <= {-1, 1}
    tot = verify.paf_sum([p, a, a, q])
    assert tot[0] == 4 * l and np.all(tot[1:] == -4), "not a bordered propus quadruple"
    return p, a, q


def certify(v):
    l = v - 1
    p, a, q = decode(l)
    s = lambda x: "".join("+" if y > 0 else "-" for y in x)
    cert = {"order": 4 * v, "construction": "propus-bordered", "symmetric": True, "group": [l],
            "sequences": [s(p), s(a), s(q)], "rowsums": [int(p.sum()), int(a.sum()), int(q.sum())],
            "K": [[1, 1, 1, 1], [1, 1, -1, -1], [1, -1, 1, -1], [1, -1, -1, 1]], "source": PAIRS[l]["source"]}
    H, _ = verify.BUILDERS[cert["construction"]](cert)
    n = H.shape[0]
    assert n == 4 * v and verify.is_pm1(H) and verify.gram_is_scalar(H, n) and (H == H.T).all()
    cert["sha256"] = hashlib.sha256(verify.canonical_csv(H)).hexdigest()
    d = ROOT / "symmetric" / str(4 * v)
    if (d / "certificate.json").exists():
        raise SystemExit(f"{d} exists")
    d.mkdir(exist_ok=True)
    (d / "certificate.json").write_text(json.dumps(cert, indent=1) + "\n")
    (d / "sequences.txt").write_text("".join(x + "\n" for x in cert["sequences"]))
    print(f"symmetric {4 * v} = 4 * {v}: H = H^T, H H^T = {n} I verified; sha256 {cert['sha256']}")


if __name__ == "__main__":
    if sys.argv[1:2] == ["check"]:
        for l in PAIRS:
            p, a, q = decode(l)
            print(f"length {l}: row sums {int(p.sum())}, {int(a.sum())}, {int(q.sum())}; PAF_p + 2 PAF_a + PAF_q = -4 off 0")
    elif sys.argv[1:2] == ["certify"]:
        for v in sys.argv[2:]:
            certify(int(v))
    else:
        raise SystemExit(__doc__)
