# Symmetric Hadamard matrix of order 1300

1300 = 4 · 325, with v = 325 = 5² · 13 = q(q + 1)/2 for q = 25.

## Construction

A Hadamard matrix H is symmetric if H = Hᵀ. This one is the case n = 26 of [CONFERENCE_PLUG.md](../CONFERENCE_PLUG.md):
a symmetric Hadamard matrix of order 2n(n − 1) from a symmetric conference matrix of order n.

- **The conference matrix.** Paley's matrix C of order 26 over GF(25) = GF(5)[x]/(x² + x + 2): index 0 is ∞, C[∞, x] = C[x, ∞] = 1,
  C[x, y] = χ(y − x) with χ the quadratic character, so C = [[0, eᵀ], [e, Q]] with Q the Paley core of order 25.
- **The plugs.** With J = J₂₅, X = [[J, J], [J, −J]] and Y = [[Q + I, Q − I], [Q − I, −Q − I]], of order 50. Both are
  symmetric, XY + YX = 0, and X² + 25Y² = 1300 I.
- **The matrix.** H = I₂₆ ⊗ X + C ⊗ Y.

It is also Turyn's symmetric image of the complex-symmetric quaternary Hadamard matrix I ⊗ J + C ⊗ (Q + iI) of order
650 (Kharaghani–Seberry).

## Search

None: the matrix is determined by q. `python3 search/sym_conference_plug.py certify 25` writes the certificate.

## Status

**Known** as a symmetric Hadamard matrix.

- Spence (An infinite family of Williamson matrices, J. Austral. Math. Soc. Ser. A 24 (1977) 252–256),
  Corollary 3 with q = 25 and r = 1, gives a symmetric Hadamard matrix of order 2qʳ(q + 1) = 1300. The matrix here
  comes from the construction of [CONFERENCE_PLUG.md](../CONFERENCE_PLUG.md), found independently; whether it is
  equivalent to Spence's was not checked.
- The Đoković–Seberry census (2025-05-28) lists v = 325 as undecided; its series omit Spence's corollary.

The sources searched are listed in [NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q (with the field polynomial for a prime power), the construction and the SHA-256 of the matrix |

The matrix itself is not stored; `verify/verify.py` rebuilds it.

## Verification

```
python3 verify/verify.py symmetric/1300
```

`verify/verify.py` (construction `symmetric-conference-plug`) builds Paley's conference matrix, checks that it is a
normalised symmetric conference matrix, that Q𝟙 = 0 and XY + YX = 0, assembles H and checks H = Hᵀ, H Hᵀ = 1300 I and
the SHA-256 of the canonical CSV.

SHA-256 of `H1300.csv`: `2c0f0eb35f3c5e640dca5312daa3477e62755821edced88839f5806fc5732740`
