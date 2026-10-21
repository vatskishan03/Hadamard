# Symmetric Hadamard matrix of order 612

612 = 4 · 153, with v = 153 = 3² · 17 = q(q + 1)/2 for q = 17.

## Construction

A Hadamard matrix H is symmetric if H = Hᵀ. This one is the case n = 18 of [CONFERENCE_PLUG.md](../CONFERENCE_PLUG.md):
a symmetric Hadamard matrix of order 2n(n − 1) from a symmetric conference matrix of order n.

- **The conference matrix.** Paley's matrix C of order 18 over GF(17) = ℤ₁₇: index 0 is ∞, C[∞, x] = C[x, ∞] = 1,
  C[x, y] = χ(y − x) with χ the quadratic character, so C = [[0, eᵀ], [e, Q]] with Q the Paley core of order 17.
- **The plugs.** With J = J₁₇, X = [[J, J], [J, −J]] and Y = [[Q + I, Q − I], [Q − I, −Q − I]], of order 34. Both are
  symmetric, XY + YX = 0, and X² + 17Y² = 612 I.
- **The matrix.** H = I₁₈ ⊗ X + C ⊗ Y.

It is also Turyn's symmetric image of the complex-symmetric quaternary Hadamard matrix I ⊗ J + C ⊗ (Q + iI) of order
306 (Kharaghani–Seberry).

## Search

None: the matrix is determined by q. `python3 search/sym_conference_plug.py certify 17` writes the certificate.

## Status

**Known** as a symmetric Hadamard matrix.

- Spence (An infinite family of Williamson matrices, J. Austral. Math. Soc. Ser. A 24 (1977) 252–256),
  Corollary 3 with q = 17 and r = 1, gives a symmetric Hadamard matrix of order 2qʳ(q + 1) = 612. The matrix here
  comes from the construction of [CONFERENCE_PLUG.md](../CONFERENCE_PLUG.md), found independently; whether it is
  equivalent to Spence's was not checked.
- The Đoković–Seberry census (2025-05-28) lists v = 153 as undecided; its series omit Spence's corollary.

The sources searched are listed in [NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q (with the field polynomial for a prime power), the construction and the SHA-256 of the matrix |

The matrix itself is not stored; `verify/verify.py` rebuilds it.

## Verification

```
python3 verify/verify.py symmetric/612
```

`verify/verify.py` (construction `symmetric-conference-plug`) builds Paley's conference matrix, checks that it is a
normalised symmetric conference matrix, that Q𝟙 = 0 and XY + YX = 0, assembles H and checks H = Hᵀ, H Hᵀ = 612 I and
the SHA-256 of the canonical CSV.

SHA-256 of `H612.csv`: `6c676695f3d4b073596bc743fff629edbd73c5a900d11a50162d9f9cda14cdaf`
