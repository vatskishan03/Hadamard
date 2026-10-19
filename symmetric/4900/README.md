# Symmetric Hadamard matrix of order 4900

4900 = 4 · 1225, with v = 1225 = 5² · 7² = q(q + 1)/2 for q = 49.

## Construction

A Hadamard matrix H is symmetric if H = Hᵀ. This one is the case n = 50 of [CONFERENCE_PLUG.md](../CONFERENCE_PLUG.md):
a symmetric Hadamard matrix of order 2n(n − 1) from a symmetric conference matrix of order n.

- **The conference matrix.** Paley's matrix C of order 50 over GF(49) = GF(7)[x]/(x² + x + 3): index 0 is ∞, C[∞, x] = C[x, ∞] = 1,
  C[x, y] = χ(y − x) with χ the quadratic character, so C = [[0, eᵀ], [e, Q]] with Q the Paley core of order 49.
- **The plugs.** With J = J₄₉, X = [[J, J], [J, −J]] and Y = [[Q + I, Q − I], [Q − I, −Q − I]], of order 98. Both are
  symmetric, XY + YX = 0, and X² + 49Y² = 4900 I.
- **The matrix.** H = I₅₀ ⊗ X + C ⊗ Y.

It is also Turyn's symmetric image of the complex-symmetric quaternary Hadamard matrix I ⊗ J + C ⊗ (Q + iI) of order
2450 (Kharaghani–Seberry).

## Search

None: the matrix is determined by q. `python3 search/sym_conference_plug.py certify 49` writes the certificate.

## Status

**New, but obtainable classically** as a symmetric Hadamard matrix.

- The Đoković–Seberry census (2025-05-28) lists v = 1225 among its undecided values for symmetric matrices.
- The ingredients are classical: Paley's conference matrix, Kharaghani and Seberry's complex Hadamard matrix of order
  n(n − 1) (Seberry–Yamada 1992, Theorem 11.14, Corollary 11.15) and Turyn's map, which is symmetric for a
  complex-symmetric input (Lemma 11.3 there). That the complex matrix is complex-symmetric, and so gives a symmetric
  real matrix, is stated in none of the sources we read ([CONFERENCE_PLUG.md](../CONFERENCE_PLUG.md#status)).
- Ordinary Hadamard matrices of order 4900 are classical. The symmetry is what is new.

The sources searched are listed in [NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q (with the field polynomial for a prime power), the construction and the SHA-256 of the matrix |

The matrix itself is not stored; `verify/verify.py` rebuilds it.

## Verification

```
python3 verify/verify.py symmetric/4900
```

`verify/verify.py` (construction `symmetric-conference-plug`) builds Paley's conference matrix, checks that it is a
normalised symmetric conference matrix, that Q𝟙 = 0 and XY + YX = 0, assembles H and checks H = Hᵀ, H Hᵀ = 4900 I and
the SHA-256 of the canonical CSV.

SHA-256 of `H4900.csv`: `f155fa50f4b5024311005c9eeb3f10bc9dac3150cd3b9621a3fea98e26ec1a23`
