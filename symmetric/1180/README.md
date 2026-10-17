# Symmetric Hadamard matrix of order 1180

1180 = 4 · 295, with 295 = 5 · 59; 295 = q + 2 with q = 293 a prime ≡ 5 (mod 8) and (q + 3)/2 = 148.

## Construction

A Hadamard matrix H is symmetric if H = Hᵀ. This one is the case q = 293 of the construction in
[YAMADA_KIYASU.md](../YAMADA_KIYASU.md): Yamada's (Kiyasu's) Lemma 26 data put in quaternary form, which is symmetric for
any skew input, and Turyn's map to a real matrix.

- **The skew input.** A skew Goethals–Seidel matrix of order 148 over ℤ₃₇ found by `search/skew_exhaustive.py 37 10` (blocks invariant under a multiplier subgroup of order 3; skew-Hadamard matrices of order 148 are classical, e.g. Đoković's supplementary difference sets of 1992). Normalised to first row all +1, its core S of order 147 has
  Sᵀ = −S, Se = 0 and S² = J − 147I.
- **Yamada's circulants.** M and N of order 147, with first rows β and α of Yamada's Lemma 3 over
  GF(293²) = GF(293)[w]/(w² − 2) and ξ = 3 + w primitive: symmetric, M of zero diagonal,
  M² + N² = 293I.
- **The quaternary matrix.** With G = ½((M + N) + i(M − N)), P = [[1, −1], [−1, 1]] and e the all-one row of length 294,
  L = [[−1, −e, −1, ie], [−eᵀ, ½(1 + i)I ⊗ J₂ + G ⊗ P, −ieᵀ, −iI ⊗ P − S ⊗ J₂], [−1, −ie, 1, e],
  [ieᵀ, −iI ⊗ P + S ⊗ J₂, eᵀ, ½(i − 1)I ⊗ J₂ + Ḡ ⊗ P]] of order 590 is symmetric with entries ±1, ±i and
  L L̄ = 590I.
- **The matrix.** With L = X + iY, H = [[X + Y, X − Y], [X − Y, −X − Y]].

## Search

None: the matrix is determined by q, r, ξ and the skew input. `python3 search/sym_yk.py certify 295` writes the
certificate.

## Status

**New** as a symmetric Hadamard matrix.

- The Đoković–Seberry census (2025-05-28) lists 295 among its 595 undecided v for symmetric matrices. Seberry–Yamada's
  Corollary 27(i) in the propus array would need a symmetric conference matrix of order 148 or a skew-Hadamard
  matrix of order 148 with an involutive anti-automorphism ([NOTES.md](../../NOTES.md) §29).
- The symmetric quaternary form of Yamada's matrix is the new step ([YAMADA_KIYASU.md](../YAMADA_KIYASU.md#status)); we
  found it in none of the sources listed there. The skew input is classical.
- Zaitsev, Kharaghani, Mohammadian and Tayfeh-Rezaie's page of quaternary Hadamard matrices (last updated 14 October
  2026) lists QH590 under "unpublished approaches". It is a quaternary Yamada–Kiyasu matrix of order 590, and a
  monomial matrix takes it to a symmetric one, from which Turyn's map gives a symmetric Hadamard matrix of order 1180
  too ([YAMADA_KIYASU.md](../YAMADA_KIYASU.md#status)). The page does not mention symmetry or real matrices.
- Ordinary Hadamard matrices of order 1180 are classical (Yamada 1989, Theorem 4, from the same skew input); the
  symmetry is what is new.

The sources searched are listed in [NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, r, ξ, the complete certificate of the skew input, and the SHA-256 of the matrix |

The matrix itself is not stored; `verify/verify.py` rebuilds it.

## Verification

```
python3 verify/verify.py symmetric/1180
```

`verify/verify.py` (construction `symmetric-yamada-kiyasu`) rebuilds the skew input from its certificate and checks it
(Q + Qᵀ = 2I), computes Yamada's circulants and checks M² + N² = 293I, builds X and Y and checks that they are symmetric
with disjoint supports covering every position, and checks H = Hᵀ, H Hᵀ = 1180 I and the SHA-256 of the canonical CSV.

SHA-256 of `H1180.csv`: `9dd048b521f99d3ef5a268adf6da698fdeb331186b5a04cae6b3d58a8e404364`
