# Symmetric Hadamard matrix of order 5156

5156 = 4 · 1289, with 1289 a prime ≡ 1 (mod 8) and (q − 1)/2 = 644.

## Construction

A Hadamard matrix H is symmetric if H = Hᵀ. This one is the case q = 1289 of the construction in
[PROPUS_MIYAMOTO.md](../PROPUS_MIYAMOTO.md): Seberry and Yamada's Williamson-type matrices W₁, W₂, W₃ = W₄ of order q
(their Corollary 29(ii), 1990), placed in the propus array.

- **The input.** A symmetric Hadamard matrix M of order 644: the symmetric member P + N of Paley's amicable pair of
  order 644 (q′ = 643 ≡ 3 mod 4).
- **The field.** GF(1289) = ℤ_1289, g = 6; χ is the quadratic character, s_i = g²ⁱ and t_j = g¹⁻²ʲ (i, j < 644).
- **Williamson-type matrices.** U₁ = −[χ(s_i − s_j)], U₂ = [χ(s_i − t_j)], U₃ = U₄ = 0 and V = (I, 0, M, M) give W₁ =
  [[1, −e], [−eᵀ, S₁]] and W_k = [[1, e], [eᵀ, S_k]] (k > 1), with S_k = U_k ⊗ J₂ + V_k ⊗ [[1, −1], [−1, 1]]. They are
  symmetric and pairwise amicable, Σ W_k W_kᵀ = 5156 I, and W₃ = W₄.
- **The matrix.** H = [[W₁, W₃, W₃, W₂], [W₃, W₂, −W₁, −W₃], [W₃, −W₁, −W₂, W₃], [W₂, −W₃, W₃, −W₁]].

## Search

None: the matrix is determined by q, g and M. `python3 search/sym_propus_miyamoto.py certify 1289` writes the
certificate.

## Status

**New** as a symmetric Hadamard matrix, but obtainable classically ([PROPUS_MIYAMOTO.md](../PROPUS_MIYAMOTO.md#status)):
Seberry and Yamada's Williamson-type matrices with a repeated block (1990), the propus array (Seberry and Balonin,
arXiv:1512.01732) and the input were all published, and we found them combined nowhere.

- The Đoković–Seberry census (2025-05-28) lists 1289 among its 595 undecided v for symmetric matrices; its
  Seberry–Balonin series uses the propus variation only for v = q + 2.
- Ordinary Hadamard matrices of order 5156 are classical (Miyamoto's theorem, from a Hadamard matrix of order q − 1 =
  1288); the symmetry is what is new.

The sources searched are listed in [NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q (and the field), g, the description of the input M, and the SHA-256 of the matrix |

The matrix itself is not stored; `verify/verify.py` rebuilds it.

## Verification

```
python3 verify/verify.py symmetric/5156
```

`verify/verify.py` (construction `symmetric-propus-miyamoto`) rebuilds M and checks that it is a symmetric Hadamard
matrix, builds U and V, checks that the W_k are ±1, symmetric, pairwise amicable, with Σ W_k W_kᵀ = 5156 I and W₃ = W₄,
assembles H and checks H = Hᵀ, H Hᵀ = 5156 I and the SHA-256 of the canonical CSV.

SHA-256 of `H5156.csv`: `150ed925aaebb008c0e0274dcada35d8170b0322ce2a3ec43cf53f5ad87cf188`
