# Symmetric Hadamard matrix of order 7012

7012 = 4 · 1753, with 1753 a prime ≡ 1 (mod 8) and (q − 1)/2 = 876.

## Construction

A Hadamard matrix H is symmetric if H = Hᵀ. This one is the case q = 1753 of the construction in
[PROPUS_MIYAMOTO.md](../PROPUS_MIYAMOTO.md): Seberry and Yamada's Williamson-type matrices W₁, W₂, W₃ = W₄ of order q
(their Corollary 29(ii), 1990), placed in the propus array.

- **The input.** A symmetric Hadamard matrix M of order 876: our symmetric Hadamard matrix of order 876
  ([symmetric/876](../876), found in 2026).
- **The field.** GF(1753) = ℤ_1753, g = 7; χ is the quadratic character, s_i = g²ⁱ and t_j = g¹⁻²ʲ (i, j < 876).
- **Williamson-type matrices.** U₁ = −[χ(s_i − s_j)], U₂ = [χ(s_i − t_j)], U₃ = U₄ = 0 and V = (I, 0, M, M) give W₁ =
  [[1, −e], [−eᵀ, S₁]] and W_k = [[1, e], [eᵀ, S_k]] (k > 1), with S_k = U_k ⊗ J₂ + V_k ⊗ [[1, −1], [−1, 1]]. They are
  symmetric and pairwise amicable, Σ W_k W_kᵀ = 7012 I, and W₃ = W₄.
- **The matrix.** H = [[W₁, W₃, W₃, W₂], [W₃, W₂, −W₁, −W₃], [W₃, −W₁, −W₂, W₃], [W₂, −W₃, W₃, −W₁]].

## Search

None: the matrix is determined by q, g and M.

**Reproduction.**

| Command | Output |
|---|---|
| `python3 search/sym_propus_miyamoto.py certify 1753` | `symmetric/7012/certificate.json`, identical to the stored file |
| `python3 verify/verify.py symmetric/7012` | `ok`, SHA-256 `ba72ed96d666f2a76f7f57a78db6ecb6e275ca8504fdd9bf9736f167c36c0ccb` |

## Status

**New** as a symmetric Hadamard matrix ([PROPUS_MIYAMOTO.md](../PROPUS_MIYAMOTO.md#status)). The construction is
published, but its input, our symmetric matrix of order 876, is from 2026.

- The Đoković–Seberry census (2025-05-28) lists 1753 among its 595 undecided v for symmetric matrices; its
  Seberry–Balonin series uses the propus variation only for v = q + 2.
- Seberry and Yamada list 1753 among the orders for which Corollary 29(ii) would give symmetric Williamson-type
  matrices "when symmetric Hadamard matrices of orders … 4.219 … are discovered" (J. Combin. Math. Combin. Comput. 7
  (1990), Remark 30). A symmetric Hadamard matrix also needs the propus array (2015).
- Miyamoto (J. Combin. Theory Ser. A 57 (1991), Theorem 5) claimed Hadamard matrices of order 4q "of symmetric block
  form" for every q ≡ 1 (mod 4) with a conference matrix of order q + 1 and a symmetric Hadamard matrix of order
  q − 1, which would include this order; the census says his proof does not give it
  ([PROPUS_MIYAMOTO.md](../PROPUS_MIYAMOTO.md#status)).
- Ordinary Hadamard matrices of order 7012 are classical (Miyamoto's theorem, from a Hadamard matrix of order q − 1 =
  1752); the symmetry is what is new.

The sources searched are listed in [NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q (and the field), g, the description of the input M, and the SHA-256 of the matrix |

The matrix itself is not stored; `verify/verify.py` rebuilds it.

## Verification

```
python3 verify/verify.py symmetric/7012
```

`verify/verify.py` (construction `symmetric-propus-miyamoto`) rebuilds M and checks that it is a symmetric Hadamard
matrix, builds U and V, checks that the W_k are ±1, symmetric, pairwise amicable, with Σ W_k W_kᵀ = 7012 I and
W₃ = W₄, assembles H and checks H = Hᵀ, H Hᵀ = 7012 I and the SHA-256 of the canonical CSV.

SHA-256 of `H7012.csv`: `ba72ed96d666f2a76f7f57a78db6ecb6e275ca8504fdd9bf9736f167c36c0ccb`
