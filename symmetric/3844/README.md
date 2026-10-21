# Symmetric Hadamard matrix of order 3844

3844 = 4 · 961, with 961 = 31² ≡ 1 (mod 8) and (q − 1)/2 = 480.

## Construction

A Hadamard matrix H is symmetric if H = Hᵀ. This one is the case q = 961 of the construction in
[PROPUS_MIYAMOTO.md](../PROPUS_MIYAMOTO.md): Seberry and Yamada's Williamson-type matrices W₁, W₂, W₃ = W₄ of order q
(their Corollary 29(ii), 1990), placed in the propus array.

- **The input.** A symmetric Hadamard matrix M of order 480: the symmetric member P + N of Paley's amicable pair of
  order 480 (q′ = 479 ≡ 3 mod 4).
- **The field.** GF(961) = GF(31)[x]/(x² + 2x + 3), with g = x primitive; χ is the quadratic character, s_i = g²ⁱ and
  t_j = g¹⁻²ʲ (i, j < 480).
- **Williamson-type matrices.** U₁ = −[χ(s_i − s_j)], U₂ = [χ(s_i − t_j)], U₃ = U₄ = 0 and V = (I, 0, M, M) give W₁ =
  [[1, −e], [−eᵀ, S₁]] and W_k = [[1, e], [eᵀ, S_k]] (k > 1), with S_k = U_k ⊗ J₂ + V_k ⊗ [[1, −1], [−1, 1]]. They are
  symmetric and pairwise amicable, Σ W_k W_kᵀ = 3844 I, and W₃ = W₄.
- **The matrix.** H = [[W₁, W₃, W₃, W₂], [W₃, W₂, −W₁, −W₃], [W₃, −W₁, −W₂, W₃], [W₂, −W₃, W₃, −W₁]].

## Search

None: the matrix is determined by q, g and M.

**Reproduction.**

| Command | Output |
|---|---|
| `python3 search/sym_propus_miyamoto.py certify 961` | `symmetric/3844/certificate.json`, identical to the stored file |
| `python3 verify/verify.py symmetric/3844` | `ok`, SHA-256 `b34c36b2b2b98c4b1ab9d8a68c4951a2a4b11651f8533b6ffe66d87bc1efffb5` |

## Status

**New, but obtainable classically** as a symmetric Hadamard matrix
([PROPUS_MIYAMOTO.md](../PROPUS_MIYAMOTO.md#status)): Seberry and Yamada's Williamson-type matrices with a repeated
block (1990), the propus array (Seberry and Balonin, arXiv:1512.01732) and the input were all published, and we found
them combined nowhere.

- The Đoković–Seberry census (2025-05-28) lists 961 among its 595 undecided v for symmetric matrices; its
  Seberry–Balonin series uses the propus variation only for v = q + 2.
- Seberry and Yamada list 961 among the orders of symmetric Williamson-type matrices that Corollary 29(ii) gives (J.
  Combin. Math. Combin. Comput. 7 (1990), Remark 30). A symmetric Hadamard matrix also needs the propus array (2015).
- Miyamoto (J. Combin. Theory Ser. A 57 (1991), Theorem 5) claimed Hadamard matrices of order 4q "of symmetric block
  form" for every q ≡ 1 (mod 4) with a conference matrix of order q + 1 and a symmetric Hadamard matrix of order
  q − 1, which would include this order; the census says his proof does not give it
  ([PROPUS_MIYAMOTO.md](../PROPUS_MIYAMOTO.md#status)).
- Ordinary Hadamard matrices of order 3844 are classical (Miyamoto's theorem, from a Hadamard matrix of order q − 1 =
  960); the symmetry is what is new.

The sources searched are listed in [NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q (and the field), g, the description of the input M, and the SHA-256 of the matrix |

The matrix itself is not stored; `verify/verify.py` rebuilds it.

## Verification

```
python3 verify/verify.py symmetric/3844
```

`verify/verify.py` (construction `symmetric-propus-miyamoto`) rebuilds M and checks that it is a symmetric Hadamard
matrix, builds U and V, checks that the W_k are ±1, symmetric, pairwise amicable, with Σ W_k W_kᵀ = 3844 I and
W₃ = W₄, assembles H and checks H = Hᵀ, H Hᵀ = 3844 I and the SHA-256 of the canonical CSV.

SHA-256 of `H3844.csv`: `b34c36b2b2b98c4b1ab9d8a68c4951a2a4b11651f8533b6ffe66d87bc1efffb5`
