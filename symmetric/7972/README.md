# Symmetric Hadamard matrix of order 7972

7972 = 4 · 1993, with 1993 a prime ≡ 1 (mod 8) and (q − 1)/2 = 996 = 4 · 249.

## Construction

A Hadamard matrix H is symmetric if H = Hᵀ. This one is the case q = 1993 of the construction in
[PROPUS_WEIGHING.md](../PROPUS_WEIGHING.md): the U-part of Seberry and Yamada's Corollary 29(ii) in Seberry and
Balonin's propus array, with its V-part replaced by a symmetric weighing matrix of order 3984 and weight 1993 built from a
Hadamard matrix K of order 996 that need not be symmetric.

- **The input.** K: our skew-Hadamard matrix of order 996 ([skew/996](../../skew/996)).
- **The field.** GF(1993) = ℤ_1993, g = 5; χ is the quadratic character, s_i = g²ⁱ and t_j = g¹⁻²ʲ (i, j < 996).
- **The U-part.** U₁ = −[χ(s_i − s_j)], U₂ = [χ(s_i − t_j)], and H₀ the propus array
  [[W₁, W₃, W₃, W₂], [W₃, W₂, −W₁, −W₃], [W₃, −W₁, −W₂, W₃], [W₂, −W₃, W₃, −W₁]] of W₁ = [[1, −e], [−eᵀ, U₁ ⊗ J₂]],
  W₂ = [[1, e], [eᵀ, U₂ ⊗ J₂]] and W₃ = [[1, e], [eᵀ, 0]], of order 7972; its cores vanish exactly where V
  (below) is nonzero.
- **The V-part.** K = [A | B] (A the first 498 columns), E = diag(1, −1) ⊗ I_498, X = [Aᵀ; Aᵀ], Y = [B, B] and
  V = [[E, X, −EX, 0], [Xᵀ, 0, I, Y], [−XᵀE, I, 0, −YE], [0, Yᵀ, −EYᵀ, E]], of order 3984: symmetric, V² = 1993I, and
  nonzero exactly where the cores of H₀ vanish.
- **The matrix.** H = H₀ + Λ(V), with V[(k, i), (l, j)] [[1, −1], [−1, 1]] placed at the core coordinates (k, i, ·),
  (l, j, ·). H₀Λ(V) = Λ(V)H₀ = 0 and H₀² + Λ(V)² = 7972I.

## Search

None: the matrix is determined by q, g and K. `python3 search/sym_weighing.py certify 1993` writes the certificate.

## Status

**New** as a symmetric Hadamard matrix.

- The Đoković–Seberry census (2025-05-28) lists 1993 among its 595 undecided v for symmetric matrices. The route of
  [PROPUS_MIYAMOTO.md](../PROPUS_MIYAMOTO.md) would need a symmetric Hadamard matrix of order 996 = 4 · 249, and 249 is
  undecided in the census too.
- The V-part from a Hadamard matrix that is not symmetric is the new step ([PROPUS_WEIGHING.md](../PROPUS_WEIGHING.md#status));
  we found it in none of the sources listed there. The input is ours (2026) and obtainable classically ([skew/996](../../skew/996)); any Hadamard matrix of order 996 would do.
- Ordinary Hadamard matrices of order 7972 are classical (Yamada's Theorem 3, 1989, from a Hadamard matrix of order 996;
  or Miyamoto's theorem); the symmetry is what is new.

The sources searched are listed in [NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, g, the complete certificate of K, and the SHA-256 of the matrix |

The matrix itself is not stored; `verify/verify.py` rebuilds it.

## Verification

```
python3 verify/verify.py symmetric/7972
```

`verify/verify.py` (construction `symmetric-propus-weighing`) rebuilds K from its certificate and checks it, builds
U₁, U₂ and V, checks that V is symmetric with V² = 1993I and the support of the zeros of H₀, assembles H and checks
H = Hᵀ, H Hᵀ = 7972 I and the SHA-256 of the canonical CSV.

SHA-256 of `H7972.csv`: `fe5f7d133f02a35fac25a32296496f59adfb0f2d1555bdeb6f664e2e18cae80d`
