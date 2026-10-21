# Symmetric Hadamard matrix of order 7652

7652 = 4 · 1913, with 1913 a prime ≡ 1 (mod 8) and (q − 1)/2 = 956 = 4 · 239.

## Construction

A Hadamard matrix H is symmetric if H = Hᵀ. This one is the case q = 1913 of the construction in
[PROPUS_WEIGHING.md](../PROPUS_WEIGHING.md): the U-part of Seberry and Yamada's Corollary 29(ii) in Seberry and
Balonin's propus array, with its V-part replaced by a symmetric weighing matrix of order 3824 and weight 1913 built
from a Hadamard matrix K of order 956 that need not be symmetric.

- **The input.** K: Đoković and Kotsireas's skew Goethals–Seidel matrix of order 956 (arXiv:1511.08734, §3;
  [search/data/skew956_dk.json](../../search/data/skew956_dk.json)).
- **The field.** GF(1913) = ℤ_1913, g = 3; χ is the quadratic character, s_i = g²ⁱ and t_j = g¹⁻²ʲ (i, j < 956).
- **The U-part.** U₁ = −[χ(s_i − s_j)], U₂ = [χ(s_i − t_j)], and H₀ the propus array
  [[W₁, W₃, W₃, W₂], [W₃, W₂, −W₁, −W₃], [W₃, −W₁, −W₂, W₃], [W₂, −W₃, W₃, −W₁]] of W₁ = [[1, −e], [−eᵀ, U₁ ⊗ J₂]],
  W₂ = [[1, e], [eᵀ, U₂ ⊗ J₂]] and W₃ = [[1, e], [eᵀ, 0]], of order 7652; its cores vanish exactly where V
  (below) is nonzero.
- **The V-part.** K = [A | B] (A the first 478 columns), E = diag(1, −1) ⊗ I_478, X = [Aᵀ; Aᵀ], Y = [B, B] and
  V = [[E, X, −EX, 0], [Xᵀ, 0, I, Y], [−XᵀE, I, 0, −YE], [0, Yᵀ, −EYᵀ, E]], of order 3824: symmetric, V² = 1913I, and
  nonzero exactly where the cores of H₀ vanish.
- **The matrix.** H = H₀ + Λ(V), with V[(k, i), (l, j)] [[1, −1], [−1, 1]] placed at the core coordinates (k, i, ·),
  (l, j, ·). H₀Λ(V) = Λ(V)H₀ = 0 and H₀² + Λ(V)² = 7652I.

## Search

None: the matrix is determined by q, g and K.

**Reproduction.**

| Command | Output |
|---|---|
| `python3 search/sym_weighing.py certify 1913` | `symmetric/7652/certificate.json`, identical to the stored file (the program refuses to overwrite it: run it in a copy of the repository without that file) |
| `python3 verify/verify.py symmetric/7652` | `ok`, SHA-256 `1cd421e82065d365f14085771c8f522aa6be382dd510ca7f38c2e7bd7528e78c` |

## Status

**New** as a symmetric Hadamard matrix.

- The Đoković–Seberry census (2025-05-28) lists 1913 among its 595 undecided v for symmetric matrices. The route of
  [PROPUS_MIYAMOTO.md](../PROPUS_MIYAMOTO.md) would need a symmetric Hadamard matrix of order 956 = 4 · 239, and 239
  is undecided in the census too.
- The V-part from a Hadamard matrix that is not symmetric is the new step
  ([PROPUS_WEIGHING.md](../PROPUS_WEIGHING.md#status)); we found it in none of the sources listed there. The input is
  classical (2015).
- Ordinary Hadamard matrices of order 7652 are classical (Yamada's Theorem 3, 1989, from a Hadamard matrix of order
  956; or Miyamoto's theorem); the symmetry is what is new.

The sources searched are listed in [NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, g, the complete certificate of K, and the SHA-256 of the matrix |

The matrix itself is not stored; `verify/verify.py` rebuilds it.

## Verification

```
python3 verify/verify.py symmetric/7652
```

`verify/verify.py` (construction `symmetric-propus-weighing`) rebuilds K from its certificate and checks it, builds
U₁, U₂ and V, checks that V is symmetric with V² = 1913I and the support of the zeros of H₀, assembles H and checks
H = Hᵀ, H Hᵀ = 7652 I and the SHA-256 of the canonical CSV.

SHA-256 of `H7652.csv`: `1cd421e82065d365f14085771c8f522aa6be382dd510ca7f38c2e7bd7528e78c`
