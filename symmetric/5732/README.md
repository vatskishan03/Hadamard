# Symmetric Hadamard matrix of order 5732

5732 = 4 · 1433, with 1433 a prime ≡ 1 (mod 8) and (q − 1)/2 = 716 = 4 · 179.

## Construction

A Hadamard matrix H is symmetric if H = Hᵀ. This one is the case q = 1433 of the construction in
[PROPUS_WEIGHING.md](../PROPUS_WEIGHING.md): the U-part of Seberry and Yamada's Corollary 29(ii) in Seberry and
Balonin's propus array, with its V-part replaced by a symmetric weighing matrix of order 2864 and weight 1433 built
from a Hadamard matrix K of order 716 that need not be symmetric.

- **The input.** K: the Hadamard matrix of order 716 of Alpöge, Voinov and Reynolds-Haertle (announced 12 August
  2026): a Goethals–Seidel array over ℤ₁₇₈ with a border of four rows and columns, its seeds and border as decoded in
  the JD-Jones-ASES Hadamard-B data ([search/data/h716_alpoge.json](../../search/data/h716_alpoge.json)), rebuilt by
  verify.py's `bordered-goethals-seidel`.
- **The field.** GF(1433) = ℤ_1433, g = 3; χ is the quadratic character, s_i = g²ⁱ and t_j = g¹⁻²ʲ (i, j < 716).
- **The U-part.** U₁ = −[χ(s_i − s_j)], U₂ = [χ(s_i − t_j)], and H₀ the propus array
  [[W₁, W₃, W₃, W₂], [W₃, W₂, −W₁, −W₃], [W₃, −W₁, −W₂, W₃], [W₂, −W₃, W₃, −W₁]] of W₁ = [[1, −e], [−eᵀ, U₁ ⊗ J₂]],
  W₂ = [[1, e], [eᵀ, U₂ ⊗ J₂]] and W₃ = [[1, e], [eᵀ, 0]], of order 5732; its cores vanish exactly where V
  (below) is nonzero.
- **The V-part.** K = [A | B] (A the first 358 columns), E = diag(1, −1) ⊗ I_358, X = [Aᵀ; Aᵀ], Y = [B, B] and
  V = [[E, X, −EX, 0], [Xᵀ, 0, I, Y], [−XᵀE, I, 0, −YE], [0, Yᵀ, −EYᵀ, E]], of order 2864: symmetric, V² = 1433I, and
  nonzero exactly where the cores of H₀ vanish.
- **The matrix.** H = H₀ + Λ(V), with V[(k, i), (l, j)] [[1, −1], [−1, 1]] placed at the core coordinates (k, i, ·),
  (l, j, ·). H₀Λ(V) = Λ(V)H₀ = 0 and H₀² + Λ(V)² = 5732I.

## Search

None: the matrix is determined by q, g and K.

**Reproduction.**

| Command | Output |
|---|---|
| `python3 search/sym_weighing.py certify 1433` | `symmetric/5732/certificate.json`, identical to the stored file (the program refuses to overwrite it: run it in a copy of the repository without that file) |
| `python3 verify/verify.py symmetric/5732` | `ok`, SHA-256 `3168182fb66f18ea1c50141ce233a59867b28c0bb22752d2e3942d31e859e9db` |

## Status

**New** as a symmetric Hadamard matrix.

- The Đoković–Seberry census (2025-05-28) lists 1433 among its 595 undecided v for symmetric matrices. The route of
  [PROPUS_MIYAMOTO.md](../PROPUS_MIYAMOTO.md) would need a symmetric Hadamard matrix of order 716 = 4 · 179, and 179
  is undecided in the census too.
- The V-part from a Hadamard matrix that is not symmetric is the new step
  ([PROPUS_WEIGHING.md](../PROPUS_WEIGHING.md#status)); we found it in none of the sources listed there. The input is
  from 2026: no Hadamard matrix of order 716 was known before Alpöge, Voinov and Reynolds-Haertle's.
- Hadamard matrices of order 5732 were already known (n = 1433 is not among the open n of Cati–Pasechnik's Table 4);
  the symmetry is what is new.

The sources searched are listed in [NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, g, the complete certificate of K, and the SHA-256 of the matrix |

The matrix itself is not stored; `verify/verify.py` rebuilds it.

## Verification

```
python3 verify/verify.py symmetric/5732
```

`verify/verify.py` (construction `symmetric-propus-weighing`) rebuilds K from its certificate and checks it, builds
U₁, U₂ and V, checks that V is symmetric with V² = 1433I and the support of the zeros of H₀, assembles H and checks
H = Hᵀ, H Hᵀ = 5732 I and the SHA-256 of the canonical CSV.

SHA-256 of `H5732.csv`: `3168182fb66f18ea1c50141ce233a59867b28c0bb22752d2e3942d31e859e9db`
