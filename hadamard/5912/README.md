# Hadamard matrix of order 5912

5912 = 8 · 739, with 739 a prime ≡ 3 (mod 4).

## Construction

**Spence's bordered construction over a group of odd order** ([search/eightp](../../search/eightp/README.md)). Let
G = ℤ₂₁₁ × F₇, of order v = 1477 = 2 · 739 − 1. Four ±1 sequences a, b, c, d on G satisfy

```
Σ PAF(g) = −4   for every g ≠ 0.
```

Then Σ r² = 4v − 4(v − 1) = 4 for their row sums r, and v is odd, so every row sum is ±1; here they are (1, 1, 1, 1).
Let S = [[a, b, c, d], [−b, a, d, −c], [−c, −d, a, b], [−d, c, −b, a]] be the Goethals–Seidel array of the row sums
(S Sᵀ = 4I), M the Goethals–Seidel array of the four group-developed blocks (order 4v), and

```
H = [ −S       S ⊗ 1ᵀ ]
    [ S ⊗ 1    M      ].
```

M Mᵀ = (4v + 4) I − 4 I₄ ⊗ J and the row sums give the borders, so H Hᵀ = 4(v + 1) I = 5912 I. Spence stated the
construction for groups of even order (row sums (2, 0, 0, 0)); the odd case needs only S in place of his border.
`verify/verify.py` (construction `spence`) accepts both.

**Symmetry: a Paley plug-in.** The sequences are constant on the orbits of Γ = ⟨19⟩ × ⟨2⟩ of order 45 (multiplication
by 19, of order 15, on ℤ₂₁₁, and by the squares of F₇*), which has 45 orbits on G. On each fibre {k} × F₇ a sequence is
either a constant on F₇* or ± the quadratic character χ₇ on F₇*, with a free sign at x = 0:
X(k, x) = e(k) + t(k) χ₇(x) + s(k) [x = 0]. The four blocks have 0, 30, 30 and 60 constant fibres.

## Search

The engine of [search/eightp/gs](../../search/eightp/README.md) ([search/NOTES.md](../../search/NOTES.md) §92, §108)
solves the K-level system of the plug-in in two levels. Level 1 joins the four blocks on their compressed
autocorrelations; level 2 chooses the signs t off the defects by a meet in the middle with row-sum and positivity
filters. The shape ranked first in the re-screen of §107. Before the production run, a near-miss calibration (level 2
with 4 or 5 of its 14 constraints left out) predicted 2.3 · 10⁻⁵ completions per distinct level-1 instance.

Two sampling processes drew level-1 instances of four allocations of constant fibres in turn, with tokens
(30, 30, 30, 30), (0, 30, 30, 60), (0, 0, 60, 60) and (0, 0, 30, 90). Duplicates under the reversal of a block, an
exact symmetry of the condition, were removed. The hit came after 75,726 distinct instances, about 4.3 core-hours,
against 1.7 expected.

**Reproduction.**

| | |
|---|---|
| Code | `search/eightp/gs/prod_5912.py` and `repro_5912.py` ([search/eightp](../../search/eightp/README.md#order-5912)) |
| Command | `cd search/eightp/gs; python3 repro_5912.py`, after the build steps of that README; about 25 min including the matrix; with `--fast --no-build`, 9 s for the sequences only |
| Randomness and caps | the sampling rounds are functions of (base seed, round index); the hit is round index 33 of base seed 1002, token (0, 30, 30, 60) |
| Result | `certificate.json`, identical; matrix SHA-256 `b4846b8e…` |

The matrix was checked by `verify/verify.py`, by the search program's circulant build after the Chinese remainder map
ℤ₂₁₁ × ℤ₇ ≅ ℤ₁₄₇₇, and by two independent verifiers with their own code (matrix products, entry formulas without
products, int32 and bit-packed Gram tests). The first verifier ran negative controls, and the second reproduced the
SHA-256 from its own build.

## Status

**New.** Cati–Pasechnik's Table 4 records 739(6). Craigen–Holzmann–Kharaghani's rule, which the table omits, gives
16 · 739 = 8(26 + 1452), so the smallest known order 2ᵐ · 739 was 16 · 739. The order 8 · 739 is one of the six orders
8p up to 10 000 that the survey of Egan, Kotsireas, London and Ó Catháin (arXiv:2610.27175, October 2026, Table 5)
lists as open. We know of no earlier construction of a Hadamard matrix of this order:

- Paley I and II: 5911 = 23 · 257 and 2955 = 3 · 5 · 197 are not prime powers.
- Doubling needs a Hadamard matrix of order 2956 = 4 · 739, which is not known ([FRONTIER.md](../../FRONTIER.md)).
- Miyamoto's theorem for order 8q needs 2q − 3 = 1475 = 5² · 59 to be a prime power. His Theorem 7 (see
  [COMPLEX_MIYAMOTO.md](../COMPLEX_MIYAMOTO.md)) needs a quaternary Hadamard matrix of order 738, which is open.
- Seberry–Yamada's order 8q needs Williamson-type matrices of order 369. Their Table A.1 gives row 369 only for
  Williamson-type matrices of even order 738.
- Craigen–Holzmann–Kharaghani's 8(x + y) needs complex Golay lengths x, y with x + y = 739. Kamali–Kharaghani's dihedral
  Golay sequences give 8(g₁g + g₂) for Golay lengths g₁, g₂ and g = 7, 9, 15. Neither reaches 739.

The Hadamard Matrix Explorer (math74.com, version 57) also lists 5912 as open, and SageMath returns Unknown. No public
list or repository of 2026 that we read gives it. The sources searched are listed in [NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | the group [211, 7] and the four sequences (index 7k + x for (k, x) ∈ ℤ₂₁₁ × ℤ₇), the SHA-256 of the matrix |

The matrix (5912² entries) is not stored; `verify/verify.py` rebuilds it (about ten minutes and 1.5 GB).

## Verification

```
python3 verify/verify.py 5912
```

SHA-256 of `H5912.csv`: `b4846b8e76f059552a232a4a6d0f815022c3b178d9103d044965894b6f9c4aa6`
