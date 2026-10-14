# Skew-Hadamard matrix of order 1860

1860 = 4 · 465, with 465 = 3 · 5 · 31.

## Construction

A skew-Hadamard matrix is a Hadamard matrix H with H + Hᵀ = 2I.

**Difference family.** Four ±1 sequences x₁, …, x₄ on ℤ₄₆₅ satisfy two conditions:

- each is constant on the orbits of the multiplier subgroup
  H = ⟨226⟩ ≤ ℤ₄₆₅*, of order 15: {1, 16, 76, 121, 196, 211, 226, 256, 286, 346, 361, 376, 391, 421, 436};
- their periodic autocorrelations sum to zero at every nonzero shift.

`certificate.json` lists, for each sequence, one representative r of each orbit rH on which the sequence is −1.

**Skew first block.** The first sequence is of skew type: x₁(0) = 1 and x₁(−t) = −x₁(t) for t ≠ 0. This is possible
because −1 ∉ H, so negation pairs the orbits.

**Row sums.** The row sums are 1, 25, 35, −3, and 1² + 25² + 35² + (−3)² = 1860.

**The matrix.** Let A, B, C, D be the circulant matrices with these first rows (`sequences.txt`) and R the
back-diagonal identity. The Goethals–Seidel array

```
    [  A     BR     CR     DR  ]
H = [ -BR    A     DᵀR   -CᵀR ]
    [ -CR   -DᵀR    A     BᵀR ]
    [ -DR    CᵀR  -BᵀR    A   ]
```

is then a Hadamard matrix of order 1860. Because A + Aᵀ = 2I, it satisfies H + Hᵀ = 2I (Wallis–Whiteman).

## Search

1. **Structure: a Paley product.** 465 = 15 · 31 and H = {h : h ≡ 1 (mod 15), h a square mod 31} = ⟨226⟩, of order 15.
   H fixes ℤ₁₅ × {0} pointwise and acts on ℤ₃₁ as the quadratic residues, so every block is u ⊗ 1 + w ⊗ J + z ⊗ χ over
   ℤ₁₅ × ℤ₃₁, with χ the Legendre sequence of ℤ₃₁, J the all-ones sequence and u, w, z sequences on ℤ₁₅ (see
   [search/NOTES.md](../../search/NOTES.md) §15).
2. **Reduction.** The identity coefficient of one of the three resulting equations forces (q + 1)N_w − 4N₋ = 4m, here
   32 N_w − 4 N₋ = 60, so exactly two positions carry w ≠ 0, one in each of two blocks, and N₋ = 1. Because 4m < 2q,
   the four sequences u_k and the four sequences z_k must then each have periodic autocorrelations summing to zero,
   and the cross-correlations Σ_k (z_k ⋆ u_k − u_k ⋆ z_k) must vanish.
3. **Exact two-stage search** (`search/paley_config.py`). Stage 1 lists all 209,536 u-quadruples (3 s). Stage 2 takes
   the first and joins all z-quadruples on the 14 remaining classes with `search/join4x.c`: 28 solutions. The first
   two were lifted to ℤ₄₆₅ and verified.

**Reproduction.**

| | |
|---|---|
| Code | commit `b3ee539` |
| Command | `python3 search/paley_config.py 465 --maxu 5` |
| Randomness and caps | none in the search (exact enumeration); which solutions join4x prints first depends on the OpenMP schedule |
| Result | 209,536 u-quadruples; the first extends to 28 z-solutions; 2 verified skew-Hadamard matrices after 16 s on 4 cores |

## Status

**New** as a skew-Hadamard matrix. A symmetric Hadamard matrix of order 1860 is known (Paley II with q = 929); it is
not skew.

- Cati–Pasechnik Table 5 lists 465(3), an entry inherited from the older tables; Koukouvinos–Stylianou (2008) also has
  465(3). The smallest known skew order of the form 2ᵐ · 465 was 2³ · 465 = 3720 (Paley I: 3719 is prime), so 1860 was
  unknown.
- The Đoković–Seberry census lists 465 among its 737 undecided values. 465 is in none of its skew series (Paley,
  Szekeres, Szekeres–Whiteman, Spence) nor among its sporadic values.
- SageMath's `skew_hadamard_matrix(1860)` has no route. Its Spence and complementary-difference-set branches need
  q = 929 ≡ 5 (mod 8), but 929 ≡ 1 (mod 8).
- Piza's deposit lists 1860 as Paley II, which is not skew.
- No infinite family reaches it:
  - Paley I: 1859 = 11 · 13² is not a prime power.
  - Szekeres, Whiteman, Spence and Momihara–Xiang 2(q + 1): q = 929 is prime but ≡ 1 (mod 8), and all these families
    need a prime power of a prime p ≡ 5 (mod 8).
  - Wallis–Whiteman 4(q + 1), with q = 464.
  - Spence 4(q² + q + 1), (p − 1)ᵘ + 1, (|t| + 1)(q + 1), mn(n − 1), products and skew × amicable products.
- The counting identity behind the search is classical in another setting. Seberry Wallis (J. Combin. Theory A, 1976),
  Lemma 8, plugs J, J − 2I and a Paley-core block into an orthogonal design of order 2ᵗ under
  a(v + 1) + b(v − 3) = 2ᵗ. With a = N_w − N₋ and b = N₋ its left side is our (q + 1)N_w − 4N₋, which here equals 4m
  rather than a power of 2. Seberry
  (Ars Combin., 1978), Theorems 6 and 7, give skew orders 2q(q − 1) and mn(n − 1) from orthogonal designs. Here the
  outer structure is a Goethals–Seidel quadruple over ℤ₁₅ instead; for 1860 those theorems would need designs that are
  not known (e.g. OD(620; 1, 155, 464)).

The sources searched are listed in [NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | the four sequences as orbit representatives, H, the row sums |
| `sequences.txt` | the four sequences, one per line |
| `H1860.csv.gz` | the matrix |

## Verification

```
python3 verify/verify.py skew/1860
gzip -dc skew/1860/H1860.csv.gz | ./check_matrix
```

SHA-256 of `H1860.csv`: `5730f7904a669ae11941d0d4b3737a038269132be9591b91fb53988866549d74`
