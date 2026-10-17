# Skew-Hadamard matrix of order 1564

1564 = 4 · 391, with 391 = 17 · 23.

## Construction

A skew-Hadamard matrix is a Hadamard matrix H with H + Hᵀ = 2I.

**Difference family.** Four ±1 sequences x₁, …, x₄ on ℤ₃₉₁ satisfy two conditions:

- each is constant on the orbits of the multiplier subgroup H = ⟨324⟩ ≤ ℤ₃₉₁*, of order 11: {1, 18, 35, 52, 154, 188,
  239, 256, 307, 324, 358};
- their periodic autocorrelations sum to zero at every nonzero shift.

`certificate.json` lists, for each sequence, one representative r of each orbit rH on which the sequence is −1.

**Skew first block.** The first sequence is of skew type: x₁(0) = 1 and x₁(−t) = −x₁(t) for t ≠ 0. This is possible
because −1 ∉ H, so negation pairs the orbits.

**Row sums.** The row sums are 1, 19, 19, 29, and 1² + 19² + 19² + 29² = 1564.

**The matrix.** Let A, B, C, D be the circulant matrices with these first rows (`sequences.txt`) and R the back-diagonal
identity. The Goethals–Seidel array

```
    [  A     BR     CR     DR  ]
H = [ -BR    A     DᵀR   -CᵀR ]
    [ -CR   -DᵀR    A     BᵀR ]
    [ -DR    CᵀR  -BᵀR    A   ]
```

is then a Hadamard matrix of order 1564. Because A + Aᵀ = 2I, it satisfies H + Hᵀ = 2I (Wallis–Whiteman).

## Search

1. **Structure: a Paley product.** 391 = 17 · 23 and H = {h : h ≡ 1 (mod 17), h a square mod 23} = ⟨324⟩, of order 11.
   Every block is u ⊗ 1 + w ⊗ J + z ⊗ χ over ℤ₁₇ × ℤ₂₃ (χ the Legendre sequence of ℤ₂₃); see
   [search/NOTES.md](../../search/NOTES.md) §15.
2. **Reduction.** The counting identity reads 24 N_w − 4 N₋ = 68, so N_w = 3 and N₋ = 1: one w-position in each of B, C
   and D, translated to 0; one of them has u = −2w there.
3. **Exact two-stage search** (`search/paley_config.py`, general mode). Stage 1: all 1,620,096 quadruples of
   compressions to ℤ₁₇ whose autocorrelations sum to zero (76 s). Stage 2 on the first: 18 z-solutions; the first two
   were lifted to ℤ₃₉₁ and verified.

**Reproduction.**

| | |
|---|---|
| Code | commit `5a713a9` |
| Command | `python3 search/paley_config.py 391 --maxu 6` |
| Randomness and caps | none in the search (exact enumeration); which solutions join4x prints first depends on the OpenMP schedule |
| Result | 1,620,096 y-quadruples; the first extends to 18 z-solutions; 2 verified skew-Hadamard matrices after 252 s on 4 cores |

## Status

**New** as a skew-Hadamard matrix. Ordinary Hadamard matrices of order 1564 have long been known; the skew property is
what is new.

- Cati–Pasechnik Table 5 lists 391(4), as do Seberry–Yamada (1992, Table 7.1) and Koukouvinos–Stylianou (2008): the
  smallest known skew order of the form 2ᵐ · 391 was 2⁴ · 391 = 6256, so 1564 and 3128 were unknown.
- The Đoković–Seberry census lists 391 among its 737 undecided values. 391 is in none of its skew series nor among its
  41 sporadic values.
- SageMath has no route to a skew matrix of order 1564: there are no difference-family data for v = 391, Paley I, Spence
  and the complementary-difference-set branches fail, and 1564 ≡ 4 (mod 8) rules out doubling.
- Piza's deposit lists 1564 as "CW(17, 23)", an ordinary Cooper–Wallis matrix.
- No infinite family reaches it:
  - Paley I: 1563 = 3 · 521 is not a prime power.
  - 2(q + 1) with q = 781 = 11 · 71, and 2ᵘ(qᵉ + 1) with qᵉ = 390.
  - Spence 4(q² + q + 1), (p − 1)ᵘ + 1 and (|t| + 1)(q + 1).
  - Seberry's theorems 6 and 7 (Ars Combin. 1978), products and skew × amicable products.
  - Good matrices of order 391 are not known.

The sources searched are listed in [NOVELTY.md](../../NOVELTY.md). Doubling (Seberry) gives a skew-Hadamard matrix of
order 3128.

## Files

| File | Contents |
|---|---|
| `certificate.json` | the four sequences as orbit representatives, H, the row sums |
| `sequences.txt` | the four sequences, one per line |
| `H1564.csv.gz` | the matrix |

## Verification

```
python3 verify/verify.py skew/1564
gzip -dc skew/1564/H1564.csv.gz | ./check_matrix
```

SHA-256 of `H1564.csv`: `7b6fb1f3f489edb7c8d1034b35e2eb4efc6122633adfdb1e45486de3c5af9c76`