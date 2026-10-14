# Skew-Hadamard matrix of order 1380

1380 = 4 · 345, with 345 = 3 · 5 · 23.

## Construction

A skew-Hadamard matrix is a Hadamard matrix H with H + Hᵀ = 2I.

**Difference family.** Four ±1 sequences x₁, …, x₄ on ℤ₃₄₅ satisfy two conditions:

- each is constant on the orbits of the multiplier subgroup
  H = ⟨301⟩ ≤ ℤ₃₄₅*, of order 11: {1, 16, 31, 121, 151, 196, 211, 256, 271, 301, 331};
- their periodic autocorrelations sum to zero at every nonzero shift.

`certificate.json` lists, for each sequence, one representative r of each orbit rH on which the sequence is −1.

**Skew first block.** The first sequence is of skew type: x₁(0) = 1 and x₁(−t) = −x₁(t) for t ≠ 0. This is possible
because −1 ∉ H, so negation pairs the orbits.

**Row sums.** The row sums are 1, 19, 17, 27, and 1² + 19² + 17² + 27² = 1380.

**The matrix.** Let A, B, C, D be the circulant matrices with these first rows (`sequences.txt`) and R the
back-diagonal identity. The Goethals–Seidel array

```
    [  A     BR     CR     DR  ]
H = [ -BR    A     DᵀR   -CᵀR ]
    [ -CR   -DᵀR    A     BᵀR ]
    [ -DR    CᵀR  -BᵀR    A   ]
```

is then a Hadamard matrix of order 1380. Because A + Aᵀ = 2I, it satisfies H + Hᵀ = 2I (Wallis–Whiteman).

## Search

1. **Structure: a Paley product.** 345 = 15 · 23 and H = {h : h ≡ 1 (mod 15), h a square mod 23} = ⟨301⟩, of order 11.
   Every block is u ⊗ 1 + w ⊗ J + z ⊗ χ over ℤ₁₅ × ℤ₂₃ (χ the Legendre sequence of ℤ₂₃, J the all-ones sequence, u, w,
   z on ℤ₁₅); see [search/NOTES.md](../../search/NOTES.md) §15.
2. **Reduction.** The counting identity (q + 1)N_w − 4N₋ = 4m reads 24 N_w − 4 N₋ = 60, so N_w = 3 and N₋ = 3: three
   positions carry w ≠ 0. Two of them in one block would leave a term of size q² in the autocorrelations of the
   compressions, so they sit one in each of the blocks B, C, D, translated to 0, each with u = −2w there.
3. **Exact two-stage search** (`search/paley_config.py`, general mode). Stage 1 lists all 2,040,192 quadruples of
   compressions y_k = u_k + q w_k to ℤ₁₅ whose autocorrelations sum to zero (5 s). Stage 2 takes the first and joins
   all z-quadruples with the prescribed autocorrelation sums and zero cross terms (`search/join4x.c`): 22 solutions.
   The first two were lifted to ℤ₃₄₅ and verified.

**Reproduction.**

| | |
|---|---|
| Code | commit `5a713a9` |
| Command | `python3 search/paley_config.py 345 --maxu 6` |
| Randomness and caps | none in the search (exact enumeration); which solutions join4x prints first depends on the OpenMP schedule. |
| Result | 2,040,192 y-quadruples; the first extends to 22 z-solutions; 2 verified skew-Hadamard matrices after 21 s on 4 cores |

## Status

**New** as a skew-Hadamard matrix. Cati–Pasechnik Table 5 lists 345(4), as does Koukouvinos–Stylianou (2008): the
smallest known skew order of the form 2ᵐ · 345 was 2⁴ · 345 = 5520 (Paley I: 5519 is prime), so 1380 and 2760 were
unknown.

- The Đoković–Seberry census lists 345 among its 737 undecided values. 345 is in none of its skew series (Paley,
  Szekeres, Szekeres–Whiteman, Spence) nor among its 41 sporadic values.
- SageMath has no route to a skew matrix of order 1380: the constructions below fail, there are no difference-family
  data for v = 345, and 1380 ≡ 4 (mod 8), so doubling does not apply.
- Piza's deposit lists 1380 as "CW(15, 23)", an ordinary Cooper–Wallis matrix.
- No infinite family reaches it:
  - Paley I: 1379 = 7 · 197 is not a prime power.
  - 2(q + 1) with q = 689 = 13 · 53, and 2ᵘ(qᵉ + 1) with qᵉ = 344.
  - Spence 4(q² + q + 1), (p − 1)ᵘ + 1 and (|t| + 1)(q + 1).
  - Products and skew × amicable products: every split leaves an odd cofactor.
  - Seberry's mn(n − 1) theorem (Ars Combin. 6, 1978, Thm 7) would need an OD(460; 1, 115, 344), which is not recorded
    anywhere we looked.

The sources searched are listed in [NOVELTY.md](../../NOVELTY.md). Doubling (Seberry) gives a skew-Hadamard matrix of
order 2760.

## Files

| File | Contents |
|---|---|
| `certificate.json` | the four sequences as orbit representatives, H, the row sums |
| `sequences.txt` | the four sequences, one per line |
| `H1380.csv.gz` | the matrix |

## Verification

```
python3 verify/verify.py skew/1380
gzip -dc skew/1380/H1380.csv.gz | ./check_matrix
```

SHA-256 of `H1380.csv`: `14993cebe74e987180ca1910e6b3ea61cfc9b33ba06bf287d92c21ec3e09dc04`
