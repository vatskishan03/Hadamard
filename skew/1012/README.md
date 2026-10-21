# Skew-Hadamard matrix of order 1012

1012 = 4 · 253, with 253 = 11 · 23.

## Construction

A skew-Hadamard matrix is a Hadamard matrix H with H + Hᵀ = 2I.

**Difference family.** Four ±1 sequences x₁, …, x₄ on ℤ₂₅₃ satisfy two conditions:

- each is constant on the orbits of the multiplier subgroup H = ⟨12⟩ ≤ ℤ₂₅₃*, of order 11: {1, 12, 78, 100, 133, 144,
  177, 188, 210, 232, 243};
- their periodic autocorrelations sum to zero at every nonzero shift.

`certificate.json` lists, for each sequence, one representative r of each orbit rH on which the sequence is −1.

**Skew first block.** The first sequence is of skew type: x₁(0) = 1 and x₁(−t) = −x₁(t) for t ≠ 0. This is possible
because −1 ∉ H, so negation pairs the orbits.

**Row sums.** The row sums are 1, 5, 19, 25, and 1² + 5² + 19² + 25² = 1012.

**The matrix.** Let A, B, C, D be the circulant matrices with these first rows (`sequences.txt`) and R the back-diagonal
identity. The Goethals–Seidel array

```
    [  A     BR     CR     DR  ]
H = [ -BR    A     DᵀR   -CᵀR ]
    [ -CR   -DᵀR    A     BᵀR ]
    [ -DR    CᵀR  -BᵀR    A   ]
```

is then a Hadamard matrix of order 1012. Because A + Aᵀ = 2I, it satisfies H + Hᵀ = 2I (Wallis–Whiteman).

## Search

1. **Structure.** H = ⟨12⟩ has order 11 and 12 ≡ 1 (mod 11), so H fixes ℤ₁₁ × {0} pointwise and acts on ℤ₂₃ as the
   quadratic residues (ℤ₂₅₃ ≅ ℤ₁₁ × ℤ₂₃). There are 33 orbits and 16 autocorrelation classes.
2. **Low-energy k-tree** (`search/skew_lowE.py`, run on the census survivors of `search/census_screen.py`):
   - all 65,536 skew first blocks;
   - 406,765, 506,562 and 872,916 low-energy candidates with row sums 5, 19 and 25 (local minima from
     `search/skew_gen.c`, three descent steps; the generator saturated far below the 8·10⁶ and 3·10⁷ requested);
   - the four-list matcher `search/match4.c` found 96 solutions at the first value of τ, after 1.7 s of matching, and
     all 96 were verified as skew-Hadamard matrices.
3. **The Gaussian models were wrong by six orders of magnitude.** `search/skew_rank.py` gave the structure 10^2.6
   solutions, and the estimate inside the lists above was 10^−1.8. The 16 classes are strongly correlated here
   (covariance eigenvalues from 31 to 4,858; for the prime 613 they lie between 494 and 633), and a full-covariance
   model still gives only 10^−2.4 for lists of 65,536, 33,230, 42,858 and 42,429 candidates that hold 19,320 solutions
   (counted exactly with `search/join4x.c`). See [search/NOTES.md](../../search/NOTES.md) §15.
4. **Why: a Paley product.** Write the blocks over ℤ₁₁ × ℤ₂₃ as u ⊗ 1 + w ⊗ J + z ⊗ χ, where χ is the Legendre sequence
   of ℤ₂₃, J the all-ones sequence and u, w, z sequences on ℤ₁₁. Then the conditions reduce to two complementary
   quadruples (u_k) and (z_k) over ℤ₁₁ and one cross condition, with w nonzero at only two places.
   `search/paley_config.py` searches this reduced problem exactly. It lists 9,520 u-quadruples, and the first of them
   already extends in 100 ways. The counting identity behind the reduction, (q + 1)N_w − 4N₋ = 4m, is the classical one
   of Seberry Wallis (J. Combin. Theory A, 1976), Lemma 8, where J, J − 2I and a Paley-core block fill an orthogonal
   design of order 2ᵗ; here the outer structure is a Goethals–Seidel quadruple over ℤ₁₁ instead.

**Reproduction.**

| | |
|---|---|
| Code | commit `b3ee539` |
| Command | `python3 search/paley_config.py 253 --maxu 1` |
| Randomness and caps | none in the search (exact enumeration); which of the solutions join4x prints first depends on the OpenMP schedule. The certificate itself is solution 0 of the earlier run `python3 search/skew_lowE.py 253 12 5,19,25 --nbc 8e6 --nd 3e7 --steps 3 --maxtau 100000 --top 10 --seed 41` (commit `ac1249a`) |
| Result | 9,520 u-quadruples; the first gives 100 z-solutions, and the lifted matrices verify; 1 s on 4 cores. The k-tree run took 21 minutes of list generation and 1.7 s of matching. |

## Status

**New** as a skew-Hadamard matrix. Cati–Pasechnik Table 5 lists 253(4): the smallest known skew order 2ᵐ · 253 was 4048,
so 1012 and 2024 were unknown.

- 1012 is among the 41 skew orders up to 1208 that the introduction of Cati–Pasechnik counts as unknown.
- The Đoković–Seberry census lists 253 among its 737 undecided v. 253 is in none of its Paley, Szekeres,
  Szekeres–Whiteman or Spence series, nor among its 41 sporadic values. Đoković's earlier lists (1994, 2008,
  arXiv:0706.1973) also have 253 as open.
- SageMath lists 1012 in `unknown_skew_hadamard`, and no branch of `skew_hadamard_matrix` constructs it.
- Piza's entry for 1012 is an ordinary, non-skew matrix.
- No infinite family reaches it:
  - Paley I (1011 = 3 · 337);
  - Szekeres, Whiteman and Spence 2(q + 1), with q = 505 = 5 · 101;
  - Wallis–Whiteman and Momihara–Xiang 2ᵘ(qᵉ + 1), with qᵉ = 252;
  - Spence 4(q² + q + 1);
  - (p − 1)ᵘ + 1, (|t| + 1)(q + 1) and mn(n − 1);
  - products 2ᵗ ∏(qᵢ + 1), and skew × amicable products.

Not checked: Koukouvinos–Stylianou (Discrete Math. 2008), which later lists supersede.

The sources searched are listed in [NOVELTY.md](../../NOVELTY.md). Doubling (Seberry) gives a skew-Hadamard matrix of
order 2024.

## Files

| File | Contents |
|---|---|
| `certificate.json` | the four sequences as orbit representatives, H, the row sums |
| `sequences.txt` | the four sequences, one per line |
| `H1012.csv.gz` | the matrix |

## Verification

```
python3 verify/verify.py skew/1012
gzip -dc skew/1012/H1012.csv.gz | ./check_matrix
```

SHA-256 of `H1012.csv`: `0cd0696e6e4f269c019aaa197ba476a3de9ccec72131bb4ba677e37cddfa58e5`
