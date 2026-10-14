# Skew-Hadamard matrix of order 3052

3052 = 4 · 763, with 763 = 7 · 109.

## Construction

A skew-Hadamard matrix is a Hadamard matrix H with H + Hᵀ = 2I.

**Difference family.** Four ±1 sequences x₁, …, x₄ on ℤ₇₆₃ satisfy two conditions:

- each is constant on the orbits of the multiplier subgroup
  H = ⟨9⟩ ≤ ℤ₇₆₃*, of order 27: {1, 9, 81, 116, 221, 233, 240, 281, 296, 298, 323, 352, 365, 372, 375, 393, 424, 457, 463, 471, 485, 550, 561, 571, 618, 634, 729};
- their periodic autocorrelations sum to zero at every nonzero shift.

`certificate.json` lists, for each sequence, one representative r of each orbit rH on which the sequence is −1.

**Skew first block.** The first sequence is of skew type: x₁(0) = 1 and x₁(−t) = −x₁(t) for t ≠ 0. This is possible because −1 ∉ H, so negation pairs the orbits.

**Row sums.** The row sums are 1, 1, 5, 55, and 1² + 1² + 5² + 55² = 3052.

**The matrix.** Let A, B, C, D be the circulant matrices with these first rows (`sequences.txt`) and R the back-diagonal identity. The Goethals–Seidel array

```
    [  A     BR     CR     DR  ]
H = [ -BR    A     DᵀR   -CᵀR ]
    [ -CR   -DᵀR    A     BᵀR ]
    [ -DR    CᵀR  -BᵀR    A   ]
```

is then a Hadamard matrix of order 3052. Because A + Aᵀ = 2I, it satisfies H + Hᵀ = 2I (Wallis–Whiteman).

## Search

1. **Exact compression screen** (`search/census_screen.py`). For every prime p | 763 the compressions of all four blocks to ℤ_p are listed exactly and joined in integers. H = ⟨9⟩ with rows (1, 1, 5, 55) survives at p = 7 (40 joint solutions) and at p = 109 (61,688).
2. **Compression chain** (`search/chain.py`, run through `search/screen_chain.py`). A dynamic program over the 31 H-orbits lists each block's joint compressions to ℤ₇ and ℤ₁₀₉. The ℤ₇ compressions are joined exactly (40 branches), and inside each branch the ℤ₁₀₉ compressions are joined. Each consistent combination is expanded by meet in the middle over the orbits, and the four lists are joined on all 15 autocorrelation classes (`search/find4b.c`).
3. The first branch held 31 combinations. Their final joins counted 3 solutions, two of which were rebuilt and verified as skew-Hadamard matrices. The search is exact and uses no randomness.

The Gaussian model of `search/skew_rank.py` predicted 10^5.5 solutions for this structure. For composite v that model is unreliable, and an earlier sampled k-tree run (400 values of τ) had found nothing here.

## Reproduction

| | |
|---|---|
| Code | commit `cd8c03a` |
| Command | `python3 search/chain.py 763 9 1,5,55` |
| Randomness and caps | none (exact search); default caps --maxfibre 3e6, --maxstates 2e7 |
| Result of that run | found in level-1 branch 1 of 40, after 31 combinations, in 86 minutes on 4 shared cores |

The tools write hit logs, temporary files and compiled helpers to `$HADAMARD_WORK` (default `<tmp>/hadamard-work`). The certificate does not depend on the search: `verify/verify.py` rebuilds and checks the matrix from `certificate.json` alone.

## Status

To our knowledge this is the first skew-Hadamard matrix of order 3052.

- Cati and Pasechnik, *A database of constructions of Hadamard matrices*, arXiv:2411.18897v2 (Aug. 2025), Table 5, lists 763(4). So the smallest known skew order of the form 2ᵐ · 763 was 12208, and 3052 was unknown.
- A literature search in October 2026 found no later construction.

The doubling construction (Seberry) then also gives a skew-Hadamard matrix of order 6104.

**Novelty checklist** (audit of 2026-10-12 by a research agent):

1. Table 5 of Cati–Pasechnik lists 763(4). So the smallest known skew order with odd part 763 was
   2⁴ · 763 = 12208. This order is above 1208, the range they checked against implemented constructions,
   so the entry was not taken on trust (items 2 and 3).
2. Đoković's papers (1992–2024, alone and with Kotsireas, Golubitsky or Balonin), Koukouvinos–Stylianou
   (Discrete Math. 2008) and the Đoković–Seberry census of skew orders 4v for odd v < 2500 (May 2025) do not
   construct it. The census lists 763 among its 737 undecided values and not among its sporadic ones.
3. The 2025–2026 papers and deposits we checked do not contain it:
   - arXiv;
   - Zenodo: Schumer, Piza, Karoui and Kline;
   - the JD-Jones-ASES repositories Hadamard-T/M/B;
   - Alpöge et al. and Schneider.

   SageMath's `unknown_skew_hadamard` list stops at 1196, below this order.
4. No infinite family reaches it:
   - Paley I (3051 = 3³ · 113);
   - Szekeres and Whiteman 2(q + 1) (q = 1525);
   - Wallis–Whiteman 4(q + 1);
   - (|t| + 1)(q + 1);
   - Spence 4(q² + q + 1);
   - Momihara–Xiang;
   - (p − 1)ᵘ + 1;
   - the product theorem 2ᵗ ∏(pᵏ + 1);
   - skew × amicable products.

This is a finding of nothing, not a proof that no earlier construction exists.

**Consequence.** By Yamada's theorem (q = 6101), this matrix gives a Hadamard matrix of order 24412, which
is probably new; see [24412](../../hadamard/24412).

## Verification

```
python3 verify/verify.py skew/3052
gzip -dc skew/3052/H3052.csv.gz | ./check_matrix
```

SHA-256 of `H3052.csv`:
`cdf0881b097a7d8b91a61329f74f36b7f3fd3ded07016bafb218ec0e664cf93d`
