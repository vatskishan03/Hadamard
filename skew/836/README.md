# Skew-Hadamard matrix of order 836

836 = 4 · 209, with 209 = 11 · 19.

## Construction

A skew-Hadamard matrix is a Hadamard matrix H with H + Hᵀ = 2I.

**Difference family.** Four ±1 sequences x₁, …, x₄ on ℤ₂₀₉ satisfy two conditions:

- each is constant on the orbits of the multiplier subgroup
  H = ⟨26⟩ = {1, 20, 26, 45, 49, 58, 64, 102, 115, 125, 144, 159, 163, 191, 201} ≤ ℤ₂₀₉*, of order 15;
- their periodic autocorrelations sum to zero at every nonzero shift.

`certificate.json` lists, for each sequence, one representative r of each orbit rH on which the sequence is −1.

**Skew first block.** The first sequence is of skew type: x₁(0) = 1 and x₁(−t) = −x₁(t) for t ≠ 0. This is possible because −1 ∉ H, so negation pairs the cosets.

**Row sums.** The row sums are 1, 9, 15, 23, and 1² + 9² + 15² + 23² = 836.

**The matrix.** Let A, B, C, D be the circulant matrices with these first rows (`sequences.txt`) and R the back-diagonal identity. The Goethals–Seidel array

```
    [  A     BR     CR     DR  ]
H = [ -BR    A     DᵀR   -CᵀR ]
    [ -CR   -DᵀR    A     BᵀR ]
    [ -DR    CᵀR  -BᵀR    A   ]
```

is then a Hadamard matrix of order 836. Because A + Aᵀ = 2I, it satisfies H + Hᵀ = 2I (Wallis–Whiteman).

## Search

1. `attack/skew_triage.py` ranks skew Goethals–Seidel multiplier structures with a Gaussian model. For H = ⟨26⟩ it predicts about 3 solutions:
   - 21 orbits and 10 autocorrelation classes;
   - 1024 skew first blocks;
   - pools of 2·10⁴ to 4·10⁴ sequences for the other three blocks.
2. `attack/skew_exhaustive.py` enumerates every admissible row-sum pattern and matches the four lists exactly.
3. The pattern (1; 9, 15, 23) gives 288 solutions, and the other two patterns give none.

The whole search takes under two minutes on four cores.

## Reproduction

| | |
|---|---|
| Code | commit `0cdf392` of this branch (the options used behave the same at later commits) |
| Command | `python3 attack/skew_exhaustive.py 209 26` |
| Randomness and caps | exact (no sampling, no seed); pair-work cap 4·10¹⁰ |
| Result of that run | pattern (1; 9, 15, 23): 288 solutions; the other two patterns: none |

The tools write their hit logs, temporary list files and compiled helpers to `$HADAMARD_WORK` (default
`<tmp>/hadamard-work`), never into the repository. `attack/skew_certify.py` (or `skew_group.py certify` for the
non-cyclic orders) turns a logged hit into the files in this directory. The certificate does not depend on
the search: `verify/verify.py` rebuilds and checks the matrix from `certificate.json` alone.

## Status

To our knowledge this is the first skew-Hadamard matrix of order 836.

- Cati and Pasechnik, *A database of constructions of Hadamard matrices*, arXiv:2411.18897v2 (Aug. 2025), Table 5, lists 209(4). So the smallest known skew order of the form 2ᵐ · 209 was 3344, and 836 was unknown.
- SageMath's Hadamard module lists 836 among its unknown skew-Hadamard orders.
- Đoković's 2007 list of odd n < 300 without a known skew-Hadamard matrix of order 4n includes 209.
- A literature search in October 2026 found no later construction.

**Novelty checklist** (audit of 2026-10-12, by independent agents that also re-verified the matrix):

1. Table 5 of Cati–Pasechnik lists 209(4): the smallest known skew order with odd part 209 was 2⁴ · 209 = 3344. This order is inside the range (≤ 1208) that they checked against implemented constructions.
2. Đoković's papers (1992–2024, alone and with Kotsireas, Golubitsky or Balonin), Koukouvinos–Stylianou
   (Discrete Math. 2008, Tables 1 and 3) and the Đoković–Seberry census of skew orders 4v for odd v < 2500
   (May 2025) do not construct it. The census lists 209 among its 737 undecided values.
3. The 2026 public deposits (Schumer; Alpöge et al.; JD-Jones-ASES; Schneider; Karoui; Piza) do not contain
   it. SageMath's `unknown_skew_hadamard` list (develop branch, fetched 2026-10-12) contains 836.
4. No infinite family reaches it: Paley I, Szekeres 2(q + 1), Wallis–Whiteman, Spence, Momihara–Xiang,
   the product theorem n(pᵏ + 1) and amicable products were checked.

This is a finding of nothing, not a proof that no earlier construction exists.

## Verification

```
python3 verify/verify.py skew/836
gzip -dc skew/836/H836.csv.gz | ./check_matrix
```

SHA-256 of `H836.csv`:
`952b2e48359cb1890c1dc7712ecb0b20187a2af495d097aa8fbb644a6418aeb8`
