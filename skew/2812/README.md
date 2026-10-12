# Skew-Hadamard matrix of order 2812

2812 = 4 · 703, with 703 = 19 · 37.

## Construction

A skew-Hadamard matrix is a Hadamard matrix H with H + Hᵀ = 2I.

**Difference family.** Four ±1 sequences x₁, …, x₄ on ℤ₇₀₃ satisfy two conditions:

- each is constant on the orbits of the multiplier subgroup
  H = ⟨16, 26⟩ ≤ ℤ₇₀₃*, of order 27: {1, 16, 26, 44, 118, 121, 157, 218, 248, 256, 271, 329, 334, 343, 367, 403, 416, 423, 441, 453, 482, 530, 567, 581, 636, 676, 682};
- their periodic autocorrelations sum to zero at every nonzero shift.

`certificate.json` lists, for each sequence, one representative r of each orbit rH on which the sequence is −1.

**Skew first block.** The first sequence is of skew type: x₁(0) = 1 and x₁(−t) = −x₁(t) for t ≠ 0. This is possible because −1 ∉ H, so negation pairs the orbits.

**Row sums.** The row sums are 1, 1, 1, 53, and 1² + 1² + 1² + 53² = 2812.

**The matrix.** Let A, B, C, D be the circulant matrices with these first rows (`sequences.txt`) and R the back-diagonal identity. The Goethals–Seidel array

```
    [  A     BR     CR     DR  ]
H = [ -BR    A     DᵀR   -CᵀR ]
    [ -CR   -DᵀR    A     BᵀR ]
    [ -DR    CᵀR  -BᵀR    A   ]
```

is then a Hadamard matrix of order 2812. Because A + Aᵀ = 2I, it satisfies H + Hᵀ = 2I (Wallis–Whiteman).

## Search

1. **Exact compression screen** (`attack/census_screen.py`). For every prime p | 703 the compressions of all four blocks to ℤ_p are listed exactly and joined in integers. H = ⟨16, 26⟩ (order 27, a non-cyclic subgroup of ℤ₇₀₃* ≅ ℤ₁₈ × ℤ₃₆) with rows (1, 1, 1, 53) survives at p = 19 (704 joint solutions) and at p = 37 (148,864).
2. **Compression chain** (`attack/chain.py`, run through `attack/screen_chain.py` in a breadth-first pass over the census survivors: at most two level-1 branches per structure). A dynamic program over the 31 H-orbits lists each block's joint compressions to ℤ₁₉ and ℤ₃₇. The ℤ₁₉ compressions are joined exactly (704 branches, heaviest first), and inside a branch the ℤ₃₇ compressions are joined. Each consistent combination is expanded by meet in the middle over the orbits, and the four lists are joined on all 15 autocorrelation classes (`attack/find4b.c`).
3. The heaviest branch held 49 combinations. Their final joins counted 3 solutions, two of which were rebuilt and verified as skew-Hadamard matrices. The search is exact and uses no randomness.

The Gaussian model of `attack/skew_rank.py` predicted 10^4.6 solutions for this structure. For composite v that model is unreliable, and an earlier sampled k-tree estimate for 703 had put the work at days.

## Reproduction

| | |
|---|---|
| Code | commit `a4cf889` |
| Command | `python3 attack/chain.py 703 16,26 1,1,53 --maxbranches 2` |
| Randomness and caps | none (exact search); default caps --maxfibre 3e6, --maxstates 2e7 |
| Result of that run | found in level-1 branch 1 of 704, after 49 combinations, in 47 minutes on 4 shared cores |

The tools write hit logs, temporary files and compiled helpers to `$HADAMARD_WORK` (default `<tmp>/hadamard-work`). The certificate does not depend on the search: `verify/verify.py` rebuilds and checks the matrix from `certificate.json` alone.

## Status

To our knowledge this is the first skew-Hadamard matrix of order 2812.

- Cati and Pasechnik, *A database of constructions of Hadamard matrices*, arXiv:2411.18897v2 (Aug. 2025), Table 5, lists 703(3). So the smallest known skew order of the form 2ᵐ · 703 was 5624, and 2812 was unknown.
- A literature search in October 2026 found no later construction.

## Verification

```
python3 verify/verify.py skew/2812
gzip -dc skew/2812/H2812.csv.gz | ./check_matrix
```

SHA-256 of `H2812.csv`:
`f2e13ba97b287a4cc1d1f51b3911d3aaca61f31a07b6b63b94cd941459269b24`
