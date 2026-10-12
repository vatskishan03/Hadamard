# Skew-Hadamard matrix of order 1852

1852 = 4 · 463 (prime).

## Construction

A skew-Hadamard matrix is a Hadamard matrix H with H + Hᵀ = 2I.

**Difference family.** Four ±1 sequences x₁, …, x₄ on ℤ₄₆₃ satisfy two conditions:

- each is constant on the orbits of the multiplier subgroup
  H = ⟨33⟩ ≤ ℤ₄₆₃*, of order 21: {1, 21, 33, 34, 118, 163, 169, 178, 182, 190, 196, 200, 230, 251, 286, 308, 318, 412, 441, 449, 450};
- their periodic autocorrelations sum to zero at every nonzero shift.

`certificate.json` lists, for each sequence, one representative r of each orbit rH on which the sequence is −1.

**Skew first block.** The first sequence is of skew type: x₁(0) = 1 and x₁(−t) = −x₁(t) for t ≠ 0. This is possible because −1 ∉ H, so negation pairs the orbits.

**Row sums.** The row sums are 1, 1, 43, 1, and 1² + 1² + 43² + 1² = 1852.

**The matrix.** Let A, B, C, D be the circulant matrices with these first rows (`sequences.txt`) and R the back-diagonal identity. The Goethals–Seidel array

```
    [  A     BR     CR     DR  ]
H = [ -BR    A     DᵀR   -CᵀR ]
    [ -CR   -DᵀR    A     BᵀR ]
    [ -DR    CᵀR  -BᵀR    A   ]
```

is then a Hadamard matrix of order 1852. Because A + Aᵀ = 2I, it satisfies H + Hᵀ = 2I (Wallis–Whiteman).

## Search

1. `attack/skew_triage.py` ranks skew Goethals–Seidel multiplier structures with a Gaussian model. For H = ⟨33⟩ it predicts about 3·10⁶ solutions:
   - 23 orbits and 11 autocorrelation classes;
   - 2048 skew first blocks;
   - pools of 6.5·10⁵ to 7.1·10⁵ sequences for the other three blocks, row sums (1, 43, 1).
2. `attack/skew_ktree.py` (four-list k-tree matcher `attack/match4.c`) enumerates the pools and matches the four lists.
3. It found 112 solutions for the most likely partial-sum value, in 0.4 s.

## Reproduction

| | |
|---|---|
| Code | commit `90ab79d` of this branch (the options used behave the same at later commits) |
| Command | `python3 attack/skew_ktree.py 463 33 --maxtau 400` |
| Randomness and caps | complete pools (no sampling, no seed); clamp S = 4; top 6; at most 400 values of τ |
| Result of that run | 112 solutions at τ #0 |

The tools write their hit logs, temporary list files and compiled helpers to `$HADAMARD_WORK` (default
`<tmp>/hadamard-work`), never into the repository. `attack/skew_certify.py` (or `skew_group.py certify` for the
non-cyclic orders) turns a logged hit into the files in this directory. The certificate does not depend on
the search: `verify/verify.py` rebuilds and checks the matrix from `certificate.json` alone.

## Status

To our knowledge this is the first skew-Hadamard matrix of order 1852.

- Cati and Pasechnik, *A database of constructions of Hadamard matrices*, arXiv:2411.18897v2 (Aug. 2025), Table 5, lists 463(7). So the smallest known skew order of the form 2ᵐ · 463 was 59264, and 1852 was unknown.
- A literature search in October 2026 found no later construction.

The doubling construction (Seberry) then also gives skew-Hadamard matrices of orders 3704, 7408, 14816, 29632.

## Verification

```
python3 verify/verify.py skew/1852
gzip -dc skew/1852/H1852.csv.gz | ./check_matrix
```

SHA-256 of `H1852.csv`:
`952cf790ef860f40f9d117a12dfaeb6cf1251980b9adc59e4bcb36452f7d066c`
