# Skew-Hadamard matrix of order 1236

1236 = 4 · 309, with 309 = 3 · 103.

## Construction

A skew-Hadamard matrix is a Hadamard matrix H with H + Hᵀ = 2I.

**Difference family.** Four ±1 sequences x₁, …, x₄ on ℤ₃₀₉ satisfy two conditions:

- each is constant on the orbits of the multiplier subgroup
  H = ⟨13⟩ ≤ ℤ₃₀₉*, of order 17: {1, 13, 34, 61, 64, 76, 79, 100, 112, 133, 169, 175, 184, 196, 214, 220, 229};
- their periodic autocorrelations sum to zero at every nonzero shift.

`certificate.json` lists, for each sequence, one representative r of each orbit rH on which the sequence is −1.

**Skew first block.** The first sequence is of skew type: x₁(0) = 1 and x₁(−t) = −x₁(t) for t ≠ 0. This is possible because −1 ∉ H, so negation pairs the orbits.

**Row sums.** The row sums are 1, 3, 35, 1, and 1² + 3² + 35² + 1² = 1236.

**The matrix.** Let A, B, C, D be the circulant matrices with these first rows (`sequences.txt`) and R the back-diagonal identity. The Goethals–Seidel array

```
    [  A     BR     CR     DR  ]
H = [ -BR    A     DᵀR   -CᵀR ]
    [ -CR   -DᵀR    A     BᵀR ]
    [ -DR    CᵀR  -BᵀR    A   ]
```

is then a Hadamard matrix of order 1236. Because A + Aᵀ = 2I, it satisfies H + Hᵀ = 2I (Wallis–Whiteman).

## Search

1. `attack/skew_triage.py` ranks skew Goethals–Seidel multiplier structures with a Gaussian model. For H = ⟨13⟩ it predicts about 16 solutions:
   - 21 orbits and 10 autocorrelation classes;
   - 1024 skew first blocks;
   - pools of 4.9·10⁴ to 1.5·10⁵ sequences for the other three blocks, row sums (3, 35, 1).
2. `attack/skew_ktree.py` (four-list k-tree matcher `attack/match4.c`) enumerates the pools and matches the four lists.
3. It found 1728 solutions for the most likely partial-sum value, in 0.4 s.

## Reproduction

| | |
|---|---|
| Code | commit `90ab79d` of this branch (the options used behave the same at later commits) |
| Command | `python3 attack/skew_ktree.py 309 13 --maxtau 400 --sols 2` |
| Randomness and caps | complete pools (no sampling, no seed); clamp S = 2; top 6; at most 400 values of τ |
| Result of that run | 1728 solutions at τ #0 |

The tools write their hit logs, temporary list files and compiled helpers to `$HADAMARD_WORK` (default
`<tmp>/hadamard-work`), never into the repository. `attack/skew_certify.py` (or `skew_group.py certify` for the
non-cyclic orders) turns a logged hit into the files in this directory. The certificate does not depend on
the search: `verify/verify.py` rebuilds and checks the matrix from `certificate.json` alone.

## Status

To our knowledge this is the first skew-Hadamard matrix of order 1236.

- Cati and Pasechnik, *A database of constructions of Hadamard matrices*, arXiv:2411.18897v2 (Aug. 2025), Table 5, lists 309(3). So the smallest known skew order of the form 2ᵐ · 309 was 2472, and 1236 was unknown.
- A literature search in October 2026 found no later construction.

## Verification

```
python3 verify/verify.py skew/1236
gzip -dc skew/1236/H1236.csv.gz | ./check_matrix
```

SHA-256 of `H1236.csv`:
`d507a23f4c879fce5dde40afdf6171ba5fcedde8ec0f1acbec7554778e726d64`
