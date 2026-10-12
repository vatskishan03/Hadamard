# Skew-Hadamard matrix of order 2404

2404 = 4 · 601 (prime).

## Construction

A skew-Hadamard matrix is a Hadamard matrix H with H + Hᵀ = 2I.

**Difference family.** Four ±1 sequences x₁, …, x₄ on ℤ₆₀₁ satisfy two conditions:

- each is constant on the orbits of the multiplier subgroup
  H = ⟨2⟩ ≤ ℤ₆₀₁*, of order 25: {1, 2, 4, 8, 16, 27, 32, 54, 64, 108, 128, 157, 216, 245, 256, 263, 301, 314, 379, 423, 432, 451, 490, 512, 526};
- their periodic autocorrelations sum to zero at every nonzero shift.

`certificate.json` lists, for each sequence, one representative r of each orbit rH on which the sequence is −1.

**Skew first block.** The first sequence is of skew type: x₁(0) = 1 and x₁(−t) = −x₁(t) for t ≠ 0. This is possible because −1 ∉ H, so negation pairs the orbits.

**Row sums.** The row sums are 1, 1, 49, 1, and 1² + 1² + 49² + 1² = 2404.

**The matrix.** Let A, B, C, D be the circulant matrices with these first rows (`sequences.txt`) and R the back-diagonal identity. The Goethals–Seidel array

```
    [  A     BR     CR     DR  ]
H = [ -BR    A     DᵀR   -CᵀR ]
    [ -CR   -DᵀR    A     BᵀR ]
    [ -DR    CᵀR  -BᵀR    A   ]
```

is then a Hadamard matrix of order 2404. Because A + Aᵀ = 2I, it satisfies H + Hᵀ = 2I (Wallis–Whiteman).

## Search

1. `attack/skew_triage.py` ranks skew Goethals–Seidel multiplier structures with a Gaussian model. For H = ⟨2⟩ it predicts about 2·10⁶ solutions:
   - 25 orbits and 12 autocorrelation classes;
   - 4096 skew first blocks;
   - pools of 2.5·10⁶ to 2.7·10⁶ sequences for the other three blocks, row sums (1, 49, 1).
2. `attack/skew_ktree.py` (four-list k-tree matcher `attack/match4.c`) enumerates the pools and matches the four lists.
3. It found 16 solutions among the 15 most likely partial-sum values, in 2 s.

## Status

To our knowledge this is the first skew-Hadamard matrix of order 2404.

- Cati and Pasechnik, *A database of constructions of Hadamard matrices*, arXiv:2411.18897v2 (Aug. 2025), Table 5, lists 601(5). So the smallest known skew order of the form 2ᵐ · 601 was 19232, and 2404 was unknown.
- A literature search in October 2026 found no later construction.

The doubling construction (Seberry) then also gives skew-Hadamard matrices of orders 4808, 9616.

## Verification

```
python3 verify/verify.py skew/2404
gzip -dc skew/2404/H2404.csv.gz | ./check_matrix
```

SHA-256 of `H2404.csv`:
`a67b40f854e9bfa1e0b80302dceb0054343baf975e15ee819637b96775d82130`
