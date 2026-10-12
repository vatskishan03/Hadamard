# Skew-Hadamard matrix of order 1412

1412 = 4 · 353 (prime).

## Construction

A skew-Hadamard matrix is a Hadamard matrix H with H + Hᵀ = 2I.

**Difference family.** Four ±1 sequences x₁, …, x₄ on ℤ₃₅₃ satisfy two conditions:

- each is constant on the orbits of the multiplier subgroup
  H = ⟨22⟩ ≤ ℤ₃₅₃*, of order 11: {1, 22, 58, 131, 140, 185, 187, 217, 231, 256, 337};
- their periodic autocorrelations sum to zero at every nonzero shift.

`certificate.json` lists, for each sequence, one representative r of each orbit rH on which the sequence is −1.

**Skew first block.** The first sequence is of skew type: x₁(0) = 1 and x₁(−t) = −x₁(t) for t ≠ 0. This is possible because −1 ∉ H, so negation pairs the orbits.

**Row sums.** The row sums are 1, 21, 23, 21, and 1² + 21² + 23² + 21² = 1412.

**The matrix.** Let A, B, C, D be the circulant matrices with these first rows (`sequences.txt`) and R the back-diagonal identity. The Goethals–Seidel array

```
    [  A     BR     CR     DR  ]
H = [ -BR    A     DᵀR   -CᵀR ]
    [ -CR   -DᵀR    A     BᵀR ]
    [ -DR    CᵀR  -BᵀR    A   ]
```

is then a Hadamard matrix of order 1412. Because A + Aᵀ = 2I, it satisfies H + Hᵀ = 2I (Wallis–Whiteman).

## Search

1. `attack/skew_rank.py` ranks skew Goethals–Seidel multiplier structures with a Gaussian model, counting only the solutions that survive when each list is cut to a sample. For H = ⟨22⟩ (order 11) the model predicts 10⁹ solutions in all:
   - 33 orbits and 16 autocorrelation classes;
   - 65536 skew first blocks;
   - pools of 5.7·10⁸ sequences for each of the other three blocks, row sums (21, 21, 23).
2. `attack/skew_ktree.py` (four-list k-tree matcher `attack/match4.c`) matched uniform samples of 1.2·10⁷ sequences from each pool.
3. It found 2 solutions among the 33 most likely partial-sum values.

## Status

To our knowledge this is the first skew-Hadamard matrix of order 1412.

- Cati and Pasechnik, *A database of constructions of Hadamard matrices*, arXiv:2411.18897v2 (Aug. 2025), Table 5, lists 353(4). So the smallest known skew order of the form 2ᵐ · 353 was 5648, and 1412 was unknown.
- A literature search in October 2026 found no later construction.

The doubling construction (Seberry) then also gives skew-Hadamard matrices of orders 2824.

## Verification

```
python3 verify/verify.py skew/1412
gzip -dc skew/1412/H1412.csv.gz | ./check_matrix
```

SHA-256 of `H1412.csv`:
`a3298cc6fed6b9c74ce477f5332a5a051124078423b24b052a589e771d7d253e`
