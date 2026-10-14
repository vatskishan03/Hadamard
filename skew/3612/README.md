# Skew-Hadamard matrix of order 3612

3612 = 4 · 903, with 903 = 3 · 7 · 43.

## Construction

A skew-Hadamard matrix is a Hadamard matrix H with H + Hᵀ = 2I.

**Difference family.** Four ±1 sequences x₁, …, x₄ on ℤ₉₀₃ satisfy two conditions:

- each is constant on the orbits of the multiplier subgroup
  H = ⟨568⟩ ≤ ℤ₉₀₃*, of order 21: {1, 64, 127, 169, 232, 253, 274, 316, 337, 358, 379, 400, 484, 526, 547, 568, 694, 778, 799, 841, 883};
- their periodic autocorrelations sum to zero at every nonzero shift.

`certificate.json` lists, for each sequence, one representative r of each orbit rH on which the sequence is −1.

**Skew first block.** The first sequence is of skew type: x₁(0) = 1 and x₁(−t) = −x₁(t) for t ≠ 0. This is possible because −1 ∉ H, so negation pairs the orbits.

**Row sums.** The row sums are 1, 43, 41, −9, and 1² + 43² + 41² + (−9)² = 3612.

**The matrix.** Let A, B, C, D be the circulant matrices with these first rows (`sequences.txt`) and R the back-diagonal identity. The Goethals–Seidel array

```
    [  A     BR     CR     DR  ]
H = [ -BR    A     DᵀR   -CᵀR ]
    [ -CR   -DᵀR    A     BᵀR ]
    [ -DR    CᵀR  -BᵀR    A   ]
```

is then a Hadamard matrix of order 3612. Because A + Aᵀ = 2I, it satisfies H + Hᵀ = 2I (Wallis–Whiteman).

## Search

1. **Structure: a Paley product.** 903 = 21 · 43 and H = {h : h ≡ 1 (mod 21), h a square mod 43} = ⟨568⟩, of order 21. H fixes ℤ₂₁ × {0} pointwise and acts on ℤ₄₃ as the quadratic residues, so every block is u ⊗ 1 + w ⊗ J + z ⊗ χ over ℤ₂₁ × ℤ₄₃, with χ the Legendre sequence of ℤ₄₃, J the all-ones sequence and u, w, z sequences on ℤ₂₁. At each position either w = 0 and u, z = ±1, or z = 0, w = ±1 and u ∈ {0, −2w}.
2. **Reduction.** The identity coefficient forces (q + 1)N_w − 4N₋ = 4m, here 44 N_w − 4 N₋ = 84, so N_w = 2 and N₋ = 1: two positions carry w ≠ 0, one in each of blocks 2 and 3 (placed at 0), and the one in block 2 has u = −2w. Because 4m < 2q, the u_k and the z_k must each have periodic autocorrelations summing to zero off the spikes' contributions, which the exact joins enforce.
3. **Two-stage search** (`search/pcg.py`, which extends `search/paley_config.py` to any configuration of the w ≠ 0 positions). With y_k = u_k + q w_k, the compression of block k to ℤ₂₁, the conditions are Σ_k y_k y_k* = 4v, Σ_k z_k z_k* = Σ_k (u_k w_k* + w_k u_k* + q w_k w_k*) and Σ_k (z_k u_k* − u_k z_k*) = 0. Stage Y joins the y-lists exactly; stage Z joins, for each y-quadruple, the z-lists on the PAF and cross classes with `search/join4x.c`.

## Reproduction

| | |
|---|---|
| Code | commit `ec11d71` |
| Command | `python3 search/pcg.py 903 --lmax 131072 --lzmax 2097152 --maxu 20 --stop --seed 2` |
| Randomness and caps | seed 2 for the sampled y-lists; the z-join is exhaustive; which solutions join4x prints first depends on the OpenMP schedule |
| Result of that run | y-lists of 1,024 (complete), 129,029, 129,037 and 131,072 (uniform samples of 2²⁰, 2²⁰ and 2²¹), 118,922 y-quadruples; for the first, the complete z-lists (2,048, 1,048,576, 1,048,576 and 2,097,152) give exactly 2 z-solutions, both lifted to ℤ₉₀₃ and verified; 10 minutes on 40 shared cores |

The tools write hit logs, temporary files and compiled helpers to `$HADAMARD_WORK` (default `<tmp>/hadamard-work`). The certificate does not depend on the search: `verify/verify.py` rebuilds and checks the matrix from `certificate.json` alone.

## Status

To our knowledge this is the first skew-Hadamard matrix of order 3612.

- Cati and Pasechnik, *A database of constructions of Hadamard matrices*, arXiv:2411.18897v2 (Aug. 2025), Table 5,
  lists 903(4). So the smallest known skew order of the form 2ᵐ · 903 was 14448, and 3612 was unknown.
- The Đoković–Seberry census (May 2025) lists 903 among its 737 undecided values.
- A literature search in October 2026 found no later construction.

## Verification

```
python3 verify/verify.py skew/3612
gzip -dc skew/3612/H3612.csv.gz | ./check_matrix
```

SHA-256 of `H3612.csv`:
`1086d13670641326f75919d40212474d6f6766deb7f53201a12993df28f06695`

**Novelty checklist** (2026-10-13):

1. Table 5: 903(4); census: undecided, and in none of its skew series (Paley, Szekeres, Szekeres–Whiteman, Spence) or sporadic values.
2. SageMath's `unknown_skew_hadamard` list stops at 1200.
3. Web searches for the order, and the 2026 deposits listed in the README (Schumer, Alpöge et al., JD-Jones-ASES, Schneider), give nothing.
4. No infinite family reaches it:
   - Paley I: 3611 = 23 · 157;
   - Szekeres, Whiteman, Spence and Momihara–Xiang 2(q + 1): q = 1805 = 5 · 19²;
   - Wallis–Whiteman 4(q + 1), with q = 902;
   - Spence 4(q² + q + 1), (p − 1)ᵘ + 1, (|t| + 1)(q + 1), mn(n − 1), products and skew × amicable products.

This is a finding of nothing, not a proof that no earlier construction exists.
