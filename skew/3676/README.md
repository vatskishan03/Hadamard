# Skew-Hadamard matrix of order 3676

3676 = 4 · 919 (prime).

## Construction

A skew-Hadamard matrix is a Hadamard matrix H with H + Hᵀ = 2I.

**Difference family.** Four ±1 sequences x₁, …, x₄ on ℤ₉₁₉ satisfy two conditions:

- each is constant on the orbits of the multiplier subgroup H = ⟨58⟩ ≤ ℤ₉₁₉*, of order 17: {1, 58, 162, 206, 229, 234,
  284, 288, 305, 338, 416, 512, 535, 607, 703, 706, 849};
- their periodic autocorrelations sum to zero at every nonzero shift.

`certificate.json` lists, for each sequence, one representative r of each orbit rH on which the sequence is −1.

**Skew first block.** The first sequence is of skew type: x₁(0) = 1 and x₁(−t) = −x₁(t) for t ≠ 0. This is possible
because −1 ∉ H, so negation pairs the orbits.

**Spin multiplier.** μ = 52 has order 3 modulo 919. It fixes x₁, and the other blocks satisfy x₃(t) = x₂(52t) and
x₄(t) = x₂(52²t).

**Row sums.** The row sums are 1, 35, 35, 35, and 1² + 35² + 35² + 35² = 3676.

**The matrix.** Let A, B, C, D be the circulant matrices with these first rows (`sequences.txt`) and R the
back-diagonal identity. The Goethals–Seidel array

```
    [  A     BR     CR     DR  ]
H = [ -BR    A     DᵀR   -CᵀR ]
    [ -CR   -DᵀR    A     BᵀR ]
    [ -DR    CᵀR  -BᵀR    A   ]
```

is then a Hadamard matrix of order 3676. Because A + Aᵀ = 2I, it satisfies H + Hᵀ = 2I (Wallis–Whiteman).

## Search

1. **Spin structure.** Following Đoković's skew matrices of order 4 · 397:
   - A is invariant under K = ⟨H, μ⟩, with μ = 52 of order 3 and |K| = 51;
   - μ permutes the other blocks: C(t) = B(μt) and D(t) = B(μ²t), so P_C(s) = P_B(μs);
   - the family is complementary iff P_A(s) + P_B(s) + P_B(μs) + P_B(μ²s) = 0 on the 9 orbits of ⟨K, −1⟩ on the
     nonzero shifts.

   The row sums are forced to (1, b, b, b) with 1 + 3b² = 4 · 919, so b = 35. `search/spin_skew.py rank` predicts
   about 2·10⁵ solutions for this structure:
   - 512 skew first blocks;
   - 2·10¹⁵ choices of B.
2. `search/spin_skew.py csearch` runs a tabu search in C (`search/spin_ls.c`) on (A, B). The moves are swaps of
   B-orbits of opposite sign, which keep the row sum, and sign changes of A-orbit pairs. The objective is the weighted
   sum of squared residuals.
3. The search found a solution in 37 s, at restart 399.

`certificate.json` records μ as `spin_multiplier`, and `verify/verify.py` checks that the last three blocks are its
images.

**Reproduction.**

| | |
|---|---|
| Code | commit `42b3ba6` (the options used behave the same at later commits) |
| Command | `python3 search/spin_skew.py csearch 919 58 52 --restarts 5000 --steps 3000 --seed 1` |
| Randomness and caps | tabu search in C (search/spin_ls.c), seed 1, at most 5000 restarts of 3000 moves |
| Result | solution at restart 399, step 1915, after 37 s |
| Certificate | `search/skew_certify.py` turns a logged hit into the files in this directory |

## Status

**Known.** A. Schumer published a skew-Hadamard matrix of order 3676 on Zenodo on 2026-09-17,
doi:10.5281/zenodo.22820997. The matrix here was found independently, by the spin tabu search above, and is a second
explicit witness.

- Schumer's witness has the same spin structure: H = ⟨703⟩ of order 17, μ = 866 of order 3, the first block skew, and
  row sums (1, 35, 35, 35).

The sources searched are listed in [NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | the four sequences as orbit representatives, H, the spin multiplier, the row sums |
| `sequences.txt` | the four sequences, one per line |
| `H3676.csv.gz` | the matrix |

## Verification

```
python3 verify/verify.py skew/3676
gzip -dc skew/3676/H3676.csv.gz | ./check_matrix
```

SHA-256 of `H3676.csv`: `5a1ac2efc4d329f95431703a59f7bfbc38b71984b5132bc7749725d3c6ecaf56`
