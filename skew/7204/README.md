# Skew-Hadamard matrix of order 7204

7204 = 4 · 1801 (prime), and 1801 = 3 · 24² + 3 · 24 + 1.

## Construction

A skew-Hadamard matrix is a Hadamard matrix H with H + Hᵀ = 2I.

**Difference family.** Four ±1 sequences x₁, …, x₄ on ℤ₁₈₀₁ satisfy two conditions:

- each is constant on the orbits of the multiplier subgroup H = ⟨2⟩ ≤ ℤ₁₈₀₁*, of order 25, the group of 72nd powers
  modulo 1801: {1, 2, 4, 8, 16, 32, 64, 128, 175, 197, 247, 256, 350, 394, 494, 512, 700, 788, 901, 988, 999, 1024,
  1351, 1400, 1576};
- their periodic autocorrelations sum to zero at every nonzero shift.

H has 72 orbits on ℤ₁₈₀₁ ∖ {0}, the cosets of H in ℤ₁₈₀₁*. `certificate.json` lists, for each sequence, one
representative r of each orbit rH on which the sequence is −1 (0 is listed when the sequence is −1 at 0).

**Skew first block.** The first sequence is of skew type: x₁(0) = 1 and x₁(−t) = −x₁(t) for t ≠ 0. This is possible
because −1 ∉ H (H has odd order), so negation pairs the orbits.

**Spin multiplier.** μ = 21 is not in H, and μ³ ∈ H. It fixes x₁, and the other blocks satisfy x₃(t) = x₂(μt) and
x₄(t) = x₂(μ²t) (so also x₂(t) = x₄(μt), since x₂ is H-invariant). This is the spin structure of Đoković
(arXiv:2404.14375), as in [skew/3676](../3676).

**Row sums.** The row sums are 1, 49, 49, 49, and 1 + 3 · 49² = 7204. A sequence constant on the orbits of H has row sum
x(0) + 25 · (an even number), that is ≡ ±1 (mod 50), and (1, 49, 49, 49) is the only solution of 4v = 1 + b² + c² + d²
in such numbers, up to order and sign.

**The matrix.** Let A, B, C, D be the circulant matrices with these first rows (`sequences.txt`) and R the
back-diagonal identity. The Goethals–Seidel array

```
    [  A     BR     CR     DR  ]
H = [ -BR    A     DᵀR   -CᵀR ]
    [ -CR   -DᵀR    A     BᵀR ]
    [ -DR    CᵀR  -BᵀR    A   ]
```

is then a Hadamard matrix of order 7204. Because A + Aᵀ = 2I, it satisfies H + Hᵀ = 2I (Wallis–Whiteman).

## Search

1. **Structure.** The family is complementary if and only if P_A(s) + P_B(s) + P_B(μs) + P_B(μ²s) = 0 on the 12
   orbits of ⟨H, μ, −1⟩ on the nonzero shifts. `search/spin_skew.py rank` predicts about 10^5.1 solutions for this
   structure ([search/NOTES.md](../../search/NOTES.md) §14).
2. **Tabu search** (`search/spin_ls2.c`, the vectorised form of `search/spin_ls.c`; [search/NOTES.md](../../search/NOTES.md)
   §33 and §115). The moves are sign changes of A on a pair of orbits {O, −O} of ⟨H, μ⟩ and swaps of two B-orbits of
   opposite sign, which keep the row sum; the objective is the weighted sum of the squared residuals; restarts of
   3000 steps.
3. The census scheduler ran the search with six workers (each `spin_ls2 input.bin 5000 3000 <seed>`, distinct seeds) on
   a 32-core machine from 1 November 2026. A solution was found on 2 November 2026 at 10:50 UTC, after 102.5 core-hours
   and 1.66 · 10⁷ restarts; it was verified exactly and certified (`search/skew_certify.py`).

`certificate.json` records μ as `spin_multiplier`, and `verify/verify.py` checks that the last three blocks are its
images.

## Status

**New** as a skew-Hadamard matrix.

- The Đoković–Seberry census lists 1801 among its 737 undecided odd v < 2500 for skew matrices.
- 1801 lies beyond Cati–Pasechnik's Table 5 (v ≤ 999) and the other skew tables we read.
- Đoković (arXiv:2404.14375) gives spin families up to v = 631, and Schumer's 919 (Zenodo, 2026) is the largest
  published one; neither gives 7204.

The sources searched are listed in [NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | the four sequences as orbit representatives, H, the spin multiplier, the row sums |
| `sequences.txt` | the four sequences, one per line |
| `H7204.csv.gz` | the matrix |

## Verification

```
python3 verify/verify.py skew/7204
gzip -dc skew/7204/H7204.csv.gz | ./check_matrix
```

SHA-256 of `H7204.csv`: `e73658c79c88839461fb178131b716fa85a21e090e2ed8d635e36c4959a19094`
