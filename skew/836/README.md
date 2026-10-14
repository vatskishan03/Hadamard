# Skew-Hadamard matrix of order 836

836 = 4 · 209, with 209 = 11 · 19.

## Construction

A skew-Hadamard matrix is a Hadamard matrix H with H + Hᵀ = 2I.

**Difference family.** Four ±1 sequences x₁, …, x₄ on ℤ₂₀₉ satisfy two conditions:

- each is constant on the orbits of the multiplier subgroup
  H = ⟨26⟩ = {1, 20, 26, 45, 49, 58, 64, 102, 115, 125, 144, 159, 163, 191, 201} ≤ ℤ₂₀₉*, of order 15;
- their periodic autocorrelations sum to zero at every nonzero shift.

`certificate.json` lists, for each sequence, one representative r of each orbit rH on which the sequence is −1.

**Skew first block.** The first sequence is of skew type: x₁(0) = 1 and x₁(−t) = −x₁(t) for t ≠ 0. This is possible
because −1 ∉ H, so negation pairs the cosets.

**Row sums.** The row sums are 1, 9, 15, 23, and 1² + 9² + 15² + 23² = 836.

**The matrix.** Let A, B, C, D be the circulant matrices with these first rows (`sequences.txt`) and R the
back-diagonal identity. The Goethals–Seidel array

```
    [  A     BR     CR     DR  ]
H = [ -BR    A     DᵀR   -CᵀR ]
    [ -CR   -DᵀR    A     BᵀR ]
    [ -DR    CᵀR  -BᵀR    A   ]
```

is then a Hadamard matrix of order 836. Because A + Aᵀ = 2I, it satisfies H + Hᵀ = 2I (Wallis–Whiteman).

## Search

1. `search/skew_triage.py` ranks skew Goethals–Seidel multiplier structures with a Gaussian model. For H = ⟨26⟩ it
   predicts about 3 solutions:
   - 21 orbits and 10 autocorrelation classes;
   - 1024 skew first blocks;
   - pools of 2·10⁴ to 4·10⁴ sequences for the other three blocks.
2. `search/skew_exhaustive.py` enumerates every admissible row-sum pattern and matches the four lists exactly.
3. The pattern (1; 9, 15, 23) gives 288 solutions, and the other two patterns give none.

The whole search takes under two minutes on four cores.

**Reproduction.**

| | |
|---|---|
| Code | commit `0cdf392` (the options used behave the same at later commits) |
| Command | `python3 search/skew_exhaustive.py 209 26` |
| Randomness and caps | exact (no sampling, no seed); pair-work cap 4·10¹⁰ |
| Result | pattern (1; 9, 15, 23): 288 solutions; the other two patterns: none |

## Status

**New** as a skew-Hadamard matrix. Cati–Pasechnik Table 5 lists 209(4): the smallest known skew order 2ᵐ · 209 was
3344, so 836 was unknown.

- 836 is inside the range (≤ 1208) that Cati and Pasechnik checked against implemented constructions.
- The Đoković–Seberry census lists 209 among its 737 undecided v. Koukouvinos–Stylianou (2008, Tables 1 and 3) and
  Đoković's papers of 1992–2024 do not construct 836.
- Đoković's 2007 list of odd n < 300 without a known skew-Hadamard matrix of order 4n includes 209.
- SageMath lists 836 in `unknown_skew_hadamard`.
- No infinite family reaches it: we checked Paley I, Szekeres 2(q + 1), Wallis–Whiteman, Spence, Momihara–Xiang, the
  product theorem n(pᵏ + 1) and amicable products.

The sources searched are listed in [NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | the four sequences as orbit representatives, H, the row sums |
| `sequences.txt` | the four sequences, one per line |
| `H836.csv.gz` | the matrix |

## Verification

```
python3 verify/verify.py skew/836
gzip -dc skew/836/H836.csv.gz | ./check_matrix
```

SHA-256 of `H836.csv`: `952b2e48359cb1890c1dc7712ecb0b20187a2af495d097aa8fbb644a6418aeb8`
