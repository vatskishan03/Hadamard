# Hadamard matrix of order 23156

23156 = 4 · 5789, with 5789 = 7 · 827.

## Construction

Paley's skew-Hadamard matrix of order 3308 (q = 3307, a prime ≡ 3 mod 8) is a 4 × 4 array of circulant blocks of
order s = 827 once its points are listed along the four orbits of a torus element of PSL(2, 3307) and rescaled by
signs. Multiplying it by T-matrices of order 7 (Seberry–Yamada 1990, Theorem 4) gives four circulant ±1 sequences of
length 5789 with Σ PAF = 0 off 0, and the Goethals–Seidel array gives the Hadamard matrix. The construction and its
proof are in [PALEY_BLOCKS.md](../PALEY_BLOCKS.md).

- Torus element: the square of x ↦ −1/(x + τ) with τ = 3; sign normalisation c = +1.
- T-matrices of order 7 (from base sequences BS(4, 3), built by the program): first rows `+000000`, `0++-000`,
  `0000+0+`, `00000+0`.
- Row sums of the four sequences: 125, 33, −41, −69 (125² + 33² + (−41)² + (−69)² = 23156).

## Search

None: the matrix is determined by q and the T-matrices. `python3 search/paley_blocks.py certify 3307 7` writes the
certificate.

## Status

**New, beyond the tables, but obtainable classically** ([PALEY_BLOCKS.md](../PALEY_BLOCKS.md#status)). In Table 1.53
of the Handbook of Combinatorial Designs (2nd ed., 2007), n = 5789 has t = 3: the smallest known Hadamard matrix with
odd part 5789 had order 2³ · 5789. The order lies above the range of Cati–Pasechnik Table 4, which stops at n = 2999.

- The entry does not appear in Đoković's updates (arXiv:0912.5091, arXiv:1008.2043).
- It follows from two published theorems that we found combined nowhere: Paley's matrix of order 3308 is equivalent
  to an array of sixteen circulant blocks of order 827 (Schmidt 1999, after Ito; Yamada 1991), and Seberry and Yamada
  (1990, Theorem 4) multiply such an array by T-matrices of order 7.
- Paley's constructions do not apply: 4n − 1 = 23155 = 5 · 11 · 421 and 2n − 1 = 11577 = 3 · 17 · 227 are not prime
  powers.
- Miyamoto's theorem does not apply: a symmetric conference matrix of order n + 1 = 5790 cannot exist (5789 is not a sum
  of two squares); Yamada–Kiyasu does not apply: n − 2 = 5787 = 3² · 643 is not a prime power.
- Cooper–Wallis does not apply: for each factorization n = a · b it needs T-matrices of order a and Williamson-type
  matrices of order b, and no such pair is known to us.
- SageMath's development version (17 October 2026) does not construct it: `hadamard_matrix(23156, existence=True)`
  returns Unknown.

We found no construction or statement of order 23156. The sources searched are listed in
[NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, τ, power, c, the T-matrices, the four sequences and the SHA-256 of the matrix |
| `sequences.txt` | the four sequences of length 5789 |

The matrix itself is not stored; `verify/verify.py` rebuilds it.

## Verification

```
python3 verify/verify.py 23156
```

SHA-256 of `H23156.csv`: `0bdf271ce07463b10951c084e9d5a4cd34b9e3ab112fa5ee60c442b38727ab05`
