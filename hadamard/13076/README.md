# Hadamard matrix of order 13076

13076 = 4 · 3269, with 3269 = 7 · 467.

## Construction

Paley's skew-Hadamard matrix of order 1868 (q = 1867, a prime ≡ 3 mod 8) is a 4 × 4 array of circulant blocks of
order s = 467 once its points are listed along the four orbits of a torus element of PSL(2, 1867) and rescaled by
signs. Multiplying it by T-matrices of order 7 (Seberry–Yamada 1990, Theorem 4) gives four circulant ±1 sequences of
length 3269 with Σ PAF = 0 off 0, and the Goethals–Seidel array gives the Hadamard matrix. The construction and its
proof are in [PALEY_BLOCKS.md](../PALEY_BLOCKS.md).

- Torus element: the square of x ↦ −1/(x + τ) with τ = 3; sign normalisation c = +1.
- T-matrices of order 7 (from base sequences BS(4, 3), built by the program): first rows `+000000`, `0++-000`,
  `0000+0+`, `00000+0`.
- Row sums of the four sequences: 87, −51, −35, 41 (87² + (−51)² + (−35)² + 41² = 13076).

## Search

None: the matrix is determined by q and the T-matrices. `python3 search/paley_blocks.py certify 1867 7` writes the
certificate.

## Status

**New, beyond the tables, but obtainable classically** ([PALEY_BLOCKS.md](../PALEY_BLOCKS.md#status)). In Table 1.53
of the Handbook of Combinatorial Designs (2nd ed., 2007), n = 3269 has t = 3: the smallest known Hadamard matrix with
odd part 3269 had order 2³ · 3269. The order lies above the range of Cati–Pasechnik Table 4, which stops at n = 2999.

- The entry does not appear in Đoković's updates (arXiv:0912.5091, arXiv:1008.2043).
- It follows from two published theorems that we found combined nowhere: Paley's matrix of order 1868 is equivalent
  to an array of sixteen circulant blocks of order 467 (Schmidt 1999, after Ito; Yamada 1991), and Seberry and Yamada
  (1990, Theorem 4) multiply such an array by T-matrices of order 7.
- Paley's constructions do not apply: 4n − 1 = 13075 = 5² · 523 and 2n − 1 = 6537 = 3 · 2179 are not prime powers.
- Miyamoto's theorem does not apply: a symmetric conference matrix of order n + 1 = 3270 cannot exist (3269 is not a sum
  of two squares); Yamada–Kiyasu does not apply: n − 2 = 3267 = 3³ · 11² is not a prime power.
- Cooper–Wallis does not apply: for each factorization n = a · b it needs T-matrices of order a and Williamson-type
  matrices of order b, and no such pair is known to us.
- SageMath's development version (October 17, 2026) does not construct it: `hadamard_matrix(13076, existence=True)`
  returns Unknown.

We found no construction or statement of order 13076. The sources searched are listed in
[NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, τ, power, c, the T-matrices, the four sequences and the SHA-256 of the matrix |
| `sequences.txt` | the four sequences of length 3269 |

The matrix itself is not stored; `verify/verify.py` rebuilds it.

## Verification

```
python3 verify/verify.py 13076
```

SHA-256 of `H13076.csv`: `a17156a57a90c7ac5e1b9cc89765979ff49f05be4084a3caa4d954f9113b9c6d`
