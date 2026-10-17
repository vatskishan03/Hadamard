# Hadamard matrix of order 35236

35236 = 4 · 8809, with 8809 = 23 · 383.

## Construction

Paley's skew-Hadamard matrix of order 1532 (q = 1531, a prime ≡ 3 mod 8) is a 4 × 4 array of circulant blocks of
order s = 383 once its points are listed along the four orbits of a torus element of PSL(2, 1531) and rescaled by
signs. Multiplying it by T-matrices of order 23 (Seberry–Yamada 1990, Theorem 4) gives four circulant ±1 sequences of
length 8809 with Σ PAF = 0 off 0, and the Goethals–Seidel array gives the Hadamard matrix. The construction and its
proof are in [PALEY_BLOCKS.md](../PALEY_BLOCKS.md).

- Torus element: x ↦ −1/(x + τ) with τ = 4; sign normalisation c = −1.
- T-matrices of order 23 ([`t-matrices/T23.txt`](../../t-matrices/T23.txt), T-sequences of length 23 as stored in
  SageMath's `t_sequences` module): first rows `+0000000000000000000000`, `0+++++--+-+-00000000000`,
  `000000000000+0+0-0+0-0+`, `0000000000000+0+0-0+0+0`.
- Row sums of the four sequences: 157, −25, −91, 41 (157² + (−25)² + (−91)² + 41² = 35236).

## Search

None: the matrix is determined by q and the T-matrices. `python3 search/paley_blocks.py certify 1531 23
t-matrices/T23.txt` writes the certificate.

## Status

**New, beyond the tables, but obtainable classically** ([PALEY_BLOCKS.md](../PALEY_BLOCKS.md#status)). In Table 1.53
of the Handbook of Combinatorial Designs (2nd ed., 2007), n = 8809 has t = 3: the smallest known Hadamard matrix with
odd part 8809 had order 2³ · 8809. The order lies above the range of Cati–Pasechnik Table 4, which stops at n = 2999.

- The entry does not appear in Đoković's updates (arXiv:0912.5091, arXiv:1008.2043).
- It follows from two published theorems that we found combined nowhere: Paley's matrix of order 1532 is equivalent
  to an array of sixteen circulant blocks of order 383 (Schmidt 1999, after Ito; Yamada 1991), and Seberry and Yamada
  (1990, Theorem 4) multiply such an array by T-matrices of order 23.
- Paley's constructions do not apply: 4n − 1 = 35235 = 3⁵ · 5 · 29 and 2n − 1 = 17617 = 79 · 223 are not prime
  powers.
- Miyamoto and Yamada–Kiyasu do not apply: 8809 is not a prime power, and n − 2 = 8807 is ≡ 7 (mod 8).
- Cooper–Wallis does not apply: for each factorization n = a · b it needs T-matrices of order a and Williamson-type
  matrices of order b, and no such pair is known to us.
- SageMath's development version (October 17, 2026) does not construct it: `hadamard_matrix(35236, existence=True)`
  returns Unknown.

We found no construction or statement of order 35236. The sources searched are listed in
[NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, τ, power, c, the T-matrices, the four sequences and the SHA-256 of the matrix |
| `sequences.txt` | the four sequences of length 8809 |

The matrix itself is not stored; `verify/verify.py` rebuilds it.

## Verification

```
python3 verify/verify.py 35236
```

SHA-256 of `H35236.csv`: `8c5dda8990a24a1467de7ebe0f83ee534d17e3bcf8e8a26b998b441552d16095`
