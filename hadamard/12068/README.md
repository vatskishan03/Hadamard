# Hadamard matrix of order 12068

12068 = 4 · 3017, with 3017 = 7 · 431.

## Construction

Paley's skew-Hadamard matrix of order 1724 (q = 1723, a prime ≡ 3 mod 8) is a 4 × 4 array of circulant blocks of
order s = 431 once its points are listed along the four orbits of a torus element of PSL(2, 1723) and rescaled by
signs. Multiplying it by T-matrices of order 7 (Seberry–Yamada 1990, Theorem 4) gives four circulant ±1 sequences of
length 3017 with Σ PAF = 0 off 0, and the Goethals–Seidel array gives the Hadamard matrix. The construction and its
proof are in [PALEY_BLOCKS.md](../PALEY_BLOCKS.md).

- Torus element: the square of x ↦ −1/(x + τ) with τ = 3; sign normalisation c = +1.
- T-matrices of order 7 (from base sequences BS(4, 3), built by the program): first rows `+000000`, `0++-000`,
  `0000+0+`, `00000+0`.
- Row sums of the four sequences: −57, −7, −11, 93 ((−57)² + (−7)² + (−11)² + 93² = 12068).

## Search

None: the matrix is determined by q and the T-matrices. `python3 search/paley_blocks.py certify 1723 7` writes the
certificate.

## Status

**New, beyond the tables, but obtainable classically** ([PALEY_BLOCKS.md](../PALEY_BLOCKS.md#status)). In Table 1.53
of the Handbook of Combinatorial Designs (2nd ed., 2007), n = 3017 has t = 3: the smallest known Hadamard matrix with
odd part 3017 had order 2³ · 3017. The order lies above the range of Cati–Pasechnik Table 4, which stops at n = 2999.

- The entry does not appear in Đoković's updates (arXiv:0912.5091, arXiv:1008.2043).
- It follows from two published theorems that we found combined nowhere: Paley's matrix of order 1724 is equivalent
  to an array of sixteen circulant blocks of order 431 (Schmidt 1999, after Ito; Yamada 1991), and Seberry and Yamada
  (1990, Theorem 4) multiply such an array by T-matrices of order 7.
- Paley's constructions do not apply: 4n − 1 = 12067 = 11 · 1097 and 2n − 1 = 6033 = 3 · 2011 are not prime powers.
- Miyamoto and Yamada–Kiyasu do not apply: 3017 is not a prime power, and neither is n − 2 = 3015 = 3² · 5 · 67.
- Cooper–Wallis does not apply: for each factorization n = a · b it needs T-matrices of order a and Williamson-type
  matrices of order b, and no such pair is known to us.
- SageMath's development version (October 17, 2026) does not construct it: `hadamard_matrix(12068, existence=True)`
  returns Unknown.

We found no construction or statement of order 12068. The sources searched are listed in
[NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, τ, power, c, the T-matrices, the four sequences and the SHA-256 of the matrix |
| `sequences.txt` | the four sequences of length 3017 |

The matrix itself is not stored; `verify/verify.py` rebuilds it.

## Verification

```
python3 verify/verify.py 12068
```

SHA-256 of `H12068.csv`: `3561359dcb52bd2608505b712d56b38af269d68f29f2b95ed0db7ff33a3ae9e1`
