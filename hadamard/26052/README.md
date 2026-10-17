# Hadamard matrix of order 26052

26052 = 4 · 6513, with 6513 = 3 · 13 · 167.

## Construction

Paley's skew-Hadamard matrix of order 2004 (q = 2003, a prime ≡ 3 mod 8) is a 4 × 4 array of circulant blocks of
order s = 501 once its points are listed along the four orbits of a torus element of PSL(2, 2003) and rescaled by
signs. Multiplying it by T-matrices of order 13 (Seberry–Yamada 1990, Theorem 4) gives four circulant ±1 sequences of
length 6513 with Σ PAF = 0 off 0, and the Goethals–Seidel array gives the Hadamard matrix. The construction and its
proof are in [PALEY_BLOCKS.md](../PALEY_BLOCKS.md).

- Torus element: the square of x ↦ −1/(x + τ) with τ = 3; sign normalisation c = +1.
- T-matrices of order 13 (from base sequences BS(7, 6), built by the program): first rows `+0+0+00000000`,
  `0+0+0-+000000`, `0000000+0000-`, `00000000++--0`.
- Row sums of the four sequences: −59, 95, 75, −89 ((−59)² + 95² + 75² + (−89)² = 26052).

## Search

None: the matrix is determined by q and the T-matrices. `python3 search/paley_blocks.py certify 2003 13` writes the
certificate.

## Status

**New, beyond the tables, but obtainable classically** ([PALEY_BLOCKS.md](../PALEY_BLOCKS.md#status)). In Table 1.53
of the Handbook of Combinatorial Designs (2nd ed., 2007), n = 6513 has t = 3: the smallest known Hadamard matrix with
odd part 6513 had order 2³ · 6513. The order lies above the range of Cati–Pasechnik Table 4, which stops at n = 2999.

- The entry does not appear in Đoković's updates (arXiv:0912.5091, arXiv:1008.2043).
- It follows from two published theorems that we found combined nowhere: Paley's matrix of order 2004 is equivalent
  to an array of sixteen circulant blocks of order 501 (Schmidt 1999, after Ito; Yamada 1991), and Seberry and Yamada
  (1990, Theorem 4) multiply such an array by T-matrices of order 13.
- Paley's constructions do not apply: 4n − 1 = 26051 = 109 · 239 and 2n − 1 = 13025 = 5² · 521 are not prime powers.
- Miyamoto and Yamada–Kiyasu do not apply: 6513 is not a prime power, and neither is n − 2 = 6511 = 17 · 383.
- Cooper–Wallis does not apply: for each factorization n = a · b it needs T-matrices of order a and Williamson-type
  matrices of order b, and no such pair is known to us.
- SageMath's development version (October 17, 2026) does not construct it: `hadamard_matrix(26052, existence=True)`
  returns Unknown.

We found no construction or statement of order 26052. The sources searched are listed in
[NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, τ, power, c, the T-matrices, the four sequences and the SHA-256 of the matrix |
| `sequences.txt` | the four sequences of length 6513 |

The matrix itself is not stored; `verify/verify.py` rebuilds it.

## Verification

```
python3 verify/verify.py 26052
```

SHA-256 of `H26052.csv`: `7b8d99f2043756bb2c67954d714f27c5cf0833bac59770e775f66d1633ce81b3`
