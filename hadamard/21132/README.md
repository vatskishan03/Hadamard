# Hadamard matrix of order 21132

21132 = 4 · 5283, with 5283 = 3² · 587.

## Construction

Paley's skew-Hadamard matrix of order 2348 (q = 2347, a prime ≡ 3 mod 8) is a 4 × 4 array of circulant blocks of
order s = 587 once its points are listed along the four orbits of a torus element of PSL(2, 2347) and rescaled by
signs. Multiplying it by T-matrices of order 9 (Seberry–Yamada 1990, Theorem 4) gives four circulant ±1 sequences of
length 5283 with Σ PAF = 0 off 0, and the Goethals–Seidel array gives the Hadamard matrix. The construction and its
proof are in [PALEY_BLOCKS.md](../PALEY_BLOCKS.md).

- Torus element: the square of x ↦ −1/(x + τ) with τ = 3; sign normalisation c = +1.
- T-matrices of order 9 ([`t-matrices/T9.txt`](../../t-matrices/T9.txt), T-sequences of length 9 from a Golay pair of
  length 8 (Turyn 1974)): first rows `+00000000`, `0+++-0000`, `00000++-+`, `000000000`.
- Row sums of the four sequences: −131, 9, 61, −13 ((−131)² + 9² + 61² + (−13)² = 21132).

## Search

None: the matrix is determined by q and the T-matrices. `python3 search/paley_blocks.py certify 2347 9
t-matrices/T9.txt` writes the certificate.

## Status

**New, beyond the tables, but obtainable classically** ([PALEY_BLOCKS.md](../PALEY_BLOCKS.md#status)). In Table 1.53
of the Handbook of Combinatorial Designs (2nd ed., 2007), n = 5283 has t = 3: the smallest known Hadamard matrix with
odd part 5283 had order 2³ · 5283. The order lies above the range of Cati–Pasechnik Table 4, which stops at n = 2999.

- The entry does not appear in Đoković's updates (arXiv:0912.5091, arXiv:1008.2043).
- It follows from two published theorems that we found combined nowhere: Paley's matrix of order 2348 is equivalent
  to an array of sixteen circulant blocks of order 587 (Schmidt 1999, after Ito; Yamada 1991), and Seberry and Yamada
  (1990, Theorem 4) multiply such an array by T-matrices of order 9.
- Paley's constructions do not apply: 4n − 1 = 21131 = 11 · 17 · 113 and 2n − 1 = 10565 = 5 · 2113 are not prime
  powers.
- Miyamoto's theorem does not apply: it needs a symmetric conference matrix of order n + 1, which cannot exist for n ≡ 3
  (mod 4). Yamada's theorem with q = n − 2 = 5281 needs a symmetric conference matrix of order 2642; 2641 = 19 · 139 is
  not a sum of two squares, so no conference matrix of order 2642 exists.
- Cooper–Wallis does not apply: for each factorization n = a · b it needs T-matrices of order a and Williamson-type
  matrices of order b, and no such pair is known to us.
- SageMath's development version (17 October 2026) does not construct it: `hadamard_matrix(21132, existence=True)`
  returns Unknown.

We found no construction or statement of order 21132. The sources searched are listed in
[NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, τ, power, c, the T-matrices, the four sequences and the SHA-256 of the matrix |
| `sequences.txt` | the four sequences of length 5283 |

The matrix itself is not stored; `verify/verify.py` rebuilds it.

## Verification

```
python3 verify/verify.py 21132
```

SHA-256 of `H21132.csv`: `294f25ffa090f6b8242fe8c98c7dc06173a92e39eae74e51d9889b9f736016f9`
