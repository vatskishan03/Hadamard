# Hadamard matrix of order 25828

25828 = 4 · 6457, with 6457 = 11 · 587.

## Construction

Paley's skew-Hadamard matrix of order 2348 (q = 2347, a prime ≡ 3 mod 8) is a 4 × 4 array of circulant blocks of
order s = 587 once its points are listed along the four orbits of a torus element of PSL(2, 2347) and rescaled by
signs. Multiplying it by T-matrices of order 11 (Seberry–Yamada 1990, Theorem 4) gives four circulant ±1 sequences of
length 6457 with Σ PAF = 0 off 0, and the Goethals–Seidel array gives the Hadamard matrix. The construction and its
proof are in [PALEY_BLOCKS.md](../PALEY_BLOCKS.md).

- Torus element: the square of x ↦ −1/(x + τ) with τ = 3; sign normalisation c = +1.
- T-matrices of order 11 ([`t-matrices/T11.txt`](../../t-matrices/T11.txt), T-sequences of length 11 from a Golay
  pair of length 10 (Turyn 1974)): first rows `+0000000000`, `0++00+0+-00`, `000++0-00-+`, `00000000000`.
- Row sums of the four sequences: −131, 23, 73, 53 ((−131)² + 23² + 73² + 53² = 25828).

## Search

None: the matrix is determined by q and the T-matrices. `python3 search/paley_blocks.py certify 2347 11
t-matrices/T11.txt` writes the certificate.

## Status

**New, beyond the tables, but obtainable classically** ([PALEY_BLOCKS.md](../PALEY_BLOCKS.md#status)). In Table 1.53
of the Handbook of Combinatorial Designs (2nd ed., 2007), n = 6457 has t = 5: the smallest known Hadamard matrix with
odd part 6457 had order 2⁵ · 6457. The order lies above the range of Cati–Pasechnik Table 4, which stops at n = 2999.

- The entry does not appear in Đoković's updates (arXiv:0912.5091, arXiv:1008.2043).
- It follows from two published theorems that we found combined nowhere: Paley's matrix of order 2348 is equivalent
  to an array of sixteen circulant blocks of order 587 (Schmidt 1999, after Ito; Yamada 1991), and Seberry and Yamada
  (1990, Theorem 4) multiply such an array by T-matrices of order 11.
- Paley's constructions do not apply: 4n − 1 = 25827 = 3 · 8609 and 2n − 1 = 12913 = 37 · 349 are not prime powers.
- Miyamoto's theorem does not apply: a symmetric conference matrix of order n + 1 = 6458 cannot exist (6457 is not a sum
  of two squares); Yamada–Kiyasu does not apply: n − 2 = 6455 = 5 · 1291 is not a prime power.
- Cooper–Wallis does not apply: for each factorization n = a · b it needs T-matrices of order a and Williamson-type
  matrices of order b, and no such pair is known to us.
- SageMath's development version (17 October 2026) does not construct it: `hadamard_matrix(25828, existence=True)`
  returns Unknown.

We found no construction or statement of order 25828. The sources searched are listed in
[NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, τ, power, c, the T-matrices, the four sequences and the SHA-256 of the matrix |
| `sequences.txt` | the four sequences of length 6457 |

The matrix itself is not stored; `verify/verify.py` rebuilds it.

## Verification

```
python3 verify/verify.py 25828
```

SHA-256 of `H25828.csv`: `9bb6b3a96f0319e32e4780f13cdc915473ad9e8d6d2e7a55791610ae76931e80`
