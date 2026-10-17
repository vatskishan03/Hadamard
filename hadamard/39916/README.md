# Hadamard matrix of order 39916

39916 = 4 · 9979, with 9979 = 17 · 587.

## Construction

Paley's skew-Hadamard matrix of order 2348 (q = 2347, a prime ≡ 3 mod 8) is a 4 × 4 array of circulant blocks of
order s = 587 once its points are listed along the four orbits of a torus element of PSL(2, 2347) and rescaled by
signs. Multiplying it by T-matrices of order 17 (Seberry–Yamada 1990, Theorem 4) gives four circulant ±1 sequences of
length 9979 with Σ PAF = 0 off 0, and the Goethals–Seidel array gives the Hadamard matrix. The construction and its
proof are in [PALEY_BLOCKS.md](../PALEY_BLOCKS.md).

- Torus element: the square of x ↦ −1/(x + τ) with τ = 3; sign normalisation c = +1.
- T-matrices of order 17 ([`t-matrices/T17.txt`](../../t-matrices/T17.txt), T-sequences of length 17 from a Golay
  pair of length 16 (Turyn 1974)): first rows `+0000000000000000`, `0+++-++-+00000000`, `000000000+++---+-`,
  `00000000000000000`.
- Row sums of the four sequences: −131, 37, 85, 119 ((−131)² + 37² + 85² + 119² = 39916).

## Search

None: the matrix is determined by q and the T-matrices. `python3 search/paley_blocks.py certify 2347 17
t-matrices/T17.txt` writes the certificate.

## Status

**New, beyond the tables, but obtainable classically** ([PALEY_BLOCKS.md](../PALEY_BLOCKS.md#status)). In Table 1.53
of the Handbook of Combinatorial Designs (2nd ed., 2007), n = 9979 has t = 3: the smallest known Hadamard matrix with
odd part 9979 had order 2³ · 9979. The order lies above the range of Cati–Pasechnik Table 4, which stops at n = 2999.

- The entry does not appear in Đoković's updates (arXiv:0912.5091, arXiv:1008.2043).
- It follows from two published theorems that we found combined nowhere: Paley's matrix of order 2348 is equivalent
  to an array of sixteen circulant blocks of order 587 (Schmidt 1999, after Ito; Yamada 1991), and Seberry and Yamada
  (1990, Theorem 4) multiply such an array by T-matrices of order 17.
- Paley's constructions do not apply: 4n − 1 = 39915 = 3² · 5 · 887 and 2n − 1 = 19957 = 7 · 2851 are not prime
  powers.
- Miyamoto and Yamada–Kiyasu do not apply: 9979 is not a prime power, and neither is n − 2 = 9977 = 11 · 907.
- Cooper–Wallis does not apply: for each factorization n = a · b it needs T-matrices of order a and Williamson-type
  matrices of order b, and no such pair is known to us.
- SageMath's development version (October 17, 2026) does not construct it: `hadamard_matrix(39916, existence=True)`
  returns Unknown.

We found no construction or statement of order 39916. The sources searched are listed in
[NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, τ, power, c, the T-matrices, the four sequences and the SHA-256 of the matrix |
| `sequences.txt` | the four sequences of length 9979 |

The matrix itself is not stored; `verify/verify.py` rebuilds it.

## Verification

```
python3 verify/verify.py 39916
```

SHA-256 of `H39916.csv`: `855909ecec0ed0d9302b0d8efd1b82c9fc982e993e2e58b6ee3341ec554ce142`
