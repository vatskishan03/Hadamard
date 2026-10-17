# Hadamard matrix of order 27524

27524 = 4 · 6881, with 6881 = 7 · 983.

## Construction

Paley's skew-Hadamard matrix of order 3932 (q = 3931, a prime ≡ 3 mod 8) is a 4 × 4 array of circulant blocks of
order s = 983 once its points are listed along the four orbits of a torus element of PSL(2, 3931) and rescaled by
signs. Multiplying it by T-matrices of order 7 (Seberry–Yamada 1990, Theorem 4) gives four circulant ±1 sequences of
length 6881 with Σ PAF = 0 off 0, and the Goethals–Seidel array gives the Hadamard matrix. The construction and its
proof are in [PALEY_BLOCKS.md](../PALEY_BLOCKS.md).

- Torus element: x ↦ −1/(x + τ) with τ = 4; sign normalisation c = −1.
- T-matrices of order 7 (from base sequences BS(4, 3), built by the program): first rows `+000000`, `0++-000`,
  `0000+0+`, `00000+0`.
- Row sums of the four sequences: −13, 45, 61, −147 ((−13)² + 45² + 61² + (−147)² = 27524).

## Search

None: the matrix is determined by q and the T-matrices. `python3 search/paley_blocks.py certify 3931 7` writes the
certificate.

## Status

**New, beyond the tables, but obtainable classically** ([PALEY_BLOCKS.md](../PALEY_BLOCKS.md#status)). In Table 1.53
of the Handbook of Combinatorial Designs (2nd ed., 2007), n = 6881 has t = 6: the smallest known Hadamard matrix with
odd part 6881 had order 2⁶ · 6881. The order lies above the range of Cati–Pasechnik Table 4, which stops at n = 2999.

- The entry does not appear in Đoković's updates (arXiv:0912.5091, arXiv:1008.2043).
- It follows from two published theorems that we found combined nowhere: Paley's matrix of order 3932 is equivalent
  to an array of sixteen circulant blocks of order 983 (Schmidt 1999, after Ito; Yamada 1991), and Seberry and Yamada
  (1990, Theorem 4) multiply such an array by T-matrices of order 7.
- Paley's constructions do not apply: 4n − 1 = 27523 = 17 · 1619 and 2n − 1 = 13761 = 3² · 11 · 139 are not prime
  powers.
- Miyamoto and Yamada–Kiyasu do not apply: 6881 is not a prime power, and neither is n − 2 = 6879 = 3 · 2293.
- Cooper–Wallis does not apply: for each factorization n = a · b it needs T-matrices of order a and Williamson-type
  matrices of order b, and no such pair is known to us.
- SageMath's development version (October 17, 2026) does not construct it: `hadamard_matrix(27524, existence=True)`
  returns Unknown.

We found no construction or statement of order 27524. The sources searched are listed in
[NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, τ, power, c, the T-matrices, the four sequences and the SHA-256 of the matrix |
| `sequences.txt` | the four sequences of length 6881 |

The matrix itself is not stored; `verify/verify.py` rebuilds it.

## Verification

```
python3 verify/verify.py 27524
```

SHA-256 of `H27524.csv`: `5aab8a0e2e82dd96bb24843fb48efd6e9f1b1d4a1205269373b01300d5a7de65`
