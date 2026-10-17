# Hadamard matrix of order 17740

17740 = 4 · 4435, with 4435 = 5 · 887.

## Construction

Paley's skew-Hadamard matrix of order 3548 (q = 3547, a prime ≡ 3 mod 8) is a 4 × 4 array of circulant blocks of
order s = 887 once its points are listed along the four orbits of a torus element of PSL(2, 3547) and rescaled by
signs. Multiplying it by T-matrices of order 5 (Seberry–Yamada 1990, Theorem 4) gives four circulant ±1 sequences of
length 4435 with Σ PAF = 0 off 0, and the Goethals–Seidel array gives the Hadamard matrix. The construction and its
proof are in [PALEY_BLOCKS.md](../PALEY_BLOCKS.md).

- Torus element: the square of x ↦ −1/(x + τ) with τ = 3; sign normalisation c = +1.
- T-matrices of order 5 (from base sequences BS(3, 2), built by the program): first rows `++000`, `00+00`, `000+-`,
  `00000`.
- Row sums of the four sequences: −31, 67, −103, −41 ((−31)² + 67² + (−103)² + (−41)² = 17740).

## Search

None: the matrix is determined by q and the T-matrices. `python3 search/paley_blocks.py certify 3547 5` writes the
certificate.

## Status

**New, beyond the tables, but obtainable classically** ([PALEY_BLOCKS.md](../PALEY_BLOCKS.md#status)). In Table 1.53
of the Handbook of Combinatorial Designs (2nd ed., 2007), n = 4435 has t = 3: the smallest known Hadamard matrix with
odd part 4435 had order 2³ · 4435. The order lies above the range of Cati–Pasechnik Table 4, which stops at n = 2999.

- The entry does not appear in Đoković's updates (arXiv:0912.5091, arXiv:1008.2043).
- It follows from two published theorems that we found combined nowhere: Paley's matrix of order 3548 is equivalent
  to an array of sixteen circulant blocks of order 887 (Schmidt 1999, after Ito; Yamada 1991), and Seberry and Yamada
  (1990, Theorem 4) multiply such an array by T-matrices of order 5.
- Paley's constructions do not apply: 4n − 1 = 17739 = 3⁵ · 73 and 2n − 1 = 8869 = 7² · 181 are not prime powers.
- Miyamoto and Yamada–Kiyasu do not apply: 4435 is not a prime power, and neither is n − 2 = 4433 = 11 · 13 · 31.
- Cooper–Wallis does not apply: for each factorization n = a · b it needs T-matrices of order a and Williamson-type
  matrices of order b, and no such pair is known to us.
- SageMath's development version (October 17, 2026) does not construct it: `hadamard_matrix(17740, existence=True)`
  returns Unknown.

We found no construction or statement of order 17740. The sources searched are listed in
[NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, τ, power, c, the T-matrices, the four sequences and the SHA-256 of the matrix |
| `sequences.txt` | the four sequences of length 4435 |

The matrix itself is not stored; `verify/verify.py` rebuilds it.

## Verification

```
python3 verify/verify.py 17740
```

SHA-256 of `H17740.csv`: `024203a68c57dfc75bdf5002c4e5c6d4d52c4079fdc2185a096f7fe6b7d53229`
