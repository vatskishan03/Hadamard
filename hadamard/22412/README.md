# Hadamard matrix of order 22412

22412 = 4 · 5603, with 5603 = 13 · 431.

## Construction

Paley's skew-Hadamard matrix of order 1724 (q = 1723, a prime ≡ 3 mod 8) is a 4 × 4 array of circulant blocks of
order s = 431 once its points are listed along the four orbits of a torus element of PSL(2, 1723) and rescaled by
signs. Multiplying it by T-matrices of order 13 (Seberry–Yamada 1990, Theorem 4) gives four circulant ±1 sequences of
length 5603 with Σ PAF = 0 off 0, and the Goethals–Seidel array gives the Hadamard matrix. The construction and its
proof are in [PALEY_BLOCKS.md](../PALEY_BLOCKS.md).

- Torus element: the square of x ↦ −1/(x + τ) with τ = 3; sign normalisation c = +1.
- T-matrices of order 13 (from base sequences BS(7, 6), built by the program): first rows `+0+0+00000000`,
  `0+0+0-+000000`, `0000000+0000-`, `00000000++--0`.
- Row sums of the four sequences: 45, −61, 45, 121 (45² + (−61)² + 45² + 121² = 22412).

## Search

None: the matrix is determined by q and the T-matrices. `python3 search/paley_blocks.py certify 1723 13` writes the
certificate.

## Status

**New, beyond the tables, but obtainable classically** ([PALEY_BLOCKS.md](../PALEY_BLOCKS.md#status)). In Table 1.53
of the Handbook of Combinatorial Designs (2nd ed., 2007), n = 5603 has t = 3: the smallest known Hadamard matrix with
odd part 5603 had order 2³ · 5603. The order lies above the range of Cati–Pasechnik Table 4, which stops at n = 2999.

- The entry does not appear in Đoković's updates (arXiv:0912.5091, arXiv:1008.2043).
- It follows from two published theorems that we found combined nowhere: Paley's matrix of order 1724 is equivalent
  to an array of sixteen circulant blocks of order 431 (Schmidt 1999, after Ito; Yamada 1991), and Seberry and Yamada
  (1990, Theorem 4) multiply such an array by T-matrices of order 13.
- Paley's constructions do not apply: 4n − 1 = 22411 = 73 · 307 and 2n − 1 = 11205 = 3³ · 5 · 83 are not prime
  powers.
- Miyamoto and Yamada–Kiyasu do not apply: 5603 is not a prime power, and neither is n − 2 = 5601 = 3 · 1867.
- Cooper–Wallis does not apply: for each factorization n = a · b it needs T-matrices of order a and Williamson-type
  matrices of order b, and no such pair is known to us.
- SageMath's development version (October 17, 2026) does not construct it: `hadamard_matrix(22412, existence=True)`
  returns Unknown.

We found no construction or statement of order 22412. The sources searched are listed in
[NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, τ, power, c, the T-matrices, the four sequences and the SHA-256 of the matrix |
| `sequences.txt` | the four sequences of length 5603 |

The matrix itself is not stored; `verify/verify.py` rebuilds it.

## Verification

```
python3 verify/verify.py 22412
```

SHA-256 of `H22412.csv`: `7e2ebf56260844fd1fafcd23cb42ac60f0aef9ce878f0719a457cbb4f6ee8acf`
