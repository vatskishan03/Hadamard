# Hadamard matrix of order 30548

30548 = 4 · 7637, with 7637 = 7 · 1091.

## Construction

Paley's skew-Hadamard matrix of order 4364 (q = 4363, a prime ≡ 3 mod 8) is a 4 × 4 array of circulant blocks of
order s = 1091 once its points are listed along the four orbits of a torus element of PSL(2, 4363) and rescaled by
signs. Multiplying it by T-matrices of order 7 (Seberry–Yamada 1990, Theorem 4) gives four circulant ±1 sequences of
length 7637 with Σ PAF = 0 off 0, and the Goethals–Seidel array gives the Hadamard matrix. The construction and its
proof are in [PALEY_BLOCKS.md](../PALEY_BLOCKS.md).

- Torus element: the square of x ↦ −1/(x + τ) with τ = 3; sign normalisation c = +1.
- T-matrices of order 7 (from base sequences BS(4, 3), built by the program): first rows `+000000`, `0++-000`,
  `0000+0+`, `00000+0`.
- Row sums of the four sequences: −117, −3, 103, −79 ((−117)² + (−3)² + 103² + (−79)² = 30548).

## Search

None: the matrix is determined by q and the T-matrices. `python3 search/paley_blocks.py certify 4363 7` writes the
certificate.

## Status

**New, beyond the tables, but obtainable classically** ([PALEY_BLOCKS.md](../PALEY_BLOCKS.md#status)). In Table 1.53
of the Handbook of Combinatorial Designs (2nd ed., 2007), n = 7637 has t = 3: the smallest known Hadamard matrix with
odd part 7637 had order 2³ · 7637. The order lies above the range of Cati–Pasechnik Table 4, which stops at n = 2999.

- The entry does not appear in Đoković's updates (arXiv:0912.5091, arXiv:1008.2043).
- It follows from two published theorems that we found combined nowhere: Paley's matrix of order 4364 is equivalent
  to an array of sixteen circulant blocks of order 1091 (Schmidt 1999, after Ito; Yamada 1991), and Seberry and
  Yamada (1990, Theorem 4) multiply such an array by T-matrices of order 7.
- Paley's constructions do not apply: 4n − 1 = 30547 = 11 · 2777 and 2n − 1 = 15273 = 3² · 1697 are not prime powers.
- Miyamoto's theorem does not apply: a symmetric conference matrix of order n + 1 = 7638 cannot exist (7637 is not a sum
  of two squares); Yamada–Kiyasu does not apply: n − 2 = 7635 = 3 · 5 · 509 is not a prime power.
- Cooper–Wallis does not apply: for each factorization n = a · b it needs T-matrices of order a and Williamson-type
  matrices of order b, and no such pair is known to us.
- SageMath's development version (October 17, 2026) does not construct it: `hadamard_matrix(30548, existence=True)`
  returns Unknown.

We found no construction or statement of order 30548. The sources searched are listed in
[NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, τ, power, c, the T-matrices, the four sequences and the SHA-256 of the matrix |
| `sequences.txt` | the four sequences of length 7637 |

The matrix itself is not stored; `verify/verify.py` rebuilds it.

## Verification

```
python3 verify/verify.py 30548
```

SHA-256 of `H30548.csv`: `40f10ffeaf0c76f5b48c4e9e50cf406256ab619b533f601265fc47f7aae0b95d`
