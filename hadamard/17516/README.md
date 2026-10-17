# Hadamard matrix of order 17516

17516 = 4 · 4379, with 4379 = 29 · 151.

## Construction

Miyamoto's Hadamard matrix of order 116 (q = 29, a prime ≡ 1 mod 4), built from Paley's conference matrix of order 30
and a Hadamard matrix of order 28 of 2N-type, is a 4 × 4 array of bordered blocks of order 29 that commute and are
reversed by one permutation. Multiplying it by T-matrices of order 151 in the Goethals–Seidel array gives the Hadamard
matrix of order 4 · 29 · 151 = 17516. The construction and its proof are in [MIYAMOTO_BLOCKS.md](../MIYAMOTO_BLOCKS.md).

- Paley's conference matrix of order 30, its points listed as ∞, 0, g^(2j), g^(2j+1) (j < 14) with the primitive root g
  = 2.
- The Hadamard matrix of order 28 is [[A, B], [−Bᵀ, Aᵀ]] with A, B negacyclic of order 14: a negaperiodic Golay pair of
  length 14, from the negacyclic Paley conference matrix of order 14 (points of PG(1, 13) along a primitive element of
  GF(13²)), whose first row (0, c₁, …, c_13) gives the pair (1, c₁, …, c_13), (−1, c₁, …, c_13) (Balonin–Đoković,
  arXiv:1508.00640, Proposition 6).
- T-matrices of order 151: [`t-matrices/T151.txt`](../../t-matrices/T151.txt) (2026).

## Search

None: the matrix is determined by q, the negaperiodic Golay pair and the T-matrices. `python3 search/miyamoto_blocks.py
certify 29 151 paley1:13 t-matrices/T151.txt` writes the certificate.

## Status

**New, beyond the tables.** In Table 1.53 of the Handbook of Combinatorial Designs (2nd ed., 2007), n = 4379 has t = 3:
the smallest known Hadamard matrix with odd part 4379 had order 2³ · 4379. The order lies above the range of
Cati–Pasechnik Table 4, which stops at n = 2999.

- The entry does not appear in Đoković's updates (arXiv:0912.5091, arXiv:1008.2043).
- It follows from published results combined in a way we found nowhere: Miyamoto's theorem (J. Combin. Theory Ser. A 57
  (1991)), the negaperiodic Golay pairs of Ito and of Balonin and Đoković, and the Goethals–Seidel array with
  T-matrices, applied to Miyamoto's blocks.
- Paley's constructions do not apply: 4n − 1 = 17515 = 5 · 31 · 113 and 2n − 1 = 8757 = 3² · 7 · 139 are not prime
  powers.
- Miyamoto and Yamada–Kiyasu do not apply: 4379 is not a prime power, and n − 2 = 4377 = 3 · 1459.
- Cooper–Wallis also gives it, from T-matrices of order 151 posted in 2026 and Williamson-type matrices of order 29
  (Seberry–Yamada 1992, Table A.1).
- The Paley-block route ([PALEY_BLOCKS.md](../PALEY_BLOCKS.md)) does not reach it.
- SageMath's development version (October 17, 2026) does not construct it: `hadamard_matrix(17516, existence=True)`
  returns Unknown.

We found no construction or statement of order 17516. The sources searched are listed in [NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, g, the negaperiodic Golay pair, the T-matrices and the SHA-256 of the matrix |

The matrix itself is not stored; `verify/verify.py` rebuilds it.

## Verification

```
python3 verify/verify.py 17516
```

SHA-256 of `H17516.csv`: `55c2d06c57976793206bf9ea9ed3a43c5960d67b0d8390870f8162b75ed7afc4`
