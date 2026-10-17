# Hadamard matrix of order 30604

30604 = 4 · 7651, with 7651 = 7 · 1093.

## Construction

Miyamoto's Hadamard matrix of order 4372 (q = 1093, a prime ≡ 1 mod 4), built from Paley's conference matrix of order
1094 and a Hadamard matrix of order 1092 of 2N-type, is a 4 × 4 array of bordered blocks of order 1093 that commute and
are reversed by one permutation. Multiplying it by T-matrices of order 7 in the Goethals–Seidel array gives the Hadamard
matrix of order 4 · 1093 · 7 = 30604. The construction and its proof are in [MIYAMOTO_BLOCKS.md](../MIYAMOTO_BLOCKS.md).

- Paley's conference matrix of order 1094, its points listed as ∞, 0, g^(2j), g^(2j+1) (j < 546) with the primitive root
  g = 5.
- The Hadamard matrix of order 1092 is [[A, B], [−Bᵀ, Aᵀ]] with A, B negacyclic of order 546: a negaperiodic Golay pair
  of length 546, from the negacyclic Paley conference matrix of order 1092 (points of PG(1, 1091) along a primitive
  element of GF(1091²)), 1091 ≡ 3 (mod 4), whose first row (0, c₁, …, c_1091) gives the pair (1, c₂, c₄, …, c_1090),
  (c₁, c₃, …, c_1091) (Balonin–Đoković, Proposition 8; Ito's series).
- T-matrices of order 7 (from base sequences BS(4, 3), built by the program): first rows `+000000`, `0++-000`,
  `0000+0+`, `00000+0`.

## Search

None: the matrix is determined by q, the negaperiodic Golay pair and the T-matrices. `python3 search/miyamoto_blocks.py
certify 1093 7 paley2:1091` writes the certificate.

## Status

**New, beyond the tables, but obtainable classically** ([MIYAMOTO_BLOCKS.md](../MIYAMOTO_BLOCKS.md#status)). In Table
1.53 of the Handbook of Combinatorial Designs (2nd ed., 2007), n = 7651 has t = 3: the smallest known Hadamard matrix
with odd part 7651 had order 2³ · 7651. The order lies above the range of Cati–Pasechnik Table 4, which stops at n =
2999.

- The entry does not appear in Đoković's updates (arXiv:0912.5091, arXiv:1008.2043).
- It follows from published results combined in a way we found nowhere: Miyamoto's theorem (J. Combin. Theory Ser. A 57
  (1991)), the negaperiodic Golay pairs of Ito and of Balonin and Đoković, and the Goethals–Seidel array with
  T-matrices, applied to Miyamoto's blocks.
- Paley's constructions do not apply: 4n − 1 = 30603 = 3 · 101² and 2n − 1 = 15301 = 11 · 13 · 107 are not prime powers.
- Miyamoto and Yamada–Kiyasu do not apply: 7651 is not a prime power, and n − 2 = 7649 is a prime ≡ 1 (mod 8).
- Cooper–Wallis does not apply: for each factorization n = a · b it needs T-matrices of order a and Williamson-type
  matrices of order b, and no such pair is known to us.
- The Paley-block route ([PALEY_BLOCKS.md](../PALEY_BLOCKS.md)) does not reach it.
- SageMath's development version (October 17, 2026) does not construct it: `hadamard_matrix(30604, existence=True)`
  returns Unknown.

We found no construction or statement of order 30604. The sources searched are listed in [NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, g, the negaperiodic Golay pair, the T-matrices and the SHA-256 of the matrix |

The matrix itself is not stored; `verify/verify.py` rebuilds it.

## Verification

```
python3 verify/verify.py 30604
```

SHA-256 of `H30604.csv`: `67892ddde17e38cba0b5dee9e99d662c0515490dcc45fea1fcd47abfa2418d0b`
