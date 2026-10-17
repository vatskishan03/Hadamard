# Hadamard matrix of order 18908

18908 = 4 · 4727, with 4727 = 29 · 163.

## Construction

Miyamoto's Hadamard matrix of order 116 (q = 29, a prime ≡ 1 mod 4), built from Paley's conference matrix of order 30
and a Hadamard matrix of order 28 of 2N-type, is a 4 × 4 array of bordered blocks of order 29 that commute and are
reversed by one permutation. Multiplying it by T-matrices of order 163 in the Goethals–Seidel array gives the Hadamard
matrix of order 4 · 29 · 163 = 18908. The construction and its proof are in [MIYAMOTO_BLOCKS.md](../MIYAMOTO_BLOCKS.md).

- Paley's conference matrix of order 30, its points listed as ∞, 0, g^(2j), g^(2j+1) (j < 14) with the primitive root g
  = 2.
- The Hadamard matrix of order 28 is [[A, B], [−Bᵀ, Aᵀ]] with A, B negacyclic of order 14: a negaperiodic Golay pair of
  length 14, from the negacyclic Paley conference matrix of order 14 (points of PG(1, 13) along a primitive element of
  GF(13²)), whose first row (0, c₁, …, c_13) gives the pair (1, c₁, …, c_13), (−1, c₁, …, c_13) (Balonin–Đoković,
  arXiv:1508.00640, Proposition 6).
- T-matrices of order 163: [`t-matrices/T163.txt`](../../t-matrices/T163.txt) (2026).

## Search

None: the matrix is determined by q, the negaperiodic Golay pair and the T-matrices. `python3 search/miyamoto_blocks.py
certify 29 163 paley1:13 t-matrices/T163.txt` writes the certificate.

## Status

**New, beyond the tables.** In Table 1.53 of the Handbook of Combinatorial Designs (2nd ed., 2007), n = 4727 has t = 3:
the smallest known Hadamard matrix with odd part 4727 had order 2³ · 4727. The order lies above the range of
Cati–Pasechnik Table 4, which stops at n = 2999.

- The entry does not appear in Đoković's updates (arXiv:0912.5091, arXiv:1008.2043).
- Seberry and Yamada (J. Combin. Math. Combin. Comput. 7 (1990), Theorem 34 and Corollary 35) write Miyamoto's matrix as
  sixteen bordered blocks, but their multiplication theorems need circulant or type-one blocks (Theorem 4) or, for
  general blocks, sixteen matrices that exist for no odd order t > 1 (Theorem 21). With a Hadamard matrix of 2N-type and
  Paley's points listed along a primitive root, the blocks commute and are transposed by one reversal, so T-matrices of
  odd order apply; we found this nowhere. Every ingredient is published (Miyamoto, J. Combin. Theory Ser. A 57 (1991);
  the negaperiodic Golay pairs of Balonin and Đoković; T-matrices).
- Paley's constructions do not apply: 4n − 1 = 18907 = 7 · 37 · 73 and 2n − 1 = 9453 = 3 · 23 · 137 are not prime
  powers.
- Miyamoto's theorem does not apply: it needs a symmetric conference matrix of order n + 1, which cannot exist for n ≡ 3
  (mod 4); Yamada–Kiyasu does not apply: n − 2 = 4725 = 3³ · 5² · 7 is not a prime power.
- Cooper–Wallis also gives it, from T-matrices of order 163 posted in 2026 and Williamson-type matrices of order 29
  (Seberry–Yamada 1992, Table A.1).
- The Paley-block route ([PALEY_BLOCKS.md](../PALEY_BLOCKS.md)) does not reach it.
- SageMath's development version (October 17, 2026) does not construct it: `hadamard_matrix(18908, existence=True)`
  returns Unknown.

We found no construction or statement of order 18908. The sources searched are listed in [NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, g, the negaperiodic Golay pair, the T-matrices and the SHA-256 of the matrix |

The matrix itself is not stored; `verify/verify.py` rebuilds it.

## Verification

```
python3 verify/verify.py 18908
```

SHA-256 of `H18908.csv`: `1f9c4213e995a4ec9ada00e018e009c108d2a9abed77809a015283b54bd70492`
