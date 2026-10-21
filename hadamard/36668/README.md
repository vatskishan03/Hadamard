# Hadamard matrix of order 36668

36668 = 4 · 9167, with 9167 = 89 · 103.

## Construction

Miyamoto's Hadamard matrix of order 356 (q = 89, a prime ≡ 1 mod 4), built from Paley's conference matrix of order 90
and a Hadamard matrix of order 88 of 2N-type, is a 4 × 4 array of bordered blocks of order 89 that commute and are
reversed by one permutation. Multiplying it by T-matrices of order 103 in the Goethals–Seidel array gives the Hadamard
matrix of order 4 · 89 · 103 = 36668. The construction and its proof are in [MIYAMOTO_BLOCKS.md](../MIYAMOTO_BLOCKS.md).

- Paley's conference matrix of order 90, its points listed as ∞, 0, g^(2j), g^(2j+1) (j < 44) with the primitive root g
  = 3.
- The Hadamard matrix of order 88 is [[A, B], [−Bᵀ, Aᵀ]] with A, B negacyclic of order 44: a negaperiodic Golay pair of
  length 44, from the negacyclic Paley conference matrix of order 44 (points of PG(1, 43) along a primitive element of
  GF(43²)), whose first row (0, c₁, …, c_43) gives the pair (1, c₁, …, c_43), (−1, c₁, …, c_43) (Balonin–Đoković,
  arXiv:1508.00640, Proposition 6).
- T-matrices of order 103: [`t-matrices/T103.txt`](../../t-matrices/T103.txt), decoded from Schneider's matrix of order
  2060 (2026).

## Search

None: the matrix is determined by q, the negaperiodic Golay pair and the T-matrices. `python3 search/miyamoto_blocks.py
certify 89 103 paley1:43 t-matrices/T103.txt` writes the certificate.

## Status

**New, beyond the tables.** In Table 1.53 of the Handbook of Combinatorial Designs (2nd ed., 2007), n = 9167 has t = 3:
the smallest known Hadamard matrix with odd part 9167 had order 2³ · 9167. The order lies above the range of
Cati–Pasechnik Table 4, which stops at n = 2999.

- The entry does not appear in Đoković's updates (arXiv:0912.5091, arXiv:1008.2043).
- Seberry and Yamada (J. Combin. Math. Combin. Comput. 7 (1990), Theorem 34 and Corollary 35) write Miyamoto's matrix as
  sixteen bordered blocks, but their multiplication theorems need circulant or type-one blocks (Theorem 4) or, for
  general blocks, sixteen matrices that exist for no odd order t > 1 (Theorem 21). With a Hadamard matrix of 2N-type and
  Paley's points listed along a primitive root, the blocks commute and are transposed by one reversal, so T-matrices of
  odd order apply; we found this nowhere. Every ingredient is published (Miyamoto, J. Combin. Theory Ser. A 57 (1991);
  the negaperiodic Golay pairs of Balonin and Đoković; T-matrices).
- Paley's constructions do not apply: 4n − 1 = 36667 = 37 · 991 and 2n − 1 = 18333 = 3³ · 7 · 97 are not prime powers.
- Miyamoto's theorem does not apply: it needs a symmetric conference matrix of order n + 1, which cannot exist for n ≡ 3
  (mod 4); Yamada–Kiyasu does not apply: n − 2 = 9165 = 3 · 5 · 13 · 47 is not a prime power.
- Cooper–Wallis also gives it, from T-matrices of order 103 posted in 2026 and Williamson-type matrices of order 89
  (Seberry–Yamada 1992, Table A.1; Table A.1 obtains them through Corollary 8.8(i), but Corollary 8.8(ii) with a
  symmetric Hadamard matrix of order 44 gives them validly).
- The Paley-block route ([PALEY_BLOCKS.md](../PALEY_BLOCKS.md)) does not reach it.
- SageMath's development version (17 October 2026) does not construct it: `hadamard_matrix(36668, existence=True)`
  returns Unknown.

We found no construction or statement of order 36668. The sources searched are listed in [NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, g, the negaperiodic Golay pair, the T-matrices and the SHA-256 of the matrix |

The matrix itself is not stored; `verify/verify.py` rebuilds it.

## Verification

```
python3 verify/verify.py 36668
```

SHA-256 of `H36668.csv`: `d2761fc973fa294baca5a6f0e727e575dda98c418e5a05faf871b4c3b68fc280`
