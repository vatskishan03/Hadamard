# Hadamard matrix of order 28092

28092 = 4 · 7023, with 7023 = 3 · 2341.

## Construction

Miyamoto's Hadamard matrix of order 9364 (q = 2341, a prime ≡ 1 mod 4), built from Paley's conference matrix of order
2342 and a Hadamard matrix of order 2340 of 2N-type, is a 4 × 4 array of bordered blocks of order 2341 that commute and
are reversed by one permutation. Multiplying it by T-matrices of order 3 in the Goethals–Seidel array gives the Hadamard
matrix of order 4 · 2341 · 3 = 28092. The construction and its proof are in [MIYAMOTO_BLOCKS.md](../MIYAMOTO_BLOCKS.md).

- Paley's conference matrix of order 2342, its points listed as ∞, 0, g^(2j), g^(2j+1) (j < 1170) with the primitive
  root g = 7.
- The Hadamard matrix of order 2340 is [[A, B], [−Bᵀ, Aᵀ]] with A, B negacyclic of order 1170: a negaperiodic Golay pair
  of length 1170, from the negacyclic Paley conference matrix of order 2340 (points of PG(1, 2339) along a primitive
  element of GF(2339²)), 2339 ≡ 3 (mod 4), whose first row (0, c₁, …, c_2339) gives the pair (1, c₂, c₄, …, c_2338),
  (c₁, c₃, …, c_2339) (Balonin–Đoković, Proposition 8: the second Paley series).
- T-matrices of order 3 (from base sequences BS(2, 1), built by the program): first rows `+00`, `0+0`, `00+`, `000`.

## Search

None: the matrix is determined by q, the negaperiodic Golay pair and the T-matrices. `python3 search/miyamoto_blocks.py
certify 2341 3 paley2:2339` writes the certificate.

## Status

**New, beyond the tables, but obtainable classically** ([MIYAMOTO_BLOCKS.md](../MIYAMOTO_BLOCKS.md#status)). In Table
1.53 of the Handbook of Combinatorial Designs (2nd ed., 2007), n = 7023 has t = 3: the smallest known Hadamard matrix
with odd part 7023 had order 2³ · 7023. The order lies above the range of Cati–Pasechnik Table 4, which stops at n =
2999.

- The entry does not appear in Đoković's updates (arXiv:0912.5091, arXiv:1008.2043).
- Seberry and Yamada (J. Combin. Math. Combin. Comput. 7 (1990), Theorem 34 and Corollary 35) write Miyamoto's matrix as
  sixteen bordered blocks, but their multiplication theorems need circulant or type-one blocks (Theorem 4) or, for
  general blocks, sixteen matrices that exist for no odd order t > 1 (Theorem 21). With a Hadamard matrix of 2N-type and
  Paley's points listed along a primitive root, the blocks commute and are transposed by one reversal, so T-matrices of
  odd order apply; we found this nowhere. Every ingredient is published (Miyamoto, J. Combin. Theory Ser. A 57 (1991);
  the negaperiodic Golay pairs of Balonin and Đoković; T-matrices).
- Paley's constructions do not apply: 4n − 1 = 28091 = 7 · 4013 and 2n − 1 = 14045 = 5 · 53² are not prime powers.
- Miyamoto's theorem does not apply: it needs a symmetric conference matrix of order n + 1, which cannot exist for n ≡ 3
  (mod 4); Yamada–Kiyasu does not apply: n − 2 = 7021 = 7 · 17 · 59 is not a prime power.
- Cooper–Wallis does not apply: for each factorization n = a · b it needs T-matrices of order a and Williamson-type
  matrices of order b, and no such pair is known to us.
- The Paley-block route ([PALEY_BLOCKS.md](../PALEY_BLOCKS.md)) does not reach it.
- SageMath's development version (17 October 2026) does not construct it: `hadamard_matrix(28092, existence=True)`
  returns Unknown.

We found no construction or statement of order 28092. The sources searched are listed in [NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, g, the negaperiodic Golay pair, the T-matrices and the SHA-256 of the matrix |

The matrix itself is not stored; `verify/verify.py` rebuilds it.

## Verification

```
python3 verify/verify.py 28092
```

SHA-256 of `H28092.csv`: `6869e348c07e8f3ae673a4c60dcaf3e0a78eeacfd16f486e3fb1c95fb346ff93`
