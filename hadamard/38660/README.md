# Hadamard matrix of order 38660

38660 = 4 · 9665, with 9665 = 5 · 1933.

## Construction

Miyamoto's Hadamard matrix of order 7732 (q = 1933, a prime ≡ 1 mod 4), built from Paley's conference matrix of order
1934 and a Hadamard matrix of order 1932 of 2N-type, is a 4 × 4 array of bordered blocks of order 1933 that commute and
are reversed by one permutation. Multiplying it by T-matrices of order 5 in the Goethals–Seidel array gives the Hadamard
matrix of order 4 · 1933 · 5 = 38660. The construction and its proof are in [MIYAMOTO_BLOCKS.md](../MIYAMOTO_BLOCKS.md).

- Paley's conference matrix of order 1934, its points listed as ∞, 0, g^(2j), g^(2j+1) (j < 966) with the primitive root
  g = 5.
- The Hadamard matrix of order 1932 is [[A, B], [−Bᵀ, Aᵀ]] with A, B negacyclic of order 966: a negaperiodic Golay pair
  of length 966, from the negacyclic Paley conference matrix of order 1932 (points of PG(1, 1931) along a primitive
  element of GF(1931²)), 1931 ≡ 3 (mod 4), whose first row (0, c₁, …, c_1931) gives the pair (1, c₂, c₄, …, c_1930),
  (c₁, c₃, …, c_1931) (Balonin–Đoković, Proposition 8: the second Paley series).
- T-matrices of order 5 (from base sequences BS(3, 2), built by the program): first rows `++000`, `00+00`, `000+-`,
  `00000`.

## Search

None: the matrix is determined by q, the negaperiodic Golay pair and the T-matrices. `python3 search/miyamoto_blocks.py
certify 1933 5 paley2:1931` writes the certificate.

## Status

**New, beyond the tables, but obtainable classically** ([MIYAMOTO_BLOCKS.md](../MIYAMOTO_BLOCKS.md#status)). In Table
1.53 of the Handbook of Combinatorial Designs (2nd ed., 2007), n = 9665 has t = 3: the smallest known Hadamard matrix
with odd part 9665 had order 2³ · 9665. The order lies above the range of Cati–Pasechnik Table 4, which stops at n =
2999.

- The entry does not appear in Đoković's updates (arXiv:0912.5091, arXiv:1008.2043).
- Seberry and Yamada (J. Combin. Math. Combin. Comput. 7 (1990), Theorem 34 and Corollary 35) write Miyamoto's matrix as
  sixteen bordered blocks, but their multiplication theorems need circulant or type-one blocks (Theorem 4) or, for
  general blocks, sixteen matrices that exist for no odd order t > 1 (Theorem 21). With a Hadamard matrix of 2N-type and
  Paley's points listed along a primitive root, the blocks commute and are transposed by one reversal, so T-matrices of
  odd order apply; we found this nowhere. Every ingredient is published (Miyamoto, J. Combin. Theory Ser. A 57 (1991);
  the negaperiodic Golay pairs of Balonin and Đoković; T-matrices).
- Paley's constructions do not apply: 4n − 1 = 38659 = 67 · 577 and 2n − 1 = 19329 = 3 · 17 · 379 are not prime powers.
- Miyamoto's theorem does not apply: it needs a symmetric conference matrix of order n + 1 = 9666, and none is known to
  us (9665 is not a prime power); Yamada–Kiyasu does not apply: n − 2 = 9663 = 3 · 3221 is not a prime power.
- Cooper–Wallis does not apply: for each factorization n = a · b it needs T-matrices of order a and Williamson-type
  matrices of order b, and no such pair is known to us.
- The Paley-block route ([PALEY_BLOCKS.md](../PALEY_BLOCKS.md)) does not reach it.
- SageMath's development version (October 17, 2026) does not construct it: `hadamard_matrix(38660, existence=True)`
  raises UnknownError in its Miyamoto test (no symmetric conference matrix of order 9666 is known to it), and returns
  Unknown when that Unknown is treated as False.

We found no construction or statement of order 38660. The sources searched are listed in [NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, g, the negaperiodic Golay pair, the T-matrices and the SHA-256 of the matrix |

The matrix itself is not stored; `verify/verify.py` rebuilds it.

## Verification

```
python3 verify/verify.py 38660
```

SHA-256 of `H38660.csv`: `ce3bffa461ecc72c63a05f172498d28911c3ff0d66ac43feb3626539c5fab9d7`
