# Hadamard matrix of order 20028

20028 = 4 · 5007, with 5007 = 3 · 1669.

## Construction

Miyamoto's Hadamard matrix of order 6676 (q = 1669, a prime ≡ 1 mod 4), built from Paley's conference matrix of order
1670 and a Hadamard matrix of order 1668 of 2N-type, is a 4 × 4 array of bordered blocks of order 1669 that commute and
are reversed by one permutation. Multiplying it by T-matrices of order 3 in the Goethals–Seidel array gives the Hadamard
matrix of order 4 · 1669 · 3 = 20028. The construction and its proof are in [MIYAMOTO_BLOCKS.md](../MIYAMOTO_BLOCKS.md).

- Paley's conference matrix of order 1670, its points listed as ∞, 0, g^(2j), g^(2j+1) (j < 834) with the primitive root
  g = 2.
- The Hadamard matrix of order 1668 is [[A, B], [−Bᵀ, Aᵀ]] with A, B negacyclic of order 834: a negaperiodic Golay pair
  of length 834, from the negacyclic Paley conference matrix of order 1668 (points of PG(1, 1667) along a primitive
  element of GF(1667²)), 1667 ≡ 3 (mod 4), whose first row (0, c₁, …, c_1667) gives the pair (1, c₂, c₄, …, c_1666),
  (c₁, c₃, …, c_1667) (Balonin–Đoković, Proposition 8; Ito's series).
- T-matrices of order 3 (from base sequences BS(2, 1), built by the program): first rows `+00`, `0+0`, `00+`, `000`.

## Search

None: the matrix is determined by q, the negaperiodic Golay pair and the T-matrices. `python3 search/miyamoto_blocks.py
certify 1669 3 paley2:1667` writes the certificate.

## Status

**New, beyond the tables, but obtainable classically** ([MIYAMOTO_BLOCKS.md](../MIYAMOTO_BLOCKS.md#status)). In Table
1.53 of the Handbook of Combinatorial Designs (2nd ed., 2007), n = 5007 has t = 3: the smallest known Hadamard matrix
with odd part 5007 had order 2³ · 5007. The order lies above the range of Cati–Pasechnik Table 4, which stops at n =
2999.

- The entry does not appear in Đoković's updates (arXiv:0912.5091, arXiv:1008.2043).
- It follows from published results combined in a way we found nowhere: Miyamoto's theorem (J. Combin. Theory Ser. A 57
  (1991)), the negaperiodic Golay pairs of Ito and of Balonin and Đoković, and the Goethals–Seidel array with
  T-matrices, applied to Miyamoto's blocks.
- Paley's constructions do not apply: 4n − 1 = 20027 = 7 · 2861 and 2n − 1 = 10013 = 17 · 19 · 31 are not prime powers.
- Miyamoto and Yamada–Kiyasu do not apply: 5007 is not a prime power, and n − 2 = 5005 = 5 · 7 · 11 · 13.
- Cooper–Wallis does not apply: for each factorization n = a · b it needs T-matrices of order a and Williamson-type
  matrices of order b, and no such pair is known to us.
- The Paley-block route ([PALEY_BLOCKS.md](../PALEY_BLOCKS.md)) does not reach it.
- SageMath's development version (October 17, 2026) does not construct it: `hadamard_matrix(20028, existence=True)`
  returns Unknown.

We found no construction or statement of order 20028. The sources searched are listed in [NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, g, the negaperiodic Golay pair, the T-matrices and the SHA-256 of the matrix |

The matrix itself is not stored; `verify/verify.py` rebuilds it.

## Verification

```
python3 verify/verify.py 20028
```

SHA-256 of `H20028.csv`: `4a75fc0d003afb5a6c4c8d62565a60be69a10ea2d4bce13990ad5a8381130398`
