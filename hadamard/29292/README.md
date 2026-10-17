# Hadamard matrix of order 29292

29292 = 4 · 7323, with 7323 = 3 · 2441.

## Construction

Miyamoto's Hadamard matrix of order 9764 (q = 2441, a prime ≡ 1 mod 4), built from Paley's conference matrix of order
2442 and a Hadamard matrix of order 2440 of 2N-type, is a 4 × 4 array of bordered blocks of order 2441 that commute and
are reversed by one permutation. Multiplying it by T-matrices of order 3 in the Goethals–Seidel array gives the Hadamard
matrix of order 4 · 2441 · 3 = 29292. The construction and its proof are in [MIYAMOTO_BLOCKS.md](../MIYAMOTO_BLOCKS.md).

- Paley's conference matrix of order 2442, its points listed as ∞, 0, g^(2j), g^(2j+1) (j < 1220) with the primitive
  root g = 6.
- The Hadamard matrix of order 2440 is [[A, B], [−Bᵀ, Aᵀ]] with A, B negacyclic of order 1220: a negaperiodic Golay pair
  of length 1220, from Turyn's multiplication (Balonin–Đoković, Proposition 3, formulas (15)–(16)) of Golay's pair of
  length 10 with an NG pair of length 122 from the negacyclic Paley conference matrix of order 122 (points of PG(1, 121)
  along a primitive element of GF(121²)), whose first row (0, c₁, …, c_121) gives the pair (1, c₁, …, c_121), (−1, c₁,
  …, c_121) (Balonin–Đoković, arXiv:1508.00640, Proposition 6).
- T-matrices of order 3 (from base sequences BS(2, 1), built by the program): first rows `+00`, `0+0`, `00+`, `000`.

## Search

None: the matrix is determined by q, the negaperiodic Golay pair and the T-matrices. `python3 search/miyamoto_blocks.py
certify 2441 3 golay10:paley1:121` writes the certificate.

## Status

**First explicit construction** ([MIYAMOTO_BLOCKS.md](../MIYAMOTO_BLOCKS.md#status)). In Table 1.53 of the Handbook of
Combinatorial Designs (2nd ed., 2007), n = 7323 has t = 3, and Đoković's update arXiv:1008.2043 lists it as settled,
through Seberry–Yamada's Corollary 8.8(i) only (below). The order lies above the range of Cati–Pasechnik Table 4, which
stops at n = 2999. Its ingredients were all known before 2026.

- Đoković's update arXiv:1008.2043 lists it (Proposition 3.2) only through Seberry–Yamada's Corollary 8.8(i), whose
  printed proof does not give Williamson-type matrices (see below).
- It follows from published results combined in a way we found nowhere: Miyamoto's theorem (J. Combin. Theory Ser. A 57
  (1991)), the negaperiodic Golay pairs of Ito and of Balonin and Đoković, and the Goethals–Seidel array with
  T-matrices, applied to Miyamoto's blocks.
- Paley's constructions do not apply: 4n − 1 = 29291 = 17 · 1723 and 2n − 1 = 14645 = 5 · 29 · 101 are not prime powers.
- Miyamoto and Yamada–Kiyasu do not apply: 7323 is not a prime power, and n − 2 = 7321 is a prime ≡ 1 (mod 8).
- Cooper–Wallis does not apply: for each factorization n = a · b it needs T-matrices of order a and Williamson-type
  matrices of order b, and no such pair is known to us.
- The Paley-block route ([PALEY_BLOCKS.md](../PALEY_BLOCKS.md)) does not reach it.
- SageMath's development version (October 17, 2026) does not construct it: `hadamard_matrix(29292, existence=True)`
  returns Unknown.

Đoković (arXiv:1008.2043, Proposition 3.2) obtained 7323 = 3 · 2441 from Williamson-type matrices of order 2441, which
he took from Corollary 8.8(i) of Seberry and Yamada (1992; Corollary 29(i) of JCMCC 7, 1990) with a Hadamard matrix of
order 1220. The printed proof of that corollary gives pairwise amicable matrices only for a symmetric input
([NOTES.md](../../NOTES.md) §7), and no symmetric Hadamard matrix of order 1220 is known to us; Miyamoto's 1991 paper,
credited with the Williamson-type case, was not accessible to us. The construction here uses Miyamoto's main theorem,
not that corollary. We found no valid construction of order 29292 before this one. The sources searched are listed in
[NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, g, the negaperiodic Golay pair, the T-matrices and the SHA-256 of the matrix |

The matrix itself is not stored; `verify/verify.py` rebuilds it.

## Verification

```
python3 verify/verify.py 29292
```

SHA-256 of `H29292.csv`: `294be44ca446e2e1d40767b907e51af0fb85ea8d9dc5226c2b5764c8bc67579b`
