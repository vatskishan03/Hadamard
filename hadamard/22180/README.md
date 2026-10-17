# Hadamard matrix of order 22180

22180 = 4 · 5545, with 5545 = 5 · 1109.

## Construction

Miyamoto's Hadamard matrix of order 4436 (q = 1109, a prime ≡ 1 mod 4) is built from Paley's conference matrix of order
1110 and a Hadamard matrix K of order 1108. Here K is the Williamson array of four symmetric, pairwise commuting ±1
matrices of order w = 277 (symmetric Williamson-type matrices) from Seberry and Yamada's constructions. Its four
quadrants then lie in a commutative algebra, so Miyamoto's matrix is a 4 × 4 array of bordered blocks of order 1109 that
commute and are transposed by one permutation. Multiplying it by T-matrices of order 5 in the Goethals–Seidel array
gives the Hadamard matrix of order 4 · 1109 · 5 = 22180. The construction and its proof are in
[MIYAMOTO_BLOCKS.md](../MIYAMOTO_BLOCKS.md#symmetric-williamson-type-matrices).

- Paley's conference matrix of order 1110, its points listed as ∞, 0, g^(2j), g^(2j+1) (j < 554) with the primitive root
  g = 2.
- Symmetric Williamson-type matrices of order 277 by Corollary 29(ii) of Seberry and Yamada (J. Combin. Math. Combin.
  Comput. 7 (1990)), with Paley's symmetric conference matrix of order 138 (p = 137); Seberry and Yamada list 277 among
  the orders this corollary gives (their Remark 30). K is their Williamson array [[X₁, X₂, X₃, X₄], [−X₂, X₁, −X₄, X₃],
  [−X₃, X₄, X₁, −X₂], [−X₄, −X₃, X₂, X₁]], relabelled so that the half swap of each quadrant becomes j ↦ 1 − j on ℤ_554.
  `search/symwt.py` builds them (recipe `cor29c(277)`) and checks that they are symmetric, commute and have X₁² + X₂² +
  X₃² + X₄² = 1108I.
- T-matrices of order 5 (first rows, from base sequences BS(3, 2), built by the program): `++000`, `00+00`, `000+-`,
  `00000`.

## Search

None: the matrix is determined by q, the recipe of the Williamson-type matrices and the T-matrices. `python3
search/miyamoto_blocks.py certify 1109 5 "wt:cor29c(277)"` writes the certificate.

## Status

**First explicit construction** ([MIYAMOTO_BLOCKS.md](../MIYAMOTO_BLOCKS.md#symmetric-williamson-type-matrices)). Table
1.53 of the Handbook of Combinatorial Designs (2nd ed., 2007) gives n = 5545 with t = 3, and Đoković (arXiv:0912.5091)
lists it in his set Δ of entries that should have had t = 2. His proof for Δ takes Williamson-type matrices from
Seberry–Yamada's Table A.1 (1992); for 5545 they are of order 1109, which Table A.1 obtains only from Corollary 8.8(i)
(key wl, with Williamson-type matrices of order 277 = (1109 − 1)/4). The printed proof of that corollary does not give
Williamson-type matrices ([NOTES.md](../../NOTES.md) §7). We found no other construction of Williamson-type matrices of
order 1109. The route here needs only symmetric ones of order (q − 1)/4 = 277, which Seberry and Yamada's corollaries
give.

- Seberry and Yamada (J. Combin. Math. Combin. Comput. 7 (1990), Theorem 34 and Corollary 35) write Miyamoto's matrix as
  sixteen bordered blocks, but their multiplication theorems need circulant or type-one blocks (Theorem 4) or, for
  general blocks, sixteen matrices that exist for no odd order t > 1 (Theorem 21). With K the Williamson array of
  symmetric Williamson-type matrices, relabelled, the blocks commute and are transposed by one permutation, so
  T-matrices of odd order apply; we found this nowhere. Every ingredient was published before 2026.
- Paley's constructions do not apply: 4n − 1 = 22179 = 3 · 7393 and 2n − 1 = 11089 = 13 · 853 are not prime powers.
- Miyamoto's theorem does not apply: 5545 is not a prime power, and a symmetric conference matrix of order 5546 is not
  known to us.
- Yamada–Kiyasu does not apply: n − 2 = 5543 = 23 · 241 is not a prime power.
- Cooper–Wallis needs Williamson-type matrices of order 1109 (only through Corollary 8.8(i), above) or T-matrices of
  order 1109; neither is known to us.
- The Paley-block route does not reach it: 4s − 1 = 4435 and 22179 are not prime powers ≡ 3 (mod 8) for the divisors s
  of 5545, and the remaining ones need T-matrices of order 5545 (s = 1) or 1109 (s = 5), which are not known to us.
- The Miyamoto-block route with a negaperiodic or periodic Golay pair of length 554 does not apply: none of that length
  is known to us.
- SageMath 10.8 (passagemath 10.8.13, October 2026) does not construct it: `hadamard_matrix(22180, existence=True)`
  raises UnknownError in its Miyamoto test (`hadamard_matrix_miyamoto_construction`), since the existence of a symmetric
  conference matrix of order 5546 is Unknown to it.

We found no explicit Hadamard matrix of order 22180 before this one. The sources searched are listed in
[NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, g, the four symmetric Williamson-type matrices (upper triangles, base64), the T-matrices and the SHA-256 of the matrix |

The matrix itself is not stored; `verify/verify.py` rebuilds it, after checking that the four matrices are symmetric ±1
matrices that commute pairwise and have X₁² + X₂² + X₃² + X₄² = 4wI.

## Verification

```
python3 verify/verify.py 22180
```

SHA-256 of `H22180.csv`: `b01773d5901de8a53f903200a391670a1525a7ed5f9f29df53ee4e0442c17994`
