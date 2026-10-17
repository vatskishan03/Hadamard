# Hadamard matrix of order 39860

39860 = 4 · 9965, with 9965 = 5 · 1993.

## Construction

Miyamoto's Hadamard matrix of order 7972 (q = 1993, a prime ≡ 1 mod 4) is built from Paley's conference matrix of order
1994 and a Hadamard matrix K of order 1992. Here K is the Williamson array of four symmetric, pairwise commuting ±1
matrices of order w = 498 (symmetric Williamson-type matrices) from Seberry and Yamada's constructions. Its four
quadrants then lie in a commutative algebra, so Miyamoto's matrix is a 4 × 4 array of bordered blocks of order 1993 that
commute and are transposed by one permutation. Multiplying it by T-matrices of order 5 in the Goethals–Seidel array
gives the Hadamard matrix of order 4 · 1993 · 5 = 39860. The construction and its proof are in
[MIYAMOTO_BLOCKS.md](../MIYAMOTO_BLOCKS.md#symmetric-williamson-type-matrices).

- Paley's conference matrix of order 1994, its points listed as ∞, 0, g^(2j), g^(2j+1) (j < 996) with the primitive root
  g = 5.
- Symmetric Williamson-type matrices of order 498 = 2 · 3 · 83 by Corollary 25 of Seberry and Yamada (J. Combin. Math.
  Combin. Comput. 7 (1990)) from Williamson matrices of order 3 (Turyn's family, q = 5) and symmetric Williamson-type
  matrices of order 83 from their Corollary 31 (q = 41, 2q − 1 = 81: Lemma 26 with Paley's symmetric core of order 41
  and Yamada's symmetric circulants M, N of order 41 with M² + N² = 81I, built over GF(81²)). K is their Williamson
  array [[X₁, X₂, X₃, X₄], [−X₂, X₁, −X₄, X₃], [−X₃, X₄, X₁, −X₂], [−X₄, −X₃, X₂, X₁]], relabelled so that the half swap
  of each quadrant becomes j ↦ 1 − j on ℤ_996. `search/symwt.py` builds them (recipe `cor25(turyn(5), cor31(41))`) and
  checks that they are symmetric, commute and have X₁² + X₂² + X₃² + X₄² = 1992I.
- T-matrices of order 5 (first rows, from base sequences BS(3, 2), built by the program): `++000`, `00+00`, `000+-`,
  `00000`.

## Search

None: the matrix is determined by q, the recipe of the Williamson-type matrices and the T-matrices. `python3
search/miyamoto_blocks.py certify 1993 5 "wt:cor25(turyn(5), cor31(41))"` writes the certificate.

## Status

**First explicit construction** ([MIYAMOTO_BLOCKS.md](../MIYAMOTO_BLOCKS.md#symmetric-williamson-type-matrices)). Table
1.53 of the Handbook of Combinatorial Designs (2nd ed., 2007) gives n = 9965 with t = 3, and Đoković (arXiv:0912.5091)
lists it in his set Δ of entries that should have had t = 2. His proof for Δ takes Williamson-type matrices from
Seberry–Yamada's Table A.1 (1992); for 9965 they are of order 1993, which Table A.1 obtains only from Corollary 8.8(i)
(key wl, which Seberry and Yamada credit to Miyamoto (1991), with Williamson-type matrices of order 498 = (1993 − 1)/4
from their Corollary 25, key wo; their Remark 30 also lists the variant with an Hadamard matrix of order 996 = (1993 −
1)/2). The printed proof of that corollary does not give Williamson-type matrices ([NOTES.md](../../NOTES.md) §7). The
route here needs only symmetric ones of order (q − 1)/4 = 498, which Seberry and Yamada's corollaries give.
Cooper–Wallis also gives this order from ingredients published before 2026 (found by the audit of October 17, 2026):
Lemma 26 needs only a symmetric quaternary Hadamard matrix of order (q − 1)/2 = 996 = 6 · 166, the Kronecker product of
Paley's conference matrix of order 6 and the repeated-block arrangement of those of order 83 from Corollary 31 is one,
and `python3 search/symwt.py "cor29q(1993, qconf(5), qwt(cor31(41)))"` gives symmetric Williamson-type matrices of order
1993 (checked).

- Seberry and Yamada (J. Combin. Math. Combin. Comput. 7 (1990), Theorem 34 and Corollary 35) write Miyamoto's matrix as
  sixteen bordered blocks, but their multiplication theorems need circulant or type-one blocks (Theorem 4) or, for
  general blocks, sixteen matrices that exist for no odd order t > 1 (Theorem 21). With K the Williamson array of
  symmetric Williamson-type matrices, relabelled, the blocks commute and are transposed by one permutation, so
  T-matrices of odd order apply; we found this nowhere. Every ingredient was published before 2026.
- Paley's constructions do not apply: 4n − 1 = 39859 = 23 · 1733 and 2n − 1 = 19929 = 3 · 7 · 13 · 73 are not prime
  powers.
- Miyamoto's theorem does not apply: 9965 is not a prime power, and a symmetric conference matrix of order 9966 is not
  known to us.
- Yamada–Kiyasu does not apply: n − 2 = 9963 = 3⁵ · 41 is not a prime power.
- Cooper–Wallis applies with the Williamson-type matrices of order 1993 of the Kronecker route above and T-matrices of
  order 5; the other split needs T-matrices of order 1993, which are not known to us.
- The Paley-block route does not reach it: 4s − 1 = 7971 and 39859 are not prime powers ≡ 3 (mod 8) for the divisors s
  of 9965, and the remaining ones need T-matrices of order 9965 (s = 1) or 1993 (s = 5), which are not known to us.
- The Miyamoto-block route with a negaperiodic or periodic Golay pair of length 996 does not apply: none of that length
  is known to us; a periodic one cannot exist, since 1992 = 2³ · 3 · 83 is not a sum of two squares.
- SageMath 10.8 (passagemath 10.8.13, October 2026) does not construct it: `hadamard_matrix(39860, existence=True)`
  raises UnknownError in its Miyamoto test (`hadamard_matrix_miyamoto_construction`), since the existence of a symmetric
  conference matrix of order 9966 is Unknown to it.

We found no explicit Hadamard matrix of order 39860 before this one. The sources searched are listed in
[NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, g, the four symmetric Williamson-type matrices (upper triangles, base64), the T-matrices and the SHA-256 of the matrix |

The matrix itself is not stored; `verify/verify.py` rebuilds it, after checking that the four matrices are symmetric ±1
matrices that commute pairwise and have X₁² + X₂² + X₃² + X₄² = 4wI.

## Verification

```
python3 verify/verify.py 39860
```

SHA-256 of `H39860.csv`: `5d01dd1fe4eadce7879099d3771f2753f196f7e12f0a0fa7f6ac8a915736f510`
