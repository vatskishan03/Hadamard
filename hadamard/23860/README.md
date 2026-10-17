# Hadamard matrix of order 23860

23860 = 4 · 5965, with 5965 = 5 · 1193.

## Construction

Miyamoto's Hadamard matrix of order 4772 (q = 1193, a prime ≡ 1 mod 4) is built from Paley's conference matrix of order
1194 and a Hadamard matrix K of order 1192. Here K is the Williamson array of four symmetric, pairwise commuting ±1
matrices of order w = 298 (symmetric Williamson-type matrices) from Seberry and Yamada's constructions. Its four
quadrants then lie in a commutative algebra, so Miyamoto's matrix is a 4 × 4 array of bordered blocks of order 1193 that
commute and are transposed by one permutation. Multiplying it by T-matrices of order 5 in the Goethals–Seidel array
gives the Hadamard matrix of order 4 · 1193 · 5 = 23860. The construction and its proof are in
[MIYAMOTO_BLOCKS.md](../MIYAMOTO_BLOCKS.md#symmetric-williamson-type-matrices).

- Paley's conference matrix of order 1194, its points listed as ∞, 0, g^(2j), g^(2j+1) (j < 596) with the primitive root
  g = 3.
- Symmetric Williamson-type matrices of order 298 = 2 · 1 · 149 by Corollary 25 of Seberry and Yamada (J. Combin. Math.
  Combin. Comput. 7 (1990)) (order 2uv from orders u and v) from the trivial ones of order 1 and symmetric
  Williamson-type matrices of order 149, which their Corollary 29(ii) gives from Paley's symmetric conference matrix of
  order 74 (p = 73). K is their Williamson array [[X₁, X₂, X₃, X₄], [−X₂, X₁, −X₄, X₃], [−X₃, X₄, X₁, −X₂], [−X₄, −X₃,
  X₂, X₁]], relabelled so that the half swap of each quadrant becomes j ↦ 1 − j on ℤ_596. `search/symwt.py` builds them
  (recipe `cor25(one(), cor29c(149))`) and checks that they are symmetric, commute and have X₁² + X₂² + X₃² + X₄² =
  1192I.
- T-matrices of order 5 (first rows, from base sequences BS(3, 2), built by the program): `++000`, `00+00`, `000+-`,
  `00000`.

## Search

None: the matrix is determined by q, the recipe of the Williamson-type matrices and the T-matrices. `python3
search/miyamoto_blocks.py certify 1193 5 "wt:cor25(one(), cor29c(149))"` writes the certificate.

## Status

**First explicit construction** ([MIYAMOTO_BLOCKS.md](../MIYAMOTO_BLOCKS.md#symmetric-williamson-type-matrices)). Table
1.53 of the Handbook of Combinatorial Designs (2nd ed., 2007) gives n = 5965 with t = 3, and Đoković (arXiv:0912.5091)
lists it in his set Δ of entries that should have had t = 2. His proof for Δ takes Williamson-type matrices from
Seberry–Yamada's Table A.1 (1992); for 5965 they are of order 1193, which Table A.1 obtains only from Corollary 8.8(i)
(key wl, which Seberry and Yamada credit to Miyamoto (1991), with Williamson-type matrices of order 298 = (1193 − 1)/4
from their Corollary 25, key wo; their Remark 30 also lists the variant with an Hadamard matrix of order 596 = (1193 −
1)/2). The printed proof of that corollary does not give Williamson-type matrices ([NOTES.md](../../NOTES.md) §7). We
found no other construction of Williamson-type matrices of order 1193: Lemma 26 would need a symmetric quaternary
Hadamard matrix of order 596 = 4 · 149 (for instance from a repeated-block Williamson-type input), and we know none. The
route here needs only symmetric ones of order (q − 1)/4 = 298, which Seberry and Yamada's corollaries give.

- Seberry and Yamada (J. Combin. Math. Combin. Comput. 7 (1990), Theorem 34 and Corollary 35) write Miyamoto's matrix as
  sixteen bordered blocks, but their multiplication theorems need circulant or type-one blocks (Theorem 4) or, for
  general blocks, sixteen matrices that exist for no odd order t > 1 (Theorem 21). With K the Williamson array of
  symmetric Williamson-type matrices, relabelled, the blocks commute and are transposed by one permutation, so
  T-matrices of odd order apply; we found this nowhere. Every ingredient was published before 2026.
- Paley's constructions do not apply: 4n − 1 = 23859 = 3² · 11 · 241 and 2n − 1 = 11929 = 79 · 151 are not prime powers.
- Miyamoto's theorem does not apply: 5965 is not a prime power, and a symmetric conference matrix of order 5966 is not
  known to us.
- Yamada–Kiyasu does not apply: n − 2 = 5963 = 67 · 89 is not a prime power.
- Cooper–Wallis needs Williamson-type matrices of order 1193 (only through Corollary 8.8(i), above) or T-matrices of
  order 1193; neither is known to us.
- The Paley-block route does not reach it: 4s − 1 = 4771 and 23859 are not prime powers ≡ 3 (mod 8) for the divisors s
  of 5965, and the remaining ones need T-matrices of order 5965 (s = 1) or 1193 (s = 5), which are not known to us.
- The Miyamoto-block route with a negaperiodic or periodic Golay pair of length 596 does not apply: none of that length
  is known to us.
- SageMath 10.8 (passagemath 10.8.13, October 2026) does not construct it: `hadamard_matrix(23860, existence=True)`
  raises UnknownError in its Miyamoto test (`hadamard_matrix_miyamoto_construction`), since the existence of a symmetric
  conference matrix of order 5966 is Unknown to it.

We found no explicit Hadamard matrix of order 23860 before this one. The sources searched are listed in
[NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, g, the four symmetric Williamson-type matrices (upper triangles, base64), the T-matrices and the SHA-256 of the matrix |

The matrix itself is not stored; `verify/verify.py` rebuilds it, after checking that the four matrices are symmetric ±1
matrices that commute pairwise and have X₁² + X₂² + X₃² + X₄² = 4wI.

## Verification

```
python3 verify/verify.py 23860
```

SHA-256 of `H23860.csv`: `e74f987d1089c87e4dceaccdb7fdc6a66599f91bac7ff70ce433a45e3ba66a59`
