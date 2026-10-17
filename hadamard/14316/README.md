# Hadamard matrix of order 14316

14316 = 4 · 3579, with 3579 = 3 · 1193.

## Construction

Miyamoto's Hadamard matrix of order 4772 (q = 1193, a prime ≡ 1 mod 4) is built from Paley's conference matrix of order
1194 and a Hadamard matrix K of order 1192. Here K is the Williamson array of four symmetric, pairwise commuting ±1
matrices of order w = 298 (symmetric Williamson-type matrices) from Seberry and Yamada's constructions. Its four
quadrants then lie in a commutative algebra, so Miyamoto's matrix is a 4 × 4 array of bordered blocks of order 1193 that
commute and are transposed by one permutation. Multiplying it by T-matrices of order 3 in the Goethals–Seidel array
gives the Hadamard matrix of order 4 · 1193 · 3 = 14316. The construction and its proof are in
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
- T-matrices of order 3 (first rows, from base sequences BS(2, 1), built by the program): `+00`, `0+0`, `00+`, `000`.

## Search

None: the matrix is determined by q, the recipe of the Williamson-type matrices and the T-matrices. `python3
search/miyamoto_blocks.py certify 1193 3 "wt:cor25(one(), cor29c(149))"` writes the certificate.

## Status

**First explicit construction** ([MIYAMOTO_BLOCKS.md](../MIYAMOTO_BLOCKS.md#symmetric-williamson-type-matrices)). Table
1.53 of the Handbook of Combinatorial Designs (2nd ed., 2007) gives n = 3579 with t = 2, as known. The only construction
we traced is Cooper–Wallis with T-matrices of order 3 and Williamson-type matrices of order 1193 from Seberry–Yamada's
Table A.1 (1992), which obtains those only from Corollary 8.8(i) (key wl, which Seberry and Yamada credit to Miyamoto
(1991), with Williamson-type matrices of order 298 = (1193 − 1)/4 from their Corollary 25, key wo; their Remark 30 also
lists the variant with an Hadamard matrix of order 596 = (1193 − 1)/2). The printed proof of that corollary does not
give Williamson-type matrices ([NOTES.md](../../NOTES.md) §7). We found no other construction of Williamson-type
matrices of order 1193: Lemma 26 would need a symmetric quaternary Hadamard matrix of order 596 = 4 · 149 (for instance
from a repeated-block Williamson-type input), and we know none. The route here needs only symmetric ones of order (q −
1)/4 = 298, which Seberry and Yamada's corollaries give.

- Seberry and Yamada (J. Combin. Math. Combin. Comput. 7 (1990), Theorem 34 and Corollary 35) write Miyamoto's matrix as
  sixteen bordered blocks, but their multiplication theorems need circulant or type-one blocks (Theorem 4) or, for
  general blocks, sixteen matrices that exist for no odd order t > 1 (Theorem 21). With K the Williamson array of
  symmetric Williamson-type matrices, relabelled, the blocks commute and are transposed by one permutation, so
  T-matrices of odd order apply; we found this nowhere. Every ingredient was published before 2026.
- Paley's constructions do not apply: 4n − 1 = 14315 = 5 · 7 · 409 and 2n − 1 = 7157 = 17 · 421 are not prime powers.
- Miyamoto's theorem does not apply: it needs a symmetric conference matrix of order n + 1, which cannot exist for n ≡ 3
  (mod 4).
- Yamada–Kiyasu does not apply: n − 2 = 3577 = 7² · 73 is not a prime power.
- Cooper–Wallis needs Williamson-type matrices of order 1193 (only through Corollary 8.8(i), above) or T-matrices of
  order 1193; neither is known to us.
- The Paley-block route does not reach it: 4s − 1 = 4771 and 14315 are not prime powers ≡ 3 (mod 8) for the divisors s
  of 3579, and the remaining ones need T-matrices of order 3579 (s = 1) or 1193 (s = 3), which are not known to us.
- The Miyamoto-block route with a negaperiodic or periodic Golay pair of length 596 does not apply: none of that length
  is known to us.
- SageMath 10.8 (passagemath 10.8.13, October 2026) does not construct it: `hadamard_matrix(14316, existence=True)`
  returns Unknown.

We found no explicit Hadamard matrix of order 14316 before this one. The sources searched are listed in
[NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, g, the four symmetric Williamson-type matrices (upper triangles, base64), the T-matrices and the SHA-256 of the matrix |

The matrix itself is not stored; `verify/verify.py` rebuilds it, after checking that the four matrices are symmetric ±1
matrices that commute pairwise and have X₁² + X₂² + X₃² + X₄² = 4wI.

## Verification

```
python3 verify/verify.py 14316
```

SHA-256 of `H14316.csv`: `e284fa4c46ed0b230e6277a22969ade30027a021a164865f6a23c0c89d9f00cd`
