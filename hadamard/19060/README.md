# Hadamard matrix of order 19060

19060 = 4 · 4765, with 4765 = 5 · 953.

## Construction

Miyamoto's Hadamard matrix of order 3812 (q = 953, a prime ≡ 1 mod 4) is built from Paley's conference matrix of order
954 and a Hadamard matrix K of order 952. Here K is the Williamson array of four symmetric, pairwise commuting ±1
matrices of order w = 238 (symmetric Williamson-type matrices) from Seberry and Yamada's constructions. Its four
quadrants then lie in a commutative algebra, so Miyamoto's matrix is a 4 × 4 array of bordered blocks of order 953 that
commute and are transposed by one permutation. Multiplying it by T-matrices of order 5 in the Goethals–Seidel array
gives the Hadamard matrix of order 4 · 953 · 5 = 19060. The construction and its proof are in
[MIYAMOTO_BLOCKS.md](../MIYAMOTO_BLOCKS.md#symmetric-williamson-type-matrices).

- Paley's conference matrix of order 954, its points listed as ∞, 0, g^(2j), g^(2j+1) (j < 476) with the primitive root
  g = 3.
- Symmetric Williamson-type matrices of order 238 = 2 · 7 · 17 by Corollary 25 of Seberry and Yamada (J. Combin. Math.
  Combin. Comput. 7 (1990)) from Williamson matrices of order 7 (Turyn's family, q = 13, built as (I + β, I − β, α, α)
  from Yamada's symmetric circulants α, β with α² + β² = 13I) and symmetric circulant Williamson matrices of order 17
  (found by the program). K is their Williamson array [[X₁, X₂, X₃, X₄], [−X₂, X₁, −X₄, X₃], [−X₃, X₄, X₁, −X₂], [−X₄,
  −X₃, X₂, X₁]], relabelled so that the half swap of each quadrant becomes j ↦ 1 − j on ℤ_476. `search/symwt.py` builds
  them (recipe `cor25(turyn(13), circ(17))`) and checks that they are symmetric, commute and have X₁² + X₂² + X₃² + X₄²
  = 952I.
- T-matrices of order 5 (first rows, from base sequences BS(3, 2), built by the program): `++000`, `00+00`, `000+-`,
  `00000`.

## Search

None: the matrix is determined by q, the recipe of the Williamson-type matrices and the T-matrices. `python3
search/miyamoto_blocks.py certify 953 5 "wt:cor25(turyn(13), circ(17))"` writes the certificate.

## Status

**First explicit construction** ([MIYAMOTO_BLOCKS.md](../MIYAMOTO_BLOCKS.md#symmetric-williamson-type-matrices)). Table
1.53 of the Handbook of Combinatorial Designs (2nd ed., 2007) gives n = 4765 with t = 3, and Đoković (arXiv:0912.5091)
lists it in his set Δ of entries that should have had t = 2. His proof for Δ takes Williamson-type matrices from
Seberry–Yamada's Table A.1 (1992); for 4765 they are of order 953, which Table A.1 obtains only from Corollary 8.8(i)
(key wl, with an Hadamard matrix of order 476 = (953 − 1)/2 (Seberry–Yamada, Remark 30)). The printed proof of that
corollary does not give Williamson-type matrices ([NOTES.md](../../NOTES.md) §7). We found no other construction of
Williamson-type matrices of order 953. The route here needs only symmetric ones of order (q − 1)/4 = 238, which Seberry
and Yamada's corollaries give.

- Seberry and Yamada (J. Combin. Math. Combin. Comput. 7 (1990), Theorem 34 and Corollary 35) write Miyamoto's matrix as
  sixteen bordered blocks, but their multiplication theorems need circulant or type-one blocks (Theorem 4) or, for
  general blocks, sixteen matrices that exist for no odd order t > 1 (Theorem 21). With K the Williamson array of
  symmetric Williamson-type matrices, relabelled, the blocks commute and are transposed by one permutation, so
  T-matrices of odd order apply; we found this nowhere. Every ingredient was published before 2026.
- Paley's constructions do not apply: 4n − 1 = 19059 = 3 · 6353 and 2n − 1 = 9529 = 13 · 733 are not prime powers.
- Miyamoto's theorem does not apply: 4765 is not a prime power, and a symmetric conference matrix of order 4766 is not
  known to us.
- Yamada–Kiyasu does not apply: n − 2 = 4763 = 11 · 433 is not a prime power.
- Cooper–Wallis needs Williamson-type matrices of order 953 (only through Corollary 8.8(i), above) or T-matrices of
  order 953; neither is known to us.
- The Paley-block route does not reach it: 4s − 1 = 3811 and 19059 are not prime powers ≡ 3 (mod 8) for the divisors s
  of 4765, and the remaining ones need T-matrices of order 4765 (s = 1) or 953 (s = 5), which are not known to us.
- The Miyamoto-block route with a negaperiodic or periodic Golay pair of length 476 does not apply: none of that length
  is known to us; a periodic one cannot exist, since 952 = 2³ · 7 · 17 is not a sum of two squares.
- SageMath 10.8 (passagemath 10.8.13, October 2026) does not construct it: `hadamard_matrix(19060, existence=True)`
  raises UnknownError in its Miyamoto test (`hadamard_matrix_miyamoto_construction`), since the existence of a symmetric
  conference matrix of order 4766 is Unknown to it.

We found no explicit Hadamard matrix of order 19060 before this one. The sources searched are listed in
[NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, g, the four symmetric Williamson-type matrices (upper triangles, base64), the T-matrices and the SHA-256 of the matrix |

The matrix itself is not stored; `verify/verify.py` rebuilds it, after checking that the four matrices are symmetric ±1
matrices that commute pairwise and have X₁² + X₂² + X₃² + X₄² = 4wI.

## Verification

```
python3 verify/verify.py 19060
```

SHA-256 of `H19060.csv`: `80fbd73e8e45e2a108f81f509179ba1c9a82c7ca7d36062d99898c68f9061216`
