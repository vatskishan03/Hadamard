# Hadamard matrix of order 33860

33860 = 4 · 8465, with 8465 = 5 · 1693.

## Construction

Miyamoto's Hadamard matrix of order 6772 (q = 1693, a prime ≡ 1 mod 4) is built from Paley's conference matrix of order
1694 and a Hadamard matrix K of order 1692. Here K is the Williamson array of four symmetric, pairwise commuting ±1
matrices of order w = 423 (symmetric Williamson-type matrices) from Seberry and Yamada's constructions. Its four
quadrants then lie in a commutative algebra, so Miyamoto's matrix is a 4 × 4 array of bordered blocks of order 1693 that
commute and are transposed by one permutation. Multiplying it by T-matrices of order 5 in the Goethals–Seidel array
gives the Hadamard matrix of order 4 · 1693 · 5 = 33860. The construction and its proof are in
[MIYAMOTO_BLOCKS.md](../MIYAMOTO_BLOCKS.md#symmetric-williamson-type-matrices).

- Paley's conference matrix of order 1694, its points listed as ∞, 0, g^(2j), g^(2j+1) (j < 846) with the primitive root
  g = 2.
- Symmetric Williamson-type matrices of order 423 by Corollary 27(i) of Seberry and Yamada (J. Combin. Math. Combin.
  Comput. 7 (1990)) with q = 421: Lemma 26 with U = (I, 0, BR, BR), V = (M, N, R, R), B Paley's skew core of order 211,
  R the back-diagonal matrix and M, N Yamada's symmetric circulants of order 211 with M² + N² = 421I (Seberry and Yamada
  list 423 in their Remark 28). K is their Williamson array [[X₁, X₂, X₃, X₄], [−X₂, X₁, −X₄, X₃], [−X₃, X₄, X₁, −X₂],
  [−X₄, −X₃, X₂, X₁]], relabelled so that the half swap of each quadrant becomes j ↦ 1 − j on ℤ_846. `search/symwt.py`
  builds them (recipe `cor27i(421)`) and checks that they are symmetric, commute and have X₁² + X₂² + X₃² + X₄² = 1692I.
- T-matrices of order 5 (first rows, from base sequences BS(3, 2), built by the program): `++000`, `00+00`, `000+-`,
  `00000`.

## Search

None: the matrix is determined by q, the recipe of the Williamson-type matrices and the T-matrices. `python3
search/miyamoto_blocks.py certify 1693 5 "wt:cor27i(421)"` writes the certificate.

## Status

**First explicit construction** ([MIYAMOTO_BLOCKS.md](../MIYAMOTO_BLOCKS.md#symmetric-williamson-type-matrices)). Table
1.53 of the Handbook of Combinatorial Designs (2nd ed., 2007) gives n = 8465 with t = 3, and Đoković (arXiv:0912.5091)
lists it in his set Δ of entries that should have had t = 2. His proof for Δ takes Williamson-type matrices from
Seberry–Yamada's Table A.1 (1992); for 8465 they are of order 1693, which Table A.1 obtains only from Corollary 8.8(i)
(key wl, with Williamson-type matrices of order 423 = (1693 − 1)/4). The printed proof of that corollary does not give
Williamson-type matrices ([NOTES.md](../../NOTES.md) §7). For q = 1693, though, the corollary's input can be chosen with
a repeated block: the symmetric Williamson-type matrices of order 423 used here have X₃ = X₄, and Lemma 26 with M₃ =
[[Y, X], [X, −Z]], M₄ = [[Z, X], [X, −Y]] (the repeated-block form of [NOTES.md](../../NOTES.md) §7) does give symmetric
Williamson-type matrices of order 1693, checked exactly by `python3 search/symwt.py "cor29r(1693, cor27i(421))"`. So
Cooper–Wallis also gives this order from ingredients published before 2026; we found that arrangement in none of the
sources we could read (Miyamoto's paper of 1991, to which Seberry and Yamada credit key wl, was not accessible).

- Seberry and Yamada (J. Combin. Math. Combin. Comput. 7 (1990), Theorem 34 and Corollary 35) write Miyamoto's matrix as
  sixteen bordered blocks, but their multiplication theorems need circulant or type-one blocks (Theorem 4) or, for
  general blocks, sixteen matrices that exist for no odd order t > 1 (Theorem 21). With K the Williamson array of
  symmetric Williamson-type matrices, relabelled, the blocks commute and are transposed by one permutation, so
  T-matrices of odd order apply; we found this nowhere. Every ingredient was published before 2026.
- Paley's constructions do not apply: 4n − 1 = 33859 = 7² · 691 and 2n − 1 = 16929 = 3⁴ · 11 · 19 are not prime powers.
- Miyamoto's theorem does not apply: 8465 is not a prime power, and a symmetric conference matrix of order 8466 is not
  known to us.
- Yamada–Kiyasu does not apply: n − 2 = 8463 = 3 · 7 · 13 · 31 is not a prime power.
- Cooper–Wallis applies with the Williamson-type matrices of order 1693 above and T-matrices of order 5; the other split
  needs T-matrices of order 1693, which are not known to us.
- The Paley-block route does not reach it: 4s − 1 = 6771 and 33859 are not prime powers ≡ 3 (mod 8) for the divisors s
  of 8465, and the remaining ones need T-matrices of order 8465 (s = 1) or 1693 (s = 5), which are not known to us.
- The Miyamoto-block route with a negaperiodic or periodic Golay pair of length 846 does not apply: none of that length
  is known to us; a periodic one cannot exist, since 1692 = 2² · 3² · 47 is not a sum of two squares.
- SageMath 10.8 (passagemath 10.8.13, October 2026) does not construct it: `hadamard_matrix(33860, existence=True)`
  raises UnknownError in its Miyamoto test (`hadamard_matrix_miyamoto_construction`), since the existence of a symmetric
  conference matrix of order 8466 is Unknown to it.

We found no explicit Hadamard matrix of order 33860 before this one. The sources searched are listed in
[NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, g, the four symmetric Williamson-type matrices (upper triangles, base64), the T-matrices and the SHA-256 of the matrix |

The matrix itself is not stored; `verify/verify.py` rebuilds it, after checking that the four matrices are symmetric ±1
matrices that commute pairwise and have X₁² + X₂² + X₃² + X₄² = 4wI.

## Verification

```
python3 verify/verify.py 33860
```

SHA-256 of `H33860.csv`: `117de87a9dc7452b78e23e99f074dc7dd2ac8f28058622a846fd0118ecc28031`
