# Hadamard matrix of order 35004

35004 = 4 · 8751, with 8751 = 3 · 2917.

## Construction

Miyamoto's Hadamard matrix of order 11668 (q = 2917, a prime ≡ 1 mod 4) is built from Paley's conference matrix of order
2918 and a Hadamard matrix K of order 2916. Here K is the Williamson array of four symmetric, pairwise commuting ±1
matrices of order w = 729 (symmetric Williamson-type matrices) from Seberry and Yamada's constructions. Its four
quadrants then lie in a commutative algebra, so Miyamoto's matrix is a 4 × 4 array of bordered blocks of order 2917 that
commute and are transposed by one permutation. Multiplying it by T-matrices of order 3 in the Goethals–Seidel array
gives the Hadamard matrix of order 4 · 2917 · 3 = 35004. The construction and its proof are in
[MIYAMOTO_BLOCKS.md](../MIYAMOTO_BLOCKS.md#symmetric-williamson-type-matrices).

- Paley's conference matrix of order 2918, its points listed as ∞, 0, g^(2j), g^(2j+1) (j < 1458) with the primitive
  root g = 5.
- Symmetric Williamson-type matrices of order 729 = 3⁶ by Corollary 29(ii) of Seberry and Yamada (J. Combin. Math.
  Combin. Comput. 7 (1990)) over GF(729), with Paley's symmetric Hadamard matrix of order 364 (second construction, q =
  181). Seberry and Yamada list 729 among the orders this corollary gives (their Remark 30). K is their Williamson array
  [[X₁, X₂, X₃, X₄], [−X₂, X₁, −X₄, X₃], [−X₃, X₄, X₁, −X₂], [−X₄, −X₃, X₂, X₁]], relabelled so that the half swap of
  each quadrant becomes j ↦ 1 − j on ℤ_1458. `search/symwt.py` builds them (recipe `cor29h(729, paley2(181))`) and
  checks that they are symmetric, commute and have X₁² + X₂² + X₃² + X₄² = 2916I.
- T-matrices of order 3 (first rows, from base sequences BS(2, 1), built by the program): `+00`, `0+0`, `00+`, `000`.

## Search

None: the matrix is determined by q, the recipe of the Williamson-type matrices and the T-matrices. `python3
search/miyamoto_blocks.py certify 2917 3 "wt:cor29h(729, paley2(181))"` writes the certificate.

## Status

**First explicit construction** ([MIYAMOTO_BLOCKS.md](../MIYAMOTO_BLOCKS.md#symmetric-williamson-type-matrices)). Table
1.53 of the Handbook of Combinatorial Designs (2nd ed., 2007) gives n = 8751 with t = 3. Đoković (arXiv:1008.2043,
Proposition 3.3; Int. Math. Forum 6 (2011) 3061–3067) settles it with Williamson-type matrices of order 2917, which he
takes from Seberry–Yamada's Corollary 8.8(i) with Williamson-type matrices of order (q − 1)/4 = 729 from their Table
A.1. The printed proof of that corollary does not give Williamson-type matrices ([NOTES.md](../../NOTES.md) §7). For q =
2917, though, the corollary's input can be chosen with a repeated block: the symmetric Williamson-type matrices of order
729 used here have X₃ = X₄, and Lemma 26 with M₃ = [[Y, X], [X, −Z]], M₄ = [[Z, X], [X, −Y]] (the repeated-block form of
[NOTES.md](../../NOTES.md) §7) does give symmetric Williamson-type matrices of order 2917, checked exactly by `python3
search/symwt.py "cor29r(2917, cor29h(729, paley2(181)))"`. So Cooper–Wallis also gives this order from ingredients
published before 2026; that arrangement is ours, and we did not find it in the literature.

- Seberry and Yamada (J. Combin. Math. Combin. Comput. 7 (1990), Theorem 34 and Corollary 35) write Miyamoto's matrix as
  sixteen bordered blocks, but their multiplication theorems need circulant or type-one blocks (Theorem 4) or, for
  general blocks, sixteen matrices that exist for no odd order t > 1 (Theorem 21). With K the Williamson array of
  symmetric Williamson-type matrices, relabelled, the blocks commute and are transposed by one permutation, so
  T-matrices of odd order apply; we found this nowhere. Every ingredient was published before 2026.
- Paley's constructions do not apply: 4n − 1 = 35003 = 17 · 29 · 71 and 2n − 1 = 17501 = 11 · 37 · 43 are not prime
  powers.
- Miyamoto's theorem does not apply: it needs a symmetric conference matrix of order n + 1, which cannot exist for n ≡ 3
  (mod 4).
- Yamada–Kiyasu does not apply: n − 2 = 8749 = 13 · 673 is not a prime power.
- Cooper–Wallis applies with the Williamson-type matrices of order 2917 above and T-matrices of order 3; the other split
  needs T-matrices of order 2917, which are not known to us.
- The Paley-block route does not reach it: 4s − 1 = 11667 and 35003 are not prime powers ≡ 3 (mod 8) for the divisors s
  of 8751, and the remaining ones need T-matrices of order 8751 (s = 1) or 2917 (s = 3), which are not known to us.
- The Miyamoto-block route with a negaperiodic or periodic Golay pair of length 1458 does not apply: none of that length
  is known to us.
- SageMath 10.8 (passagemath 10.8.13, October 2026) does not construct it: `hadamard_matrix(35004, existence=True)`
  returns Unknown.

We found no explicit Hadamard matrix of order 35004 before this one. The sources searched are listed in
[NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, g, the four symmetric Williamson-type matrices (upper triangles, base64), the T-matrices and the SHA-256 of the matrix |

The matrix itself is not stored; `verify/verify.py` rebuilds it, after checking that the four matrices are symmetric ±1
matrices that commute pairwise and have X₁² + X₂² + X₃² + X₄² = 4wI.

## Verification

```
python3 verify/verify.py 35004
```

SHA-256 of `H35004.csv`: `b58d6fb5b0619a7f2611ffd8ac7150eca1af95ab8355b2a6082855593acb3590`
