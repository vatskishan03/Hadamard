# Hadamard matrix of order 32748

32748 = 4 · 8187, with 8187 = 3 · 2729.

## Construction

Miyamoto's Hadamard matrix of order 10916 (q = 2729, a prime ≡ 1 mod 4) is built from Paley's conference matrix of order
2730 and a Hadamard matrix K of order 2728. Here K is the Williamson array of four symmetric, pairwise commuting ±1
matrices of order w = 682 (symmetric Williamson-type matrices) from Seberry and Yamada's constructions. Its four
quadrants then lie in a commutative algebra, so Miyamoto's matrix is a 4 × 4 array of bordered blocks of order 2729 that
commute and are transposed by one permutation. Multiplying it by T-matrices of order 3 in the Goethals–Seidel array
gives the Hadamard matrix of order 4 · 2729 · 3 = 32748. The construction and its proof are in
[MIYAMOTO_BLOCKS.md](../MIYAMOTO_BLOCKS.md#symmetric-williamson-type-matrices).

- Paley's conference matrix of order 2730, its points listed as ∞, 0, g^(2j), g^(2j+1) (j < 1364) with the primitive
  root g = 3.
- Symmetric Williamson-type matrices of order 682 = 2 · 11 · 31 by Corollary 25 of Seberry and Yamada (J. Combin. Math.
  Combin. Comput. 7 (1990)) from symmetric circulant Williamson matrices of order 11 (found by the program) and
  Williamson matrices of order 31 (Turyn's family, q = 61, built as (I + β, I − β, α, α) from Yamada's symmetric
  circulants α, β with α² + β² = 61I). K is their Williamson array [[X₁, X₂, X₃, X₄], [−X₂, X₁, −X₄, X₃], [−X₃, X₄, X₁,
  −X₂], [−X₄, −X₃, X₂, X₁]], relabelled so that the half swap of each quadrant becomes j ↦ 1 − j on ℤ_1364.
  `search/symwt.py` builds them (recipe `cor25(circ(11), turyn(61))`) and checks that they are symmetric, commute and
  have X₁² + X₂² + X₃² + X₄² = 2728I.
- T-matrices of order 3 (first rows, from base sequences BS(2, 1), built by the program): `+00`, `0+0`, `00+`, `000`.

## Search

None: the matrix is determined by q, the recipe of the Williamson-type matrices and the T-matrices. `python3
search/miyamoto_blocks.py certify 2729 3 "wt:cor25(circ(11), turyn(61))"` writes the certificate.

## Status

**First explicit construction** ([MIYAMOTO_BLOCKS.md](../MIYAMOTO_BLOCKS.md#symmetric-williamson-type-matrices)). Table
1.53 of the Handbook of Combinatorial Designs (2nd ed., 2007) gives n = 8187 with t = 3. Đoković (arXiv:1008.2043,
Proposition 3.2; Int. Math. Forum 6 (2011) 3061–3067) settles it with Williamson-type matrices of order 2729, which he
takes from Seberry–Yamada's Corollary 8.8(i) with an Hadamard matrix of order (q − 1)/2 = 1364. The printed proof of
that corollary gives pairwise amicable matrices only for a symmetric input ([NOTES.md](../../NOTES.md) §7), and no
symmetric Hadamard matrix of order 1364 is known to us. Cooper–Wallis also gives this order from ingredients published
before 2026: Lemma 26 needs only a symmetric quaternary Hadamard matrix of order (q − 1)/2 = 1364 = 22 · 62, the
Kronecker product of the repeated-block arrangement of the symmetric Williamson-type matrices of order 11 from Corollary
31 and Paley's conference matrix of order 62 is one, and
`python3 search/symwt.py "cor29q(2729, qwt(cor31(5)), qconf(61))"` gives symmetric Williamson-type matrices of order
2729 (checked). The route here needs only symmetric ones of order (q − 1)/4 = 682, which Seberry and Yamada's
corollaries give.

- Seberry and Yamada (J. Combin. Math. Combin. Comput. 7 (1990), Theorem 34 and Corollary 35) write Miyamoto's matrix as
  sixteen bordered blocks, but their multiplication theorems need circulant or type-one blocks (Theorem 4) or, for
  general blocks, sixteen matrices that exist for no odd order t > 1 (Theorem 21). With K the Williamson array of
  symmetric Williamson-type matrices, relabelled, the blocks commute and are transposed by one permutation, so
  T-matrices of odd order apply; we found this nowhere. Every ingredient was published before 2026.
- Paley's constructions do not apply: 4n − 1 = 32747 = 11 · 13 · 229 and 2n − 1 = 16373 = 7 · 2339 are not prime powers.
- Miyamoto's theorem does not apply: it needs a symmetric conference matrix of order n + 1, which cannot exist for n ≡ 3
  (mod 4).
- Yamada–Kiyasu does not apply: n − 2 = 8185 = 5 · 1637 is not a prime power.
- Cooper–Wallis applies with the Williamson-type matrices of order 2729 of the Kronecker route above and T-matrices of
  order 3; the other split needs T-matrices of order 2729, which are not known to us.
- The Paley-block route does not reach it: 4s − 1 = 10915 and 32747 are not prime powers ≡ 3 (mod 8) for the divisors s
  of 8187, and the remaining ones need T-matrices of order 8187 (s = 1) or 2729 (s = 3), which are not known to us.
- The Miyamoto-block route with a negaperiodic or periodic Golay pair of length 1364 does not apply: none of that length
  is known to us; a periodic one cannot exist, since 2728 = 2³ · 11 · 31 is not a sum of two squares.
- SageMath 10.8 (passagemath 10.8.13, October 2026) does not construct it: `hadamard_matrix(32748, existence=True)`
  returns Unknown.

We found no explicit Hadamard matrix of order 32748 before this one. The sources searched are listed in
[NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, g, the four symmetric Williamson-type matrices (upper triangles, base64), the T-matrices and the SHA-256 of the matrix |

The matrix itself is not stored; `verify/verify.py` rebuilds it, after checking that the four matrices are symmetric ±1
matrices that commute pairwise and have X₁² + X₂² + X₃² + X₄² = 4wI.

## Verification

```
python3 verify/verify.py 32748
```

SHA-256 of `H32748.csv`: `983b1417a597807f55f4f75b940f7480523292ecae98b37959a4edafe797896d`
