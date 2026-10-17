# Hadamard matrix of order 35436

35436 = 4 · 8859, with 8859 = 3 · 2953.

## Construction

Miyamoto's Hadamard matrix of order 11812 (q = 2953, a prime ≡ 1 mod 4) is built from Paley's conference matrix of order
2954 and a Hadamard matrix K of order 2952. Here K is the Williamson array of four symmetric, pairwise commuting ±1
matrices of order w = 738 (symmetric Williamson-type matrices) from Seberry and Yamada's constructions. Its four
quadrants then lie in a commutative algebra, so Miyamoto's matrix is a 4 × 4 array of bordered blocks of order 2953 that
commute and are transposed by one permutation. Multiplying it by T-matrices of order 3 in the Goethals–Seidel array
gives the Hadamard matrix of order 4 · 2953 · 3 = 35436. The construction and its proof are in
[MIYAMOTO_BLOCKS.md](../MIYAMOTO_BLOCKS.md#symmetric-williamson-type-matrices).

- Paley's conference matrix of order 2954, its points listed as ∞, 0, g^(2j), g^(2j+1) (j < 1476) with the primitive
  root g = 13.
- Symmetric Williamson-type matrices of order 738 = 2 · 9 · 41 by Corollary 25 of Seberry and Yamada (J. Combin. Math.
  Combin. Comput. 7 (1990)) from Williamson matrices of order 9 (Turyn's family, q = 17) and symmetric circulant
  Williamson matrices of order 41 constant on the orbits of ⟨9, −1⟩ (found by the program). K is their Williamson array
  [[X₁, X₂, X₃, X₄], [−X₂, X₁, −X₄, X₃], [−X₃, X₄, X₁, −X₂], [−X₄, −X₃, X₂, X₁]], relabelled so that the half swap of
  each quadrant becomes j ↦ 1 − j on ℤ_1476. `search/symwt.py` builds them (recipe `cor25(turyn(17), circm(41, 9))`) and
  checks that they are symmetric, commute and have X₁² + X₂² + X₃² + X₄² = 2952I.
- T-matrices of order 3 (first rows, from base sequences BS(2, 1), built by the program): `+00`, `0+0`, `00+`, `000`.

## Search

None: the matrix is determined by q, the recipe of the Williamson-type matrices and the T-matrices. `python3
search/miyamoto_blocks.py certify 2953 3 "wt:cor25(turyn(17), circm(41, 9))"` writes the certificate.

## Status

**First explicit construction** ([MIYAMOTO_BLOCKS.md](../MIYAMOTO_BLOCKS.md#symmetric-williamson-type-matrices)). Table
1.53 of the Handbook of Combinatorial Designs (2nd ed., 2007) gives n = 8859 with t = 3. Đoković (arXiv:1008.2043,
Proposition 3.2; Int. Math. Forum 6 (2011) 3061–3067) settles it with Williamson-type matrices of order 2953, which he
takes from Seberry–Yamada's Corollary 8.8(i) with an Hadamard matrix of order (q − 1)/2 = 1476. The printed proof of
that corollary gives pairwise amicable matrices only for a symmetric input ([NOTES.md](../../NOTES.md) §7), and no
symmetric Hadamard matrix of order 1476 is known to us. Cooper–Wallis also gives this order from ingredients published
before 2026 (found by the audit of October 17, 2026): Lemma 26 needs only a symmetric quaternary Hadamard matrix of
order (q − 1)/2 = 1476 = 18 · 82, the Kronecker product of Paley's conference matrices of orders 18 and 82 (over GF(81))
is one, and `python3 search/symwt.py "cor29q(2953, qconf(17), qconf(81))"` gives symmetric Williamson-type matrices of
order 2953 (checked). The route here needs only symmetric ones of order (q − 1)/4 = 738, which Seberry and Yamada's
corollaries give.

- Seberry and Yamada (J. Combin. Math. Combin. Comput. 7 (1990), Theorem 34 and Corollary 35) write Miyamoto's matrix as
  sixteen bordered blocks, but their multiplication theorems need circulant or type-one blocks (Theorem 4) or, for
  general blocks, sixteen matrices that exist for no odd order t > 1 (Theorem 21). With K the Williamson array of
  symmetric Williamson-type matrices, relabelled, the blocks commute and are transposed by one permutation, so
  T-matrices of odd order apply; we found this nowhere. Every ingredient was published before 2026.
- Paley's constructions do not apply: 4n − 1 = 35435 = 5 · 19 · 373 and 2n − 1 = 17717 = 7 · 2531 are not prime powers.
- Miyamoto's theorem does not apply: it needs a symmetric conference matrix of order n + 1, which cannot exist for n ≡ 3
  (mod 4).
- Yamada–Kiyasu does not apply: n − 2 = 8857 = 17 · 521 is not a prime power.
- Cooper–Wallis applies with the Williamson-type matrices of order 2953 of the Kronecker route above and T-matrices of
  order 3; the other split needs T-matrices of order 2953, which are not known to us.
- The Paley-block route does not reach it: 4s − 1 = 11811 and 35435 are not prime powers ≡ 3 (mod 8) for the divisors s
  of 8859, and the remaining ones need T-matrices of order 8859 (s = 1) or 2953 (s = 3), which are not known to us.
- The Miyamoto-block route with a negaperiodic or periodic Golay pair of length 1476 does not apply: none of that length
  is known to us.
- SageMath 10.8 (passagemath 10.8.13, October 2026) does not construct it: `hadamard_matrix(35436, existence=True)`
  returns Unknown.

We found no explicit Hadamard matrix of order 35436 before this one. The sources searched are listed in
[NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, g, the four symmetric Williamson-type matrices (upper triangles, base64), the T-matrices and the SHA-256 of the matrix |

The matrix itself is not stored; `verify/verify.py` rebuilds it, after checking that the four matrices are symmetric ±1
matrices that commute pairwise and have X₁² + X₂² + X₃² + X₄² = 4wI.

## Verification

```
python3 verify/verify.py 35436
```

SHA-256 of `H35436.csv`: `4f5489fbddf4cfccd6fb5bbca40e40e5c9422210c158ea05936a7615d1d936dd`
