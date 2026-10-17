# Hadamard matrix of order 20796

20796 = 4 · 5199, with 5199 = 3 · 1733.

## Construction

Miyamoto's Hadamard matrix of order 6932 (q = 1733, a prime ≡ 1 mod 4) is built from Paley's conference matrix of order
1734 and a Hadamard matrix K of order 1732. Here K is the Williamson array of four symmetric, pairwise commuting ±1
matrices of order w = 433 (symmetric Williamson-type matrices) from Seberry and Yamada's constructions. Its four
quadrants then lie in a commutative algebra, so Miyamoto's matrix is a 4 × 4 array of bordered blocks of order 1733 that
commute and are transposed by one permutation. Multiplying it by T-matrices of order 3 in the Goethals–Seidel array
gives the Hadamard matrix of order 4 · 1733 · 3 = 20796. The construction and its proof are in
[MIYAMOTO_BLOCKS.md](../MIYAMOTO_BLOCKS.md#symmetric-williamson-type-matrices).

- Paley's conference matrix of order 1734, its points listed as ∞, 0, g^(2j), g^(2j+1) (j < 866) with the primitive root
  g = 2.
- Symmetric Williamson-type matrices of order 433 by Corollary 29(ii) of Seberry and Yamada (J. Combin. Math. Combin.
  Comput. 7 (1990)), with the symmetric Hadamard matrix H₁₀₈ ⊗ H₂ of order 216 (H₁₀₈ Paley's second construction for q =
  53). Seberry and Yamada list 433 among the orders this corollary gives (their Remark 30). K is their Williamson array
  [[X₁, X₂, X₃, X₄], [−X₂, X₁, −X₄, X₃], [−X₃, X₄, X₁, −X₂], [−X₄, −X₃, X₂, X₁]], relabelled so that the half swap of
  each quadrant becomes j ↦ 1 − j on ℤ_866. `search/symwt.py` builds them (recipe `cor29h(433, kron(paley2(53), h2()))`)
  and checks that they are symmetric, commute and have X₁² + X₂² + X₃² + X₄² = 1732I.
- T-matrices of order 3 (first rows, from base sequences BS(2, 1), built by the program): `+00`, `0+0`, `00+`, `000`.

## Search

None: the matrix is determined by q, the recipe of the Williamson-type matrices and the T-matrices. `python3
search/miyamoto_blocks.py certify 1733 3 "wt:cor29h(433, kron(paley2(53), h2()))"` writes the certificate.

## Status

**First explicit construction** ([MIYAMOTO_BLOCKS.md](../MIYAMOTO_BLOCKS.md#symmetric-williamson-type-matrices)). Table
1.53 of the Handbook of Combinatorial Designs (2nd ed., 2007) gives n = 5199 with t = 3, and Đoković (arXiv:0912.5091)
lists it in his set Δ of entries that should have had t = 2. His proof for Δ takes Williamson-type matrices from
Seberry–Yamada's Table A.1 (1992); for 5199 they are of order 1733, which Table A.1 obtains only from Corollary 8.8(i)
(key wl, with Williamson-type matrices of order 433 = (1733 − 1)/4). The printed proof of that corollary does not give
Williamson-type matrices ([NOTES.md](../../NOTES.md) §7). For q = 1733, though, the corollary's input can be chosen with
a repeated block: the symmetric Williamson-type matrices of order 433 used here have X₃ = X₄, and Lemma 26 with M₃ =
[[Y, X], [X, −Z]], M₄ = [[Z, X], [X, −Y]] (the repeated-block form of [NOTES.md](../../NOTES.md) §7) does give symmetric
Williamson-type matrices of order 1733, checked exactly by `python3 search/symwt.py "cor29r(1733, cor29h(433,
kron(paley2(53), h2())))"`. So Cooper–Wallis also gives this order from ingredients published before 2026; that
arrangement is ours, and we did not find it in the literature.

- Seberry and Yamada (J. Combin. Math. Combin. Comput. 7 (1990), Theorem 34 and Corollary 35) write Miyamoto's matrix as
  sixteen bordered blocks, but their multiplication theorems need circulant or type-one blocks (Theorem 4) or, for
  general blocks, sixteen matrices that exist for no odd order t > 1 (Theorem 21). With K the Williamson array of
  symmetric Williamson-type matrices, relabelled, the blocks commute and are transposed by one permutation, so
  T-matrices of odd order apply; we found this nowhere. Every ingredient was published before 2026.
- Paley's constructions do not apply: 4n − 1 = 20795 = 5 · 4159 and 2n − 1 = 10397 = 37 · 281 are not prime powers.
- Miyamoto's theorem does not apply: it needs a symmetric conference matrix of order n + 1, which cannot exist for n ≡ 3
  (mod 4).
- Yamada–Kiyasu would need a skew-Hadamard matrix of order 2600, which is not known to us.
- Cooper–Wallis applies with the Williamson-type matrices of order 1733 above and T-matrices of order 3; the other split
  needs T-matrices of order 1733, which are not known to us.
- The Paley-block route does not reach it: 4s − 1 = 6931 and 20795 are not prime powers ≡ 3 (mod 8) for the divisors s
  of 5199, and the remaining ones need T-matrices of order 5199 (s = 1) or 1733 (s = 3), which are not known to us.
- The Miyamoto-block route with a negaperiodic or periodic Golay pair of length 866 does not apply: none of that length
  is known to us.
- SageMath 10.8 (passagemath 10.8.13, October 2026) does not construct it: `hadamard_matrix(20796, existence=True)`
  returns Unknown.

We found no explicit Hadamard matrix of order 20796 before this one. The sources searched are listed in
[NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, g, the four symmetric Williamson-type matrices (upper triangles, base64), the T-matrices and the SHA-256 of the matrix |

The matrix itself is not stored; `verify/verify.py` rebuilds it, after checking that the four matrices are symmetric ±1
matrices that commute pairwise and have X₁² + X₂² + X₃² + X₄² = 4wI.

## Verification

```
python3 verify/verify.py 20796
```

SHA-256 of `H20796.csv`: `4cfd33ac62b77cbf3bad331849413d13aaf2b0ec118f2df8fd68c29aaa5eaa83`
