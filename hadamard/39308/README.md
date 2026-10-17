# Hadamard matrix of order 39308

39308 = 4 · 9827, with 9827 = 31 · 317.

## Construction

Miyamoto's Hadamard matrix of order 1268 (q = 317, a prime ≡ 1 mod 4) is built from Paley's conference matrix of order
318 and a Hadamard matrix K of order 316. Here K is the Williamson array of four pairwise commuting symmetric matrices
of order w = 79, obtained from near-Williamson matrices. Its four quadrants then lie in a commutative algebra, so
Miyamoto's matrix is a 4 × 4 array of bordered blocks of order 317 that commute and are transposed by one permutation.
Multiplying it by T-matrices of order 31 in the Goethals–Seidel array gives the Hadamard matrix of order 4 · 317 · 31 =
39308. The construction and its proof are in [MIYAMOTO_BLOCKS.md](../MIYAMOTO_BLOCKS.md#near-williamson-matrices).

- Paley's conference matrix of order 318, its points listed as ∞, 0, g^(2j), g^(2j+1) (j < 158) with the primitive root
  g = 2.
- Near-Williamson matrices of order 79 from the database of arXiv:2605.08661v2 (cs.uleth.ca/~zaitsevv/nw/NW79.txt),
  [near-williamson/NW79.txt](../../near-williamson/NW79.txt): circulant ±1 matrices A, B, C, D with B, C, D symmetric
  and AAᵀ + BBᵀ + CCᵀ + DDᵀ = 316I. With R the back-identity, AR, B, C, D are symmetric and commute pairwise, and K is
  their Williamson array [[X₁, X₂, X₃, X₄], [−X₂, X₁, −X₄, X₃], [−X₃, X₄, X₁, −X₂], [−X₄, −X₃, X₂, X₁]], relabelled so
  that the half swap of each quadrant becomes j ↦ 1 − j on ℤ_158.
- T-matrices of order 31 (first rows, [t-matrices/T31.txt](../../t-matrices/T31.txt), from `search/tmat_mult.py 31 5`):
  `+0000000--0000-0000000000000000`, `0+000+00000+00000-0000--++0+000`, `00+-000000+00-0--0-+00000000-00`,
  `0000+0-+0000+0000000++0000-00+-`.

## Search

None: the matrix is determined by q, the near-Williamson matrices and the T-matrices. `python3 search/miyamoto_blocks.py
certify 317 31 nw:near-williamson/NW79.txt t-matrices/T31.txt` writes the certificate.

## Status

**Known** ([MIYAMOTO_BLOCKS.md](../MIYAMOTO_BLOCKS.md#near-williamson-matrices)); the matrix here is independent. Table
1.53 of the Handbook of Combinatorial Designs (2nd ed., 2007) gives n = 9827 with t = 3, but Đoković (arXiv:0912.5091)
lists it in his set Δ of entries that should have had t = 2: Cooper–Wallis with T-matrices of order 31 and
Williamson-type matrices of order 317 from Seberry–Yamada's Table A.1 (1992), which Corollary 8.8(ii) gives validly
from Paley's symmetric conference matrix of order 158 (157 is a prime ≡ 1 mod 4). The route here was open classically
too: Turyn's Williamson matrices of order 79 = (157 + 1)/2 (J. Combin. Theory Ser. A 12 (1972)) satisfy the same
hypotheses as the near-Williamson matrices, and Paley's negaperiodic Golay pair of length 158 (Balonin–Đoković,
Proposition 6, p = 157) gives the Golay-pair version of [MIYAMOTO_BLOCKS.md](../MIYAMOTO_BLOCKS.md). An earlier version
of this README labelled the order New, beyond the tables; the audit of October 17, 2026 found it in Δ.

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, g, the near-Williamson matrices, the T-matrices and the SHA-256 of the matrix |

The matrix itself is not stored; `verify/verify.py` rebuilds it.

## Verification

```
python3 verify/verify.py 39308
```

SHA-256 of `H39308.csv`: `06c84a1f171821310315cba5e7421b9c34882728a8a43f93038e8ec82f7b32bd`
