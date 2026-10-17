# Hadamard matrix of order 18292

18292 = 4 · 4573, with 4573 = 17 · 269.

## Construction

Miyamoto's Hadamard matrix of order 1076 (q = 269, a prime ≡ 1 mod 4) is built from Paley's conference matrix of order
270 and a Hadamard matrix K of order 268. Here K is the Williamson array of four pairwise commuting symmetric matrices
of order w = 67, obtained from near-Williamson matrices. Its four quadrants then lie in a commutative algebra, so
Miyamoto's matrix is a 4 × 4 array of bordered blocks of order 269 that commute and are transposed by one permutation.
Multiplying it by T-matrices of order 17 in the Goethals–Seidel array gives the Hadamard matrix of order 4 · 269 · 17 =
18292. The construction and its proof are in [MIYAMOTO_BLOCKS.md](../MIYAMOTO_BLOCKS.md#near-williamson-matrices).

- Paley's conference matrix of order 270, its points listed as ∞, 0, g^(2j), g^(2j+1) (j < 134) with the primitive root
  g = 2.
- Near-Williamson matrices of order 67 printed in Table 7 of arXiv:2605.08661v2 (September 2026),
  [near-williamson/NW67.txt](../../near-williamson/NW67.txt): circulant ±1 matrices A, B, C, D with B, C, D symmetric
  and AAᵀ + BBᵀ + CCᵀ + DDᵀ = 268I. With R the back-identity, AR, B, C, D are symmetric and commute pairwise, and K is
  their Williamson array [[X₁, X₂, X₃, X₄], [−X₂, X₁, −X₄, X₃], [−X₃, X₄, X₁, −X₂], [−X₄, −X₃, X₂, X₁]], relabelled so
  that the half swap of each quadrant becomes j ↦ 1 − j on ℤ_134.
- T-matrices of order 17 (first rows, [t-matrices/T17.txt](../../t-matrices/T17.txt), from a Golay pair of length 16
  (Turyn)): `+0000000000000000`, `0+++-++-+00000000`, `000000000+++---+-`, `00000000000000000`.

## Search

None: the matrix is determined by q, the near-Williamson matrices and the T-matrices. `python3 search/miyamoto_blocks.py
certify 269 17 nw:near-williamson/NW67.txt t-matrices/T17.txt` writes the certificate.

## Status

**New, beyond the tables** ([MIYAMOTO_BLOCKS.md](../MIYAMOTO_BLOCKS.md#near-williamson-matrices)). Table 1.53 of the
Handbook of Combinatorial Designs (2nd ed., 2007) gives n = 4573 with t = 3: no Hadamard matrix of order 4n was known.
Đoković's updates (arXiv:0912.5091, arXiv:1008.2043) do not cover it, and we found no later construction. The
near-Williamson matrices of order 67 are from 2026, so the order is not obtainable classically by this route.

- Paley's constructions do not apply: 4n − 1 = 18291 = 3 · 7 · 13 · 67 and 2n − 1 = 9145 = 5 · 31 · 59 are not prime
  powers.
- Miyamoto's theorem does not apply: 4573 is not a prime power, and a symmetric conference matrix of order 4574 is not
  known to us.
- Yamada–Kiyasu does not apply: n − 2 = 4571 = 7 · 653 is not a prime power ≡ 1 (mod 4).
- Cooper–Wallis needs Williamson-type matrices of order 269 or T-matrices of order 269; neither is known to us.
- The Paley-block route does not reach it: 4s − 1 = 1075 and 18291 are not prime powers ≡ 3 (mod 8) for the other
  divisors s of 4573, and the remaining ones need T-matrices of order 4573 (s = 1) or 269 (s = 17), which are not known
  to us.
- The Miyamoto-block route with a negaperiodic or periodic Golay pair of length 134 does not apply: no negaperiodic pair
  of that length is known (Balonin–Đoković, list (31), quasi-Williamson matrices of order 67), and a periodic one cannot
  exist (268 = 2² · 67 is not a sum of two squares).
- SageMath 10.8 (passagemath 10.8.13, October 2026) does not construct it: `hadamard_matrix(18292, existence=True)`
  raises UnknownError in its Miyamoto test (`hadamard_matrix_miyamoto_construction`), since the existence of a symmetric
  conference matrix of order 4574 is Unknown to it.

We found no valid construction of order 18292 before this one. The sources searched are listed in
[NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, g, the near-Williamson matrices, the T-matrices and the SHA-256 of the matrix |

The matrix itself is not stored; `verify/verify.py` rebuilds it.

## Verification

```
python3 verify/verify.py 18292
```

SHA-256 of `H18292.csv`: `a60d2e0b4df3a9262dd02ef31361d05d0da1a923c8699fda36cc1d1a03cb5fcd`
