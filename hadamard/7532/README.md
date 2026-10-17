# Hadamard matrix of order 7532

7532 = 4 · 1883, with 1883 = 7 · 269.

## Construction

Miyamoto's Hadamard matrix of order 1076 (q = 269, a prime ≡ 1 mod 4) is built from Paley's conference matrix of order
270 and a Hadamard matrix K of order 268. Here K is the Williamson array of four pairwise commuting symmetric matrices
of order w = 67, obtained from near-Williamson matrices. Its four quadrants then lie in a commutative algebra, so
Miyamoto's matrix is a 4 × 4 array of bordered blocks of order 269 that commute and are transposed by one permutation.
Multiplying it by T-matrices of order 7 in the Goethals–Seidel array gives the Hadamard matrix of order 4 · 269 · 7 =
7532. The construction and its proof are in [MIYAMOTO_BLOCKS.md](../MIYAMOTO_BLOCKS.md#near-williamson-matrices).

- Paley's conference matrix of order 270, its points listed as ∞, 0, g^(2j), g^(2j+1) (j < 134) with the primitive root
  g = 2.
- Near-Williamson matrices of order 67 printed in Table 7 of arXiv:2605.08661v2 (September 2026),
  [near-williamson/NW67.txt](../../near-williamson/NW67.txt): circulant ±1 matrices A, B, C, D with B, C, D symmetric
  and AAᵀ + BBᵀ + CCᵀ + DDᵀ = 268I. With R the back-identity, AR, B, C, D are symmetric and commute pairwise, and K is
  their Williamson array [[X₁, X₂, X₃, X₄], [−X₂, X₁, −X₄, X₃], [−X₃, X₄, X₁, −X₂], [−X₄, −X₃, X₂, X₁]], relabelled so
  that the half swap of each quadrant becomes j ↦ 1 − j on ℤ_134.
- T-matrices of order 7 (first rows, from base sequences BS(4, 3)): `+000000`, `0++-000`, `0000+0+`, `00000+0`.

## Search

None: the matrix is determined by q, the near-Williamson matrices and the T-matrices. `python3 search/miyamoto_blocks.py
certify 269 7 nw:near-williamson/NW67.txt` writes the certificate.

## Status

**New** ([MIYAMOTO_BLOCKS.md](../MIYAMOTO_BLOCKS.md#near-williamson-matrices)). Cati–Pasechnik's Table 4 lists n = 1883
as open, and [FRONTIER.md](../../FRONTIER.md) listed it as unresolved. The near-Williamson matrices of order 67 are from
September 2026. Seberry–Yamada's Corollary 8.8(i) would turn them into Williamson-type matrices of order 269, and
Cooper–Wallis would then give this order. But the printed proof of that corollary does not give Williamson-type matrices
([NOTES.md](../../NOTES.md) §7), and we found no claim of n = 1883 through it.

- Paley's constructions do not apply: 4n − 1 = 7531 = 17 · 443 and 2n − 1 = 3765 = 3 · 5 · 251 are not prime powers.
- Miyamoto's theorem does not apply: it needs a symmetric conference matrix of order n + 1, which cannot exist for n ≡ 3
  (mod 4).
- Yamada–Kiyasu does not apply: n − 2 = 1881 = 3² · 11 · 19 is not a prime power ≡ 1 (mod 4).
- Cooper–Wallis needs Williamson-type matrices of order 269 or T-matrices of order 269; neither is known to us.
- The Paley-block route does not reach it: 4s − 1 = 1075 and 7531 are not prime powers ≡ 3 (mod 8) for the other
  divisors s of 1883, and the remaining ones need T-matrices of order 1883 (s = 1) or 269 (s = 7), which are not known
  to us.
- The Miyamoto-block route with a negaperiodic or periodic Golay pair of length 134 does not apply: no negaperiodic pair
  of that length is known (Balonin–Đoković, list (31), quasi-Williamson matrices of order 67), and a periodic one cannot
  exist (268 = 2² · 67 is not a sum of two squares).
- SageMath 10.8 (passagemath 10.8.13, October 2026) does not construct it: `hadamard_matrix(7532, existence=True)`
  returns Unknown.

We found no valid construction of order 7532 before this one. The sources searched are listed in
[NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, g, the near-Williamson matrices, the T-matrices and the SHA-256 of the matrix |

The matrix itself is not stored; `verify/verify.py` rebuilds it.

## Verification

```
python3 verify/verify.py 7532
```

SHA-256 of `H7532.csv`: `7dc4eaeb49f9846cc4839ffc9b9c8c2c3c8b57121f1ad4497da0b69062340e95`
