# Hadamard matrix of order 10444

10444 = 4 · 2611, with 2611 = 7 · 373.

## Construction

Miyamoto's Hadamard matrix of order 1492 (q = 373, a prime ≡ 1 mod 4) is built from Paley's conference matrix of order
374 and a Hadamard matrix K of order 372. Here K is the Williamson array of four pairwise commuting symmetric matrices
of order w = 93, obtained from near-Williamson matrices. Its four quadrants then lie in a commutative algebra, so
Miyamoto's matrix is a 4 × 4 array of bordered blocks of order 373 that commute and are transposed by one permutation.
Multiplying it by T-matrices of order 7 in the Goethals–Seidel array gives the Hadamard matrix of order 4 · 373 · 7 =
10444. The construction and its proof are in [MIYAMOTO_BLOCKS.md](../MIYAMOTO_BLOCKS.md#near-williamson-matrices).

- Paley's conference matrix of order 374, its points listed as ∞, 0, g^(2j), g^(2j+1) (j < 186) with the primitive root
  g = 2.
- Near-Williamson matrices of order 93 from the database of arXiv:2605.08661v2 (cs.uleth.ca/~zaitsevv/nw/NW93.txt),
  [near-williamson/NW93.txt](../../near-williamson/NW93.txt): circulant ±1 matrices A, B, C, D with B, C, D symmetric
  and AAᵀ + BBᵀ + CCᵀ + DDᵀ = 372I. With R the back-identity, AR, B, C, D are symmetric and commute pairwise, and K is
  their Williamson array [[X₁, X₂, X₃, X₄], [−X₂, X₁, −X₄, X₃], [−X₃, X₄, X₁, −X₂], [−X₄, −X₃, X₂, X₁]], relabelled so
  that the half swap of each quadrant becomes j ↦ 1 − j on ℤ_186.
- T-matrices of order 7 (first rows, from base sequences BS(4, 3)): `+000000`, `0++-000`, `0000+0+`, `00000+0`.

## Search

None: the matrix is determined by q, the near-Williamson matrices and the T-matrices. `python3 search/miyamoto_blocks.py
certify 373 7 nw:near-williamson/NW93.txt` writes the certificate.

## Status

**First explicit construction** ([MIYAMOTO_BLOCKS.md](../MIYAMOTO_BLOCKS.md#near-williamson-matrices)). Cati–Pasechnik's
Table 4 omits n = 2611: it inherits m = 2 from Seberry–Yamada's 2020 Table A.17, and their 1992 Table A.2 builds the
order (code o2) from Williamson-type matrices of order 373. Table A.1 obtains those only from Corollary 8.8(i), with
Williamson-type matrices of order 93. The printed proof of that corollary does not give Williamson-type matrices
([NOTES.md](../../NOTES.md) §7), and [FRONTIER.md](../../FRONTIER.md) listed the order as unresolved.

- Paley's constructions do not apply: 4n − 1 = 10443 = 3 · 59² and 2n − 1 = 5221 = 23 · 227 are not prime powers.
- Miyamoto's theorem does not apply: it needs a symmetric conference matrix of order n + 1, which cannot exist for n ≡ 3
  (mod 4).
- Yamada–Kiyasu does not apply: n − 2 = 2609 is a prime ≡ 1 (mod 8), and Yamada's theorem for that case needs a
  symmetric conference matrix of order 1306, which is not known to us.
- Cooper–Wallis needs Williamson-type matrices of order 373 or T-matrices of order 373; neither is known to us.
- The Paley-block route does not reach it: 4s − 1 = 1491 and 10443 are not prime powers ≡ 3 (mod 8) for the other
  divisors s of 2611, and the remaining ones need T-matrices of order 2611 (s = 1) or 373 (s = 7), which are not known
  to us.
- The Miyamoto-block route with a negaperiodic or periodic Golay pair of length 186 does not apply: no negaperiodic pair
  of that length is known (Balonin–Đoković, list (31), quasi-Williamson matrices of order 93), and a periodic one cannot
  exist (372 = 2² · 3 · 31 is not a sum of two squares).
- SageMath (passagemath 10.8) does not construct it: `hadamard_matrix(10444, existence=True)` returns Unknown.

We found no valid construction of order 10444 before this one. The sources searched are listed in
[NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, g, the near-Williamson matrices, the T-matrices and the SHA-256 of the matrix |

The matrix itself is not stored; `verify/verify.py` rebuilds it.

## Verification

```
python3 verify/verify.py 10444
```

SHA-256 of `H10444.csv`: `3070fdd4473e3a915139462c73506f3bfc3b2e9d6c495c6c060fad99b5e28aa0`
