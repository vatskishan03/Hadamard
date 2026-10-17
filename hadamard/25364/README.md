# Hadamard matrix of order 25364

25364 = 4 · 6341, with 6341 = 17 · 373.

## Construction

Miyamoto's Hadamard matrix of order 1492 (q = 373, a prime ≡ 1 mod 4) is built from Paley's conference matrix of order
374 and a Hadamard matrix K of order 372. Here K is the Williamson array of four pairwise commuting symmetric matrices
of order w = 93, obtained from near-Williamson matrices. Its four quadrants then lie in a commutative algebra, so
Miyamoto's matrix is a 4 × 4 array of bordered blocks of order 373 that commute and are transposed by one permutation.
Multiplying it by T-matrices of order 17 in the Goethals–Seidel array gives the Hadamard matrix of order 4 · 373 · 17 =
25364. The construction and its proof are in [MIYAMOTO_BLOCKS.md](../MIYAMOTO_BLOCKS.md#near-williamson-matrices).

- Paley's conference matrix of order 374, its points listed as ∞, 0, g^(2j), g^(2j+1) (j < 186) with the primitive root
  g = 2.
- Near-Williamson matrices of order 93 from the database of arXiv:2605.08661v2 (cs.uleth.ca/~zaitsevv/nw/NW93.txt),
  [near-williamson/NW93.txt](../../near-williamson/NW93.txt): circulant ±1 matrices A, B, C, D with B, C, D symmetric
  and AAᵀ + BBᵀ + CCᵀ + DDᵀ = 372I. With R the back-identity, AR, B, C, D are symmetric and commute pairwise, and K is
  their Williamson array [[X₁, X₂, X₃, X₄], [−X₂, X₁, −X₄, X₃], [−X₃, X₄, X₁, −X₂], [−X₄, −X₃, X₂, X₁]], relabelled so
  that the half swap of each quadrant becomes j ↦ 1 − j on ℤ_186.
- T-matrices of order 17 (first rows, [t-matrices/T17.txt](../../t-matrices/T17.txt), from a Golay pair of length 16
  (Turyn)): `+0000000000000000`, `0+++-++-+00000000`, `000000000+++---+-`, `00000000000000000`.

## Search

None: the matrix is determined by q, the near-Williamson matrices and the T-matrices. `python3 search/miyamoto_blocks.py
certify 373 17 nw:near-williamson/NW93.txt t-matrices/T17.txt` writes the certificate.

## Status

**First explicit construction** ([MIYAMOTO_BLOCKS.md](../MIYAMOTO_BLOCKS.md#near-williamson-matrices)). Table 1.53 of
the Handbook of Combinatorial Designs (2nd ed., 2007) gives n = 6341 with t = 8, and Đoković (arXiv:0912.5091) lists it
in his set Δ of entries that should have had t = 2. His proof for Δ takes Williamson-type matrices from Seberry–Yamada's
Table A.1 (1992); for 6341 they are of order 373, which Table A.1 obtains only from Corollary 8.8(i) (key wl,
Williamson-type matrices of order 93). The printed proof of that corollary does not give Williamson-type matrices
([NOTES.md](../../NOTES.md) §7), and we found no other construction of the order.

- Paley's constructions do not apply: 4n − 1 = 25363 = 13 · 1951 and 2n − 1 = 12681 = 3² · 1409 are not prime powers.
- Miyamoto's theorem does not apply: 6341 is not a prime power, and a symmetric conference matrix of order 6342 is not
  known to us.
- Yamada–Kiyasu does not apply: n − 2 = 6339 = 3 · 2113 is not a prime power ≡ 1 (mod 4).
- Cooper–Wallis needs Williamson-type matrices of order 373 or T-matrices of order 373; neither is known to us.
- The Paley-block route does not reach it: 4s − 1 = 1491 and 25363 are not prime powers ≡ 3 (mod 8) for the other
  divisors s of 6341, and the remaining ones need T-matrices of order 6341 (s = 1) or 373 (s = 17), which are not known
  to us.
- The Miyamoto-block route with a negaperiodic or periodic Golay pair of length 186 does not apply: no negaperiodic pair
  of that length is known (Balonin–Đoković, list (31), quasi-Williamson matrices of order 93), and a periodic one cannot
  exist (372 = 2² · 3 · 31 is not a sum of two squares).
- SageMath 10.8 (passagemath 10.8.13, October 2026) does not construct it: `hadamard_matrix(25364, existence=True)`
  raises UnknownError in its Miyamoto test (`hadamard_matrix_miyamoto_construction`), since the existence of a symmetric
  conference matrix of order 6342 is Unknown to it.

We found no valid construction of order 25364 before this one. The sources searched are listed in
[NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, g, the near-Williamson matrices, the T-matrices and the SHA-256 of the matrix |

The matrix itself is not stored; `verify/verify.py` rebuilds it.

## Verification

```
python3 verify/verify.py 25364
```

SHA-256 of `H25364.csv`: `ad33cd96edb0065bdbbb2db04cd0ee4853f57688719040533e04093f7886abfb`
