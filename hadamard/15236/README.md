# Hadamard matrix of order 15236

15236 = 4 · 3809, with 3809 = 13 · 293.

## Construction

Paley's skew-Hadamard matrix of order 1172 (q = 1171, a prime ≡ 3 mod 8) is a 4 × 4 array of circulant blocks of order
s = 293 once its points are listed along the four orbits of a torus element of PSL(2, 1171) and rescaled by signs.
Multiplying it by T-matrices of order 13 (Seberry–Yamada 1990, Theorem 4) gives four circulant ±1 sequences of length
3809 with Σ PAF = 0 off 0, and the Goethals–Seidel array gives the Hadamard matrix. The construction and its proof are
in [PALEY_BLOCKS.md](../PALEY_BLOCKS.md).

- Torus element: x ↦ −1/(x + τ) with τ = 4; sign normalisation c = −1.
- T-matrices of order 13 (first rows, from base sequences BS(7, 6)): `+0+0+00000000`, `0+0+0-+000000`, `0000000+0000-`,
  `00000000++--0`.
- Row sums of the four sequences: −31, 53, −21, −105 ((−31)² + 53² + (−21)² + (−105)² = 15236).

## Search

None: the matrix is determined by q and the T-matrices. `python3 search/paley_blocks.py certify 1171 13` writes the
certificate.

## Status

**First explicit construction** ([PALEY_BLOCKS.md](../PALEY_BLOCKS.md#orders-recorded-only-through-corollary-88i)).
Table 1.53 of the Handbook of Combinatorial Designs (2nd ed., 2007) gives n = 3809 with t = 2, as known. The only
construction we traced uses Williamson-type matrices of order 293 from Seberry–Yamada's Table A.1 (1992). Table A.1
obtains Williamson-type matrices of order 293 only from Corollary 8.8(i) (key wl, with Williamson-type matrices of order
73 = (293 − 1)/4), whose printed proof does not give Williamson-type matrices ([NOTES.md](../../NOTES.md) §7). Table A.1
credits that case to Miyamoto (J. Combin. Theory Ser. A 57 (1991)), whose paper was not accessible to us; we found no
valid proof of it and no other construction of Williamson-type matrices of order 293.

- Paley's constructions do not apply: 4n − 1 = 15235 = 5 · 11 · 277 and 2n − 1 = 7617 = 3 · 2539 are not prime powers.
- Miyamoto's theorem does not apply: it needs a symmetric conference matrix of order n + 1 = 3810, and none is known to
  us (3809 is not a prime power); Yamada–Kiyasu does not apply: n − 2 = 3807 = 3⁴ · 47 is not a prime power.
- Cooper–Wallis does not apply otherwise: the other splits of n need Williamson-type matrices or T-matrices that are not
  known to us.
- SageMath's development version (October 17, 2026) does not construct it: `hadamard_matrix(15236, existence=True)`
  raises UnknownError in its Miyamoto test (no symmetric conference matrix of order 3810 is known to it), and returns
  Unknown when that Unknown is treated as False.

We found no valid construction of order 15236 before this one. The sources searched are listed in
[NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, τ, power, c, the T-matrices, the four sequences and the SHA-256 of the matrix |
| `sequences.txt` | the four sequences of length 3809 |

The matrix itself is not stored; `verify/verify.py` rebuilds it.

## Verification

```
python3 verify/verify.py 15236
```

SHA-256 of `H15236.csv`: `76f4668c0a44df9724a6a91428585b55648ef4f0caf925b37474ed5c10445612`
