# Hadamard matrix of order 8996

8996 = 4 · 2249, with 2249 = 13 · 173.

## Construction

Paley's skew-Hadamard matrix of order 692 (q = 691, a prime ≡ 3 mod 8) is a 4 × 4 array of circulant blocks of order
s = 173 once its points are listed along the four orbits of a torus element of PSL(2, 691) and rescaled by signs.
Multiplying it by T-matrices of order 13 (Seberry–Yamada 1990, Theorem 4) gives four circulant ±1 sequences of length
2249 with Σ PAF = 0 off 0, and the Goethals–Seidel array gives the Hadamard matrix. The construction and its proof are
in [PALEY_BLOCKS.md](../PALEY_BLOCKS.md).

- Torus element: x ↦ −1/(x + τ) with τ = 4; sign normalisation c = +1.
- T-matrices of order 13 (first rows, from base sequences BS(7, 6)): `+0+0+00000000`, `0+0+0-+000000`, `0000000+0000-`,
  `00000000++--0`.
- Row sums of the four sequences: −15, 29, 73, −51 ((−15)² + 29² + 73² + (−51)² = 8996).

## Search

None: the matrix is determined by q and the T-matrices. `python3 search/paley_blocks.py certify 691 13` writes the
certificate.

## Status

**First explicit construction** ([PALEY_BLOCKS.md](../PALEY_BLOCKS.md#orders-recorded-only-through-corollary-88i)).
Cati–Pasechnik's Table 4 omits n = 2249: it inherits m = 2 from Seberry–Yamada's 2020 Table A.17, and their 1992 Table
A.2 builds the order (code o2) from Williamson-type matrices of order 173. Table A.1 obtains Williamson-type matrices of
order 173 only from Corollary 8.8(i) (key wl, with Williamson-type matrices of order 43 = (173 − 1)/4), whose printed
proof does not give Williamson-type matrices ([NOTES.md](../../NOTES.md) §7). Table A.1 credits that case to Miyamoto
(J. Combin. Theory Ser. A 57 (1991)), whose paper was not accessible to us; we found no valid proof of it and no other
construction of Williamson-type matrices of order 173.

- Paley's constructions do not apply: 4n − 1 = 8995 = 5 · 7 · 257 and 2n − 1 = 4497 = 3 · 1499 are not prime powers.
- Miyamoto's theorem does not apply: it needs a symmetric conference matrix of order n + 1 = 2250, and none is known to
  us (2249 is not a prime power); Yamada–Kiyasu does not apply: n − 2 = 2247 = 3 · 7 · 107 is not a prime power.
- Cooper–Wallis does not apply otherwise: the other splits of n need Williamson-type matrices or T-matrices that are not
  known to us.
- SageMath's development version (October 17, 2026) does not construct it: `hadamard_matrix(8996, existence=True)`
  raises UnknownError in its Miyamoto test (no symmetric conference matrix of order 2250 is known to it), and returns
  Unknown when that Unknown is treated as False.

We found no valid construction of order 8996 before this one. The sources searched are listed in
[NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, τ, power, c, the T-matrices, the four sequences and the SHA-256 of the matrix |
| `sequences.txt` | the four sequences of length 2249 |

The matrix itself is not stored; `verify/verify.py` rebuilds it.

## Verification

```
python3 verify/verify.py 8996
```

SHA-256 of `H8996.csv`: `46e630fd9eb4fe5a9f14ae7dfb4e6899baa2a3dccd48f45823bbaf05e305845a`
