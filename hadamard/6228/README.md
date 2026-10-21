# Hadamard matrix of order 6228

6228 = 4 · 1557, with 1557 = 3² · 173.

## Construction

Paley's skew-Hadamard matrix of order 692 (q = 691, a prime ≡ 3 mod 8) is a 4 × 4 array of circulant blocks of order
s = 173 once its points are listed along the four orbits of a torus element of PSL(2, 691) and rescaled by signs.
Multiplying it by T-matrices of order 9 (Seberry–Yamada 1990, Theorem 4) gives four circulant ±1 sequences of length
1557 with Σ PAF = 0 off 0, and the Goethals–Seidel array gives the Hadamard matrix. The construction and its proof are
in [PALEY_BLOCKS.md](../PALEY_BLOCKS.md).

- Torus element: x ↦ −1/(x + τ) with τ = 4; sign normalisation c = +1.
- T-matrices of order 9 (first rows, from base sequences BS(5, 4)): `++0+00000`, `00+0-0000`, `00000+00-`, `000000+-0`.
- Row sums of the four sequences: 3, 27, 27, −69 (3² + 27² + 27² + (−69)² = 6228).

## Search

None: the matrix is determined by q and the T-matrices. `python3 search/paley_blocks.py certify 691 9` writes the
certificate.

## Status

**First explicit construction** ([PALEY_BLOCKS.md](../PALEY_BLOCKS.md#orders-recorded-only-through-corollary-88i)).
Cati–Pasechnik's Table 4 omits n = 1557 (m = 2, which by Cati–Pasechnik §6 and their footnote 11 comes from
Seberry–Yamada's 2020 Table A.17; we did not read that table), and their 1992 Table A.2 builds the order (code o2) from
Williamson-type matrices of order 173. Table A.1 obtains Williamson-type matrices of order 173 only from Corollary
8.8(i) (key wl, with Williamson-type matrices of order 43 = (173 − 1)/4), whose printed proof does not give
Williamson-type matrices ([NOTES.md](../../NOTES.md) §7). Table A.1 credits that case to Miyamoto (J. Combin. Theory
Ser. A 57 (1991)), whose paper was not accessible to us; we found no valid proof of it. Since 6 October 2026 there is
another construction: Verma's symmetric conference matrix of order 86 (arXiv:2610.08893, checked exactly here) gives
symmetric Williamson-type matrices of order 173 by Corollary 8.8(ii), and Cooper–Wallis then gives this order; we found
no statement of this.

- Paley's constructions do not apply: 4n − 1 = 6227 = 13 · 479 and 2n − 1 = 3113 = 11 · 283 are not prime powers.
- Miyamoto's theorem does not apply: it needs a symmetric conference matrix of order n + 1 = 1558, and none is known to
  us (1557 is not a prime power); Yamada–Kiyasu does not apply: n − 2 = 1555 = 5 · 311 is not a prime power.
- Cooper–Wallis does not apply otherwise: the other splits of n need Williamson-type matrices or T-matrices that are not
  known to us.
- SageMath's development version (17 October 2026) does not construct it: `hadamard_matrix(6228, existence=True)`
  raises UnknownError in its Miyamoto test (no symmetric conference matrix of order 1558 is known to it), and returns
  Unknown when that Unknown is treated as False.

We found no valid construction of order 6228 before this one. The sources searched are listed in
[NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, τ, power, c, the T-matrices, the four sequences and the SHA-256 of the matrix |
| `sequences.txt` | the four sequences of length 1557 |

The matrix itself is not stored; `verify/verify.py` rebuilds it.

## Verification

```
python3 verify/verify.py 6228
```

SHA-256 of `H6228.csv`: `fe2d880c2f21beed09e294a07190c8db2f1678fe366de364a71b7683aace75c7`
