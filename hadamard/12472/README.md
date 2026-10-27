# Hadamard matrix of order 12472

12472 = 8 · 1559, with 1559 a prime ≡ 3 (mod 4).

## Construction

Miyamoto's construction in quaternary form ([COMPLEX_MIYAMOTO.md](../COMPLEX_MIYAMOTO.md)), then Turyn's doubling.

1. **A quaternary Hadamard matrix of order 1558 = q − 1** from this repository
   ([complex-hadamard/C1558.json](../../complex-hadamard/C1558.json)). It is Szöllősi's array over ℤ₇₇₉ with a multiplier
   subgroup of order 60; Zaitsev's page has a matrix of this order too.
2. **A quaternary Hadamard matrix of order 4q = 6236.** Paley's skew conference matrix of order 1560, multiplied by i,
   is Hermitian. Its core on GF(1559)* goes into Miyamoto's array of sixteen bordered blocks together with the matrix
   of order 1558, and the border is chosen to match the row sums of the Hermitian core.
3. **Turyn's doubling** H = [[Hr + Hi, Hr − Hi], [−(Hr − Hi), Hr + Hi]] gives the real matrix of order 12472.

`verify/verify.py` rebuilds the matrix of order 1558 from its four sequences with `verify/verify_complex.py` and checks
it. It then builds the quaternary matrix of order 6236, checking its border identities, and finally checks
H Hᵀ = 12472 I.

## Search

None. The route is a theorem; the only search was the 4 × 4 border, solved once for all q
([COMPLEX_MIYAMOTO.md](../COMPLEX_MIYAMOTO.md)).

## Status

The literature check for this order and for the theorem is recorded in [NOVELTY.md](../../NOVELTY.md).
Cati–Pasechnik's Table 4 records 1559(4): order 16 · 1559 known, 8 · 1559 not. A Hadamard matrix of order
6236 = 4 · 1559 is unknown ([FRONTIER.md](../../FRONTIER.md)), so there is no doubling. Other routes:

- Paley: 12471 = 3 · 4157 and 6235 = 5 · 29 · 43 are not prime powers.
- Seberry–Yamada's order 8q needs Williamson-type matrices of order 779 = 19 · 41.
- Miyamoto's order 8q needs 2q − 3 = 3115 = 5 · 7 · 89 to be a prime power.
- The complex Golay decomposition 8(x + y), x + y = 1559, needs an odd complex Golay length x ∈ {1, 3, 5, 11, 13}, and
  1559 − x is not a complex Golay length.
- Kamali and Kharaghani's 8(g₁g + g₂) with their dihedral Golay lengths g = 7, 9, 15 does not reach 1559
  ([search/NOTES.md](../../search/NOTES.md) §104).

The quaternary matrix of order 6236 is an order where a quaternary Hadamard matrix exists but no real one is known.

## Files

| File | Contents |
|---|---|
| `certificate.json` | q and the complete certificate of the quaternary matrix of order 1558 |

The matrix (12472² entries) is not stored; `verify/verify.py` rebuilds it in about 40 seconds.

## Verification

```
python3 verify/verify.py 12472
```

SHA-256 of `H12472.csv`: `7be2f93a2f42a331dff0112674d677d2f0c3cbd064cde11e78b9f8d21e643c91`
