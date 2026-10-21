# Hadamard matrix of order 45044

45044 = 4 · 11261, with 11261 prime and 11261 ≡ 1 (mod 4).

## Construction

Miyamoto's theorem: if q ≡ 1 (mod 4) is a prime power and a Hadamard matrix of order q − 1 exists, then a Hadamard
matrix of order 4q exists (M. Miyamoto, A construction of Hadamard matrices, J. Combin. Theory Ser. A 57 (1991)
86–108).

Here q = 11261. The inputs are

- Paley's symmetric conference matrix of order 11262;
- the Hadamard matrix of order 11260 = q − 1 in [11260](../11260), from Paley's matrix of order 2252 in sixteen
  circulant blocks and T-matrices of order 5.

The blocks are arranged as in SageMath's `hadamard_matrix_miyamoto_construction`; `build_miyamoto` in
`verify/verify.py` spells out the layout. `certificate.json` contains q, the choice of conference matrix, the complete
certificate of the order-11260 matrix and the SHA-256 of the result.

## Status

**New, beyond the tables, but obtainable classically**. The order lies above the range of every table we know: the
Handbook's Table 1.53 stops at n = 9999.

The input of order 11260 = 4 · 2815 is one of the fourteen orders of [PALEY_BLOCKS.md](../PALEY_BLOCKS.md): Paley's
matrix of order 2252 in sixteen circulant blocks of order 563, times T-matrices of order 5. Cati–Pasechnik Table 4
lists 2815(3), but that matrix was obtainable from published theorems (Seberry–Yamada 1990, Theorem 4, with
Seberry–Yamada 1993, Theorem 1), so this order was too. We found no construction or statement of order 45044.

No other construction we know of reaches 45044:

- Paley I and II fail (45043 = 31 · 1453, 22521 = 3 · 7507). A conference matrix of order 22522 cannot exist, since
  22521 is not a sum of two squares.
- Yamada and Yamada–Kiyasu at q + 2 = 11261 need 11259 = 3⁴ · 139 to be a prime power.
- 11261 is prime, so a product theorem needs T-matrices or Williamson-type matrices of order 11261, which are not
  known to us; 11261 is not of Spence's form q² + q + 1.
- The arXiv sources of Hadamard-related papers that we collected (to 8 October 2026) do not mention 11261 or 45044.

The sources searched are listed in [NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, conference matrix, certificate of the order-11260 matrix, SHA-256 |

The matrix is not stored.

## Verification

```
python3 verify/verify.py 45044
```

This rebuilds the matrix, checks H Hᵀ = 45044 I and compares the SHA-256.

SHA-256 of `H45044.csv`: `dee9d1fe40f5974966b3007501d2e14649fb15afac09cf93af8cadf32a941a4f`
