# Hadamard matrix of order 21044

21044 = 4 · 5261, with 5261 prime and 5261 ≡ 1 (mod 4).

## Construction

Miyamoto's theorem: if q ≡ 1 (mod 4) is a prime power and a Hadamard matrix of order q − 1 exists, then a Hadamard
matrix of order 4q exists (M. Miyamoto, A construction of Hadamard matrices, J. Combin. Theory Ser. A 57 (1991)
86–108).

Here q = 5261. The inputs are

- Paley's symmetric conference matrix of order 5262;
- the Hadamard matrix of order 5260 = q − 1 in [5260](../5260), from Paley's matrix of order 1052 in sixteen
  circulant blocks and T-matrices of order 5.

The blocks are arranged as in SageMath's `hadamard_matrix_miyamoto_construction`; `build_miyamoto` in
`verify/verify.py` spells out the layout. `certificate.json` contains q, the choice of conference matrix, the complete
certificate of the order-5260 matrix and the SHA-256 of the result.

## Status

**New, beyond the tables, but obtainable classically**. In Table 1.53 of the Handbook of Combinatorial Designs (2nd
ed., 2007), n = 5261 has t = 6: the smallest known Hadamard matrix with odd part 5261 had order 2⁶ · 5261. The order
lies above the range of Cati–Pasechnik Table 4, which stops at n = 2999. Doubling gives the orders 2³ · 5261, 2⁴ ·
5261 and 2⁵ · 5261, which the Handbook also lists as unknown.

The input of order 5260 = 4 · 1315 is one of the fourteen orders of [PALEY_BLOCKS.md](../PALEY_BLOCKS.md): Paley's
matrix of order 1052 in sixteen circulant blocks of order 263, times T-matrices of order 5. Cati–Pasechnik Table 4
lists 1315(3), but that matrix was obtainable from published theorems (Seberry–Yamada 1990, Theorem 4, with
Seberry–Yamada 1993, Theorem 1), so this order was too. We found no construction or statement of order 21044.

No other construction we know of reaches 21044:

- Paley I and II fail (21043 = 11 · 1913, 10521 = 3² · 7 · 167). A conference matrix of order 10522 cannot exist,
  since 10521 is not a sum of two squares.
- Yamada and Yamada–Kiyasu at q + 2 = 5261 need 5259 = 3 · 1753 to be a prime power.
- 5261 is prime, so a product theorem needs T-matrices or Williamson-type matrices of order 5261, which are not known
  to us; 5261 is not of Spence's form q² + q + 1.
- The arXiv sources of Hadamard-related papers that we collected (to October 8, 2026) do not mention 5261 or 21044.

The sources searched are listed in [NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, conference matrix, certificate of the order-5260 matrix, SHA-256 |

The matrix is not stored.

## Verification

```
python3 verify/verify.py 21044
```

This rebuilds the matrix, checks H Hᵀ = 21044 I and compares the SHA-256.

SHA-256 of `H21044.csv`: `c723f88c80f676399687d297664f3efa1db9a61c0258dda02c6923e954af5f5d`
