# Hadamard matrix of order 40244

40244 = 4 · 10061, with 10061 prime and 10061 ≡ 1 (mod 4).

## Construction

Miyamoto's theorem: if q ≡ 1 (mod 4) is a prime power and a Hadamard matrix of order q − 1 exists, then a Hadamard
matrix of order 4q exists (M. Miyamoto, A construction of Hadamard matrices, J. Combin. Theory Ser. A 57 (1991)
86–108).

Here q = 10061. The inputs are

- Paley's symmetric conference matrix of order 10062;
- the Hadamard matrix of order 10060 = q − 1 in [10060](../10060), from Paley's matrix of order 2012 in sixteen
  circulant blocks and T-matrices of order 5.

The blocks are arranged as in SageMath's `hadamard_matrix_miyamoto_construction`; `build_miyamoto` in
`verify/verify.py` spells out the layout. `certificate.json` contains q, the choice of conference matrix, the complete
certificate of the order-10060 matrix and the SHA-256 of the result.

## Status

**New, beyond the tables**, but obtainable classically. The order lies above the range of every table we know: the
Handbook's Table 1.53 stops at n = 9999.

The input of order 10060 = 4 · 2515 is one of the fourteen orders of [PALEY_BLOCKS.md](../PALEY_BLOCKS.md): Paley's
matrix of order 2012 in sixteen circulant blocks of order 503, times T-matrices of order 5. Cati–Pasechnik Table 4
lists 2515(3), but that matrix was obtainable from published theorems (Seberry–Yamada 1990, Theorem 4, with
Seberry–Yamada 1993, Theorem 1), so this order was too. We found no construction or statement of order 40244.

No other construction we know of reaches 40244:

- Paley I and II fail (40243 = 7 · 5749, 20121 = 3 · 19 · 353). A conference matrix of order 20122 cannot exist,
  since 20121 is not a sum of two squares.
- Yamada and Yamada–Kiyasu at q + 2 = 10061 need 10059 = 3 · 7 · 479 to be a prime power.
- 10061 is prime, so a product theorem needs T-matrices or Williamson-type matrices of order 10061, which are not
  known to us; 10061 is not of Spence's form q² + q + 1.
- The arXiv sources of Hadamard-related papers that we collected (to October 8, 2026) do not mention 10061 or 40244.

The sources searched are listed in [NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, conference matrix, certificate of the order-10060 matrix, SHA-256 |

The matrix is not stored.

## Verification

```
python3 verify/verify.py 40244
```

This rebuilds the matrix, checks H Hᵀ = 40244 I and compares the SHA-256.

SHA-256 of `H40244.csv`: `a6f052289c4714b1f3564452adc113af6d77e8896d4a33a9069399c3657d7e3e`
