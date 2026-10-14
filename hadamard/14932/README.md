# Order 14932

14932 = 4 · 3733, with 3733 prime and 3733 ≡ 1 (mod 4).

## Construction

Miyamoto's theorem: if q ≡ 1 (mod 4) is a prime power and a Hadamard matrix of order q − 1 exists, then a Hadamard
matrix of order 4q exists (M. Miyamoto, A construction of Hadamard matrices, J. Combin. Theory Ser. A 57 (1991)
86–108).

Here q = 3733. The inputs are

- Paley's symmetric conference matrix of order 3734: the matrix χ(y − x), x, y ∈ ℤ₃₇₃₃, with χ the quadratic
  character, bordered by a row and a column of ones;
- the Hadamard matrix of order 3732 = q − 1 in [3732](../3732), a skew Goethals–Seidel matrix over ℤ₉₃₃.

The blocks are arranged as in SageMath's `hadamard_matrix_miyamoto_construction`; `build_miyamoto` in
[`verify/verify.py`](../../verify/verify.py) spells out the layout. `certificate.json` contains q, the choice of
conference matrix, the complete certificate of the order-3732 matrix and the SHA-256 of the result.

The matrix (about 450 MB as CSV) is not stored. `python3 verify/verify.py 14932` rebuilds it, checks H Hᵀ = 14932 I
and compares the SHA-256, in about 2 minutes on 4 cores.

## Status

First explicit construction we know of, with one caveat.

- The order follows from Miyamoto's theorem only now that a Hadamard matrix of order 3732 exists; Cati–Pasechnik
  Table 4 lists 933(4), and Table 1.53 of the Handbook of Combinatorial Designs (2007) has 933 → 4.
- **Caveat.** Handbook Table 1.53 marks n = 3733 with "." (t = 2), that is, as already known. We found no construction
  behind that entry: Paley I and II fail (14931 = 3³ · 7 · 79, 7465 = 5 · 1493), Miyamoto needed H(3732), which the
  same table lists as unknown, Yamada needs q ≡ 1 (mod 8), Yamada–Kiyasu at q + 2 = 3733 needs q = 3731 = 7 · 13 · 41,
  and no Williamson-type or T-matrix construction of the needed orders is known. For five other primes q ≡ 1 (mod 4)
  below 10000 that the table marks "." without a visible route (3821, 4157, 5021, 5693, 7517), 2q − 3 is a prime
  ≡ 3 (mod 4), so Paley I gives H(2(q − 1)) rather than H(q − 1); that pattern does not cover 3733
  (2 · 3733 − 3 = 7463 = 17 · 439). We read the table from an unofficial copy, so the entry may also be a
  transcription error.
- Cati–Pasechnik and Seberry–Yamada (1992) stop at n ≤ 2999, so they do not cover 3733. No 2026 deposit, arXiv
  paper or Zenodo record we checked mentions 14932 or 3733 (audit of 2026-10-13 by a research agent).

## Verification

```
python3 verify/verify.py 14932
```

SHA-256 of `H14932.csv` (not stored): `1a0ce328dc4ca1013d286976b84610480bb09459ab3531184bf7042bf83fca9c`
