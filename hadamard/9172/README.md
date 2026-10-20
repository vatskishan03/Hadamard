# Hadamard matrix of order 9172

9172 = 4 · 2293, with 2293 prime and 2293 ≡ 1 (mod 4).

## Construction

Miyamoto's theorem: if q ≡ 1 (mod 4) is a prime power and a Hadamard matrix of order q − 1 exists, then a Hadamard
matrix of order 4q exists (M. Miyamoto, A construction of Hadamard matrices, J. Combin. Theory Ser. A 57 (1991) 86–108).

Here q = 2293. The inputs are

- Paley's symmetric conference matrix of order 2294: the matrix χ(y − x), x, y ∈ ℤ₂₂₉₃, with χ the quadratic character,
  bordered by a row and a column of ones;
- the Hadamard matrix of order 2292 = 4 · 573 from [2292](../2292) (Spence's bordered construction over ℤ₂ × ℤ₂ × ℤ₁₄₃
  with a spin structure).

`certificate.json` contains q, the choice of conference matrix, the complete certificate of the matrix of order 2292 and
the SHA-256 of the result.

## Status

**New.** The entry 2293(6) is unresolved in Cati–Pasechnik Table 4. So doubling this matrix also settles the orders
2³ · 2293 = 18344, 2⁴ · 2293 = 36688 and 2⁵ · 2293 = 73376.

- Paley's constructions do not apply: 4 · 2293 − 1 = 9171 = 3² · 1019 and 2 · 2293 − 1 = 4585 = 5 · 7 · 131.
- Miyamoto's theorem needed a Hadamard matrix of order 2292. Middletownbooks/Hadamard published one on 9 October 2026,
  before ours, and we learned of it on 19 October ([2292](../2292/README.md#status)). We found no matrix of order 9172
  published before ours of 14 October 2026.

The sources searched are listed in [NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, conference matrix, the certificate of the matrix of order 2292, SHA-256 |

The matrix is not stored (about 210 MB as CSV); `verify/verify.py 9172` rebuilds it.

## Verification

```
python3 verify/verify.py 9172
```

SHA-256 of `H9172.csv` (not stored): `25011d2ed180a0ba9ac2f98d381789657560450b2bd361eedc47f6b782b52915`