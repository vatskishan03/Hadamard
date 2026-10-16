# Hadamard matrix of order 6356

6356 = 4 · 1589, with 1589 = 7 · 227.

## Construction

Paley's skew-Hadamard matrix of order 908 (q = 907, a prime ≡ 3 mod 8) is a 4 × 4 array of circulant blocks of
order s = 227 once its points are listed along the four orbits of a torus element of PSL(2, 907) and rescaled by
signs. Multiplying it by T-matrices of order 7 (Seberry–Yamada 1990, Theorem 4) gives four circulant ±1 sequences
of length 1589 with Σ PAF = 0 off 0, and the Goethals–Seidel array gives the Hadamard matrix. The construction and
its proof are in [PALEY_BLOCKS.md](../PALEY_BLOCKS.md).

- Torus element: the square of x ↦ −1/(x + τ) with τ = 3; sign normalisation c = +1.
- T-matrices of order 7 (first rows, from base sequences BS(4, 3)): `+000000`, `0++-000`, `0000+0+`, `00000+0`.
- Row sums of the four sequences: 33, −67, 7, 27 (33² + (−67)² + 7² + 27² = 6356).

## Search

None: the matrix is determined by q and the T-matrices. `python3 search/paley_blocks.py certify 907 7` writes the
certificate.

## Status

**New**, but obtainable classically ([PALEY_BLOCKS.md](../PALEY_BLOCKS.md#status)). Cati–Pasechnik Table 4 lists
1589(3): the smallest known order 2ᵐ · 1589 was 8 · 1589. We found no construction of order 6356 before this one. It
follows from two published theorems that we found combined nowhere: Paley's matrix of order 908 is equivalent to an
array of sixteen circulant blocks of order 227 (Schmidt 1999, after Ito; Yamada 1991), and Seberry and Yamada (1990,
Theorem 4) multiply such an array by T-matrices of order 7.

No other route we know of reaches 6356: Paley I and II fail (4n − 1 = 6355 = 5 · 31 · 41 and 2n − 1 = 3177 = 3² · 353
are not prime powers), and Cooper–Wallis would need Williamson-type matrices or T-matrices of order 227 or 1589,
which are not known to us. Paley's quaternion-type blocks are not amicable ([NOTES.md](../../NOTES.md) §21).

The sources searched are listed in [NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, τ, power, c, the T-matrices, the four sequences and the SHA-256 of the matrix |
| `sequences.txt` | the four sequences of length 1589 |

The matrix itself is not stored; `verify/verify.py` rebuilds it.

## Verification

```
python3 verify/verify.py 6356
```

SHA-256 of `H6356.csv`: `de25e0148c551c2a211c68b1c0b79148cc1f6675d4dc06e707d6f5327d9ea5d9`
