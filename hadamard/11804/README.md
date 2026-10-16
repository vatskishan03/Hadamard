# Hadamard matrix of order 11804

11804 = 4 · 2951, with 2951 = 13 · 227.

## Construction

Paley's skew-Hadamard matrix of order 908 (q = 907, a prime ≡ 3 mod 8) is a 4 × 4 array of circulant blocks of
order s = 227 once its points are listed along the four orbits of a torus element of PSL(2, 907) and rescaled by
signs. Multiplying it by T-matrices of order 13 (Seberry–Yamada 1990, Theorem 4) gives four circulant ±1 sequences
of length 2951 with Σ PAF = 0 off 0, and the Goethals–Seidel array gives the Hadamard matrix. The construction and
its proof are in [PALEY_BLOCKS.md](../PALEY_BLOCKS.md).

- Torus element: the square of x ↦ −1/(x + τ) with τ = 3; sign normalisation c = +1.
- T-matrices of order 13 (first rows, from base sequences BS(7, 6)): `+0+0+00000000`, `0+0+0-+000000`, `0000000+0000-`, `00000000++--0`.
- Row sums of the four sequences: 45, −61, 57, −53 (45² + (−61)² + 57² + (−53)² = 11804).

## Search

None: the matrix is determined by q and the T-matrices. `python3 search/paley_blocks.py certify 907 13` writes the
certificate.

## Status

**New**, but obtainable classically ([PALEY_BLOCKS.md](../PALEY_BLOCKS.md#status)). Cati–Pasechnik Table 4 lists
2951(3): the smallest known order 2ᵐ · 2951 was 8 · 2951. We found no construction of order 11804 before this one. It
follows from two published theorems that we found combined nowhere: Paley's matrix of order 908 is equivalent to an
array of sixteen circulant blocks of order 227 (Schmidt 1999, after Ito; Yamada 1991), and Seberry and Yamada (1990,
Theorem 4) multiply such an array by T-matrices of order 13.

No other route we know of reaches 11804: Paley I and II fail (4n − 1 = 11803 = 11 · 29 · 37 and 2n − 1 = 5901 = 3 · 7
· 281 are not prime powers), and Cooper–Wallis would need Williamson-type matrices or T-matrices of order 227 or
2951, which are not known to us. Paley's quaternion-type blocks are not amicable ([NOTES.md](../../NOTES.md) §21).

The sources searched are listed in [NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, τ, power, c, the T-matrices, the four sequences and the SHA-256 of the matrix |
| `sequences.txt` | the four sequences of length 2951 |

The matrix itself is not stored; `verify/verify.py` rebuilds it.

## Verification

```
python3 verify/verify.py 11804
```

SHA-256 of `H11804.csv`: `45fcff69e8d7c4512782baac3c1daea80ae7a71b37038e1c3526f0dee323e8e7`
