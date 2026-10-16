# Hadamard matrix of order 7364

7364 = 4 · 1841, with 1841 = 7 · 263.

## Construction

Paley's skew-Hadamard matrix of order 1052 (q = 1051, a prime ≡ 3 mod 8) is a 4 × 4 array of circulant blocks of
order s = 263 once its points are listed along the four orbits of a torus element of PSL(2, 1051) and rescaled by
signs. Multiplying it by T-matrices of order 7 (Seberry–Yamada 1990, Theorem 4) gives four circulant ±1 sequences
of length 1841 with Σ PAF = 0 off 0, and the Goethals–Seidel array gives the Hadamard matrix. The construction and
its proof are in [PALEY_BLOCKS.md](../PALEY_BLOCKS.md).

- Torus element: x ↦ −1/(x + τ) with τ = 4; sign normalisation c = −1.
- T-matrices of order 7 (first rows, from base sequences BS(4, 3)): `+000000`, `0++-000`, `0000+0+`, `00000+0`.
- Row sums of the four sequences: −49, 9, −11, 69 ((−49)² + 9² + (−11)² + 69² = 7364).

## Search

None: the matrix is determined by q and the T-matrices. `python3 search/paley_blocks.py certify 1051 7` writes the
certificate.

## Status

**New**, but obtainable classically ([PALEY_BLOCKS.md](../PALEY_BLOCKS.md#status)). Cati–Pasechnik Table 4 lists
1841(3): the smallest known order 2ᵐ · 1841 was 8 · 1841. We found no construction of order 7364 before this one. It
follows from two published theorems that we found combined nowhere: Paley's matrix of order 1052 is equivalent to an
array of sixteen circulant blocks of order 263 (Schmidt 1999, after Ito; Yamada 1991), and Seberry and Yamada (1990,
Theorem 4) multiply such an array by T-matrices of order 7.

No other route we know of reaches 7364: Paley I and II fail (4n − 1 = 7363 = 37 · 199 and 2n − 1 = 3681 = 3² · 409
are not prime powers), and Cooper–Wallis would need Williamson-type matrices or T-matrices of order 263 or 1841,
which are not known to us. Paley's quaternion-type blocks are not amicable ([NOTES.md](../../NOTES.md) §21).

The sources searched are listed in [NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, τ, power, c, the T-matrices, the four sequences and the SHA-256 of the matrix |
| `sequences.txt` | the four sequences of length 1841 |

The matrix itself is not stored; `verify/verify.py` rebuilds it.

## Verification

```
python3 verify/verify.py 7364
```

SHA-256 of `H7364.csv`: `8d7fc4e418f63672f0000483e45d7a66c999b961c8ce7306353b9ab1c6a49144`
