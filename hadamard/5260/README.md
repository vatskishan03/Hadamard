# Hadamard matrix of order 5260

5260 = 4 · 1315, with 1315 = 5 · 263.

## Construction

Paley's skew-Hadamard matrix of order 1052 (q = 1051, a prime ≡ 3 mod 8) is a 4 × 4 array of circulant blocks of
order s = 263 once its points are listed along the four orbits of a torus element of PSL(2, 1051) and rescaled by
signs. Multiplying it by T-matrices of order 5 (Seberry–Yamada 1990, Theorem 4) gives four circulant ±1 sequences
of length 1315 with Σ PAF = 0 off 0, and the Goethals–Seidel array gives the Hadamard matrix. The construction and
its proof are in [PALEY_BLOCKS.md](../PALEY_BLOCKS.md).

- Torus element: x ↦ −1/(x + τ) with τ = 4; sign normalisation c = −1.
- T-matrices of order 5 (first rows, from base sequences BS(3, 2)): `++000`, `00+00`, `000+-`, `00000`.
- Row sums of the four sequences: 15, −25, 21, 63 (15² + (−25)² + 21² + 63² = 5260).

## Search

None: the matrix is determined by q and the T-matrices. `python3 search/paley_blocks.py certify 1051 5` writes the
certificate.

## Status

**New**, but obtainable classically ([PALEY_BLOCKS.md](../PALEY_BLOCKS.md#status)). Cati–Pasechnik Table 4 lists
1315(3): the smallest known order 2ᵐ · 1315 was 8 · 1315. We found no construction of order 5260 before this one. It
follows from two published theorems that we found combined nowhere: Paley's matrix of order 1052 is equivalent to an
array of sixteen circulant blocks of order 263 (Schmidt 1999, after Ito; Yamada 1991), and Seberry and Yamada (1990,
Theorem 4) multiply such an array by T-matrices of order 5.

No other route we know of reaches 5260: Paley I and II fail (4n − 1 = 5259 = 3 · 1753 and 2n − 1 = 2629 = 11 · 239
are not prime powers), and Cooper–Wallis would need Williamson-type matrices or T-matrices of order 263 or 1315,
which are not known to us. Paley's quaternion-type blocks are not amicable ([NOTES.md](../../NOTES.md) §21).

The sources searched are listed in [NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, τ, power, c, the T-matrices, the four sequences and the SHA-256 of the matrix |
| `sequences.txt` | the four sequences of length 1315 |

The matrix itself is not stored; `verify/verify.py` rebuilds it.

## Verification

```
python3 verify/verify.py 5260
```

SHA-256 of `H5260.csv`: `48d441537747d813674ba855c335a9de91a3388df3f739b9b978b2cb63884d3b`
