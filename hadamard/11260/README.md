# Hadamard matrix of order 11260

11260 = 4 · 2815, with 2815 = 5 · 563.

## Construction

Paley's skew-Hadamard matrix of order 2252 (q = 2251, a prime ≡ 3 mod 8) is a 4 × 4 array of circulant blocks of
order s = 563 once its points are listed along the four orbits of a torus element of PSL(2, 2251) and rescaled by
signs. Multiplying it by T-matrices of order 5 (Seberry–Yamada 1990, Theorem 4) gives four circulant ±1 sequences
of length 2815 with Σ PAF = 0 off 0, and the Goethals–Seidel array gives the Hadamard matrix. The construction and
its proof are in [PALEY_BLOCKS.md](../PALEY_BLOCKS.md).

- Torus element: x ↦ −1/(x + τ) with τ = 4; sign normalisation c = −1.
- T-matrices of order 5 (first rows, from base sequences BS(3, 2)): `++000`, `00+00`, `000+-`, `00000`.
- Row sums of the four sequences: −35, 75, −63, −21 ((−35)² + 75² + (−63)² + (−21)² = 11260).

## Search

None: the matrix is determined by q and the T-matrices. `python3 search/paley_blocks.py certify 2251 5` writes the
certificate.

## Status

**New**, but obtainable classically ([PALEY_BLOCKS.md](../PALEY_BLOCKS.md#status)). Cati–Pasechnik Table 4 lists
2815(3): the smallest known order 2ᵐ · 2815 was 8 · 2815. We found no construction of order 11260 before this one. It
follows from two published theorems that we found combined nowhere: Paley's matrix of order 2252 is equivalent to an
array of sixteen circulant blocks of order 563 (Schmidt 1999, after Ito; Yamada 1991), and Seberry and Yamada (1990,
Theorem 4) multiply such an array by T-matrices of order 5.

No other route we know of reaches 11260: Paley I and II fail (4n − 1 = 11259 = 3⁴ · 139 and 2n − 1 = 5629 = 13 · 433
are not prime powers), and Cooper–Wallis would need Williamson-type matrices or T-matrices of order 563 or 2815,
which are not known to us. Paley's quaternion-type blocks are not amicable ([NOTES.md](../../NOTES.md) §21).

The sources searched are listed in [NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, τ, power, c, the T-matrices, the four sequences and the SHA-256 of the matrix |
| `sequences.txt` | the four sequences of length 2815 |

The matrix itself is not stored; `verify/verify.py` rebuilds it.

## Verification

```
python3 verify/verify.py 11260
```

SHA-256 of `H11260.csv`: `a2d9ab4a5da3ed797f5dd6edae28a862bee544baeed9d9070ef483a4490ab5f2`
