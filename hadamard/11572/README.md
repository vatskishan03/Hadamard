# Hadamard matrix of order 11572

11572 = 4 · 2893, with 2893 = 11 · 263.

## Construction

Paley's skew-Hadamard matrix of order 1052 (q = 1051, a prime ≡ 3 mod 8) is a 4 × 4 array of circulant blocks of
order s = 263 once its points are listed along the four orbits of a torus element of PSL(2, 1051) and rescaled by
signs. Multiplying it by T-matrices of order 11 (Seberry–Yamada 1990, Theorem 4) gives four circulant ±1 sequences
of length 2893 with Σ PAF = 0 off 0, and the Goethals–Seidel array gives the Hadamard matrix. The construction and
its proof are in [PALEY_BLOCKS.md](../PALEY_BLOCKS.md).

- Torus element: x ↦ −1/(x + τ) with τ = 4; sign normalisation c = −1.
- T-matrices of order 11 (first rows, from base sequences BS(6, 5)): `++00+000000`, `00++0-00000`, `000000++-+-`,
  `00000000000`.
- Row sums of the four sequences: −5, −17, 43, 97 ((−5)² + (−17)² + 43² + 97² = 11572).

## Search

None: the matrix is determined by q and the T-matrices. `python3 search/paley_blocks.py certify 1051 11` writes the
certificate.

## Status

**New, but obtainable classically** ([PALEY_BLOCKS.md](../PALEY_BLOCKS.md#status)). Cati–Pasechnik Table 4 lists
2893(3): the smallest known order 2ᵐ · 2893 was 8 · 2893. We found no construction of order 11572 before this one. It
follows from two published theorems that we found combined nowhere: Paley's matrix of order 1052 is equivalent to an
array of sixteen circulant blocks of order 263 (Schmidt 1999, after Ito; Yamada 1991), and Seberry and Yamada (1990,
Theorem 4) multiply such an array by T-matrices of order 11.

No other route we know of reaches 11572: Paley I and II fail (4n − 1 = 11571 = 3 · 7 · 19 · 29 and 2n − 1 = 5785 = 5
· 13 · 89 are not prime powers), and Cooper–Wallis would need Williamson-type matrices or T-matrices of order 263 or
2893, which are not known to us. Paley's quaternion-type blocks are not amicable ([NOTES.md](../../NOTES.md) §21).

The sources searched are listed in [NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, τ, power, c, the T-matrices, the four sequences and the SHA-256 of the matrix |
| `sequences.txt` | the four sequences of length 2893 |

The matrix itself is not stored; `verify/verify.py` rebuilds it.

## Verification

```
python3 verify/verify.py 11572
```

SHA-256 of `H11572.csv`: `88de65af512ebd085f128f7a2b16184f8c4a58b69535e599a14d2fafcddc51a3`
