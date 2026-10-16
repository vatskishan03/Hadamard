# Hadamard matrix of order 10060

10060 = 4 · 2515, with 2515 = 5 · 503.

## Construction

Paley's skew-Hadamard matrix of order 2012 (q = 2011, a prime ≡ 3 mod 8) is a 4 × 4 array of circulant blocks of
order s = 503 once its points are listed along the four orbits of a torus element of PSL(2, 2011) and rescaled by
signs. Multiplying it by T-matrices of order 5 (Seberry–Yamada 1990, Theorem 4) gives four circulant ±1 sequences
of length 2515 with Σ PAF = 0 off 0, and the Goethals–Seidel array gives the Hadamard matrix. The construction and
its proof are in [PALEY_BLOCKS.md](../PALEY_BLOCKS.md).

- Torus element: x ↦ −1/(x + τ) with τ = 4; sign normalisation c = +1.
- T-matrices of order 5 (first rows, from base sequences BS(3, 2)): `++000`, `00+00`, `000+-`, `00000`.
- Row sums of the four sequences: 45, −85, 27, −9 (45² + (−85)² + 27² + (−9)² = 10060).

## Search

None: the matrix is determined by q and the T-matrices. `python3 search/paley_blocks.py certify 2011 5` writes the
certificate.

## Status

**New**, but obtainable classically ([PALEY_BLOCKS.md](../PALEY_BLOCKS.md#status)). Cati–Pasechnik Table 4 lists
2515(3): the smallest known order 2ᵐ · 2515 was 8 · 2515. We found no construction of order 10060 before this one. It
follows from two published theorems that we found combined nowhere: Paley's matrix of order 2012 is equivalent to an
array of sixteen circulant blocks of order 503 (Schmidt 1999, after Ito; Yamada 1991), and Seberry and Yamada (1990,
Theorem 4) multiply such an array by T-matrices of order 5.

No other route we know of reaches 10060: Paley I and II fail (4n − 1 = 10059 = 3 · 7 · 479 and 2n − 1 = 5029 = 47 ·
107 are not prime powers), and Cooper–Wallis would need Williamson-type matrices or T-matrices of order 503 or 2515,
which are not known to us. Paley's quaternion-type blocks are not amicable ([NOTES.md](../../NOTES.md) §21).

The sources searched are listed in [NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, τ, power, c, the T-matrices, the four sequences and the SHA-256 of the matrix |
| `sequences.txt` | the four sequences of length 2515 |

The matrix itself is not stored; `verify/verify.py` rebuilds it.

## Verification

```
python3 verify/verify.py 10060
```

SHA-256 of `H10060.csv`: `e8372df9d7434993d39d3881c0a5eae788ae1d8f2ef47fcc845a296c4e2538fb`
