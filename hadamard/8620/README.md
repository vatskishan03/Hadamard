# Hadamard matrix of order 8620

8620 = 4 · 2155, with 2155 = 5 · 431.

## Construction

Paley's skew-Hadamard matrix of order 1724 (q = 1723, a prime ≡ 3 mod 8) is a 4 × 4 array of circulant blocks of
order s = 431 once its points are listed along the four orbits of a torus element of PSL(2, 1723) and rescaled by
signs. Multiplying it by T-matrices of order 5 (Seberry–Yamada 1990, Theorem 4) gives four circulant ±1 sequences
of length 2155 with Σ PAF = 0 off 0, and the Goethals–Seidel array gives the Hadamard matrix. The construction and
its proof are in [PALEY_BLOCKS.md](../PALEY_BLOCKS.md).

- Torus element: the square of x ↦ −1/(x + τ) with τ = 3; sign normalisation c = +1.
- T-matrices of order 5 (first rows, from base sequences BS(3, 2)): `++000`, `00+00`, `000+-`, `00000`.
- Row sums of the four sequences: 23, −41, 37, 71 (23² + (−41)² + 37² + 71² = 8620).

## Search

None: the matrix is determined by q and the T-matrices. `python3 search/paley_blocks.py certify 1723 5` writes the
certificate.

## Status

**New, but obtainable classically** ([PALEY_BLOCKS.md](../PALEY_BLOCKS.md#status)). Cati–Pasechnik Table 4 lists
2155(3): the smallest known order 2ᵐ · 2155 was 8 · 2155. We found no construction of order 8620 before this one. It
follows from two published theorems that we found combined nowhere: Paley's matrix of order 1724 is equivalent to an
array of sixteen circulant blocks of order 431 (Schmidt 1999, after Ito; Yamada 1991), and Seberry and Yamada (1990,
Theorem 4) multiply such an array by T-matrices of order 5.

No other route we know of reaches 8620: Paley I and II fail (4n − 1 = 8619 = 3 · 13² · 17 and 2n − 1 = 4309 = 31 ·
139 are not prime powers), and Cooper–Wallis would need Williamson-type matrices or T-matrices of order 431 or 2155,
which are not known to us. Paley's quaternion-type blocks are not amicable ([NOTES.md](../../NOTES.md) §21).

The sources searched are listed in [NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, τ, power, c, the T-matrices, the four sequences and the SHA-256 of the matrix |
| `sequences.txt` | the four sequences of length 2155 |

The matrix itself is not stored; `verify/verify.py` rebuilds it.

## Verification

```
python3 verify/verify.py 8620
```

SHA-256 of `H8620.csv`: `70c2b88676ec7e387ae55445ae578bfe26540a46c03344e04a8cd9cc3930b33f`
