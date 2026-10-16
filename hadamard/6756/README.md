# Hadamard matrix of order 6756

6756 = 4 · 1689, with 1689 = 3 · 563.

## Construction

Paley's skew-Hadamard matrix of order 2252 (q = 2251, a prime ≡ 3 mod 8) is a 4 × 4 array of circulant blocks of
order s = 563 once its points are listed along the four orbits of a torus element of PSL(2, 2251) and rescaled by
signs. Multiplying it by T-matrices of order 3 (Seberry–Yamada 1990, Theorem 4) gives four circulant ±1 sequences
of length 1689 with Σ PAF = 0 off 0, and the Goethals–Seidel array gives the Hadamard matrix. The construction and
its proof are in [PALEY_BLOCKS.md](../PALEY_BLOCKS.md).

- Torus element: x ↦ −1/(x + τ) with τ = 4; sign normalisation c = −1.
- T-matrices of order 3 (first rows, from base sequences BS(2, 1)): `+00`, `0+0`, `00+`, `000`.
- Row sums of the four sequences: −15, 59, −41, 37 ((−15)² + 59² + (−41)² + 37² = 6756).

## Search

None: the matrix is determined by q and the T-matrices. `python3 search/paley_blocks.py certify 2251 3` writes the
certificate.

## Status

**New**, but obtainable classically ([PALEY_BLOCKS.md](../PALEY_BLOCKS.md#status)). Cati–Pasechnik Table 4 lists
1689(3): the smallest known order 2ᵐ · 1689 was 8 · 1689. We found no construction of order 6756 before this one. It
follows from two published theorems that we found combined nowhere: Paley's matrix of order 2252 is equivalent to an
array of sixteen circulant blocks of order 563 (Schmidt 1999, after Ito; Yamada 1991), and Seberry and Yamada (1990,
Theorem 4) multiply such an array by T-matrices of order 3.

No other route we know of reaches 6756: Paley I and II fail (4n − 1 = 6755 = 5 · 7 · 193 and 2n − 1 = 3377 = 11 · 307
are not prime powers), and Cooper–Wallis would need Williamson-type matrices or T-matrices of order 563 or 1689,
which are not known to us. Paley's quaternion-type blocks are not amicable ([NOTES.md](../../NOTES.md) §21).

The sources searched are listed in [NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, τ, power, c, the T-matrices, the four sequences and the SHA-256 of the matrix |
| `sequences.txt` | the four sequences of length 1689 |

The matrix itself is not stored; `verify/verify.py` rebuilds it.

## Verification

```
python3 verify/verify.py 6756
```

SHA-256 of `H6756.csv`: `a44c19f88c1a50b671030a5407e7b65575c6cf1a40d00600741b30d733721a1f`
