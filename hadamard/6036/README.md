# Hadamard matrix of order 6036

6036 = 4 · 1509, with 1509 = 3 · 503.

## Construction

Paley's skew-Hadamard matrix of order 2012 (q = 2011, a prime ≡ 3 mod 8) is a 4 × 4 array of circulant blocks of
order s = 503 once its points are listed along the four orbits of a torus element of PSL(2, 2011) and rescaled by
signs. Multiplying it by T-matrices of order 3 (Seberry–Yamada 1990, Theorem 4) gives four circulant ±1 sequences
of length 1509 with Σ PAF = 0 off 0, and the Goethals–Seidel array gives the Hadamard matrix. The construction and
its proof are in [PALEY_BLOCKS.md](../PALEY_BLOCKS.md).

- Torus element: x ↦ −1/(x + τ) with τ = 4; sign normalisation c = +1.
- T-matrices of order 3 (first rows, from base sequences BS(2, 1)): `+00`, `0+0`, `00+`, `000`.
- Row sums of the four sequences: 35, −51, 19, 43 (35² + (−51)² + 19² + 43² = 6036).

## Search

None: the matrix is determined by q and the T-matrices. `python3 search/paley_blocks.py certify 2011 3` writes the
certificate.

## Status

**New, but obtainable classically** ([PALEY_BLOCKS.md](../PALEY_BLOCKS.md#status)). Cati–Pasechnik Table 4 lists
1509(3): the smallest known order 2ᵐ · 1509 was 8 · 1509. We found no construction of order 6036 before this one. It
follows from two published theorems that we found combined nowhere: Paley's matrix of order 2012 is equivalent to an
array of sixteen circulant blocks of order 503 (Schmidt 1999, after Ito; Yamada 1991), and Seberry and Yamada (1990,
Theorem 4) multiply such an array by T-matrices of order 3.

No other route we know of reaches 6036: Paley I and II fail (4n − 1 = 6035 = 5 · 17 · 71 and 2n − 1 = 3017 = 7 · 431
are not prime powers), and Cooper–Wallis would need Williamson-type matrices or T-matrices of order 503 or 1509,
which are not known to us. Paley's quaternion-type blocks are not amicable ([NOTES.md](../../NOTES.md) §21).

The sources searched are listed in [NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, τ, power, c, the T-matrices, the four sequences and the SHA-256 of the matrix |
| `sequences.txt` | the four sequences of length 1509 |

The matrix itself is not stored; `verify/verify.py` rebuilds it.

## Verification

```
python3 verify/verify.py 6036
```

SHA-256 of `H6036.csv`: `064aaf52eed6f4487eeca658aa84d0429415ae8620a17de4c83e4c93914b33bf`
