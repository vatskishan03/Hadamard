# Hadamard matrix of order 10644

10644 = 4 · 2661, with 2661 = 3 · 887.

## Construction

Paley's skew-Hadamard matrix of order 3548 (q = 3547, a prime ≡ 3 mod 8) is a 4 × 4 array of circulant blocks of
order s = 887 once its points are listed along the four orbits of a torus element of PSL(2, 3547) and rescaled by
signs. Multiplying it by T-matrices of order 3 (Seberry–Yamada 1990, Theorem 4) gives four circulant ±1 sequences
of length 2661 with Σ PAF = 0 off 0, and the Goethals–Seidel array gives the Hadamard matrix. The construction and
its proof are in [PALEY_BLOCKS.md](../PALEY_BLOCKS.md).

- Torus element: the square of x ↦ −1/(x + τ) with τ = 3; sign normalisation c = +1.
- T-matrices of order 3 (first rows, from base sequences BS(2, 1)): `+00`, `0+0`, `00+`, `000`.
- Row sums of the four sequences: 1, 71, −69, 29 (1² + 71² + (−69)² + 29² = 10644).

## Search

None: the matrix is determined by q and the T-matrices. `python3 search/paley_blocks.py certify 3547 3` writes the
certificate.

## Status

**New, but obtainable classically** ([PALEY_BLOCKS.md](../PALEY_BLOCKS.md#status)). Cati–Pasechnik Table 4 lists
2661(3): the smallest known order 2ᵐ · 2661 was 8 · 2661. We found no construction of order 10644 before this one. It
follows from two published theorems that we found combined nowhere: Paley's matrix of order 3548 is equivalent to an
array of sixteen circulant blocks of order 887 (Schmidt 1999, after Ito; Yamada 1991), and Seberry and Yamada (1990,
Theorem 4) multiply such an array by T-matrices of order 3.

No other route we know of reaches 10644: Paley I and II fail (4n − 1 = 10643 = 29 · 367 and 2n − 1 = 5321 = 17 · 313
are not prime powers), and Cooper–Wallis would need Williamson-type matrices or T-matrices of order 887 or 2661,
which are not known to us. Paley's quaternion-type blocks are not amicable ([NOTES.md](../../NOTES.md) §21).

The sources searched are listed in [NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, τ, power, c, the T-matrices, the four sequences and the SHA-256 of the matrix |
| `sequences.txt` | the four sequences of length 2661 |

The matrix itself is not stored; `verify/verify.py` rebuilds it.

## Verification

```
python3 verify/verify.py 10644
```

SHA-256 of `H10644.csv`: `241afedc477ab1823dbbb44d35ab21d8e60046560e052b7ed8a259c908a4ddc8`
