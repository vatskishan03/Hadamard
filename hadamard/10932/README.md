# Hadamard matrix of order 10932

10932 = 4 · 2733, with 2733 = 3 · 911.

## Construction

Paley's skew-Hadamard matrix of order 3644 (q = 3643, a prime ≡ 3 mod 8) is a 4 × 4 array of circulant blocks of
order s = 911 once its points are listed along the four orbits of a torus element of PSL(2, 3643) and rescaled by
signs. Multiplying it by T-matrices of order 3 (Seberry–Yamada 1990, Theorem 4) gives four circulant ±1 sequences
of length 2733 with Σ PAF = 0 off 0, and the Goethals–Seidel array gives the Hadamard matrix. The construction and
its proof are in [PALEY_BLOCKS.md](../PALEY_BLOCKS.md).

- Torus element: the square of x ↦ −1/(x + τ) with τ = 3; sign normalisation c = +1.
- T-matrices of order 3 (first rows, from base sequences BS(2, 1)): `+00`, `0+0`, `00+`, `000`.
- Row sums of the four sequences: −49, −17, 69, −59 ((−49)² + (−17)² + 69² + (−59)² = 10932).

## Search

None: the matrix is determined by q and the T-matrices. `python3 search/paley_blocks.py certify 3643 3` writes the
certificate.

## Status

**New**, but obtainable classically ([PALEY_BLOCKS.md](../PALEY_BLOCKS.md#status)). Cati–Pasechnik Table 4 lists
2733(3): the smallest known order 2ᵐ · 2733 was 8 · 2733. We found no construction of order 10932 before this one. It
follows from two published theorems that we found combined nowhere: Paley's matrix of order 3644 is equivalent to an
array of sixteen circulant blocks of order 911 (Schmidt 1999, after Ito; Yamada 1991), and Seberry and Yamada (1990,
Theorem 4) multiply such an array by T-matrices of order 3.

No other route we know of reaches 10932: Paley I and II fail (4n − 1 = 10931 = 17 · 643 and 2n − 1 = 5465 = 5 · 1093
are not prime powers), and Cooper–Wallis would need Williamson-type matrices or T-matrices of order 911 or 2733,
which are not known to us. Paley's quaternion-type blocks are not amicable ([NOTES.md](../../NOTES.md) §21).

The sources searched are listed in [NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, τ, power, c, the T-matrices, the four sequences and the SHA-256 of the matrix |
| `sequences.txt` | the four sequences of length 2733 |

The matrix itself is not stored; `verify/verify.py` rebuilds it.

## Verification

```
python3 verify/verify.py 10932
```

SHA-256 of `H10932.csv`: `633a6c813cdb5afcd9be3a5ec57e7bc82daf0de5f840fa9833aa6e57e3fe98c8`
