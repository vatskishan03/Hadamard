# Hadamard matrix of order 7660

7660 = 4 · 1915, with 1915 = 5 · 383.

## Construction

Paley's skew-Hadamard matrix of order 1532 (q = 1531, a prime ≡ 3 mod 8) is a 4 × 4 array of circulant blocks of
order s = 383 once its points are listed along the four orbits of a torus element of PSL(2, 1531) and rescaled by
signs. Multiplying it by T-matrices of order 5 (Seberry–Yamada 1990, Theorem 4) gives four circulant ±1 sequences
of length 1915 with Σ PAF = 0 off 0, and the Goethals–Seidel array gives the Hadamard matrix. The construction and
its proof are in [PALEY_BLOCKS.md](../PALEY_BLOCKS.md).

- Torus element: x ↦ −1/(x + τ) with τ = 4; sign normalisation c = −1.
- T-matrices of order 5 (first rows, from base sequences BS(3, 2)): `++000`, `00+00`, `000+-`, `00000`.
- Row sums of the four sequences: 39, −73, −9, −27 (39² + (−73)² + (−9)² + (−27)² = 7660).

## Search

None: the matrix is determined by q and the T-matrices. `python3 search/paley_blocks.py certify 1531 5` writes the
certificate.

## Status

**New, but obtainable classically** ([PALEY_BLOCKS.md](../PALEY_BLOCKS.md#status)). Cati–Pasechnik Table 4 lists
1915(3): the smallest known order 2ᵐ · 1915 was 8 · 1915. We found no construction of order 7660 before this one. It
follows from two published theorems that we found combined nowhere: Paley's matrix of order 1532 is equivalent to an
array of sixteen circulant blocks of order 383 (Schmidt 1999, after Ito; Yamada 1991), and Seberry and Yamada (1990,
Theorem 4) multiply such an array by T-matrices of order 5.

No other route we know of reaches 7660: Paley I and II fail (4n − 1 = 7659 = 3² · 23 · 37 and 2n − 1 = 3829 = 7 · 547
are not prime powers), and Cooper–Wallis would need Williamson-type matrices or T-matrices of order 383 or 1915,
which are not known to us. Paley's quaternion-type blocks are not amicable ([NOTES.md](../../NOTES.md) §21).

The sources searched are listed in [NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, τ, power, c, the T-matrices, the four sequences and the SHA-256 of the matrix |
| `sequences.txt` | the four sequences of length 1915 |

The matrix itself is not stored; `verify/verify.py` rebuilds it.

## Verification

```
python3 verify/verify.py 7660
```

SHA-256 of `H7660.csv`: `fbc10cd7287ac21d9c79d804f84e94bc63c511098d4b6948ad6698ad464772ea`
