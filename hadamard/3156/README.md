# Hadamard matrix of order 3156

3156 = 4 · 789, with 789 = 3 · 263.

## Construction

Paley's skew-Hadamard matrix of order 1052 (q = 1051, a prime ≡ 3 mod 8) is a 4 × 4 array of circulant blocks of
order s = 263 once its points are listed along the four orbits of a torus element of PSL(2, 1051) and rescaled by
signs. Multiplying it by T-matrices of order 3 (Seberry–Yamada 1990, Theorem 4) gives four circulant ±1 sequences
of length 789 with Σ PAF = 0 off 0, and the Goethals–Seidel array gives the Hadamard matrix. The construction and
its proof are in [PALEY_BLOCKS.md](../PALEY_BLOCKS.md).

- Torus element: x ↦ −1/(x + τ) with τ = 4; sign normalisation c = −1.
- T-matrices of order 3 (first rows, from base sequences BS(2, 1)): `+00`, `0+0`, `00+`, `000`.
- Row sums of the four sequences: −7, 9, 1, 55 ((−7)² + 9² + 1² + 55² = 3156).

## Search

None: the matrix is determined by q and the T-matrices. `python3 search/paley_blocks.py certify 1051 3` writes the
certificate.

## Status

**New, but obtainable classically** ([PALEY_BLOCKS.md](../PALEY_BLOCKS.md#status)). Cati–Pasechnik Table 4 lists
789(3): the smallest known order 2ᵐ · 789 was 8 · 789. We found no construction of order 3156 before this one. It
follows from two published theorems that we found combined nowhere: Paley's matrix of order 1052 is equivalent to an
array of sixteen circulant blocks of order 263 (Schmidt 1999, after Ito; Yamada 1991), and Seberry and Yamada (1990,
Theorem 4) multiply such an array by T-matrices of order 3.

No other route we know of reaches 3156: Paley I and II fail (4n − 1 = 3155 = 5 · 631 and 2n − 1 = 1577 = 19 · 83 are
not prime powers), and Cooper–Wallis would need Williamson-type matrices or T-matrices of order 263 or 789, which are
not known to us. Paley's quaternion-type blocks are not amicable ([NOTES.md](../../NOTES.md) §21).

The sources searched are listed in [NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, τ, power, c, the T-matrices, the four sequences and the SHA-256 of the matrix |
| `sequences.txt` | the four sequences of length 789 |

The matrix itself is not stored; `verify/verify.py` rebuilds it.

## Verification

```
python3 verify/verify.py 3156
```

SHA-256 of `H3156.csv`: `c318b8b6d037347dc5ea01f008141833ae12a1719e10a80d563764b872147834`
