# Hadamard matrix of order 3460

3460 = 4 · 865, with 865 = 5 · 173.

## Construction

Paley's skew-Hadamard matrix of order 692 (q = 691, a prime ≡ 3 mod 8) is a 4 × 4 array of circulant blocks of order
s = 173 once its points are listed along the four orbits of a torus element of PSL(2, 691) and rescaled by signs.
Multiplying it by T-matrices of order 5 (Seberry–Yamada 1990, Theorem 4) gives four circulant ±1 sequences of length 865
with Σ PAF = 0 off 0, and the Goethals–Seidel array gives the Hadamard matrix. The construction and its proof are in
[PALEY_BLOCKS.md](../PALEY_BLOCKS.md).

- Torus element: x ↦ −1/(x + τ) with τ = 4; sign normalisation c = +1.
- T-matrices of order 5 (first rows, from base sequences BS(3, 2)): `++000`, `00+00`, `000+-`, `00000`.
- Row sums of the four sequences: −7, 19, 41, −37 ((−7)² + 19² + 41² + (−37)² = 3460).

## Search

None: the matrix is determined by q and the T-matrices. `python3 search/paley_blocks.py certify 691 5` writes the
certificate.

## Status

**First explicit construction** ([PALEY_BLOCKS.md](../PALEY_BLOCKS.md#orders-recorded-only-through-corollary-88i)).
Cati–Pasechnik's Table 4 omits n = 865: it inherits m = 2 from Seberry–Yamada's 2020 Table A.17, and their 1992 Table
A.2 builds the order (code o2) from Williamson-type matrices of order 173. Table A.1 obtains Williamson-type matrices of
order 173 only from Corollary 8.8(i) (key wl, with Williamson-type matrices of order 43 = (173 − 1)/4), whose printed
proof does not give Williamson-type matrices ([NOTES.md](../../NOTES.md) §7). Table A.1 credits that case to Miyamoto
(J. Combin. Theory Ser. A 57 (1991)), whose paper was not accessible to us; we found no valid proof of it and no other
construction of Williamson-type matrices of order 173.

- Paley's constructions do not apply: 4n − 1 = 3459 = 3 · 1153 and 2n − 1 = 1729 = 7 · 13 · 19 are not prime powers.
- Miyamoto and Yamada–Kiyasu do not apply: 865 is not a prime power, and n − 2 = 863 is ≡ 7 (mod 8).
- Cooper–Wallis does not apply otherwise: the other splits of n need Williamson-type matrices or T-matrices that are not
  known to us.
- SageMath's development version (October 17, 2026) does not construct it: `hadamard_matrix(3460, existence=True)`
  returns Unknown.

We found no valid construction of order 3460 before this one. The sources searched are listed in
[NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, τ, power, c, the T-matrices, the four sequences and the SHA-256 of the matrix |
| `sequences.txt` | the four sequences of length 865 |

The matrix itself is not stored; `verify/verify.py` rebuilds it.

## Verification

```
python3 verify/verify.py 3460
```

SHA-256 of `H3460.csv`: `06aa16c919a8b23b993c19d579104018c65500bca13e475947013a7c4f35f7ae`
