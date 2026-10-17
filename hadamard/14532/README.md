# Hadamard matrix of order 14532

14532 = 4 · 3633, with 3633 = 3 · 7 · 173.

## Construction

Paley's skew-Hadamard matrix of order 692 (q = 691, a prime ≡ 3 mod 8) is a 4 × 4 array of circulant blocks of order
s = 173 once its points are listed along the four orbits of a torus element of PSL(2, 691) and rescaled by signs.
Multiplying it by T-matrices of order 21 (Seberry–Yamada 1990, Theorem 4) gives four circulant ±1 sequences of length
3633 with Σ PAF = 0 off 0, and the Goethals–Seidel array gives the Hadamard matrix. The construction and its proof are
in [PALEY_BLOCKS.md](../PALEY_BLOCKS.md).

- Torus element: x ↦ −1/(x + τ) with τ = 4; sign normalisation c = +1.
- T-matrices of order 21 ([`t-matrices/T21.txt`](../../t-matrices/T21.txt), from a Golay pair of length 20 (Turyn: (1,
  0, …, 0), (0, (A + B)/2), (0, (A − B)/2), 0)).
- Row sums of the four sequences: −53, −33, 103, −5 ((−53)² + (−33)² + 103² + (−5)² = 14532).

## Search

None: the matrix is determined by q and the T-matrices. `python3 search/paley_blocks.py certify 691 21
t-matrices/T21.txt` writes the certificate.

## Status

**First explicit construction** ([PALEY_BLOCKS.md](../PALEY_BLOCKS.md#orders-recorded-only-through-corollary-88i)).
Table 1.53 of the Handbook of Combinatorial Designs (2nd ed., 2007) gives n = 3633 with t = 2, as known. The only
construction we traced uses Williamson-type matrices of order 173 from Seberry–Yamada's Table A.1 (1992). Table A.1
obtains Williamson-type matrices of order 173 only from Corollary 8.8(i) (key wl, with Williamson-type matrices of order
43 = (173 − 1)/4), whose printed proof does not give Williamson-type matrices ([NOTES.md](../../NOTES.md) §7). Table A.1
credits that case to Miyamoto (J. Combin. Theory Ser. A 57 (1991)), whose paper was not accessible to us; we found no
valid proof of it and no other construction of Williamson-type matrices of order 173.

- Paley's constructions do not apply: 4n − 1 = 14531 = 11 · 1321 and 2n − 1 = 7265 = 5 · 1453 are not prime powers.
- Miyamoto and Yamada–Kiyasu do not apply: 3633 is not a prime power, and n − 2 = 3631 is ≡ 7 (mod 8).
- Cooper–Wallis does not apply otherwise: the other splits of n need Williamson-type matrices or T-matrices that are not
  known to us.
- SageMath's development version (October 17, 2026) does not construct it: `hadamard_matrix(14532, existence=True)`
  returns Unknown.

We found no valid construction of order 14532 before this one. The sources searched are listed in
[NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, τ, power, c, the T-matrices, the four sequences and the SHA-256 of the matrix |
| `sequences.txt` | the four sequences of length 3633 |

The matrix itself is not stored; `verify/verify.py` rebuilds it.

## Verification

```
python3 verify/verify.py 14532
```

SHA-256 of `H14532.csv`: `bfbc3b8841286e9fdd9c09511419bccccd6a099dc1f9fd89d086ac0df1dea90c`
