# Hadamard matrix of order 22836

22836 = 4 · 5709, with 5709 = 3 · 11 · 173.

## Construction

Paley's skew-Hadamard matrix of order 692 (q = 691, a prime ≡ 3 mod 8) is a 4 × 4 array of circulant blocks of order
s = 173 once its points are listed along the four orbits of a torus element of PSL(2, 691) and rescaled by signs.
Multiplying it by T-matrices of order 33 (Seberry–Yamada 1990, Theorem 4) gives four circulant ±1 sequences of length
5709 with Σ PAF = 0 off 0, and the Goethals–Seidel array gives the Hadamard matrix. The construction and its proof are
in [PALEY_BLOCKS.md](../PALEY_BLOCKS.md).

- Torus element: x ↦ −1/(x + τ) with τ = 4; sign normalisation c = +1.
- T-matrices of order 33 ([`t-matrices/T33.txt`](../../t-matrices/T33.txt), from a Golay pair of length 32 (Turyn: (1,
  0, …, 0), (0, (A + B)/2), (0, (A − B)/2), 0)).
- Row sums of the four sequences: −71, −79, 105, −23 ((−71)² + (−79)² + 105² + (−23)² = 22836).

## Search

None: the matrix is determined by q and the T-matrices. `python3 search/paley_blocks.py certify 691 33
t-matrices/T33.txt` writes the certificate.

## Status

**First explicit construction** ([PALEY_BLOCKS.md](../PALEY_BLOCKS.md#orders-recorded-only-through-corollary-88i)).
Table 1.53 of the Handbook of Combinatorial Designs (2nd ed., 2007) gives n = 5709 with t = 3, and Đoković
(arXiv:0912.5091) lists it in his set Δ of entries that should have had t = 2. His proof for Δ takes Williamson-type
matrices from Seberry–Yamada's Table A.1 (1992); for 5709 they are of order 173. Table A.1 obtains Williamson-type
matrices of order 173 only from Corollary 8.8(i) (key wl, with Williamson-type matrices of order 43 = (173 − 1)/4),
whose printed proof does not give Williamson-type matrices ([NOTES.md](../../NOTES.md) §7). Table A.1 credits that case
to Miyamoto (J. Combin. Theory Ser. A 57 (1991)), whose paper was not accessible to us; we found no valid proof of it
and no other construction of Williamson-type matrices of order 173.

- Paley's constructions do not apply: 4n − 1 = 22835 = 5 · 4567 and 2n − 1 = 11417 = 7² · 233 are not prime powers.
- Miyamoto and Yamada–Kiyasu do not apply: 5709 is not a prime power, and n − 2 = 5707 = 13 · 439 is not a prime power.
- Cooper–Wallis does not apply otherwise: the other splits of n need Williamson-type matrices or T-matrices that are not
  known to us.
- SageMath's development version (October 17, 2026) does not construct it: `hadamard_matrix(22836, existence=True)`
  returns Unknown.

We found no valid construction of order 22836 before this one. The sources searched are listed in
[NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, τ, power, c, the T-matrices, the four sequences and the SHA-256 of the matrix |
| `sequences.txt` | the four sequences of length 5709 |

The matrix itself is not stored; `verify/verify.py` rebuilds it.

## Verification

```
python3 verify/verify.py 22836
```

SHA-256 of `H22836.csv`: `d91351f71e5b17e78fded84270a9cb00b933d5e714b303645859f4acff1a9c78`
