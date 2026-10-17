# Hadamard matrix of order 28372

28372 = 4 · 7093, with 7093 = 41 · 173.

## Construction

Paley's skew-Hadamard matrix of order 692 (q = 691, a prime ≡ 3 mod 8) is a 4 × 4 array of circulant blocks of order
s = 173 once its points are listed along the four orbits of a torus element of PSL(2, 691) and rescaled by signs.
Multiplying it by T-matrices of order 41 (Seberry–Yamada 1990, Theorem 4) gives four circulant ±1 sequences of length
7093 with Σ PAF = 0 off 0, and the Goethals–Seidel array gives the Hadamard matrix. The construction and its proof are
in [PALEY_BLOCKS.md](../PALEY_BLOCKS.md).

- Torus element: x ↦ −1/(x + τ) with τ = 4; sign normalisation c = +1.
- T-matrices of order 41 ([`t-matrices/T41.txt`](../../t-matrices/T41.txt), sources in
  [t-matrices/README.md](../../t-matrices/README.md)).
- Row sums of the four sequences: −49, −151, −31, 47 ((−49)² + (−151)² + (−31)² + 47² = 28372).

## Search

None: the matrix is determined by q and the T-matrices. `python3 search/paley_blocks.py certify 691 41
t-matrices/T41.txt` writes the certificate.

## Status

**First explicit construction** ([PALEY_BLOCKS.md](../PALEY_BLOCKS.md#orders-recorded-only-through-corollary-88i)).
Table 1.53 of the Handbook of Combinatorial Designs (2nd ed., 2007) gives n = 7093 with t = 3, and Đoković
(arXiv:0912.5091) lists it in his set Δ of entries that should have had t = 2. His proof for Δ takes Williamson-type
matrices from Seberry–Yamada's Table A.1 (1992); for 7093 they are of order 173. Table A.1 obtains Williamson-type
matrices of order 173 only from Corollary 8.8(i) (key wl, with Williamson-type matrices of order 43 = (173 − 1)/4),
whose printed proof does not give Williamson-type matrices ([NOTES.md](../../NOTES.md) §7). Table A.1 credits that case
to Miyamoto (J. Combin. Theory Ser. A 57 (1991)), whose paper was not accessible to us; we found no valid proof of it
and no other construction of Williamson-type matrices of order 173.

- Paley's constructions do not apply: 4n − 1 = 28371 = 3 · 7² · 193 and 2n − 1 = 14185 = 5 · 2837 are not prime powers.
- Miyamoto and Yamada–Kiyasu do not apply: 7093 is not a prime power, and n − 2 = 7091 = 7 · 1013 is not a prime power.
- Cooper–Wallis does not apply otherwise: the other splits of n need Williamson-type matrices or T-matrices that are not
  known to us.
- SageMath's development version (October 17, 2026) does not construct it: `hadamard_matrix(28372, existence=True)`
  returns Unknown.

We found no valid construction of order 28372 before this one. The sources searched are listed in
[NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, τ, power, c, the T-matrices, the four sequences and the SHA-256 of the matrix |
| `sequences.txt` | the four sequences of length 7093 |

The matrix itself is not stored; `verify/verify.py` rebuilds it.

## Verification

```
python3 verify/verify.py 28372
```

SHA-256 of `H28372.csv`: `f4d14020f51bce0c0500f7bff092e5c0142f95c2dc04a2bcac309319a4aae5ed`
