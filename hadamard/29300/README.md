# Hadamard matrix of order 29300

29300 = 4 · 7325, with 7325 = 5² · 293.

## Construction

Paley's skew-Hadamard matrix of order 1172 (q = 1171, a prime ≡ 3 mod 8) is a 4 × 4 array of circulant blocks of order
s = 293 once its points are listed along the four orbits of a torus element of PSL(2, 1171) and rescaled by signs.
Multiplying it by T-matrices of order 25 (Seberry–Yamada 1990, Theorem 4) gives four circulant ±1 sequences of length
7325 with Σ PAF = 0 off 0, and the Goethals–Seidel array gives the Hadamard matrix. The construction and its proof are
in [PALEY_BLOCKS.md](../PALEY_BLOCKS.md).

- Torus element: x ↦ −1/(x + τ) with τ = 4; sign normalisation c = −1.
- T-matrices of order 25 ([`t-matrices/T25.txt`](../../t-matrices/T25.txt), sources in
  [t-matrices/README.md](../../t-matrices/README.md)).
- Row sums of the four sequences: −47, 71, −21, −147 ((−47)² + 71² + (−21)² + (−147)² = 29300).

## Search

None: the matrix is determined by q and the T-matrices. `python3 search/paley_blocks.py certify 1171 25
t-matrices/T25.txt` writes the certificate.

## Status

**First explicit construction** ([PALEY_BLOCKS.md](../PALEY_BLOCKS.md#orders-recorded-only-through-corollary-88i)).
Table 1.53 of the Handbook of Combinatorial Designs (2nd ed., 2007) gives n = 7325 with t = 3, and Đoković
(arXiv:0912.5091) lists it in his set Δ of entries that should have had t = 2. His proof for Δ takes Williamson-type
matrices from Seberry–Yamada's Table A.1 (1992); for 7325 they are of order 293. Table A.1 obtains Williamson-type
matrices of order 293 only from Corollary 8.8(i) (key wl, with Williamson-type matrices of order 73 = (293 − 1)/4),
whose printed proof does not give Williamson-type matrices ([NOTES.md](../../NOTES.md) §7). Table A.1 credits that case
to Miyamoto (J. Combin. Theory Ser. A 57 (1991)), whose paper was not accessible to us; we found no valid proof of it
and no other construction of Williamson-type matrices of order 293.

- Paley's constructions do not apply: 4n − 1 = 29299 = 83 · 353 and 2n − 1 = 14649 = 3 · 19 · 257 are not prime powers.
- Miyamoto and Yamada–Kiyasu do not apply: 7325 is not a prime power, and n − 2 = 7323 = 3 · 2441 is not a prime power.
- Cooper–Wallis does not apply otherwise: the other splits of n need Williamson-type matrices or T-matrices that are not
  known to us.
- SageMath's development version (October 17, 2026) does not construct it: `hadamard_matrix(29300, existence=True)`
  returns Unknown.

We found no valid construction of order 29300 before this one. The sources searched are listed in
[NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, τ, power, c, the T-matrices, the four sequences and the SHA-256 of the matrix |
| `sequences.txt` | the four sequences of length 7325 |

The matrix itself is not stored; `verify/verify.py` rebuilds it.

## Verification

```
python3 verify/verify.py 29300
```

SHA-256 of `H29300.csv`: `fcc62dba526e15358a4870b62c5fc4c9f7b7f85889b9d6829f144d93e66ea9cf`
