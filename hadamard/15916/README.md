# Hadamard matrix of order 15916

15916 = 4 · 3979, with 3979 = 23 · 173.

## Construction

Paley's skew-Hadamard matrix of order 692 (q = 691, a prime ≡ 3 mod 8) is a 4 × 4 array of circulant blocks of order
s = 173 once its points are listed along the four orbits of a torus element of PSL(2, 691) and rescaled by signs.
Multiplying it by T-matrices of order 23 (Seberry–Yamada 1990, Theorem 4) gives four circulant ±1 sequences of length
3979 with Σ PAF = 0 off 0, and the Goethals–Seidel array gives the Hadamard matrix. The construction and its proof are
in [PALEY_BLOCKS.md](../PALEY_BLOCKS.md).

- Torus element: x ↦ −1/(x + τ) with τ = 4; sign normalisation c = +1.
- T-matrices of order 23 ([`t-matrices/T23.txt`](../../t-matrices/T23.txt), sources in
  [t-matrices/README.md](../../t-matrices/README.md)).
- Row sums of the four sequences: 25, −61, 107, −11 (25² + (−61)² + 107² + (−11)² = 15916).

## Search

None: the matrix is determined by q and the T-matrices. `python3 search/paley_blocks.py certify 691 23
t-matrices/T23.txt` writes the certificate.

## Status

**First explicit construction** ([PALEY_BLOCKS.md](../PALEY_BLOCKS.md#orders-recorded-only-through-corollary-88i)).
Table 1.53 of the Handbook of Combinatorial Designs (2nd ed., 2007) gives n = 3979 with t = 2, as known. The only
construction we traced uses Williamson-type matrices of order 173 from Seberry–Yamada's Table A.1 (1992). Table A.1
obtains Williamson-type matrices of order 173 only from Corollary 8.8(i) (key wl, with Williamson-type matrices of order
43 = (173 − 1)/4), whose printed proof does not give Williamson-type matrices ([NOTES.md](../../NOTES.md) §7). Table A.1
credits that case to Miyamoto (J. Combin. Theory Ser. A 57 (1991)), whose paper was not accessible to us; we found no
valid proof of it. Since 6 October 2026 there is another construction: Verma's symmetric conference matrix of order 86
(arXiv:2610.08893, checked exactly here) gives symmetric Williamson-type matrices of order 173 by Corollary 8.8(ii), and
Cooper–Wallis then gives this order; we found no statement of this.

- Paley's constructions do not apply: 4n − 1 = 15915 = 3 · 5 · 1061 and 2n − 1 = 7957 = 73 · 109 are not prime powers.
- Miyamoto's theorem does not apply: it needs a symmetric conference matrix of order n + 1, which cannot exist for n ≡ 3
  (mod 4); Yamada–Kiyasu does not apply: n − 2 = 3977 = 41 · 97 is not a prime power.
- Cooper–Wallis does not apply otherwise: the other splits of n need Williamson-type matrices or T-matrices that are not
  known to us.
- SageMath's development version (October 17, 2026) does not construct it: `hadamard_matrix(15916, existence=True)`
  returns Unknown.

We found no valid construction of order 15916 before this one. The sources searched are listed in
[NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, τ, power, c, the T-matrices, the four sequences and the SHA-256 of the matrix |
| `sequences.txt` | the four sequences of length 3979 |

The matrix itself is not stored; `verify/verify.py` rebuilds it.

## Verification

```
python3 verify/verify.py 15916
```

SHA-256 of `H15916.csv`: `c608719bf3127d88d5b8ae16babf7a4e45008c794fe4a2a08894605572dd8a4b`
