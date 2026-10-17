# Hadamard matrix of order 36332

36332 = 4 · 9083, with 9083 = 31 · 293.

## Construction

Paley's skew-Hadamard matrix of order 1172 (q = 1171, a prime ≡ 3 mod 8) is a 4 × 4 array of circulant blocks of order
s = 293 once its points are listed along the four orbits of a torus element of PSL(2, 1171) and rescaled by signs.
Multiplying it by T-matrices of order 31 (Seberry–Yamada 1990, Theorem 4) gives four circulant ±1 sequences of length
9083 with Σ PAF = 0 off 0, and the Goethals–Seidel array gives the Hadamard matrix. The construction and its proof are
in [PALEY_BLOCKS.md](../PALEY_BLOCKS.md).

- Torus element: x ↦ −1/(x + τ) with τ = 4; sign normalisation c = −1.
- T-matrices of order 31 ([`t-matrices/T31.txt`](../../t-matrices/T31.txt), sources in
  [t-matrices/README.md](../../t-matrices/README.md)).
- Row sums of the four sequences: −53, 95, 153, 33 ((−53)² + 95² + 153² + 33² = 36332).

## Search

None: the matrix is determined by q and the T-matrices. `python3 search/paley_blocks.py certify 1171 31
t-matrices/T31.txt` writes the certificate.

## Status

**First explicit construction** ([PALEY_BLOCKS.md](../PALEY_BLOCKS.md#orders-recorded-only-through-corollary-88i)).
Table 1.53 of the Handbook of Combinatorial Designs (2nd ed., 2007) gives n = 9083 with t = 3, and Đoković
(arXiv:0912.5091) lists it in his set Δ of entries that should have had t = 2. His proof for Δ takes Williamson-type
matrices from Seberry–Yamada's Table A.1 (1992); for 9083 they are of order 293. Table A.1 obtains Williamson-type
matrices of order 293 only from Corollary 8.8(i) (key wl, with Williamson-type matrices of order 73 = (293 − 1)/4),
whose printed proof does not give Williamson-type matrices ([NOTES.md](../../NOTES.md) §7). Table A.1 credits that case
to Miyamoto (J. Combin. Theory Ser. A 57 (1991)), whose paper was not accessible to us; we found no valid proof of it.
For q = 293 the corollary can be repaired: Seberry–Yamada's symmetric Williamson-type matrices of order 73 (Corollary
8.8(ii) with Paley's symmetric Hadamard matrix of order 36) have a repeated block, and Lemma 26 in the repeated-block
form of [NOTES.md](../../NOTES.md) §7 turns them into Williamson-type matrices of order 293 (checked exactly: `python3
search/symwt.py "cor29r(293, cor29h(73, paley2(17)))"`). So the order was obtainable from published lemmas, by an
arrangement we did not find in the literature.

- Paley's constructions do not apply: 4n − 1 = 36331 = 47 · 773 and 2n − 1 = 18165 = 3 · 5 · 7 · 173 are not prime
  powers.
- Miyamoto's theorem does not apply: it needs a symmetric conference matrix of order n + 1, which cannot exist for n ≡ 3
  (mod 4); Yamada–Kiyasu does not apply: n − 2 = 9081 = 3² · 1009 is not a prime power.
- Cooper–Wallis does not apply otherwise: the other splits of n need Williamson-type matrices or T-matrices that are not
  known to us.
- SageMath's development version (October 17, 2026) does not construct it: `hadamard_matrix(36332, existence=True)`
  returns Unknown.

We found no valid construction of order 36332 before this one. The sources searched are listed in
[NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, τ, power, c, the T-matrices, the four sequences and the SHA-256 of the matrix |
| `sequences.txt` | the four sequences of length 9083 |

The matrix itself is not stored; `verify/verify.py` rebuilds it.

## Verification

```
python3 verify/verify.py 36332
```

SHA-256 of `H36332.csv`: `7f8c0b916928d26f88597662c9b90afdd0a64cf427c2bcefd185d97bf8297933`
