# Hadamard matrix of order 37300

37300 = 4 · 9325, with 9325 = 5² · 373.

## Construction

Paley's skew-Hadamard matrix of order 7460 (q = 7459, a prime ≡ 3 mod 8) is a 4 × 4 array of circulant blocks of order
s = 1865 once its points are listed along the four orbits of a torus element of PSL(2, 7459) and rescaled by signs.
Multiplying it by T-matrices of order 5 (Seberry–Yamada 1990, Theorem 4) gives four ±1 sequences on the group ℤ₅ × ℤ₁₈₆₅
of order 9325 with Σ PAF = 0 off 0, and the Goethals–Seidel array gives the Hadamard matrix. The construction and its
proof are in [PALEY_BLOCKS.md](../PALEY_BLOCKS.md).

- Torus element: x ↦ −1/(x + τ) with τ = 4; sign normalisation c = −1.
- T-matrices of order 5 (first rows, from base sequences BS(3, 2)): `++000`, `00+00`, `000+-`, `00000`.
- Row sums of the four sequences: 33, −61, −57, −171 (33² + (−61)² + (−57)² + (−171)² = 37300).
- Since gcd(s, t) = 5, the sequences live on ℤ_5 × ℤ_1865, not on ℤ_9325. `sequences.txt` lists the entry at (i, a) in
  position i · 1865 + a, and the Goethals–Seidel array uses the group-developed matrices and the reversal x ↦ −x of the
  group.

## Search

None: the matrix is determined by q and the T-matrices. `python3 search/paley_blocks.py certify 7459 5` writes the
certificate.

## Status

**First explicit construction** ([PALEY_BLOCKS.md](../PALEY_BLOCKS.md#orders-recorded-only-through-corollary-88i)).
Table 1.53 of the Handbook of Combinatorial Designs (2nd ed., 2007) gives n = 9325 with t = 3, and Đoković
(arXiv:0912.5091) lists it in his set Δ of entries that should have had t = 2. His proof for Δ takes Williamson-type
matrices from Seberry–Yamada's Table A.1 (1992); for 9325 they are of order 373. Table A.1 obtains Williamson-type
matrices of order 373 only from Corollary 8.8(i) (key wl, with Williamson-type matrices of order 93 = (373 − 1)/4),
whose printed proof does not give Williamson-type matrices ([NOTES.md](../../NOTES.md) §7). Table A.1 credits that case
to Miyamoto (J. Combin. Theory Ser. A 57 (1991)), whose paper was not accessible to us; we found no valid proof of it
and no other construction of Williamson-type matrices of order 373.

- Paley's constructions do not apply: 4n − 1 = 37299 = 3 · 12433 and 2n − 1 = 18649 = 17 · 1097 are not prime powers.
- Miyamoto's theorem does not apply: it needs a symmetric conference matrix of order n + 1 = 9326, and none is known to
  us (9325 is not a prime power); Yamada–Kiyasu does not apply: n − 2 = 9323 is ≡ 3 (mod 8).
- Cooper–Wallis does not apply otherwise: the other splits of n need Williamson-type matrices or T-matrices that are not
  known to us.
- SageMath's development version (October 17, 2026) does not construct it: `hadamard_matrix(37300, existence=True)`
  raises UnknownError in its Miyamoto test (no symmetric conference matrix of order 9326 is known to it), and returns
  Unknown when that Unknown is treated as False.

We found no valid construction of order 37300 before this one. The sources searched are listed in
[NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, τ, power, c, the T-matrices, the four sequences and the SHA-256 of the matrix |
| `sequences.txt` | the four sequences of length 9325 |

The matrix itself is not stored; `verify/verify.py` rebuilds it.

## Verification

```
python3 verify/verify.py 37300
```

SHA-256 of `H37300.csv`: `bb9b910bf524d7580529dea7e9f1ed8189d72b8ff884827c7fcef8d23281d200`
