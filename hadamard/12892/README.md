# Hadamard matrix of order 12892

12892 = 4 · 3223, with 3223 = 11 · 293.

## Construction

Paley's skew-Hadamard matrix of order 1172 (q = 1171, a prime ≡ 3 mod 8) is a 4 × 4 array of circulant blocks of order
s = 293 once its points are listed along the four orbits of a torus element of PSL(2, 1171) and rescaled by signs.
Multiplying it by T-matrices of order 11 (Seberry–Yamada 1990, Theorem 4) gives four circulant ±1 sequences of length
3223 with Σ PAF = 0 off 0, and the Goethals–Seidel array gives the Hadamard matrix. The construction and its proof are
in [PALEY_BLOCKS.md](../PALEY_BLOCKS.md).

- Torus element: x ↦ −1/(x + τ) with τ = 4; sign normalisation c = −1.
- T-matrices of order 11 (first rows, from base sequences BS(6, 5)): `++00+000000`, `00++0-00000`, `000000++-+-`,
  `00000000000`.
- Row sums of the four sequences: 7, 31, −41, −101 (7² + 31² + (−41)² + (−101)² = 12892).

## Search

None: the matrix is determined by q and the T-matrices. `python3 search/paley_blocks.py certify 1171 11` writes the
certificate.

## Status

**First explicit construction** ([PALEY_BLOCKS.md](../PALEY_BLOCKS.md#orders-recorded-only-through-corollary-88i)).
Table 1.53 of the Handbook of Combinatorial Designs (2nd ed., 2007) gives n = 3223 with t = 2, as known. The only
construction we traced uses Williamson-type matrices of order 293 from Seberry–Yamada's Table A.1 (1992). Table A.1
obtains Williamson-type matrices of order 293 only from Corollary 8.8(i) (key wl, with Williamson-type matrices of order
73 = (293 − 1)/4), whose printed proof does not give Williamson-type matrices ([NOTES.md](../../NOTES.md) §7). Table A.1
credits that case to Miyamoto (J. Combin. Theory Ser. A 57 (1991)), whose paper was not accessible to us; we found no
valid proof of it. For q = 293 the corollary can be repaired: Seberry–Yamada's symmetric Williamson-type matrices of
order 73 (Corollary 8.8(ii) with Paley's symmetric Hadamard matrix of order 36) have a repeated block, and Lemma 26 in
the repeated-block form of [NOTES.md](../../NOTES.md) §7 turns them into Williamson-type matrices of order 293 (checked
exactly: `python3 search/symwt.py "cor29r(293, cor29h(73, paley2(17)))"`). So the order was obtainable from published
lemmas, by an arrangement we did not find in the literature.

- Paley's constructions do not apply: 4n − 1 = 12891 = 3 · 4297 and 2n − 1 = 6445 = 5 · 1289 are not prime powers.
- Miyamoto's theorem does not apply: it needs a symmetric conference matrix of order n + 1, which cannot exist for n ≡ 3
  (mod 4). Yamada–Kiyasu (n − 2 = 3221 a prime ≡ 5 mod 8) also gives it, from the skew-Hadamard matrix of order 1612 of
  this repository ([skew/1612](../../skew/1612)), found in October 2026.
- Cooper–Wallis does not apply otherwise: the other splits of n need Williamson-type matrices or T-matrices that are not
  known to us.
- SageMath's development version (17 October 2026) does not construct it: `hadamard_matrix(12892, existence=True)`
  returns Unknown.

We found no valid construction of order 12892 before this one. The sources searched are listed in
[NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, τ, power, c, the T-matrices, the four sequences and the SHA-256 of the matrix |
| `sequences.txt` | the four sequences of length 3223 |

The matrix itself is not stored; `verify/verify.py` rebuilds it.

## Verification

```
python3 verify/verify.py 12892
```

SHA-256 of `H12892.csv`: `60ad9080c13be3ca41b7a7bcd3c67746535f04416f1429d5fd1bf276e338aa71`
