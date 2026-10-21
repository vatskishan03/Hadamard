# Hadamard matrix of order 8204

8204 = 4 · 2051, with 2051 = 7 · 293.

## Construction

Paley's skew-Hadamard matrix of order 1172 (q = 1171, a prime ≡ 3 mod 8) is a 4 × 4 array of circulant blocks of order
s = 293 once its points are listed along the four orbits of a torus element of PSL(2, 1171) and rescaled by signs.
Multiplying it by T-matrices of order 7 (Seberry–Yamada 1990, Theorem 4) gives four circulant ±1 sequences of length
2051 with Σ PAF = 0 off 0, and the Goethals–Seidel array gives the Hadamard matrix. The construction and its proof are
in [PALEY_BLOCKS.md](../PALEY_BLOCKS.md).

- Torus element: x ↦ −1/(x + τ) with τ = 4; sign normalisation c = −1.
- T-matrices of order 7 (first rows, from base sequences BS(4, 3)): `+000000`, `0++-000`, `0000+0+`, `00000+0`.
- Row sums of the four sequences: 47, −3, 19, −75 (47² + (−3)² + 19² + (−75)² = 8204).

## Search

None: the matrix is determined by q and the T-matrices. `python3 search/paley_blocks.py certify 1171 7` writes the
certificate.

## Status

**First explicit construction** ([PALEY_BLOCKS.md](../PALEY_BLOCKS.md#orders-recorded-only-through-corollary-88i)).
Cati–Pasechnik's Table 4 omits n = 2051 (m = 2, which by Cati–Pasechnik §6 and their footnote 11 comes from
Seberry–Yamada's 2020 Table A.17; we did not read that table), and their 1992 Table A.2 builds the order (code o2) from
Williamson-type matrices of order 293. Table A.1 obtains Williamson-type matrices of order 293 only from Corollary
8.8(i) (key wl, with Williamson-type matrices of order 73 = (293 − 1)/4), whose printed proof does not give
Williamson-type matrices ([NOTES.md](../../NOTES.md) §7). Table A.1 credits that case to Miyamoto (J. Combin. Theory
Ser. A 57 (1991)), whose paper was not accessible to us; we found no valid proof of it. For q = 293 the corollary can be
repaired: Seberry–Yamada's symmetric Williamson-type matrices of order 73 (Corollary 8.8(ii) with Paley's symmetric
Hadamard matrix of order 36) have a repeated block, and Lemma 26 in the repeated-block form of
[NOTES.md](../../NOTES.md) §7 turns them into Williamson-type matrices of order 293 (checked exactly: `python3
search/symwt.py "cor29r(293, cor29h(73, paley2(17)))"`). So the order was obtainable from published lemmas, by an
arrangement we did not find in the literature.

- Paley's constructions do not apply: 4n − 1 = 8203 = 13 · 631 and 2n − 1 = 4101 = 3 · 1367 are not prime powers.
- Miyamoto's theorem does not apply: it needs a symmetric conference matrix of order n + 1, which cannot exist for n ≡ 3
  (mod 4); Yamada–Kiyasu does not apply: n − 2 = 2049 = 3 · 683 is not a prime power.
- Cooper–Wallis does not apply otherwise: the other splits of n need Williamson-type matrices or T-matrices that are not
  known to us.
- SageMath's development version (17 October 2026) does not construct it: `hadamard_matrix(8204, existence=True)`
  returns Unknown.

We found no valid construction of order 8204 before this one. The sources searched are listed in
[NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, τ, power, c, the T-matrices, the four sequences and the SHA-256 of the matrix |
| `sequences.txt` | the four sequences of length 2051 |

The matrix itself is not stored; `verify/verify.py` rebuilds it.

## Verification

```
python3 verify/verify.py 8204
```

SHA-256 of `H8204.csv`: `94f5e4fe453486627c03ff778d0665eca999b1c020103a13f72cc5d4b5bc7a02`
