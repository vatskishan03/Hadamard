# Hadamard matrix of order 7612

7612 = 4 · 1903, with 1903 = 11 · 173.

## Construction

Paley's skew-Hadamard matrix of order 692 (q = 691, a prime ≡ 3 mod 8) is a 4 × 4 array of circulant blocks of order
s = 173 once its points are listed along the four orbits of a torus element of PSL(2, 691) and rescaled by signs.
Multiplying it by T-matrices of order 11 (Seberry–Yamada 1990, Theorem 4) gives four circulant ±1 sequences of length
1903 with Σ PAF = 0 off 0, and the Goethals–Seidel array gives the Hadamard matrix. The construction and its proof are
in [PALEY_BLOCKS.md](../PALEY_BLOCKS.md).

- Torus element: x ↦ −1/(x + τ) with τ = 4; sign normalisation c = +1.
- T-matrices of order 11 (first rows, from base sequences BS(6, 5)): `++00+000000`, `00++0-00000`, `000000++-+-`,
  `00000000000`.
- Row sums of the four sequences: −15, 5, 51, −69 ((−15)² + 5² + 51² + (−69)² = 7612).

## Search

None: the matrix is determined by q and the T-matrices. `python3 search/paley_blocks.py certify 691 11` writes the
certificate.

## Status

**First explicit construction** ([PALEY_BLOCKS.md](../PALEY_BLOCKS.md#orders-recorded-only-through-corollary-88i)).
Cati–Pasechnik's Table 4 omits n = 1903: it inherits m = 2 from Seberry–Yamada's 2020 Table A.17, and their 1992 Table
A.2 builds the order (code o2) from Williamson-type matrices of order 173. Table A.1 obtains Williamson-type matrices of
order 173 only from Corollary 8.8(i) (key wl, with Williamson-type matrices of order 43 = (173 − 1)/4), whose printed
proof does not give Williamson-type matrices ([NOTES.md](../../NOTES.md) §7). Table A.1 credits that case to Miyamoto
(J. Combin. Theory Ser. A 57 (1991)), whose paper was not accessible to us; we found no valid proof of it and no other
construction of Williamson-type matrices of order 173.

- Paley's constructions do not apply: 4n − 1 = 7611 = 3 · 43 · 59 and 2n − 1 = 3805 = 5 · 761 are not prime powers.
- Miyamoto's theorem does not apply: it needs a symmetric conference matrix of order n + 1, which cannot exist for n ≡ 3
  (mod 4). Yamada–Kiyasu (n − 2 = 1901 a prime ≡ 5 mod 8) would need a skew-Hadamard matrix of order 952 = 8 · 119; none
  is known to us (119 is undecided in the Đoković–Seberry census of skew-Hadamard matrices of order 4v, so doubling does
  not give one either).
- Cooper–Wallis does not apply otherwise: the other splits of n need Williamson-type matrices or T-matrices that are not
  known to us.
- SageMath's development version (October 17, 2026) does not construct it: `hadamard_matrix(7612, existence=True)`
  returns Unknown.

We found no valid construction of order 7612 before this one. The sources searched are listed in
[NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, τ, power, c, the T-matrices, the four sequences and the SHA-256 of the matrix |
| `sequences.txt` | the four sequences of length 1903 |

The matrix itself is not stored; `verify/verify.py` rebuilds it.

## Verification

```
python3 verify/verify.py 7612
```

SHA-256 of `H7612.csv`: `b1ea8cea14548ebc0af1da5161007f60a995a37859c94c7be184920324a3f618`
