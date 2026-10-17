# Hadamard matrix of order 9708

9708 = 4 · 2427, with 2427 = 3 · 809.

## Construction

Miyamoto's Hadamard matrix of order 3236 (q = 809, a prime ≡ 1 mod 4), built from Paley's conference matrix of order 810
and a Hadamard matrix of order 808 made of two circulant blocks, is a 4 × 4 array of bordered blocks of order 809 that
commute and are transposed by the reversal j ↦ −j. Multiplying it by T-matrices of order 3 in the Goethals–Seidel array
gives the Hadamard matrix of order 4 · 809 · 3 = 9708. The construction and its proof are in
[MIYAMOTO_BLOCKS.md](../MIYAMOTO_BLOCKS.md#periodic-golay-pairs).

- Paley's conference matrix of order 810, its points listed as ∞, 0, g^(2j), g^(2j+1) (j < 404) with the primitive root
  g = 3.
- The Hadamard matrix of order 808 is [[A, B], [−Bᵀ, Aᵀ]] with A, B circulant of order 404: a periodic Golay pair of
  length 404 = 2 · 202, the product (formulas (15) and (16) of Balonin–Đoković, which hold for periodic pairs as well;
  Georgiou, Stylianou, Drosou and Koukouvinos 2014, Theorem 1) of the Golay pair of length 2 and the periodic Golay pair
  of length 202 of Đoković and Kotsireas (Numer. Algorithms 69 (2015); arXiv:1310.5773, orbits of H = {1, 87, 95, 137,
  185} on ℤ₂₀₂).
- T-matrices of order 3 (first rows, from base sequences BS(2, 1)): `+00`, `0+0`, `00+`, `000`.

## Search

None beyond the Williamson matrices: the matrix is determined by q, the pair and the T-matrices. `python3
search/miyamoto_blocks.py certify 809 3 golay2:dk202` writes the certificate.

## Status

**First explicit construction** ([MIYAMOTO_BLOCKS.md](../MIYAMOTO_BLOCKS.md#periodic-golay-pairs)). Cati–Pasechnik's
Table 4 omits n = 2427 (m = 2, which by Cati–Pasechnik §6 and their footnote 11 comes from Seberry–Yamada's 2020 Table
A.17; we did not read that table), and their 1992 Table A.2 builds the order (code o2) from Williamson-type matrices of
order 809. Table A.1 obtains Williamson-type matrices of order 809 only from Corollary 8.8(i), whose printed proof does
not give Williamson-type matrices ([NOTES.md](../../NOTES.md) §7). Table A.1 credits that case to Miyamoto (J. Combin.
Theory Ser. A 57 (1991)), whose paper was not accessible to us; we found no valid proof of it and no other construction
of Williamson-type matrices of order 809.

- Paley's constructions do not apply: 4n − 1 = 9707 = 17 · 571 and 2n − 1 = 4853 = 23 · 211 are not prime powers.
- Miyamoto's theorem does not apply: it needs a symmetric conference matrix of order n + 1, which cannot exist for n ≡ 3
  (mod 4); Yamada–Kiyasu does not apply: n − 2 = 2425 = 5² · 97 is not a prime power.
- Cooper–Wallis needs Williamson-type matrices of order 809 or T-matrices of order 809; neither is known to us.
- The Paley-block route does not reach it: 4s − 1 = 3235 and 9707 are not prime powers ≡ 3 (mod 8) for the other
  divisors s of 2427, and the remaining ones need T-matrices of order 2427 (s = 1) or 809 (s = 3), which are not known
  to us. Nor does the Miyamoto-block route with a negaperiodic Golay pair: none of length 404 is known (Balonin–Đoković,
  list (31)).
- SageMath's development version (October 17, 2026) does not construct it: `hadamard_matrix(9708, existence=True)`
  returns Unknown.

We found no valid construction of order 9708 before this one. The sources searched are listed in
[NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, g, the periodic Golay pair, the T-matrices and the SHA-256 of the matrix |

The matrix itself is not stored; `verify/verify.py` rebuilds it.

## Verification

```
python3 verify/verify.py 9708
```

SHA-256 of `H9708.csv`: `6152893e9822fd46cab53ff457625d2b8c9686b716e63efb8b14ccf67c90e212`
