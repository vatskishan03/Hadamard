# Hadamard matrix of order 35596

35596 = 4 · 8899, with 8899 = 11 · 809.

## Construction

Miyamoto's Hadamard matrix of order 3236 (q = 809, a prime ≡ 1 mod 4), built from Paley's conference matrix of order 810
and a Hadamard matrix of order 808 made of two circulant blocks, is a 4 × 4 array of bordered blocks of order 809 that
commute and are transposed by the reversal j ↦ −j. Multiplying it by T-matrices of order 11 in the Goethals–Seidel array
gives the Hadamard matrix of order 4 · 809 · 11 = 35596. The construction and its proof are in
[MIYAMOTO_BLOCKS.md](../MIYAMOTO_BLOCKS.md#periodic-golay-pairs).

- Paley's conference matrix of order 810, its points listed as ∞, 0, g^(2j), g^(2j+1) (j < 404) with the primitive root
  g = 3.
- The Hadamard matrix of order 808 is [[A, B], [−Bᵀ, Aᵀ]] with A, B circulant of order 404: a periodic Golay pair of
  length 404 = 2 · 202, the product (formulas (15) and (16) of Balonin–Đoković, which hold for periodic pairs as well;
  Georgiou, Stylianou, Drosou and Koukouvinos 2014, Theorem 1) of the Golay pair of length 2 and the periodic Golay pair
  of length 202 of Đoković and Kotsireas (Numer. Algorithms 69 (2015); arXiv:1310.5773, orbits of H = {1, 87, 95, 137,
  185} on ℤ₂₀₂).
- T-matrices of order 11 (first rows, from base sequences BS(6, 5)): `++00+000000`, `00++0-00000`, `000000++-+-`,
  `00000000000`.

## Search

None: the matrix is determined by q, the periodic Golay pair and the T-matrices. `python3 search/miyamoto_blocks.py
certify 809 11 golay2:dk202` writes the certificate.

## Status

**First explicit construction** ([MIYAMOTO_BLOCKS.md](../MIYAMOTO_BLOCKS.md#periodic-golay-pairs)). Table 1.53 of the
Handbook of Combinatorial Designs (2nd ed., 2007) gives n = 8899 with t = 3, and Đoković (arXiv:0912.5091) lists it in
his set Δ of entries that should have had t = 2. His proof for Δ takes Williamson-type matrices from Seberry–Yamada's
Table A.1 (1992); for 8899 they are of order 809. Table A.1 obtains Williamson-type matrices of order 809 only from
Corollary 8.8(i), whose printed proof does not give Williamson-type matrices ([NOTES.md](../../NOTES.md) §7). Table A.1
credits that case to Miyamoto (J. Combin. Theory Ser. A 57 (1991)), whose paper was not accessible to us; we found no
valid proof of it and no other construction of Williamson-type matrices of order 809.

- Paley's constructions do not apply: 4n − 1 = 35595 = 3² · 5 · 7 · 113 and 2n − 1 = 17797 = 13 · 37² are not prime
  powers.
- Miyamoto and Yamada–Kiyasu do not apply: 8899 is not a prime power, and n − 2 = 8897 = 7 · 31 · 41 is not a prime
  power.
- Cooper–Wallis needs Williamson-type matrices of order 809 or T-matrices of order 809; neither is known to us.
- The Paley-block route does not reach it: 4s − 1 = 3235 and 35595 are not prime powers ≡ 3 (mod 8) for the other
  divisors s of 8899, and the remaining ones need T-matrices of order 8899 (s = 1) or 809 (s = 11), which are not known
  to us. Nor does the Miyamoto-block route with a negaperiodic Golay pair: none of length 404 is known (Balonin–Đoković,
  list (31)).
- SageMath's development version (October 17, 2026) does not construct it: `hadamard_matrix(35596, existence=True)`
  returns Unknown.

We found no valid construction of order 35596 before this one. The sources searched are listed in
[NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, g, the periodic Golay pair, the T-matrices and the SHA-256 of the matrix |

The matrix itself is not stored; `verify/verify.py` rebuilds it.

## Verification

```
python3 verify/verify.py 35596
```

SHA-256 of `H35596.csv`: `94bde853b23b5468e4945d8a58d7103fa6bf4ac9e7afef8377dc26de05bb0ab1`
