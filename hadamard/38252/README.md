# Hadamard matrix of order 38252

38252 = 4 · 9563, with 9563 = 73 · 131.

## Construction

Paley's skew-Hadamard matrix of order 524 (q = 523, a prime ≡ 3 mod 8) is a 4 × 4 array of circulant blocks of order
s = 131 once its points are listed along the four orbits of a torus element of PSL(2, 523) and rescaled by signs.
Multiplying it by T-matrices of order 73 (Seberry–Yamada 1990, Theorem 4) gives four circulant ±1 sequences of length
9563 with Σ PAF = 0 off 0, and the Goethals–Seidel array gives the Hadamard matrix. The construction and its proof
are in [PALEY_BLOCKS.md](../PALEY_BLOCKS.md).

- Torus element: the square of x ↦ −1/(x + τ) with τ = 3; sign normalisation c = +1.
- T-matrices of order 73 ([`t-matrices/T73.txt`](../../t-matrices/T73.txt), found by `search/tmat_mult.py 73 8`):
  first rows in `certificate.json`.
- Row sums of the four sequences: −79, 39, −167, −51 ((−79)² + 39² + (−167)² + (−51)² = 38252).

## Search

None: the matrix is determined by q and the T-matrices. `python3 search/paley_blocks.py certify 523 73
t-matrices/T73.txt` writes the certificate.

## Status

**New, beyond the tables, but obtainable classically** ([PALEY_BLOCKS.md](../PALEY_BLOCKS.md#status)). In Table 1.53
of the Handbook of Combinatorial Designs (2nd ed., 2007), n = 9563 has t = 3: the smallest known Hadamard matrix with
odd part 9563 had order 2³ · 9563. The order lies above the range of Cati–Pasechnik Table 4, which stops at n = 2999.

- The entry does not appear in Đoković's updates (arXiv:0912.5091, arXiv:1008.2043).
- It follows from two published theorems that we found combined nowhere: Paley's matrix of order 524 is equivalent to
  an array of sixteen circulant blocks of order 131 (Schmidt 1999, after Ito; Yamada 1991), and Seberry and Yamada
  (1990, Theorem 4) multiply such an array by T-matrices of order 73.
- Paley's constructions do not apply: 4n − 1 = 38251 = 29 · 1319 and 2n − 1 = 19125 = 3² · 5³ · 17 are not prime
  powers.
- Miyamoto's theorem does not apply: it needs a symmetric conference matrix of order n + 1, which cannot exist for n ≡ 3
  (mod 4); Yamada–Kiyasu does not apply: n − 2 = 9561 = 3 · 3187 is not a prime power.
- SageMath's development version (17 October 2026) does not construct it: `hadamard_matrix(38252, existence=True)`
  returns Unknown.

A second route reaches the same order: Cooper–Wallis with T-matrices of order 131 and Williamson-type matrices of
order 73 (Seberry–Yamada 1992, Table A.1). T-sequences of length 131 follow from the Turyn-type sequences TT(44) of
S. London and I. Kotsireas (Cryptogr. Commun. 17 (2025) 1601–1610), which we could not read; if that paper lists
Hadamard orders, this one may be among them. Our T-matrices of order 131 are in
[`t-matrices/T131.txt`](../../t-matrices/T131.txt).

We found no construction or statement of order 38252. The sources searched are listed in
[NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, τ, power, c, the T-matrices, the four sequences and the SHA-256 of the matrix |
| `sequences.txt` | the four sequences of length 9563 |

The matrix itself is not stored; `verify/verify.py` rebuilds it.

## Verification

```
python3 verify/verify.py 38252
```

SHA-256 of `H38252.csv`: `0bdb4c312c909baf6e02c7359c785d1a17953ff9173ce827d91bf992e509d46d`
