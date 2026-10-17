# Hadamard matrix of order 21484

21484 = 4 · 5371, with 5371 = 41 · 131.

## Construction

Paley's skew-Hadamard matrix of order 524 (q = 523, a prime ≡ 3 mod 8) is a 4 × 4 array of circulant blocks of order
s = 131 once its points are listed along the four orbits of a torus element of PSL(2, 523) and rescaled by signs.
Multiplying it by T-matrices of order 41 (Seberry–Yamada 1990, Theorem 4) gives four circulant ±1 sequences of length
5371 with Σ PAF = 0 off 0, and the Goethals–Seidel array gives the Hadamard matrix. The construction and its proof
are in [PALEY_BLOCKS.md](../PALEY_BLOCKS.md).

- Torus element: the square of x ↦ −1/(x + τ) with τ = 3; sign normalisation c = +1.
- T-matrices of order 41 ([`t-matrices/T41.txt`](../../t-matrices/T41.txt), found by `search/tmat_mult.py 41 10`):
  first rows `+-00000000-00000-0-000000000000000000-000`, `0000+000000--0000000000+0+00-00+00-000-0+`,
  `00+000-0000000-+0-0-+0+0+0-+0000++0++0000`, `000+0-0+--000+0000000-0000000++00000000-0`.
- Row sums of the four sequences: −99, −9, −71, 81 ((−99)² + (−9)² + (−71)² + 81² = 21484).

## Search

None: the matrix is determined by q and the T-matrices. `python3 search/paley_blocks.py certify 523 41
t-matrices/T41.txt` writes the certificate.

## Status

**New, beyond the tables, but obtainable classically** ([PALEY_BLOCKS.md](../PALEY_BLOCKS.md#status)). In Table 1.53
of the Handbook of Combinatorial Designs (2nd ed., 2007), n = 5371 has t = 3: the smallest known Hadamard matrix with
odd part 5371 had order 2³ · 5371. The order lies above the range of Cati–Pasechnik Table 4, which stops at n = 2999.

- The entry does not appear in Đoković's updates (arXiv:0912.5091, arXiv:1008.2043).
- It follows from two published theorems that we found combined nowhere: Paley's matrix of order 524 is equivalent to
  an array of sixteen circulant blocks of order 131 (Schmidt 1999, after Ito; Yamada 1991), and Seberry and Yamada
  (1990, Theorem 4) multiply such an array by T-matrices of order 41.
- Paley's constructions do not apply: 4n − 1 = 21483 = 3² · 7 · 11 · 31 and 2n − 1 = 10741 = 23 · 467 are not prime
  powers.
- Miyamoto and Yamada–Kiyasu do not apply: 5371 is not a prime power, and neither is n − 2 = 5369 = 7 · 13 · 59.
- SageMath's development version (October 17, 2026) does not construct it: `hadamard_matrix(21484, existence=True)`
  returns Unknown.

A second route reaches the same order: Cooper–Wallis with T-matrices of order 131 and Williamson-type matrices of
order 41 (Seberry–Yamada 1992, Table A.1). T-sequences of length 131 follow from the Turyn-type sequences TT(44) of
S. London and I. Kotsireas (Cryptogr. Commun. 17 (2025) 1601–1610), which we could not read; if that paper lists
Hadamard orders, this one may be among them. Our T-matrices of order 131 are in
[`t-matrices/T131.txt`](../../t-matrices/T131.txt).

We found no construction or statement of order 21484. The sources searched are listed in
[NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, τ, power, c, the T-matrices, the four sequences and the SHA-256 of the matrix |
| `sequences.txt` | the four sequences of length 5371 |

The matrix itself is not stored; `verify/verify.py` rebuilds it.

## Verification

```
python3 verify/verify.py 21484
```

SHA-256 of `H21484.csv`: `300e57eb6ed9581d8ec9cd12a479dd39801f1c852f558f36cd6098fc25ca2111`
