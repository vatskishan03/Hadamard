# Hadamard matrix of order 37228

37228 = 4 · 9307, with 9307 = 41 · 227.

## Construction

Paley's skew-Hadamard matrix of order 908 (q = 907, a prime ≡ 3 mod 8) is a 4 × 4 array of circulant blocks of order
s = 227 once its points are listed along the four orbits of a torus element of PSL(2, 907) and rescaled by signs.
Multiplying it by T-matrices of order 41 (Seberry–Yamada 1990, Theorem 4) gives four circulant ±1 sequences of length
9307 with Σ PAF = 0 off 0, and the Goethals–Seidel array gives the Hadamard matrix. The construction and its proof
are in [PALEY_BLOCKS.md](../PALEY_BLOCKS.md).

- Torus element: the square of x ↦ −1/(x + τ) with τ = 3; sign normalisation c = +1.
- T-matrices of order 41 ([`t-matrices/T41.txt`](../../t-matrices/T41.txt), found by `search/tmat_mult.py 41 10`):
  first rows `+-00000000-00000-0-000000000000000000-000`, `0000+000000--0000000000+0+00-00+00-000-0+`,
  `00+000-0000000-+0-0-+0+0+0-+0000++0++0000`, `000+0-0+--000+0000000-0000000++00000000-0`.
- Row sums of the four sequences: −29, −21, −15, 189 ((−29)² + (−21)² + (−15)² + 189² = 37228).

## Search

None: the matrix is determined by q and the T-matrices. `python3 search/paley_blocks.py certify 907 41
t-matrices/T41.txt` writes the certificate.

## Status

**New, beyond the tables, but obtainable classically** ([PALEY_BLOCKS.md](../PALEY_BLOCKS.md#status)). In Table 1.53
of the Handbook of Combinatorial Designs (2nd ed., 2007), n = 9307 has t = 3: the smallest known Hadamard matrix with
odd part 9307 had order 2³ · 9307. The order lies above the range of Cati–Pasechnik Table 4, which stops at n = 2999.

- The entry does not appear in Đoković's updates (arXiv:0912.5091, arXiv:1008.2043).
- It follows from two published theorems that we found combined nowhere: Paley's matrix of order 908 is equivalent to
  an array of sixteen circulant blocks of order 227 (Schmidt 1999, after Ito; Yamada 1991), and Seberry and Yamada
  (1990, Theorem 4) multiply such an array by T-matrices of order 41.
- Paley's constructions do not apply: 4n − 1 = 37227 = 3 · 12409 and 2n − 1 = 18613 = 7 · 2659 are not prime powers.
- Miyamoto and Yamada–Kiyasu do not apply: 9307 is not a prime power, and neither is n − 2 = 9305 = 5 · 1861.
- Cooper–Wallis does not apply: for each factorization n = a · b it needs T-matrices of order a and Williamson-type
  matrices of order b, and no such pair is known to us.
- SageMath's development version (October 17, 2026) does not construct it: `hadamard_matrix(37228, existence=True)`
  returns Unknown.

We found no construction or statement of order 37228. The sources searched are listed in
[NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, τ, power, c, the T-matrices, the four sequences and the SHA-256 of the matrix |
| `sequences.txt` | the four sequences of length 9307 |

The matrix itself is not stored; `verify/verify.py` rebuilds it.

## Verification

```
python3 verify/verify.py 37228
```

SHA-256 of `H37228.csv`: `a9ed5b29aa4b29b0035b091da0c885b82de4303765718d589c5d905705d9ae1c`
