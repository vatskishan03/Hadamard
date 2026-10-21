# Hadamard matrix of order 26332

26332 = 4 · 6583, with 6583 = 29 · 227.

## Construction

Paley's skew-Hadamard matrix of order 908 (q = 907, a prime ≡ 3 mod 8) is a 4 × 4 array of circulant blocks of order
s = 227 once its points are listed along the four orbits of a torus element of PSL(2, 907) and rescaled by signs.
Multiplying it by T-matrices of order 29 (Seberry–Yamada 1990, Theorem 4) gives four circulant ±1 sequences of length
6583 with Σ PAF = 0 off 0, and the Goethals–Seidel array gives the Hadamard matrix. The construction and its proof
are in [PALEY_BLOCKS.md](../PALEY_BLOCKS.md).

- Torus element: the square of x ↦ −1/(x + τ) with τ = 3; sign normalisation c = +1.
- T-matrices of order 29 ([`t-matrices/T29.txt`](../../t-matrices/T29.txt), T-sequences of length 29 as stored in
  SageMath's `t_sequences` module): first rows `++00++-0-++00+000000000000000`, `00-+000+000+-0+00000000000000`,
  `000000000000000+000-0000+000-`, `0000000000000000+++0-+-+0---0`.
- Row sums of the four sequences: 47, −103, 67, −95 (47² + (−103)² + 67² + (−95)² = 26332).

## Search

None: the matrix is determined by q and the T-matrices. `python3 search/paley_blocks.py certify 907 29
t-matrices/T29.txt` writes the certificate.

## Status

**New, beyond the tables, but obtainable classically** ([PALEY_BLOCKS.md](../PALEY_BLOCKS.md#status)). In Table 1.53
of the Handbook of Combinatorial Designs (2nd ed., 2007), n = 6583 has t = 4: the smallest known Hadamard matrix with
odd part 6583 had order 2⁴ · 6583. The order lies above the range of Cati–Pasechnik Table 4, which stops at n = 2999.

- The entry does not appear in Đoković's updates (arXiv:0912.5091, arXiv:1008.2043).
- It follows from two published theorems that we found combined nowhere: Paley's matrix of order 908 is equivalent to
  an array of sixteen circulant blocks of order 227 (Schmidt 1999, after Ito; Yamada 1991), and Seberry and Yamada
  (1990, Theorem 4) multiply such an array by T-matrices of order 29.
- Paley's constructions do not apply: 4n − 1 = 26331 = 3 · 67 · 131 and 2n − 1 = 13165 = 5 · 2633 are not prime
  powers.
- Miyamoto's theorem does not apply: it needs a symmetric conference matrix of order n + 1, which cannot exist for n ≡ 3
  (mod 4). Kiyasu's theorem with q = n − 2 = 6581 needs a skew-Hadamard matrix of order 3292 = 4 · 823, and v = 823 is
  undecided in the Đoković–Seberry census.
- Cooper–Wallis does not apply: for each factorization n = a · b it needs T-matrices of order a and Williamson-type
  matrices of order b, and no such pair is known to us.
- SageMath's development version (17 October 2026) does not construct it: `hadamard_matrix(26332, existence=True)`
  returns Unknown.

We found no construction or statement of order 26332. The sources searched are listed in
[NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, τ, power, c, the T-matrices, the four sequences and the SHA-256 of the matrix |
| `sequences.txt` | the four sequences of length 6583 |

The matrix itself is not stored; `verify/verify.py` rebuilds it.

## Verification

```
python3 verify/verify.py 26332
```

SHA-256 of `H26332.csv`: `12e3d8cc89ca056806dcd17c5f8dab8a96d13bb421b6a10365d5d4357a6b1fea`
