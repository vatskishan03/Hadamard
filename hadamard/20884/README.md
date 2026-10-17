# Hadamard matrix of order 20884

20884 = 4 · 5221, with 5221 = 23 · 227.

## Construction

Paley's skew-Hadamard matrix of order 908 (q = 907, a prime ≡ 3 mod 8) is a 4 × 4 array of circulant blocks of order
s = 227 once its points are listed along the four orbits of a torus element of PSL(2, 907) and rescaled by signs.
Multiplying it by T-matrices of order 23 (Seberry–Yamada 1990, Theorem 4) gives four circulant ±1 sequences of length
5221 with Σ PAF = 0 off 0, and the Goethals–Seidel array gives the Hadamard matrix. The construction and its proof
are in [PALEY_BLOCKS.md](../PALEY_BLOCKS.md).

- Torus element: the square of x ↦ −1/(x + τ) with τ = 3; sign normalisation c = +1.
- T-matrices of order 23 ([`t-matrices/T23.txt`](../../t-matrices/T23.txt), T-sequences of length 23 as stored in
  SageMath's `t_sequences` module): first rows `+0000000000000000000000`, `0+++++--+-+-00000000000`,
  `000000000000+0+0-0+0-0+`, `0000000000000+0+0-0+0+0`.
- Row sums of the four sequences: 117, −75, 7, 39 (117² + (−75)² + 7² + 39² = 20884).

## Search

None: the matrix is determined by q and the T-matrices. `python3 search/paley_blocks.py certify 907 23
t-matrices/T23.txt` writes the certificate.

## Status

**New, beyond the tables, but obtainable classically** ([PALEY_BLOCKS.md](../PALEY_BLOCKS.md#status)). In Table 1.53
of the Handbook of Combinatorial Designs (2nd ed., 2007), n = 5221 has t = 3: the smallest known Hadamard matrix with
odd part 5221 had order 2³ · 5221. The order lies above the range of Cati–Pasechnik Table 4, which stops at n = 2999.

- The entry does not appear in Đoković's updates (arXiv:0912.5091, arXiv:1008.2043).
- It follows from two published theorems that we found combined nowhere: Paley's matrix of order 908 is equivalent to
  an array of sixteen circulant blocks of order 227 (Schmidt 1999, after Ito; Yamada 1991), and Seberry and Yamada
  (1990, Theorem 4) multiply such an array by T-matrices of order 23.
- Paley's constructions do not apply: 4n − 1 = 20883 = 3 · 6961 and 2n − 1 = 10441 = 53 · 197 are not prime powers.
- Miyamoto's theorem does not apply: a symmetric conference matrix of order n + 1 = 5222 cannot exist (5221 is not a sum
  of two squares); Yamada–Kiyasu does not apply: n − 2 = 5219 = 17 · 307 is not a prime power.
- Cooper–Wallis does not apply: for each factorization n = a · b it needs T-matrices of order a and Williamson-type
  matrices of order b, and no such pair is known to us.
- SageMath's development version (October 17, 2026) does not construct it: `hadamard_matrix(20884, existence=True)`
  returns Unknown.

We found no construction or statement of order 20884. The sources searched are listed in
[NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, τ, power, c, the T-matrices, the four sequences and the SHA-256 of the matrix |
| `sequences.txt` | the four sequences of length 5221 |

The matrix itself is not stored; `verify/verify.py` rebuilds it.

## Verification

```
python3 verify/verify.py 20884
```

SHA-256 of `H20884.csv`: `8dd4fa20c7bde31fdb6bd11c202004e9a329790a33bb016c9c22b0a6407ce7b5`
