# Hadamard matrix of order 30916

30916 = 4 · 7729, with 7729 = 59 · 131.

## Construction

Paley's skew-Hadamard matrix of order 524 (q = 523, a prime ≡ 3 mod 8) is a 4 × 4 array of circulant blocks of order
s = 131 once its points are listed along the four orbits of a torus element of PSL(2, 523) and rescaled by signs.
Multiplying it by T-matrices of order 59 (Seberry–Yamada 1990, Theorem 4) gives four circulant ±1 sequences of length
7729 with Σ PAF = 0 off 0, and the Goethals–Seidel array gives the Hadamard matrix. The construction and its proof
are in [PALEY_BLOCKS.md](../PALEY_BLOCKS.md).

- Torus element: the square of x ↦ −1/(x + τ) with τ = 3; sign normalisation c = +1.
- T-matrices of order 59 ([`t-matrices/T59.txt`](../../t-matrices/T59.txt), T-sequences of length 59 as stored in
  SageMath's `t_sequences` module): first rows in `certificate.json`.
- Row sums of the four sequences: 7, −29, 85, 151 (7² + (−29)² + 85² + 151² = 30916).

## Search

None: the matrix is determined by q and the T-matrices. `python3 search/paley_blocks.py certify 523 59
t-matrices/T59.txt` writes the certificate.

## Status

**New, beyond the tables, but obtainable classically** ([PALEY_BLOCKS.md](../PALEY_BLOCKS.md#status)). In Table 1.53
of the Handbook of Combinatorial Designs (2nd ed., 2007), n = 7729 has t = 3: the smallest known Hadamard matrix with
odd part 7729 had order 2³ · 7729. The order lies above the range of Cati–Pasechnik Table 4, which stops at n = 2999.

- The entry does not appear in Đoković's updates (arXiv:0912.5091, arXiv:1008.2043).
- It follows from two published theorems that we found combined nowhere: Paley's matrix of order 524 is equivalent to
  an array of sixteen circulant blocks of order 131 (Schmidt 1999, after Ito; Yamada 1991), and Seberry and Yamada
  (1990, Theorem 4) multiply such an array by T-matrices of order 59.
- Paley's constructions do not apply: 4n − 1 = 30915 = 3³ · 5 · 229 and 2n − 1 = 15457 = 13 · 29 · 41 are not prime
  powers.
- Miyamoto and Yamada–Kiyasu do not apply: 7729 is not a prime power, and n − 2 = 7727 is ≡ 7 (mod 8).
- Cooper–Wallis does not apply: for each factorization n = a · b it needs T-matrices of order a and Williamson-type
  matrices of order b, and no such pair is known to us.
- SageMath's development version (October 17, 2026) does not construct it: `hadamard_matrix(30916, existence=True)`
  returns Unknown.

We found no construction or statement of order 30916. The sources searched are listed in
[NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, τ, power, c, the T-matrices, the four sequences and the SHA-256 of the matrix |
| `sequences.txt` | the four sequences of length 7729 |

The matrix itself is not stored; `verify/verify.py` rebuilds it.

## Verification

```
python3 verify/verify.py 30916
```

SHA-256 of `H30916.csv`: `22b2f84da7acd2a0a73672c69be0045fa7be6e662117dabb3c3c307612ade9eb`
