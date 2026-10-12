# Order 21748

21748 = 4 · 5437, with 5437 prime and 5437 ≡ 1 (mod 4).

## Construction

Miyamoto's theorem: if q ≡ 1 (mod 4) is a prime power and a Hadamard matrix of
order q − 1 exists, then a Hadamard matrix of order 4q exists (M. Miyamoto,
A construction of Hadamard matrices, J. Combin. Theory Ser. A 57 (1991) 86–108).
The construction combines a symmetric conference matrix of order q + 1 with the
Hadamard matrix of order q − 1.

Here q = 5437. The inputs are

- Paley's symmetric conference matrix of order 5438: the matrix χ(y − x),
  x, y ∈ ℤ₅₄₃₇, with χ the quadratic character, bordered by a row and a column
  of ones;
- a Hadamard matrix K of order 5436 = 4 · 151 · 9. It is the Cooper–Wallis
  product of the T-matrices of order 151 in
  [`t-matrices/T151.txt`](../t-matrices/T151.txt) and the Williamson matrices of
  order 9 listed in SageMath, assembled as for [order 2060](../2060).

The blocks are arranged as in SageMath's `hadamard_matrix_miyamoto_construction`.
`build_miyamoto` in [`verify/verify.py`](../verify/verify.py) spells out the
layout, including the reordering of the conference matrix.

`certificate.json` contains q, the choice of conference matrix, the complete
certificate of K (T-matrices, Williamson matrices and the SHA-256 of K) and the
SHA-256 of the result.

## Status

Probably new. In Table 1.53 of the Handbook of Combinatorial Designs (2nd ed.,
2007), n = 5437 has t = 8: the smallest known Hadamard matrix with odd part 5437 had
order 2^8 · 5437. The input order 4 · 1359 has t = 3 there and in Cati–Pasechnik
Table 4. Neither entry appears in Đoković's updates (arXiv:0912.5091,
arXiv:1008.2043), and no classical family gives order 4 · 5437. The order has
followed from Miyamoto's theorem since Schumer published a matrix of order 5436 on
September 17, 2026, but we found no construction or statement of it. We read Table
1.53 from an unofficial copy, so treat this claim with that caveat.

The order is above the range of Cati–Pasechnik Table 4, which stops at n = 2999.

The input K is not new. Schumer built a Hadamard matrix of order 5436 the same
way, from T-matrices of order 151 and Williamson matrices of order 9 (Zenodo,
September 17, 2026, doi:10.5281/zenodo.22820997).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, conference matrix, certificate of K, SHA-256 |

The matrix is not stored. Its CSV encoding takes 1.2 GB (80 MB gzipped).

```
python3 verify/verify.py --write 21748
```

rebuilds and checks it and writes `21748/H21748.csv`. On six cores of an
AMD EPYC server this takes 30 seconds and 1 GB of memory.

SHA-256 of `H21748.csv`:
`9aec7ca4e8cf18716d2b4d25e9ee4e75cbb1ca0d94590f88b69fe9f8a442546c`
