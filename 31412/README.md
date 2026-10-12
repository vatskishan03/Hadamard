# Order 31412

31412 = 4 · 7853, with 7853 prime and 7853 ≡ 1 (mod 4).

## Construction

Miyamoto's theorem with q = 7853, as described for [order 21748](../21748). The
inputs are

- Paley's symmetric conference matrix of order 7854 over ℤ₇₈₅₃;
- a Hadamard matrix K of order 7852 = 4 · 151 · 13. It is the Cooper–Wallis
  product of the T-matrices of order 151 in
  [`t-matrices/T151.txt`](../t-matrices/T151.txt) and the Williamson matrices of
  order 13 listed in SageMath.

`certificate.json` contains q, the choice of conference matrix, the complete
certificate of K and the SHA-256 of the result.

## Status

Probably new. In Table 1.53 of the Handbook of Combinatorial Designs (2nd ed.,
2007), n = 7853 has t = 6: the smallest known Hadamard matrix with odd part 7853 had
order 2^6 · 7853. The input order 4 · 1963 has t = 3 there and in Cati–Pasechnik
Table 4. Neither entry appears in Đoković's updates (arXiv:0912.5091,
arXiv:1008.2043), and no classical family gives order 4 · 7853. The order has
followed from Miyamoto's theorem since Schumer published a matrix of order 7852 on
September 17, 2026, but we found no construction or statement of it. We read Table
1.53 from an unofficial copy, so treat this claim with that caveat.

The order is above the range of Cati–Pasechnik Table 4, which stops at n = 2999.

The input K is not new. Schumer built a Hadamard matrix of order 7852 the same
way, from T-matrices of order 151 and Williamson matrices of order 13 (Zenodo,
September 17, 2026, doi:10.5281/zenodo.22820997).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, conference matrix, certificate of K, SHA-256 |

The matrix is not stored. Its CSV encoding takes 2.5 GB (148 MB gzipped).

```
python3 verify/verify.py --write 31412
```

rebuilds and checks it and writes `31412/H31412.csv`. On six cores of an
AMD EPYC server this takes 80 seconds and 2 GB of memory.

SHA-256 of `H31412.csv`:
`1ce3cadb49a121c9ec9b4a0f62080babfb80309b84ec0a8c94ae8c60c3680da5`
