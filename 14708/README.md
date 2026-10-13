# Order 14708

14708 = 4 · 3677, with 3677 prime and 3677 ≡ 1 (mod 4).

## Construction

Miyamoto's theorem: if q ≡ 1 (mod 4) is a prime power and a Hadamard matrix of
order q − 1 exists, then a Hadamard matrix of order 4q exists (M. Miyamoto,
A construction of Hadamard matrices, J. Combin. Theory Ser. A 57 (1991) 86–108).
The construction combines a symmetric conference matrix of order q + 1 with the
Hadamard matrix of order q − 1.

Here q = 3677. The inputs are

- Paley's symmetric conference matrix of order 3678: the matrix χ(y − x),
  x, y ∈ ℤ₃₆₇₇, with χ the quadratic character, bordered by a row and a column
  of ones;
- a Hadamard matrix K of order 3676 = 4 · 919. K is a skew-Hadamard matrix: a
  Goethals–Seidel matrix over ℤ₉₁₉ whose four blocks are constant on the cosets
  of the multiplier subgroup ⟨58⟩ of order 17. The first block is of skew type,
  the other three are related by the multiplier 52 (a spin structure), and the
  row sums are 1, 35, 35, 35. It came out of our skew-Hadamard searches.

The blocks are arranged as in SageMath's `hadamard_matrix_miyamoto_construction`.
`build_miyamoto` in [`verify/verify.py`](../verify/verify.py) spells out the
layout, including the reordering of the conference matrix.

`certificate.json` contains q, the choice of conference matrix, the complete
certificate of K and the SHA-256 of the result.

## Status

Probably new. In Table 1.53 of the Handbook of Combinatorial Designs (2nd ed.,
2007), n = 3677 has t = 4: the smallest known Hadamard matrix with odd part 3677 had
order 2^4 · 3677. The entry does not appear in Đoković's updates (arXiv:0912.5091,
arXiv:1008.2043), and Paley's constructions do not apply (4 · 3677 − 1 =
7 · 11 · 191 and 2 · 3677 − 1 = 3² · 19 · 43). The order has followed from
Miyamoto's theorem since Schumer published a skew-Hadamard matrix of order 3676 on
September 17, 2026, but we found no construction or statement of it; that record
lists only the orders 3676, 5436, 6644 and 7852. We read Table 1.53 from an
unofficial copy, so treat this claim with that caveat.

The order is above the range of Cati–Pasechnik Table 4, which stops at n = 2999.

The input order is not new. 3676 was unresolved in Cati–Pasechnik Table 4 until
Schumer's construction (Zenodo, September 17, 2026, doi:10.5281/zenodo.22820997),
which has the same spin structure. Ours was found independently.

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, conference matrix, certificate of K, SHA-256 |

The matrix is not stored. Its CSV encoding takes 541 MB (42 MB gzipped).

```
python3 verify/verify.py --write 14708
```

rebuilds and checks it and writes `14708/H14708.csv`. On a 48-core server this
takes about 20 seconds.

SHA-256 of `H14708.csv`:
`35d1bd3b6c2abfc33757943b617a8b37a8d7fa99f6f5dcc2c49ea5f9bbc83bf9`
