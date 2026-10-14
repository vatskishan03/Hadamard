# Order 27604

27604 = 4 · 6901, with 6901 = 67 · 103.

## Construction

Cooper–Wallis with Williamson-type matrices that are not Williamson matrices. Four T-matrices of order 103, in
[`t-matrices/T103.txt`](../../t-matrices/T103.txt), are circulant matrices T₁, …, T₄ with entries 0, ±1 whose supports
partition ℤ_103 and with Σ TᵢTᵢᵀ = 103 I. Putting them into the Goethals–Seidel array in the usual way gives an
orthogonal design OD(412; 103, 103, 103, 103). Substituting for its four variables four Williamson-type matrices of
order 67, that is four ±1 matrices W₁, …, W₄ with WₖWₗᵀ = WₗWₖᵀ and Σ WₖWₖᵀ = 268 I, gives a Hadamard matrix of
order 27604 (J. Cooper and J. Wallis, Bull. Austral. Math. Soc. 7 (1972); the amicability of the Wₖ is what cancels the
cross terms of H Hᵀ).

The Williamson-type matrices come from the near-Williamson matrices of order 67 in
[`near-williamson/NW67.txt`](../../near-williamson/NW67.txt): circulant ±1 matrices A, B, C, D with B, C, D symmetric
and A Aᵀ + B Bᵀ + C Cᵀ + D Dᵀ = 268 I, taken from Table 7 of arXiv:2605.08661v2 (Kharaghani, Mohammadian, Tayfeh-Rezaie, Zaitsev). With R the back-identity, (AR, B, C, D) are
Williamson-type matrices: AR is symmetric (it is back-circulant) and commutes with the symmetric circulants B, C, D.
A itself is not symmetric. Williamson matrices of order 67 are not known.

`certificate.json` contains the first rows of the T-matrices (`0`, `+`, `-`) and of A, B, C, D (`+`, `-`).
`build_cooper_wallis_near_williamson` in [`verify/verify.py`](../../verify/verify.py) spells out the assembly; it works in
int8, block by block, because a Kronecker product in int64 would not fit in memory.

## Status

Probably new. In Table 1.53 of the Handbook of Combinatorial Designs (2nd ed., 2007), n = 6901 has t = 3: the
smallest known Hadamard matrix with odd part 6901 had order 2^3 · 6901. The entry does not appear in
Đoković's updates (arXiv:0912.5091, arXiv:1008.2043), and no classical family gives order 4 · 6901: Paley's first construction would need 4 · 6901 − 1 = 27603 = 3² · 3067 and the second 2 · 6901 − 1 = 13801 = 37 · 373 to be prime powers; 6901 is not a prime power (Miyamoto), and 6901 − 2 = 6899 is prime but ≡ 3 (mod 4), so the theorems of Yamada and Kiyasu do not apply; neither
the Mathon, Seberry–Whiteman, Spence, Whiteman nor Leung–Momihara–Xiang families contain it. Since 6901 = 67 · 103, a product theorem needs T-sequences or T-matrices of length t and Williamson-type matrices of order w with t · w = 6901. With t = 67, w = 103: T-sequences of length 67 are known, but Williamson-type matrices of order 103 are not in Table A.1 of Seberry and Yamada (1992) or in Đoković's lists. With t = 103, w = 67: T-sequences of length 103 are not known, and Williamson-type matrices of order 67 are absent from Table A.1 and appeared only in 2026.
Every route therefore uses an object from 2026. T-matrices of order 103 were decoded from the Hadamard matrix of order 2060 that Schneider posted publicly on August 23, 2026 (JD-Jones-ASES/Hadamard-T, August 29, 2026); ours were found independently, see [order 2060](../2060). Near-Williamson matrices of order 67
are from Table 7 of arXiv:2605.08661v2 (Kharaghani, Mohammadian, Tayfeh-Rezaie, Zaitsev). That paper lists only quaternary Hadamard matrices (of orders 94, 118, 130, 134, 206 and 266) as
applications, and neither JD-Jones-ASES/Hadamard-T nor Schumer's deposit pairs T-matrices with anything but Williamson matrices of
order at most 29. We found no construction or statement of order 27604. We read Table 1.53 from an unofficial copy, so
treat this claim with that caveat.

The order is above the range of Cati–Pasechnik Table 4, which stops at n = 2999.

The ingredients are not new, only their combination.

## Files

| File | Contents |
|---|---|
| `certificate.json` | T-matrices of order 103, near-Williamson matrices of order 67, SHA-256 |

The matrix is not stored. Its CSV encoding takes 1.9 GB (107 MB gzipped).

```
python3 verify/verify.py --write 27604
```

rebuilds and checks it and writes `hadamard/27604/H27604.csv`. On this machine (a 10-core laptop shared with other jobs) this takes 2 minutes 16 seconds and 1.5 GB of memory.

SHA-256 of `H27604.csv`:
`29182c796b25084283d3b39bb6fcbdb4e72e9aacfb9b9877b03a6fd45e393252`
