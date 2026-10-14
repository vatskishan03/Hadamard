# Order 38468

38468 = 4 · 9617, with 9617 = 59 · 163.

## Construction

Cooper–Wallis with Williamson-type matrices that are not Williamson matrices. Four T-matrices of order 163, in
[`t-matrices/T163.txt`](../../t-matrices/T163.txt), are circulant matrices T₁, …, T₄ with entries 0, ±1 whose supports
partition ℤ_163 and with Σ TᵢTᵢᵀ = 163 I. Putting them into the Goethals–Seidel array in the usual way gives an
orthogonal design OD(652; 163, 163, 163, 163). Substituting for its four variables four Williamson-type matrices of
order 59, that is four ±1 matrices W₁, …, W₄ with WₖWₗᵀ = WₗWₖᵀ and Σ WₖWₖᵀ = 236 I, gives a Hadamard matrix of
order 38468 (J. Cooper and J. Wallis, Bull. Austral. Math. Soc. 7 (1972); the amicability of the Wₖ is what cancels the
cross terms of H Hᵀ).

The Williamson-type matrices come from the near-Williamson matrices of order 59 in
[`near-williamson/NW59.txt`](../../near-williamson/NW59.txt): circulant ±1 matrices A, B, C, D with B, C, D symmetric
and A Aᵀ + B Bᵀ + C Cᵀ + D Dᵀ = 236 I, taken from Table 5 of H. Kharaghani, A. Mohammadian, B. Tayfeh-Rezaie and V. Zaitsev, arXiv:2605.08661v2 (September 2026). With R the back-identity, (AR, B, C, D) are
Williamson-type matrices: AR is symmetric (it is back-circulant) and commutes with the symmetric circulants B, C, D.
A itself is not symmetric. No Williamson matrices of order 59 exist (Holzmann, Kharaghani and Tayfeh-Rezaie 2008).

`certificate.json` contains the first rows of the T-matrices (`0`, `+`, `-`) and of A, B, C, D (`+`, `-`).
`build_cooper_wallis_near_williamson` in [`verify/verify.py`](../../verify/verify.py) spells out the assembly; it works in
int8, block by block, because a Kronecker product in int64 would not fit in memory.

## Status

Probably new. In Table 1.53 of the Handbook of Combinatorial Designs (2nd ed., 2007), n = 9617 has t = 4: the
smallest known Hadamard matrix with odd part 9617 had order 2^4 · 9617. The entry does not appear in
Đoković's updates (arXiv:0912.5091, arXiv:1008.2043), and no classical family gives order 4 · 9617: Paley's first construction would need 4 · 9617 − 1 = 38467 = 11 · 13 · 269 and the second 2 · 9617 − 1 = 19233 = 3² · 2137 to be prime powers; 9617 is not a prime power (Miyamoto), and 9617 − 2 = 9615 = 3 · 5 · 641 is not one either (Yamada, Kiyasu); neither
the Mathon, Seberry–Whiteman, Spence, Whiteman nor Leung–Momihara–Xiang families contain it. Since 9617 = 59 · 163, a product theorem needs T-sequences or T-matrices of length t and Williamson-type matrices of order w with t · w = 9617. With t = 59, w = 163: T-sequences of length 59 are known, but Williamson-type matrices of order 163 are not in Table A.1 of Seberry and Yamada (1992) or in Đoković's lists. With t = 163, w = 59: T-sequences of length 163 are not known (only periodic T-matrices, since August 2026), and Williamson-type matrices of order 59 appeared only in 2026.
Every route therefore uses an object from 2026. T-matrices of order 163 were first made public in JD-Jones-ASES/Hadamard-T (August 29, 2026). Near-Williamson matrices of order 59
are from Table 5 of H. Kharaghani, A. Mohammadian, B. Tayfeh-Rezaie and V. Zaitsev, arXiv:2605.08661v2 (September 2026). That paper lists only quaternary Hadamard matrices (of orders 94, 118, 130, 134, 206 and 266) as
applications, and neither JD-Jones-ASES/Hadamard-T nor Schumer's deposit pairs T-matrices with anything but Williamson matrices of
order at most 29. We found no construction or statement of order 38468. We read Table 1.53 from an unofficial copy, so
treat this claim with that caveat.

The order is above the range of Cati–Pasechnik Table 4, which stops at n = 2999.

The ingredients are not new, only their combination.

## Files

| File | Contents |
|---|---|
| `certificate.json` | T-matrices of order 163, near-Williamson matrices of order 59, SHA-256 |

The matrix is not stored. Its CSV encoding takes 3.7 GB (199 MB gzipped).

```
python3 verify/verify.py --write 38468
```

rebuilds and checks it and writes `hadamard/38468/H38468.csv`. On this machine (a 10-core laptop shared with other jobs) this takes 5 minutes 13 seconds and 2.7 GB of memory.

SHA-256 of `H38468.csv`:
`7efe5e75738d1e14e0f2891ec066a4c58d81ba1ba90bf249ad7c146877c6f071`
