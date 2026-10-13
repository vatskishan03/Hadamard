# Order 39260

39260 = 4 · 9815, with 9815 = 5 · 13 · 151.

## Construction

Cooper–Wallis with Williamson-type matrices that are not Williamson matrices. Four T-matrices of order 151, in
[`t-matrices/T151.txt`](../t-matrices/T151.txt), are circulant matrices T₁, …, T₄ with entries 0, ±1 whose supports
partition ℤ_151 and with Σ TᵢTᵢᵀ = 151 I. Putting them into the Goethals–Seidel array in the usual way gives an
orthogonal design OD(604; 151, 151, 151, 151). Substituting for its four variables four Williamson-type matrices of
order 65, that is four ±1 matrices W₁, …, W₄ with WₖWₗᵀ = WₗWₖᵀ and Σ WₖWₖᵀ = 260 I, gives a Hadamard matrix of
order 39260 (J. Cooper and J. Wallis, Bull. Austral. Math. Soc. 7 (1972); the amicability of the Wₖ is what cancels the
cross terms of H Hᵀ).

The Williamson-type matrices come from the near-Williamson matrices of order 65 in
[`near-williamson/NW65.txt`](../near-williamson/NW65.txt): circulant ±1 matrices A, B, C, D with B, C, D symmetric
and A Aᵀ + B Bᵀ + C Cᵀ + D Dᵀ = 260 I, taken from Table 6 of arXiv:2605.08661v2 (Kharaghani, Mohammadian, Tayfeh-Rezaie, Zaitsev). With R the back-identity, (AR, B, C, D) are
Williamson-type matrices: AR is symmetric (it is back-circulant) and commutes with the symmetric circulants B, C, D.
A itself is not symmetric. Williamson matrices of order 65 are not known.

`certificate.json` contains the first rows of the T-matrices (`0`, `+`, `-`) and of A, B, C, D (`+`, `-`).
`build_cooper_wallis_near_williamson` in [`verify/verify.py`](../verify/verify.py) spells out the assembly; it works in
int8, block by block, because a Kronecker product in int64 would not fit in memory.

## Status

Probably new. In Table 1.53 of the Handbook of Combinatorial Designs (2nd ed., 2007), n = 9815 has t = 3: the
smallest known Hadamard matrix with odd part 9815 had order 2^3 · 9815. The entry does not appear in
Đoković's updates (arXiv:0912.5091, arXiv:1008.2043), and no classical family gives order 4 · 9815: Paley's first construction would need 4 · 9815 − 1 = 39259 = 11 · 43 · 83 and the second 2 · 9815 − 1 = 19629 = 3³ · 727 to be prime powers; 9815 is not a prime power (Miyamoto), and 9815 − 2 = 9813 = 3 · 3271 is not one either (Yamada, Kiyasu); neither
the Mathon, Seberry–Whiteman, Spence, Whiteman nor Leung–Momihara–Xiang families contain it. Since 9815 = 5 · 13 · 151, a product theorem needs T-sequences or T-matrices of length t and Williamson-type matrices of order w with t · w = 9815. With t = 65, w = 151: T-sequences of length 65 are known, but Williamson-type matrices of order 151 are not in Table A.1 of Seberry and Yamada (1992) or in Đoković's lists. With t = 151, w = 65: T-sequences of length 151 are not known (only periodic T-matrices, since September 2026), and Williamson-type matrices of order 65 appeared only in 2026. The other divisors would need T-sequences of length 1963 or 755, which are not known.
Every route therefore uses an object from 2026. T-matrices of order 151 were first made public by Schumer (Zenodo, September 17, 2026, doi:10.5281/zenodo.22820997). Near-Williamson matrices of order 65
are from Table 6 of arXiv:2605.08661v2 (Kharaghani, Mohammadian, Tayfeh-Rezaie, Zaitsev). That paper lists only quaternary Hadamard matrices (of orders 94, 118, 130, 134, 206 and 266) as
applications, and neither JD-Jones-ASES/Hadamard-T nor Schumer's deposit pairs T-matrices with anything but Williamson matrices of
order at most 29. We found no construction or statement of order 39260. We read Table 1.53 from an unofficial copy, so
treat this claim with that caveat.

The order is above the range of Cati–Pasechnik Table 4, which stops at n = 2999.

The ingredients are not new, only their combination.

## Files

| File | Contents |
|---|---|
| `certificate.json` | T-matrices of order 151, near-Williamson matrices of order 65, SHA-256 |

The matrix is not stored. Its CSV encoding takes 3.9 GB (201 MB gzipped).

```
python3 verify/verify.py --write 39260
```

rebuilds and checks it and writes `39260/H39260.csv`. On this machine (a 10-core laptop shared with other jobs) this takes 3 minutes 25 seconds and 2.8 GB of memory.

SHA-256 of `H39260.csv`:
`cbb53ae69f72073e0ef3a2103509711b3ffb73385a04d26dfb152aaaac8a07a1`
