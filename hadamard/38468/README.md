# Hadamard matrix of order 38468

38468 = 4 · 9617, with 9617 = 59 · 163.

## Construction

Cooper–Wallis with Williamson-type matrices that are not Williamson matrices. Four T-matrices of order 163, in
[`t-matrices/T163.txt`](../../t-matrices/T163.txt), are circulant matrices T₁, …, T₄ with entries 0, ±1 whose supports
partition ℤ₁₆₃ and with Σ TᵢTᵢᵀ = 163 I. Putting them into the Goethals–Seidel array in the usual way gives an
orthogonal design OD(652; 163, 163, 163, 163). Substituting for its four variables four Williamson-type matrices of
order 59, that is four ±1 matrices W₁, …, W₄ with WₖWₗᵀ = WₗWₖᵀ and Σ WₖWₖᵀ = 236 I, gives a Hadamard matrix of order
38468 (J. Cooper and J. Wallis, Bull. Austral. Math. Soc. 7 (1972); the amicability of the Wₖ is what cancels the
cross terms of H Hᵀ).

The Williamson-type matrices come from the near-Williamson matrices of order 59 in
[`near-williamson/NW59.txt`](../../near-williamson/NW59.txt): circulant ±1 matrices A, B, C, D with B, C, D symmetric
and A Aᵀ + B Bᵀ + C Cᵀ + D Dᵀ = 236 I, taken from Table 5 of H. Kharaghani, A. Mohammadian, B. Tayfeh-Rezaie and
V. Zaitsev, arXiv:2605.08661v2 (September 2026). With R the back-identity, (AR, B, C, D) are Williamson-type matrices:
AR is symmetric (it is back-circulant) and commutes with the symmetric circulants B, C, D. A itself is not symmetric.
No Williamson matrices of order 59 exist (Holzmann, Kharaghani and Tayfeh-Rezaie 2008).

`certificate.json` contains the first rows of the T-matrices (`0`, `+`, `-`) and of A, B, C, D (`+`, `-`).
`build_cooper_wallis_near_williamson` in [`verify/verify.py`](../../verify/verify.py) spells out the assembly; it
works in int8, block by block, because a Kronecker product in int64 would not fit in memory.

## Status

**New, beyond the tables.** The Handbook (2007, Table 1.53) gives n = 9617 as t = 4: the smallest known Hadamard
matrix with odd part 9617 had order 2⁴ · 9617. The order lies above Cati–Pasechnik Table 4, which stops at n = 2999.

No construction we know of reaches 38468:

- Paley I and II need 38467 = 11 · 13 · 269 and 19233 = 3² · 2137 to be prime powers.
- Miyamoto's theorem does not apply: a symmetric conference matrix of order n + 1 = 9618 cannot exist (9617 is not a sum
  of two squares). Neither is 9617 − 2 = 9615 = 3 · 5 · 641 (Yamada, Kiyasu).
- The Mathon, Seberry–Whiteman, Spence, Whiteman and Leung–Momihara–Xiang families do not contain it.
- Since 9617 = 59 · 163, a product theorem needs T-sequences or T-matrices of length t and Williamson-type matrices of
  order w with t · w = 9617:
  - t = 59, w = 163: T-sequences of length 59 are known, but Williamson-type matrices of order 163 are not in
    Table A.1 of Seberry and Yamada (1992) or in Đoković's lists.
  - t = 163, w = 59: T-sequences of length 163 are not known (only periodic T-matrices, since August 2026), and
    Williamson-type matrices of order 59 appeared only in 2026.

Every route therefore uses an object from 2026:

- T-matrices of order 163 were first made public in JD-Jones-ASES/Hadamard-T (August 29, 2026).
- Near-Williamson matrices of order 59 are from Table 5 of Kharaghani, Mohammadian, Tayfeh-Rezaie and Zaitsev,
  arXiv:2605.08661v2 (September 2026). That paper lists only quaternary Hadamard matrices (of orders 94, 118, 130,
  134, 206 and 266) as applications.

Neither JD-Jones-ASES/Hadamard-T nor Schumer's deposit pairs T-matrices with anything but Williamson matrices of
order at most 29. We found no construction or statement of order 38468. The ingredients are not new, only their
combination. The sources searched are listed in [NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | T-matrices of order 163, near-Williamson matrices of order 59, SHA-256 |

The matrix is not stored. Its CSV encoding takes 3.7 GB (199 MB gzipped).

## Verification

```
python3 verify/verify.py --write 38468
```

This rebuilds and checks the matrix and writes `hadamard/38468/H38468.csv`. It takes 5 minutes 13 seconds and 2.7 GB
of memory.

SHA-256 of `H38468.csv`: `7efe5e75738d1e14e0f2891ec066a4c58d81ba1ba90bf249ad7c146877c6f071`
