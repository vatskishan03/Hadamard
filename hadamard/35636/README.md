# Hadamard matrix of order 35636

35636 = 4 · 8909, with 8909 = 59 · 151.

## Construction

Cooper–Wallis with Williamson-type matrices that are not Williamson matrices. Four T-matrices of order 151, in
[`t-matrices/T151.txt`](../../t-matrices/T151.txt), are circulant matrices T₁, …, T₄ with entries 0, ±1 whose supports
partition ℤ₁₅₁ and with Σ TᵢTᵢᵀ = 151 I. Putting them into the Goethals–Seidel array in the usual way gives an
orthogonal design OD(604; 151, 151, 151, 151). Substituting for its four variables four Williamson-type matrices of
order 59, that is four ±1 matrices W₁, …, W₄ with WₖWₗᵀ = WₗWₖᵀ and Σ WₖWₖᵀ = 236 I, gives a Hadamard matrix of order
35636 (J. Cooper and J. Wallis, Bull. Austral. Math. Soc. 7 (1972); the amicability of the Wₖ is what cancels the
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

**New, beyond the tables.** The Handbook (2007, Table 1.53) gives n = 8909 as t = 3: the smallest known Hadamard
matrix with odd part 8909 had order 2³ · 8909. The order lies above Cati–Pasechnik Table 4, which stops at n = 2999.

No construction we know of reaches 35636:

- Paley I and II need 35635 = 5 · 7127 and 17817 = 3 · 5939 to be prime powers.
- Miyamoto's theorem does not apply: a symmetric conference matrix of order n + 1 = 8910 cannot exist (8909 is not a sum
  of two squares). Neither is 8909 − 2 = 8907 = 3 · 2969 (Yamada, Kiyasu).
- The Mathon, Seberry–Whiteman, Spence, Whiteman and Leung–Momihara–Xiang families do not contain it.
- Since 8909 = 59 · 151, a product theorem needs T-sequences or T-matrices of length t and Williamson-type matrices of
  order w with t · w = 8909:
  - t = 59, w = 151: T-sequences of length 59 are known, but Williamson-type matrices of order 151 are not in
    Table A.1 of Seberry and Yamada (1992) or in Đoković's lists.
  - t = 151, w = 59: T-sequences of length 151 are not known (only periodic T-matrices, since September 2026), and
    Williamson-type matrices of order 59 appeared only in 2026.

Every route therefore uses an object from 2026:

- T-matrices of order 151 were first made public by Schumer (Zenodo, September 17, 2026, doi:10.5281/zenodo.22820997).
- Near-Williamson matrices of order 59 are from Table 5 of Kharaghani, Mohammadian, Tayfeh-Rezaie and Zaitsev,
  arXiv:2605.08661v2 (September 2026). That paper lists only quaternary Hadamard matrices (of orders 94, 118, 130,
  134, 206 and 266) as applications.

Neither JD-Jones-ASES/Hadamard-T nor Schumer's deposit pairs T-matrices with anything but Williamson matrices of
order at most 29. We found no construction or statement of order 35636. The ingredients are not new, only their
combination. The sources searched are listed in [NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | T-matrices of order 151, near-Williamson matrices of order 59, SHA-256 |

The matrix is not stored. Its CSV encoding takes 3.2 GB (161 MB gzipped).

## Verification

```
python3 verify/verify.py --write 35636
```

This rebuilds and checks the matrix and writes `hadamard/35636/H35636.csv`. It takes 3 minutes 44 seconds and 2.3 GB
of memory.

SHA-256 of `H35636.csv`: `33f584fcbfa4896f49384aa13e858d5216f377f7ebe92cfa115ee6e0634cc77a`
