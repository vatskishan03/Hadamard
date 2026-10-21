# Hadamard matrix of order 24308

24308 = 4 · 6077, with 6077 = 59 · 103.

## Construction

Cooper–Wallis with Williamson-type matrices that are not Williamson matrices. Four T-matrices of order 103, in
[`t-matrices/T103.txt`](../../t-matrices/T103.txt), are circulant matrices T₁, …, T₄ with entries 0, ±1 whose supports
partition ℤ₁₀₃ and with Σ TᵢTᵢᵀ = 103 I. Putting them into the Goethals–Seidel array in the usual way gives an
orthogonal design OD(412; 103, 103, 103, 103). Substituting for its four variables four Williamson-type matrices of
order 59, that is four ±1 matrices W₁, …, W₄ with WₖWₗᵀ = WₗWₖᵀ and Σ WₖWₖᵀ = 236 I, gives a Hadamard matrix of order
24308 (J. Cooper and J. Wallis, Bull. Austral. Math. Soc. 7 (1972); the amicability of the Wₖ is what cancels the
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

**New, beyond the tables.** The Handbook (2007, Table 1.53) gives n = 6077 as t = 4: the smallest known Hadamard
matrix with odd part 6077 had order 2⁴ · 6077. The order lies above Cati–Pasechnik Table 4, which stops at n = 2999.

No construction we know of reaches 24308:

- Paley I and II need 24307 = 109 · 223 and 12153 = 3 · 4051 to be prime powers.
- Miyamoto's theorem does not apply: a symmetric conference matrix of order n + 1 = 6078 cannot exist (6077 is not a sum
  of two squares). Neither is 6077 − 2 = 6075 = 3⁵ · 5² (Yamada, Kiyasu).
- The Mathon, Seberry–Whiteman, Spence, Whiteman and Leung–Momihara–Xiang families do not contain it.
- Since 6077 = 59 · 103, a product theorem needs T-sequences or T-matrices of length t and Williamson-type matrices of
  order w with t · w = 6077:
  - t = 59, w = 103: T-sequences of length 59 are known, but Williamson-type matrices of order 103 are not in
    Table A.1 of Seberry and Yamada (1992) or in Đoković's lists.
  - t = 103, w = 59: T-sequences of length 103 are not known (they are known up to length 100 except possibly 97),
    and Williamson-type matrices of order 59 appeared only in 2026.

Every route therefore uses an object from 2026:

- T-matrices of order 103 were decoded from the Hadamard matrix of order 2060 that Schneider posted publicly on
  23 August 2026 (JD-Jones-ASES/Hadamard-T, 29 August 2026). Ours were found independently, see
  [order 2060](../2060).
- Near-Williamson matrices of order 59 are from Table 5 of Kharaghani, Mohammadian, Tayfeh-Rezaie and Zaitsev,
  arXiv:2605.08661v2 (September 2026). That paper lists only quaternary Hadamard matrices (of orders 94, 118, 130,
  134, 206 and 266) as applications.

Neither JD-Jones-ASES/Hadamard-T nor Schumer's deposit pairs T-matrices with anything but Williamson matrices of
order at most 29. We found no construction or statement of order 24308. The ingredients are not new, only their
combination. The sources searched are listed in [NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | T-matrices of order 103, near-Williamson matrices of order 59, SHA-256 |

The matrix is not stored. Its CSV encoding takes 1.5 GB (76 MB gzipped).

## Verification

```
python3 verify/verify.py --write 24308
```

This rebuilds and checks the matrix and writes `hadamard/24308/H24308.csv`. It takes 1 minute 20 seconds and 1.4 GB of
memory.

SHA-256 of `H24308.csv`: `bf9d8a5578a4cf241d9475401cf51d7993f4bed19d3013089e7d59002f72f29e`
