# Hadamard matrix of order 39260

39260 = 4 · 9815, with 9815 = 5 · 13 · 151.

## Construction

Cooper–Wallis with Williamson-type matrices that are not Williamson matrices. Four T-matrices of order 151, in
[`t-matrices/T151.txt`](../../t-matrices/T151.txt), are circulant matrices T₁, …, T₄ with entries 0, ±1 whose supports
partition ℤ₁₅₁ and with Σ TᵢTᵢᵀ = 151 I. Putting them into the Goethals–Seidel array in the usual way gives an
orthogonal design OD(604; 151, 151, 151, 151). Substituting for its four variables four Williamson-type matrices of
order 65, that is four ±1 matrices W₁, …, W₄ with WₖWₗᵀ = WₗWₖᵀ and Σ WₖWₖᵀ = 260 I, gives a Hadamard matrix of order
39260 (J. Cooper and J. Wallis, Bull. Austral. Math. Soc. 7 (1972); the amicability of the Wₖ is what cancels the
cross terms of H Hᵀ).

The Williamson-type matrices come from the near-Williamson matrices of order 65 in
[`near-williamson/NW65.txt`](../../near-williamson/NW65.txt): circulant ±1 matrices A, B, C, D with B, C, D symmetric
and A Aᵀ + B Bᵀ + C Cᵀ + D Dᵀ = 260 I, taken from Table 6 of arXiv:2605.08661v2 (Kharaghani, Mohammadian,
Tayfeh-Rezaie, Zaitsev). With R the back-identity, (AR, B, C, D) are Williamson-type matrices: AR is symmetric (it is
back-circulant) and commutes with the symmetric circulants B, C, D. A itself is not symmetric. Williamson matrices of
order 65 are not known.

`certificate.json` contains the first rows of the T-matrices (`0`, `+`, `-`) and of A, B, C, D (`+`, `-`).
`build_cooper_wallis_near_williamson` in [`verify/verify.py`](../../verify/verify.py) spells out the assembly; it
works in int8, block by block, because a Kronecker product in int64 would not fit in memory.

## Status

**New, beyond the tables.** The Handbook (2007, Table 1.53) gives n = 9815 as t = 3: the smallest known Hadamard
matrix with odd part 9815 had order 2³ · 9815. The order lies above Cati–Pasechnik Table 4, which stops at n = 2999.

No construction we know of reaches 39260:

- Paley I and II need 39259 = 11 · 43 · 83 and 19629 = 3³ · 727 to be prime powers.
- Miyamoto does not apply: 9815 is not a prime power. Neither is 9815 − 2 = 9813 = 3 · 3271 (Yamada, Kiyasu).
- The Mathon, Seberry–Whiteman, Spence, Whiteman and Leung–Momihara–Xiang families do not contain it.
- Since 9815 = 5 · 13 · 151, a product theorem needs T-sequences or T-matrices of length t and Williamson-type
  matrices of order w with t · w = 9815:
  - t = 65, w = 151: T-sequences of length 65 are known, but Williamson-type matrices of order 151 are not in
    Table A.1 of Seberry and Yamada (1992) or in Đoković's lists.
  - t = 151, w = 65: T-sequences of length 151 are not known (only periodic T-matrices, since September 2026), and
    Williamson-type matrices of order 65 appeared only in 2026.
  - The other divisors would need T-sequences of length 1963 or 755, which are not known.

Every route therefore uses an object from 2026:

- T-matrices of order 151 were first made public by Schumer (Zenodo, September 17, 2026, doi:10.5281/zenodo.22820997).
- Near-Williamson matrices of order 65 are from Table 6 of arXiv:2605.08661v2 (Kharaghani, Mohammadian,
  Tayfeh-Rezaie, Zaitsev). That paper lists only quaternary Hadamard matrices (of orders 94, 118, 130, 134, 206 and
  266) as applications.

Neither JD-Jones-ASES/Hadamard-T nor Schumer's deposit pairs T-matrices with anything but Williamson matrices of
order at most 29. We found no construction or statement of order 39260. The ingredients are not new, only their
combination. The sources searched are listed in [NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | T-matrices of order 151, near-Williamson matrices of order 65, SHA-256 |

The matrix is not stored. Its CSV encoding takes 3.9 GB (201 MB gzipped).

## Verification

```
python3 verify/verify.py --write 39260
```

This rebuilds and checks the matrix and writes `hadamard/39260/H39260.csv`. It takes 3 minutes 25 seconds and 2.8 GB
of memory.

SHA-256 of `H39260.csv`: `cbb53ae69f72073e0ef3a2103509711b3ffb73385a04d26dfb152aaaac8a07a1`
