# Hadamard matrix of order 25972

25972 = 4 · 6493, with 6493 = 43 · 151.

## Construction

Cooper–Wallis. Four T-matrices of order 151, in [`t-matrices/T151.txt`](../../t-matrices/T151.txt), are circulant
matrices T₁, …, T₄ with entries 0, ±1 whose supports partition ℤ₁₅₁ and with Σ TᵢTᵢᵀ = 151 I. Putting them into the
Goethals–Seidel array in the usual way gives an orthogonal design OD(604; 151, 151, 151, 151). Substituting for its four
variables four Williamson matrices of order 43, symmetric circulant ±1 matrices W₁, …, W₄ with Σ Wₖ² = 172 I, gives a
Hadamard matrix of order 25972 (J. Cooper and J. Wallis, Bull. Austral. Math. Soc. 7 (1972)). The Wₖ are symmetric and
commute, so the cross terms of H Hᵀ cancel.

The Williamson matrices are those of SageMath's database (`williamson_type_quadruples_smallcases`, after Hall (1983)),
copied as data into [`search/data/williamson_sage.json`](../../search/data/williamson_sage.json). `certificate.json`
contains the first rows of the T-matrices (`0`, `+`, `-`) and of W₁, …, W₄ (`+`, `-`). `build_cooper_wallis` and
`assemble_cooper_wallis` in [`verify/verify.py`](../../verify/verify.py) spell out the assembly; they work in int8,
block by block.

## Search

None: both ingredients were known. `python3 search/cw_certify.py t-matrices/T151.txt 43 hadamard/25972` rebuilds the
matrix, checks it and writes the certificate.

## Status

**New, beyond the tables.** The Handbook (2007, Table 1.53) gives n = 6493 as t = 8: the smallest known Hadamard matrix
with odd part 6493 had order 2⁸ · 6493. The order lies above Cati–Pasechnik Table 4, which stops at n = 2999.

No construction we know of reaches it with inputs published before 2026:

- Paley I and II need 25971 = 3 · 11 · 787 and 12985 = 5 · 7² · 53 to be prime powers.
- Miyamoto's theorem does not apply: a symmetric conference matrix of order n + 1 = 6494 cannot exist (6493 is not a sum
  of two squares). Yamada–Kiyasu does not apply: n − 2 = 6491 is a prime ≡ 3 (mod 8), and the theorems need 1 or 5 (mod
  8).
- The entry is not among the Handbook entries that Đoković's updates settle (arXiv:0912.5091, arXiv:1008.2043).
- Since 6493 = 43 · 151, a product theorem needs T-sequences or T-matrices of order a and Williamson-type matrices of
  order b with a · b = 6493:
  - a = 43, b = 151: T-matrices of order 43 are classical, but Seberry–Yamada's Table A.1 (1992) marks 151 with t = 3
    (only 8-Williamson matrices), and no later list has Williamson-type matrices of order 151.
  - a = 151, b = 43: T-matrices of order 151 are not among the orders of Seberry–Yamada (1992, §3), and no base
    sequences give T-sequences of length 151: BS(76, 75) would need l = 75, beyond the known l ≤ 43 (Wang and Zhu,
    arXiv:2506.20296), and 3l − 1 = 151 has no integer solution. They appeared in 2026.
- The closure of the published theorems (Cooper–Wallis with Yang's and Golay products of base sequences, Miyamoto,
  Yamada, Kiyasu, conference matrices, Turyn's doubling) does not reach 4 · 6493 ([NOTES.md](../../NOTES.md) §28).
- SageMath's development version (October 17, 2026) does not construct it: `hadamard_matrix(25972, existence=True)`
  returns Unknown.

The route therefore uses an object from 2026: T-matrices of order 151 were posted by Schumer (Zenodo, September 17,
2026, doi:10.5281/zenodo.22820997). Lou, Xia and Xia's paper on T-matrices of order 6m + 1 (151 = 6 · 25 + 1), which we
could not read, might contain order 151; if it does, this order was obtainable classically. We found no construction or
statement of order 25972. The ingredients are not new, only their combination.

The order has followed from published objects since September 17, 2026, but we found it stated nowhere. Schumer's
deposit, which posted the T-matrices of order 151, says that Cooper–Wallis gives a Hadamard matrix of order 604w
whenever Williamson-type matrices of order w are available, but it builds only w = 9, 11 and 13 and names only values
with 151w ≤ 2999; w = 43 is not among them.

The sources searched are listed in [NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | T-matrices of order 151, Williamson matrices of order 43, SHA-256 |

The matrix is not stored. Its CSV encoding takes about 1.7 GB.

## Verification

```
python3 verify/verify.py --write 25972
```

SHA-256 of `H25972.csv`: `073135e9aadc51e204e2f471b9503fdb9cb9cde5f85027857b3d1f012feac4fa`
