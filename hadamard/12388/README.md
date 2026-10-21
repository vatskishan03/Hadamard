# Hadamard matrix of order 12388

12388 = 4 · 3097, with 3097 = 19 · 163.

## Construction

Cooper–Wallis. Four T-matrices of order 163, in [`t-matrices/T163.txt`](../../t-matrices/T163.txt), are circulant
matrices T₁, …, T₄ with entries 0, ±1 whose supports partition ℤ₁₆₃ and with Σ TᵢTᵢᵀ = 163 I. Putting them into the
Goethals–Seidel array in the usual way gives an orthogonal design OD(652; 163, 163, 163, 163). Substituting for its four
variables four Williamson matrices of order 19, symmetric circulant ±1 matrices W₁, …, W₄ with Σ Wₖ² = 76 I, gives a
Hadamard matrix of order 12388 (J. Cooper and J. Wallis, Bull. Austral. Math. Soc. 7 (1972)). The Wₖ are symmetric and
commute, so the cross terms of H Hᵀ cancel.

The Williamson matrices are those of SageMath's database (`williamson_type_quadruples_smallcases`, after Hall (1983)),
copied as data into [`search/data/williamson_sage.json`](../../search/data/williamson_sage.json). `certificate.json`
contains the first rows of the T-matrices (`0`, `+`, `-`) and of W₁, …, W₄ (`+`, `-`). `build_cooper_wallis` and
`assemble_cooper_wallis` in [`verify/verify.py`](../../verify/verify.py) spell out the assembly; they work in int8,
block by block.

## Search

None: both ingredients were known. `python3 search/cw_certify.py t-matrices/T163.txt 19 hadamard/12388` rebuilds the
matrix, checks it and writes the certificate.

## Status

**New, beyond the tables.** The Handbook (2007, Table 1.53) gives n = 3097 as t = 3: the smallest known Hadamard matrix
with odd part 3097 had order 2³ · 3097. The order lies above Cati–Pasechnik Table 4, which stops at n = 2999.

No construction we know of reaches it with inputs published before 2026:

- Paley I and II need 12387 = 3 · 4129 and 6193 = 11 · 563 to be prime powers.
- Miyamoto's theorem does not apply: a symmetric conference matrix of order n + 1 = 3098 cannot exist (3097 is not a sum
  of two squares). Yamada–Kiyasu does not apply: n − 2 = 3095 = 5 · 619 is not a prime power.
- The entry is not among the Handbook entries that Đoković's updates settle (arXiv:0912.5091, arXiv:1008.2043).
- Since 3097 = 19 · 163, a product theorem needs T-sequences or T-matrices of order a and Williamson-type matrices of
  order b with a · b = 3097:
  - a = 19, b = 163: T-matrices of order 19 are classical, but Seberry–Yamada's Table A.1 (1992) marks 163 with t = 3
    (only 8-Williamson matrices), and no later list has Williamson-type matrices of order 163.
  - a = 163, b = 19: T-matrices of order 163 are not among the orders of Seberry–Yamada (1992, §3), and no base
    sequences give T-sequences of length 163: BS(82, 81) would need l = 81, beyond the known l ≤ 43 (Wang and Zhu,
    arXiv:2506.20296), and 3l − 1 = 163 has no integer solution. They appeared in 2026.
- The closure of the published theorems (Cooper–Wallis with Yang's and Golay products of base sequences, Miyamoto,
  Yamada, Kiyasu, conference matrices, Turyn's doubling) does not reach 4 · 3097 ([NOTES.md](../../NOTES.md) §28).
- SageMath's development version (October 17, 2026) does not construct it: `hadamard_matrix(12388, existence=True)`
  returns Unknown.

The route therefore uses an object from 2026: T-matrices of order 163 were posted by JD-Jones-ASES/Hadamard-T on August
29, 2026. Lou, Xia and Xia's paper on T-matrices of order 6m + 1 (163 = 6 · 27 + 1), which we could not read, might
contain order 163; if it does, this order was obtainable classically. We found no construction or statement of order
12388. The ingredients are not new, only their combination.

The order has followed from published objects since August 29, 2026, but we found it stated nowhere.
JD-Jones-ASES/Hadamard-T pairs its T-matrices of order 163 only with Williamson matrices of orders 11 and 13 (orders
7172 and 8476). JD-Jones-ASES/Hadamard-formal (August 30, 2026) holds kernel-checked T-matrices of order 163, Williamson
matrices of orders 19 and 23 and a generic Cooper–Wallis theorem, so 12388 is one instantiation away there, but its
results state only eight orders, from 2884 to 11948. Schumer's deposit (September 17, 2026) contains an independent
T(163) as a control, with no product.

The sources searched are listed in [NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | T-matrices of order 163, Williamson matrices of order 19, SHA-256 |

The matrix is not stored. Its CSV encoding takes about 0.4 GB.

## Verification

```
python3 verify/verify.py --write 12388
```

SHA-256 of `H12388.csv`: `7fa0cc0e75c3be4e2d4e0b7d19cf40b6bddd2082368beec00aaceeea00f9acce`
