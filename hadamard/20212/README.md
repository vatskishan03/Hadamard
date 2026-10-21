# Hadamard matrix of order 20212

20212 = 4 · 5053, with 5053 = 31 · 163.

## Construction

Cooper–Wallis. Four T-matrices of order 163, in [`t-matrices/T163.txt`](../../t-matrices/T163.txt), are circulant
matrices T₁, …, T₄ with entries 0, ±1 whose supports partition ℤ₁₆₃ and with Σ TᵢTᵢᵀ = 163 I. Putting them into the
Goethals–Seidel array in the usual way gives an orthogonal design OD(652; 163, 163, 163, 163). Substituting for its four
variables four Williamson matrices of order 31, symmetric circulant ±1 matrices W₁, …, W₄ with Σ Wₖ² = 124 I, gives a
Hadamard matrix of order 20212 (J. Cooper and J. Wallis, Bull. Austral. Math. Soc. 7 (1972)). The Wₖ are symmetric and
commute, so the cross terms of H Hᵀ cancel.

The Williamson matrices are those of SageMath's database (`williamson_type_quadruples_smallcases`, after London (2013)),
copied as data into [`search/data/williamson_sage.json`](../../search/data/williamson_sage.json). `certificate.json`
contains the first rows of the T-matrices (`0`, `+`, `-`) and of W₁, …, W₄ (`+`, `-`). `build_cooper_wallis` and
`assemble_cooper_wallis` in [`verify/verify.py`](../../verify/verify.py) spell out the assembly; they work in int8,
block by block.

## Search

None: both ingredients were known. `python3 search/cw_certify.py t-matrices/T163.txt 31 hadamard/20212` rebuilds the
matrix, checks it and writes the certificate.

## Status

**New, beyond the tables.** The Handbook (2007, Table 1.53) gives n = 5053 as t = 3: the smallest known Hadamard matrix
with odd part 5053 had order 2³ · 5053. The order lies above Cati–Pasechnik Table 4, which stops at n = 2999.

No construction we know of reaches it with inputs published before 2026:

- Paley I and II need 20211 = 3 · 6737 and 10105 = 5 · 43 · 47 to be prime powers.
- Miyamoto's theorem does not apply: a symmetric conference matrix of order n + 1 = 5054 cannot exist (5053 is not a sum
  of two squares). Yamada–Kiyasu does not apply: n − 2 = 5051 is a prime ≡ 3 (mod 8), and the theorems need 1 or 5 (mod
  8).
- The entry is not among the Handbook entries that Đoković's updates settle (arXiv:0912.5091, arXiv:1008.2043).
- Since 5053 = 31 · 163, a product theorem needs T-sequences or T-matrices of order a and Williamson-type matrices of
  order b with a · b = 5053:
  - a = 31, b = 163: T-matrices of order 31 are classical, but Seberry–Yamada's Table A.1 (1992) marks 163 with t = 3
    (only 8-Williamson matrices), and no later list has Williamson-type matrices of order 163.
  - a = 163, b = 31: T-matrices of order 163 are not among the orders of Seberry–Yamada (1992, §3), and no base
    sequences give T-sequences of length 163: BS(82, 81) would need l = 81, beyond the known l ≤ 43 (Wang and Zhu,
    arXiv:2506.20296), and 3l − 1 = 163 has no integer solution. They appeared in 2026.
- The closure of the published theorems (Cooper–Wallis with Yang's and Golay products of base sequences, Miyamoto,
  Yamada, Kiyasu, conference matrices, Turyn's doubling) does not reach 4 · 5053 ([NOTES.md](../../NOTES.md) §28).
- SageMath's development version (17 October 2026) does not construct it: `hadamard_matrix(20212, existence=True)`
  returns Unknown.

The route therefore uses an object from 2026: T-matrices of order 163 were posted by JD-Jones-ASES/Hadamard-T on 29
August 2026. Lou, Xia and Xia's paper on T-matrices of order 6m + 1 (163 = 6 · 27 + 1), which we could not read, might
contain order 163; if it does, this order was obtainable classically. We found no construction or statement of order
20212. The ingredients are not new, only their combination.

The order has followed from published objects since 29 August 2026, but we found it stated nowhere.
JD-Jones-ASES/Hadamard-T pairs its T-matrices of order 163 only with Williamson matrices of orders 11 and 13 (orders
7172 and 8476). JD-Jones-ASES/Hadamard-formal (30 August 2026) has a generic Cooper–Wallis theorem and T-matrices of
order 163, but no Williamson matrices of order 31, and states only eight orders, from 2884 to 11948. Schumer's deposit
(17 September 2026) contains an independent T(163) as a control, with no product.

The sources searched are listed in [NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | T-matrices of order 163, Williamson matrices of order 31, SHA-256 |

The matrix is not stored. Its CSV encoding takes about 1.0 GB.

## Verification

```
python3 verify/verify.py --write 20212
```

SHA-256 of `H20212.csv`: `451a7a8f51a07958d2a477e9afe5178ed7bb2aba5c85fdd32f5b552cce7a3641`
