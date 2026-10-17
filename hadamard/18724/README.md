# Hadamard matrix of order 18724

18724 = 4 · 4681, with 4681 = 31 · 151.

## Construction

Cooper–Wallis. Four T-matrices of order 151, in [`t-matrices/T151.txt`](../../t-matrices/T151.txt), are circulant
matrices T₁, …, T₄ with entries 0, ±1 whose supports partition ℤ₁₅₁ and with Σ TᵢTᵢᵀ = 151 I. Putting them into the
Goethals–Seidel array in the usual way gives an orthogonal design OD(604; 151, 151, 151, 151). Substituting for its four
variables four Williamson matrices of order 31, symmetric circulant ±1 matrices W₁, …, W₄ with Σ Wₖ² = 124 I, gives a
Hadamard matrix of order 18724 (J. Cooper and J. Wallis, Bull. Austral. Math. Soc. 7 (1972)). The Wₖ are symmetric and
commute, so the cross terms of H Hᵀ cancel.

The Williamson matrices are those of SageMath's database (`williamson_type_quadruples_smallcases`, after London (2013)),
copied as data into [`search/data/williamson_sage.json`](../../search/data/williamson_sage.json). `certificate.json`
contains the first rows of the T-matrices (`0`, `+`, `-`) and of W₁, …, W₄ (`+`, `-`). `build_cooper_wallis` and
`assemble_cooper_wallis` in [`verify/verify.py`](../../verify/verify.py) spell out the assembly; they work in int8,
block by block.

## Search

None: both ingredients were known. `python3 search/cw_certify.py t-matrices/T151.txt 31 hadamard/18724` rebuilds the
matrix, checks it and writes the certificate.

## Status

**New, beyond the tables.** The Handbook (2007, Table 1.53) gives n = 4681 as t = 3: the smallest known Hadamard matrix
with odd part 4681 had order 2³ · 4681. The order lies above Cati–Pasechnik Table 4, which stops at n = 2999.

No construction we know of reaches it with inputs published before 2026:

- Paley I and II need 18723 = 3 · 79² and 9361 = 11 · 23 · 37 to be prime powers.
- Miyamoto's theorem does not apply: a symmetric conference matrix of order n + 1 = 4682 cannot exist (4681 is not a sum
  of two squares). Yamada–Kiyasu does not apply: n − 2 = 4679 is a prime ≡ 7 (mod 8), and the theorems need 1 or 5 (mod
  8).
- The entry is not among the Handbook entries that Đoković's updates settle (arXiv:0912.5091, arXiv:1008.2043).
- Since 4681 = 31 · 151, a product theorem needs T-sequences or T-matrices of order a and Williamson-type matrices of
  order b with a · b = 4681:
  - a = 31, b = 151: T-matrices of order 31 are classical, but Seberry–Yamada's Table A.1 (1992) marks 151 with t = 3
    (only 8-Williamson matrices), and no later list has Williamson-type matrices of order 151.
  - a = 151, b = 31: T-matrices of order 151 are not among the orders of Seberry–Yamada (1992, §3), and no base
    sequences give T-sequences of length 151: BS(76, 75) would need l = 75, beyond the known l ≤ 43 (Wang and Zhu,
    arXiv:2506.20296), and 3l − 1 = 151 has no integer solution. They appeared in 2026.
- The closure of the published theorems (Cooper–Wallis with Yang's and Golay products of base sequences, Miyamoto,
  Yamada, Kiyasu, conference matrices, Turyn's doubling) does not reach 4 · 4681 ([NOTES.md](../../NOTES.md) §28).
- SageMath's development version (October 17, 2026) does not construct it: `hadamard_matrix(18724, existence=True)`
  returns Unknown.

The route therefore uses an object from 2026: T-matrices of order 151 were posted by Schumer (Zenodo, September 17,
2026, doi:10.5281/zenodo.22820997). Lou, Xia and Xia's paper on T-matrices of order 6m + 1 (151 = 6 · 25 + 1), which we
could not read, might contain order 151; if it does, this order was obtainable classically. We found no construction or
statement of order 18724. The ingredients are not new, only their combination. The sources searched are listed in
[NOVELTY.md](../../NOVELTY.md).

The order has followed from published objects since September 17, 2026, but we found it stated nowhere. Schumer's
deposit, which posted the T-matrices of order 151, says that Cooper–Wallis gives a Hadamard matrix of order 604w
whenever Williamson-type matrices of order w are available, but it builds only w = 9, 11 and 13 and names only values
with 151w ≤ 2999; w = 31 is not among them.

## Files

| File | Contents |
|---|---|
| `certificate.json` | T-matrices of order 151, Williamson matrices of order 31, SHA-256 |

The matrix is not stored. Its CSV encoding takes about 0.9 GB.

## Verification

```
python3 verify/verify.py --write 18724
```

SHA-256 of `H18724.csv`: `8e17f33c0783dc3b5956740fbb044e967a2923877e44e0dacbeb1d2eec706aa7`
