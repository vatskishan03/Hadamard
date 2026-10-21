# Hadamard matrix of order 20116

20116 = 4 · 5029, with 5029 = 47 · 107.

## Construction

Cooper–Wallis with Williamson-type matrices that are not Williamson matrices. Four T-matrices of order 107, in
[`t-matrices/T107.txt`](../../t-matrices/T107.txt), give an orthogonal design OD(428; 107, 107, 107, 107) in the
Goethals–Seidel array. They are T-sequences of length 107 obtained from the Turyn-type sequences TT(36) of Kharaghani
and Tayfeh-Rezaie (2005), through the base sequences BS(71, 36), as SageMath's `t_sequences` module stores and builds
them (`T_sequences_smallcases(107)`). Substituting for the four variables Williamson-type matrices of order 47, that is
four ±1 matrices W₁, …, W₄ with WₖWₗᵀ = WₗWₖᵀ and Σ WₖWₖᵀ = 188 I, gives a Hadamard matrix of order 20116 (J. Cooper and
J. Wallis, Bull. Austral. Math. Soc. 7 (1972)). The amicability of the Wₖ is what cancels the cross terms of H Hᵀ.

The Williamson-type matrices of order 47 are the bordered quadruple of
[`williamson-type/WT47.json`](../../williamson-type/WT47.json), which also gives [12596](../12596) and
[19364](../19364). `certificate.json` contains the T-matrices and the four sequences on ℤ₄₆;
`build_cooper_wallis_bordered` in [`verify/verify.py`](../../verify/verify.py) rebuilds the bordered matrices and checks
them.

## Search

None: both ingredients were known. `python3 search/cw_certify.py t-matrices/T107.txt williamson-type/WT47.json
hadamard/20116` rebuilds the matrix, checks it and writes the certificate.

## Status

**New, beyond the tables.** The Handbook (2007, Table 1.53) gives n = 5029 as t = 3: the smallest known Hadamard matrix
with odd part 5029 had order 2³ · 5029. The order lies above Cati–Pasechnik Table 4, which stops at n = 2999.

No construction we know of reaches it with inputs published before 2026:

- Paley I and II need 20115 = 3³ · 5 · 149 and 10057 = 89 · 113 to be prime powers.
- Miyamoto's theorem does not apply: a symmetric conference matrix of order n + 1 = 5030 cannot exist (5029 is not a sum
  of two squares). Yamada–Kiyasu does not apply: n − 2 = 5027 = 11 · 457 is not a prime power.
- The entry is not among the Handbook entries that Đoković's updates settle (arXiv:0912.5091, arXiv:1008.2043).
- Since 5029 = 47 · 107, a product theorem needs T-sequences or T-matrices of order a and Williamson-type matrices of
  order b with a · b = 5029:
  - a = 47, b = 107: T-matrices of order 47 are classical, but Seberry–Yamada's Table A.1 (1992) marks 107 with t = 3
    (only 8-Williamson matrices), and no later list has Williamson-type matrices of order 107.
  - a = 107, b = 47: T-sequences of length 107 are classical (2005), but Williamson matrices of order 47 do not exist
    (Holzmann, Kharaghani and Tayfeh-Rezaie 2008), and Williamson-type matrices of order 47 appeared only in 2026
    (near-Williamson matrices of Kharaghani, Mohammadian, Tayfeh-Rezaie and Zaitsev, arXiv:2605.08661).
- The closure of the published theorems (Cooper–Wallis with Yang's and Golay products of base sequences, Miyamoto,
  Yamada, Kiyasu, conference matrices, Turyn's doubling) does not reach 4 · 5029 ([NOTES.md](../../NOTES.md) §28).
- SageMath's development version (October 17, 2026) does not construct it: `hadamard_matrix(20116, existence=True)`
  returns Unknown.

The route therefore uses an object from 2026: the Williamson-type matrices of order 47. We found no construction or
statement of order 20116. The ingredients are not new, only their combination.

The order has followed from published objects since May 9, 2026 (the near-Williamson matrices of order 47 in the first
version of arXiv:2605.08661, with the T-sequences of length 107 of 2005), but we found it stated nowhere: that paper
applies its matrices to quaternary Hadamard matrices only, and SageMath's database has no Williamson-type matrices of
order 47.

The sources searched are listed in [NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | T-matrices of order 107, the bordered Williamson-type sequences of order 47, SHA-256 |

The matrix is not stored. Its CSV encoding takes about 1.0 GB.

## Verification

```
python3 verify/verify.py --write 20116
```

SHA-256 of `H20116.csv`: `8c02bff1b8fd297131d20c5aa51b6e4872f2ffbe9bbf03ab3dfec434b7129615`
