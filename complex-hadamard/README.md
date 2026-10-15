# Complex Hadamard matrices of orders 2N, N odd

A quaternary (complex) Hadamard matrix of order n has entries 1, −1, i, −i and satisfies K K* = n I. Its order is 1
or even, and it is conjectured that one exists for every even order (Turyn). A real Hadamard matrix is a quaternary
one, and real Hadamard matrices are known for every multiple of 4 below 2000 since August 2026 ([README](../README.md)).
The open cases below 2000 are therefore orders 2N with N odd.

Recent progress at the bottom of that range:
- Szöllősi constructed order 94 (arXiv:2603.09572, March 2026), noting that "the next outstanding order seems to be"
  118.
- Kharaghani, Mohammadian, Tayfeh-Rezaie and Zaitsev (arXiv:2605.08661v2, September 2026) gave the first examples of
  orders 118, 130, 134, 206 and 266, from near-Williamson matrices.
- V. Zaitsev's web page [Quaternary Hadamard matrices](https://www.cs.uleth.ca/~zaitsevv/quaternary.html) (last
  updated 14 October 2026) gives matrices of 117 orders that are blank in the 1992 table, 23 of them "of the
  Szollosi type".

This directory gives 22 orders 2N that are blank in the 1992 table. We found them independently of Zaitsev's page,
which we learned of on 15 October 2026, after they were committed here; the page also has all 22. See
[Status](#status).

## Construction

Szöllősi, Theorem 4 of arXiv:2603.09572. Let A, B, C, D be circulant ±1 matrices of odd order N with

    A Aᵀ + B Bᵀ + C Cᵀ + D Dᵀ = 4N I,          A and B symmetric,

and let R be the back-diagonal identity. Then

    K = [ (A+B)/2 + i(A−B)/2           (C+DR)/2 + i(C−DR)/2 ]
        [ (Cᵀ+DR)/2 − i(Cᵀ−DR)/2      −(A+B)/2 + i(A−B)/2  ]

is a quaternary Hadamard matrix of order 2N.

**Proof sketch.** Put X = (A+B)/2 + i(A−B)/2 and Y = (C+DR)/2 + i(C−DR)/2. Their entries lie in {±1, ±i}. X is
symmetric and circulant, so it commutes with C, D and R, and K = [[X, Y], [Y*, −X*]].
- Since A, B and C, DR are amicable pairs, X X* + Y Y* = (A Aᵀ + B Bᵀ + C Cᵀ + D Dᵀ)/2 = 2N I.
- X and Y are normal, so X* X + Y* Y = 2N I as well.
- XY = YX makes the two block rows orthogonal.

Williamson and near-Williamson matrices are the special cases with three or four symmetric blocks. Allowing two
arbitrary circulant blocks enlarges the search space by far.

## Search

All four blocks are constant on the orbits of a multiplier group H ≤ ℤ_N* with −1 ∉ H, so C and D are H-invariant.
A and B are symmetric, hence constant on the orbits of ⟨H, −1⟩. The condition Σ PAF = 0 off zero then holds on the
⟨H, −1⟩-classes of ℤ_N \ {0}.

For each row-sum pattern (a, b, c, d) with a² + b² + c² + d² = 4N:
- the four blocks are listed completely (`search/bnw_enum.c`);
- the lists are joined exactly (`search/join4x.c`, pairing each symmetric list with an arbitrary one);
- every hit is turned into K and checked exactly, and the real Goethals–Seidel matrix of order 4N is checked as well.

Code: `search/ssxx.py` (triage, search), `search/ssxx_sweep.py` (the candidate orders and the sweep),
`search/ssxx_sample.py` (for structures too large to list: uniform samples of the large lists, joined exactly; for
prime N the Gaussian estimate of the hits inside the samples is accurate: 1445 expected and 1430 found at N = 71,
10^4.0 and 12618 at N = 131), `search/ssxx_group.py` (the same array over the non-cyclic groups GF(p^k) × ℤ_m and ℤ_p²
× ℤ_m, with diagonal field multiplications or odd-order subgroups of GL(2, p): the theorem holds verbatim for
group-developed matrices over any abelian group, and `verify/verify_complex.py` checks such certificates),
`search/ssxx_certify.py` (certificates), `search/cchain4.py` and `search/cchain_sweep.py` (structures whose complete
lists are too large to join: an exact compression chain over the subgroups of ℤ_N, as for order 3964; order 310 came
from ⟨16⟩ on ℤ₁₅₅ with lists of 10⁴, 10⁴, 7 · 10⁸ and 6 · 10⁸, compressed to ℤ₅ and split on ℤ₃₁).

Structures whose blocks C and D have 10⁸–10¹¹ members are searched by the low-energy k-tree `search/lowe_sz.py`: each
big block is enumerated in full and only its members of low energy (the most likely to complete a solution under a
Gaussian model of the other three blocks) are kept, and `search/kt4.c` joins the kept lists (orders 342, 366, 426, 430
and 454). The structures are ranked by `search/gauss_mc.py`, a hybrid estimate of the number of solutions: for composite
N (or N − 1) the per-class Gaussian estimate of `ssxx.py` can be low by several orders of magnitude, and 430 and 454
came from structures it put at 10^−1.2 and 10^0.6 solutions (the hybrid estimate: 10^4.1 and 10^4.4). See
[search/NOTES.md](../search/NOTES.md) §32–34.

Two exact necessary conditions now come first. For each prime power d dividing the group order, the four blocks'
compressions to ℤ_d must join exactly (`search/compress_check.py`). On 291 structures with exact counts this rules out
244 of the 264 without solutions and none of the 27 with solutions; joining pairs of coprime levels together
(`--joint`) rules out 15 more. Each block then keeps only members whose compressions occur
in some level-1 solution (fibre filters). For order 366 the filters left 2.1% of C, all of which was kept, and the
search found 1536 solutions where earlier runs without them had found none. See
[search/NOTES.md](../search/NOTES.md) §36–39.

A second construction, the bordered Szöllősi array on a group of order N − 1, is stated and proved in
[BORDERED.md](BORDERED.md). Two twisted specialisations of both arrays (B = A∘μ, D = C∘μ, and a half twist), which
turn the four-sequence search into a two-list search, are in [TWISTED.md](TWISTED.md).

## Matrices

| order 2N | N | multiplier group / source | row sums (a, b; c, d) | file |
|---:|---|---|---|---|
| 142 | 71 | ⟨5⟩, order 5 | 9, 9; 1, 11 | [C142.json](C142.json) |
| 154 | 77 = 7 · 11 | ⟨15⟩, order 5 | 3, 3; 1, 17 | [C154.json](C154.json) |
| 222 | 111 = 3 · 37 | near-Williamson matrices of order 111 ([NW111.txt](../near-williamson/NW111.txt)), reordered (B, C, A, D) | 1, 1; 21, 1 | [C222.json](C222.json) |
| 238 | 119 = 7 · 17 | ⟨8⟩, order 8 | 1, 5; 15, 15 | [C238.json](C238.json) |
| 262 | 131 | ⟨53⟩, order 5 | 9, 9; 1, 19 | [C262.json](C262.json) |
| 286 | 143 = 11 · 13 | ⟨3⟩, order 15 | 3, 9; 11, 19 | [C286.json](C286.json) |
| 302 | 151 | ⟨8⟩, order 5 | 9, 9; 9, 19 | [C302.json](C302.json) |
| 310 | 155 = 5 · 31 | ⟨16⟩, order 5 (compression chain, `search/cchain4.py`) | 9, 9; 13, 17 | [C310.json](C310.json) |
| 326 | 163 | ⟨38⟩, order 9 | 17, 19; 1, 1 | [C326.json](C326.json) |
| 342 | 171 = 9 · 19 | ⟨11⟩, order 6 (low-energy k-tree, `search/lowe_sz.py`) | 1, 11; 11, 21 | [C342.json](C342.json) |
| 366 | 183 = 3 · 61 | bordered array over ℤ₁₈₂ ([BORDERED.md](BORDERED.md)), ⟨3⟩, order 6 (low-energy k-tree with fibre filters) | 0, 0; 2, 0 | [C366.json](C366.json) |
| 370 | 185 = 5 · 37 | ⟨9⟩, order 18 | 3, 3; 19, 19 | [C370.json](C370.json) |
| 382 | 191 | ⟨39⟩, order 5 (streaming k-tree, `search/kstream.py`) | 11, 11; 9, 21 | [C382.json](C382.json) |
| 426 | 213 = 3 · 71 | ⟨37⟩, order 7 (low-energy k-tree) | 13, 17; 13, 15 | [C426.json](C426.json) |
| 430 | 215 = 5 · 43 | ⟨11⟩, order 7 (low-energy k-tree) | 11, 17; 15, 15 | [C430.json](C430.json) |
| 454 | 227 | bordered array over ℤ₂₂₆ ([BORDERED.md](BORDERED.md)), ⟨49⟩, order 7 (low-energy k-tree) | 0, 0; 2, 0 | [C454.json](C454.json) |
| 494 | 247 = 13 · 19 | ⟨9, 35⟩, order 27 | 17, 19; 7, 17 | [C494.json](C494.json) |
| 518 | 259 = 7 · 37 | ⟨44⟩, order 9 (compression chain) | 5, 29; 1, 13 | [C518.json](C518.json) |
| 534 | 267 = 3 · 89 | ⟨4⟩, order 11 | 1, 1; 21, 25 | [C534.json](C534.json) |
| 742 | 371 = 7 · 53 | ⟨13⟩, order 26 | 3, 3; 25, 29 | [C742.json](C742.json) |
| 754 | 377 = 13 · 29 | ⟨9⟩, order 42 | 27, 27; 5, 5 | [C754.json](C754.json) |
| 1558 | 779 = 19 · 41 | ⟨49⟩, order 60 | 5, 7; 39, 39 | [C1558.json](C1558.json) |

Each certificate gives N, the multiplier group, the row sums, the first rows of A, B, C, D (`+`, `-`) and the
SHA-256 of the canonical text of K: 2N lines of 2N digits, digit e standing for the entry iᵉ.

```
python3 verify/verify_complex.py            # all certificates
python3 verify/verify_complex.py --write 142 # also writes complex-hadamard/C142.txt
```

This rebuilds K from the four sequences alone and checks every entry and K K* = 2N I exactly, in integer arithmetic
on the real and imaginary parts.

## Status

**Primary source.** Table 11.2 of J. Seberry and M. Yamada, "Hadamard matrices, sequences, and block designs"
(Contemporary Design Theory, Wiley, 1992) records, for every odd q < 1000, whether a complex Hadamard matrix of
order 2q was known. Its constructions are conference matrices (codes c1–c6) and Williamson-type matrices (Turyn;
Kharaghani–Seberry, codes x2 and x6). Every order in this directory is **blank** in that table. The parsed table is
[`search/data/complex_known_1992.txt`](../search/data/complex_known_1992.txt).

**Results since 1992** that we checked, none of which gives these orders:
- Đoković's good matrices of orders 35 and 127 (1993): 70 and 254.
- Đoković's Williamson-type matrices of orders 509, 857 and 883 (arXiv:1008.2043, 2010, from theorems of Miyamoto
  and Yamada): 1018, 1714 and 1766. Our first status list missed these three.
- Szöllősi (2026): 94.
- Kharaghani, Mohammadian, Tayfeh-Rezaie and Zaitsev (2026): 118, 130, 134, 206 and 266.
- Quaternary Legendre pairs, giving order 2ℓ + 2 (Kotsireas–Winterhof 2022; Kotsireas–Koutschan–Winterhof 2025;
  Jedwab–Pender 2025; Lebedev and Verma 2026): every even length ℓ up to 100 except 70, 76, 88, 92, 94 and 100,
  plus two infinite families. Order 142 needs ℓ = 70 and 154 needs ℓ = 76, both open. Our larger orders fall
  outside both families.
- Mathon's conference matrices (q²(q ± 2) + 1). We also count them as known where the 1992 table is blank
  (N = 221, 787).

The sweep's candidate list (`python3 search/ssxx_sweep.py list`) is the blank entries of the 1992 table minus these
(and, since 15 October 2026, minus the orders on Zaitsev's page). Szöllősi's paper and arXiv:2605.08661 also take
their open orders from Table 11.2 of the 1992 survey. We did not consult the tables of Seberry and Yamada's 2020 book
(Wiley), so a construction published between 1992 and 2020 and recorded only there is the remaining caveat for the
earlier literature. At q = 163 our copy of the table is blank; Zaitsev's page calls the entry ambiguous and treats
the order as unknown, because the construction claimed there would need 2q − 1 = 325 to be a prime power.

**Concurrent work.** V. Zaitsev's page (University of Lethbridge; "Last updated: 14 October 2026"; read on
15 October 2026) gives matrices of 117 orders blank in the 1992 table: 5 from arXiv:2605.08661, 23 "of the Szollosi
type, with two of the four blocks symmetric", 27 by "unpublished approaches" and 62 from Paley-type matrices. All 22
orders of this directory are among them: 21 in its Szöllősi-type list and 454 among the unpublished approaches.
- Our certificates were first committed on 14 October 2026 (142 and 154 at 05:24 UTC, the last, 366, at 20:09 UTC).
  We do not know when each order first appeared on the page.
- We downloaded the page's matrices of our 22 orders and checked K K* = nI exactly for each.
- They are different matrices. Where the page's matrix has the layout of our array (orders 142, 262, 302, 326 and
  382), none of its four sequences has the PAF spectrum of any of ours (the spectrum is unchanged by shifts,
  multipliers, negation and reversal).
- So these 22 orders are not new to the literature. What this directory adds is a second, independent construction
  for each, with a certificate that `verify/verify_complex.py` checks from the four sequences alone.

**What remains open** (`python3 search/complex_status.py`): 115 orders 2N below 2000 with N odd have no quaternary
Hadamard matrix known to us. This is exactly the list of unknown orders on Zaitsev's page. The smallest are 334, 358,
418, 438, 446, 478 and 498. Order 334 (N = 167) is out of reach of every structure we have: ℤ₁₆₇ has no usable
multiplier group, and on ℤ₁₆₆ = ℤ₂ × ℤ₈₃ the only one has order 41 and six orbits. Order 358 (N = 179) is similar: the
cyclic array has no structure in our range of orbit counts, and the bordered array's best has an estimate of 10⁻⁵
solutions. For 418 the best structures have estimates of one to three solutions in pools of 10¹¹–10¹⁵ members; for 438,
446, 478 and 498 every estimate is below one. The one open order where the array has strong structures is 682 (N = 341
= 11 · 31): eight multiplier groups of order 10 with estimates of 10^1.7–10^3.0 solutions. The two best joins
found no solution in 97% and 95% of their collision mass, which the estimate gave a chance of 2–5% each
([search/NOTES.md](../search/NOTES.md) §39–41).
