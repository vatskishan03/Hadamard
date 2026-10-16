# Complex Hadamard matrices of orders 2N, N odd

A quaternary (complex) Hadamard matrix of order n has entries 1, −1, i, −i and satisfies K K* = n I. Its order is 1
or even, and it is conjectured that one exists for every even order (Turyn). Orders divisible by 4 are covered by real
Hadamard matrices. The open cases are therefore the orders 2N with N odd.

Recent progress at the bottom of that range:
- Szöllősi constructed order 94 (arXiv:2603.09572, March 2026), noting that "the next outstanding order seems to be"
  118.
- Kharaghani, Mohammadian, Tayfeh-Rezaie and Zaitsev (arXiv:2605.08661v2, September 2026) gave the first examples of
  orders 118, 130, 134, 206 and 266, from near-Williamson matrices.
- V. Zaitsev's page [Quaternary Hadamard matrices](https://www.cs.uleth.ca/~zaitsevv/quaternary.html) (last updated
  14 October 2026; work with the same coauthors) gives explicit matrices of 122 orders 2N, most of them from Szöllősi's
  array with two symmetric blocks, the construction used here.

This directory constructs 11 orders that had no construction in the published work we checked: 214, 238, 302, 326, 366,
382, 426, 494, 518, 534 and 754. All 11 are also on Zaitsev's page, found concurrently and independently (we saw the
page on 16 October 2026, after our searches); we claim no priority for them. It also gives explicit matrices of 76
orders that are blank in the 1992 table but follow from a result we had missed: Paley's matrix of order q + 1, q ≡ 3
(mod 4) a prime power, makes the dicyclic group of order 2(q + 1) a Hadamard group (Ito), and Egan (Des. Codes Cryptogr.
87 (2019) 67–74, Theorem 4.5) turns this into a complex Hadamard matrix of order (q + 1)/2. So every order 2N with N odd
and 4N − 1 a prime power is known. We found the construction independently, with a direct proof ([PALEY.md](PALEY.md)),
and had at first labelled 13 of these orders New; they are not.

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
from ⟨16⟩ on ℤ₁₅₅ with lists of 10⁴, 10⁴, 7 · 10⁸ and 6 · 10⁸, compressed to ℤ₅ and split on ℤ₃₁),
`search/ssxx_paley.py` (Paley products N = m · q with H the squares of GF(q) acting trivially on ℤ_m: every block is
u ⊗ δ + w ⊗ J + z ⊗ χ, and the search is exact over the w-positions and two joins over ℤ_m; order 330 came from
ℤ₁₅ × GF(11), [search/NOTES.md](../search/NOTES.md) §41), `search/lowe_sz.py` (structures whose arbitrary blocks
have 10⁸–10¹¹ members: every member is visited, those of low energy under the Gaussian model are kept, and the kept
lists are joined by a k-tree; order 426 came from ⟨37⟩ on ℤ₂₁₃ with lists of 6435, 6435, 3 · 10⁷ and 3 · 10⁷).
`search/comp_est.py`, `search/comp_triage.py` and `search/lowe_comp.py` condition the estimates and the low-energy
lists on the block compressions to a quotient ℤ_d: a row pattern with no feasible compression has no solution, and
the live patterns have far more than the plain Gaussian says ([search/NOTES.md](../search/NOTES.md) §42).
Order 366 came from ⟨34⟩ on ℤ₁₈₃: its pattern 1, 23; 9, 11 lives on a single class of ℤ₃ compressions
(A ≡ (−1, 1, 1), B ≡ (1, 11, 11), C ≡ (−1, −1, 11), D ≡ (−11, 11, 11)), which the unconditioned low-energy lists
had missed. `search/comp_scan.py plain` runs this conditioning over every structure ℤ_d × ℤ_q with H ≤ ℤ_q* fixing
ℤ_d, for all open N < 1000; only three structures stay alive, and the first, N = 215 = 5 · 43 with ⟨41⟩ of order 7,
gave order 430: its pattern 15, 15; 11, 17 has a compression class with pools of only 1200, 1200, 2.4 · 10⁶ and
2.4 · 10⁶ and about 200 predicted solutions, and the exact join of the complete lists found 5 in nine minutes. The
same conditioning over every multiplier group of ℤ_N* and every listable quotient (`comp_scan.py plainall`) gave order
342: N = 171 = 9 · 19 with ⟨11⟩ of order 6, which permutes the units of ℤ₉ transitively and has order 3 on ℤ₁₉, rows
11, 11; 9, 19, conditioned on ℤ₉ (13 feasible tuples); the first tuple (pools 1134, 378, 2.8 · 10⁷, 2.8 · 10⁷) gave 12
solutions in four minutes. (A ℤ₃ estimate alone can mislead: for ⟨31⟩ the best ℤ₃ tuple is dead on ℤ₉.)

A second construction, the bordered Szöllősi array on a group of order N − 1, is stated and proved in
[BORDERED.md](BORDERED.md). Order 454 comes from it: 227 is a safe prime, so ℤ₂₂₇ has no usable multiplier group,
while ℤ₂₂₆ = ℤ₂ × ℤ₁₁₃ has ⟨49⟩ of order 7. Its fibre sums over ℤ₂ are ≡ ±1 (mod 14), which the Gaussian estimate
ignores (10^0.6); conditioned on the compressions to ℤ₂ (`search/comp_est.py`) the estimate is 10^6.0, and the
low-energy k-tree on one compression class (pools 4900, 4900, 1.7 · 10⁸, 1.3 · 10⁸, lists of 2 · 10⁷) found 64
solutions in six minutes ([search/NOTES.md](../search/NOTES.md) §42).

A third construction, the character-bordered Szöllősi array, is stated and proved in [TBORDERED.md](TBORDERED.md).
The core lives on a group G of order N − t with a quotient φ: G → Q, the PAF sums are prescribed by a set T of t
characters of Q (Σ_k PAF_k(g) = −4 Σ_{χ∈T} χ(φ(g)) off 0), and a border of size 2t completes it; t = 1 is the
bordered array. With Q = ℤ₄ and T = {1, χ, χ̄} it applies to every N ≡ 3 (mod 4) through a group of order N − 3.
Order 214 comes from it: N = 107 is a safe prime and 53 = (107 − 1)/2 is prime, so neither ℤ₁₀₇ nor ℤ₁₀₆ has a
usable multiplier group, but ℤ₁₀₄ = ℤ₈ × ℤ₁₃ with H of order 3 on ℤ₁₃ (40 orbits) rates 10^5.7, and one compression
class on ℤ₄ (pools 5040, 2520, 2 · 10⁷ and 1.9 · 10⁸, about 10^3.8 predicted solutions) gave 5 solutions in six minutes
(`search/tbsz.py`). With Q = ℤ₆ and T all characters but the sign (t = 5) it applies to every N ≡ 5 (mod 6); its
border is found by an exhaustive search (`search/border_exh.c`). Order 322 comes from it: on ℤ₂ × ℤ₂ × ℤ₃ × ℤ₁₃
with H of order 6 (a shear on ℤ₂², −1 on ℤ₃, order 3 on ℤ₁₃; 35 orbits) the pattern 2, 0; 4, 0 rates 10^2.5 on ℤ₆,
and one compression class with complete lists (864, 432, 8.2 · 10⁶ and 1.3 · 10⁶) gave solutions after 80 seconds of
the join.

Two twisted specialisations of both arrays (B = A∘μ, D = C∘μ, and a half twist), which turn the four-sequence search
into a two-list search, are in [TWISTED.md](TWISTED.md).

## Matrices

**Constructed here and concurrently by Zaitsev et al.** (Concurrent). For none of them is 4N − 1 or 2N − 1 a prime
power, and none is covered by the published sources we checked (Status section). Zaitsev's page has matrices of all
11. For 10 of them ours is inequivalent to theirs: the multiset over row pairs of the counts
#{c : K[a,c] K̄[b,c] = iᵏ}, k = 0..3, taken up to rotation and reversal, is invariant under row and column
permutations, unit scalings, conjugation and transposition, and it differs. For 238 it agrees, and so does a finer
row-wise version; the two may be equivalent ([NOTES.md](../NOTES.md) §20).

| order 2N | N | multiplier group / source | row sums (a, b; c, d) | file |
|---:|---|---|---|---|
| 214 | 107 | character-bordered array ([TBORDERED.md](TBORDERED.md)) on ℤ₁₀₄ = ℤ₈ × ℤ₁₃, Q = ℤ₄, T = {1, χ, χ̄}, ⟨(1, 3)⟩ of order 3 (low-energy k-tree on one compression class, `search/tbsz.py`) | 2, 0; 2, 2 | [C214.json](C214.json) |
| 238 | 119 = 7 · 17 | ⟨8⟩, order 8 | 1, 5; 15, 15 | [C238.json](C238.json) |
| 302 | 151 | ⟨8⟩, order 5 | 9, 9; 9, 19 | [C302.json](C302.json) |
| 326 | 163 | ⟨38⟩, order 9 | 17, 19; 1, 1 | [C326.json](C326.json) |
| 366 | 183 = 3 · 61 | ⟨34⟩, order 5 (low-energy k-tree on one compression class, `search/lowe_comp.py`) | 1, 23; 9, 11 | [C366.json](C366.json) |
| 382 | 191 | ⟨39⟩, order 5 (streaming k-tree, `search/kstream.py`) | 11, 11; 9, 21 | [C382.json](C382.json) |
| 426 | 213 = 3 · 71 | ⟨37⟩, order 7 (low-energy k-tree, `search/lowe_sz.py`) | 13, 17; 13, 15 | [C426.json](C426.json) |
| 494 | 247 = 13 · 19 | ⟨9, 35⟩, order 27 | 17, 19; 7, 17 | [C494.json](C494.json) |
| 518 | 259 = 7 · 37 | ⟨44⟩, order 9 (compression chain) | 5, 29; 1, 13 | [C518.json](C518.json) |
| 534 | 267 = 3 · 89 | ⟨4⟩, order 11 | 1, 1; 21, 25 | [C534.json](C534.json) |
| 754 | 377 = 13 · 29 | ⟨9⟩, order 42 | 27, 27; 5, 5 | [C754.json](C754.json) |

**Known orders, independent matrices** (Known: Ito; Egan 2019, since 4N − 1 is a prime power). We found 13 of them
(142, 154, 222, 262, 286, 310, 322, 330, 342, 370, 430, 454, 742) by the searches before we knew this; the other 63
are halves of Paley's matrix.

| order 2N | N | multiplier group / source | row sums (a, b; c, d) | file |
|---:|---|---|---|---|
| 142 | 71 | ⟨5⟩, order 5 | 9, 9; 1, 11 | [C142.json](C142.json) |
| 154 | 77 = 7 · 11 | ⟨15⟩, order 5 | 3, 3; 1, 17 | [C154.json](C154.json) |
| 222 | 111 = 3 · 37 | near-Williamson matrices of order 111 ([NW111.txt](../near-williamson/NW111.txt)), reordered (B, C, A, D) | 1, 1; 21, 1 | [C222.json](C222.json) |
| 262 | 131 | ⟨53⟩, order 5 | 9, 9; 1, 19 | [C262.json](C262.json) |
| 286 | 143 = 11 · 13 | ⟨3⟩, order 15 | 3, 9; 11, 19 | [C286.json](C286.json) |
| 310 | 155 = 5 · 31 | ⟨16⟩, order 5 (compression chain, `search/cchain4.py`) | 9, 9; 13, 17 | [C310.json](C310.json) |
| 322 | 161 = 7 · 23 | character-bordered array ([TBORDERED.md](TBORDERED.md)) on ℤ₂ × ℤ₂ × ℤ₃ × ℤ₁₃, Q = ℤ₆, T = all characters but the sign, H of order 6 (exact k-tree on one compression class, `search/tbsz.py`) | 2, 0; 4, 0 | [C322.json](C322.json) |
| 330 | 165 = 3 · 5 · 11 | ⟨136⟩, order 5 (Paley product ℤ₁₅ × GF(11), `search/ssxx_paley.py`) | 7, 9; 1, 23 | [C330.json](C330.json) |
| 342 | 171 = 9 · 19 | ⟨11⟩, order 6 (exact k-tree on one ℤ₉ compression class, `search/lowe_comp.py`) | 11, 11; 9, 19 | [C342.json](C342.json) |
| 370 | 185 = 5 · 37 | ⟨9⟩, order 18 | 3, 3; 19, 19 | [C370.json](C370.json) |
| 406 | 203 | Paley's matrix of order 812, q = 811 ([PALEY.md](PALEY.md)) | – | [C406.json](C406.json) |
| 414 | 207 | Paley's matrix of order 828, q = 827 ([PALEY.md](PALEY.md)) | – | [C414.json](C414.json) |
| 430 | 215 = 5 · 43 | ⟨41⟩, order 7 (exact k-tree on one compression class, `search/lowe_comp.py`) | 15, 15; 11, 17 | [C430.json](C430.json) |
| 454 | 227 | bordered array ([BORDERED.md](BORDERED.md)) on ℤ₂₂₆, ⟨49⟩ of order 7 (low-energy k-tree on one compression class, `search/lowe_comp.py`) | 0, 0; 2, 0 | [C454.json](C454.json) |
| 474 | 237 | Paley's matrix of order 948, q = 947 ([PALEY.md](PALEY.md)) | – | [C474.json](C474.json) |
| 526 | 263 | Paley's matrix of order 1052, q = 1051 ([PALEY.md](PALEY.md)) | – | [C526.json](C526.json) |
| 546 | 273 | Paley's matrix of order 1092, q = 1091 ([PALEY.md](PALEY.md)) | – | [C546.json](C546.json) |
| 582 | 291 | Paley's matrix of order 1164, q = 1163 ([PALEY.md](PALEY.md)) | – | [C582.json](C582.json) |
| 646 | 323 | Paley's matrix of order 1292, q = 1291 ([PALEY.md](PALEY.md)) | – | [C646.json](C646.json) |
| 666 | 333 | Paley's matrix of order 1332, q = 1331 = 11³ ([PALEY.md](PALEY.md)) | – | [C666.json](C666.json) |
| 714 | 357 | Paley's matrix of order 1428, q = 1427 ([PALEY.md](PALEY.md)) | – | [C714.json](C714.json) |
| 742 | 371 = 7 · 53 | ⟨13⟩, order 26 | 3, 3; 25, 29 | [C742.json](C742.json) |
| 766 | 383 | Paley's matrix of order 1532, q = 1531 ([PALEY.md](PALEY.md)) | – | [C766.json](C766.json) |
| 786 | 393 | Paley's matrix of order 1572, q = 1571 ([PALEY.md](PALEY.md)) | – | [C786.json](C786.json) |
| 790 | 395 | Paley's matrix of order 1580, q = 1579 ([PALEY.md](PALEY.md)) | – | [C790.json](C790.json) |
| 814 | 407 | Paley's matrix of order 1628, q = 1627 ([PALEY.md](PALEY.md)) | – | [C814.json](C814.json) |
| 834 | 417 | Paley's matrix of order 1668, q = 1667 ([PALEY.md](PALEY.md)) | – | [C834.json](C834.json) |
| 850 | 425 | Paley's matrix of order 1700, q = 1699 ([PALEY.md](PALEY.md)) | – | [C850.json](C850.json) |
| 862 | 431 | Paley's matrix of order 1724, q = 1723 ([PALEY.md](PALEY.md)) | – | [C862.json](C862.json) |
| 874 | 437 | Paley's matrix of order 1748, q = 1747 ([PALEY.md](PALEY.md)) | – | [C874.json](C874.json) |
| 894 | 447 | Paley's matrix of order 1788, q = 1787 ([PALEY.md](PALEY.md)) | – | [C894.json](C894.json) |
| 906 | 453 | Paley's matrix of order 1812, q = 1811 ([PALEY.md](PALEY.md)) | – | [C906.json](C906.json) |
| 934 | 467 | Paley's matrix of order 1868, q = 1867 ([PALEY.md](PALEY.md)) | – | [C934.json](C934.json) |
| 966 | 483 | Paley's matrix of order 1932, q = 1931 ([PALEY.md](PALEY.md)) | – | [C966.json](C966.json) |
| 990 | 495 | Paley's matrix of order 1980, q = 1979 ([PALEY.md](PALEY.md)) | – | [C990.json](C990.json) |
| 994 | 497 | Paley's matrix of order 1988, q = 1987 ([PALEY.md](PALEY.md)) | – | [C994.json](C994.json) |
| 1002 | 501 | Paley's matrix of order 2004, q = 2003 ([PALEY.md](PALEY.md)) | – | [C1002.json](C1002.json) |
| 1006 | 503 | Paley's matrix of order 2012, q = 2011 ([PALEY.md](PALEY.md)) | – | [C1006.json](C1006.json) |
| 1066 | 533 | Paley's matrix of order 2132, q = 2131 ([PALEY.md](PALEY.md)) | – | [C1066.json](C1066.json) |
| 1102 | 551 | Paley's matrix of order 2204, q = 2203 ([PALEY.md](PALEY.md)) | – | [C1102.json](C1102.json) |
| 1126 | 563 | Paley's matrix of order 2252, q = 2251 ([PALEY.md](PALEY.md)) | – | [C1126.json](C1126.json) |
| 1134 | 567 | Paley's matrix of order 2268, q = 2267 ([PALEY.md](PALEY.md)) | – | [C1134.json](C1134.json) |
| 1170 | 585 | Paley's matrix of order 2340, q = 2339 ([PALEY.md](PALEY.md)) | – | [C1170.json](C1170.json) |
| 1174 | 587 | Paley's matrix of order 2348, q = 2347 ([PALEY.md](PALEY.md)) | – | [C1174.json](C1174.json) |
| 1206 | 603 | Paley's matrix of order 2412, q = 2411 ([PALEY.md](PALEY.md)) | – | [C1206.json](C1206.json) |
| 1266 | 633 | Paley's matrix of order 2532, q = 2531 ([PALEY.md](PALEY.md)) | – | [C1266.json](C1266.json) |
| 1270 | 635 | Paley's matrix of order 2540, q = 2539 ([PALEY.md](PALEY.md)) | – | [C1270.json](C1270.json) |
| 1330 | 665 | Paley's matrix of order 2660, q = 2659 ([PALEY.md](PALEY.md)) | – | [C1330.json](C1330.json) |
| 1342 | 671 | Paley's matrix of order 2684, q = 2683 ([PALEY.md](PALEY.md)) | – | [C1342.json](C1342.json) |
| 1366 | 683 | Paley's matrix of order 2732, q = 2731 ([PALEY.md](PALEY.md)) | – | [C1366.json](C1366.json) |
| 1422 | 711 | Paley's matrix of order 2844, q = 2843 ([PALEY.md](PALEY.md)) | – | [C1422.json](C1422.json) |
| 1426 | 713 | Paley's matrix of order 2852, q = 2851 ([PALEY.md](PALEY.md)) | – | [C1426.json](C1426.json) |
| 1486 | 743 | Paley's matrix of order 2972, q = 2971 ([PALEY.md](PALEY.md)) | – | [C1486.json](C1486.json) |
| 1506 | 753 | Paley's matrix of order 3012, q = 3011 ([PALEY.md](PALEY.md)) | – | [C1506.json](C1506.json) |
| 1510 | 755 | Paley's matrix of order 3020, q = 3019 ([PALEY.md](PALEY.md)) | – | [C1510.json](C1510.json) |
| 1534 | 767 | Paley's matrix of order 3068, q = 3067 ([PALEY.md](PALEY.md)) | – | [C1534.json](C1534.json) |
| 1542 | 771 | Paley's matrix of order 3084, q = 3083 ([PALEY.md](PALEY.md)) | – | [C1542.json](C1542.json) |
| 1582 | 791 | Paley's matrix of order 3164, q = 3163 ([PALEY.md](PALEY.md)) | – | [C1582.json](C1582.json) |
| 1626 | 813 | Paley's matrix of order 3252, q = 3251 ([PALEY.md](PALEY.md)) | – | [C1626.json](C1626.json) |
| 1630 | 815 | Paley's matrix of order 3260, q = 3259 ([PALEY.md](PALEY.md)) | – | [C1630.json](C1630.json) |
| 1650 | 825 | Paley's matrix of order 3300, q = 3299 ([PALEY.md](PALEY.md)) | – | [C1650.json](C1650.json) |
| 1654 | 827 | Paley's matrix of order 3308, q = 3307 ([PALEY.md](PALEY.md)) | – | [C1654.json](C1654.json) |
| 1662 | 831 | Paley's matrix of order 3324, q = 3323 ([PALEY.md](PALEY.md)) | – | [C1662.json](C1662.json) |
| 1666 | 833 | Paley's matrix of order 3332, q = 3331 ([PALEY.md](PALEY.md)) | – | [C1666.json](C1666.json) |
| 1674 | 837 | Paley's matrix of order 3348, q = 3347 ([PALEY.md](PALEY.md)) | – | [C1674.json](C1674.json) |
| 1746 | 873 | Paley's matrix of order 3492, q = 3491 ([PALEY.md](PALEY.md)) | – | [C1746.json](C1746.json) |
| 1750 | 875 | Paley's matrix of order 3500, q = 3499 ([PALEY.md](PALEY.md)) | – | [C1750.json](C1750.json) |
| 1774 | 887 | Paley's matrix of order 3548, q = 3547 ([PALEY.md](PALEY.md)) | – | [C1774.json](C1774.json) |
| 1786 | 893 | Paley's matrix of order 3572, q = 3571 ([PALEY.md](PALEY.md)) | – | [C1786.json](C1786.json) |
| 1822 | 911 | Paley's matrix of order 3644, q = 3643 ([PALEY.md](PALEY.md)) | – | [C1822.json](C1822.json) |
| 1830 | 915 | Paley's matrix of order 3660, q = 3659 ([PALEY.md](PALEY.md)) | – | [C1830.json](C1830.json) |
| 1846 | 923 | Paley's matrix of order 3692, q = 3691 ([PALEY.md](PALEY.md)) | – | [C1846.json](C1846.json) |
| 1870 | 935 | Paley's matrix of order 3740, q = 3739 ([PALEY.md](PALEY.md)) | – | [C1870.json](C1870.json) |
| 1926 | 963 | Paley's matrix of order 3852, q = 3851 ([PALEY.md](PALEY.md)) | – | [C1926.json](C1926.json) |
| 1962 | 981 | Paley's matrix of order 3924, q = 3923 ([PALEY.md](PALEY.md)) | – | [C1962.json](C1962.json) |
| 1966 | 983 | Paley's matrix of order 3932, q = 3931 ([PALEY.md](PALEY.md)) | – | [C1966.json](C1966.json) |

Each certificate gives N, the multiplier group, the row sums, the first rows of A, B, C, D (`+`, `-`) (for the
bordered array also the group of order N − 1 and the border; for the character-bordered array the group, the fibre
labels of φ and the three border matrices; for Paley's matrix only q, and the modulus of GF(q) when q is not prime)
and the SHA-256 of the canonical text of K: 2N lines of 2N digits, digit e standing for the entry iᵉ.

```
python3 verify/verify_complex.py            # all certificates
python3 verify/verify_complex.py --write 142 # also writes complex-hadamard/C142.txt
```

This rebuilds K from the four sequences alone (or from q, with its own arithmetic in GF(q)) and checks every entry
and K K* = 2N I exactly, in integer arithmetic on the real and imaginary parts. All 87 certificates take about 13
minutes, most of it the orders above 1000.

## Status

**Primary source.** Table 11.2 of J. Seberry and M. Yamada, "Hadamard matrices, sequences, and block designs"
(Contemporary Design Theory, Wiley, 1992) records, for every odd q < 1000, whether a complex Hadamard matrix of
order 2q was known. Its constructions are conference matrices (codes c1–c6) and Williamson-type matrices (Turyn;
Kharaghani–Seberry, codes x2 and x6). Every order in this directory is **blank** in that table, although the half-order
lemma behind PALEY.md is in the same survey (Lemma 11.16). The parsed table is
[`search/data/complex_known_1992.txt`](../search/data/complex_known_1992.txt).

**Results since 1992** that we checked. The first gives the 76 known orders of this directory, the last the 11
concurrent ones:
- Ito's dicyclic Hadamard groups (Ito, On Hadamard groups, 1994–2000; Schmidt, Des. Codes Cryptogr. 17 (1999) 61–68,
  Corollary 3.6) and R. Egan, Phased unitary Golay pairs, Butson Hadamard matrices and a conjecture of Ito's,
  Des. Codes Cryptogr. 87 (2019) 67–74, Theorem 4.5: a complex Hadamard matrix of order 2v whenever v = 2^a 10^b 26^c m
  with 4m − 1 or 2m − 1 a prime power, or m odd with a Williamson matrix of order m, and for every v ≤ 46 (Ito's
  conjecture is verified there). For N odd this covers 4N − 1 a prime power; 2N − 1 and Williamson matrices were in
  the 1992 table already.
- Đoković's good matrices of orders 35 and 127 (1993): 70 and 254.
- Szöllősi (2026): 94.
- Kharaghani, Mohammadian, Tayfeh-Rezaie and Zaitsev (2026): 118, 130, 134, 206 and 266.
- Quaternary Legendre pairs, giving order 2ℓ + 2 (Kotsireas–Winterhof 2022; Kotsireas–Koutschan–Winterhof 2025;
  Jedwab–Pender 2025; Lebedev and Verma 2026): every even length ℓ up to 100 except 70, 76, 88, 92, 94 and 100,
  plus two infinite families. Order 142 needs ℓ = 70 and 154 needs ℓ = 76, both open. Our larger orders fall
  outside both families.
- Skew-regular quaternary Hadamard matrices of order 1 + p², p an odd prime power (arXiv:2510.06489, 2025): these
  orders are known already from Turyn's theorem, since 2N − 1 = p² is a prime power.
- The corollaries of Miyamoto's lemma (Seberry–Yamada 1992, Corollary 8.8; the M-structure paper, Corollaries 16 and
  18) with inputs found after 1992, which give Williamson-type matrices and so complex matrices of order 2N:
  N = 509 (good matrices of order 127), 857 (H(428)) and 883 (q = 441, Mathon's conference matrix of order 442),
  stated by Đoković (arXiv:1008.2043, Propositions 3.2–3.4), and N = 269 (near-Williamson matrices of order 67, 2026;
  this combination is ours). So 538, 1018, 1714 and 1766 are known. We could not rebuild the matrices behind
  Corollary 8.8(i) ([NOTES.md](../NOTES.md) §19) and follow the literature in counting them.
- Mathon's conference matrices (q²(q ± 2) + 1). We also count them as known where the 1992 table is blank
  (N = 221, 787).
- Zaitsev's page (above, read on 16 October 2026). Its 122 matrices were downloaded and checked exactly
  ([`search/data/complex_zaitsev_2026.txt`](../search/data/complex_zaitsev_2026.txt)). Besides the 11 orders of the
  first table it covers 24 orders that we had listed as open, among them 470, 490, 686, 1306, 1466, 1706 and 1882.

The sweep's candidate list (`python3 search/ssxx_sweep.py list`) is exactly the blank entries of the 1992 table
minus these. The 2026 papers cite the updated Table 11.2 of Seberry and Yamada's 2020 book (Wiley), which we could
not consult. A construction published between 1992 and 2020 and recorded only there is the remaining caveat.

**What remains open** (`python3 search/complex_status.py`): after the orders in this directory and on Zaitsev's page,
115 values of N < 1000 have no complex Hadamard matrix of order 2N known to us. The smallest open orders are 334, 358,
418, 438, 446, 478, 498, 502, 506, 550, 574, 598, …; for none of them is 4N − 1 a prime power. (The earlier count of 143
also listed 538, 1018, 1714 and 1766, which follow from Miyamoto's corollaries, and that of 139 the 24 orders of
Zaitsev's page; see above.) 214, 238, 302 and 326 were successively the
smallest open order; they are in the table above. No quaternary Legendre pair family covers the length 106 of order
214 (Jedwab–Pender need 213 or 105 to be a prime power), and the 2026 quaternary Legendre pairs (arXiv:2609.04589,
arXiv:2609.25573) have lengths below 100.

**Correction.** An earlier version of this directory called 76 orders New that follow from Ito and Egan: the 63 halves
of Paley's matrix and 13 orders found by search, among them 142 and 154 (which we had called the two smallest open
orders after the 2026 papers) and 322. Neither 2026 paper cites Egan, and both take their smallest open order from the
1992 table; we did the same. The record of how the error arose is in [PALEY.md](PALEY.md#novelty).

A second correction: the 11 orders of the first table were labelled New until we found Zaitsev's page on 16 October
2026. They are Concurrent: the page lists all 11, and its last modification (15 October 2026, 02:50 UTC) precedes our
commits of 214 and 366.
