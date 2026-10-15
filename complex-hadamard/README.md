# Complex Hadamard matrices of orders 2N, N odd

A quaternary (complex) Hadamard matrix of order n has entries 1, −1, i, −i and satisfies K K* = n I. Its order is 1
or even, and it is conjectured that one exists for every even order (Turyn). Orders divisible by 4 are covered by real
Hadamard matrices. The open cases are therefore the orders 2N with N odd.

Recent progress at the bottom of that range:
- Szöllősi constructed order 94 (arXiv:2603.09572, March 2026), noting that "the next outstanding order seems to be"
  118.
- Kharaghani, Mohammadian, Tayfeh-Rezaie and Zaitsev (arXiv:2605.08661v2, September 2026) gave the first examples of
  orders 118, 130, 134, 206 and 266, from near-Williamson matrices.

This directory adds orders that, to our knowledge, had no construction.

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
the live patterns have far more than the plain Gaussian says ([search/NOTES.md](../search/NOTES.md) §42). Order 366 came
from ⟨34⟩ on ℤ₁₈₃: its pattern 1, 23; 9, 11 lives on a single class of ℤ₃ compressions (A ≡ (−1, 1, 1), B ≡ (1, 11,
11), C ≡ (−1, −1, 11), D ≡ (−11, 11, 11)), which the unconditioned low-energy lists had missed.

A second construction, the bordered Szöllősi array on a group of order N − 1, is stated and proved in
[BORDERED.md](BORDERED.md). Order 454 comes from it: 227 is a safe prime, so ℤ₂₂₇ has no usable multiplier group,
while ℤ₂₂₆ = ℤ₂ × ℤ₁₁₃ has ⟨49⟩ of order 7. Its fibre sums over ℤ₂ are ≡ ±1 (mod 14), which the Gaussian estimate
ignores (10^0.6); conditioned on the compressions to ℤ₂ (`search/comp_est.py`) the estimate is 10^6.0, and the
low-energy k-tree on one compression class (pools 4900, 4900, 1.7 · 10⁸, 1.3 · 10⁸, lists of 2 · 10⁷) found 64
solutions in six minutes ([search/NOTES.md](../search/NOTES.md) §42). Two twisted specialisations of both arrays (B = A∘μ, D = C∘μ, and a half twist), which
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
| 330 | 165 = 3 · 5 · 11 | ⟨136⟩, order 5 (Paley product ℤ₁₅ × GF(11), `search/ssxx_paley.py`) | 7, 9; 1, 23 | [C330.json](C330.json) |
| 366 | 183 = 3 · 61 | ⟨34⟩, order 5 (low-energy k-tree on one compression class, `search/lowe_comp.py`) | 1, 23; 9, 11 | [C366.json](C366.json) |
| 370 | 185 = 5 · 37 | ⟨9⟩, order 18 | 3, 3; 19, 19 | [C370.json](C370.json) |
| 382 | 191 | ⟨39⟩, order 5 (streaming k-tree, `search/kstream.py`) | 11, 11; 9, 21 | [C382.json](C382.json) |
| 426 | 213 = 3 · 71 | ⟨37⟩, order 7 (low-energy k-tree, `search/lowe_sz.py`) | 13, 17; 13, 15 | [C426.json](C426.json) |
| 454 | 227 | bordered array ([BORDERED.md](BORDERED.md)) on ℤ₂₂₆, ⟨49⟩ of order 7 (low-energy k-tree on one compression class, `search/lowe_comp.py`) | 0, 0; 2, 0 | [C454.json](C454.json) |
| 494 | 247 = 13 · 19 | ⟨9, 35⟩, order 27 | 17, 19; 7, 17 | [C494.json](C494.json) |
| 518 | 259 = 7 · 37 | ⟨44⟩, order 9 (compression chain) | 5, 29; 1, 13 | [C518.json](C518.json) |
| 534 | 267 = 3 · 89 | ⟨4⟩, order 11 | 1, 1; 21, 25 | [C534.json](C534.json) |
| 742 | 371 = 7 · 53 | ⟨13⟩, order 26 | 3, 3; 25, 29 | [C742.json](C742.json) |
| 754 | 377 = 13 · 29 | ⟨9⟩, order 42 | 27, 27; 5, 5 | [C754.json](C754.json) |

Each certificate gives N, the multiplier group, the row sums, the first rows of A, B, C, D (`+`, `-`) (for the
bordered array also the group of order N − 1 and the border) and the SHA-256 of the canonical text of K: 2N lines of 2N digits, digit e standing for the entry iᵉ.

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
- Szöllősi (2026): 94.
- Kharaghani, Mohammadian, Tayfeh-Rezaie and Zaitsev (2026): 118, 130, 134, 206 and 266.
- Quaternary Legendre pairs, giving order 2ℓ + 2 (Kotsireas–Winterhof 2022; Kotsireas–Koutschan–Winterhof 2025;
  Jedwab–Pender 2025; Lebedev and Verma 2026): every even length ℓ up to 100 except 70, 76, 88, 92, 94 and 100,
  plus two infinite families. Order 142 needs ℓ = 70 and 154 needs ℓ = 76, both open. Our larger orders fall
  outside both families.
- Mathon's conference matrices (q²(q ± 2) + 1). We also count them as known where the 1992 table is blank
  (N = 221, 787).

The sweep's candidate list (`python3 search/ssxx_sweep.py list`) is exactly the blank entries of the 1992 table
minus these. The 2026 papers cite the updated Table 11.2 of Seberry and Yamada's 2020 book (Wiley), which we could
not consult. A construction published between 1992 and 2020 and recorded only there is the remaining caveat.

**What remains open** (`python3 search/complex_status.py`): after the orders in this directory, 210 values of N < 1000
have no complex Hadamard matrix of order 2N known to us. The smallest open orders are 214, 322, 334, 342, 358, 406,
414, 418, 430, 438, 446, 470, …. Order 214 (N = 107) is out of reach of every structure we have: 107 and 53 = (107 −
1)/2 are both prime, so neither ℤ₁₀₇ nor ℤ₁₀₆ has a usable multiplier group, and no quaternary Legendre pair family
covers length 106.
