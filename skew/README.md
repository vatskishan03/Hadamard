# Skew-Hadamard matrices

A Hadamard matrix H is skew if H + Hᵀ = 2I. Skew-Hadamard matrices are conjectured to exist for every order divisible
by 4. Each directory here holds one order: a certificate, the matrix (where it is small enough to store) and a README
with the construction, the search and the evidence that the order is new. The order 3732 is in
[hadamard/3732](../hadamard/3732), since it was also an unresolved Hadamard order.

**Before** gives the entry of Cati–Pasechnik Table 5 (v ≤ 999): v(t) means that the smallest skew-Hadamard order known
with odd part v was 2ᵗv. **Census** is the status of v in the Đoković–Seberry census of skew orders 4v, odd v < 2500.
Doubling (Seberry) turns a skew-Hadamard matrix of order 4v into ones of every order 2ᵏ · 4v. Status labels are
defined in [NOVELTY.md](../NOVELTY.md).

## Goethals–Seidel difference families

Four ±1 sequences on a group of order v, invariant under a group H of automorphisms with −1 ∉ H, the first of skew
type, with periodic autocorrelations summing to zero off the identity. The Goethals–Seidel array then gives a
skew-Hadamard matrix of order 4v. The group is ℤᵥ unless stated.

| Order | v | Symmetry H | Row sums | Before | Status |
|------:|---|---|---|---|---|
| [820](820) | 205 = 5 · 41 | ⟨51⟩, order 5 (cyclotomic product of index 8) | 1, 1, −17, −23 | 205(3) | New |
| [836](836) | 209 = 11 · 19 | ⟨26⟩, order 15 | 1, 9, 15, 23 | 209(4) | New |
| [1012](1012) | 253 = 11 · 23 | ⟨12⟩, order 11 (Paley product) | 1, 5, 19, 25 | 253(4) | New |
| [1036](1036) | 259 = 7 · 37 | ⟨9⟩, order 9 | 1, 7, 31, 5 | 259(4) | New |
| [1044](1044) | 261 = 3² · 29 | ⟨190⟩, order 7 (cyclotomic product of index 4) | 1, 9, −1, 31 | 261(3) | New |
| [1100](1100) | 275 = 5² · 11, group ℤ₅ × ℤ₅ × ℤ₁₁ | order 15: an element of order 3 of the Singer cycle of ℤ₅², and ⟨3⟩ on ℤ₁₁ | 1, 3, 33, 1 | 275(4) | New |
| [1140](1140) | 285 = 3 · 5 · 19 | ⟨61⟩, order 9 (Paley product) | 1, 17, 15, 25 | 285(3) | New; see its README |
| [1148](1148) | 287 = 7 · 41 | ⟨16⟩, order 15 | 1, 3, 33, 7 | 287(4) | New |
| [1204](1204) | 301 = 7 · 43 | ⟨127⟩, order 7 (cyclotomic product of index 6) | 1, −11, −11, −31 | 301(3) | New |
| [1236](1236) | 309 = 3 · 103 | ⟨13⟩, order 17 | 1, 3, 35, 1 | 309(3) | New |
| [1276](1276) | 319 = 11 · 29 | ⟨45⟩, order 7 (cyclotomic product of index 4) | 1, 25, 19, 17 | 319(3) | New |
| [1364](1364) | 341 = 11 · 31 | ⟨14⟩, order 15 | 1, 9, 21, 29 | 341(4) | New |
| [1372](1372) | 343 = 7³, group GF(343) | the subgroup of order 19 of GF(343)* | 1, 37, 1, 1 | 343(4) | New |
| [1380](1380) | 345 = 3 · 5 · 23 | ⟨301⟩, order 11 (Paley product) | 1, 19, 17, 27 | 345(4) | New |
| [1412](1412) | 353, prime | ⟨22⟩, order 11 | 1, 21, 23, 21 | 353(4) | New |
| [1508](1508) | 377 = 13 · 29 | ⟨16⟩, order 21 | 1, 27, 27, 7 | 377(5) | New |
| [1564](1564) | 391 = 17 · 23 | ⟨324⟩, order 11 (Paley product) | 1, 19, 19, 29 | 391(4) | New |
| [1612](1612) | 403 = 13 · 31 | ⟨100⟩, order 15 | 1, 5, 25, 31 | 403(5) | New |
| [1804](1804) | 451 = 11 · 41 | ⟨16, 42⟩, order 25 | 1, 1, 11, 41 | none; see its README | New |
| [1852](1852) | 463, prime | ⟨33⟩, order 21 | 1, 1, 43, 1 | 463(7), out of date | Known (Đoković 1994) |
| [1860](1860) | 465 = 3 · 5 · 31 | ⟨226⟩, order 15 (Paley product) | 1, 25, 35, −3 | 465(3) | New |
| [1876](1876) | 469 = 7 · 67 | ⟨9⟩, order 33 | 1, 1, 5, 43 | 469(3) | New |
| [2052](2052) | 513 = 3³ · 19, group ℤ₁₉ × GF(27) | the squares of GF(27)*, order 13 (Paley product) | 1, 27, 19, 31 | 513(4) | New |
| [2404](2404) | 601, prime | ⟨2⟩, order 25 | 1, 1, 49, 1 | 601(5) | New |
| [2452](2452) | 613, prime | ⟨37⟩, order 17 | 1, 1, 35, 35 | 613(3) | New |
| [2812](2812) | 703 = 19 · 37 | ⟨16, 26⟩, order 27 | 1, 1, 1, 53 | 703(3) | New |
| [3052](3052) | 763 = 7 · 109 | ⟨9⟩, order 27 | 1, 1, 5, 55 | 763(4) | New |
| [3380](3380) | 845 = 5 · 13², group ℤ₅ × GF(169) | ⟨w⁸⟩ on GF(169), order 21 (cyclotomic product of index 8) | 1, 43, −3, −39 | 845(4) | New |
| [3556](3556) | 889 = 7 · 127 | ⟨2, 25⟩, order 63 | 1, 5, 41, 43 | 889(4) | New |
| [3612](3612) | 903 = 3 · 7 · 43 | ⟨568⟩, order 21 (Paley product) | 1, 43, 41, −9 | 903(4) | New |
| [3676](3676) | 919, prime | ⟨58⟩, order 17, blocks related by a multiplier of order 3 | 1, 35, 35, 35 | — | Known (Schumer 2026) |
| [3732](../hadamard/3732) | 933 = 3 · 311 | ⟨7⟩, order 31 (cyclotomic product of index 10) | 1, 61, 3, 1 | 933(10) | New |
| [5668](5668) | 1417 = 13 · 109 | ⟨1405, 328⟩, order 81 (cyclotomic product of index 4 with ⟨3⟩ on ℤ₁₃) | 1, 7, −53, −53 | beyond Table 5; 45344 = 32 · 1417 | New |
| [6004](6004) | 1501 = 19 · 79 | ⟨799⟩, order 39 (Paley product) | 1, 77, 7, −5 | beyond Table 5; 12008 = 8 · 1501 | New |

## Seberry's theorem with new ingredients

Seberry (1978): amicable Hadamard matrices of order n and a skew-type orthogonal design OD(mn; 1, m, mn − m − 1) give
a skew-Hadamard matrix of order mn(n − 1). The orders below use amicable pairs or orthogonal designs that were not
known. The amicable pairs of orders 52, 76 and 100 come from skew-propus quadruples, that is from symmetric orthogonal
designs of type (1, n − 1) ([search/NOTES.md](../search/NOTES.md) §31). Seberry's tables of amicable orders 2ᵗq (1978,
1991, 1992) give t = 3 for q = 13, 19 and 25, so these pairs are new as far as we know. For m = 1 the skew-type
designs come from G-matrices, which are known for these orders. The designs for m ≥ 3 were found here.

| Order | v | Amicable pair | Skew-type orthogonal design | Census | Status |
|------:|---|---|---|---|---|
| [1716](1716) | 429 = 13 · 3 · 11 | order 12, Paley (q = 11) | OD(156; 1, 13, 142) over ℤ₃₉, new | undecided; Table 5: 429(3) | New |
| [1900](1900) | 475 = 5 · 5 · 19 | order 20, Paley (q = 19) | OD(100; 1, 5, 94) over ℤ₂₅, new | undecided; Table 5: 475(4) | New |
| [2652](2652) | 663 = 13 · 51 | order 52, new (skew-propus quadruple of length 13) | G-matrices of order 13 | undecided; Table 5 implies known | First explicit; see its README |
| [2772](2772) | 693 = 21 · 3 · 11 | order 12, Paley (q = 11) | OD(252; 1, 21, 230) over ℤ₆₃, new | undecided; Table 5: 693(4) | New |
| [5700](5700) | 1425 = 19 · 75 | order 76, new (the design of [symmetric/5700](../symmetric/5700)) | G-matrices of order 19 | undecided | New |
| [7956](7956) | 1989 = 3 · 13 · 51 | order 52, as for 2652 | OD(156; 1, 3, 152) over ℤ₃₉, new | undecided | New |
| [9828](9828) | 2457 = 13 · 7 · 27 | order 28, Paley (q = 27) | OD(364; 1, 13, 350) over ℤ₉₁, new | undecided | New |
| [9900](9900) | 2475 = 25 · 99 | order 100, new (skew-propus quadruple of length 25) | G-matrices of order 25 | undecided | New |

The matrix of order 9900 is not stored; `verify/verify.py skew/9900` rebuilds it.

## How they were found

The Goethals–Seidel families come from searches over multiplier structures: exact joins of complete block lists for
small structures, exact compression chains (`search/chain.py`, run over the census by `search/census_screen.py` and
`search/screen_chain.py`), Paley-product and cyclotomic-product configurations (`search/paley_config.py`,
`search/cyclo_product.py`), low-energy k-trees on lists too large to store (`search/cp_lowe.py`,
`search/skew_lowE.py`) and, for 3676, a spin structure (`search/spin_skew.py`). The orthogonal designs come from
`search/odm_engine.py` and `search/odm_join.py`. [search/NOTES.md](../search/NOTES.md) §10–§16, §20, §23 and §31
describe the methods, and each README records the command, code commit, seed and caps of its search.

`python3 verify/verify.py skew/<order>` rebuilds a matrix from its certificate and checks H Hᵀ = nI and H + Hᵀ = 2I.
