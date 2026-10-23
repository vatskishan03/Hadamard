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
| [476](476) | 119 = 7 · 17 | ⟨86⟩, order 3 (Paley product with nine w-positions) | 1, −5, 15, −15 | 119(4) | New |
| [820](820) | 205 = 5 · 41 | ⟨51⟩, order 5 (cyclotomic product of index 8) | 1, 1, −17, −23 | 205(3) | New |
| [836](836) | 209 = 11 · 19 | ⟨26⟩, order 15 | 1, 9, 15, 23 | 209(4) | New |
| [1012](1012) | 253 = 11 · 23 | ⟨12⟩, order 11 (Paley product) | 1, 5, 19, 25 | 253(4) | New |
| [1036](1036) | 259 = 7 · 37 | ⟨9⟩, order 9 | 1, 7, 31, 5 | 259(4) | New |
| [1044](1044) | 261 = 3² · 29 | ⟨190⟩, order 7 (cyclotomic product of index 4) | 1, 9, −1, 31 | 261(3) | New |
| [1100](1100) | 275 = 5² · 11, group ℤ₅ × ℤ₅ × ℤ₁₁ | order 15: an element of order 3 of the Singer cycle of ℤ₅², and ⟨3⟩ on ℤ₁₁ | 1, 3, 33, 1 | 275(4) | New |
| [1140](1140) | 285 = 3 · 5 · 19 | ⟨61⟩, order 9 (Paley product) | 1, 17, 15, 25 | 285(3) | New |
| [1148](1148) | 287 = 7 · 41 | ⟨16⟩, order 15 | 1, 3, 33, 7 | 287(4) | New |
| [1204](1204) | 301 = 7 · 43 | ⟨127⟩, order 7 (cyclotomic product of index 6) | 1, −11, −11, −31 | 301(3) | New |
| [1236](1236) | 309 = 3 · 103 | ⟨13⟩, order 17 | 1, 3, 35, 1 | 309(3) | New |
| [1276](1276) | 319 = 11 · 29 | ⟨45⟩, order 7 (cyclotomic product of index 4) | 1, 25, 19, 17 | 319(3) | New |
| [1364](1364) | 341 = 11 · 31 | ⟨14⟩, order 15 | 1, 9, 21, 29 | 341(4) | New |
| [1372](1372) | 343 = 7³, group GF(343) | the subgroup of order 19 of GF(343)* | 1, 37, 1, 1 | 343(4) | New |
| [1380](1380) | 345 = 3 · 5 · 23 | ⟨301⟩, order 11 (Paley product) | 1, 19, 17, 27 | 345(4) | New |
| [1412](1412) | 353, prime | ⟨22⟩, order 11 | 1, 21, 23, 21 | 353(4) | New |
| [1508](1508) | 377 = 13 · 29 | ⟨16⟩, order 21 | 1, 27, 27, 7 | 377(5) | New |
| [1540](1540) | 385 = 5 · 7 · 11 | ⟨246, 221⟩, order 15 (Paley product with a multiplier on ℤ₃₅) | 1, −7, 31, 23 | 385(3) | New |
| [1564](1564) | 391 = 17 · 23 | ⟨324⟩, order 11 (Paley product) | 1, 19, 19, 29 | 391(4) | New |
| [1612](1612) | 403 = 13 · 31 | ⟨100⟩, order 15 | 1, 5, 25, 31 | 403(5) | New |
| [1804](1804) | 451 = 11 · 41 | ⟨16, 42⟩, order 25 | 1, 1, 11, 41 | none; see its README | New |
| [1852](1852) | 463, prime | ⟨33⟩, order 21 | 1, 1, 43, 1 | 463(7), out of date | Known (Đoković 1994) |
| [1860](1860) | 465 = 3 · 5 · 31 | ⟨226⟩, order 15 (Paley product) | 1, 25, 35, −3 | 465(3) | New |
| [1876](1876) | 469 = 7 · 67 | ⟨9⟩, order 33 | 1, 1, 5, 43 | 469(3) | New |
| [1892](1892) | 473 = 11 · 43 | ⟨345, 221⟩, order 15 (Paley product with a multiplier of order 3 on ℤ₄₃) | 1, −9, 17, −39 | 473(4) | New |
| [2052](2052) | 513 = 3³ · 19, group ℤ₁₉ × GF(27) | the squares of GF(27)*, order 13 (Paley product) | 1, 27, 19, 31 | 513(4) | New |
| [2156](2156) | 539 = 7² · 11, group ℤ₇ × ℤ₇ × ℤ₁₁ | order 15: the scalar 2 (order 3) on ℤ₇², and the squares ⟨4⟩ (order 5) on ℤ₁₁ (Paley product with a multiplier on ℤ₇²) | 1, −29, 15, 33 | 539(4) | New |
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

## Negacyclic Goethals–Seidel arrays

Four ±1 sequences of length 2v (v odd) whose negaperiodic autocorrelations sum to 8v δ, the first of skew type, in
the Goethals–Seidel array with negacyclic blocks; searched as quaternary sequences on ℤᵥ constant on the orbits of a
multiplier group H of odd order (for 2392 and 4712 as Paley plug-ins over ℤ_m × F_q). The array gives a
skew-Hadamard matrix of order 8v; see [NEGACYCLIC.md](NEGACYCLIC.md).

| Order | v | Symmetry H | Sums of the sequences | Before | Status |
|------:|---|---|---|---|---|
| [1528](1528) | 191, prime | ⟨180⟩, order 19 | 6, −6, −22, 18 | 191(226) | New |
| [2216](2216) | 277, prime | ⟨273⟩, order 23 | −28, −34, −14, 32 | 277(5) | New |
| [2392](2392) | 299 = 13 · 23 | ⟨185, 209⟩, order 33 (Paley plug-in over ℤ₁₃ × F₂₃) | −60, −8, 12, 68 | 299(4) | New |
| [4712](4712) | 589 = 19 · 31 | ⟨311, 381⟩, order 45 (Paley plug-in over ℤ₁₉ × F₃₁) | −10, 34, −46, 12 | 589(5) | New |

## Seberry's theorem with new ingredients

Seberry (1978): amicable Hadamard matrices of order n and a skew-type orthogonal design OD(mn; 1, m, mn − m − 1) give a
skew-Hadamard matrix of order mn(n − 1). The orders below use amicable pairs or orthogonal designs that no source
states. The amicable pairs of orders 52, 76 and 100 come from skew-propus quadruples, that is from symmetric orthogonal
designs of type (1, n − 1) ([search/NOTES.md](../search/NOTES.md) §31). Seberry's tables of amicable orders 2ᵗq (1978,
1991, 1992) give t = 3 for q = 13, 19 and 25, and no source we know states amicable pairs of these orders, but the
skew-propus difference families behind them were published by Abuzin, Balonin, Đoković and Kotsireas (Inform. Control
Syst. 2019, no. 5, 2–9, Table 1; doi:10.31799/1684-8853-2019-5-2-9), so the pairs are first explicit, not new objects.
For m = 1 the skew-type designs come from G-matrices, which are known for these orders. The designs for m ≥ 3 were found
here.

| Order | v | Amicable pair | Skew-type orthogonal design | Census | Status |
|------:|---|---|---|---|---|
| [1716](1716) | 429 = 13 · 3 · 11 | order 12, Paley (q = 11) | OD(156; 1, 13, 142) over ℤ₃₉, new | undecided; Table 5: 429(3) | New |
| [1900](1900) | 475 = 5 · 5 · 19 | order 20, Paley (q = 19) | OD(100; 1, 5, 94) over ℤ₂₅, new | undecided; Table 5: 475(4) | New |
| [2652](2652) | 663 = 13 · 51 | order 52, first explicit (skew-propus family of length 13) | G-matrices of order 13 | undecided; Table 5 implies known | First explicit construction |
| [2772](2772) | 693 = 21 · 3 · 11 | order 12, Paley (q = 11) | OD(252; 1, 21, 230) over ℤ₆₃, new | undecided; Table 5: 693(4) | New |
| [5700](5700) | 1425 = 19 · 75 | order 76, first explicit (the design of [symmetric/5700](../symmetric/5700)) | G-matrices of order 19 | undecided | New, but obtainable classically |
| [7956](7956) | 1989 = 3 · 13 · 51 | order 52, as for 2652 | OD(156; 1, 3, 152) over ℤ₃₉, new | undecided | New |
| [9828](9828) | 2457 = 13 · 7 · 27 | order 28, Paley (q = 27) | OD(364; 1, 13, 350) over ℤ₉₁, new | undecided | New |
| [9900](9900) | 2475 = 25 · 99 | order 100, first explicit (skew-propus family of length 25) | G-matrices of order 25 | undecided | New, but obtainable classically |

The matrix of order 9900 is not stored; `verify/verify.py skew/9900` rebuilds it.

## Wallis's product with the amicable pairs

J. Wallis (1976; Seberry–Yamada 1992, SXI): a skew-Hadamard matrix I + S of order h and an amicable pair (A, M) of order
m give the skew-Hadamard matrix I ⊗ A + S ⊗ M of order hm. With the pairs of orders 52 and 76 above, and a pair of order
124 from the skew-propus quadruple of length 31 ([search/NOTES.md](../search/NOTES.md) §31; classically the smallest
amicable order with odd part 31 is 248), four orders 16v fall below Cati–Pasechnik's Table 5 ([NOTES.md](../NOTES.md)
§28). The skew factors are classical or published, and so are the skew-propus families of lengths 13, 19 and 31 behind
the amicable pairs (Abuzin–Balonin–Đoković–Kotsireas 2019); the one step we found nowhere in print is the passage from
such a family to the amicable pair, so the orders are labelled New, but obtainable classically.

| Order | v | Skew factor | Amicable pair | Before | Status |
|------:|---|---|---|---|---|
| [9424](9424) | 589 = 19 · 31 | 124 = 4 · 31: Paley plus three blocks under ⟨2⟩ (Szekeres 1969 has this order) | 76 | 589(5) | New, but obtainable classically |
| [9776](9776) | 611 = 13 · 47 | 188 = 4 · 47: Đoković (2008) | 52 | 611(5) | New, but obtainable classically |
| [13072](13072) | 817 = 19 · 43 | 172 = 4 · 43: Paley plus three blocks under ⟨4⟩ (Đoković 1992 has this order) | 76 | 817(5) | New, but obtainable classically |
| [14384](14384) | 899 = 29 · 31 | 116 = 4 · 29: good matrices of order 29 (Szekeres 1988) | 124 | 899(6) | New, but obtainable classically |

The matrices are not stored; `verify/verify.py` rebuilds them.

## Paley's matrix in sixteen blocks, times three

A skew-Hadamard matrix of order 3(q + 1) for every prime power q ≡ 3 (mod 4) ([PALEY_SKEW.md](PALEY_SKEW.md)): Paley's
matrix of order q + 1 = 4s is a 4 × 4 array of negacyclic s × s blocks, Frobenius pairs two blocks of one block column,
and T-matrices of order 3 in the Goethals–Seidel array then give a skew-type first block. No search is involved. The
orders below are the 52 values v = 3(q + 1)/4 < 2500 that the census lists as undecided and that were not already in
this repository, and six orders 8v and 16v (q ≡ 7 mod 8), five of which improve Cati–Pasechnik Table 5; 7656 is Known
(Wallis–Whiteman). The block form and its pairing are Yamada's (1991) and the multiplication is Seberry and Yamada's
(1990), so the new orders were obtainable classically ([status](PALEY_SKEW.md#status)). The family also contains
the census values 261, 285, 429, 465, 513, 663 and 2475, which are above with other constructions. The matrices are not
stored; `verify/verify.py skew/<order>` rebuilds them.

| Order | v | q | s | Census; Table 5 | Status |
|------:|---|---:|---:|---|---|
| [996](996) | 249 = 3 · 83 | 331 | 83 | undecided; Table 5: 249(4) | New, but obtainable classically |
| [1476](1476) | 369 = 3² · 41 | 491 | 123 | undecided; Table 5: 369(4) | New, but obtainable classically |
| [1692](1692) | 423 = 3² · 47 | 563 | 141 | undecided; Table 5: 423(4) | New, but obtainable classically |
| [1764](1764) | 441 = 3² · 7² | 587 | 147 | undecided; Table 5: 441(3) | New, but obtainable classically |
| [1800](1800) | 225 = 3² · 5² (order 8 · v) | 599 | 150 | Table 5: 225(4) | New, but obtainable classically |
| [2076](2076) | 519 = 3 · 173 | 691 | 173 | undecided; Table 5: 519(4) | New, but obtainable classically |
| [2436](2436) | 609 = 3 · 7 · 29 | 811 | 203 | undecided; Table 5: 609(3) | New, but obtainable classically |
| [2484](2484) | 621 = 3³ · 23 | 827 | 207 | undecided; Table 5: 621(3) | New, but obtainable classically |
| [2724](2724) | 681 = 3 · 227 | 907 | 227 | undecided; Table 5: 681(4) | New, but obtainable classically |
| [2916](2916) | 729 = 3⁶ | 971 | 243 | undecided; Table 5: 729(4) | New, but obtainable classically |
| [3060](3060) | 765 = 3² · 5 · 17 | 1019 | 255 | undecided; Table 5: 765(4) | New, but obtainable classically |
| [3156](3156) | 789 = 3 · 263 | 1051 | 263 | undecided; Table 5: 789(3) | New, but obtainable classically |
| [3516](3516) | 879 = 3 · 293 | 1171 | 293 | undecided; not in Table 5 | First explicit construction |
| [3564](3564) | 891 = 3⁴ · 11 | 1187 | 297 | undecided; Table 5: 891(3) | New, but obtainable classically |
| [3876](3876) | 969 = 3 · 17 · 19 | 1291 | 323 | undecided; Table 5: 969(4) | New, but obtainable classically |
| [4344](4344) | 543 = 3 · 181 (order 8 · v) | 1447 | 362 | Table 5: 543(4) | New, but obtainable classically |
| [4356](4356) | 1089 = 3² · 11² | 1451 | 363 | undecided | New, but obtainable classically |
| [4380](4380) | 1095 = 3 · 5 · 73 | 1459 | 365 | undecided | New, but obtainable classically |
| [4500](4500) | 1125 = 3² · 5³ | 1499 | 375 | undecided | New, but obtainable classically |
| [4572](4572) | 1143 = 3² · 127 | 1523 | 381 | undecided | New, but obtainable classically |
| [4596](4596) | 1149 = 3 · 383 | 1531 | 383 | undecided | New, but obtainable classically |
| [4632](4632) | 579 = 3 · 193 (order 8 · v) | 1543 | 386 | Table 5: 579(5) | New, but obtainable classically |
| [4740](4740) | 1185 = 3 · 5 · 79 | 1579 | 395 | undecided | New, but obtainable classically |
| [4860](4860) | 1215 = 3⁵ · 5 | 1619 | 405 | undecided | New, but obtainable classically |
| [4884](4884) | 1221 = 3 · 11 · 37 | 1627 | 407 | undecided | New, but obtainable classically |
| [5364](5364) | 1341 = 3² · 149 | 1787 | 447 | undecided | New, but obtainable classically |
| [5436](5436) | 1359 = 3² · 151 | 1811 | 453 | undecided | New, but obtainable classically |
| [5604](5604) | 1401 = 3 · 467 | 1867 | 467 | undecided | New, but obtainable classically |
| [5796](5796) | 1449 = 3² · 7 · 23 | 1931 | 483 | undecided | New, but obtainable classically |
| [5964](5964) | 1491 = 3 · 7 · 71 | 1987 | 497 | undecided | New, but obtainable classically |
| [6036](6036) | 1509 = 3 · 503 | 2011 | 503 | undecided | New, but obtainable classically |
| [6084](6084) | 1521 = 3² · 13² | 2027 | 507 | undecided | New, but obtainable classically |
| [6396](6396) | 1599 = 3 · 13 · 41 | 2131 | 533 | undecided | New, but obtainable classically |
| [6540](6540) | 1635 = 3 · 5 · 109 | 2179 | 545 | undecided | New, but obtainable classically |
| [6612](6612) | 1653 = 3 · 19 · 29 | 2203 | 551 | undecided | New, but obtainable classically |
| [6732](6732) | 1683 = 3² · 11 · 17 | 2243 | 561 | undecided | New, but obtainable classically |
| [6756](6756) | 1689 = 3 · 563 | 2251 | 563 | undecided | New, but obtainable classically |
| [7236](7236) | 1809 = 3³ · 67 | 2411 | 603 | undecided | New, but obtainable classically |
| [7272](7272) | 909 = 3² · 101 (order 8 · v) | 2423 | 606 | Table 5: 909(4) | New, but obtainable classically |
| [7380](7380) | 1845 = 3² · 5 · 41 | 2459 | 615 | undecided | New, but obtainable classically |
| [7620](7620) | 1905 = 3 · 5 · 127 | 2539 | 635 | undecided | New, but obtainable classically |
| [7656](7656) | 957 = 3 · 11 · 29 (order 8 · v) | 2551 | 638 | Table 5: 957(4) | Known (Wallis–Whiteman) |
| [7740](7740) | 1935 = 3² · 5 · 43 | 2579 | 645 | undecided | New, but obtainable classically |
| [8052](8052) | 2013 = 3 · 11 · 61 | 2683 | 671 | undecided | New, but obtainable classically |
| [8100](8100) | 2025 = 3⁴ · 5² | 2699 | 675 | undecided | New, but obtainable classically |
| [8196](8196) | 2049 = 3 · 683 | 2731 | 683 | undecided | New, but obtainable classically |
| [8412](8412) | 2103 = 3 · 701 | 2803 | 701 | undecided | New, but obtainable classically |
| [8532](8532) | 2133 = 3³ · 79 | 2843 | 711 | undecided | New, but obtainable classically |
| [8556](8556) | 2139 = 3 · 23 · 31 | 2851 | 713 | undecided | New, but obtainable classically |
| [8892](8892) | 2223 = 3² · 13 · 19 | 2963 | 741 | undecided | New, but obtainable classically |
| [8916](8916) | 2229 = 3 · 743 | 2971 | 743 | undecided | New, but obtainable classically |
| [9252](9252) | 2313 = 3² · 257 | 3083 | 771 | undecided | New, but obtainable classically |
| [9564](9564) | 2391 = 3 · 797 | 3187 | 797 | undecided | New, but obtainable classically |
| [9612](9612) | 2403 = 3³ · 89 | 3203 | 801 | undecided | New, but obtainable classically |
| [9780](9780) | 2445 = 3 · 5 · 163 | 3259 | 815 | undecided | New, but obtainable classically |
| [9972](9972) | 2493 = 3² · 277 | 3323 | 831 | undecided | New, but obtainable classically |
| [9996](9996) | 2499 = 3 · 7² · 17 | 3331 | 833 | undecided | New, but obtainable classically |
| [14352](14352) | 897 = 3 · 13 · 23 (order 16 · v) | 4783 | 1196 | Table 5: 897(5) | New, but obtainable classically |

## How they were found

The orders 3(q + 1) need no search (`search/paley_blocks.py skew q`). The Goethals–Seidel families come from searches
over multiplier structures: exact joins of complete block lists for small structures, exact compression chains
(`search/chain.py`, run over the census by `search/census_screen.py` and `search/screen_chain.py`), Paley-product and
cyclotomic-product configurations (`search/paley_config.py`, `search/cyclo_product.py`), the two-level sampler for
Paley products with many w-positions (`search/skewplug`, for 476, 1540, 1892 and 2156), the quaternary search over negacyclic
arrays and its Paley plug-ins (`search/negaskew`, for 1528, 2216, 2392 and 4712), low-energy k-trees on lists too large to store
(`search/cp_lowe.py`, `search/skew_lowE.py`) and, for 3676, a spin structure (`search/spin_skew.py`). The orthogonal
designs come from `search/odm_engine.py` and `search/odm_join.py`. [search/NOTES.md](../search/NOTES.md) §10–§16, §20,
§23, §31 and §87–90 describe the methods, and each README records the command, code commit, seed and caps of its
search.

`python3 verify/verify.py skew/<order>` rebuilds a matrix from its certificate and checks H Hᵀ = nI and H + Hᵀ = 2I.
