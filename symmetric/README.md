# Symmetric Hadamard matrices

A Hadamard matrix H is symmetric if H = Hᵀ. The Đoković–Seberry census (2025-05-28) lists 595 odd v < 2500 for which
no symmetric Hadamard matrix of order 4v was known; Đoković settled three of them, 189, 473 and 1139, in August 2025
(arXiv:2508.17141). The 102 orders here decide 102 more, leaving 490. For 84 of them no published construction gives a
symmetric matrix; 9 are Concurrent: 8 also follow from quaternary matrices published on Zaitsev et al.'s page of
14 October 2026, and 940 (v = 235) from the conference graph on 117 vertices that J. Y. Z. Ho reported in October 2026; 9 follow from published theorems that the census omits: 1804 and 3532 from Seberry–Balonin's
Corollary 1, and 612, 1300, 1740, 2812, 4732, 4900 and 7564 from Spence's Corollary 3 (1977). Among them is 260, the
census's smallest value v = 65.

Each directory holds one order with its certificate, its README and, for the orders found by search, the matrix.
Status labels are defined in [NOVELTY.md](../NOVELTY.md). The orders by construction:

| Family | Document or program | Orders | Count | Status |
|---|---|---|---:|---|
| Propus quadruples found by search | `search/propus.py`, `search/prow.c` | [532](532), [732](732), [876](876) | 3 | 3 New |
| Quaternary Legendre pair | `search/sym_qlp.py` | [260](260) | 1 | 1 New |
| Miyamoto's lemma in the propus array | [PROPUS_MIYAMOTO.md](PROPUS_MIYAMOTO.md) | 48 orders from [356](356) to [9892](9892) | 48 | 1 New; 47 New, but obtainable classically |
| Any Hadamard matrix in the propus array | [PROPUS_WEIGHING.md](PROPUS_WEIGHING.md) | 10 orders from [3236](3236) to [9764](9764) | 10 | 10 New |
| Yamada–Kiyasu in quaternary form | [YAMADA_KIYASU.md](YAMADA_KIYASU.md) | 28 orders from [1180](1180) to [9340](9340) | 28 | 18 New; 8 Concurrent; 2 Known |
| [[P, Q], [Q, P]] with Paley's matrices | [PERIODIC_GOLAY.md](PERIODIC_GOLAY.md) | [404](404), [596](596) | 2 | 2 New |
| Spence's recursion | [CONFERENCE_PLUG.md](CONFERENCE_PLUG.md) | [612](612), [1300](1300), [1740](1740), [4900](4900), [7564](7564) | 5 | 5 Known |
| Seberry and Balonin's Corollary 1 | `search/sym_propus_sb.py` | [940](940), [1804](1804), [3532](3532) | 3 | 1 Concurrent; 2 Known |
| Orders m(m − 1) and m(m + 3) | `search/sym_mm1.py` | [4828](4828), [5700](5700) | 2 | 1 New; 1 New, but obtainable classically |
| Total | | | 102 | 36 New; 48 New, but obtainable classically; 9 Concurrent; 9 Known |

| Order | v | Construction | Census | Status |
|------:|---|---|---|---|
| [260](260) | 65 = 5 · 13 | bordered propus quadruple over ℤ₆₄ decoded from Verma's quaternary Legendre pair of length 64 (A real, B symmetric); row sums 2, 0, 0 | undecided | New |
| [356](356) | 89 | Miyamoto's lemma in the propus array ([PROPUS_MIYAMOTO.md](PROPUS_MIYAMOTO.md)), M = Paley I, order 44 (q′ = 43) | undecided | New, but obtainable classically |
| [404](404) | 101 | a Hadamard matrix [[P, Q], [Q, P]] of order q − 1 with Paley's matrices of order q + 1 ([PERIODIC_GOLAY.md](PERIODIC_GOLAY.md)), K from a periodic Golay pair of length 50 (published; ours over ℤ₂ × ℤ₅²) | undecided | New |
| [532](532) | 133 = 7 · 19 | propus quadruple over ℤ₁₃₃, H = ⟨26⟩ of order 6; row sums 11, 5, 19 | undecided | New |
| [596](596) | 149 | a Hadamard matrix [[P, Q], [Q, P]] of order q − 1 with Paley's matrices of order q + 1 ([PERIODIC_GOLAY.md](PERIODIC_GOLAY.md)), K from Đoković–Kotsireas's periodic Golay pair of length 74 | undecided | New |
| [612](612) | 153 = 3² · 17 | symmetric conference matrix of order n in I ⊗ X + C ⊗ Y ([CONFERENCE_PLUG.md](CONFERENCE_PLUG.md)), Paley's matrix of order n = 18 (q = 17) | undecided (omitted) | Known |
| [732](732) | 183 = 3 · 61 | propus quadruple over ℤ₁₈₃, H = ⟨34⟩ of order 5 (cyclotomic product ℤ₃ × GF(61) of index 12); row sums 21, 1, 17 | undecided | New |
| [772](772) | 193 | Miyamoto's lemma in the propus array ([PROPUS_MIYAMOTO.md](PROPUS_MIYAMOTO.md)), M = H₂ ⊗ Paley I, order 48 (q′ = 47) | undecided | New, but obtainable classically |
| [876](876) | 219 = 3 · 73 | propus quadruple over ℤ₂₁₉, H = ⟨4⟩ of order 9 (cyclotomic product ℤ₃ × GF(73) of index 8); row sums 3, 17, 17 | undecided | New |
| [932](932) | 233 | Miyamoto's lemma in the propus array ([PROPUS_MIYAMOTO.md](PROPUS_MIYAMOTO.md)), M = propus, order 116 (Di Matteo–Đoković–Kotsireas) | undecided | New, but obtainable classically |
| [940](940) | 235 = 5 · 47 | Seberry–Balonin's Corollary 1 with q = 233 and the core of a symmetric conference matrix of order 118 found here by an exact cover (`search/conference118`) | undecided | Concurrent (J. Y. Z. Ho's conference graph on 117 vertices) |
| [964](964) | 241 | Miyamoto's lemma in the propus array ([PROPUS_MIYAMOTO.md](PROPUS_MIYAMOTO.md)), M = H₂ ⊗ Paley I, order 60 (q′ = 59) | undecided | New, but obtainable classically |
| [1028](1028) | 257 | Miyamoto's lemma in the propus array ([PROPUS_MIYAMOTO.md](PROPUS_MIYAMOTO.md)), M = Paley I, order 128 (q′ = 127) | undecided | New, but obtainable classically |
| [1180](1180) | 295 = 5 · 59 | Yamada–Kiyasu in symmetric quaternary form ([YAMADA_KIYASU.md](YAMADA_KIYASU.md)), skew input of order 148: Goethals–Seidel over ℤ₃₇ | undecided | Concurrent |
| [1276](1276) | 319 = 11 · 29 | Yamada–Kiyasu in symmetric quaternary form ([YAMADA_KIYASU.md](YAMADA_KIYASU.md)), skew input of order 160: Paley (p = 19), doubled three times | undecided | Concurrent |
| [1300](1300) | 325 = 5² · 13 | symmetric conference matrix of order n in I ⊗ X + C ⊗ Y ([CONFERENCE_PLUG.md](CONFERENCE_PLUG.md)), Paley's matrix of order n = 26 (q = 25) | undecided (omitted) | Known |
| [1412](1412) | 353 | Miyamoto's lemma in the propus array ([PROPUS_MIYAMOTO.md](PROPUS_MIYAMOTO.md)), M = H₂^⊗2 ⊗ Paley I, order 44 (q′ = 43) | undecided | New, but obtainable classically |
| [1444](1444) | 361 = 19² | Miyamoto's lemma in the propus array ([PROPUS_MIYAMOTO.md](PROPUS_MIYAMOTO.md)), M = Paley I, order 180 (q′ = 179) | undecided | New, but obtainable classically |
| [1564](1564) | 391 = 17 · 23 | Yamada–Kiyasu in symmetric quaternary form ([YAMADA_KIYASU.md](YAMADA_KIYASU.md)), skew input of order 196: Goethals–Seidel over ℤ₄₉ | undecided | Concurrent |
| [1604](1604) | 401 | Miyamoto's lemma in the propus array ([PROPUS_MIYAMOTO.md](PROPUS_MIYAMOTO.md)), M = Paley I, order 200 (q′ = 199) | undecided | New, but obtainable classically |
| [1636](1636) | 409 | Miyamoto's lemma in the propus array ([PROPUS_MIYAMOTO.md](PROPUS_MIYAMOTO.md)), M = Paley II, order 204 (q′ = 101) | undecided | New, but obtainable classically |
| [1732](1732) | 433 | Miyamoto's lemma in the propus array ([PROPUS_MIYAMOTO.md](PROPUS_MIYAMOTO.md)), M = H₂ ⊗ Paley I, order 108 (q′ = 107) | undecided | New, but obtainable classically |
| [1740](1740) | 435 = 3 · 5 · 29 | symmetric conference matrix of order n in I ⊗ X + C ⊗ Y ([CONFERENCE_PLUG.md](CONFERENCE_PLUG.md)), Paley's matrix of order n = 30 (q = 29) | undecided (omitted) | Known |
| [1796](1796) | 449 | Miyamoto's lemma in the propus array ([PROPUS_MIYAMOTO.md](PROPUS_MIYAMOTO.md)), M = Paley I, order 224 (q′ = 223) | undecided | New, but obtainable classically |
| [1804](1804) | 451 = 11 · 41 | Seberry–Balonin's Corollary 1 with q = 449 and the core of a symmetric conference matrix of order 226 (from a skew-Hadamard matrix of order 16) | undecided (omitted) | Known |
| [1828](1828) | 457 | Miyamoto's lemma in the propus array ([PROPUS_MIYAMOTO.md](PROPUS_MIYAMOTO.md)), M = Paley I, order 228 (q′ = 227) | undecided | New, but obtainable classically |
| [1852](1852) | 463 | Yamada–Kiyasu in symmetric quaternary form ([YAMADA_KIYASU.md](YAMADA_KIYASU.md)), skew input of order 232: good matrices of order 29 (Szekeres), doubled once | undecided | Concurrent |
| [2116](2116) | 529 = 23² | Miyamoto's lemma in the propus array ([PROPUS_MIYAMOTO.md](PROPUS_MIYAMOTO.md)), M = Paley I, order 264 (q′ = 263) | undecided | New, but obtainable classically |
| [2276](2276) | 569 | Miyamoto's lemma in the propus array ([PROPUS_MIYAMOTO.md](PROPUS_MIYAMOTO.md)), M = Paley I, order 284 (q′ = 283) | undecided | New, but obtainable classically |
| [2564](2564) | 641 | Miyamoto's lemma in the propus array ([PROPUS_MIYAMOTO.md](PROPUS_MIYAMOTO.md)), M = H₂^⊗2 ⊗ Paley I, order 80 (q′ = 79) | undecided | New, but obtainable classically |
| [2620](2620) | 655 = 5 · 131 | Yamada–Kiyasu in symmetric quaternary form ([YAMADA_KIYASU.md](YAMADA_KIYASU.md)), skew input of order 328: Paley (p = 163), doubled once | undecided | Concurrent |
| [2692](2692) | 673 | Miyamoto's lemma in the propus array ([PROPUS_MIYAMOTO.md](PROPUS_MIYAMOTO.md)), M = H₂ ⊗ Paley I, order 168 (q′ = 167) | undecided | New, but obtainable classically |
| [2716](2716) | 679 = 7 · 97 | Yamada–Kiyasu in symmetric quaternary form ([YAMADA_KIYASU.md](YAMADA_KIYASU.md)), skew input of order 340: Whiteman, GF(13²) | undecided | Concurrent |
| [2812](2812) | 703 = 19 · 37 | Yamada–Kiyasu in symmetric quaternary form ([YAMADA_KIYASU.md](YAMADA_KIYASU.md)), skew input of order 352: Paley (p = 43), doubled three times | undecided (omitted) | Known |
| [3076](3076) | 769 | Miyamoto's lemma in the propus array ([PROPUS_MIYAMOTO.md](PROPUS_MIYAMOTO.md)), M = Paley I, order 384 (q′ = 383) | undecided | New, but obtainable classically |
| [3236](3236) | 809 | any Hadamard matrix of order (q − 1)/2 in the propus array ([PROPUS_WEIGHING.md](PROPUS_WEIGHING.md)), K = Cooper–Wallis, T-matrices of order 101, order 404 | undecided | New |
| [3292](3292) | 823 | Yamada–Kiyasu in symmetric quaternary form ([YAMADA_KIYASU.md](YAMADA_KIYASU.md)), skew input of order 412: Goethals–Seidel over ℤ₁₀₃ (Đoković 1994) | undecided | Concurrent |
| [3428](3428) | 857 | any Hadamard matrix of order (q − 1)/2 in the propus array ([PROPUS_WEIGHING.md](PROPUS_WEIGHING.md)), K = Cooper–Wallis, T-matrices of order 107, order 428 | undecided | New |
| [3524](3524) | 881 | Miyamoto's lemma in the propus array ([PROPUS_MIYAMOTO.md](PROPUS_MIYAMOTO.md)), M = Paley I, order 440 (q′ = 439) | undecided | New, but obtainable classically |
| [3532](3532) | 883 | Seberry–Balonin's Corollary 1 with q = 881 and the core of Mathon's conference matrix of order 442 | undecided (omitted) | Known |
| [3716](3716) | 929 | Miyamoto's lemma in the propus array ([PROPUS_MIYAMOTO.md](PROPUS_MIYAMOTO.md)), M = Paley I, order 464 (q′ = 463) | undecided | New, but obtainable classically |
| [3772](3772) | 943 = 23 · 41 | Yamada–Kiyasu in symmetric quaternary form ([YAMADA_KIYASU.md](YAMADA_KIYASU.md)), skew input of order 472: Goethals–Seidel over ℤ₅₉ (Fletcher–Koukouvinos–Seberry), doubled once | undecided | Concurrent |
| [3812](3812) | 953 | any Hadamard matrix of order (q − 1)/2 in the propus array ([PROPUS_WEIGHING.md](PROPUS_WEIGHING.md)), K = Cooper–Wallis, T-matrices of order 17 and Turyn's Williamson matrices of order 7, order 476 | undecided | New |
| [3844](3844) | 961 = 31² | Miyamoto's lemma in the propus array ([PROPUS_MIYAMOTO.md](PROPUS_MIYAMOTO.md)), M = Paley I, order 480 (q′ = 479) | undecided | New, but obtainable classically |
| [4132](4132) | 1033 | Miyamoto's lemma in the propus array ([PROPUS_MIYAMOTO.md](PROPUS_MIYAMOTO.md)), M = Paley II, order 516 (q′ = 257) | undecided | New, but obtainable classically |
| [4196](4196) | 1049 | Miyamoto's lemma in the propus array ([PROPUS_MIYAMOTO.md](PROPUS_MIYAMOTO.md)), M = Paley I, order 524 (q′ = 523) | undecided | New, but obtainable classically |
| [4252](4252) | 1063 | Yamada–Kiyasu in symmetric quaternary form ([YAMADA_KIYASU.md](YAMADA_KIYASU.md)), skew input of order 532: Goethals–Seidel over ℤ₁₃₃ | undecided | New |
| [4388](4388) | 1097 | Miyamoto's lemma in the propus array ([PROPUS_MIYAMOTO.md](PROPUS_MIYAMOTO.md)), M = Paley I, order 548 (q′ = 547) | undecided | New, but obtainable classically |
| [4516](4516) | 1129 | Miyamoto's lemma in the propus array ([PROPUS_MIYAMOTO.md](PROPUS_MIYAMOTO.md)), M = Paley I, order 564 (q′ = 563) | undecided | New, but obtainable classically |
| [4612](4612) | 1153 | Miyamoto's lemma in the propus array ([PROPUS_MIYAMOTO.md](PROPUS_MIYAMOTO.md)), M = H₂^⊗3 ⊗ Paley I, order 72 (q′ = 71) | undecided | New, but obtainable classically |
| [4732](4732) | 1183 = 7 · 13² | Yamada–Kiyasu in symmetric quaternary form ([YAMADA_KIYASU.md](YAMADA_KIYASU.md)), skew input of order 592: Goethals–Seidel over ℤ₃₇, doubled twice | undecided (omitted) | Known |
| [4772](4772) | 1193 | any Hadamard matrix of order (q − 1)/2 in the propus array ([PROPUS_WEIGHING.md](PROPUS_WEIGHING.md)), K = Miyamoto (q′ = 149) with T-matrices of order 37, order 596 | undecided | New |
| [4828](4828) | 1207 = 17 · 71 | m(m + 3) with m = 68: Paley's symmetric OD(68; 1, 67), x ↦ J − 2I, y ↦ the core of a symmetric Hadamard matrix of order 72 | undecided | New |
| [4868](4868) | 1217 | Miyamoto's lemma in the propus array ([PROPUS_MIYAMOTO.md](PROPUS_MIYAMOTO.md)), M = Paley I, order 608 (q′ = 607) | undecided | New, but obtainable classically |
| [4900](4900) | 1225 = 5² · 7² | symmetric conference matrix of order n in I ⊗ X + C ⊗ Y ([CONFERENCE_PLUG.md](CONFERENCE_PLUG.md)), Paley's matrix of order n = 50 (q = 49) | undecided (omitted) | Known |
| [4924](4924) | 1231 | Yamada–Kiyasu in symmetric quaternary form ([YAMADA_KIYASU.md](YAMADA_KIYASU.md)), skew input of order 616: Paley (p = 307), doubled once | undecided | New |
| [4996](4996) | 1249 | Miyamoto's lemma in the propus array ([PROPUS_MIYAMOTO.md](PROPUS_MIYAMOTO.md)), M = H₂ ⊗ Paley I, order 312 (q′ = 311) | undecided | New, but obtainable classically |
| [5156](5156) | 1289 | Miyamoto's lemma in the propus array ([PROPUS_MIYAMOTO.md](PROPUS_MIYAMOTO.md)), M = Paley I, order 644 (q′ = 643) | undecided | New, but obtainable classically |
| [5212](5212) | 1303 | Yamada–Kiyasu in symmetric quaternary form ([YAMADA_KIYASU.md](YAMADA_KIYASU.md)), skew input of order 652: Goethals–Seidel over ℤ₁₆₃ (Đoković 1992) | undecided | New |
| [5284](5284) | 1321 | Miyamoto's lemma in the propus array ([PROPUS_MIYAMOTO.md](PROPUS_MIYAMOTO.md)), M = Paley I, order 660 (q′ = 659) | undecided | New, but obtainable classically |
| [5476](5476) | 1369 = 37² | Miyamoto's lemma in the propus array ([PROPUS_MIYAMOTO.md](PROPUS_MIYAMOTO.md)), M = Paley I, order 684 (q′ = 683) | undecided | New, but obtainable classically |
| [5636](5636) | 1409 | Miyamoto's lemma in the propus array ([PROPUS_MIYAMOTO.md](PROPUS_MIYAMOTO.md)), M = H₂^⊗4 ⊗ Paley I, order 44 (q′ = 43) | undecided | New, but obtainable classically |
| [5700](5700) | 1425 = 3 · 5² · 19 | m(m − 1) with m = 76: a skew-propus quadruple of length 19 (row sums 5, 5, 1, invariant under the multiplier 7) gives a symmetric OD(76; 1, 75); x ↦ J, y ↦ the core of P + N | undecided | New, but obtainable classically |
| [5732](5732) | 1433 | any Hadamard matrix of order (q − 1)/2 in the propus array ([PROPUS_WEIGHING.md](PROPUS_WEIGHING.md)), K = Alpöge–Voinov–Reynolds-Haertle's matrix of order 716 (2026), order 716 | undecided | New |
| [5956](5956) | 1489 | Miyamoto's lemma in the propus array ([PROPUS_MIYAMOTO.md](PROPUS_MIYAMOTO.md)), M = Paley I, order 744 (q′ = 743) | undecided | New, but obtainable classically |
| [5980](5980) | 1495 = 5 · 13 · 23 | Yamada–Kiyasu in symmetric quaternary form ([YAMADA_KIYASU.md](YAMADA_KIYASU.md)), skew input of order 748: Szekeres, GF(373) | undecided | New |
| [6396](6396) | 1599 = 3 · 13 · 41 | Yamada–Kiyasu in symmetric quaternary form ([YAMADA_KIYASU.md](YAMADA_KIYASU.md)), skew input of order 800: Paley (p = 199), doubled twice | undecided | New |
| [6404](6404) | 1601 | Miyamoto's lemma in the propus array ([PROPUS_MIYAMOTO.md](PROPUS_MIYAMOTO.md)), M = H₂^⊗2 ⊗ Paley I, order 200 (q′ = 199) | undecided | New, but obtainable classically |
| [6556](6556) | 1639 = 11 · 149 | Yamada–Kiyasu in symmetric quaternary form ([YAMADA_KIYASU.md](YAMADA_KIYASU.md)), skew input of order 820: our skew matrix of order 820 | undecided | New |
| [6684](6684) | 1671 = 3 · 557 | Yamada–Kiyasu in symmetric quaternary form ([YAMADA_KIYASU.md](YAMADA_KIYASU.md)), skew input of order 836: our skew matrix of order 836 | undecided | New |
| [6788](6788) | 1697 | Miyamoto's lemma in the propus array ([PROPUS_MIYAMOTO.md](PROPUS_MIYAMOTO.md)), M = H₂^⊗2 ⊗ Paley I, order 212 (q′ = 211) | undecided | New, but obtainable classically |
| [7012](7012) | 1753 | Miyamoto's lemma in the propus array ([PROPUS_MIYAMOTO.md](PROPUS_MIYAMOTO.md)), M = our symmetric matrix of order 876 | undecided | New |
| [7108](7108) | 1777 | Miyamoto's lemma in the propus array ([PROPUS_MIYAMOTO.md](PROPUS_MIYAMOTO.md)), M = Paley I, order 888 (q′ = 887) | undecided | New, but obtainable classically |
| [7204](7204) | 1801 | Miyamoto's lemma in the propus array ([PROPUS_MIYAMOTO.md](PROPUS_MIYAMOTO.md)), M = Paley II, order 900 (q′ = 449) | undecided | New, but obtainable classically |
| [7492](7492) | 1873 | Miyamoto's lemma in the propus array ([PROPUS_MIYAMOTO.md](PROPUS_MIYAMOTO.md)), M = H₂ ⊗ Paley I, order 468 (q′ = 467) | undecided | New, but obtainable classically |
| [7556](7556) | 1889 | Miyamoto's lemma in the propus array ([PROPUS_MIYAMOTO.md](PROPUS_MIYAMOTO.md)), M = H₂^⊗2 ⊗ propus, order 236 (Đoković) | undecided | New, but obtainable classically |
| [7564](7564) | 1891 = 31 · 61 | symmetric conference matrix of order n in I ⊗ X + C ⊗ Y ([CONFERENCE_PLUG.md](CONFERENCE_PLUG.md)), Paley's matrix of order n = 62 (q = 61) | undecided (omitted) | Known |
| [7612](7612) | 1903 = 11 · 173 | Yamada–Kiyasu in symmetric quaternary form ([YAMADA_KIYASU.md](YAMADA_KIYASU.md)), skew input of order 952: our skew matrix of order 476, doubled once | undecided | New |
| [7652](7652) | 1913 | any Hadamard matrix of order (q − 1)/2 in the propus array ([PROPUS_WEIGHING.md](PROPUS_WEIGHING.md)), K = Đoković–Kotsireas's skew matrix of order 956, order 956 | undecided | New |
| [7804](7804) | 1951 | Yamada–Kiyasu in symmetric quaternary form ([YAMADA_KIYASU.md](YAMADA_KIYASU.md)), skew input of order 976: Paley (p = 487), doubled once | undecided | New |
| [7900](7900) | 1975 = 5² · 79 | Yamada–Kiyasu in symmetric quaternary form ([YAMADA_KIYASU.md](YAMADA_KIYASU.md)), skew input of order 988: Goethals–Seidel over ℤ₂₄₇ (Đoković 2008) | undecided | New |
| [7972](7972) | 1993 | any Hadamard matrix of order (q − 1)/2 in the propus array ([PROPUS_WEIGHING.md](PROPUS_WEIGHING.md)), K = our skew matrix of order 996, order 996 | undecided | New |
| [7996](7996) | 1999 | Yamada–Kiyasu in symmetric quaternary form ([YAMADA_KIYASU.md](YAMADA_KIYASU.md)), skew input of order 1000: Paley (p = 499), doubled once | undecided | New |
| [8068](8068) | 2017 | Miyamoto's lemma in the propus array ([PROPUS_MIYAMOTO.md](PROPUS_MIYAMOTO.md)), M = H₂ ⊗ Paley I, order 504 (q′ = 503) | undecided | New, but obtainable classically |
| [8284](8284) | 2071 = 19 · 109 | Yamada–Kiyasu in symmetric quaternary form ([YAMADA_KIYASU.md](YAMADA_KIYASU.md)), skew input of order 1036: our skew matrix of order 1036 | undecided | New |
| [8324](8324) | 2081 | Miyamoto's lemma in the propus array ([PROPUS_MIYAMOTO.md](PROPUS_MIYAMOTO.md)), M = Paley I, order 1040 (q′ = 1039) | undecided | New, but obtainable classically |
| [8452](8452) | 2113 | Miyamoto's lemma in the propus array ([PROPUS_MIYAMOTO.md](PROPUS_MIYAMOTO.md)), M = H₂^⊗2 ⊗ Paley I, order 264 (q′ = 263) | undecided | New, but obtainable classically |
| [8516](8516) | 2129 | Miyamoto's lemma in the propus array ([PROPUS_MIYAMOTO.md](PROPUS_MIYAMOTO.md)), M = Paley I, order 1064 (q′ = 1063) | undecided | New, but obtainable classically |
| [8572](8572) | 2143 | Yamada–Kiyasu in symmetric quaternary form ([YAMADA_KIYASU.md](YAMADA_KIYASU.md)), skew input of order 1072: Goethals–Seidel over ℤ₆₇ (Đoković 1992), doubled twice | undecided | New |
| [8612](8612) | 2153 | any Hadamard matrix of order (q − 1)/2 in the propus array ([PROPUS_WEIGHING.md](PROPUS_WEIGHING.md)), K = Miyamoto (q′ = 269) with T-matrices of order 67, order 1076 | undecided | New |
| [8644](8644) | 2161 | Miyamoto's lemma in the propus array ([PROPUS_MIYAMOTO.md](PROPUS_MIYAMOTO.md)), M = H₂ ⊗ Paley II, order 540 (q′ = 269) | undecided | New, but obtainable classically |
| [8836](8836) | 2209 = 47² | Miyamoto's lemma in the propus array ([PROPUS_MIYAMOTO.md](PROPUS_MIYAMOTO.md)), M = Paley I, order 1104 (q′ = 1103) | undecided | New, but obtainable classically |
| [8892](8892) | 2223 = 3² · 13 · 19 | Yamada–Kiyasu in symmetric quaternary form ([YAMADA_KIYASU.md](YAMADA_KIYASU.md)), skew input of order 1112: Szekeres, GF(277), doubled once | undecided | New |
| [8956](8956) | 2239 | Yamada–Kiyasu in symmetric quaternary form ([YAMADA_KIYASU.md](YAMADA_KIYASU.md)), skew input of order 1120: Paley (p = 139), doubled three times | undecided | New |
| [9084](9084) | 2271 = 3 · 757 | Yamada–Kiyasu in symmetric quaternary form ([YAMADA_KIYASU.md](YAMADA_KIYASU.md)), skew input of order 1136: Paley (p = 283), doubled twice | undecided | New |
| [9180](9180) | 2295 = 3³ · 5 · 17 | Yamada–Kiyasu in symmetric quaternary form ([YAMADA_KIYASU.md](YAMADA_KIYASU.md)), skew input of order 1148: our skew matrix of order 1148 | undecided | New |
| [9340](9340) | 2335 = 5 · 467 | Yamada–Kiyasu in symmetric quaternary form ([YAMADA_KIYASU.md](YAMADA_KIYASU.md)), skew input of order 1168: Goethals–Seidel over ℤ₇₃, doubled twice | undecided | New |
| [9508](9508) | 2377 | Miyamoto's lemma in the propus array ([PROPUS_MIYAMOTO.md](PROPUS_MIYAMOTO.md)), M = Paley I, order 1188 (q′ = 1187) | undecided | New, but obtainable classically |
| [9572](9572) | 2393 | any Hadamard matrix of order (q − 1)/2 in the propus array ([PROPUS_WEIGHING.md](PROPUS_WEIGHING.md)), K = Cooper–Wallis, T-matrices of order 23 and Turyn's Williamson matrices of order 13, order 1196 | undecided | New |
| [9668](9668) | 2417 | Miyamoto's lemma in the propus array ([PROPUS_MIYAMOTO.md](PROPUS_MIYAMOTO.md)), M = H₂ ⊗ propus, order 604 (Balonin–Đoković) | undecided | New, but obtainable classically |
| [9764](9764) | 2441 | any Hadamard matrix of order (q − 1)/2 in the propus array ([PROPUS_WEIGHING.md](PROPUS_WEIGHING.md)), K = Cooper–Wallis, T-matrices of order 5 and Turyn's Williamson matrices of order 61, order 1220 | undecided | New |
| [9892](9892) | 2473 | Miyamoto's lemma in the propus array ([PROPUS_MIYAMOTO.md](PROPUS_MIYAMOTO.md)), M = Paley II, order 1236 (q′ = 617) | undecided | New, but obtainable classically |

## Constructions

**Propus quadruples.** Circulant ±1 matrices (A, B, B, D) of order v with A symmetric and
A Aᵀ + 2 B Bᵀ + D Dᵀ = 4v I fill a symmetric array of order 4v (Seberry and Balonin, Australas. J. Combin. 69 (2017);
Balonin and Đoković, arXiv:1803.08787). The quadruples here are invariant under a multiplier group H and were found by
exact joins (`search/propus.py`, [search/NOTES.md](../search/NOTES.md) §17), through compressions for 876
(`search/propus_chain.py`, §25 and §28), and by the row search `search/prow.c` for 732 (§30 and §31).

**Quaternary Legendre pair.** Verma's quaternary Legendre pair (A, B) of length 64 (arXiv:2609.25573) has A real and B
symmetric. Writing B = ((1 + i)p + (1 − i)q)/2 with p, q = ±1 gives a bordered propus quadruple (p, A, A, q) over ℤ₆₄,
which fills the bordered propus array (`search/sym_qlp.py`, `verify.py` construction `propus-bordered`,
[NOTES.md](../NOTES.md) §32): order 260 ([260](260)).

**Miyamoto's lemma in the propus array.** For a prime power q ≡ 1 (mod 8) and a symmetric Hadamard matrix M of order
(q − 1)/2, Seberry and Yamada's Williamson-type matrices of order q (Corollary 29(ii), 1990) have two equal blocks and
fill Seberry and Balonin's propus array, which gives a symmetric Hadamard matrix of order 4q
([PROPUS_MIYAMOTO.md](PROPUS_MIYAMOTO.md), `search/sym_propus_miyamoto.py`). This settles 48 census values; for 47 of
them every ingredient was published before 2026, and 7012 = 4 · 1753 uses our matrix of order 876.

**Any Hadamard matrix in the propus array.** In the construction above only the square and the support of the
V-part matter, so it can be replaced by a symmetric weighing matrix of order 2(q − 1) and weight q with the same
support, and one is built from any Hadamard matrix of order (q − 1)/2: symmetric Hadamard matrices of order 4q for
every prime power q ≡ 1 (mod 8) with a Hadamard matrix of order (q − 1)/2 ([PROPUS_WEIGHING.md](PROPUS_WEIGHING.md),
`search/sym_weighing.py`, `verify.py` construction `symmetric-propus-weighing`). This settles the ten census values
q ≡ 1 (mod 8) that the symmetric input left open, among them 1433 with Alpöge, Voinov and Reynolds-Haertle's matrix of
order 716 (2026).

**Yamada–Kiyasu in quaternary form.** Yamada's Lemma 26 data for his (Kiyasu's) Hadamard matrices of order 4(q + 2)
(q ≡ 5 mod 8 a prime power, a skew-Hadamard matrix of order (q + 3)/2) form a symmetric quaternary Hadamard matrix of
order 2(q + 2) for every skew input, and Turyn's map [[X + Y, X − Y], [X − Y, −X − Y]] makes it a symmetric Hadamard
matrix of order 4(q + 2) ([YAMADA_KIYASU.md](YAMADA_KIYASU.md), `search/sym_yk.py`, `verify.py` construction
`symmetric-yamada-kiyasu`). This settles 28 census values, two of them (2812, 4732) also covered by Spence's
corollary; six more need a skew order that is undecided.

**[[P, Q], [Q, P]] with Paley's matrices.** For a prime power q ≡ 1 (mod 4) and a Hadamard matrix K = [[P, Q], [Q, P]]
of order q − 1, for example one built from a periodic Golay pair of length (q − 1)/2, Paley's conference matrix of
order q + 1 gives a symmetric Hadamard matrix of order 4q with an involutory automorphism with 2q − 2 fixed points
([PERIODIC_GOLAY.md](PERIODIC_GOLAY.md), `search/sym_involution.py`, `verify.py` construction `symmetric-involution`).
This settles v = 101 and 149 (q = 101, 149).

**Spence's recursion.** From a symmetric conference matrix C of order n, H = I ⊗ X + C ⊗ Y is a symmetric Hadamard
matrix of order 2n(n − 1) ([CONFERENCE_PLUG.md](CONFERENCE_PLUG.md), `search/sym_conference_plug.py`, `verify.py`
construction `symmetric-conference-plug`). For Paley's matrices (n = q + 1) this is the case r = 1 of Spence's
Corollary 3 (J. Austral. Math. Soc. Ser. A 24 (1977) 252–256), which gives order 2qʳ(q + 1) for every r ≥ 0. The five
orders here, 612, 1300, 1740, 4900 and 7564 (q = 17, 25, 29, 49, 61), are Known: the census omits the corollary.

**Seberry and Balonin's Corollary 1.** For a prime power q ≡ 1 (mod 4) such that (q + 1)/2 is a prime power or the
order of the core of a symmetric conference matrix, Lemma 26 with Paley's or that core and Yamada's circulants gives
Williamson-type matrices of order q + 2 with a repeated block, hence a symmetric Hadamard matrix of order 4(q + 2).
The census lists neither of the orders that the cores of orders 225 and 441 give, 1804 and 3532
(`search/sym_propus_sb.py`, `verify.py` construction `symmetric-propus-sb`). With the symmetric conference matrix of order 118 found
here (`search/conference118`, [search/NOTES.md](../search/NOTES.md) §95) it gives 940 (q = 233), Concurrent with the
conference graph on 117 vertices that J. Y. Z. Ho reported in October 2026.

**Order m(m − 1).** A skew-Hadamard matrix S = I + C of order m with an involutive signed permutation P such that
PCP = −C is the same as a symmetric orthogonal design xP + yN of type (1, m − 1), with N = CP. Replacing x by J and y
by the core D of the symmetric Hadamard matrix P + N gives the symmetric Hadamard matrix P ⊗ J + N ⊗ D of order
m(m − 1) ([search/NOTES.md](../search/NOTES.md) §31, `search/sym_mm1.py`). Đoković's theorem (arXiv:2508.17141) is the
case m − 1 = q ≡ 3 (mod 8) a prime power. The design can also come from a skew-propus quadruple of any odd length n
(m = 4n), found here by exhaustive search, or from a Legendre pair with one skew and one symmetric block.

**Order m(m + 3).** Replacing x by J − 2I instead of J, with the core of a symmetric Hadamard matrix of order m + 4,
gives order m(m + 3) ([search/NOTES.md](../search/NOTES.md) §31). This is the symmetric form of Wallis's construction
of Hadamard and Williamson-type matrices of order m(m + 3) (1970, 1975).

`python3 verify/verify.py symmetric/<order>` rebuilds a matrix from its certificate and checks H Hᵀ = nI and H = Hᵀ.
