# Unresolved orders

This page covers the 153 entries of Cati–Pasechnik Table 4 (n ≤ 2999) that are
still unresolved. The other 42 were settled as follows:

- 12 by Alpöge et al.,
- 18 by other public constructions from 2026 (Schneider; JD-Jones-ASES/Hadamard-T
  and Hadamard-M; Schumer),
- 9 by Đoković in 2010,
- 3 in this repository.

**Columns.** For each order the table lists the cheapest structure of each kind
found by `search/triage.py` that is expected to contain at least ten solutions:

- multiplier subgroups of ℤ_n* for Goethals–Seidel,
- groups of order n − 1 for Spence's bordered construction.

log₁₀ E is the expected number of solutions in the structure. log₁₀ W is the
predicted matching work per solution with lists of 10⁹ candidates (see
[METHOD.md](METHOD.md)). The searches so far needed about a tenth of W. One
machine performs about 10¹¹ list operations per hour, so orders up to
log₁₀ W ≈ 13 are within reach.

**Smaller structures.** Smaller structures, which the Gaussian model does not
describe reliably, were searched exhaustively:

- 1,319 structures in total: Spence structures with pools of at most 4 · 10⁴
  sequences, and Goethals–Seidel structures with at most 18 cosets.
- For primes above about 400 there are no Goethals–Seidel structures of the latter
  kind, by the √p bound in METHOD.md.
- 1,159 of these were searched completely and contain no solution. 80 were searched
  partially, and 80 have more than 22 orbits and were not searched.

**Route column.** The last column names a recursive theorem that would settle the
order from a smaller object, where one applies:

- Miyamoto: an Hadamard matrix of order n − 1.
- Kiyasu: a skew Hadamard matrix of order (n + 1)/2.
- Yamada: a symmetric conference matrix of order (n + 1)/2.

A conference matrix of that order cannot exist when (n − 1)/2 is not a sum of two
squares.

| 4n | n | GS: \|H\| | cosets | classes | log₁₀ E | log₁₀ W | Spence: group | log₁₀ W | Route |
|---:|---:|---:|---:|---:|---:|---:|---|---:|---|
| 2148 | 537 | – | – | – | – | – | ℤ8 × ℤ67 | 24.1 |  |
| 2284 | 571 | – | – | – | – | – | ℤ2 × ℤ3 × ℤ5 × ℤ19 | 33.4 | Yamada from a symmetric conference matrix of order 286 (cannot exist) |
| 2292 | 573 | 10 | 78 | 39 | 20.5 | 30.0 | ℤ2xℤ2 × ℤ11 × ℤ13 | 25.3 |  |
| 2396 | 599 | – | – | – | – | – | ℤ2 × ℤ13 × ℤ23 | 36.4 |  |
| 2572 | 643 | – | – | – | – | – | – | – | Yamada from a symmetric conference matrix of order 322 (cannot exist) |
| 2588 | 647 | – | – | – | – | – | ℤ2 × ℤ17 × ℤ19 | 35.7 |  |
| 2636 | 659 | 7 | 95 | 47 | 39.8 | 26.0 | ℤ2 × ℤ7 × ℤ47 | 19.4 |  |
| 2676 | 669 | – | – | – | – | – | – | – |  |
| 2876 | 719 | – | – | – | – | – | – | – |  |
| 2956 | 739 | 9 | 83 | 41 | 33.5 | 22.7 | ℤ2 × ℤ3xℤ3 × ℤ41 | 42.4 |  |
| 3004 | 751 | – | – | – | – | – | ℤ2 × ℤ3 × ℤ5xℤ5xℤ5 | 46.0 |  |
| 3156 | 789 | – | – | – | – | – | – | – |  |
| 3356 | 839 | – | – | – | – | – | ℤ2 × ℤ419 | 15.0 |  |
| 3436 | 859 | 13 | 67 | 33 | 25.7 | 17.2 | ℤ2 × ℤ3 × ℤ11 × ℤ13 | 61.8 | Yamada from a symmetric conference matrix of order 430 (cannot exist) |
| 3452 | 863 | – | – | – | – | – | – | – |  |
| 3628 | 907 | – | – | – | – | – | ℤ2 × ℤ3 × ℤ151 | 33.3 |  |
| 3676 | 919 | 17 | 55 | 27 | 20.5 | 12.7 | ℤ2 × ℤ3xℤ3xℤ3 × ℤ17 | 31.8 |  |
| 3732 | 933 | – | – | – | – | – | – | – |  |
| 3788 | 947 | – | – | – | – | – | ℤ2 × ℤ11 × ℤ43 | 37.4 |  |
| 3820 | 955 | 19 | 55 | 27 | 2.9 | 28.8 | ℤ2 × ℤ3xℤ3 × ℤ53 | 52.8 | Yamada from a symmetric conference matrix of order 478 |
| 3884 | 971 | – | – | – | – | – | – | – |  |
| 3964 | 991 | – | – | – | – | – | – | – |  |
| 4076 | 1019 | – | – | – | – | – | – | – |  |
| 4124 | 1031 | – | – | – | – | – | ℤ2 × ℤ5 × ℤ103 | 40.8 |  |
| 4156 | 1039 | – | – | – | – | – | – | – |  |
| 4204 | 1051 | – | – | – | – | – | ℤ2 × ℤ3 × ℤ5xℤ5 × ℤ7 | 50.4 | Yamada from a symmetric conference matrix of order 526 (cannot exist) |
| 4348 | 1087 | – | – | – | – | – | ℤ2 × ℤ3 × ℤ181 | 41.9 |  |
| 4412 | 1103 | 19 | 59 | 29 | 21.1 | 15.5 | ℤ2 × ℤ19 × ℤ29 | 40.1 |  |
| 4460 | 1115 | – | – | – | – | – | – | – |  |
| 4676 | 1169 | – | – | – | – | – | ℤ16 × ℤ73 | 50.2 |  |
| 4748 | 1187 | – | – | – | – | – | ℤ2 × ℤ593 | 15.2 |  |
| 4892 | 1223 | 13 | 95 | 47 | 33.7 | 32.4 | ℤ2 × ℤ13 × ℤ47 | 48.8 |  |
| 5020 | 1255 | – | – | – | – | – | – | – |  |
| 5028 | 1257 | 22 | 78 | 39 | 12.7 | 37.8 | ℤ8 × ℤ157 | 58.9 |  |
| 5036 | 1259 | – | – | – | – | – | – | – |  |
| 5132 | 1283 | – | – | – | – | – | – | – |  |
| 5164 | 1291 | 15 | 87 | 43 | 30.2 | 29.3 | ℤ2 × ℤ3 × ℤ5 × ℤ43 | 57.1 | Yamada from a symmetric conference matrix of order 646 (cannot exist) |
| 5260 | 1315 | – | – | – | – | – | ℤ2 × ℤ9 × ℤ73 | 48.9 |  |
| 5276 | 1319 | – | – | – | – | – | – | – |  |
| 5308 | 1327 | – | – | – | – | – | ℤ2 × ℤ3 × ℤ13 × ℤ17 | 58.6 |  |
| 5468 | 1367 | – | – | – | – | – | ℤ2 × ℤ683 | 19.6 |  |
| 5692 | 1423 | – | – | – | – | – | – | – |  |
| 5708 | 1427 | – | – | – | – | – | – | – |  |
| 5748 | 1437 | – | – | – | – | – | – | – |  |
| 5756 | 1439 | – | – | – | – | – | – | – |  |
| 5788 | 1447 | – | – | – | – | – | ℤ2 × ℤ3 × ℤ241 | 58.3 |  |
| 5804 | 1451 | – | – | – | – | – | ℤ2 × ℤ5xℤ5 × ℤ29 | 55.2 |  |
| 5884 | 1471 | – | – | – | – | – | ℤ2 × ℤ3 × ℤ5 × ℤ7xℤ7 | 61.1 |  |
| 5892 | 1473 | – | – | – | – | – | ℤ8xℤ8 × ℤ23 | 35.4 |  |
| 5948 | 1487 | – | – | – | – | – | – | – |  |
| 5996 | 1499 | – | – | – | – | – | – | – |  |
| 6036 | 1509 | – | – | – | – | – | ℤ2xℤ2 × ℤ13 × ℤ29 | 50.8 |  |
| 6172 | 1543 | – | – | – | – | – | – | – |  |
| 6236 | 1559 | – | – | – | – | – | – | – |  |
| 6268 | 1567 | – | – | – | – | – | – | – |  |
| 6284 | 1571 | – | – | – | – | – | – | – |  |
| 6316 | 1579 | – | – | – | – | – | – | – |  |
| 6332 | 1583 | – | – | – | – | – | ℤ2 × ℤ7 × ℤ113 | 64.9 |  |
| 6356 | 1589 | – | – | – | – | – | – | – |  |
| 6476 | 1619 | – | – | – | – | – | – | – |  |
| 6652 | 1663 | – | – | – | – | – | ℤ2 × ℤ3 × ℤ277 | 49.7 |  |
| 6668 | 1667 | 17 | 99 | 49 | 32.0 | 37.5 | ℤ2 × ℤ7xℤ7 × ℤ17 | 59.0 |  |
| 6756 | 1689 | – | – | – | – | – | ℤ8 × ℤ211 | 57.6 |  |
| 6796 | 1699 | – | – | – | – | – | – | – | Yamada from a symmetric conference matrix of order 850 (cannot exist) |
| 6852 | 1713 | 38 | 62 | 31 | 5.4 | 32.1 | – | – |  |
| 6876 | 1719 | – | – | – | – | – | ℤ2 × ℤ859 | 22.1 |  |
| 6988 | 1747 | – | – | – | – | – | – | – |  |
| 7132 | 1783 | – | – | – | – | – | ℤ2 × ℤ3xℤ3xℤ3xℤ3 × ℤ11 | 54.1 |  |
| 7148 | 1787 | 19 | 95 | 47 | 30.0 | 36.1 | – | – |  |
| 7180 | 1795 | – | – | – | – | – | – | – |  |
| 7292 | 1823 | – | – | – | – | – | ℤ2 × ℤ911 | 26.5 |  |
| 7324 | 1831 | – | – | – | – | – | – | – |  |
| 7364 | 1841 | – | – | – | – | – | – | – |  |
| 7388 | 1847 | – | – | – | – | – | ℤ2 × ℤ13 × ℤ71 | 72.5 |  |
| 7484 | 1871 | – | – | – | – | – | – | – |  |
| 7516 | 1879 | – | – | – | – | – | – | – | Kiyasu from a skew Hadamard matrix of order 940 |
| 7532 | 1883 | – | – | – | – | – | ℤ2 × ℤ941 | 20.5 |  |
| 7572 | 1893 | 42 | 62 | 31 | 4.4 | 33.1 | ℤ2xℤ2 × ℤ11 × ℤ43 | 70.7 |  |
| 7628 | 1907 | – | – | – | – | – | – | – |  |
| 7660 | 1915 | – | – | – | – | – | – | – | Yamada from a symmetric conference matrix of order 958 (cannot exist) |
| 7716 | 1929 | – | – | – | – | – | – | – |  |
| 7876 | 1969 | – | – | – | – | – | – | – |  |
| 7916 | 1979 | – | – | – | – | – | – | – |  |
| 7924 | 1981 | – | – | – | – | – | – | – |  |
| 7948 | 1987 | – | – | – | – | – | – | – |  |
| 8108 | 2027 | – | – | – | – | – | ℤ2 × ℤ1013 | 44.7 |  |
| 8156 | 2039 | – | – | – | – | – | – | – |  |
| 8252 | 2063 | – | – | – | – | – | – | – |  |
| 8332 | 2083 | – | – | – | – | – | – | – | Yamada from a symmetric conference matrix of order 1042 (cannot exist) |
| 8348 | 2087 | – | – | – | – | – | – | – |  |
| 8380 | 2095 | – | – | – | – | – | – | – |  |
| 8396 | 2099 | – | – | – | – | – | – | – |  |
| 8620 | 2155 | – | – | – | – | – | – | – | Yamada from a symmetric conference matrix of order 1078 (cannot exist) |
| 8684 | 2171 | – | – | – | – | – | – | – |  |
| 8812 | 2203 | – | – | – | – | – | – | – |  |
| 8828 | 2207 | – | – | – | – | – | ℤ2 × ℤ1103 | 40.2 |  |
| 8860 | 2215 | 34 | 81 | 41 | 1.9 | 50.1 | – | – | Kiyasu from a skew Hadamard matrix of order 1108 |
| 9004 | 2251 | – | – | – | – | – | – | – |  |
| 9148 | 2287 | – | – | – | – | – | – | – |  |
| 9172 | 2293 | – | – | – | – | – | – | – | Miyamoto from H(2292) |
| 9260 | 2315 | – | – | – | – | – | – | – |  |
| 9308 | 2327 | – | – | – | – | – | – | – |  |
| 9356 | 2339 | – | – | – | – | – | – | – |  |
| 9388 | 2347 | – | – | – | – | – | – | – |  |
| 9484 | 2371 | – | – | – | – | – | – | – |  |
| 9532 | 2383 | – | – | – | – | – | – | – | Kiyasu from a skew Hadamard matrix of order 1192 |
| 9596 | 2399 | – | – | – | – | – | – | – |  |
| 9692 | 2423 | – | – | – | – | – | – | – |  |
| 9716 | 2429 | – | – | – | – | – | – | – |  |
| 9836 | 2459 | – | – | – | – | – | – | – |  |
| 10012 | 2503 | – | – | – | – | – | – | – |  |
| 10052 | 2513 | – | – | – | – | – | – | – |  |
| 10060 | 2515 | – | – | – | – | – | – | – |  |
| 10124 | 2531 | – | – | – | – | – | – | – |  |
| 10172 | 2543 | – | – | – | – | – | ℤ2 × ℤ31 × ℤ41 | 46.5 |  |
| 10316 | 2579 | – | – | – | – | – | – | – |  |
| 10356 | 2589 | – | – | – | – | – | – | – |  |
| 10364 | 2591 | – | – | – | – | – | – | – |  |
| 10516 | 2629 | – | – | – | – | – | – | – |  |
| 10588 | 2647 | – | – | – | – | – | – | – |  |
| 10636 | 2659 | – | – | – | – | – | – | – | Yamada from a symmetric conference matrix of order 1330 (cannot exist) |
| 10644 | 2661 | – | – | – | – | – | – | – |  |
| 10684 | 2671 | – | – | – | – | – | – | – |  |
| 10708 | 2677 | – | – | – | – | – | – | – | Miyamoto from H(2676) |
| 10732 | 2683 | – | – | – | – | – | – | – |  |
| 10748 | 2687 | – | – | – | – | – | – | – |  |
| 10796 | 2699 | – | – | – | – | – | – | – |  |
| 10844 | 2711 | – | – | – | – | – | – | – |  |
| 10924 | 2731 | – | – | – | – | – | – | – | Yamada from a symmetric conference matrix of order 1366 (cannot exist) |
| 10932 | 2733 | – | – | – | – | – | – | – |  |
| 11068 | 2767 | – | – | – | – | – | – | – |  |
| 11212 | 2803 | – | – | – | – | – | – | – | Yamada from a symmetric conference matrix of order 1402 (cannot exist) |
| 11260 | 2815 | – | – | – | – | – | – | – |  |
| 11276 | 2819 | – | – | – | – | – | – | – |  |
| 11292 | 2823 | – | – | – | – | – | – | – |  |
| 11364 | 2841 | – | – | – | – | – | – | – |  |
| 11372 | 2843 | – | – | – | – | – | – | – |  |
| 11420 | 2855 | – | – | – | – | – | ℤ2 × ℤ1427 | 52.1 |  |
| 11460 | 2865 | – | – | – | – | – | – | – |  |
| 11516 | 2879 | – | – | – | – | – | – | – |  |
| 11548 | 2887 | – | – | – | – | – | – | – |  |
| 11572 | 2893 | – | – | – | – | – | – | – |  |
| 11596 | 2899 | – | – | – | – | – | – | – | Yamada from a symmetric conference matrix of order 1450 (cannot exist) |
| 11612 | 2903 | – | – | – | – | – | ℤ2 × ℤ1451 | 56.7 |  |
| 11652 | 2913 | – | – | – | – | – | – | – |  |
| 11708 | 2927 | – | – | – | – | – | – | – |  |
| 11756 | 2939 | – | – | – | – | – | – | – |  |
| 11804 | 2951 | – | – | – | – | – | – | – |  |
| 11852 | 2963 | – | – | – | – | – | ℤ2 × ℤ1481 | 47.1 |  |
| 11884 | 2971 | 55 | 55 | 27 | 14.0 | 19.2 | – | – | Yamada from a symmetric conference matrix of order 1486 (cannot exist) |
| 11892 | 2973 | – | – | – | – | – | – | – |  |
| 11980 | 2995 | – | – | – | – | – | – | – |  |
| 11996 | 2999 | – | – | – | – | – | – | – |  |
