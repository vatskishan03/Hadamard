# Unresolved orders

Table 4 of Cati and Pasechnik lists 195 odd n ≤ 2999 for which no Hadamard matrix of order 4n was known in 2025. 70
of them are settled now, and the 125 that remain are listed below with the cheapest search structures known for each.
The smallest is n = 599, order 2396. Of the 125, 87 have 4n < 10000. Both counts treat n = 573 and 647 (orders 2292
and 2588, Middletownbooks/Hadamard) as settled; if they are not, 127 remain and 89 have 4n < 10000.

Settled:

- 12 by Alpöge, Voinov and Reynolds-Haertle (2026);
- 19 by other public constructions of 2026 (Schneider; the JD-Jones-ASES repositories Hadamard-T and Hadamard-M;
  Schumer). Five of them, n = 131 · t with t = 7, 11, 13, 17, 19, were already reachable in 2025: Cooper–Wallis with
  the T-matrices of order 131 that London and Kotsireas's Turyn-type sequences TT(44) give ([NOTES.md](NOTES.md) §28);
- 6 by Đoković in 2010: n = 823, 883, 1063, 1303, 2143 and 2335 (arXiv:1008.2043, Propositions 3.1 and 3.4);
- 1 in 1988: n = 1823, from Seberry and Whiteman's symmetric conference matrix of order 3646
  ([NOTES.md](NOTES.md) §7);
- 2 by the GitHub repository Middletownbooks/Hadamard: n = 573 (order 2292, 9 October 2026; our independent matrix is
  [hadamard/2292](hadamard/2292)) and n = 647 (order 2588, 14 October 2026);
- 30 in this repository:
  - n = 523, 537, 571, 659, 859, 933, 991, 1051, 1123, 1255, 1483, 1699 and 2293 by searches and recursive theorems
    (orders [2092](hadamard/2092), [2148](hadamard/2148), [2284](hadamard/2284), [2636](hadamard/2636),
    [3436](hadamard/3436), [3732](hadamard/3732), [3964](hadamard/3964), [4204](hadamard/4204),
    [4492](hadamard/4492), [5020](hadamard/5020), [5932](hadamard/5932), [6796](hadamard/6796) and
    [9172](hadamard/9172)); 571 concurrently with Zaitsev et al.;
  - n = 789, 1315, 1509, 1589, 1689, 1841, 1915, 2155, 2515, 2661, 2733, 2815, 2893 and 2951 by Paley's matrix in
    sixteen circulant blocks times T-matrices, two published theorems that the tables did not combine
    ([hadamard/PALEY_BLOCKS.md](hadamard/PALEY_BLOCKS.md));
  - n = 1527, 1883 and 2545 by Miyamoto's matrix in sixteen bordered blocks with near-Williamson matrices of 2026,
    times T-matrices ([hadamard/MIYAMOTO_BLOCKS.md](hadamard/MIYAMOTO_BLOCKS.md#near-williamson-matrices)).

**Corollary 8.8(i).** Đoković (arXiv:1008.2043, Propositions 3.2–3.4) also claimed n = 1527, 2545 and 2571, through
Williamson-type matrices of orders 509 and 857 that Seberry–Yamada's Corollary 8.8(i) would give. Its printed proof
does not give Williamson-type matrices ([NOTES.md](NOTES.md) §7, §19), and we found no other proof. 1527 and 2545 are
constructed here without it; 2571 = 3 · 857 is listed below. For the same reason, ten n that Table 4 omits as known
had no valid construction: 519, 865, 1211, 1465, 1557, 1903, 2051, 2249, 2427 and 2611. All ten are constructed here
([hadamard](hadamard/README.md#orders-recorded-only-through-corollary-88i)).

**Columns.** For each order the table lists the cheapest structure of each kind found by `search/triage.py` that is
expected to contain at least ten solutions:

- multiplier subgroups of ℤ_n* for Goethals–Seidel,
- groups of order n − 1 for Spence's bordered construction.

log₁₀ E is the expected number of solutions in the structure. log₁₀ W is the predicted matching work per solution with
lists of 10⁹ candidates (see [METHOD.md](METHOD.md)). The searches so far needed about a tenth of W. One machine
performs about 10¹¹ list operations per hour, so orders up to log₁₀ W ≈ 13 are within reach. A dash means that no
structure of that kind is expected to contain ten solutions; n.c. means not computed.

**Smaller structures.** Structures too small for the Gaussian model to describe reliably were searched exhaustively:

- 1,319 structures in total: Spence structures with pools of at most 4 · 10⁴ sequences, and Goethals–Seidel structures
  with at most 18 cosets.
- For primes above about 400 there are no Goethals–Seidel structures of the latter kind, by the √p bound in METHOD.md.
- 1,159 of these were searched completely and contain no solution. 80 were searched partially, and 80 have more than
  22 orbits and were not searched.

**Route column.** The last column names a recursive theorem that would settle the order from a smaller object, where
one applies:

- Miyamoto: an Hadamard matrix of order n − 1.
- Kiyasu: a skew-Hadamard matrix of order (n + 1)/2.
- Yamada: a symmetric conference matrix of order (n + 1)/2, which cannot exist when (n − 1)/2 is not a sum of two
  squares.
- complex Yamada–Kiyasu: a quaternary conference matrix of order (n + 1)/2 (zero diagonal, entries ±1, ±i elsewhere),
  when n − 2 is a prime power ≡ 1 (mod 4) (the builder takes primes). It has no sum-of-two-squares condition; it
  settled 859, 1051 and 1699 ([hadamard/3436](hadamard/3436), [hadamard/4204](hadamard/4204),
  [hadamard/6796](hadamard/6796), [NOTES.md](NOTES.md) §31, §34).
- Cooper–Wallis: Williamson-type matrices of order w, for n = t·w with T-matrices of order t known. Listed for w = 191
  (t = 5, 9, 15), where neither T-matrices of order 191 nor near-Williamson matrices of order 191 with a multiplier
  symmetry exist ([search/NOTES.md](search/NOTES.md) §6, §37). Also listed for w = 857, whose Williamson-type matrices
  are claimed only through Seberry–Yamada's Corollary 8.8(i).
- Miyamoto blocks ([hadamard/MIYAMOTO_BLOCKS.md](hadamard/MIYAMOTO_BLOCKS.md)): for n = q·t with q = 4w + 1 a prime
  and T-matrices of order t known, symmetric Williamson-type matrices of order w (near-Williamson matrices give them;
  Kharaghani et al. define those for odd w only) or a negaperiodic Golay pair of length 2w. Not listed in the table;
  the open cases are 2571 (w = 214), 2823 (w = 235) and 2677 (q = 2677, t = 1, which needs a Hadamard matrix of order
  2676 = 4 · 669). 2823 would follow from a symmetric conference matrix of order 118 (Brouwer's table lists
  SRG(117, 58, 28, 29) as existing; the matrix is not public; see
  [hadamard/MIYAMOTO_BLOCKS.md](hadamard/MIYAMOTO_BLOCKS.md)).

| 4n | n | GS: \|H\| | cosets | classes | log₁₀ E | log₁₀ W | Spence: group | log₁₀ W | Route |
|---:|---:|---:|---:|---:|---:|---:|---|---:|---|
| 2396 | 599 | – | – | – | – | – | ℤ₂ × ℤ₁₃ × ℤ₂₃ | 36.4 |  |
| 2572 | 643 | – | – | – | – | – | – | – | Yamada from a symmetric conference matrix of order 322 (cannot exist); complex Yamada–Kiyasu from a quaternary one |
| 2676 | 669 | – | – | – | – | – | – | – |  |
| 2876 | 719 | – | – | – | – | – | – | – |  |
| 2956 | 739 | 9 | 83 | 41 | 33.5 | 22.7 | ℤ₂ × ℤ₃ × ℤ₃ × ℤ₄₁ | 42.4 |  |
| 3004 | 751 | – | – | – | – | – | ℤ₂ × ℤ₃ × ℤ₅ × ℤ₅ × ℤ₅ | 46.0 |  |
| 3356 | 839 | – | – | – | – | – | ℤ₂ × ℤ₄₁₉ | 15.0 |  |
| 3452 | 863 | – | – | – | – | – | – | – |  |
| 3628 | 907 | – | – | – | – | – | ℤ₂ × ℤ₃ × ℤ₁₅₁ | 33.3 |  |
| 3788 | 947 | – | – | – | – | – | ℤ₂ × ℤ₁₁ × ℤ₄₃ | 37.4 |  |
| 3820 | 955 | 19 | 55 | 27 | 2.9 | 28.8 | ℤ₂ × ℤ₃ × ℤ₃ × ℤ₅₃ | 52.8 | Yamada from a symmetric conference matrix of order 478 (or a quaternary one, complex Yamada–Kiyasu); Cooper–Wallis from Williamson-type matrices of order 191 |
| 3884 | 971 | – | – | – | – | – | – | – |  |
| 4076 | 1019 | – | – | – | – | – | – | – |  |
| 4124 | 1031 | – | – | – | – | – | ℤ₂ × ℤ₅ × ℤ₁₀₃ | 40.8 |  |
| 4156 | 1039 | – | – | – | – | – | – | – |  |
| 4348 | 1087 | – | – | – | – | – | ℤ₂ × ℤ₃ × ℤ₁₈₁ | 41.9 |  |
| 4412 | 1103 | 19 | 59 | 29 | 21.1 | 15.5 | ℤ₂ × ℤ₁₉ × ℤ₂₉ | 40.1 |  |
| 4460 | 1115 | – | – | – | – | – | – | – |  |
| 4676 | 1169 | – | – | – | – | – | ℤ₁₆ × ℤ₇₃ | 50.2 |  |
| 4748 | 1187 | – | – | – | – | – | ℤ₂ × ℤ₅₉₃ | 15.2 |  |
| 4892 | 1223 | 13 | 95 | 47 | 33.7 | 32.4 | ℤ₂ × ℤ₁₃ × ℤ₄₇ | 48.8 |  |
| 5028 | 1257 | 22 | 78 | 39 | 12.7 | 37.8 | ℤ₈ × ℤ₁₅₇ | 58.9 |  |
| 5036 | 1259 | – | – | – | – | – | – | – |  |
| 5132 | 1283 | – | – | – | – | – | – | – |  |
| 5164 | 1291 | 15 | 87 | 43 | 30.2 | 29.3 | ℤ₂ × ℤ₃ × ℤ₅ × ℤ₄₃ | 57.1 | Yamada from a symmetric conference matrix of order 646 (cannot exist); complex Yamada–Kiyasu from a quaternary one |
| 5276 | 1319 | – | – | – | – | – | – | – |  |
| 5308 | 1327 | – | – | – | – | – | ℤ₂ × ℤ₃ × ℤ₁₃ × ℤ₁₇ | 58.6 |  |
| 5468 | 1367 | – | – | – | – | – | ℤ₂ × ℤ₆₈₃ | 19.6 |  |
| 5692 | 1423 | – | – | – | – | – | – | – |  |
| 5708 | 1427 | – | – | – | – | – | – | – |  |
| 5748 | 1437 | – | – | – | – | – | – | – |  |
| 5756 | 1439 | – | – | – | – | – | – | – |  |
| 5788 | 1447 | – | – | – | – | – | ℤ₂ × ℤ₃ × ℤ₂₄₁ | 58.3 |  |
| 5804 | 1451 | – | – | – | – | – | ℤ₂ × ℤ₅ × ℤ₅ × ℤ₂₉ | 55.2 |  |
| 5884 | 1471 | – | – | – | – | – | ℤ₂ × ℤ₃ × ℤ₅ × ℤ₇ × ℤ₇ | 61.1 |  |
| 5892 | 1473 | – | – | – | – | – | ℤ₈ × ℤ₈ × ℤ₂₃ | 35.4 |  |
| 5948 | 1487 | – | – | – | – | – | – | – |  |
| 5996 | 1499 | – | – | – | – | – | – | – |  |
| 6172 | 1543 | – | – | – | – | – | – | – |  |
| 6236 | 1559 | – | – | – | – | – | – | – |  |
| 6268 | 1567 | – | – | – | – | – | – | – |  |
| 6284 | 1571 | – | – | – | – | – | – | – |  |
| 6316 | 1579 | – | – | – | – | – | – | – |  |
| 6332 | 1583 | – | – | – | – | – | ℤ₂ × ℤ₇ × ℤ₁₁₃ | 64.9 |  |
| 6476 | 1619 | – | – | – | – | – | – | – |  |
| 6652 | 1663 | – | – | – | – | – | ℤ₂ × ℤ₃ × ℤ₂₇₇ | 49.7 |  |
| 6668 | 1667 | 17 | 99 | 49 | 32.0 | 37.5 | ℤ₂ × ℤ₇ × ℤ₇ × ℤ₁₇ | 59.0 |  |
| 6852 | 1713 | 38 | 62 | 31 | 5.4 | 32.1 | – | – |  |
| 6876 | 1719 | – | – | – | – | – | ℤ₂ × ℤ₈₅₉ | 22.1 | Cooper–Wallis from Williamson-type matrices of order 191 |
| 6988 | 1747 | – | – | – | – | – | – | – |  |
| 7132 | 1783 | – | – | – | – | – | ℤ₂ × ℤ₃ × ℤ₃ × ℤ₃ × ℤ₃ × ℤ₁₁ | 54.1 |  |
| 7148 | 1787 | 19 | 95 | 47 | 30.0 | 36.1 | – | – |  |
| 7180 | 1795 | – | – | – | – | – | – | – |  |
| 7324 | 1831 | – | – | – | – | – | – | – |  |
| 7388 | 1847 | – | – | – | – | – | ℤ₂ × ℤ₁₃ × ℤ₇₁ | 72.5 |  |
| 7484 | 1871 | – | – | – | – | – | – | – |  |
| 7516 | 1879 | – | – | – | – | – | – | – | Kiyasu from a skew Hadamard matrix of order 940, or complex Yamada–Kiyasu from a quaternary conference matrix of that order |
| 7572 | 1893 | 42 | 62 | 31 | 4.4 | 33.1 | ℤ₂ × ℤ₂ × ℤ₁₁ × ℤ₄₃ | 70.7 |  |
| 7628 | 1907 | – | – | – | – | – | – | – |  |
| 7716 | 1929 | – | – | – | – | – | – | – |  |
| 7876 | 1969 | – | – | – | – | – | – | – |  |
| 7916 | 1979 | – | – | – | – | – | – | – |  |
| 7924 | 1981 | – | – | – | – | – | – | – |  |
| 7948 | 1987 | – | – | – | – | – | – | – |  |
| 8108 | 2027 | – | – | – | – | – | ℤ₂ × ℤ₁₀₁₃ | 44.7 |  |
| 8156 | 2039 | – | – | – | – | – | – | – |  |
| 8252 | 2063 | – | – | – | – | – | – | – |  |
| 8332 | 2083 | – | – | – | – | – | – | – | Yamada from a symmetric conference matrix of order 1042 (cannot exist); complex Yamada–Kiyasu from a quaternary one |
| 8348 | 2087 | – | – | – | – | – | – | – |  |
| 8380 | 2095 | – | – | – | – | – | – | – |  |
| 8396 | 2099 | – | – | – | – | – | – | – |  |
| 8684 | 2171 | – | – | – | – | – | – | – |  |
| 8812 | 2203 | – | – | – | – | – | – | – |  |
| 8828 | 2207 | – | – | – | – | – | ℤ₂ × ℤ₁₁₀₃ | 40.2 |  |
| 8860 | 2215 | 34 | 81 | 41 | 1.9 | 50.1 | – | – | Kiyasu from a skew Hadamard matrix of order 1108, or complex Yamada–Kiyasu from a quaternary conference matrix of that order |
| 9004 | 2251 | – | – | – | – | – | – | – |  |
| 9148 | 2287 | – | – | – | – | – | – | – |  |
| 9260 | 2315 | – | – | – | – | – | – | – |  |
| 9308 | 2327 | – | – | – | – | – | – | – |  |
| 9356 | 2339 | – | – | – | – | – | – | – |  |
| 9388 | 2347 | – | – | – | – | – | – | – |  |
| 9484 | 2371 | – | – | – | – | – | – | – |  |
| 9532 | 2383 | – | – | – | – | – | – | – | Kiyasu from a skew Hadamard matrix of order 1192, or complex Yamada–Kiyasu from a quaternary conference matrix of that order |
| 9596 | 2399 | – | – | – | – | – | – | – |  |
| 9692 | 2423 | – | – | – | – | – | – | – |  |
| 9716 | 2429 | – | – | – | – | – | – | – |  |
| 9836 | 2459 | – | – | – | – | – | – | – |  |
| 10012 | 2503 | – | – | – | – | – | – | – |  |
| 10052 | 2513 | – | – | – | – | – | – | – |  |
| 10124 | 2531 | – | – | – | – | – | – | – |  |
| 10172 | 2543 | – | – | – | – | – | ℤ₂ × ℤ₃₁ × ℤ₄₁ | 46.5 |  |
| 10284 | 2571 | – | – | – | – | – | n.c. | n.c. | Cooper–Wallis from Williamson-type matrices of order 857 (claimed only through Seberry–Yamada Cor. 8.8(i)) |
| 10316 | 2579 | – | – | – | – | – | – | – |  |
| 10356 | 2589 | – | – | – | – | – | – | – |  |
| 10364 | 2591 | – | – | – | – | – | – | – |  |
| 10516 | 2629 | – | – | – | – | – | – | – |  |
| 10588 | 2647 | – | – | – | – | – | – | – |  |
| 10636 | 2659 | – | – | – | – | – | – | – | Yamada from a symmetric conference matrix of order 1330 (cannot exist); complex Yamada–Kiyasu from a quaternary one |
| 10684 | 2671 | – | – | – | – | – | – | – |  |
| 10708 | 2677 | – | – | – | – | – | – | – | Miyamoto from H(2676) |
| 10732 | 2683 | – | – | – | – | – | – | – |  |
| 10748 | 2687 | – | – | – | – | – | – | – |  |
| 10796 | 2699 | – | – | – | – | – | – | – |  |
| 10844 | 2711 | – | – | – | – | – | – | – |  |
| 10924 | 2731 | – | – | – | – | – | – | – | Yamada from a symmetric conference matrix of order 1366 (cannot exist); complex Yamada–Kiyasu from a quaternary one |
| 11068 | 2767 | – | – | – | – | – | – | – |  |
| 11212 | 2803 | – | – | – | – | – | – | – | Yamada from a symmetric conference matrix of order 1402 (cannot exist); complex Yamada–Kiyasu from a quaternary one |
| 11276 | 2819 | – | – | – | – | – | – | – |  |
| 11292 | 2823 | – | – | – | – | – | – | – |  |
| 11364 | 2841 | – | – | – | – | – | – | – |  |
| 11372 | 2843 | – | – | – | – | – | – | – |  |
| 11420 | 2855 | – | – | – | – | – | ℤ₂ × ℤ₁₄₂₇ | 52.1 |  |
| 11460 | 2865 | – | – | – | – | – | – | – | Cooper–Wallis from Williamson-type matrices of order 191 |
| 11516 | 2879 | – | – | – | – | – | – | – |  |
| 11548 | 2887 | – | – | – | – | – | – | – |  |
| 11596 | 2899 | – | – | – | – | – | – | – | Yamada from a symmetric conference matrix of order 1450 (cannot exist); complex Yamada–Kiyasu from a quaternary one |
| 11612 | 2903 | – | – | – | – | – | ℤ₂ × ℤ₁₄₅₁ | 56.7 |  |
| 11652 | 2913 | – | – | – | – | – | – | – |  |
| 11708 | 2927 | – | – | – | – | – | – | – |  |
| 11756 | 2939 | – | – | – | – | – | – | – |  |
| 11852 | 2963 | – | – | – | – | – | ℤ₂ × ℤ₁₄₈₁ | 47.1 |  |
| 11884 | 2971 | 55 | 55 | 27 | 14.0 | 19.2 | – | – | Yamada from a symmetric conference matrix of order 1486 (cannot exist); complex Yamada–Kiyasu from a quaternary one |
| 11892 | 2973 | – | – | – | – | – | – | – |  |
| 11980 | 2995 | – | – | – | – | – | – | – |  |
| 11996 | 2999 | – | – | – | – | – | – | – |  |
