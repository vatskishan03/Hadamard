# Novelty

How we decided whether an order is new, and what we searched. Every claim of novelty in this repository means "new
to our knowledge" in the sense below. A literature search cannot show that no construction exists elsewhere;
corrections are welcome.

## Status labels

| Label | Meaning |
|---|---|
| **New** | The order is open in the reference table for its family (below), and no later construction was found. |
| **New, but obtainable classically** | As New, but every ingredient (theorems and objects) was published before 2026; only their combination is new. |
| **New, beyond the tables** | The order lies above the range of the reference tables. Where Table 1.53 of the Handbook covers it, the table gives t ≥ 3 (no Hadamard matrix of order 4n known), and no later construction was found. |
| **New, beyond the tables, but obtainable classically** | As New, beyond the tables, but every ingredient (theorems and objects) was published before 2026; only their combination is new. |
| **First explicit construction** | A table records the order as known, but we could not trace any construction behind the entry. |
| **Concurrent** | No published construction was found, but another group constructed the order independently at about the same time. We claim no priority. |
| **Known** | An earlier construction exists. The matrix here is independent of it. |

For skew and symmetric matrices the property is what is new: an ordinary Hadamard matrix of the same order may have
been known. Each order's README says so where it applies.

When a table's entry rests on a printed route that fails (Corollary 8.8(i) of Seberry–Yamada, [NOTES.md](NOTES.md) §7),
the order is labelled First explicit construction even if a repaired route from ingredients published before 2026
exists; the README names that route where we found one.

## Reference tables

| Family | Table | Range |
|---|---|---|
| Hadamard | M. Cati, D. V. Pasechnik, *A database of constructions of Hadamard matrices*, arXiv:2411.18897 (v1 of 2024-11-28 and v2 of 2025-08-30 have the same rows), Table 4 | odd n ≤ 2999 |
| Hadamard | *Handbook of Combinatorial Designs*, 2nd ed. (2007), Table 1.53 | odd n < 10000 |
| Skew | D. Ž. Đoković, J. Seberry, census of skew and symmetric Hadamard matrices of order 4v (2025-05-28), and Cati–Pasechnik Table 5 | odd v < 2500 |
| Symmetric | the same census: 595 undecided v (three of them settled by Đoković, arXiv:2508.17141); Handbook Table V.1.52 | odd v < 2500 |
| Complex | J. Seberry, M. Yamada, *Hadamard matrices, sequences, and block designs*, in *Contemporary Design Theory* (Wiley, 1992), Table 11.2 | odd N < 1000 |
| Williamson-type | Seberry–Yamada (1992), Table A.1 | |
| T-matrices | Seberry–Yamada (1992), §3 | t ≤ 210 |

We also read the older tables these replace: Seberry–Yamada (1992), Tables 7.1 and A.2, and, for most skew orders,
Koukouvinos–Stylianou, Discrete Math. 308 (2008), Tables 1 and 3. A README says so where one of these was not read.
The Handbook tables were read from a copy available online and not checked against the printed edition.

## Sources searched

The searches were made between 12 and 17 October 2026.

**For every order.**

- Đoković's updates of the Hadamard tables (arXiv:0912.5091, arXiv:1008.2043), Đoković–Golubitsky–Kotsireas
  (arXiv:1301.3671), and the papers of Đoković and coauthors on skew and symmetric matrices, among them
  arXiv:1511.08734, arXiv:2404.14375, arXiv:2508.17141, Đoković–Kotsireas, Math. Comput. Sci. 12 (2018) 373–388, and
  Abuzin–Balonin–Đoković–Kotsireas, Inform. Control Syst. 2019, no. 5, 2–9 (Table 1, skew-propus families).
- SageMath, development version of October 2026: `hadamard_matrix(n, existence=True)`,
  `skew_hadamard_matrix(n, existence=True)` and the list `unknown_skew_hadamard`.
- The constructions published in 2026: Alpöge, Voinov and Reynolds-Haertle (2026-08-12); Schneider (2026-08-23); the
  JD-Jones-ASES repositories Hadamard-T, Hadamard-M, Hadamard-B and Hadamard-formal; Schumer, orders 3676, 5436, 6644
  and 7852 (Zenodo, 2026-09-17, doi:10.5281/zenodo.22820997); the deposits of Piza, Kline, Bedani, Karoui and
  Blackthorn; Szöllősi, arXiv:2603.09572; Kharaghani, Mohammadian, Tayfeh-Rezaie and Zaitsev, arXiv:2605.08661
  (near-Williamson matrices); H. Verma, A symmetric conference matrix of order 86, arXiv:2610.08893 (2026-10-06),
  which gives Williamson-type matrices of order 173 through Seberry–Yamada's Corollary 8.8(ii) ([NOTES.md](NOTES.md)
  §7).
- Zenodo (records matching "hadamard"), GitHub, the e-print sources of Hadamard-related arXiv papers (to 2026-10-08)
  and the math.CO listings (to 2026-10-12), the OEIS, MathOverflow, Balonin's online catalogue and the Palomar
  registry.
- Found after these searches, on 19 October 2026: the GitHub repository Middletownbooks/Hadamard (see Corrections).

**Complex orders.** Ito's papers on Hadamard groups (1994–2000); B. Schmidt, Des. Codes Cryptogr. 17 (1999) 61–68;
R. Egan, Des. Codes Cryptogr. 87 (2019) 67–74, Theorem 4.5, which gives every order 2N with 4N − 1 a prime power;
V. Zaitsev's page Quaternary Hadamard matrices (https://www.cs.uleth.ca/~zaitsevv/quaternary.html, with Kharaghani,
Mohammadian and Tayfeh-Rezaie, last updated 2026-10-14), with matrices of 122 orders 2N; Czifra, Matolcsi and
Szöllősi, arXiv:2511.11230 (Butson matrices BH(n, q) for q = 6, 10 and 14, none quaternary); the quaternary Legendre
pairs of arXiv:2609.04589 and arXiv:2609.25573 (lengths 42 to 80); the skew-regular quaternary Hadamard matrices of
order 1 + p² of arXiv:2510.06489; Đoković's Williamson-type matrices of orders 509, 857 and 883 (arXiv:1008.2043).

**T-matrices.** S. London and I. Kotsireas, New Turyn-type sequences, Cryptogr. Commun. 17 (2025) 1601–1610 (TT(40),
TT(42), TT(44)), which give T-sequences of length 131.

**Paley's matrix in circulant blocks** ([hadamard/PALEY_BLOCKS.md](hadamard/PALEY_BLOCKS.md)). Seberry–Yamada, JCMCC 7
(1990) 97–137 (Theorem 4) and its Kobe proceedings version; the 1992 survey (Theorems 3.10 and 8.1); Seberry–Yamada,
JCMCC 13 (1993) 97–106 (Theorems 1 and 6, Corollary 2); Schmidt (1999); Schmidt–Tan, Des. Codes Cryptogr. 73 (2014)
105–119; Balonin–Đoković, arXiv:1508.00640; Barrera Acevedo–Ó Catháin–Dietrich, arXiv:1904.11460; Xia–Xia–Seberry–Wu,
Discrete Appl. Math. 145 (2005) 498–504; the literature on cocyclic and dicyclic Hadamard matrices. Both steps of the
construction are published; we found no source that combines them. Yamada's paper of 1991 and Ito's papers of 1994
were read through Seberry–Yamada (1993) and Schmidt (1999). For the orders beyond Table 4 also, for each order,
SageMath's `hadamard_matrix(4n, existence=True)` (Unknown except 16892 and 25132), the e-print texts for every n and
4n, and Cooper–Wallis, Miyamoto, Yamada–Kiyasu and Paley with the T-matrix and Williamson-type orders known in 2026.

**Miyamoto's matrix in bordered blocks** ([hadamard/MIYAMOTO_BLOCKS.md](hadamard/MIYAMOTO_BLOCKS.md)) and the orders
recorded only through Corollary 8.8(i). Seberry–Yamada, JCMCC 7 (1990), Theorems 4, 21, 34 and Corollaries 22, 35,
its Kobe version (Theorems 1, 8, 21, Corollaries 9, 22) and the 1992 survey (Theorems 3.10, 8.1, 8.2, 8.10, Corollary
8.11; Tables A.1, A.2; §5, Yang's theorems): Miyamoto's matrix is an M-structure there but is never multiplied by
T-matrices, and Theorem 21 admits only t = 1 or even t. Seberry–Yamada, JCMCC 13 (1993), Theorems 6 and 7; Handbook
§V.2, 2.45–2.49; Đoković, arXiv:0912.5091 and arXiv:1008.2043; Balonin–Đoković, arXiv:1508.00640 (negaperiodic Golay
pairs, and their list (31) of open quasi-Williamson orders); Đoković–Kotsireas, arXiv:1310.5773 (periodic Golay
pairs); SageMath for every order (Unknown for all); Paley I and II, Miyamoto, Yamada–Kiyasu (with the skew-Hadamard
matrices of the census and of this repository) and the T × WT splits of Cooper–Wallis for each n; the e-print texts for
Miyamoto's theorem combined with T-matrices or block multiplication, and for every n and 4n (only Đoković's claims of
7179 and 7323 in arXiv:1008.2043). Yang's multiplication is used as restated in the 1992 survey (Theorem 5.23).

**Skew-Hadamard matrices of order 3(q + 1)** ([skew/PALEY_SKEW.md](skew/PALEY_SKEW.md)). Seberry, Ars Combin. 6
(1978) (Summary 2, Theorems 2–7, Table 2); Seberry–Yamada, JCMCC 7 (1990) and its Kobe version, Utilitas Math. 40
(1991), the 1992 survey (§7, the families SI–SXI, Table 7.1, Theorem 3.10, §8) and JCMCC 13 (1993); Wallis's lists
of 1971 and 1973; Xia–Xia–Seberry, Discrete Appl. Math. 145 (2005), and the group's papers of 2005–2006; Handbook §V.1
with Table 1.51 (online copy); Koukouvinos–Stylianou (2008), first pages only; Cati–Pasechnik Table 5 and footnote
13; the census; Đoković, arXiv:2508.17141 (XXSW at order q + 1 only) and his earlier skew papers; SageMath
(`unknown_skew_hadamard`, and the branches of `skew_hadamard_matrix`, Unknown for all 58 orders); the arXiv
skew-Hadamard records and the Hadamard records of 2025–2026, Zenodo, the JD-Jones-ASES repositories and Schumer's
deposit. None gives skew matrices of order 3(q + 1) in general or of any of the 58 orders; Seberry–Yamada's Table 7.1
and Cati–Pasechnik's Table 5 record 663 and 879 as known without a construction. JCMCC 13 (1993), Theorem 1 (after
Yamada 1991) contains the sixteen-block form with the transposed pair, so the orders are New, but obtainable
classically. Wallis–Whiteman (1972), as stated by Momihara–Xiang (arXiv:1801.08776), gives 7656.

**Symmetric Hadamard matrices of order 4q** ([symmetric/PROPUS_MIYAMOTO.md](symmetric/PROPUS_MIYAMOTO.md)).
Seberry–Yamada, JCMCC 7 (1990), Lemma 26, Corollary 29 and Remarks 19 and 30, and its Kobe version (Corollary 16,
Remark 17); Seberry–Balonin, arXiv:1512.01732 and Australas. J. Combin. 69 (2017) (Theorem 3, Corollary 1, the
unresolved list); the census (its series, sporadic values and remark on Miyamoto's Theorem 5); Handbook Table V.1.52
as quoted by Đoković; Balonin–Đoković, arXiv:1803.08787; Đoković, arXiv:2101.05429, 2301.02751, 2404.14375 and
2508.17141; Balonin–Đoković–Karbovskiy, arXiv:1708.05098 and 1710.03037; Di Matteo–Đoković–Kotsireas,
arXiv:1503.04226; chapter 1 and the index of Seberry–Yamada (2020); Yamada, J. Austral. Math. Soc. A 46 (1989);
Seberry, Australas. J. Combin. 55 (2013); SageMath; Balonin's online catalogue (propus pages); arXiv, Semantic Scholar
and OpenAlex for Miyamoto 1991 and its citations. None states the family; Miyamoto's claim of 1991 is discussed in
PROPUS_MIYAMOTO.md.

**Symmetric Hadamard matrices without a symmetric input**
([symmetric/PROPUS_WEIGHING.md](symmetric/PROPUS_WEIGHING.md),
[symmetric/YAMADA_KIYASU.md](symmetric/YAMADA_KIYASU.md)). Checked on 17 October 2026: Yamada, J. Austral. Math. Soc.
A 46 (1989), read in full (Theorems 3–5, ordinary matrices only); Seberry–Yamada, JCMCC 7 (1990) and its Kobe version
(Lemma 26, Corollaries 27, 29 and 31, Remarks 17 and 30, Theorem 34 with its Case B, a V-part coupled across sixteen
blocks for ordinary M-structure matrices, and Corollary 35), the 1992 survey (Lemma 8.6, Corollaries 8.7, 8.8, 8.11,
Theorem 8.10, Lemma 11.3) and JCMCC 13 (1993), Corollary 5; Google Books snippets of Seberry–Yamada (2020), §6.5;
Seberry–Balonin (2015, 2017); Balonin–Đoković, arXiv:1803.08787; the census; Đoković, arXiv:1008.2043, 2101.05429,
2404.14375 and 2508.17141; Valtonen, arXiv:2608.06920; Goryainov–Haemers–Konstantinova–Li, arXiv:2512.16653; Barrera
Acevedo–Dietrich, Australas. J. Combin. 73 (2019); Cati–Pasechnik; Balonin's online catalogue; Kiyasu, RIMS Kōkyūroku
587 (1986); arXiv, Crossref, OpenAlex and CiNii searches for symmetric Hadamard matrices, propus matrices and
symmetric weighing matrices; and Zaitsev's page, whose quaternary matrices QH590, QH638, QH782, QH926, QH1310, QH1358,
QH1646 and QH1886 are of the Yamada–Kiyasu kind and monomially equivalent to symmetric ones, so 1180, 1276, 1564,
1852, 2620, 2716, 3292 and 3772 are labelled Concurrent. None states either construction or a symmetric matrix of any
of the 37 orders.

**Symmetric orders 404 and 596** ([symmetric/PERIODIC_GOLAY.md](symmetric/PERIODIC_GOLAY.md)). Checked on 18 October
2026: the census (both v undecided); Miyamoto's claim of 1991 as the census reports it (a symmetric input of order q −
1, rejected by the census); Seberry–Balonin (2015, 2017); Balonin–Đoković, arXiv:1803.08787; Di Matteo, Đoković and
Kotsireas, arXiv:1503.04226; Đoković, arXiv:2101.05429 and 2508.17141; Đoković–Kotsireas, arXiv:1310.5773 (the pair of
length 74); Momihara–Xiang, arXiv:1801.08776; Zaitsev's page (no quaternary matrix of order 202 or 298); web searches
for symmetric Hadamard matrices of orders 404 and 596 and for symmetric constructions from periodic Golay pairs
(Balonin and Đoković's 2015 paper on symmetric two-circulant Hadamard matrices and periodic Golay pairs was seen only
by title; two-circulant matrices from a pair of length ℓ have order 2ℓ, not 4(2ℓ + 1)). None gives a symmetric matrix
of either order or the construction.

**Symmetric orders 612, 1300, 1740, 4900 and 7564** ([symmetric/CONFERENCE_PLUG.md](symmetric/CONFERENCE_PLUG.md)).
Checked on 19 October 2026: the census (all five v undecided); Seberry–Yamada (1992), §11 (Kharaghani–Seberry's
Theorem 11.14 and Corollaries 11.15 and 11.20, complex Hadamard matrices of order w(w − 1), stated without symmetry,
and Lemma 11.3, Turyn's map); Seberry, J. Austral. Math. Soc. A 29 (1980) 235–242 (Mukhopadhyay's recursion with skew
outer designs, ordinary matrices, read from page images); Seberry–Balonin (2017); Balonin–Đoković, arXiv:1803.08787;
Đoković, arXiv:2101.05429 and 2508.17141; web searches for symmetric Hadamard matrices of these orders and of order
2n(n − 1) or 2q(q + 1); read by the round-13 judges: Farouk–Wang, Filomat 34 (2020) 815–834 (Scarpis–Đoković matrices
of order 2q(q + 1), not symmetric) and Đoković, arXiv:2508.17141 in full (the dual scheme, order q(q + 1)). Not read:
Kharaghani–Seberry, Congr. Numer. 75 (1990); Wallis, Linear Multilinear Algebra 1 (1973); Seberry–Wallis, LNM 292, pp.
349–353; Seberry–Zhang (1993); the Seberry–Yamada book of 2020. None of the sources read states that the complex
matrix is complex-symmetric or gives a symmetric matrix of any of the five orders. Every ingredient is older than
2026, hence the label New, but obtainable classically.

**Symmetric order 260** ([symmetric/260](symmetric/260)). Checked on 17 October 2026: the census; Đoković,
arXiv:2101.05429 (65 first among the open v < 200, and the Handbook's Table 1.52 list as he quotes it);
Seberry–Balonin (2015, 2017); Balonin–Đoković, arXiv:1803.08787; Verma, arXiv:2609.25573 (the pair only); Lebedev,
arXiv:2609.04589; Kotsireas–Winterhof and Jedwab–Pender as cited by Verma; Zaitsev's page (QH130, not equivalent to a
symmetric matrix); web searches for a symmetric Hadamard matrix of order 260.

**Order 6796** ([hadamard/6796](hadamard/6796)). Checked on 20 October 2026:
- *Tables and recursions.* Cati–Pasechnik, arXiv:2411.18897 v1 and v2 (Table 4: 1699(3)); Handbook Table 1.53;
  Đoković, arXiv:0912.5091 (the set Δ and four further values) and arXiv:1008.2043 (42 values);
  Đoković–Golubitsky–Kotsireas, arXiv:1301.3671; SageMath development source of 18 October 2026, read branch by
  branch, not run (Unknown).
- *Classical families.* Yamada 1989 (Theorem 5) and Seberry–Yamada 1992 (key y3, §11). For conference matrices of
  order 850: Seberry–Whiteman 1980 (via Egan, arXiv:2309.07522), Handbook §V.6 (BGW and GC families) and Et-Taoui,
  arXiv:1409.5720. Fallon–Iverson, arXiv:2610.09169: a complex conference matrix of order 850, not quaternary in any
  equivalent form.
- *Recent work.* The 2026 sources of "For every order" and Middletownbooks/Hadamard (main and its branch of 18
  October); Zaitsev's page (unchanged since 15 October; orders below 2000); the October 2026 arXiv listings of
  math.CO, cs.DM and cs.IT, with the e-print sources of 20 related papers; Zenodo, GitHub (repositories and code)
  and web searches for 6796, 3398 and 1699.
- *Finding.* None gives a Hadamard matrix of order 6796, a quaternary Hadamard matrix of order 3398 or a quaternary
  conference matrix of order 850, or states the complex Yamada–Kiyasu template.
- *Not read.* The Seberry–Yamada book of 2020 and Seberry–Whiteman 1980 itself.

**Orders 3436 and 4204** ([hadamard/3436](hadamard/3436), [hadamard/4204](hadamard/4204)). Checked on 20 October
2026:
- *Tables and recursions.* Cati–Pasechnik v1 and v2 (Table 4: 859(3), 1051(3)); Handbook Table 1.53 (t = 3 for
  both); Đoković, arXiv:0912.5091 and arXiv:1008.2043; Đoković–Golubitsky–Kotsireas, arXiv:1301.3671; SageMath
  development source of 20 October 2026, read branch by branch, not run (Unknown for both); Seberry–Yamada 1992,
  Tables A.1 and 11.2.
- *Conference matrices of orders 430 and 526.* Belevitch's condition; the generalized Paley, Berman, BGW and GC,
  Seberry–Whiteman (via Egan, arXiv:2309.07522) and Et-Taoui families. Fallon–Iverson's ancillary library
  (arXiv:2610.09169) has complex conference matrices of both orders, not quaternary in any equivalent form.
  Miyamoto's Corollary 4 (through the JD-Jones-ASES Hadamard-M note) needs Cᵀ = C.
- *Recent work.* Zaitsev's page (unchanged since 15 October; 1718 on its unknown list); Middletownbooks/Hadamard,
  both branches (last commit 18 October; 859 listed as work in progress); the JD-Jones-ASES repositories
  (unchanged); Schumer's deposit; Zenodo; the arXiv listings of 20 October and the arXiv API to 19 October; the
  e-print sources of 20 papers of 2026; GitHub code and repository search; Palomar; Balonin's catalogue; web searches
  for 3436, 4204, 1718, 2102, 859 and 1051.
- *Finding.* None gives a Hadamard matrix of order 3436 or 4204, a quaternary Hadamard matrix of order 1718 or 2102,
  or a quaternary conference matrix of order 430 or 526.
- *Not read.* As for 6796, and Miyamoto's §7 lists. ProofWatch could not be reached.

**Order 2572** ([hadamard/2572](hadamard/2572)). Checked on 21 October 2026, before the core was found, by an audit
and an adversarial re-search:
- *Tables and recursions.* Cati–Pasechnik v1 and v2 (Table 4: 643(3); its §6 sentence naming 643 concerns the skew
  table); Đoković, arXiv:0912.5091 and arXiv:1008.2043; Đoković–Golubitsky–Kotsireas, arXiv:1301.3671; SageMath
  development source of 21 October 2026, read branch by branch, not run (Unknown); Seberry–Yamada 1992, Tables A.1 and
  11.2; Brouwer's table (srg(321, 160, 79, 80) does not exist, srg(1285, 642, 320, 321) undecided).
- *Recent work.* The survey of Egan, Kotsireas, London and Ó Catháin, arXiv:2610.27175 (19 October; 4 · 643 among its
  95 open orders, which it lists before our results of 20 October); the math74 explorer v51; Zaitsev's page (last
  modified 15 October; 1286 on its unknown list); Middletownbooks/Hadamard (643 open, screened on 6 October); the
  arXiv listings to 21 October; GitHub and web searches for 2572, 1286 and 643.
- *Conference matrices of order 322.* As for 430 and 526; Fallon–Iverson's library has a complex one, not
  quaternary in any equivalent form (invariants up to 0.765 from {±1, ±i}).
- *Finding.* None gives a Hadamard matrix of order 2572, a quaternary Hadamard matrix of order 1286 or a quaternary
  conference matrix of order 322.
- *Not read.* As for 6796; the Handbook's Table 1.53 row for 643 was not read directly; Miyamoto's §7 lists.

**Complex order 418** ([complex-hadamard](complex-hadamard/README.md)). Checked on 22 October 2026:
- *Tables and recursions.* Seberry–Yamada 1992, Table 11.2 (N = 209 blank). Turyn's route needs 417 = 3 · 139, and
  Ito's dicyclic groups with Egan's Theorem 4.5 need 835 = 5 · 167, to be a prime power. A Kronecker product would
  need a quaternary factor of odd order 11, 19 or 209.
- *Recent work.* Zaitsev's page (last modified 15 October, read on 22 October: 418 on its unknown list, no QH418);
  Kharaghani, Mohammadian, Tayfeh-Rezaie and Zaitsev, arXiv:2605.08661; Szöllősi's talk at IWDMA 2026 (17 March
  2026: complex Hadamard matrices up to order 94, next open case 118); Munim's status report of 14 October
  2026 on dev.to (real orders only); web searches for a quaternary Hadamard matrix of order 418.
- *Finding.* None gives a quaternary Hadamard matrix of order 418.
- *Not read.* As for 6796. Zaitsev et al. search these orders at present; their page lists 418 as unknown.

## Not checked

- J. Seberry, M. Yamada, *Hadamard Matrices: Constructions using Number Theory and Linear Algebra* (Wiley, 2020):
  its updated tables, and chapter 6 (§6.4 Multiplication Theorems Using M-structures, §6.5 Miyamoto's Theorem and
  Corollaries via M-structures), the most likely place for a prior statement of the Miyamoto-block route; only its
  contents, chapter 1 and index were seen. Its index puts "symmetric propus-type" on p. 154 and symmetric Hadamard
  matrices on pp. 211–214, the most likely places for a prior statement of the symmetric 4q family.
- M. Miyamoto, A construction of Hadamard matrices, J. Combin. Theory Ser. A 57 (1991) 86–108, including its
  Corollaries 2–5 on matrices of symmetric block form, to which Seberry and Yamada credit Table A.1's key wl, and the
  lists of its §7, which Cati and Pasechnik could not verify either. We know it through the census's remark, Seberry
  and Yamada's restatements, the zbMATH review (Zbl 0733.05017) and the transcription in the JD-Jones-ASES Hadamard-M
  note, by which list (2) contains 127 (unsupported by Miyamoto's Corollary 4, per that note), so that list (3) may
  record 1527 = 3 · 509.
- Yang's 1989 paper (Proc. Amer. Math. Soc. 107); Whiteman, Discrete Math. 2 (1972); Haviland's thesis (2023) and
  Anderson–Haviland–Holmes–Humphries–Magland, Graphs Combin. 39 (2023), on skew relative difference sets in dicyclic
  groups (by its abstract one infinite family from Ito's and Schmidt's results; our matrices of orders 60 and 132 have
  automorphism groups of order 2, so they are not group-developed of that kind); Xia–Zuo–Lou–Xia, Trans. Comb. 13
  (2024), read only through a summary (its Theorems 2.6 and 2.9 and the sets S1–S3 reach none of our orders).
- J. Seberry, *Orthogonal Designs* (Springer, 2017), and Geramita–Seberry (1979).
- The skew table of the Handbook (Table 1.51) beyond the online copy.
- MathSciNet, zbMATH, and the full texts of paywalled journals from 2025–2026.

A construction recorded only in one of these is the remaining risk for every claim.

## Corrections

Labels of the first release (14 October 2026) that have changed:

- Complex orders 142, 154, 222, 262, 286, 310, 370 and 742: **Known**. For each, 4N − 1 is prime, and Ito's dicyclic
  Hadamard groups with Egan's Theorem 4.5 give the order ([complex-hadamard/PALEY.md](complex-hadamard/PALEY.md)).
- Complex orders 238, 302, 326, 382, 494, 518, 534 and 754: **Concurrent**, with Zaitsev et al.
- 2284: **Concurrent**, through Zaitsev et al.'s quaternary matrix of order 1142 and Turyn's doubling.
- 2292: **Known**. The GitHub repository Middletownbooks/Hadamard published an independent matrix of this order on
  9 October 2026, before ours of 14 October; we learned of it on 19 October.
- 31316: **Known**; SageMath's development version builds it since 3 October 2026.
- 19028, 21748 and 39260: **New, beyond the tables, but obtainable classically**.
- Skew 5700 and 9900, symmetric 5700: **New, but obtainable classically**. The amicable Hadamard matrices of orders 52,
  76 and 100 behind them are first explicit, not new: they come from the skew-propus families of Abuzin, Balonin,
  Đoković and Kotsireas (2019).
- T-matrices of order 131: **Known**, from London and Kotsireas's TT(44) (2025).
