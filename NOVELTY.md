# Novelty

How we decided whether an order is new, and what we searched. Every claim of novelty in this repository means "new
to our knowledge" in the sense below. A literature search cannot show that no construction exists elsewhere;
corrections are welcome.

## Status labels

| Label | Meaning |
|---|---|
| **New** | The order is open in the reference table for its family (below), and no later construction was found. |
| **New, but obtainable classically** | As New, but every ingredient (theorems and objects) was published before 2026; only their combination is new. |
| **New, beyond the tables** | The order lies above the range of Table 4 (for complex orders, Table 11.2). Where Table 1.53 of the Handbook covers it, the table gives t ≥ 3 (no Hadamard matrix of order 4n known), and no later construction was found. |
| **New, beyond the tables, but obtainable classically** | As New, beyond the tables, but every ingredient (theorems and objects) was published before 2026; only their combination is new. |
| **First explicit construction** | A table records the order as known, but we could not trace any construction behind the entry. |
| **Concurrent** | No published construction was found, but another group constructed the order independently at about the same time. We claim no priority. Also used when published or posted objects of another group give the order by a short, explicit step (a monomial equivalence and a classical map); the order's README says which. |
| **Known** | An earlier construction exists. Where the matrix here is independent of it, the README says so. |

New is used when the construction that gives the symmetric (skew, complex) matrix needs a proof not found in any
source read, even if its inputs are classical; New, but obtainable classically when a published theorem gives the
matrix, or gives it after a substitution that needs no new proof.

For skew and symmetric matrices the property is what is new: an ordinary Hadamard matrix of the same order may have
been known. Each order's README says so where it applies.

When a table's entry rests on a printed route that fails (Corollary 8.8(i) of Seberry–Yamada,
[NOTES.md](NOTES.md) §7), the order is labelled First explicit construction even if a repaired route from ingredients
published before 2026 exists; the README names that route where we found one.

## Reference tables

| Family | Table | Range |
|---|---|---|
| Hadamard | M. Cati, D. V. Pasechnik, *A database of constructions of Hadamard matrices*, arXiv:2411.18897 (v1 of 2024-11-28 and v2 of 2025-08-30 have the same rows), Table 4 | odd n ≤ 2999 |
| Hadamard | *Handbook of Combinatorial Designs*, 2nd ed. (2007), Table 1.53 | odd n < 10000 |
| Skew | D. Ž. Đoković, J. Seberry, census of skew and symmetric Hadamard matrices of order 4v (2025-05-28), and Cati–Pasechnik Table 5 | odd v < 2500 |
| Symmetric | the same census: 595 undecided v (three of them settled by Đoković, arXiv:2508.17141, and nine by published theorems that its series omit: Seberry–Balonin's Corollary 1 and Spence's Corollary 3); the census itself warns that it may have missed some orders; Handbook Table V.1.52 | odd v < 2500 |
| Complex | J. Seberry, M. Yamada, *Hadamard matrices, sequences, and block designs*, in *Contemporary Design Theory* (Wiley, 1992), Table 11.2 | odd N < 1000 |
| Williamson-type | Seberry–Yamada (1992), Table A.1 | |
| T-matrices | Seberry–Yamada (1992), §3 | t ≤ 210 |

We also read the older tables these replace: Seberry–Yamada (1992), Tables 7.1 and A.2, and, for most skew orders,
Koukouvinos–Stylianou, Discrete Math. 308 (2008), Tables 1 and 3. A README says so where one of these was not read.
The Handbook tables were read from a copy available online and not checked against the printed edition.

## Sources searched

The searches were made between 12 and 22 October 2026; each section below gives its own dates.

**For every order.**

- Đoković's updates of the Hadamard tables (arXiv:0912.5091, arXiv:1008.2043), Đoković–Golubitsky–Kotsireas
  (arXiv:1301.3671), and the papers of Đoković and coauthors on skew and symmetric matrices, among them
  arXiv:1511.08734, arXiv:2404.14375, arXiv:2508.17141, Đoković–Kotsireas, Math. Comput. Sci. 12 (2018) 373–388, and
  Abuzin–Balonin–Đoković–Kotsireas, Inform. Control Syst. 2019, no. 5, 2–9 (Table 1, skew-propus families).
- SageMath, development version of October 2026: `hadamard_matrix(n, existence=True)`,
  `skew_hadamard_matrix(n, existence=True)` and the list `unknown_skew_hadamard`.
- The constructions published in 2026: Alpöge, Voinov and Reynolds-Haertle (2026-08-12); Schneider (2026-08-23); the
  JD-Jones-ASES repositories Hadamard-T, Hadamard-M (https://github.com/JD-Jones-ASES/Hadamard-M), Hadamard-B and
  Hadamard-formal; Schumer, orders 3676, 5436, 6644 and 7852 (Zenodo, 2026-09-17, doi:10.5281/zenodo.22820997); the
  deposits of Piza, Kline, Bedani, Karoui and Blackthorn; Szöllősi, arXiv:2603.09572; Kharaghani, Mohammadian,
  Tayfeh-Rezaie and Zaitsev, arXiv:2605.08661 (near-Williamson matrices); H. Verma, A symmetric conference matrix of
  order 86, arXiv:2610.08893 (2026-10-06), which gives Williamson-type matrices of order 173 through Seberry–Yamada's
  Corollary 8.8(ii) ([NOTES.md](NOTES.md) §7); R. Egan, I. S. Kotsireas, S. London and P. Ó Catháin, Recent progress
  in the construction of Hadamard matrices, arXiv:2610.27175 (2026-10-19), whose Table 5 lists the 95 orders up to
  10000 at which no matrix was known to its authors.
- Zenodo (records matching "hadamard"), GitHub, the e-print sources of Hadamard-related arXiv papers (to 2026-10-08)
  and the math.CO listings (to 2026-10-12), the OEIS, MathOverflow, Balonin's online catalogue and the Palomar
  registry.
- The GitHub repository Middletownbooks/Hadamard (https://github.com/Middletownbooks/Hadamard), found on 19 October
  2026 and read again on 20 October (see Corrections).

**Complex orders.** Ito's papers on Hadamard groups (1994–2000); B. Schmidt, Des. Codes Cryptogr. 17 (1999) 61–68;
R. Egan, Des. Codes Cryptogr. 87 (2019) 67–74, Theorem 4.5, which gives every order 2N with 4N − 1 a prime power;
V. Zaitsev's page Quaternary Hadamard matrices (https://www.cs.uleth.ca/~zaitsevv/quaternary.html, with Kharaghani,
Mohammadian and Tayfeh-Rezaie, last updated 14 October 2026), with matrices of 122 orders 2N; Czifra, Matolcsi and
Szöllősi, arXiv:2511.11230 (Butson matrices BH(n, q) for q = 6, 10 and 14, none quaternary); the quaternary Legendre
pairs of arXiv:2609.04589 and arXiv:2609.25573 (lengths 42 to 80); the skew-regular quaternary Hadamard matrices of
order 1 + p² of arXiv:2510.06489; Đoković's Williamson-type matrices of orders 509, 857 and 883 (arXiv:1008.2043).

**T-matrices.** S. London and I. Kotsireas, New Turyn-type sequences, Cryptogr. Commun. 17 (2025) 1601–1610 (TT(40),
TT(42), TT(44)), which give T-sequences of length 131. For the orders 97 and 109 of the form 6m + 1, Xia, Xia, Seberry
and Qin, Far East J. Math. Sci. 101 (2017) 1731–1749 (abstract only: concrete constructions for m ≤ 13, m ≠ 9, that is
up to order 79); Lou, Xia and Xia, Adv. Appl. Discrete Math. 21 (2019) 163–181, could not be read (see Not checked).

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
recorded only through Corollary 8.8(i). Seberry–Yamada, JCMCC 7 (1990), Theorems 4, 21, 34 and Corollaries 22, 35, its
Kobe version (Theorems 1, 8, 21, Corollaries 9, 22) and the 1992 survey (Theorems 3.10, 8.1, 8.2, 8.10, Corollary
8.11; Tables A.1, A.2; §5, Yang's theorems): Miyamoto's matrix is an M-structure there but is never multiplied by
T-matrices, and Theorem 21 admits only t = 1 or even t. Seberry–Yamada, JCMCC 13 (1993), Theorems 6 and 7;
Handbook §V.2, 2.45–2.49; Đoković, arXiv:0912.5091 and arXiv:1008.2043; Balonin–Đoković, arXiv:1508.00640
(negaperiodic Golay pairs, and their list (31) of open quasi-Williamson orders); Đoković–Kotsireas, arXiv:1310.5773
(periodic Golay pairs); SageMath for every order (Unknown for all); Paley I and II, Miyamoto, Yamada–Kiyasu (with the
skew-Hadamard matrices of the census and of this repository) and the T × WT splits of Cooper–Wallis for each n; the
e-print texts for Miyamoto's theorem combined with T-matrices or block multiplication, and for every n and 4n (only
Đoković's claims of 7179 and 7323 in arXiv:1008.2043). Yang's multiplication is used as restated in the 1992 survey
(Theorem 5.23).

**Skew-Hadamard matrices of order 3(q + 1)** ([skew/PALEY_SKEW.md](skew/PALEY_SKEW.md)). Seberry, Ars Combin. 6
(1978) (Summary 2, Theorems 2–7, Table 2); Seberry–Yamada, JCMCC 7 (1990) and its Kobe version, Utilitas Math. 40
(1991), the 1992 survey (§7, the families SI–SXI, Table 7.1, Theorem 3.10, §8) and JCMCC 13 (1993); Wallis's lists
of 1971 and 1973; Xia–Xia–Seberry, Discrete Appl. Math. 145 (2005), and the group's papers of 2005–2006; Handbook §V.1
with Table 1.51 (online copy); Koukouvinos–Stylianou (2008), first pages only; Cati–Pasechnik Table 5 and footnote
13; the census; Đoković, arXiv:2508.17141 (XXSW at order q + 1 only) and his earlier skew papers; SageMath
(`unknown_skew_hadamard`, and the cases of `skew_hadamard_matrix`, Unknown for all 58 orders); the arXiv
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
of the 37 orders of PROPUS_WEIGHING.md and YAMADA_KIYASU.md, except that Spence (1977), Corollary 3, gives symmetric
Hadamard matrices of orders 2812 and 4732 by another construction, so those two are labelled Known. The eight
equivalences were found with nauty and checked exactly; no script for them is included in this repository. Of the 27
matrices on the page whose order 2v has v among the census values still undecided when the check was made, nauty finds
no symmetric form for 20: QH130, QH214, QH238, QH326, QH426, QH470, QH490, QH494, QH518, QH534, QH538, QH566, QH610,
QH686, QH754, QH806, QH890, QH1142, QH1462 and QH1526 (QH538 has 16 side-swaps, none of which gives a symmetric
matrix). QH1018 has 16 side-swaps, none of which gives one either, checked without conjugation only. QH1306 and QH1466
timed out, and QH1558, QH1690, QH1706 and QH1882 were not tested. QH266, QH366 and QH1714 (v = 133, 183 and 857, which
this repository settles by the orders 532, 732 and 3428) were not tested.

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
**Known**: Spence, An infinite family of Williamson matrices, J. Austral. Math. Soc. Ser. A 24 (1977) 252–256,
Corollary 3, gives symmetric Hadamard matrices of order 2qʳ(q + 1) for every prime power q ≡ 1 (mod 4) and r ≥ 0;
these five orders are r = 1 with q = 17, 25, 29, 49 and 61, and it also gives 2812 (q = 37) and 4732 (q = 13, r = 2).
The census lists all seven v as undecided; its series omit the corollary. Also read: Seberry–Yamada (1992), §11
(Kharaghani–Seberry's Theorem 11.14 and Corollaries 11.15 and 11.20, and Lemma 11.3, Turyn's map); Seberry,
J. Austral. Math. Soc. A 29 (1980) 235–242; Seberry–Balonin (2017); Balonin–Đoković, arXiv:1803.08787; Đoković,
arXiv:2101.05429 and 2508.17141; Farouk–Wang, Filomat 34 (2020) 815–834 (Scarpis–Đoković matrices of order
2q(q + 1), not symmetric).

**Symmetric order 260** ([symmetric/260](symmetric/260)). Checked on 17 October 2026: the census; Đoković,
arXiv:2101.05429 (65 first among the open v < 200, and the Handbook's Table 1.52 list as he quotes it);
Seberry–Balonin (2015, 2017); Balonin–Đoković, arXiv:1803.08787; Verma, arXiv:2609.25573 (the pair only); Lebedev,
arXiv:2609.04589; Kotsireas–Winterhof and Jedwab–Pender as cited by Verma; Zaitsev's page (QH130, not equivalent to a
symmetric matrix); web searches for a symmetric Hadamard matrix of order 260.

**Symmetric conference matrices of order 66** ([search/NOTES.md](search/NOTES.md) §72). Checked on 19 October 2026:
Gritsenko, On strongly regular graph with parameters (65, 32, 15, 16), arXiv:2102.05432 (2021), constructs an
srg(65, 32, 15, 16) and so a symmetric conference matrix of order 66; Brouwer's table of strongly regular graphs lists
(65, 32, 15, 16) as existing, and OEIS A000952 contains 66. Goryainov, Haemers, Konstantinova and Li
(arXiv:2512.16653) name 86 as the smallest open order, which Verma's matrix (arXiv:2610.08893, 2026-10-06) settles.
Order 66 is therefore not open, and the order-66 case of Spence's recursion (v = 2145) is not a census value. Our
enumeration under two groups gives 7 classes of these matrices; two of them, and their 16 descendants
srg(65, 32, 15, 16), are not in McKay's database of srg(65, 32, 15, 16) (32 graphs, marked as possibly incomplete); no
other enumeration was located.

**Near-Williamson matrices of order 111** ([near-williamson](near-williamson)). Checked in October 2026: Kharaghani,
Mohammadian, Tayfeh-Rezaie and Zaitsev, arXiv:2605.08661v2, whose orders are every odd w ≤ 69 and 73, 79, 91, 93, 95,
97, 103, 109, 127, 133 and 181 (the last group in their database); Seberry–Yamada (1992), Table A.1, which records
only Williamson-type matrices of order 222 for 111 (t = 3). None gives near-Williamson or Williamson-type matrices of
order 111, which is therefore New.

**Order 6796** ([hadamard/6796](hadamard/6796)). Checked on 20 October 2026, Zaitsev's page again on 22 October:
- *Tables and recursions.* Cati–Pasechnik, arXiv:2411.18897 v1 and v2 (Table 4: 1699(3)); Handbook Table 1.53;
  Đoković, arXiv:0912.5091 (the set Δ and four further values) and arXiv:1008.2043 (42 values);
  Đoković–Golubitsky–Kotsireas, arXiv:1301.3671; SageMath development source of 18 October 2026, read case by case,
  not run (Unknown).
- *Classical families.* Yamada 1989 (Theorem 5) and Seberry–Yamada 1992 (key y3, §11). For conference matrices of
  order 850: Seberry–Whiteman 1980 (via Egan, arXiv:2309.07522), Handbook §V.6 (BGW and GC families) and Et-Taoui,
  arXiv:1409.5720. Fallon–Iverson, arXiv:2610.09169: a complex conference matrix of order 850, not quaternary in any
  equivalent form.
- *Recent work.* The 2026 sources of "For every order" and Middletownbooks/Hadamard (default branch, last commit 20
  October, a change to its README only; a second branch holds search logs); Egan, Kotsireas, London and Ó Catháin,
  arXiv:2610.27175 (19 October 2026), Table 5, which lists 1699 as open; Zaitsev's page (last updated 14 October,
  unchanged on 22 October; orders below 2000); the October 2026 arXiv listings of math.CO, cs.DM and cs.IT, with the
  e-print sources of 20 related papers; Zenodo, GitHub (repositories and code) and web searches for 6796, 3398 and
  1699.
- *Finding.* None gives a Hadamard matrix of order 6796, a quaternary Hadamard matrix of order 3398 or a quaternary
  conference matrix of order 850, or states the complex Yamada–Kiyasu template.
- *Not read.* The Seberry–Yamada book of 2020, Seberry–Whiteman 1980 itself and Miyamoto's §7 lists.
- *Labels.* Hadamard order 6796: New. Complex order 3398 (N = 1699, above the range N < 1000 of the reference table):
  New, beyond the tables.

**Orders 3436 and 4204** ([hadamard/3436](hadamard/3436), [hadamard/4204](hadamard/4204)). Checked on 20 October
2026, Zaitsev's page again on 22 October:
- *Tables and recursions.* Cati–Pasechnik v1 and v2 (Table 4: 859(3), 1051(3)); Handbook Table 1.53 (t = 3 for
  both); Đoković, arXiv:0912.5091 and arXiv:1008.2043; Đoković–Golubitsky–Kotsireas, arXiv:1301.3671; SageMath
  development source of 20 October 2026, read case by case, not run (Unknown for both); Seberry–Yamada 1992,
  Tables A.1 and 11.2.
- *Conference matrices of orders 430 and 526.* Belevitch's condition; the generalized Paley, Berman, BGW and GC,
  Seberry–Whiteman (via Egan, arXiv:2309.07522) and Et-Taoui families. Fallon–Iverson's ancillary library
  (arXiv:2610.09169) has complex conference matrices of both orders, not quaternary in any equivalent form.
  Miyamoto's Corollary 4 (through the JD-Jones-ASES Hadamard-M note) needs Cᵀ = C.
- *Recent work.* Zaitsev's page (last updated 14 October, unchanged on 22 October; 1718 on its unknown list);
  Middletownbooks/Hadamard (default branch, last commit 20 October, a change to its README only, which lists 859 as
  open; its list of open orders below 4000 ends at v = 971; a second branch holds search logs); Egan, Kotsireas,
  London and Ó Catháin, arXiv:2610.27175 (19 October 2026), Table 5, which lists 859 and 1051 as open; the
  JD-Jones-ASES repositories (unchanged); Schumer's deposit; Zenodo; the arXiv listings of 20 October and the arXiv
  API to 19 October; the e-print sources of 20 papers of 2026; GitHub code and repository search; Palomar; Balonin's
  catalogue; web searches for 3436, 4204, 1718, 2102, 859 and 1051.
- *Finding.* None gives a Hadamard matrix of order 3436 or 4204, a quaternary Hadamard matrix of order 1718 or 2102,
  or a quaternary conference matrix of order 430 or 526.
- *Not read.* As for 6796.
- *Labels.* Hadamard orders 3436 and 4204: New. Complex order 1718 (N = 859): New; complex order 2102 (N = 1051, above
  the range N < 1000 of the reference table): New, beyond the tables.

**Order 2572** ([hadamard/2572](hadamard/2572)). Checked on 21 October 2026:
- *Tables and recursions.* Cati–Pasechnik v1 and v2 (Table 4: 643(3); its §6 sentence naming 643 concerns the skew
  table); Đoković, arXiv:0912.5091 and arXiv:1008.2043; Đoković–Golubitsky–Kotsireas, arXiv:1301.3671; SageMath
  development source of 21 October 2026, read case by case, not run (Unknown); Seberry–Yamada 1992, Tables A.1 and
  11.2; Brouwer's table (srg(321, 160, 79, 80) does not exist, srg(1285, 642, 320, 321) undecided).
- *Conference matrices of order 322.* As for orders 430 and 526. Fallon–Iverson's ancillary library has a complex
  one, not quaternary in any equivalent form.
- *Recent work.* Egan, Kotsireas, London and Ó Catháin, arXiv:2610.27175 (19 October 2026), Table 5, which lists 643
  as open; Zaitsev's page (last updated 14 October, unchanged on 21 October; 1286 on its unknown list);
  Middletownbooks/Hadamard (its open list of 20 October has 643; its screen of 6 October covers Goethals–Seidel and
  Legendre-quadruple structures over ℤ₆₄₃); the arXiv listings to 21 October; GitHub and web searches for 2572, 1286
  and 643.
- *Finding.* None gives a Hadamard matrix of order 2572, a quaternary Hadamard matrix of order 1286 or a quaternary
  conference matrix of order 322.
- *Not read.* As for 6796; the Handbook's Table 1.53 row for 643, which was not read directly.
- *Labels.* Hadamard order 2572: New. Complex order 1286 (N = 643): New.

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

**Skew order 476 and symmetric order 7612** ([skew/476](skew/476), [symmetric/7612](symmetric/7612)). Checked on
22–23 October 2026:
- *Tables.* The Đoković–Seberry census (119 and 1903 undecided, for skew and symmetric matrices); Cati–Pasechnik
  Table 5 (119(4), so the smallest skew order known with odd part 119 was 1904) and Table 3 (no skew construction at
  119); SageMath's `unknown_skew_hadamard` (it contains 476) and the development thread of SageMath issue 34848.
- *Recent work.* The survey of Egan, Kotsireas, London and Ó Catháin, arXiv:2610.27175 (smallest open skew order 356,
  no construction named at 476 or 7612; its entry “4n for odd n < 250” under known skew orders contradicts its own
  356 and the census); the 2026 sources of "For every order"; web searches for skew 476 and symmetric 7612.
- *Finding.* None gives a skew-Hadamard matrix of order 476 or 952 or a symmetric Hadamard matrix of order 7612.
- *Not read.* As for 6796; work after the census was found only through the survey and web searches.

**Skew orders 1528, 1892, 2216, 2392 and 4712** ([skew/1528](skew/1528), [skew/1892](skew/1892),
[skew/2216](skew/2216), [skew/2392](skew/2392), [skew/4712](skew/4712)). Checked on 23 October 2026:
- *Tables.* Cati–Pasechnik Table 5, v1 and v2 (191(226), 277(5), 299(4), 473(4), 589(5)); the Đoković–Seberry census
  (191, 277, 299, 473 and 589 undecided, which concerns orders 4v and is not settled here except for 473); SageMath's
  development source of 21 October 2026, `skew_hadamard_matrix` read branch by branch, not run (Paley, Spence,
  complementary difference sets, Spence 1975, the orthogonal-design branch, doubling, Williamson's product and the
  stored skew difference families give none of the five orders, so they are Unknown).
- *Classical routes* (each README lists them): Paley I needs 1527 = 3 · 509, 1891 = 31 · 61, 2215 = 5 · 443,
  2391 = 3 · 797 or 4711 = 7 · 673 to be a prime power; Szekeres and Momihara–Xiang 2(qᵉ + 1) need 763, 945, 1107, 1195
  or 2355; Wallis–Whiteman and Momihara–Xiang 2ᵘ(qᵉ + 1) need 381 or 190, 472, 553 or 276, 597, 1177; Spence
  4(q² + q + 1) and Turyn (p − 1)ᵘ + 1 do not apply; Seberry's mn(n − 1) needs a skew-type OD(44; 1, 1, 42) for 1892
  (there are no G-matrices of order 11) and an OD(mn; 1, m, mn − m − 1) with n = 2 otherwise; doubling needs skew orders 764,
  946 (not a multiple of 4), 1108, 1196 or 2356; skew × amicable products have a factor that is not a Hadamard order,
  apart from the doubling case 2 · (n/2).
- *Negacyclic arrays.* Finlayson and Seberry, Orthogonal designs from negacyclic matrices, Australas. J. Combin. 30
  (2004) 319–330 (the negacyclic Goethals–Seidel array, Theorem 4, with blocks of odd order 3, 5, 7: orders 12, 20,
  28); Balonin and Đoković, arXiv:1508.00640 (negaperiodic Golay pairs, two blocks, not skew); Kotsireas and
  Koukouvinos, New skew-Hadamard matrices via computational algebra, Australas. J. Combin. 41 (2008) 235–248 (orders
  60 to 416).
- *Recent work.* The survey of Egan, Kotsireas, London and Ó Catháin, arXiv:2610.27175 (no skew construction at
  these orders; 191 appears only in its table of open ordinary orders 4 · 5p and 4 · 9p); Đoković, arXiv:2301.02751
  (191 among the undecided skew v); the GitHub repositories of the 2026 Hadamard announcements (munimtechnologies,
  Middletownbooks: no skew orders); web searches for the five orders and for negacyclic or quaternary skew
  Goethals–Seidel arrays.
- *Finding.* None gives a skew-Hadamard matrix of order 1528, 1892, 2216, 2392 or 4712, or a multiplier search over
  negacyclic Goethals–Seidel arrays.
- *Not read.* As for 6796.

**Skew order 2156** ([skew/2156](skew/2156)). Checked on 23 October 2026:
- *Tables.* The Đoković–Seberry census (539 undecided); Cati–Pasechnik Table 5, v1 and v2 (539(4)) and Table 4 (539 not
  listed, so an ordinary Hadamard matrix of order 2156 is known); SageMath's development source of 21 October 2026,
  `skew_hadamard_matrix` read branch by branch, not run (none of its branches gives 2156).
- *Classical routes.* Paley I needs 2155 = 5 · 431, Szekeres and Momihara–Xiang 2(qᵉ + 1) need 1077 = 3 · 359,
  Wallis–Whiteman 4(q + 1) needs 538; Spence 4(q² + q + 1), Turyn (p − 1)ᵘ + 1 and Seberry's mn(n − 1) (only n = 2)
  do not apply; 2156 ≡ 4 (mod 8) rules out doubling, and every product of a skew and an amicable order has a factor
  that is not a Hadamard order. Momihara and Xiang's skew Hadamard difference families (arXiv:1801.08776) live on
  F_{qᵉ}, and 539 is not a prime power.
- *Recent work.* The survey of Egan, Kotsireas, London and Ó Catháin, arXiv:2610.27175; Đoković, arXiv:2301.02751 and
  arXiv:1511.08734; web searches for the order.
- *Finding.* None gives a skew-Hadamard matrix of order 2156.
- *Not read.* As for 6796.

**Skew order 1540** ([skew/1540](skew/1540)). Checked on 23 October 2026:
- *Tables.* The Đoković–Seberry census (385 undecided); Cati–Pasechnik Table 5 (385(3)) and Table 4 (385 not listed,
  so an ordinary Hadamard matrix of order 1540 is known); SageMath's development source of 21 October 2026,
  `skew_hadamard_matrix` and its difference-set tables read, not run (no entry for 385 or 1540).
- *Classical routes.* Paley I needs 1539 = 3⁴ · 19; the two-block families of order 2(q + 1) need q = 769 ≡ 5 (mod 8)
  (769 ≡ 1); Momihara–Xiang's families (arXiv:1801.08776) give 2ᵘ⁻¹(qᵉ + 1) with q ≡ 2ᵘ + 1 (mod 2ᵘ⁺¹), which would need
  qᵉ = 769 with u = 2 or qᵉ = 384 with u = 3; Wallis–Whiteman 4(q + 1), Spence 4(q² + q + 1), Turyn and the family
  3(q + 1) do not apply; every product of a skew and an amicable order has a factor that is not a Hadamard order.
- *Recent work.* The survey of Egan, Kotsireas, London and Ó Catháin, arXiv:2610.27175 (no skew matrix of this order);
  Đoković, arXiv:2301.02751 and arXiv:1511.08734; Middletownbooks/Hadamard; web searches for the order.
- *Finding.* None gives a skew-Hadamard matrix of order 1540.
- *Not read.* As for 6796.

**Skew orders 980 and 3124** ([skew/980](skew/980), [skew/3124](skew/3124)). Checked on 24 October 2026:
- *Tables.* The Đoković–Seberry census (245 and 781 undecided; 245 also on Đoković's list of 2007, arXiv:0706.1973);
  Cati–Pasechnik Table 5 (245(4), 781(3)) and Table 4 (neither listed, so ordinary Hadamard matrices of orders 980 and
  3124 are known); SageMath's development source of 21 October 2026 lists 980 among its unknown skew orders.
- *Classical routes.* Paley I needs 979 = 11 · 89 and 3123 = 3² · 347; the orders 2(q + 1) need q = 489 = 3 · 163 and
  1561 = 7 · 223; Momihara–Xiang's 2ᵘ⁻¹(qᵉ + 1) would need qᵉ = 489 or 244, and 1561 or 780; Wallis–Whiteman 4(q + 1),
  Spence 4(q² + q + 1), Turyn and the family 3(q + 1) do not apply; every product of a skew and an amicable order has a
  factor that is not a Hadamard order.
- *Recent work.* The survey of Egan, Kotsireas, London and Ó Catháin, arXiv:2610.27175 (neither order); Middletownbooks/
  Hadamard; web searches for both orders.
- *Finding.* None gives a skew-Hadamard matrix of order 980 or 3124.
- *Not read.* As for 6796.

**Orders 8696 and 12952** (checked on 24 October 2026). Seberry–Yamada (1992), Tables A.1 (keys wi and wl), A.2,
Lemma 8.12 and Corollaries 8.8, 8.14 and 8.18; Cati–Pasechnik's Table 4 (1087(3), 1619(3)); the survey of Egan,
Kotsireas, London and Ó Catháin, arXiv:2610.27175 (its Table 5 lists 8696 as open, after excluding Williamson-type
matrices of order 543); Đoković, arXiv:0912.5091 and arXiv:1008.2043; SageMath. 8696 is **Known** (A.1's key wi is a
valid route and A.2 records the order); 12952 is a **First explicit construction** (key wl is Corollary 8.8(i)).

**The conference matrix of order 118, symmetric 940 and order 11292** (checked on 24–25 October 2026). Brouwer's table
of strongly regular graphs (read on 17 October 2026: SRG(117, 58, 28, 29) exists, J. Y. Z. Ho, email, October 2026,
not public); Mathon, Canad. J. Math. 30 (1978) (order pq² + 1 for p = q + 2, which excludes 117 = 13 · 3²);
Seberry–Whiteman (1988); the census (235 undecided) and Seberry–Balonin (2017); Cati–Pasechnik's Table 4 (2823 open);
the survey of Egan et al. (2026); SageMath; the 2026 sources listed above. No published conference matrix of order 118,
symmetric Hadamard matrix of order 940 or Hadamard matrix of order 11292 was found. Ho's graph gives 940 through
Seberry–Balonin's Corollary 1, so 940 is **Concurrent**, like the conference matrix itself; 11292 needs the Miyamoto
blocks (or the repeated-block arrangement of [NOTES.md](NOTES.md) §7) on top of it and is **New**.

**Orders 2092, 3732, 4492, 5020 and 5932** (rechecked on 23 October 2026). The survey of Egan, Kotsireas, London and
Ó Catháin (arXiv:2610.27175, 19 October 2026) states that 3732, 4492, 5020 and 5932 "were found by the present
authors", by the Goethals–Seidel array with multiplier groups of orders 31, 33, 50 and 39, the groups our matrices use,
and lists 2092 under Cooley (Middletownbooks/Hadamard). That repository committed matrices of orders 2092 and 3732 on
14 October 2026 and calls ours independent. Ours were committed on 11 October (2092, 4492, 5932) and 13 October (3732,
5020). The survey cites this repository at commit 2d45dbb (18 October), which already held all five, for 2148, 2284,
2636 and 3964. With three groups working within days of each other, the five orders are **Concurrent**; neither other
matrix of order 3732 is stated to be skew, so the skew-Hadamard matrix of order 3732 stays **New**. The same survey
settles 4 · 1447 = 5788 (its Theorem 7) and lists six open orders 8p ([FRONTIER.md](FRONTIER.md)).

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
- L. Lou, M. Xia, T. Xia, On T-matrices of order 6m + 1, Adv. Appl. Discrete Math. 21 (2019) 163–181: it could not be
  read and may contain T-matrices of orders 97 = 6 · 16 + 1 and 109 = 6 · 18 + 1. Xia, Xia, Seberry and Qin, Far East
  J. Math. Sci. 101 (2017) 1731–1749, give concrete constructions of T-matrices of order 6m + 1 only for m ≤ 13,
  m ≠ 9. 197 is not of the form 6m + 1.
- Kharaghani and Seberry, Congr. Numer. 75 (1990); Seberry and Wallis, LNM 292, pp. 349–353; Seberry and Zhang (1993):
  sources for symmetric matrices of order 2n(n − 1), not read.
- The skew table of the Handbook (Table 1.51) beyond the online copy.
- MathSciNet, zbMATH, and the full texts of paywalled journals from 2025–2026.

A construction recorded only in one of these is the remaining risk for every claim.

## Corrections

Statements and labels of the first release (14 October 2026) that have changed:

- Complex orders 142, 154, 222, 262, 286, 310, 370 and 742: **Known**. For each, 4N − 1 is prime, and Ito's dicyclic
  Hadamard groups with Egan's Theorem 4.5 give the order ([complex-hadamard/PALEY.md](complex-hadamard/PALEY.md)).
- Complex orders 238, 302, 326, 382, 494, 518, 534 and 754: **Concurrent**, with Zaitsev et al.
- 2284: **Concurrent**, through Zaitsev et al.'s quaternary matrix of order 1142 and Turyn's doubling.
- 2092, 3732, 4492, 5020 and 5932: **Concurrent** (23 October), with Cooley (2092, 3732) and with Egan, Kotsireas,
  London and Ó Catháin (3732, 4492, 5020, 5932); the skew-Hadamard matrix of order 3732 stays New.
- 2292: **Known**. The GitHub repository Middletownbooks/Hadamard published an independent matrix of this order on
  9 October 2026, before ours of 14 October; we learned of it on 19 October.
- 31316: **Known**; SageMath's development version builds it since 3 October 2026.
- 19028, 21748 and 39260: **New, beyond the tables, but obtainable classically**.
- Skew 5700 and 9900, symmetric 5700: **New, but obtainable classically**. The amicable Hadamard matrices of orders
  52, 76 and 100 behind them are first explicit, not new: they come from the skew-propus families of Abuzin, Balonin,
  Đoković and Kotsireas (2019).
- T-matrices of order 131: **Known**, from London and Kotsireas's TT(44) (2025).
- 5272, 7464, 18344, 36688 and 73376: the first release called these orders open and settled by doubling the matrices
  of orders 2636, 3732 and 9172. They were not open. Complex Golay pairs of lengths x and y give a Hadamard matrix of
  order 8(x + y) (Craigen, Holzmann and Kharaghani, Discrete Math. 252 (2002) 73–89); 659 = 11 + 648 and
  2293 = 5 + 2288 give 5272 and 18344, and doubling gives 36688 and 73376. Agaian's product (Lecture Notes in Math.
  1168, 1985) of Hadamard matrices of orders 12 and 1244, the latter by Alpöge, Voinov and Reynolds-Haertle (2026),
  gives 7464 = 8 · 933. Egan, Kotsireas, London and Ó Catháin (arXiv:2610.27175) list none of them as open. The
  matrices of orders 2636, 3732 and 9172 are unaffected.
