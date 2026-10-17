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
the order keeps the label First explicit construction even if a repaired route from ingredients published before 2026
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

## Later work

For every order we checked:

- Đoković's updates of the Hadamard tables (arXiv:0912.5091, arXiv:1008.2043), Đoković–Golubitsky–Kotsireas
  (arXiv:1301.3671), and the later papers of Đoković and coauthors on skew and symmetric matrices (among them
  arXiv:1511.08734, arXiv:2404.14375 and arXiv:2508.17141);
- SageMath, development version of October 2026: `hadamard_matrix(n, existence=True)`,
  `skew_hadamard_matrix(n, existence=True)` and the list `unknown_skew_hadamard` (we missed that it builds 31316 since
  the merge of October 3; corrected in [hadamard/31316](hadamard/31316/README.md#status));
- the constructions published in 2026:
  - L. Alpöge, P. Voinov and S. Reynolds-Haertle, twelve orders below 2000 (2026-08-12);
  - Schneider, order 2060 (2026-08-23);
  - the JD-Jones-ASES repositories Hadamard-T, Hadamard-M, Hadamard-B and Hadamard-formal;
  - Schumer, orders 3676, 5436, 6644 and 7852 (Zenodo, 2026-09-17, doi:10.5281/zenodo.22820997);
  - the deposits of Piza, Kline, Bedani, Karoui and Blackthorn;
  - Szöllősi, arXiv:2603.09572 (complex Hadamard matrices);
  - for the complex orders also Ito's papers on Hadamard groups (1994–2000), B. Schmidt, Williamson matrices and a
    conjecture of Ito's, Des. Codes Cryptogr. 17 (1999) 61–68, and R. Egan, Des. Codes Cryptogr. 87 (2019) 67–74,
    Theorem 4.5, which give every order 2N with 4N − 1 a prime power (found on 16 October 2026, after we had
    labelled 76 such orders New; corrected in [complex-hadamard](complex-hadamard/README.md#status));
  - Kharaghani, Mohammadian, Tayfeh-Rezaie and Zaitsev, arXiv:2605.08661 (near-Williamson matrices);
  - H. Verma, A symmetric conference matrix of order 86, arXiv:2610.08893 (6 October 2026; found by the audit of 17
    October 2026 and checked exactly here), which gives Williamson-type matrices of order 173 through Seberry–Yamada's
    Corollary 8.8(ii) ([NOTES.md](NOTES.md) §7);
  - for the T-matrices, S. London and I. Kotsireas, New Turyn-type sequences, Cryptogr. Commun. 17 (2025)
    1601–1610 (TT(40), TT(42), TT(44)), which give T-sequences of length 131 (found on 16 October 2026, after we
    had listed 131 as new; corrected in [t-matrices](t-matrices/README.md));
  - V. Zaitsev's page Quaternary Hadamard matrices (https://www.cs.uleth.ca/~zaitsevv/quaternary.html, last updated
    14 October 2026, read on 16 October 2026; work with Kharaghani, Mohammadian and Tayfeh-Rezaie), with matrices of
    122 complex orders, among them all 11 that we had labelled New; they are now labelled Concurrent. Its matrix of
    order 1142 also gives the real order 2284 by Turyn's doubling, so [2284](hadamard/2284) is Concurrent as well (no
    other order of the page doubles to an order of FRONTIER.md or of this repository);
  - Czifra, Matolcsi and Szöllősi, arXiv:2511.11230 (Butson matrices BH(n, q) for q = 6, 10 and 14, none quaternary),
    the quaternary Legendre pairs of arXiv:2609.04589 and arXiv:2609.25573 (lengths 42 to 80), and the skew-regular
    quaternary Hadamard matrices of order 1 + p² of arXiv:2510.06489 (complex orders,
    checked on 15 October 2026);
- for the fourteen orders from Paley's matrix in sixteen circulant blocks ([hadamard/PALEY_BLOCKS.md](hadamard/PALEY_BLOCKS.md)),
  checked on 16 October 2026: Seberry–Yamada, JCMCC 7 (1990) 97–137 (Theorem 4) and its Kobe proceedings version, the
  1992 survey (Theorems 3.10 and 8.1), Seberry–Yamada, JCMCC 13 (1993) 97–106 (Theorems 1 and 6, Corollary 2),
  Schmidt, Des. Codes Cryptogr. 17 (1999) 61–68, Schmidt–Tan, Des. Codes Cryptogr. 73 (2014) 105–119,
  Balonin–Đoković, arXiv:1508.00640, Barrera Acevedo–Ó Catháin–Dietrich, arXiv:1904.11460, Xia–Xia–Seberry–Wu,
  Discrete Appl. Math. 145 (2005) 498–504, and the literature on cocyclic and dicyclic Hadamard matrices. Both steps of
  the construction are published; we found no source that combines them. Yamada's paper of 1991 and Ito's papers of
  1994 were not accessible; we read their results through Seberry–Yamada (1993) and Schmidt (1999);
- for the 115 Paley-block orders beyond Table 4
  ([hadamard/PALEY_BLOCKS.md](hadamard/PALEY_BLOCKS.md#beyond-table-4)), checked on 17 October 2026:
  `hadamard_matrix(4n, existence=True)` of SageMath's development version for every order (Unknown except 16892 and
  25132); the Hadamard-related e-print texts of our arXiv cache for every n and 4n (no relevant mention);
  Cooper–Wallis, Miyamoto, Yamada–Kiyasu and Paley for each n with the T-matrix and Williamson-type orders known in
  2026 (the routes found are stated in the READMEs). London and Kotsireas, Cryptogr. Commun. 17 (2025), whose
  T-sequences of length 131 give eight of the orders by Cooper–Wallis, was not accessible;
- for the Miyamoto-block orders ([hadamard/MIYAMOTO_BLOCKS.md](hadamard/MIYAMOTO_BLOCKS.md)) and the orders recorded
  only through Corollary 8.8(i)
  ([hadamard/PALEY_BLOCKS.md](hadamard/PALEY_BLOCKS.md#orders-recorded-only-through-corollary-88i)), checked on 17
  October 2026: SageMath's development version for every order (Unknown for all; its Miyamoto test raises UnknownError
  for some, rerun with Unknown treated as False); Paley I and II, Miyamoto, Yamada–Kiyasu (with the skew-Hadamard
  matrices of the Đoković–Seberry census and of this repository) and the T × WT splits of Cooper–Wallis for each n;
  Seberry–Yamada, JCMCC 7 (1990), Theorems 4, 21, 34 and Corollaries 22, 35, its Kobe version (Theorems 1, 8, 21,
  Corollaries 9, 22) and the 1992 survey (Theorems 3.10, 8.1, 8.2, 8.10, Corollary 8.11; Tables A.1, A.2; §5, Yang's
  theorems): Miyamoto's matrix is an M-structure there, but is never multiplied by T-matrices, and Theorem 21 admits
  only t = 1 or even t; Seberry–Yamada, JCMCC 13 (1993), Theorems 6 and 7; Handbook §V.2, 2.45–2.49; Đoković,
  arXiv:0912.5091 and arXiv:1008.2043; Balonin–Đoković, arXiv:1508.00640 (negaperiodic Golay pairs and their list (31)
  of open quasi-Williamson orders); Đoković–Kotsireas, arXiv:1310.5773 (periodic Golay pairs); the e-print texts of our
  arXiv cache for Miyamoto's theorem combined with T-matrices or block multiplication (no source found), and for every n
  and 4n with a pattern that also matches numbers followed by a comma or period (the first grep missed comma-separated
  lists; rerun on 17 October 2026, finding only Đoković's claims of 7179 and 7323 in arXiv:1008.2043). Miyamoto's 1991
  paper (J. Combin. Theory Ser. A 57) and Yang's 1989 paper (Proc. Amer. Math. Soc. 107) were not accessible; we used
  Yang's multiplication as restated in the 1992 survey (Theorem 5.23);
- for the skew-Hadamard matrices of order 3(q + 1) ([skew/PALEY_SKEW.md](skew/PALEY_SKEW.md)), checked on 17 October
  2026: Seberry, Ars Combin. 6 (1978) (Summary 2, Theorems 2–7, Table 2); Seberry–Yamada, JCMCC 7 (1990) and its Kobe
  version, Utilitas Math. 40 (1991), the 1992 survey (§7, the families SI–SXI, Table 7.1 on the scan, Theorem 3.10, §8)
  and JCMCC 13 (1993); Wallis's lists of 1971 and 1973; Xia–Xia–Seberry, Discrete Appl. Math. 145 (2005), and the
  group's papers of 2005–2006; Handbook §V.1 with Table 1.51 (online copy); Koukouvinos–Stylianou (2008), first pages
  only; Cati–Pasechnik Table 5 and footnote 13; the census; Đoković, arXiv:2508.17141 (XXSW at order q + 1 only) and his
  earlier skew papers; SageMath's development source (`unknown_skew_hadamard`, and an emulation of the branches of
  `skew_hadamard_matrix`, Unknown for all 58 orders); the arXiv skew-Hadamard records (39) and Hadamard records of
  2025–2026, Zenodo, the JD-Jones-ASES repositories and Schumer's deposit. None gives skew matrices of order 3(q + 1) in
  general or of any of the 58 orders; Seberry–Yamada's Table 7.1 and Cati–Pasechnik's Table 5 record 663 and 879 as
  known without a construction. JCMCC 13 (1993), Theorem 1 (after Yamada 1991) contains the sixteen-block form with the
  transposed pair, so the orders are labelled New, but obtainable classically. Wallis–Whiteman (1972), as stated by
  Momihara–Xiang (arXiv:1801.08776), gives 7656 (relabelled Known). Not accessible: Whiteman, Discrete Math. 2 (1972);
  Haviland's BYU thesis (2023) and Anderson–Haviland–Holmes–Humphries–Magland, Graphs Combin. 39 (2023), on skew
  relative difference sets in dicyclic groups, whose abstract describes one infinite family built from Ito's and
  Schmidt's results (our matrices of orders 60 and 132 have automorphism groups of order 2, so they are not
  group-developed of that kind); Xia–Zuo–Lou–Xia, Trans. Comb. 13 (2024);
- searches of Zenodo (records matching "hadamard"), GitHub, the e-print sources of Hadamard-related arXiv papers
  (to 2026-10-08) and the math.CO listings (to 2026-10-12), the OEIS, MathOverflow, Balonin's catalogue and the Palomar
  registry.

The searches were made between 12 and 14 October 2026.

## Not checked

- J. Seberry, M. Yamada, *Hadamard Matrices: Constructions using Number Theory and Linear Algebra* (Wiley, 2020):
  its updated tables, and chapter 6 (§6.4 Multiplication Theorems Using M-structures, §6.5 Miyamoto's Theorem and
  Corollaries via M-structures), the most likely place for a prior statement of the Miyamoto-block route; only its
  contents and index were seen.
- J. Seberry, *Orthogonal Designs* (Springer, 2017), and Geramita–Seberry (1979).
- Miyamoto, J. Combin. Theory Ser. A 57 (1991): its Corollaries 2–5 on matrices of symmetric block form, to which
  Seberry and Yamada credit Table A.1's key wl (the Williamson-type matrices behind the Miyamoto-block orders of
  [hadamard/MIYAMOTO_BLOCKS.md](hadamard/MIYAMOTO_BLOCKS.md)), and the lists in its §7, which Cati and Pasechnik could
  not verify either. We know them only through secondary sources: the zbMATH review (Zbl 0733.05017) and the
  transcription in the JD-Jones-ASES Hadamard-M note, by which list (2) contains 127 (unsupported by Miyamoto's
  Corollary 4, per that note), so that list (3) may record 1527 = 3 · 509.
- The skew table of the Handbook (Table 1.51).
- MathSciNet, zbMATH, and the full texts of paywalled journals from 2025–2026.

A construction recorded only in one of these is the remaining risk for every claim.
