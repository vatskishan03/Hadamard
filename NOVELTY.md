# Novelty

How we decided whether an order is new, and what we searched. Every claim of novelty in this repository means "new
to our knowledge" in the sense below. A literature search cannot show that no construction exists elsewhere;
corrections are welcome.

## Status labels

| Label | Meaning |
|---|---|
| **New** | The order is open in the reference table for its family (below), and no later construction was found. |
| **New, beyond the tables** | The order lies above the range of the reference tables. Where Table 1.53 of the Handbook covers it, the table gives t ≥ 3 (no Hadamard matrix of order 4n known), and no later construction was found. |
| **First explicit construction** | A table records the order as known, but we could not trace any construction behind the entry. |
| **Concurrent** | No published construction was found, but another group constructed the order independently at about the same time. We claim no priority. |
| **Known** | An earlier construction exists. The matrix here is independent of it. |

For skew and symmetric matrices the property is what is new: an ordinary Hadamard matrix of the same order may have
been known. Each order's README says so where it applies.

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
- searches of Zenodo (records matching "hadamard"), GitHub, the e-print sources of Hadamard-related arXiv papers
  (to 2026-10-08) and the math.CO listings (to 2026-10-12), the OEIS, MathOverflow, Balonin's catalogue and the Palomar
  registry.

The searches were made between 12 and 14 October 2026.

## Not checked

- The updated tables in J. Seberry, M. Yamada, *Hadamard Matrices: Constructions using Number Theory and Linear
  Algebra* (Wiley, 2020), seen only through the papers that cite them.
- J. Seberry, *Orthogonal Designs* (Springer, 2017), and Geramita–Seberry (1979).
- The lists in §7 of Miyamoto, J. Combin. Theory Ser. A 57 (1991), which Cati and Pasechnik could not verify either.
- The skew table of the Handbook (Table 1.51).
- MathSciNet, zbMATH, and the full texts of paywalled journals from 2025–2026.

A construction recorded only in one of these is the remaining risk for every claim.
