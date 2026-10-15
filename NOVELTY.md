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
  `skew_hadamard_matrix(n, existence=True)` and the list `unknown_skew_hadamard`;
- the constructions published in 2026:
  - L. Alpöge, P. Voinov and S. Reynolds-Haertle, twelve orders below 2000 (2026-08-12);
  - Schneider, order 2060 (2026-08-23);
  - the JD-Jones-ASES repositories Hadamard-T, Hadamard-M, Hadamard-B and Hadamard-formal;
  - Schumer, orders 3676, 5436, 6644 and 7852 (Zenodo, 2026-09-17, doi:10.5281/zenodo.22820997);
  - the deposits of Piza, Kline, Bedani, Karoui and Blackthorn;
  - Szöllősi, arXiv:2603.09572 (complex Hadamard matrices);
  - Kharaghani, Mohammadian, Tayfeh-Rezaie and Zaitsev, arXiv:2605.08661 (near-Williamson matrices);
  - Czifra, Matolcsi and Szöllősi, arXiv:2511.11230 (Butson matrices BH(n, q) for q = 6, 10 and 14, none quaternary),
    and the quaternary Legendre pairs of arXiv:2609.04589 and arXiv:2609.25573 (lengths 42 to 80) (complex orders,
    checked on 15 October 2026);
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
