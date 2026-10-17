# New Hadamard matrices

Explicit Hadamard, skew-Hadamard and symmetric Hadamard matrices of orders that had no known construction, and complex
Hadamard matrices of orders found here and, concurrently, by another group. Every matrix comes with a certificate from
which `verify/` rebuilds it and checks it in exact arithmetic.

A Hadamard matrix of order n is an n × n matrix H with entries ±1 and H Hᵀ = n I. It can exist only for n = 1, 2 or n ≡
0 (mod 4), and the Hadamard conjecture asserts that it exists for every multiple of 4. In August 2026 Alpöge, Voinov and
Reynolds-Haertle constructed the twelve unresolved orders below 2000 [1], and Schneider constructed 2060 [2]. This
repository constructs the next four, 2092, 2148, 2284 and 2292, and twenty-three more open orders up to 11804. Fourteen
of these follow from two published theorems that the tables did not combine: Paley's matrix of order q + 1 is an array
of sixteen circulants when q ≡ 3 (mod 8), and Seberry and Yamada multiply such arrays by T-matrices
([hadamard/PALEY_BLOCKS.md](hadamard/PALEY_BLOCKS.md)). Three more, 6108, 7532 and 10180, come from Miyamoto's matrix of
order 4q in sixteen commuting bordered blocks, built with the near-Williamson matrices of 2026 [11] and multiplied by
T-matrices ([hadamard/MIYAMOTO_BLOCKS.md](hadamard/MIYAMOTO_BLOCKS.md#near-williamson-matrices)). 2284 also follows from
a quaternary matrix of order 1142 posted concurrently by Zaitsev et al. [13].
The smallest order we know to be unresolved is now 2396 = 4 · 599 ([FRONTIER.md](FRONTIER.md)).

For skew-Hadamard matrices (H + Hᵀ = 2I) we prove that one of order 3(q + 1) exists for every prime power q ≡ 3
(mod 4): Paley's matrix is an array of sixteen negacyclic blocks, Frobenius pairs two of them, and T-matrices of order
3 then keep the first Goethals–Seidel block skew-type ([skew/PALEY_SKEW.md](skew/PALEY_SKEW.md)). This settles 52 values
v < 2500 that the Đoković–Seberry census lists as undecided, from 996 = 4 · 249 to 9996 = 4 · 2499. The block form and
the pairing are Yamada's (1991) and the multiplication is Seberry and Yamada's (1990); we found them combined nowhere,
so these orders were obtainable classically.

## Results

**Unresolved Hadamard orders.** For each order 4n below, n is listed as unresolved in Table 4 of Cati and
Pasechnik [3], and we found no later construction, except that 2284 follows from a quaternary matrix of order 1142
posted concurrently by Zaitsev et al. [13] (Turyn's doubling; [hadamard/2284](hadamard/2284/README.md#status)).
The fourteen orders built from Paley's matrix were within reach of two published theorems [14, 15], which the tables did
not combine; we found no earlier construction of them. 6108 and 10180 (n = 1527 and 2545) are first explicit
constructions: Đoković claimed them in 2010 through Seberry–Yamada's Corollary 8.8(i), whose printed proof does not give
Williamson-type matrices ([NOTES.md](NOTES.md) §7).

| Order | n | Construction |
|------:|---|---|
| [2092](hadamard/2092) | 523 | Goethals–Seidel array [6] over ℤ₅₂₃, multiplier subgroup of order 9 |
| [2148](hadamard/2148) | 537 = 3 · 179 | bordered Goethals–Seidel array over ℤ₃ × ℤ₁₇₈, border of size 12 |
| [2284](hadamard/2284) | 571 | Goethals–Seidel array over ℤ₅₇₁, three blocks related by a multiplier of order 3 |
| [2292](hadamard/2292) | 573 = 3 · 191 | Spence's construction over ℤ₂ × ℤ₂ × ℤ₁₄₃, three blocks related by an automorphism of order 15 |
| [3156](hadamard/3156) | 789 = 3 · 263 | Paley's matrix of order 1052 in sixteen circulant blocks [14], T-matrices of order 3 [15] |
| [3732](hadamard/3732) | 933 = 3 · 311 | Goethals–Seidel array over ℤ₃ × GF(311), multiplier subgroup of order 31; skew |
| [3964](hadamard/3964) | 991 | bordered Goethals–Seidel array over ℤ₂ × ℤ₂ × ℤ₂₄₇, border of size 12 |
| [4492](hadamard/4492) | 1123 | Goethals–Seidel array over ℤ₁₁₂₃, multiplier subgroup of order 33 |
| [5020](hadamard/5020) | 1255 = 5 · 251 | Goethals–Seidel array over ℤ₅ × GF(251), multiplier subgroup of order 50 |
| [5260](hadamard/5260) | 1315 = 5 · 263 | Paley's matrix of order 1052 in sixteen circulant blocks [14], T-matrices of order 5 [15] |
| [5932](hadamard/5932) | 1483 | Goethals–Seidel array over ℤ₁₄₈₃, multiplier subgroup of order 39 |
| [6036](hadamard/6036) | 1509 = 3 · 503 | Paley's matrix of order 2012 in sixteen circulant blocks [14], T-matrices of order 3 [15] |
| [6108](hadamard/6108) | 1527 = 3 · 509 | Miyamoto's matrix [8] in sixteen bordered blocks with near-Williamson matrices of order 127 [11], T-matrices of order 3 [15] |
| [6356](hadamard/6356) | 1589 = 7 · 227 | Paley's matrix of order 908 in sixteen circulant blocks [14], T-matrices of order 7 [15] |
| [6756](hadamard/6756) | 1689 = 3 · 563 | Paley's matrix of order 2252 in sixteen circulant blocks [14], T-matrices of order 3 [15] |
| [7364](hadamard/7364) | 1841 = 7 · 263 | Paley's matrix of order 1052 in sixteen circulant blocks [14], T-matrices of order 7 [15] |
| [7532](hadamard/7532) | 1883 = 7 · 269 | Miyamoto's matrix [8] in sixteen bordered blocks with near-Williamson matrices of order 67 [11], T-matrices of order 7 [15] |
| [7660](hadamard/7660) | 1915 = 5 · 383 | Paley's matrix of order 1532 in sixteen circulant blocks [14], T-matrices of order 5 [15] |
| [8620](hadamard/8620) | 2155 = 5 · 431 | Paley's matrix of order 1724 in sixteen circulant blocks [14], T-matrices of order 5 [15] |
| [9172](hadamard/9172) | 2293 | Miyamoto's theorem [8] from the matrix of order 2292 |
| [10060](hadamard/10060) | 2515 = 5 · 503 | Paley's matrix of order 2012 in sixteen circulant blocks [14], T-matrices of order 5 [15] |
| [10180](hadamard/10180) | 2545 = 5 · 509 | Miyamoto's matrix [8] in sixteen bordered blocks with near-Williamson matrices of order 127 [11], T-matrices of order 5 [15] |
| [10644](hadamard/10644) | 2661 = 3 · 887 | Paley's matrix of order 3548 in sixteen circulant blocks [14], T-matrices of order 3 [15] |
| [10932](hadamard/10932) | 2733 = 3 · 911 | Paley's matrix of order 3644 in sixteen circulant blocks [14], T-matrices of order 3 [15] |
| [11260](hadamard/11260) | 2815 = 5 · 563 | Paley's matrix of order 2252 in sixteen circulant blocks [14], T-matrices of order 5 [15] |
| [11572](hadamard/11572) | 2893 = 11 · 263 | Paley's matrix of order 1052 in sixteen circulant blocks [14], T-matrices of order 11 [15] |
| [11804](hadamard/11804) | 2951 = 13 · 227 | Paley's matrix of order 908 in sixteen circulant blocks [14], T-matrices of order 13 [15] |

**Other results.**

- **Hadamard matrices of 30 larger orders** from 12388 to 51580, above the range of Table 4, for which we found no
  earlier construction. They come from recursive theorems (Miyamoto [8], Yamada–Kiyasu [9], Cooper–Wallis [7]) applied
  to the matrices here and to objects published in 2026. For three of them the Handbook records the order as known, with
  no construction we could trace, and six (19028, 21044, 21748, 39260, 40244, 45044) were obtainable classically by the
  Paley-block route [14, 15]. Nine come from the closure audit of October 17 ([NOTES.md](NOTES.md) §28): Cooper–Wallis
  with the T-matrices of orders 151 and 163 of 2026 and Williamson matrices of orders 19 to 43 (12388, 13892, 14996,
  18724, 20212, 25972, 28036), with T-sequences of length 107 and the Williamson-type matrices of order 47 of 2026
  (20116), and Yamada–Kiyasu from Đoković and Kotsireas's skew matrix of order 956 (15292, obtainable classically). Also
  independent matrices of the known orders 2060, 3532 and 31316 (which SageMath's development version builds since
  October 3, 2026). See [hadamard](hadamard).
- **Hadamard matrices of 117 orders 4n with 3001 ≤ n ≤ 9999** (12052 to 39916) that Table 1.53 of the Handbook, with
  Đoković's updates of 2009–2010, lists as unknown (one of them, 32316, Đoković claimed only through a corollary whose
  printed proof fails). They come from the same Paley-block route [14, 15]. For 109 of them the T-matrices (orders 3 to
  73; for 28820 those of order 55, which we built by Yang's multiplication) were known before 2026, so these orders were
  obtainable classically; six need T-matrices of order 103, 151 or 163 from 2026; and two (16892, 25132) are built by
  SageMath's development version since October 3, 2026. See
  [hadamard/PALEY_BLOCKS.md](hadamard/PALEY_BLOCKS.md#beyond-table-4). Also 21036 = 4 · 3 · 1753, which Đoković claimed
  through Williamson-type matrices of order 1753 that rested on the same corollary: part (ii) of it with our symmetric
  matrix of order 876 gives them validly ([hadamard/21036](hadamard/21036)).
- **Hadamard matrices of 15 more orders 4n with 3279 ≤ n ≤ 9991** that Handbook Table 1.53 lists as unknown. Seberry and
  Yamada write Miyamoto's matrix of order 4q as an array of sixteen bordered blocks, but their multiplication theorems
  do not apply to it; when its Hadamard matrix of order q − 1 is of 2N-type (a negaperiodic Golay pair) and Paley's
  points are listed along a primitive root, the blocks commute and are transposed by one permutation, so T-matrices
  multiply it ([hadamard/MIYAMOTO_BLOCKS.md](hadamard/MIYAMOTO_BLOCKS.md)). Nine are obtainable classically, four need
  the T-matrices of 2026, and two (7179 and 7323) Đoković claimed only through the failing corollary below. With
  near-Williamson matrices of 2026 in the same bordered blocks, four more Handbook entries, 13988, 18292, 33356 and
  39812 (n = 3497, 4573, 8339 and 9953), are new beyond the tables, and 39308 (n = 9827) is a new matrix of an order
  that Đoković's set Δ already settles
  ([hadamard/MIYAMOTO_BLOCKS.md](hadamard/MIYAMOTO_BLOCKS.md#near-williamson-matrices)).
- **Orders recorded only through Seberry–Yamada's Corollary 8.8(i)**, whose printed proof does not give Williamson-type
  matrices ([NOTES.md](NOTES.md) §7): first valid constructions of 52 orders whose only traced construction used it. Ten
  are omitted by Table 4 (n = 519, 865, 1211, 1465, 1557, 1903, 2051, 2249, 2427 and 2611), among them 2427 and 2611
  (orders [9708](hadamard/9708) and [10444](hadamard/10444)), which [FRONTIER.md](FRONTIER.md) listed as unresolved; 26
  are Handbook entries from Đoković's set Δ, his 2010 claims or with t = 2, orders 12260 to 39860. 26 come from the
  Paley-block route, and 26 from Miyamoto's matrix in bordered blocks: two with a periodic Golay pair of length 404, two
  with Williamson matrices of order 153, six with near-Williamson matrices of order 93 or 127, and sixteen with
  Seberry–Yamada's symmetric Williamson-type matrices of orders 193 to 802
  ([hadamard/PALEY_BLOCKS.md](hadamard/PALEY_BLOCKS.md#orders-recorded-only-through-corollary-88i),
  [hadamard/MIYAMOTO_BLOCKS.md](hadamard/MIYAMOTO_BLOCKS.md#symmetric-williamson-type-matrices)). For 15 of the 52 the
  corollary can be repaired after all: Seberry–Yamada's symmetric inputs of orders 73, 193, 423, 433 and 729 have a
  repeated block, and Lemma 26 then gives Williamson-type matrices of orders 293, 773, 1693, 1733 and 2917 (checked), so
  those orders were obtainable classically; and the 173 · t orders have a published route since Verma's symmetric
  conference matrix of order 86 (arXiv:2610.08893, October 2026). Table 4 entries that Đoković settled through the
  corollary are in the table above (6108, 10180) or still open (2571).
- **Skew-Hadamard matrices of 101 orders** from 820 to 14352 for which we found no earlier construction: 32 from
  Goethals–Seidel difference families with a skew block, 8 from Seberry's theorem with new amicable Hadamard matrices or
  new orthogonal designs, 4 orders 16v (9424, 9776, 13072, 14384) from Wallis's product of a published skew matrix with
  our amicable pairs of orders 52, 76 and 124, and 57 from the theorem above, of order 3(q + 1), obtainable classically:
  52 census values 4v (one of them, 3516 = 4 · 879, recorded as known by two older tables without a construction) and
  five orders 8v and 16v below the smallest power of two in Cati–Pasechnik's Table 5. Also independent matrices of the
  known skew orders 1852, 3676 and 7656 (Wallis–Whiteman, which Table 5 missed). See [skew](skew).
- **Symmetric Hadamard matrices** of orders 532, 732, 876, 4828 and 5700, all undecided in the Đoković–Seberry
  census [4]. See [symmetric](symmetric).
- **Complex Hadamard matrices** of 11 orders 2N with N odd, blank in Seberry and Yamada's table of known orders [5]
  and in the papers since [10, 11, 12]: 214, 238, 302, 326, 366, 382, 426, 494, 518, 534 and 754, from searches over
  Szöllősi's array [10] and its bordered, character-bordered and twisted forms. They were found concurrently and
  independently by Zaitsev et al. [13], whose page has matrices of all 11; we claim no priority. 334 is the smallest
  open order now. Also explicit matrices of 76 known orders from 142 to 1966 with 4N − 1 a prime power, among them
  halves of Paley's matrix of order q + 1 ([complex-hadamard/PALEY.md](complex-hadamard/PALEY.md)); these follow from
  Ito's dicyclic Hadamard groups and Egan's Theorem 4.5 [12], which we had missed at first. First explicit matrices
  of 1618 and 1906, which the 1992 table records only through Seberry–Yamada's Corollary 8.8(i): its printed proof
  does not give Williamson-type matrices, but its matrices give a complex Hadamard matrix of order 2q from any
  quaternary Hadamard matrix of order (q − 1)/2, which we prove and use ([NOTES.md](NOTES.md) §19; Zaitsev's matrices of
  orders 538, 1018 and 1714 have this form, unstated). See [complex-hadamard](complex-hadamard).
- **Auxiliary objects:** Williamson-type matrices of order 111 ([near-williamson](near-williamson)); T-matrices of
  orders 97, 109 and 197, for which we found no earlier construction (97 and 109 may appear in Lou–Xia–Xia's work on
  order 6m + 1, which we could not read), and an independent instance of the known order 131 ([t-matrices](t-matrices));
  amicable Hadamard matrices of orders 52, 76, 100 and 124; and the orthogonal designs OD(100; 1, 5, 94), OD(156; 1, 3,
  152), OD(156; 1, 13, 142), OD(252; 1, 21, 230) and OD(364; 1, 13, 350)
  ([skew](skew#seberrys-theorem-with-new-ingredients)).

"New" means new to our knowledge after the literature search described in [NOVELTY.md](NOVELTY.md), which also
defines the status labels used in each README.

## Verification

```
python3 verify/verify.py                  # every real matrix and auxiliary object
python3 verify/verify.py 2092 skew/820    # selected orders
python3 verify/verify_complex.py          # the complex Hadamard matrices
```

`verify.py` rebuilds each matrix from its `certificate.json`, checks the defining identities of the ingredients
(complementary sequences, T-matrices, Williamson-type matrices, conference matrices, skew or symmetric inputs), checks
H Hᵀ = n I exactly, checks H + Hᵀ = 2I or H = Hᵀ where claimed, and compares the SHA-256 digest of the canonical CSV
encoding with the recorded one. It needs Python 3 and NumPy. The large orders dominate the running time: 51580 takes
about 40 minutes and 3.5 GB of memory. `verify_complex.py` works in exact integer arithmetic on the four sequences.

An independent checker in C reads a stored matrix and compares every pair of rows:

```
cc -O2 -fopenmp verify/check_matrix.c -o check_matrix
gzip -dc hadamard/2092/H2092.csv.gz | ./check_matrix
```

## Methods

| Construction | Orders | Description |
|---|---|---|
| Goethals–Seidel difference families invariant under a multiplier subgroup, found by an extended k-tree search | 2092, 3532, 4492, 5932 | [METHOD.md](METHOD.md) |
| The same with blocks permuted by a multiplier of order 3 | 2284 | [METHOD.md](METHOD.md), Block symmetries |
| Bordered Goethals–Seidel arrays: exact compression chains, and a split along ℤ₂ | 2148, 3964 | [METHOD.md](METHOD.md), [search/NOTES.md](search/NOTES.md) §18–19 |
| Cyclotomic products ℤ_m × GF(q) of index 4 to 10 | 3732, 5020; skew 820, 1044, 1204, 1276, 3380, 5668 | [search/NOTES.md](search/NOTES.md) §16, §20 |
| Spence's construction with blocks (A, B, B∘μ, B∘μ²) and a Legendre lift | 2292 | [NOTES.md](NOTES.md) §13, [search/spence_spin](search/spence_spin) |
| Skew difference families: compression chains, Paley products, low-energy k-trees | skew orders | [skew](skew), [search/NOTES.md](search/NOTES.md) §10–16, §20 |
| Paley's matrix in sixteen negacyclic blocks, Frobenius pairing, T-matrices of order 3 | skew 3(q + 1): 996 to 14352 | [skew/PALEY_SKEW.md](skew/PALEY_SKEW.md), [NOTES.md](NOTES.md) §27 |
| Propus quadruples and the symmetric m(m − 1), m(m + 3) constructions | symmetric orders | [symmetric](symmetric), [search/NOTES.md](search/NOTES.md) §17, §27–31 |
| Szöllősi's array, its bordered, character-bordered and twisted forms | complex orders | [complex-hadamard](complex-hadamard) |
| Cooper–Wallis with T-matrices, near-Williamson and bordered Williamson-type matrices | 12596 to 39260 | [williamson-type](williamson-type), [t-matrices](t-matrices) |

[NOTES.md](NOTES.md) and [search/NOTES.md](search/NOTES.md) also record the structures that turned out to be empty,
the obstructions found on the way, and calibrations of the solution-count estimates. [FRONTIER.md](FRONTIER.md) lists
the unresolved orders up to 4 · 2999 with the cheapest structure known for each.

## Layout

| Path | Contents |
|---|---|
| `hadamard/<order>/`, `skew/<order>/`, `symmetric/<order>/` | `certificate.json`, `README.md`, and where stored `sequences.txt` and the matrix `H<order>.csv.gz` (one row per line, entries separated by commas) |
| `complex-hadamard/` | certificates `C<order>.json` and the constructions |
| `t-matrices/`, `near-williamson/`, `williamson-type/` | auxiliary objects used by the Cooper–Wallis orders |
| `verify/` | the verifiers |
| `search/` | the search programs, with [search/README.md](search/README.md) as an index |

## References

1. L. Alpöge, P. Voinov, S. Reynolds-Haertle, Hadamard matrices of orders 668, 716, 892, 1132, 1244, 1388, 1436,
   1676, 1772, 1916, 1948 and 1964, announcement, August 12, 2026.
2. Schneider, Hadamard matrix of order 2060, public gist, August 23, 2026.
3. M. Cati, D. V. Pasechnik, A database of constructions of Hadamard matrices, arXiv:2411.18897v2 (2025).
4. D. Ž. Đoković, J. Seberry, census of skew and symmetric Hadamard matrices of order 4v for odd v < 2500,
   May 28, 2025.
5. J. Seberry, M. Yamada, Hadamard matrices, sequences, and block designs, in J. H. Dinitz, D. R. Stinson (eds.),
   Contemporary Design Theory: A Collection of Surveys, Wiley (1992) 431–560.
6. J.-M. Goethals, J. J. Seidel, Orthogonal matrices with zero diagonal, Canad. J. Math. 19 (1967) 1001–1010.
7. J. Cooper, J. Wallis, A construction for Hadamard arrays, Bull. Austral. Math. Soc. 7 (1972) 269–278.
8. M. Miyamoto, A construction of Hadamard matrices, J. Combin. Theory Ser. A 57 (1991) 86–108.
9. M. Yamada, Some new series of Hadamard matrices, J. Austral. Math. Soc. Ser. A 46 (1989) 371–383.
10. Szöllősi, A complex Hadamard matrix of order 94, arXiv:2603.09572 (2026).
11. H. Kharaghani, A. Mohammadian, B. Tayfeh-Rezaie, V. Zaitsev, arXiv:2605.08661v2 (2026).
12. R. Egan, Phased unitary Golay pairs, Butson Hadamard matrices and a conjecture of Ito's, Des. Codes Cryptogr. 87
    (2019) 67–74.
13. V. Zaitsev (with H. Kharaghani, A. Mohammadian, B. Tayfeh-Rezaie), Quaternary Hadamard matrices,
    https://www.cs.uleth.ca/~zaitsevv/quaternary.html, last updated October 14, 2026.
14. B. Schmidt, Williamson matrices and a conjecture of Ito's, Des. Codes Cryptogr. 17 (1999) 61–68, Theorems 2.2
    and 3.3 (the latter due to Ito, 1994); M. Yamada, Hadamard matrices of generalized quaternion type, Discrete Math.
    87 (1991) 187–196.
15. J. Seberry, M. Yamada, On the products of Hadamard matrices, Williamson matrices and other orthogonal matrices
    using M-structures, J. Combin. Math. Combin. Comput. 7 (1990) 97–137, Theorem 4.

## Citation

If you use these matrices or programs, please cite the repository ([CITATION.cff](CITATION.cff)):

> Kishan, *New Hadamard matrices*, https://github.com/vatskishan03/Hadamard, 2026.

## License

Apache License 2.0; see [LICENSE](LICENSE).
