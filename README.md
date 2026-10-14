# Hadamard

Explicit Hadamard matrices for orders that had no published construction, with
certificates that can be checked independently.

A Hadamard matrix of order n is an n × n matrix H with entries ±1 and HHᵀ = nI.
Such a matrix can exist only if n = 1, 2 or n ≡ 0 (mod 4), and the Hadamard
conjecture states that it exists for every multiple of 4. After the constructions
announced in August 2026 for the twelve unresolved orders below 2000 [1], the
smallest unresolved orders were 2060 and 2092. This repository constructs 2092 and 2148; with 2060 (Schneider,
August 2026), the smallest order we know to be unresolved is now 2292 = 4 · 573.

## Results

| Order | n | Construction | Status |
|------:|---|---|---|
| [2092](2092) | 4 · 523 | Goethals–Seidel array, supplementary difference set in ℤ₅₂₃ invariant under a multiplier subgroup of order 9 | New |
| [2148](2148) | 4 · 537 | Bordered Goethals–Seidel array over ℤ₅₃₄ ≅ ℤ₆ × GF(89) with a border of size 12 (G/K ≅ ℤ₃), sequences invariant under a multiplier subgroup of order 11; found by splitting the conditions along ℤ₂ ([attack/NOTES.md](attack/NOTES.md) §19) | New |
| [2284](2284) | 4 · 571 | Goethals–Seidel, ℤ₅₇₁, multiplier subgroup of order 5; three blocks related by a multiplier of order 3 | New |
| [3732](3732) | 4 · 933 | Goethals–Seidel array, skew; ℤ₉₃₃ ≅ ℤ₃ × GF(311), sequences invariant under the multiplier subgroup of order 31 (h ≡ 1 mod 3, h a 10th power mod 311): a cyclotomic product of index 10 | New; also gives 7464 by doubling |
| [3964](3964) | 4 · 991 | Bordered Goethals–Seidel array over ℤ₂ × ℤ₂ × ℤ₂₄₇ with a border of size 12; sequences invariant under a group of order 27 | New |
| [4492](4492) | 4 · 1123 | Goethals–Seidel, ℤ₁₁₂₃, multiplier subgroup of order 33 | New |
| [5020](5020) | 4 · 1255 | Goethals–Seidel; ℤ₁₂₅₅ ≅ ℤ₅ × GF(251), multiplier subgroup of order 50 (h ≡ ±1 mod 5, h a 10th power mod 251) | New |
| [5932](5932) | 4 · 1483 | Goethals–Seidel, ℤ₁₄₈₃, multiplier subgroup of order 39 | New |
| [12596](12596) | 4 · 3149 | Cooper–Wallis: T-matrices of order 67 and bordered Williamson-type matrices of order 47 ([williamson-type](williamson-type)) | Probably new; above the range of [3] |
| [14708](14708) | 4 · 3677 | Miyamoto's theorem from a skew-Hadamard matrix of order 3676 (Goethals–Seidel over ℤ₉₁₉, multiplier subgroup of order 17, spin structure) | Probably new; above the range of [3] |
| [14932](14932) | 4 · 3733 | Miyamoto's theorem (q = 3733) from the matrix of order 3732 | First explicit construction we know of; Handbook Table 1.53 lists 3733 with t = 2, with no construction we could trace |
| [19364](19364) | 4 · 4841 | Cooper–Wallis: T-matrices of order 103 and bordered Williamson-type matrices of order 47 | Probably new; above the range of [3] |
| [20084](20084) | 4 · 5021 | Miyamoto's theorem (q = 5021) from the matrix of order 5020 | First explicit construction we know of; Handbook Table 1.53 lists 5021 with t = 2, with no construction we could trace |
| [21748](21748) | 4 · 5437 | Miyamoto's theorem from a Cooper–Wallis matrix of order 5436 (T-matrices of order 151, Williamson order 9) | Probably new; above the range of [3] |
| [24124](24124) | 4 · 6031 | Yamada's theorem (q = 6029) from a skew-Hadamard matrix of order 3016, the double of a skew Goethals–Seidel matrix of order 1508 over ℤ₃₇₇ | Probably new; above the range of [3] |
| [24308](24308) | 4 · 6077 | Cooper–Wallis: T-matrices of order 103 and Williamson-type matrices of order 59 from near-Williamson matrices (Kharaghani, Mohammadian, Tayfeh-Rezaie, Zaitsev) | Probably new; above the range of [3] |
| [24412](24412) | 4 · 6103 | Yamada's theorem (Kiyasu), q = 6101, from this repository's [skew 3052](skew/3052) (no doubling) | Probably new; Handbook Table 1.53 lists 6103 with t = 3 |
| [27604](27604) | 4 · 6901 | Cooper–Wallis: T-matrices of order 103 and Williamson-type matrices of order 67 from near-Williamson matrices (Kharaghani, Mohammadian, Tayfeh-Rezaie, Zaitsev) | Probably new; above the range of [3] |
| [29404](29404) | 4 · 7351 | Yamada's theorem (q = 7349) from this repository's [skew 3676](skew/3676) (an order first published by Schumer in September 2026) | First explicit construction we know of; Handbook Table 1.53 lists 7351 with t = 2, with no construction we could trace |
| [31316](31316) | 4 · 7829 | Miyamoto's theorem from a Cooper–Wallis matrix of order 7828 (T-matrices of order 103, Williamson order 19) | Probably new; above the range of [3] |
| [31412](31412) | 4 · 7853 | Miyamoto's theorem from a Cooper–Wallis matrix of order 7852 (T-matrices of order 151, Williamson order 13) | Probably new; above the range of [3] |
| [34196](34196) | 4 · 8549 | Cooper–Wallis: T-matrices of order 103 and bordered Williamson-type matrices of order 83 (four sequences on ℤ₈₂, [williamson-type](williamson-type)) | Probably new; above the range of [3] |
| [35636](35636) | 4 · 8909 | Cooper–Wallis: T-matrices of order 151 and Williamson-type matrices of order 59 from near-Williamson matrices (Kharaghani, Mohammadian, Tayfeh-Rezaie, Zaitsev) | Probably new; above the range of [3] |
| [38468](38468) | 4 · 9617 | Cooper–Wallis: T-matrices of order 163 and Williamson-type matrices of order 59 from near-Williamson matrices (Kharaghani, Mohammadian, Tayfeh-Rezaie, Zaitsev) | Probably new; above the range of [3] |
| [39260](39260) | 4 · 9815 | Cooper–Wallis: T-matrices of order 151 and Williamson-type matrices of order 65 from near-Williamson matrices (Kharaghani, Mohammadian, Tayfeh-Rezaie, Zaitsev) | Probably new; above the range of [3] |
| [51580](51580) | 4 · 12895 | Yamada's theorem (Kiyasu), q = 12893, from the skew-Hadamard matrix of order 6448 = 16 · 403: this repository's [skew 1612](skew/1612), doubled twice | Probably new; above the range of [3] |
| [2060](2060) | 4 · 515 | Cooper–Wallis: T-matrices of order 103 and Williamson matrices of order 5 | Independent construction; the order was first settled publicly on 2026-08-23 [2] |
| [3532](3532) | 4 · 883 | Goethals–Seidel, ℤ₈₈₃, multiplier subgroup of order 21 | Order known [4]; new explicit matrix |

An order is marked new if it is listed as unresolved in the table of Cati and
Pasechnik [3], is not covered by the later updates of Đoković [4, 5] or by the
recursive constructions of Paley, Miyamoto, Yamada, Mathon and Spence, and no
construction was found in a literature search in October 2026. The table in [3]
is not complete (nine of its entries were settled in [4]), so these claims carry
that caveat.

The eighteen orders above 4 · 2999 lie beyond [3]. For the seventeen with n < 10000 we used Table 1.53 of the
Handbook of Combinatorial Designs (2007) together with the updates in [4, 5]. For 14932 (n = 3733), 20084
(n = 5021) and 29404 (n = 7351) that table claims t = 2, but we could not trace a construction behind any of these entries; see their READMEs. For 51580 (n = 12895), which is
beyond that table too, we checked Paley, the Williamson-type families of Seberry–Yamada's 1992 survey and the
2026 deposits. The README of each order gives the details. Their matrices are not stored (up to 7 GB as CSV);
`verify/verify.py` rebuilds them from the certificates.

Four of these orders come from this repository's skew-Hadamard matrices through Yamada's theorem:
- 24124 from skew 1508, doubled once;
- 24412 from skew 3052;
- 29404 from skew 3676;
- 51580 from skew 1612, doubled twice.

14932 and 20084 come from this repository's matrices of orders 3732 and 5020 through Miyamoto's theorem.

24124 is also reachable by Cooper–Wallis from T-matrices of order 163 and Williamson matrices of order 37;
see its README.

Orders under investigation: [4412](4412). [2148](2148), listed here before as the smallest order we knew to be
unresolved, is now settled above. [3676](3676) was settled by Schumer [8] while our search for it was running.

## Skew-Hadamard results

A Hadamard matrix H is skew if H + Hᵀ = 2I. Skew-Hadamard matrices are conjectured to exist for every
order divisible by 4. Each matrix below comes from a Goethals–Seidel difference family with a skew
first block. The group is ℤᵥ, except for 1100 (ℤ₅ × ℤ₅ × ℤ₁₁), 1372 (GF(343) = ℤ₇³),
2052 (ℤ₁₉ × GF(27) = ℤ₁₉ × ℤ₃³) and 3380 (ℤ₅ × GF(169) = ℤ₅ × ℤ₁₃²). All four blocks are invariant under a group H of automorphisms with −1 ∉ H.

**New orders** (to our knowledge; see the checklist below):

| Order | v | H | Row sums | Table 5 entry |
|------:|---|---|---|---|
| [820](skew/820) | 205 (5 · 41) | ⟨51⟩, order 5 (cyclotomic product of index 8) | 1, 1, −17, −23 | 205(3) |
| [836](skew/836) | 209 (11 · 19) | ⟨26⟩, order 15 | 1, 9, 15, 23 | 209(4) |
| [1012](skew/1012) | 253 (11 · 23) | ⟨12⟩, order 11 (Paley product) | 1, 5, 19, 25 | 253(4) |
| [1036](skew/1036) | 259 (7 · 37) | ⟨9⟩, order 9 | 1, 5, 7, 31 | 259(4) |
| [1044](skew/1044) | 261 (3² · 29) | ⟨190⟩, order 7 (cyclotomic product of index 4) | 1, 9, −1, 31 | 261(3) |
| [1100](skew/1100) | 275 (5² · 11), group ℤ₅ × ℤ₅ × ℤ₁₁ | Singer-3 on ℤ₅² × ⟨3⟩ on ℤ₁₁, order 15 | 1, 1, 3, 33 | 275(4) |
| [1140](skew/1140) | 285 (3 · 5 · 19) | ⟨61⟩, order 9 (Paley product) | 1, 17, 15, 25 | 285(3); first recorded, but Seberry's 1978 theorem gives it from an OD(60; 1, 3, 56), see its README |
| [1148](skew/1148) | 287 (7 · 41) | ⟨16⟩, order 15 | 1, 3, 7, 33 | 287(4) |
| [1204](skew/1204) | 301 (7 · 43) | ⟨127⟩, order 7 (cyclotomic product of index 6) | 1, −11, −11, −31 | 301(3) |
| [1236](skew/1236) | 309 (3 · 103) | ⟨13⟩, order 17 | 1, 1, 3, 35 | 309(3) |
| [1276](skew/1276) | 319 (11 · 29) | ⟨45⟩, order 7 (cyclotomic product of index 4) | 1, 25, 19, 17 | 319(3) |
| [1364](skew/1364) | 341 (11 · 31) | ⟨14⟩, order 15 | 1, 9, 21, 29 | 341(4) |
| [1372](skew/1372) | 343 (7³), group GF(343) = ℤ₇³ | cyclotomic subgroup of order 19 of GF(343)* | 1, 1, 1, 37 | 343(4) |
| [1380](skew/1380) | 345 (3 · 5 · 23) | ⟨301⟩, order 11 (Paley product) | 1, 19, 17, 27 | 345(4) |
| [1412](skew/1412) | 353 (prime) | ⟨22⟩, order 11 | 1, 21, 21, 23 | 353(4) |
| [1508](skew/1508) | 377 (13 · 29) | ⟨16⟩, order 21 | 1, 7, 27, 27 | 377(5) |
| [1564](skew/1564) | 391 (17 · 23) | ⟨324⟩, order 11 (Paley product) | 1, 19, 19, 29 | 391(4) |
| [1612](skew/1612) | 403 (13 · 31) | ⟨100⟩, order 15 | 1, 5, 25, 31 | 403(5) |
| [1804](skew/1804) | 451 (11 · 41) | ⟨16, 42⟩, order 25 | 1, 1, 11, 41 | none (implicit 2, unexplained); 451(3) in Seberry–Yamada 1992 and Koukouvinos–Stylianou 2008; census: undecided |
| [1860](skew/1860) | 465 (3 · 5 · 31) | ⟨226⟩, order 15 (Paley product) | 1, 25, 35, 3 | 465(3) |
| [1876](skew/1876) | 469 (7 · 67) | ⟨9⟩, order 33 | 1, 1, 5, 43 | 469(3) |
| [2052](skew/2052) | 513 (3³ · 19), group ℤ₁₉ × GF(27) | squares of GF(27)*, order 13 (Paley product) | 1, 27, 19, 31 | 513(4) |
| [2404](skew/2404) | 601 (prime) | ⟨2⟩, order 25 | 1, 1, 1, 49 | 601(5) |
| [2452](skew/2452) | 613 (prime) | ⟨37⟩, order 17 | 1, 1, 35, 35 | 613(3) |
| [2812](skew/2812) | 703 (19 · 37) | ⟨16, 26⟩, order 27 | 1, 1, 1, 53 | 703(3) |
| [3052](skew/3052) | 763 (7 · 109) | ⟨9⟩, order 27 | 1, 1, 5, 55 | 763(4) |
| [3380](skew/3380) | 845 (5 · 13²), group ℤ₅ × GF(169) | ⟨w⁸⟩ on GF(169), order 21 (cyclotomic product of index 8) | 1, 43, −3, −39 | 845(4); census: undecided; smallest known before was 16 · 845 = 13520 |
| [3556](skew/3556) | 889 (7 · 127) | ⟨2,25⟩, order 63 | 1, 5, 41, 43 | 889(4) |
| [3732](3732) | 933 (3 · 311) | ⟨7⟩, order 31 (cyclotomic product of index 10) | 1, 61, 3, 1 | 933(10); census: undecided; also an unresolved Hadamard order |
| [3612](skew/3612) | 903 (3 · 7 · 43) | ⟨568⟩, order 21 (Paley product) | 1, 43, 41, 9 | 903(4) |
| [5668](skew/5668) | 1417 (13 · 109) | ⟨1405, 328⟩, order 81 (cyclotomic product of index 4 with ⟨3⟩ on ℤ₁₃) | 1, 7, −53, −53 | beyond Table 5 (n ≤ 999); census: undecided; smallest known before was 32 · 1417 = 45344 |
| [6004](skew/6004) | 1501 (19 · 79) | ⟨799⟩, order 39 (Paley product) | 1, 77, 7, −5 | beyond Table 5 (n ≤ 999); census: undecided; smallest known before was 8 · 1501 = 12008 |

**Known orders, independent matrices.** These orders were claimed as new at first. Both claims are
withdrawn:

| Order | v | H | Row sums | Prior construction |
|------:|---|---|---|---|
| [1852](skew/1852) | 463 (prime) | ⟨33⟩, order 21 | 1, 1, 1, 43 | Đoković, Australas. J. Combin. 10 (1994) 259–264 |
| [3676](skew/3676) | 919 (prime) | ⟨58⟩, order 17; spin μ = 52 | 1, 35, 35, 35 | Schumer, Zenodo 2026-09-17, doi:10.5281/zenodo.22820997 (same spin structure) |

Doubling (Seberry) then gives skew matrices of orders 2ᵏ · 4v for every k ≥ 0.

**Novelty checklist.** An order is called new only if all of the following hold. Each order's README
records the result of each check.

1. **Cati–Pasechnik Table 5** (arXiv:2411.18897v2) lists v with m ≥ 3. That table is verified by
   implemented constructions only up to order 1208. Beyond that it copies the older Seberry–Yamada
   Table 9.2 and can be outdated: 463 is an example.
2. **Đoković's papers** (1992–2024, alone and with Kotsireas and Golubitsky), Koukouvinos–Stylianou
   (Discrete Math. 2008, Tables 1 and 3), and Đoković–Seberry's census of skew orders 4v for v < 2500
   do not list the order.
3. **2026 public deposits** (Schumer; Alpöge et al.; JD-Jones-ASES; Schneider; Karoui; Piza) do not
   contain the order. SageMath's `unknown_skew_hadamard` list is consulted as well.
4. **No infinite family reaches it.** The families checked are Paley I, Szekeres 2(q + 1),
   Wallis–Whiteman, Spence, Momihara–Xiang, the product theorem n(pᵏ + 1), and amicable products.

The searches use `attack/skew_triage.py`, `attack/skew_rank.py`, `attack/skew_ktree.py` with `attack/match4.c`,
`attack/skew_exhaustive.py`, `attack/skew_group.py`, `attack/skew_field.py`, the exact compression-chain search
`attack/chain.py` (driven over the census by `attack/census_screen.py` and `attack/screen_chain.py`, which found
2812 and 3052), the unequal-list low-energy search `attack/skew_lowE.py` (which found 2452 and 1012), the Paley-product configuration search `attack/paley_config.py` (§15 of `attack/NOTES.md`), the cyclotomic-product search `attack/cyclo_product.py` with `attack/join3.c` (driven over the census by `attack/skew_kmult_sweep.py`, which found 5668), the low-energy four-list k-tree `attack/cp_lowe.py` with `attack/cp_enum.c` and `attack/kt4.c` (which found 820, 1044, 1204, 1276 and 3380), and, for the spin structure, `attack/spin_skew.py` with `attack/spin_ls.c`. Each directory's README
records the command, code commit, seed and caps of its search. `verify/verify.py` checks H + Hᵀ = 2I for
these certificates.

## Symmetric Hadamard results

A Hadamard matrix H is symmetric if H = Hᵀ. The Đoković–Seberry census lists 595 odd v < 2500 for which no
symmetric Hadamard matrix of order 4v was known. The matrices below come from propus-type quadruples (A, B, B, D) of
circulants with A symmetric and A Aᵀ + 2 B Bᵀ + D Dᵀ = 4v I, invariant under a multiplier group H
(`attack/propus.py`, `attack/NOTES.md` §17; 876 through compressions, `attack/propus_chain.py`, §25 and §28).

| Order | v | H | Row sums (A, B, D) | Census |
|------:|---|---|---|---|
| [532](symmetric/532) | 133 (7 · 19) | ⟨26⟩, order 6 | 11, 5, 19 | undecided |
| [876](symmetric/876) | 219 (3 · 73) | ⟨4⟩, order 9 (cyclotomic product ℤ₃ × GF(73) of index 8) | 3, 17, 17 | undecided |

## Verification

```
python3 verify/verify.py
```

rebuilds every matrix from its `certificate.json`, checks the defining identity,
checks HHᵀ = nI exactly and compares the SHA-256 digest of the canonical CSV
encoding with the recorded value. It needs NumPy and takes about five minutes.

An independent checker in C verifies the shipped matrices directly:

```
cc -O2 -fopenmp verify/check_matrix.c -o check_matrix
gzip -dc 4492/H4492.csv.gz | ./check_matrix
```

## Layout

Each order has its own directory containing

- `certificate.json`: the construction data from which the matrix is rebuilt,
- `sequences.txt`: for Goethals–Seidel orders, the first rows of the four circulant blocks,
- `H<n>.csv.gz`: the matrix, one row per line, entries separated by commas,
- `README.md`: construction details and status.

Other top-level contents:

- [`skew`](skew): skew-Hadamard matrices, one directory per order, in the same layout.
- [`symmetric`](symmetric): symmetric Hadamard matrices, one directory per order, in the same layout.
- [`t-matrices`](t-matrices): the T-matrices of orders 97 to 197 found during this work.
- [`williamson-type`](williamson-type): bordered Williamson-type matrices of order w from four sequences on ℤ_(w−1)
  (construction, proof, calibration and limits).
- [`near-williamson`](near-williamson): near-Williamson matrices of Kharaghani, Mohammadian, Tayfeh-Rezaie and Zaitsev
  (2026), used as Williamson-type matrices for the orders 24308 to 39260.
- [`search`](search): the programs that found the Goethals–Seidel orders.
- [`METHOD.md`](METHOD.md): how the search works.
- [`FRONTIER.md`](FRONTIER.md): the unresolved orders up to 4 · 2999, with the cheapest structure for each and its predicted cost.
- [`NOTES.md`](NOTES.md): what the October 2026 searches showed, including the structures that turned out empty.

## References

1. L. Alpöge, P. Voinov, S. Reynolds-Haertle, Hadamard matrices of orders 668, 716,
   892, 1132, 1244, 1388, 1436, 1676, 1772, 1916, 1948 and 1964, announcement,
   August 12, 2026.
2. Hadamard matrix of order 2060, public gist by Schneider, August 23, 2026; decoded
   in JD-Jones-ASES/Hadamard-T.
3. M. Cati, D. V. Pasechnik, A database of constructions of Hadamard matrices,
   arXiv:2411.18897v2 (2025), Table 4.
4. D. Ž. Đoković, Small orders of Hadamard matrices and base sequences,
   Int. Math. Forum 6 (2011) 3061–3067, arXiv:1008.2043.
5. D. Ž. Đoković, Hadamard matrices of small order and Yang conjecture,
   J. Combin. Des. 18 (2010) 254–259, arXiv:0912.5091.
6. J.-M. Goethals, J. J. Seidel, Orthogonal matrices with zero diagonal,
   Canad. J. Math. 19 (1967) 1001–1010.
7. J. Cooper, J. Wallis, A construction for Hadamard arrays,
   Bull. Austral. Math. Soc. 7 (1972) 269–278.
8. Schumer, Hadamard matrices of orders 3676, 5436, 6644 and 7852, Zenodo,
   September 17, 2026, doi:10.5281/zenodo.22820997.
