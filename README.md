# New Hadamard matrices

Explicit Hadamard, skew-Hadamard and symmetric Hadamard matrices of orders for which no construction was known, and
complex Hadamard matrices found concurrently with another group. Every matrix has a certificate from which `verify/`
rebuilds it and checks it in exact arithmetic.

A Hadamard matrix of order n is an n × n matrix H with entries ±1 and H Hᵀ = n I. It can exist only for n = 1, 2 or
n ≡ 0 (mod 4), and the Hadamard conjecture asserts that it exists for every multiple of 4. In August 2026 Alpöge,
Voinov and Reynolds-Haertle constructed the twelve unresolved orders below 2000 [1], and Schneider constructed 2060 [2].
This repository constructs the next three, 2092, 2148 and 2284, and 24 more of the orders that Table 4 of Cati and
Pasechnik [3] lists as unresolved, up to 11804. Our matrix of order 2292 turned out to have been preceded by an
independent one [16]. The smallest unresolved order is now 2396 = 4 · 599 ([FRONTIER.md](FRONTIER.md)).

Most of the larger orders come from two block decompositions. Paley's matrix of order q + 1, q ≡ 3 (mod 8), is an
array of sixteen circulant blocks [14], and Seberry and Yamada multiply such arrays by T-matrices [15]; we found the two
results combined nowhere ([hadamard/PALEY_BLOCKS.md](hadamard/PALEY_BLOCKS.md)). Miyamoto's matrix of order 4q [8] is
an array of sixteen bordered blocks. When its input is a negaperiodic or periodic Golay pair, or symmetric
Williamson-type matrices, and Paley's points are listed along a primitive root, these blocks commute and one
permutation transposes them all, so T-matrices multiply this array as well; we found this stated nowhere
([hadamard/MIYAMOTO_BLOCKS.md](hadamard/MIYAMOTO_BLOCKS.md)). For skew-Hadamard matrices the blocks of Paley's matrix
give a theorem: a skew-Hadamard matrix of order 3(q + 1) exists for every prime power q ≡ 3 (mod 4)
([skew/PALEY_SKEW.md](skew/PALEY_SKEW.md)).

## Results

| Family | Orders | Status |
|---|---|---|
| Hadamard, unresolved in Table 4 [3] | [27](#unresolved-orders) | 10 New, 14 New but obtainable classically, 2 First explicit construction, 1 Concurrent |
| Hadamard, omitted by Table 4 through Seberry–Yamada's Corollary 8.8(i) | [10](hadamard#orders-recorded-only-through-corollary-88i) | First explicit construction |
| Hadamard, beyond Table 4 (12052 to 51580) | [207](hadamard#orders-beyond-the-tables) | 34 New, 124 New but obtainable classically, 49 First explicit construction |
| Skew-Hadamard (820 to 14384) | [101](skew) | 37 New, 62 New but obtainable classically, 2 First explicit construction |
| Symmetric Hadamard (356 to 9892) | [53](symmetric) | 5 New, 48 New but obtainable classically |
| Complex Hadamard (214 to 1906) | [14](complex-hadamard) | 12 Concurrent, 2 First explicit construction |
| T-matrices | [97, 109, 197](t-matrices) | New |

The labels are defined in [NOVELTY.md](NOVELTY.md), which also lists the literature searched. "Obtainable
classically" means that every ingredient was published before 2026 and only their combination is new. The repository
also holds independent matrices of known orders (7 Hadamard, 3 skew, 2 symmetric and 76 complex), the first explicit
amicable Hadamard matrices of orders 52, 76, 100 and 124 and new orthogonal designs
([skew](skew#seberrys-theorem-with-new-ingredients)), and Williamson-type matrices of order 111
([near-williamson](near-williamson)).

### Unresolved orders

For each order 4n, n is listed as unresolved in Table 4 of Cati and Pasechnik [3].

| Order | n | Construction | Status |
|------:|---|---|---|
| [2092](hadamard/2092) | 523 | Goethals–Seidel array [6] over ℤ₅₂₃, multiplier subgroup of order 9 | New |
| [2148](hadamard/2148) | 537 = 3 · 179 | bordered Goethals–Seidel array over ℤ₃ × ℤ₁₇₈, border of size 12 | New |
| [2284](hadamard/2284) | 571 | Goethals–Seidel array over ℤ₅₇₁, three blocks related by a multiplier of order 3 | Concurrent [13] |
| [2636](hadamard/2636) | 659 | bordered Goethals–Seidel array over ℤ₄ × ℤ₄ × ℤ₄₁, border of size 12 | New |
| [3156](hadamard/3156) | 789 = 3 · 263 | Paley's matrix of order 1052 in sixteen circulant blocks [14], T-matrices of order 3 [15] | New, but obtainable classically |
| [3732](hadamard/3732) | 933 = 3 · 311 | Goethals–Seidel array over ℤ₃ × GF(311), multiplier subgroup of order 31; skew | New |
| [3964](hadamard/3964) | 991 | bordered Goethals–Seidel array over ℤ₂ × ℤ₂ × ℤ₂₄₇, border of size 12 | New |
| [4492](hadamard/4492) | 1123 | Goethals–Seidel array over ℤ₁₁₂₃, multiplier subgroup of order 33 | New |
| [5020](hadamard/5020) | 1255 = 5 · 251 | Goethals–Seidel array over ℤ₅ × GF(251), multiplier subgroup of order 50 | New |
| [5260](hadamard/5260) | 1315 = 5 · 263 | Paley's matrix of order 1052 in sixteen circulant blocks, T-matrices of order 5 | New, but obtainable classically |
| [5932](hadamard/5932) | 1483 | Goethals–Seidel array over ℤ₁₄₈₃, multiplier subgroup of order 39 | New |
| [6036](hadamard/6036) | 1509 = 3 · 503 | Paley's matrix of order 2012 in sixteen circulant blocks, T-matrices of order 3 | New, but obtainable classically |
| [6108](hadamard/6108) | 1527 = 3 · 509 | Miyamoto's matrix [8] in sixteen bordered blocks with near-Williamson matrices of order 127 [11], T-matrices of order 3 | First explicit construction |
| [6356](hadamard/6356) | 1589 = 7 · 227 | Paley's matrix of order 908 in sixteen circulant blocks, T-matrices of order 7 | New, but obtainable classically |
| [6756](hadamard/6756) | 1689 = 3 · 563 | Paley's matrix of order 2252 in sixteen circulant blocks, T-matrices of order 3 | New, but obtainable classically |
| [7364](hadamard/7364) | 1841 = 7 · 263 | Paley's matrix of order 1052 in sixteen circulant blocks, T-matrices of order 7 | New, but obtainable classically |
| [7532](hadamard/7532) | 1883 = 7 · 269 | Miyamoto's matrix in sixteen bordered blocks with near-Williamson matrices of order 67, T-matrices of order 7 | New |
| [7660](hadamard/7660) | 1915 = 5 · 383 | Paley's matrix of order 1532 in sixteen circulant blocks, T-matrices of order 5 | New, but obtainable classically |
| [8620](hadamard/8620) | 2155 = 5 · 431 | Paley's matrix of order 1724 in sixteen circulant blocks, T-matrices of order 5 | New, but obtainable classically |
| [9172](hadamard/9172) | 2293 | Miyamoto's theorem [8] from the matrix of order 2292 | New |
| [10060](hadamard/10060) | 2515 = 5 · 503 | Paley's matrix of order 2012 in sixteen circulant blocks, T-matrices of order 5 | New, but obtainable classically |
| [10180](hadamard/10180) | 2545 = 5 · 509 | Miyamoto's matrix in sixteen bordered blocks with near-Williamson matrices of order 127, T-matrices of order 5 | First explicit construction |
| [10644](hadamard/10644) | 2661 = 3 · 887 | Paley's matrix of order 3548 in sixteen circulant blocks, T-matrices of order 3 | New, but obtainable classically |
| [10932](hadamard/10932) | 2733 = 3 · 911 | Paley's matrix of order 3644 in sixteen circulant blocks, T-matrices of order 3 | New, but obtainable classically |
| [11260](hadamard/11260) | 2815 = 5 · 563 | Paley's matrix of order 2252 in sixteen circulant blocks, T-matrices of order 5 | New, but obtainable classically |
| [11572](hadamard/11572) | 2893 = 11 · 263 | Paley's matrix of order 1052 in sixteen circulant blocks, T-matrices of order 11 | New, but obtainable classically |
| [11804](hadamard/11804) | 2951 = 13 · 227 | Paley's matrix of order 908 in sixteen circulant blocks, T-matrices of order 13 | New, but obtainable classically |

2284 also follows from the quaternary Hadamard matrix of order 1142 that Zaitsev et al. posted concurrently [13], by
Turyn's doubling. Our matrix of order 2292 = 4 · 573 (Spence's construction over ℤ₂ × ℤ₂ × ℤ₁₄₃) is not in the table:
Middletownbooks/Hadamard had published an independent matrix of that order on 9 October 2026, before ours of 14 October,
and we learned of it on 19 October [16] ([hadamard/2292](hadamard/2292/README.md#status)). The same repository settled
2588 = 4 · 647 on 14 October 2026 [16]. Đoković claimed 1527 and 2545 in 2010 through Seberry–Yamada's Corollary 8.8(i),
whose printed proof does not give Williamson-type matrices ([NOTES.md](NOTES.md) §7). Doubling 2636, 3732 and 9172
settles the further open orders 5272, 7464, 18344, 36688 and 73376.

### Other results

- **Hadamard matrices beyond Table 4**: 207 orders from 12052 to 51580. 158 of them are listed as unknown in Table 1.53
  of the Handbook or lie above its range, and 49 are recorded as known without a valid or traceable construction. Most
  come from the two block decompositions, the rest from the recursive theorems of Miyamoto [8], Yamada–Kiyasu [9] and
  Cooper–Wallis [7] applied to matrices here and to objects published in 2026
  ([hadamard](hadamard#orders-beyond-the-tables)).
- **Corollary 8.8(i).** 58 orders were claimed or recorded as known only through Seberry–Yamada's Corollary 8.8(i):
  ten that Table 4 omits (2076 to 10444), 6108 and 10180, and 46 beyond Table 4. We give the first valid constructions of all 58
  ([hadamard](hadamard#orders-recorded-only-through-corollary-88i)). The corollary's matrices still give a complex
  Hadamard matrix of order 2q from any quaternary Hadamard matrix of order (q − 1)/2 ([NOTES.md](NOTES.md) §19).
- **Skew-Hadamard matrices** of 101 orders: Goethals–Seidel difference families with a skew block found by search,
  Seberry's theorem with new amicable pairs and orthogonal designs, Wallis's product, and the 3(q + 1) theorem, which
  settles 52 values v < 2500 that the Đoković–Seberry census [4] lists as undecided ([skew](skew)).
- **Symmetric Hadamard matrices** of 53 orders that the census lists as undecided: 532, 732, 876, 4828 and 5700 from
  searches and plug-in constructions, and 48 orders 4q (q ≡ 1 mod 8 a prime power) from Seberry–Yamada's
  Williamson-type matrices with a repeated block in Seberry and Balonin's propus array, among them 7012 = 4 · 1753 from
  our matrix of order 876 ([symmetric/PROPUS_MIYAMOTO.md](symmetric/PROPUS_MIYAMOTO.md)).
- **Complex Hadamard matrices** of the 12 orders 214, 238, 302, 326, 366, 382, 426, 494, 518, 534, 754 and 1558,
  blank in Seberry and Yamada's table [5] and in the papers since [10, 11], from searches over Szöllősi's array [10]
  and its bordered, character-bordered and twisted forms. Zaitsev et al. found all twelve concurrently and
  independently [13]; we claim no priority. The smallest open order is 334. Also the first explicit matrices of 1618
  and 1906, and explicit matrices of 76 orders 2N with 4N − 1 a prime power, which Ito's dicyclic Hadamard groups and
  Egan's Theorem 4.5 [12] give ([complex-hadamard](complex-hadamard)).
- **T-matrices** of orders 97, 109 and 197 ([t-matrices](t-matrices)). 97 and 109 have the form 6m + 1 and may appear
  in Lou, Xia and Xia's work on T-matrices of that order, which we could not read.

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

| Construction | Orders | Details |
|---|---|---|
| Goethals–Seidel difference families invariant under a multiplier subgroup, found by an extended k-tree search | 2092, 3532, 4492, 5932 | [METHOD.md](METHOD.md) |
| The same with blocks permuted by a multiplier of order 3 | 2284 | [METHOD.md](METHOD.md#block-symmetries) |
| Bordered Goethals–Seidel arrays with exact compression chains | 2148, 2636, 3964 | [METHOD.md](METHOD.md#bordered-structures), [search/NOTES.md](search/NOTES.md) §18–19 |
| Cyclotomic products ℤ_m × GF(q) of index 4 to 10 | 3732, 5020; skew 820, 1044, 1204, 1276, 3380, 5668 | [search/NOTES.md](search/NOTES.md) §16, §20 |
| Spence's construction with blocks (A, B, B∘μ, B∘μ²) and a Legendre lift | 2292 | [NOTES.md](NOTES.md) §13, [search/spence_spin](search/spence_spin) |
| Paley's matrix in sixteen circulant blocks, times T-matrices | 157 orders, 14 of them unresolved | [hadamard/PALEY_BLOCKS.md](hadamard/PALEY_BLOCKS.md) |
| Miyamoto's matrix in sixteen commuting bordered blocks, times T-matrices | 49 orders, 3 of them unresolved | [hadamard/MIYAMOTO_BLOCKS.md](hadamard/MIYAMOTO_BLOCKS.md) |
| Recursive theorems of Miyamoto, Yamada–Kiyasu and Cooper–Wallis | 9172 and orders beyond Table 4 | [hadamard](hadamard), [williamson-type](williamson-type), [t-matrices](t-matrices) |
| Skew difference families: compression chains, Paley products, low-energy k-trees | skew orders | [skew](skew), [search/NOTES.md](search/NOTES.md) §10–16, §20 |
| Paley's matrix in sixteen negacyclic blocks, Frobenius pairing, T-matrices of order 3 | skew orders 3(q + 1) | [skew/PALEY_SKEW.md](skew/PALEY_SKEW.md) |
| Propus quadruples, the symmetric m(m − 1) and m(m + 3) constructions, Miyamoto's lemma in the propus array | symmetric orders | [symmetric](symmetric), [search/NOTES.md](search/NOTES.md) §17, §27–31 |
| Szöllősi's array and its bordered, character-bordered and twisted forms | complex orders | [complex-hadamard](complex-hadamard) |

[NOTES.md](NOTES.md) and [search/NOTES.md](search/NOTES.md) are the research notebooks: the constructions, the
structures that turned out to be empty, the obstructions, and calibrations of the solution-count estimates.
[FRONTIER.md](FRONTIER.md) lists the unresolved orders up to 4 · 2999 with the cheapest search structures known for
each.

## Layout

| Path | Contents |
|---|---|
| `hadamard/<order>/`, `skew/<order>/`, `symmetric/<order>/` | `certificate.json`, `README.md`, and where stored `sequences.txt` and the matrix `H<order>.csv.gz` (one row per line, entries separated by commas) |
| `complex-hadamard/` | certificates `C<order>.json` and the constructions |
| `t-matrices/`, `near-williamson/`, `williamson-type/` | auxiliary objects |
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
16. Middletownbooks, Hadamard matrices of previously unknown orders, https://github.com/Middletownbooks/Hadamard
    (order 2292 committed October 9, 2026; order 2588 October 14, 2026).

## Citation

If you use these matrices or programs, please cite the repository ([CITATION.cff](CITATION.cff)):

> Kishan, *New Hadamard matrices*, https://github.com/vatskishan03/Hadamard, 2026.

## License

Apache License 2.0; see [LICENSE](LICENSE).
