# New Hadamard matrices

Explicit Hadamard, skew-Hadamard, symmetric Hadamard and complex Hadamard matrices of orders that had no known
construction. Every matrix comes with a certificate from which `verify/` rebuilds it and checks it in exact
arithmetic.

A Hadamard matrix of order n is an n × n matrix H with entries ±1 and H Hᵀ = n I. It can exist only for n = 1, 2 or
n ≡ 0 (mod 4), and the Hadamard conjecture asserts that it exists for every multiple of 4. In August 2026 Alpöge,
Voinov and Reynolds-Haertle constructed the twelve unresolved orders below 2000 [1], and Schneider constructed
2060 [2]. This repository constructs the next four, 2092, 2148, 2284 and 2292, and seven more open orders up to 9172.
The smallest order we know to be unresolved is now 2396 = 4 · 599 ([FRONTIER.md](FRONTIER.md)).

## Results

**Unresolved Hadamard orders.** For each order 4n below, n is listed as unresolved in Table 4 of Cati and
Pasechnik [3], and we found no later construction.

| Order | n | Construction |
|------:|---|---|
| [2092](hadamard/2092) | 523 | Goethals–Seidel array [6] over ℤ₅₂₃, multiplier subgroup of order 9 |
| [2148](hadamard/2148) | 537 = 3 · 179 | bordered Goethals–Seidel array over ℤ₃ × ℤ₁₇₈, border of size 12 |
| [2284](hadamard/2284) | 571 | Goethals–Seidel array over ℤ₅₇₁, three blocks related by a multiplier of order 3 |
| [2292](hadamard/2292) | 573 = 3 · 191 | Spence's construction over ℤ₂ × ℤ₂ × ℤ₁₄₃, three blocks related by an automorphism of order 15 |
| [2636](hadamard/2636) | 659 | bordered Goethals–Seidel array over ℤ₄ × ℤ₄ × ℤ₄₁, border of size 12 |
| [3732](hadamard/3732) | 933 = 3 · 311 | Goethals–Seidel array over ℤ₃ × GF(311), multiplier subgroup of order 31; skew |
| [3964](hadamard/3964) | 991 | bordered Goethals–Seidel array over ℤ₂ × ℤ₂ × ℤ₂₄₇, border of size 12 |
| [4492](hadamard/4492) | 1123 | Goethals–Seidel array over ℤ₁₁₂₃, multiplier subgroup of order 33 |
| [5020](hadamard/5020) | 1255 = 5 · 251 | Goethals–Seidel array over ℤ₅ × GF(251), multiplier subgroup of order 50 |
| [5932](hadamard/5932) | 1483 | Goethals–Seidel array over ℤ₁₄₈₃, multiplier subgroup of order 39 |
| [9172](hadamard/9172) | 2293 | Miyamoto's theorem [8] from the matrix of order 2292 |

**Other results.**

- **Hadamard matrices of 19 larger orders** from 12596 to 51580, above the range of Table 4, for which we found no
  earlier construction. They come from recursive theorems (Miyamoto [8], Yamada–Kiyasu [9], Cooper–Wallis [7])
  applied to the matrices here and to objects published in 2026. For three of them the Handbook records the order as
  known, with no construction we could trace. Also independent matrices of the known orders 2060 and 3532. See
  [hadamard](hadamard).
- **Skew-Hadamard matrices of 40 orders** from 820 to 9900 for which we found no earlier construction: 32 from
  Goethals–Seidel difference families with a skew block, and 8 from Seberry's theorem with new amicable Hadamard
  matrices or new orthogonal designs. Also independent matrices of the known skew orders 1852 and 3676.
  See [skew](skew).
- **Symmetric Hadamard matrices** of orders 532, 732, 876, 4828 and 5700, all undecided in the Đoković–Seberry
  census [4]. See [symmetric](symmetric).
- **Complex Hadamard matrices** of 16 orders 2N with N odd, from 142 to 754, all blank in Seberry and Yamada's table
  of known orders [5]. 142 and 154 were the two smallest open orders after the work of 2026 [10, 11]. See
  [complex-hadamard](complex-hadamard).
- **Auxiliary objects:** Williamson-type matrices of order 111 ([near-williamson](near-williamson)); T-matrices of
  orders 97, 109, 131 and 197, for which we found no earlier construction ([t-matrices](t-matrices)); amicable
  Hadamard matrices of orders 52, 76 and 100; and the orthogonal designs OD(100; 1, 5, 94), OD(156; 1, 3, 152),
  OD(156; 1, 13, 142), OD(252; 1, 21, 230) and OD(364; 1, 13, 350) ([skew](skew#seberrys-theorem-with-new-ingredients)).

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
| Bordered Goethals–Seidel arrays: exact compression chains, and a split along ℤ₂ | 2148, 2636, 3964 | [METHOD.md](METHOD.md), [search/NOTES.md](search/NOTES.md) §18–19 |
| Cyclotomic products ℤ_m × GF(q) of index 4 to 10 | 3732, 5020; skew 820, 1044, 1204, 1276, 3380, 5668 | [search/NOTES.md](search/NOTES.md) §16, §20 |
| Spence's construction with blocks (A, B, B∘μ, B∘μ²) and a Legendre lift | 2292 | [NOTES.md](NOTES.md) §13, [search/spence_spin](search/spence_spin) |
| Skew difference families: compression chains, Paley products, low-energy k-trees | skew orders | [skew](skew), [search/NOTES.md](search/NOTES.md) §10–16, §20 |
| Propus quadruples and the symmetric m(m − 1), m(m + 3) constructions | symmetric orders | [symmetric](symmetric), [search/NOTES.md](search/NOTES.md) §17, §27–31 |
| Szöllősi's array, its bordered and twisted forms | complex orders | [complex-hadamard](complex-hadamard) |
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
8. Y. Miyamoto, A construction of Hadamard matrices, J. Combin. Theory Ser. A 57 (1991) 86–108.
9. M. Yamada, Some new series of Hadamard matrices, J. Austral. Math. Soc. Ser. A 46 (1989) 371–383.
10. Szöllősi, A complex Hadamard matrix of order 94, arXiv:2603.09572 (2026).
11. H. Kharaghani, A. Mohammadian, B. Tayfeh-Rezaie, V. Zaitsev, arXiv:2605.08661v2 (2026).

## Citation

If you use these matrices or programs, please cite the repository ([CITATION.cff](CITATION.cff)):

> Kishan, *New Hadamard matrices*, https://github.com/vatskishan03/Hadamard, 2026.

## License

Apache License 2.0; see [LICENSE](LICENSE).
