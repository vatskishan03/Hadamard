# Skew-Hadamard matrix of order 3044

3044 = 4 · 761 (prime).

## Construction

A skew-Hadamard matrix is a Hadamard matrix H with H + Hᵀ = 2I.

**Difference family.** Four ±1 sequences x₁, …, x₄ on ℤ₇₆₁ satisfy two conditions:

- each is constant on the orbits of the multiplier subgroup H = ⟨498⟩ ≤ ℤ₇₆₁*, of order 19, which is the group of 40th
  powers modulo 761: {1, 25, 152, 232, 233, 258, 274, 357, 362, 405, 410, 473, 498, 554, 625, 636, 679, 680, 756};
- their periodic autocorrelations sum to zero at every nonzero shift.

H has 41 orbits on ℤ₇₆₁: the fixed point 0 and the 40 cosets of H in ℤ₇₆₁*, each of size 19. `certificate.json` lists,
for each sequence, one representative r of each orbit rH on which the sequence is −1. All four sequences are +1 at 0,
and −1 on 20, 20, 19 and 19 of the 40 nonzero orbits (380, 380, 361 and 361 positions), a 4-{761; 380, 380, 361, 361;
721} supplementary difference set.

**Skew first block.** The first sequence is of skew type: x₁(0) = 1 and x₁(−t) = −x₁(t) for t ≠ 0. This is possible
because −1 ∉ H (H has odd order; −1 is a 20th power modulo 761 but not a 40th power), so negation pairs the 40 nonzero
orbits into 20 pairs (the coset of index j with the coset of index j + 20) and x₁ is −1 on exactly one orbit of each
pair.

**Row sums.** The row sums are 1, 1, 39, 39, and 1² + 1² + 39² + 39² = 3044. A sequence constant on the orbits of H has
row sum x(0) + 19 · (an even number), that is ≡ x(0) (mod 38), and 3044 = 1 + b² + c² + d² has only one solution in
numbers ≡ ±1 (mod 38), up to order and sign, (b, c, d) = (1, 39, 39). With the row sums taken positive (a block may be
negated in the array), 1 ≡ 39 ≡ 1 (mod 38) forces x(0) = +1 in every block.

**The matrix.** Let A, B, C, D be the circulant matrices with these first rows (`sequences.txt`) and R the
back-diagonal identity. The Goethals–Seidel array

```
    [  A     BR     CR     DR  ]
H = [ -BR    A     DᵀR   -CᵀR ]
    [ -CR   -DᵀR    A     BᵀR ]
    [ -DR    CᵀR  -BᵀR    A   ]
```

is then a Hadamard matrix of order 3044. Because A + Aᵀ = 2I, it satisfies H + Hᵀ = 2I (Wallis–Whiteman).

## Search

1. **Structure.** Over the prime field ℤ₇₆₁ each block is a choice of sign on the 40 orbits of H, with the value +1 at
   the point 0. The row sums (1, 1, 39, 39) are forced, which leaves one class of block shapes:
   - 2²⁰ = 1,048,576 skew first blocks;
   - C(40, 20) = 137,846,528,820 choices for B;
   - C(40, 19) = 131,282,408,400 choices each for C and D.

   The complementarity condition is Σₖ Pₖ(s) = 0 on the 20 orbits of ⟨H, −1⟩ on the nonzero shifts, of rank 19 once the
   row sums are fixed. The Gaussian estimate of the number of solutions of the class is 10^9.50.
2. **Reversal twins.** Reversing a sequence (x(t) → x(−t)) leaves its autocorrelations unchanged, so reversing any of
   the four blocks of a solution gives another solution: every solution lies in a class of up to 16 twins. On the 40
   class signs, reversal is the rotation by 20 (for the skew block it negates all nonzero signs and keeps x(0) = 1).
   Every list below keeps one member of each twin: for B, C and D the member whose class mask is at most its rotation by
   20, for the skew block the 524,288 members with sign + on the first class.
3. **Low-energy k-tree** (`search/listgen/cp_lowe_dd.py`, `search/cp_enum.c`, `search/kt4.c`;
   [search/NOTES.md](../../search/NOTES.md) §20–22, §110–111,
   [search/listgen/README.md](../../search/listgen/README.md)).
   - `cp_enum.c` enumerates the class-sign vectors of B, C and D (the 40 classes split into 5 groups of 8,
     `--groups 5`). It keeps the members whose autocorrelation vector is closest to its conditional mean given that the
     four vectors sum to zero (whitened energy at most τ, rank 19): τ = 6.57, 4.63 and 4.57 keep 590,647,863,
     33,937,598 and 30,532,180 members of B, C and D, which the twin reduction cuts to 295,323,935, 16,968,798 and
     15,266,090. With the 524,288 skew blocks (all kept) the Gaussian estimate inside these lists is 10^1.21 solutions
     (counted once per twin class; 10^2.41 before the reduction); it is a heuristic, not a rate.
   - `kt4.c` joins the four lists over 10⁵ split values of five merge classes (2, 7, 8, 12, 14), chosen from 1,087,935
     sampled values and ranked by expected solutions per unit of work (they hold 0.442 of the expected solutions;
     predicted 4.02 · 10¹² pairs on the AB side and 3.86 · 10¹³ on the CD side in all), with on average 1.76 · 10⁷ pairs
     on the AB side (8.5 · 10⁶ to 2.3 · 10⁷) and 6.3 · 10⁸ on the CD side (7.7 · 10⁷ to 1.4 · 10⁹) of each. The first
     2702 split values gave nothing. The 2703rd gave two solutions; kt4 printed the first (its limit is 1), which was
     verified. The run stopped after 4.77 · 10¹⁰ AB pairs and 1.70 · 10¹² CD pairs, 1.2 % and 4.4 % of the predicted
     totals.
4. **Earlier attempts.** Both runs below used the same class without the twin reduction and found no solution before
   they were stopped.
   - Lists of about 10⁸ members per non-skew block (100,007,217, 92,730,345 and 102,496,968 kept, all 1,048,576 skew
     blocks, estimate 10^2.60 inside): 44 split values in 4614 s on 12 threads.
   - The lists above before the reduction (1,048,576, 590,647,863, 33,937,598 and 30,532,180 members, estimate 10^2.41):
     377 split values in 11,116 s on 12 threads, 2.77 · 10¹⁰ AB pairs and 1.09 · 10¹² CD pairs, about 46 GB of memory.
     Per split value that is 29.5 s against 3.2 s with the reduction (4.2 times the AB pairs and 4.6 times the CD
     pairs).

**Reproduction.**

| | |
|---|---|
| Code | `search/listgen/cp_lowe_dd.py` (md5 `2b9f668d4a99969bebc8ea4e2b3bc317`), which is `search/cp_lowe.py` plus the option `--dedup` (the twin reduction, above) and `--only-sigma` (a validation option, not used here), run through the resumable wrapper `search/listgen/dd_run.sh` (md5 `be9c5fafeb317108cde6f8ca2f5a58af`; the file that ran, md5 `f9a321939f884403ede551e9a0355432`, differs from it only in the line that gave `HERE`, the directory of `cp_lowe_dd.py`, as an absolute path on the machine of the run, where the published file takes it from the directory of the script); both ran from a working copy and are added with this order (`cp_lowe_dd.py` byte for byte) and are not in commit `919e73d`. The other programs are those of commit `919e73d`, byte for byte: `search/cp_enum.c` (`4c2094d92d03faff5b970b136c837958`), `search/kt4.c` (`d22c49d986b9b3ba9ee6d3e62ded6fea`), `search/cyclo_product.py` (`98be16555ffe42590a0cd2185add9f3e`), `search/workdir.py` (`b845fde4e2f09a184b42ec79d7776f4b`), `verify/verify.py` (`970e5d9aba849fd38a7a4109a85a091d`), and `search/skew_certify.py` (`813b89691dad146f8d68b217b292659a`, with `search/certify.py` `f27dd2834f970c8424a1d6f45a94e896`) for the certificate. On the machine of the run `cp_lowe_dd.py` stood in the directory of `cp_enum.c`, `kt4.c`, `cyclo_product.py` and `workdir.py`, with `verify/verify.py` in the sibling directory `verify`, as `search/cp_lowe.py` does. Here it is in `search/listgen/`, where the two directories it puts on `sys.path` do not hold those modules: `PYTHONPATH=search` supplies them ([search/listgen/README.md](../../search/listgen/README.md)). |
| Command | `PYTHONPATH=search python3 search/listgen/cp_lowe_dd.py 761 --q 761 --e 40 --config 0 --L 6e8,3e7,3e7 --dedup --stream --groups 5 --samples 2000000 --mcpairs 10000000 --maxpairs 300000000 --nsig 100000 --K 1 --seed 1 --threads 12`, one round, with `HADAMARD_WORK` the run directory and the list files on a RAM disk (`search/listgen/dd_run.sh WORKROOT` followed by the arguments above, with WORKROOT and `PYTHONPATH` absolute paths) |
| Randomness and caps | `--seed 1` (numpy `default_rng`) draws the 2 · 10⁶ samples per block that set the energy thresholds and the 10⁷ pairs per side that rank the split values; `--L 6e8,3e7,3e7` the target number of members kept per non-skew block before the twin reduction (the threshold is a quantile of the sampled energies, so the lists hold 590,647,863, 33,937,598 and 30,532,180, and 295,323,935, 16,968,798 and 15,266,090 after it); `--nsig 100000` split values; `--maxpairs 3e8` skips a split value with more AB pairs (none did: the largest had 2.3 · 10⁷); `--K 1` stops kt4 at the first solution; `--kgen 1 --twist 0` (the defaults) |
| Result | one verified solution after 9955 s (2.77 h) on 12 threads, 30.7 GB peak resident memory. Elapsed times from the start: key tables and the skew block 53 s; the list of B 480 s; C 782 s; D 1076 s; the ranking of the split values 1083 s; kt4 needed about 115 s more to map and index the lists; its 2703 split values then took 8755 s (3.2 s each) and the first solution was logged at 9951 s. |
| Machine | 32-core AMD EPYC 9R14 shared with other jobs (12 threads used), 246 GB, Python 3.12.3, numpy 1.26.4, gcc 13.3.0 (`-O3 -march=native -fopenmp`; the binaries that ran have md5 `fd442e1ee5f538eda5e8dd8be1ae524e` for `cp_enum` and `ffe1b35a47b8bc1fb8d2cbb5fda2e986` for `kt4`) |
| Certificate | `search/skew_certify.py` turns the logged hit into the files in this directory |

## Status

**New** as a skew-Hadamard matrix. An ordinary Hadamard matrix of order 3044 is classical (see below); the skew property
is what is new.

- The Đoković–Seberry census lists 761 among its 737 undecided odd v < 2500 for skew matrices, and not in its Paley,
  Szekeres, Szekeres–Whiteman or Spence series or among its sporadic values (its Spence series has 757 = 27² + 27 + 1
  and not 761). 761 occurs in its Turyn C-series of symmetric matrices, 2v = 1 + (h − 1)ᵏ with k even and h the order of
  a skew-Hadamard matrix: 2 · 761 = 1 + 39², h = 40.
- Cati–Pasechnik Table 5 (v2, 30 August 2025; entries n(m) give the smallest m with a skew matrix of order 2ᵐn known,
  and m = 2 is not shown) lists 761(22): the smallest known skew order with odd part 761 was 2²² · 761 =
  3,191,865,344, from Paley I (3,191,865,343 is prime and ≡ 3 (mod 4); 2ᵐ · 761 − 1 is not a prime power for m = 2, …,
  21, nor a perfect power). Table 5 updates Table 9.2 of Seberry–Yamada's 2020 book (by its §6 it lowers 67 entries
  and fills in 14 that the book lists as unknown: 191, 449, 491, 577, 643, 659, 709, 733, 739, 829, 877, 929, 933, 954
  (sic: so printed in §6, but 954 is even and Table 5 has 953(14))), and 761 is not among the 14; we infer that the
  book's entry for 761 is a number at least 22 (we did not see the book's table). Koukouvinos–Stylianou (2008), Table
  1 (odd v up to 493) and Seberry–Yamada (Utilitas Math. 1991; its abstract promises tables to q ≤ 2000, but its
  printed skew table stops at 499) stop before 761; Seberry–Yamada (1992) Table 7.1 (to 999) is built from the
  families SI–SXI, none of which gives 3044 (see below).
- The survey of Egan, Kotsireas, London and Ó Catháin (arXiv:2610.27175, 19 October 2026) lists in its Table 3 as
  known skew orders Paley's q + 1, doubling, Spence's family, 4n for odd n < 250, 276 and 292, with 356 the smallest
  open order; neither 3044 nor 761 occurs in its text (a search of the PDF's extracted text, 2 November 2026), and we
  found no statement in it that a skew matrix of order 3044 is known. SageMath's `skew_hadamard_matrix` implements the
  known orders up to 1200 and only some above (Paley I, Spence 1975 and 1977, Whiteman, complementary difference sets,
  orthogonal designs, doubling, Williamson, stored small cases); none of these gives 3044 (see below).
- The ordinary matrix of order 3044 is classical, in two ways. First, J. Seberry (Wallis), Some (1, −1) matrices, J.
  Combin. Theory Ser. B 10 (1971) 1–11 (doi:10.1016/0095-8956(71)90064-5), builds from a skew-Hadamard matrix of order h
  an n-type matrix N = I + R (R a symmetric conference matrix, R² = (n − 1)I) of order (h − 1)² + 1, here 39² + 1 = 1522
  for h = 40 (skew orders 16, 36 and 40 give n-type orders 226, 1226 and 1522), and from an n-type matrix of order n a
  symmetric Hadamard matrix of order 2n; its abstract lists 3044, and its text says "symmetric Hadamard matrices of
  orders 452, 2452, and 3044". Turyn's theorem on C-matrices (Canad. J. Math. 23 (1971) 531–535) gives the same
  symmetric conference matrix of order 1522 from the skew-Hadamard matrix of order 40, and the census's Turyn C-series
  lists v = 761 for it. Second, Miyamoto's theorem (J. Combin. Theory Ser. A 57 (1991) 86–108): for a prime power q ≡ 1
  (mod 4) and a Hadamard matrix of order q − 1 there is one of order 4q; here q = 761 and q − 1 = 760 = 2 · (379 + 1),
  the Paley I matrix of order 380 doubled. Cati–Pasechnik Table 4 does not list 761, so m = 2 there (their Theorem 6 is
  this theorem; Egan et al.'s Table 4 calls it Miyamoto I), and the Hadamard Matrix Explorer (math74.com, orders
  3001–4000) builds order 3044 exactly so, as `M(K(S(2), P1(380)))`; it prefers its Goethals–Seidel routes (`G(A(n))`,
  and `G(SK(n))` for skew-type sequences) to Miyamoto's, so by its stated order of routes it holds no such sequences for
  761.
- None of the known constructions (Seberry–Yamada's families SI–SXI, the census series, Cati–Pasechnik, SageMath) gives
  a skew matrix of order 3044:
  - Paley I and its products 2ᵗ ∏(qᵢ + 1): 3043 = 17 · 179 is not a prime power, and the divisors of 3044 that are
    ≡ 0 (mod 4) are 4 and 3044, which leaves the odd factor 761.
  - Szekeres, Whiteman, Spence and Momihara–Xiang 2(qᵉ + 1): 1521 = 3² · 13² is not a prime power.
  - Wallis–Whiteman and Momihara–Xiang 4(qᵉ + 1): 760 = 2³ · 5 · 19 is not a prime power.
  - Spence 4(q² + q + 1): 761 = q² + q + 1 has no integer solution (the discriminant 1 + 4 · 760 = 3041 is not a square;
    27² + 27 + 1 = 757 and 28² + 28 + 1 = 813).
  - Turyn (h − 1)ᵘ + 1 from a skew matrix of order h, u odd: 3043 = 17 · 179 is not a perfect power, so only u = 1 (the
    order itself) remains.
  - Doubling needs a skew matrix of order 1522, which is not a multiple of 4.
  - Skew × amicable products need 3044 = mn with m and n Hadamard orders, and the only such factorisation is 1 · 3044
    (the others have a factor 761 or 1522). Williamson's h(pᵏ + 1), with pᵏ + 1 ≡ 0 (mod 4), needs pᵏ + 1 = 4, that is a
    skew matrix of the odd order h = 761, which cannot exist (or pᵏ + 1 = 3044, and 3043 = 17 · 179 is not a prime
    power).
  - Seberry's mn(n − 1) is possible only with n = 2, which needs an OD(3044; 1, 1522, 1521).
  - Seberry–Yamada's remaining families: 2ᵗq needs q ≡ 5 (mod 8), and 761 ≡ 1 (mod 8); (|t| + 1)(q + 1) with q ≡ 5
    (mod 8) a prime power needs q + 1 to divide 3044, and none of 1, 3, 760, 1521, 3043 is such a prime power; the
    sporadic 4m (m the odd integers from 3 to 31 and 37, 39, 43, 49, 65, 67, 93, 113, 121, 127, 129, 133, 157, 163,
    181, 217, 219, 241, 267) does not contain 761.
  - The 3(q + 1) theorem needs 3 | 3044. Blatt–Szekeres needs complementary difference sets in a group of order 1521 = 2
    · 760 + 1 = 3² · 13²; the known families need 1521 or 2 · 1521 + 1 = 3043 = 17 · 179 to be a prime power, and
    neither is.
  - Good matrices and skew Williamson-type matrices of order 761 are not recorded.

The sources searched are listed in [NOVELTY.md](../../NOVELTY.md). Doubling (Seberry) gives skew-Hadamard matrices of
orders 6088, 12176, …, 2²¹ · 761 = 1,595,932,672 (nineteen orders), which are new as far as we know; 2²² · 761 was
already known.

## Files

| File | Contents |
|---|---|
| `certificate.json` | the four sequences as orbit representatives, H, the row sums |
| `sequences.txt` | the four sequences, one per line |
| `H3044.csv.gz` | the matrix |

## Verification

```
python3 verify/verify.py skew/3044
gzip -dc skew/3044/H3044.csv.gz | ./check_matrix
```

SHA-256 of `H3044.csv`: `751d175eb8a03dc995975e13f6169f485e6f0338f130e8cc4ba55314d6ed2965`
