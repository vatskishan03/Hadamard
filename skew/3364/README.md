# Skew-Hadamard matrix of order 3364

3364 = 4 · 841, with 841 = 29².

## Construction

A skew-Hadamard matrix is a Hadamard matrix H with H + Hᵀ = 2I.

**Difference family.** Four ±1 sequences x₁, …, x₄ on the group G = ℤ₂₉ × ℤ₂₉ satisfy two conditions:

- each is constant on the orbits of the multiplier group H = ⟨w⁴⁰⟩ ≤ GF(29²)*, of order 21, the group of 40th powers,
  where G is the additive group of GF(29²) = 𝔽₂₉[x]/(x² + 5x + 2), w is the class of x (a primitive element) and the
  element (c₀, c₁) of G stands for c₀ + c₁w, at position 29c₀ + c₁ of a sequence. On the coordinates the generator of H is
  the matrix [[10, 6], [26, 25]], which is the 40th power of the matrix [[0, 27], [1, 24]] of multiplication by w; it is
  an automorphism of G, and the certificate states it;
- their periodic autocorrelations over G sum to zero at every nonzero shift.

H has 41 orbits on G: the fixed point 0 and the 40 cosets of H in GF(29²)*, each of size 21; the coset of wⁱ is the
class i mod 40. The four sequences are −1 at 0 in x₃ and x₄ and +1 at 0 in x₁ and x₂, and on the classes below (each
class is an orbit of size 21):

| | x(0) | classes j (0 ≤ j < 40) on which the sequence is −1 |
|---|---:|---|
| x₁ | +1 | 2, 6, 7, 8, 10, 11, 12, 13, 14, 15, 16, 17, 18, 20, 21, 23, 24, 25, 29, 39 |
| x₂ | +1 | 2, 4, 5, 6, 9, 14, 15, 16, 17, 18, 19, 21, 23, 26, 27, 29, 36, 37, 38, 39 |
| x₃ | −1 | 1, 2, 7, 9, 10, 18, 22, 23, 25, 26, 27, 29, 31, 32, 35, 36, 37, 38, 39 |
| x₄ | −1 | 0, 6, 8, 9, 10, 13, 15, 17, 22, 26, 27, 28, 30, 31, 33, 34, 35, 36, 39 |

So the sequences are −1 on 20, 20, 19 and 19 of the 40 nonzero orbits, and on 420, 420, 400 and 400 positions (x₃ and x₄
include the point 0): a 4-{841; 420, 420, 400, 400; 799} supplementary difference set in G (λ = Σ kᵢ − 841 = 799, and
λ · 840 = 671,160 = Σ kᵢ(kᵢ − 1)). `certificate.json` lists the four sequences in full, `sequences.txt` one per line.

**Skew first block.** The first sequence is of skew type: x₁(0) = 1 and x₁(−t) = −x₁(t) for t ≠ 0. This is possible
because −1 ∉ H (H has odd order): −1 = w⁴²⁰ is a 20th power but not a 40th power, so negation maps the class j to the
class j + 20 and pairs the 40 nonzero orbits into 20 pairs, and x₁ is −1 on exactly one orbit of each pair.

**Row sums.** The row sums are 1, 1, 41, 41, and 1² + 1² + 41² + 41² = 3364. A sequence constant on the orbits of H has
row sum x(0) + 21 · (an even number), that is ≡ x(0) (mod 42), and 3364 = 1 + b² + c² + d² has only one solution in
numbers ≡ ±1 (mod 42), up to order and sign, (b, c, d) = (1, 41, 41). With the row sums taken positive (a block may be
negated in the array), 1 ≡ 1 forces x(0) = +1 in the first two blocks and 41 ≡ −1 (mod 42) forces x(0) = −1 in the other
two.

**The matrix.** Let A, B, C, D be the group-developed matrices X[g, h] = x(g − h) of the four sequences and R the
permutation matrix with R[g, h] = 1 when g + h = 0. The Goethals–Seidel array

```
    [  A     BR     CR     DR  ]
H = [ -BR    A     DᵀR   -CᵀR ]
    [ -CR   -DᵀR    A     BᵀR ]
    [ -DR    CᵀR  -BᵀR    A   ]
```

is then a Hadamard matrix of order 4 · 841 = 3364. Because A + Aᵀ = 2I, it satisfies H + Hᵀ = 2I (Wallis–Whiteman).

## Search

Numbers below are those printed by the log of the run, except those marked (recomputed).

1. **Structure.** The group is not cyclic, so this is the prime-power case of the cyclotomic product of
   [search/NOTES.md](../../search/NOTES.md) §21 with no factor ℤ_m (`cp_lowe.py --field 29,2`; `search/cp_field.py`
   takes the first primitive polynomial over 𝔽₂₉ of degree 2, x² + 5x + 2). Each block is a choice of sign on the 40
   classes and at the point 0. The row sums (1, 1, 41, 41) are forced (above), which leaves one class of block shapes:
   - 2²⁰ = 1,048,576 skew first blocks;
   - C(40, 20) = 137,846,528,820 choices for B (x(0) = +1, 20 classes −1);
   - C(40, 19) = 131,282,408,400 choices each for C and D (x(0) = −1, 19 classes −1).

   The complementarity condition is Σₖ Pₖ(s) = 0 on the 20 orbits of ⟨H, −1⟩ on the nonzero shifts, of rank 19 once the
   row sums are fixed (the log: key dimension 20, rank 19). The Gaussian estimate of the number of solutions of the
   class is 10^9.02.
2. **Low-energy k-tree** (`search/cp_lowe.py`, `search/cp_enum.c`, `search/kt4.c`;
   [search/NOTES.md](../../search/NOTES.md) §20–22, §110–111, [search/listgen/README.md](../../search/listgen/README.md)).
   - `cp_enum.c` enumerates the class-sign vectors of B, C and D (the 40 classes split into 5 groups of 8, `--groups 5`).
     It keeps the members whose autocorrelation vector is closest to its conditional mean given that the four vectors
     sum to zero (whitened energy at most τ, rank 19): τ = 6.51, 4.51 and 4.50 keep 588,899,837, 28,751,728 and
     28,469,572 members of B, C and D, and all 1,048,576 skew blocks are kept. The Gaussian estimate inside these lists
     is 10^1.85 solutions; it is a heuristic, not a rate. The twin reduction of `search/listgen/cp_lowe_dd.py` (reversing
     a block, x(t) → x(−t), leaves its autocorrelations unchanged) was not used, so the lists hold the reversal twins.
   - `kt4.c` joins the four lists over 10⁵ split values of five merge classes (1, 7, 14, 15, 18), chosen from 1,125,425
     sampled values and ranked by expected solutions per unit of work (they hold 0.411 of the expected solutions;
     predicted 1.57 · 10¹³ pairs on the AB side and 1.02 · 10¹⁴ on the CD side in all). The first 1128 split values
     gave nothing: 7.6 · 10¹⁰ AB pairs (6.7 · 10⁷ on average, 3.6 · 10⁷ to 8.7 · 10⁷) and 2.0 · 10¹² CD pairs (1.8 · 10⁹
     on average, 3.0 · 10⁸ to 3.5 · 10⁹), 0.48 % and 2.0 % of the predicted totals (recomputed from the 1128 lines of
     the log); they took 15,009 s on kt4's own clock, 13.3 s each (3 to 27 s). The next split value was not finished:
     the log ends with three solutions (the limit was K = 3), printed at 16,310 s from the start, and with no line for
     the split value that held them, so the number of pairs of that split value and the number of its solutions are not
     known.
   - The three printed solutions were each rebuilt and checked (autocorrelations, row sums (1, 1, 41, 41), invariance
     under H). The second is the first with the block x₂ reversed; the third is not a reversal twin of either (recomputed
     from the hit records). The certificate is the third (`skew_certify_group.py 841 2`).
3. **Earlier attempts.** None recorded for this structure: neither the repository nor the run directories of the
   machine show an earlier run of v = 841 as a skew Goethals–Seidel family. The only earlier search on ℤ₂₉ × ℤ₂₉ is the
   sweep of negacyclic Goethals–Seidel arrays of order 8v = 6728 ([search/NOTES.md](../../search/NOTES.md) §88, a
   different construction), in which the GF(841) shapes are empty.

**Reproduction.**

| | |
|---|---|
| Code | the programs of commit `919e73d`, byte for byte, run unchanged: `search/cp_lowe.py` (md5 `4e2f4c952a56d1d0942b5c1c1a7d4d70`), `search/cp_field.py` (`3e7d729e5dc9abc18876bc5b9dbef55f`), `search/cp_enum.c` (`4c2094d92d03faff5b970b136c837958`), `search/kt4.c` (`d22c49d986b9b3ba9ee6d3e62ded6fea`), `search/cyclo_product.py` (`98be16555ffe42590a0cd2185add9f3e`), `search/workdir.py` (`b845fde4e2f09a184b42ec79d7776f4b`), `verify/verify.py` (`970e5d9aba849fd38a7a4109a85a091d`), and `search/skew_certify_group.py` (`2759b46f5e366a4d26785c2c5f9268c1`, with `search/certify.py` `f27dd2834f970c8424a1d6f45a94e896`) for the certificate. The md5 of each of these files on the machine of the run is that of the file in the repository. Nothing is added to `search/` with this order |
| Command | `python3 search/cp_lowe.py 841 --q 841 --e 40 --config 0 --L 6e8,3e7,3e7 --nsig 100000 --threads 12 --K 3 --seed 1 --groups 5 --samples 2000000 --mcpairs 10000000 --maxpairs 300000000 --stream --field 29,2`, one round, with `HADAMARD_WORK` the run directory and the list files on a RAM disk |
| Randomness and caps | `--seed 1` (numpy `default_rng`) draws the 2 · 10⁶ samples per block that set the energy thresholds and the 10⁷ pairs per side that rank the split values; `--L 6e8,3e7,3e7` the target number of members kept per non-skew block (the threshold is a quantile of the sampled energies, so the lists hold 588,899,837, 28,751,728 and 28,469,572); `--nsig 100000` split values; `--maxpairs 3e8` skips a split value with more AB pairs (none of the 1128 logged did: the largest had 8.7 · 10⁷); `--K 3` stops kt4 at the third solution; `--kgen 1 --twist 0` (the defaults) |
| Result | three verified solutions, logged 16,310 s (4.53 h) after the start on 12 threads. Elapsed times from the start: key tables and the skew block 53 s and 54 s; the list of B 459 s; C 758 s; D 1051 s; kt4 started at 1058 s. The scheduler of the machine recorded 16,315 s of wall time, 195,779 core-seconds (54.4 core-hours counted as 12 times the wall time) and a peak resident memory of 45.1 GB (its periodic sample of the process group) |
| Machine | 32-core AMD EPYC 9R14 shared with other jobs (12 threads used), 246 GB, Python 3.12.3, numpy 1.26.4, gcc 13.3.0 (`-O3 -march=native -fopenmp`; the binaries that ran have md5 `fd442e1ee5f538eda5e8dd8be1ae524e` for `cp_enum` and `ffe1b35a47b8bc1fb8d2cbb5fda2e986` for `kt4`) |
| Certificate | `search/skew_certify_group.py 841 2` turns the third logged hit into the files in this directory |

## Status

**New** as a skew-Hadamard matrix. An ordinary Hadamard matrix of order 3364 is classical (see below); the skew
property is what is new.

- The Đoković–Seberry census lists 841 among its 737 undecided odd v < 2500 for skew matrices, and not in its Paley,
  Szekeres, Szekeres–Whiteman or Spence series or among its sporadic values. For symmetric matrices it is decided
  there, in the Paley series v = (q + 1)/2 with q = 1681. The census says that it may have missed some orders.
- Cati–Pasechnik Table 5 (v2, 30 August 2025; entries n(m) give the smallest m with a skew matrix of order 2ᵐn known,
  and m = 2 is not shown) prints 841(6): the smallest skew order with odd part 841 that the table records is
  2⁶ · 841 = 53824. Among the odd squares of the table 841 is the outlier (289(3), 361(3), 441(3), 529(3), 625(3),
  729(4), 961(3)), and the table omits a published family: Momihara and Xiang (arXiv:1801.08776, Corollary 1.7) give a
  skew-Hadamard matrix of order 4(qᵉ + 1) for every prime power q ≡ 9 (mod 16) and every e ≥ 1, and q = 41, e = 2 gives
  4 · 1682 = 6728 = 2³ · 841. Neither Momihara–Xiang nor Wallis–Whiteman is among the 60 references of Table 5, and
  [NOTES.md](../../NOTES.md) §28 records the omission. So the smallest order with odd part 841 known before this one
  was 2³ · 841 = 6728, and this matrix lowers the exponent from 3 (6 as printed) to 2.
- The other skew tables do not settle it: Koukouvinos–Stylianou (2008, Table 1) covers odd m < 500, and the rows for 841
  of Seberry–Yamada's Table 7.1 (1992) and Table 9.2 (2020) were not read. None of the families SI–SXI of the 1992
  survey reaches 3364 (below). SageMath's list `unknown_skew_hadamard` stops at 1196, below this order, and no branch of
  `skew_hadamard_matrix` reaches 3364 as far as its docstrings and source say (not run).
- The ordinary matrix of order 3364 is classical: Paley II with q = 1681 = 41² (the Hadamard Matrix Explorer builds 3364
  as `P2(3364)`), Turyn's Williamson matrices of order (q + 1)/2 = 841 for q = 1681 (1972), and the Williamson-type
  matrices of order 4q² for prime powers q ≡ 1 (mod 4) of Xia and Liu (1991), over the same group GF(q²)⁺ (here q = 29).
  All of these are symmetric, or have four symmetric blocks, so none is skew; the skew property is what is new.
- No classical construction gives a skew matrix of order 3364:
  - Paley I: 3363 = 3 · 19 · 59 is not a prime power.
  - Szekeres, Whiteman, Spence and Momihara–Xiang 2(qᵉ + 1): qᵉ = 1681 = 41², and neither 41 nor 1681 is ≡ 5 (mod 8).
  - Wallis–Whiteman and Momihara–Xiang 4(qᵉ + 1): 840 is not a prime power (it does give 4 · 1682 = 6728).
  - Spence 4(q² + q + 1): q² + q + 1 = 841 has no integer solution (the discriminant 3361 is not a square).
  - Turyn (h − 1)ᵘ + 1 from a skew matrix of order h, u odd: 3363 is not a perfect power, so only u = 1 (the order
    itself) remains.
  - Doubling needs a skew matrix of order 1682, which is not a multiple of 4.
  - Skew × amicable products need 3364 = mn with m and n Hadamard orders, and the only such factorisation is 1 · 3364
    (2 · 1682, 4 · 841, 29 · 116 and 58 · 58 have a factor 1682, 841, 29 or 58). Williamson's h(pᵏ + 1), with
    pᵏ + 1 ≡ 0 (mod 4), needs pᵏ + 1 = 4 and a skew matrix of the odd order h = 841, which cannot exist, or
    pᵏ + 1 = 116 or 3364, and 115 and 3363 are not prime powers.
  - Seberry's mn(n − 1) is possible only with n = 2, which needs an OD(3364; 1, 1682, 1681).
  - Yamada's (|t| + 1)(q + 1), with q = s² + 4t² ≡ 5 (mod 8) a prime power, has q + 1 ≡ 6 (mod 8), and no divisor of
    3364 is. Seberry's 2ᵗq with q = s² + 4r² ≡ 5 (mod 8) a prime power needs q = 841, which is ≡ 1 (mod 8).
    (Families SVIII and SX of the list of known skew orders in Seberry–Yamada 1992.)
  - The 3(q + 1) theorem needs 3 | 3364. Blatt–Szekeres needs complementary difference sets in a group of order
    2 · 840 + 1 = 1681 = 41²; the known families need 1681 = pᵗ with p ≡ 5 (mod 8) (SageMath's
    `complementary_difference_sets`) or 2 · 1681 + 1 = 3363 to be a prime power, and 41 ≡ 1 (mod 8) and
    3363 = 3 · 19 · 59.
  - The two nearest cyclotomic searches, with unions of 16th-power classes of GF(1681) in the bordered two-block family
    of order 2(1681 + 1) and with unions of 8th-power classes and the point 0 of GF(841) in the four-block family, find
    nothing ([NOVELTY.md](../../NOVELTY.md)). The family here has index 40.

The sources searched are listed in [NOVELTY.md](../../NOVELTY.md). Doubling (Seberry) gives skew-Hadamard matrices of
orders 6728, 13456, …, which were already known, since 6728 = 8 · 841 is Momihara and Xiang's.

## Files

| File | Contents |
|---|---|
| `certificate.json` | the group ℤ₂₉ × ℤ₂₉, the automorphism matrix of H, the four sequences, the row sums |
| `sequences.txt` | the four sequences, one per line |
| `H3364.csv.gz` | the matrix |

## Verification

```
python3 verify/verify.py skew/3364
gzip -dc skew/3364/H3364.csv.gz | ./check_matrix
```

`check_matrix` read the stored matrix (5,659,930 row pairs, no violation; run on 2 November 2026).

SHA-256 of `H3364.csv`: `50243e1b4cf345ea27de6393f13d8621f8eb222ce697b16832d15a85db68054c`
