# Skew-Hadamard matrix of order 1004

1004 = 4 · 251 (prime).

## Construction

A skew-Hadamard matrix is a Hadamard matrix H with H + Hᵀ = 2I.

**Difference family.** Four ±1 sequences x₁, …, x₄ on ℤ₂₅₁ satisfy two conditions:

- each is constant on the orbits of the multiplier subgroup H = ⟨219⟩ ≤ ℤ₂₅₁*, of order 5, which is the group of 50th
  powers modulo 251: {1, 20, 113, 149, 219};
- their periodic autocorrelations sum to zero at every nonzero shift.

H has 51 orbits on ℤ₂₅₁: the fixed point 0 and 50 orbits of size 5. `certificate.json` lists, for each sequence, one
representative r of each orbit rH on which the sequence is −1. The four sequences are −1 on 25, 25, 25 and 23 of the
51 orbits (125, 121, 121 and 111 positions), a 4-{251; 125, 121, 121, 111; 227} supplementary difference set.

**Skew first block.** The first sequence is of skew type: x₁(0) = 1 and x₁(−t) = −x₁(t) for t ≠ 0. This is possible
because −1 ∉ H (H has odd order), so negation pairs the 50 nonzero orbits into 25 pairs and x₁ is −1 on exactly one
orbit of each pair. Here x₁ is the Paley–Todd block: it is −1 exactly on the 125 nonzero squares modulo 251 (25
orbits, because H has odd order and so lies in the squares), and its autocorrelation is −1 at every nonzero shift.
The other three sequences therefore form, on their own, a 3-{251; 121, 121, 111; 165} supplementary difference set
(165 = 121 + 121 + 111 − (3 · 251 − 1)/4), with autocorrelations summing to +1 at every nonzero shift.

**Row sums.** The row sums are 1, 9, 9, 29, and 1² + 9² + 9² + 29² = 1004. A sequence constant on the orbits of H has
row sum ±1 + 5w, where w is the sum of 50 orbit signs and so even: the row sum is ≡ ±1 (mod 10). A skew first block
has row sum 1, and 1004 = 1 + b² + c² + d² has three solutions in odd numbers, up to order and sign: (b, c, d) =
(7, 15, 27), (9, 9, 29) and (11, 21, 21). The first is not ≡ ±1 (mod 10), which leaves two classes: (1, 9, 9, 29),
the one of this matrix, and (1, 11, 21, 21), the class of the non-skew matrix of Đoković, Golubitsky and Kotsireas
(see Status).

**The matrix.** Let A, B, C, D be the circulant matrices with these first rows (`sequences.txt`) and R the
back-diagonal identity. The Goethals–Seidel array

```
    [  A     BR     CR     DR  ]
H = [ -BR    A     DᵀR   -CᵀR ]
    [ -CR   -DᵀR    A     BᵀR ]
    [ -DR    CᵀR  -BᵀR    A   ]
```

is then a Hadamard matrix of order 1004. Because A + Aᵀ = 2I, it satisfies H + Hᵀ = 2I (Wallis–Whiteman).

## Search

1. **Structure.** Over the prime field ℤ₂₅₁ each block is a choice of sign on the 50 orbits of H and on the point 0. The
   row sums (1, 9, 9, 29) fix the shapes of the blocks:
   - 2²⁵ = 33,554,432 skew first blocks, 2²⁴ = 16,777,216 of them up to reversal (reversal negates the nonzero part
     of a skew block and leaves its autocorrelation unchanged);
   - C(50, 24) = 1.22 · 10¹⁴ choices each for B and C (x(0) = −1 and 24 orbits −1);
   - C(50, 22) = 8.87 · 10¹³ choices for D (x(0) = −1 and 22 orbits −1).

   The complementarity condition is Σₖ Pₖ(s) = 0 on the 25 orbits of ⟨H, −1⟩ on the nonzero shifts, of rank 24 once
   the row sums are fixed. The key of a member is its vector of 25 autocorrelations, one per orbit.
2. **Generated lists** (`search/listgen/lg_lowe.py`, `lg_lib.py`, `gen_run.py`, `lg_pt.c`;
   [search/listgen/README.md](../../search/listgen/README.md), [search/NOTES.md](../../search/NOTES.md) §111). The
   lists for B, C and D cannot be enumerated (10¹³ to 10¹⁴ members), so each holds low-energy members found by
   sampling. The energy of a member is the one of `search/cp_lowe.py`: the squared distance of its key from the
   conditional mean of its block given that the four keys sum to zero, in the metric of the conditional covariance (24
   dimensions), estimated from 4 · 10⁵ random members per block (`--seed 2`). Of a member and its reversal only one is
   kept.
   - A: the 2²⁴ skew first blocks up to reversal are enumerated by Gray code and the 343,764 of energy ≤ 10.69 are
     kept (the threshold at which the nominal Gaussian ball holds 3 · 10⁵ members). 12 s.
   - B: 12 Metropolis chains at the fixed temperature T = 0.42, one per thread, each move swapping a +1 and a −1 orbit
     sign, collect distinct members of energy ≤ 3.677 (nominal ball 7 · 10⁷) until 14,000,000 are found: 452 s, 1.62 ·
     10¹¹ proposals, 2.24 · 10⁹ accepted.
   - C: B and C have the same row sum, so C is taken from the list of B: its 4,000,000 members of lowest energy (cut
     3.289).
   - D: the same chains, energy ≤ 4.504 (nominal ball 4 · 10⁸): 10,000,008 members in 27 s, of which the 3,000,000 of
     lowest energy are kept (cut 3.884).

   Generation took 5,617 CPU seconds, 557 s of wall time from the start including 55 s for the moment estimates. The
   Gaussian estimate of the number of key-quadruples inside the lists [343,764; 14,000,000; 4,000,000; 3,000,000] is
   10^1.18 (15); it is a heuristic, not a rate. A Monte-Carlo test of it on lists of the same kind but smaller
   (A 340,468, B 2 · 10⁶, C 6 · 10⁵, D 5 · 10⁵; 7.5 · 10⁹ random quadruples) found the number of quadruples within
   squared whitened distance 4 of the target at 0.066 of the Gaussian count (0.07 to 0.08 for squared distances 2.5
   to 5). The lists are truncated balls, not Gaussian samples; an isotropic-ball calculation of the central density
   for these thresholds gives 0.029 of the Gaussian value. With the Monte-Carlo factor the expected number of
   key-quadruples inside the lists is about 1, with the isotropic one 0.4.
3. **k-tree join** (`search/cp_lowe.py`, `search/kt4.c`; [search/NOTES.md](../../search/NOTES.md) §20, §110, §111).
   `kt4.c` joins the four lists over split values of five merge classes (0, 7, 11, 14, 22 of the 25 autocorrelation
   classes). 6,000 split values were chosen from 143,967 candidates, ranked by expected solutions per unit of work
   (4 · 10⁶ sampled pairs per side); they hold 0.579 of the estimated solutions (8.8 of the 15 expected
   key-quadruples). The first 52 split values, in the order of that ranking, took 1323 s from the start of the join
   (about 25 s each) with 2.8 · 10⁹ pairs on the AB side (5.3 · 10⁷ on average) and 1.07 · 10¹¹ on the CD side
   (2.1 · 10⁹ on average, 7.7 · 10⁸ to 3.0 · 10⁹) and gave nothing. The 53rd gave one solution, verified, and the run
   stopped at its limit of 1 solution after 29 s more. Its split value, (−2, −6, −6, −2, −6) at the five merge
   classes, is the 53rd of the ranking: rebuilding the ranking from the saved lists reproduces the 52 logged split
   values and places it there. The first 53 carry 1.8 % of the estimated share, that is 0.27 of the Gaussian
   expectation (about 0.02 with the factor above).
4. **Earlier attempt.** A first run on the same class, with seed 1 and lists of about the same sizes (A 342,527, then
   14,000,000, 4,000,000 and 3,000,000), joined 17 of its 100 chosen split values in 810 s (other merge classes: 1, 6,
   7, 10, 11) and found nothing.

**Reproduction.**

| | |
|---|---|
| Code | commit `919e73d`: `search/cp_lowe.py`, `search/cyclo_product.py`, `search/workdir.py`, `search/kt4.c`, `verify/verify.py` (the files that ran) and `search/skew_certify.py` (certificate). Not in that commit: `lg_lowe.py`, `lg_lib.py`, `gen_run.py` and `lg_pt.c`, which ran from a working copy and are added as `search/listgen/lg_lowe.py`, `search/listgen/lg_lib.py`, `search/listgen/gen_run.py` and `search/listgen/lg_pt.c` ([search/listgen/README.md](../../search/listgen/README.md)). `lg_lowe.py` and `lg_pt.c` are byte for byte the files that ran. `lg_lib.py` and `gen_run.py` differ from them only in the one line each that put directories of the working copy on the machine of the run on `sys.path` (published md5: `lg_lib.py` 2662b4539b3cb78670bc6d9cae013633, `gen_run.py` d0f14e225bf572f7669c0101f4159335); `PYTHONPATH=search` supplies the modules |
| md5 of the files that ran | `lg_lowe.py` 415f7f57091a6e61ad4ba423c3f1f2a7; `lg_lib.py` 92721f7d60d12fe8b874b0479d632402; `gen_run.py` 72aa3b2851863e1bb05115fcc175e7bc; `lg_pt.c` 8e9fab9877257775e78d0dd42e455d18; `cp_lowe.py` 4e2f4c952a56d1d0942b5c1c1a7d4d70; `cyclo_product.py` 98be16555ffe42590a0cd2185add9f3e; `workdir.py` b845fde4e2f09a184b42ec79d7776f4b; `kt4.c` d22c49d986b9b3ba9ee6d3e62ded6fea; `verify.py` 970e5d9aba849fd38a7a4109a85a091d |
| Build | `cc -O3 -march=native -ffast-math -fopenmp search/listgen/lg_pt.c -o search/listgen/lg_pt -lm` (`lg_lowe.py` calls the binary next to it); `kt4` is compiled by `search/workdir.py` on first use (`-O3 -march=native -fopenmp -lm`). gcc 13.3.0 (Ubuntu 24.04), Python 3.12.3, NumPy 1.26.4. Rebuilt from these sources on the machine of the run, the binaries have the md5 of the ones that ran: `lg_pt` 84a3b7bed9e13368d64e6eca4fbd5435, `kt4` ffe1b35a47b8bc1fb8d2cbb5fda2e986 |
| Command | `PYTHONPATH=search OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 python3 search/listgen/lg_lowe.py 251 50 0 --LA 3e5 --LB 7e7 --LC 7e7 --LD 4e8 --nB 1.4e7 --nD 1e7 --kC 4e6 --kD 3e6 --share-bc --resume --keep-lists --secs-gen 1800 --threads 12 --T 0.42 --nsig 6000 --mcpairs 4000000 --maxpairs 3e8 --join-core-h 120 --K 1 --seed 2 --workdir W` with `HADAMARD_WORK` set to a work directory (the hit is appended to its `skew_hits.jsonl`) |
| Randomness and caps | `--seed 2` seeds the 4 · 10⁵ moment samples per block (`--samples`, the default), the chains (10 · seed + k = 20 + k for block k, and 1000003 times that plus t for thread t) and the 4 · 10⁶ pairs per side that rank the split values; `--LA/--LB/--LC/--LD` set the energy thresholds through the nominal ball sizes; `--nB`, `--nD` the number of distinct members the chains collect (block A: all below its threshold); `--kC`, `--kD` the members of lowest energy kept; `--secs-gen` a time cap per block (not reached); `--nsig` the number of split values; `--maxpairs` the cap on pairs per side of a split value; `--join-core-h 120` stops the join after 120 core-hours; `--K 1` stops it at the first solution. The lists are not reproducible bit for bit (the threads race into one hash table and the chains stop when the target is reached), so a rerun gives lists of the same sizes and thresholds but other members; the ranking of the split values is determined by the lists |
| Result | one verified solution after 1918 s (0.53 h) wall on 12 threads, 23,017 core-seconds (6.4 core-hours; 12 threads times the wall time): lists 557 s, ranking of the split values 7 s, join 1352 s; about 4.2 GB peak resident memory (periodic sample of the process group) |
| Certificate | `search/skew_certify.py` turns the logged hit into the files in this directory |

## Status

**New** as a skew-Hadamard matrix. An ordinary Hadamard matrix of order 1004 is known (Đoković, Golubitsky and
Kotsireas, J. Combin. Des. 22 (2014) 270–277); the skew property is what is new.

- The Đoković–Seberry census lists 251 among its 737 undecided odd v < 2500 for skew matrices, and not in its Paley,
  Szekeres, Szekeres–Whiteman or Spence series or among its sporadic values.
- Cati–Pasechnik Table 5 (v2, 30 August 2025; entries n(m) give the smallest m with a skew matrix of order 2ᵐn known,
  and m = 2 is not shown) lists 251(6): the smallest known skew order with odd part 251 was 2⁶ · 251 = 16064, from
  Paley I (16063 is prime and ≡ 3 (mod 4)). For 2ᵘ · 251 with u = 2, 3, 4, 5 the numbers 1003 = 17 · 59, 2007 = 3² ·
  223, 4015 = 5 · 11 · 73 and 8031 = 3 · 2677 are not prime powers. Their Table 3, the constructions of skew matrices
  of order 4v, has no entry for 251, while their Table 2 lists 251 for ordinary matrices (SDS). Koukouvinos–Stylianou
  (2008), Table 1 also gives 251(6) (there m(t) means that a skew matrix of order 2ᵗm exists and m(·) that none is
  known; odd m < 500 missing from the table have 4m known). Seberry–Yamada's tables (1992, Table 7.1; Utilitas Math.
  1991; 2020, Table 9.2) reach 251; we did not re-read their entries for it, but Cati–Pasechnik's Table 5 updates the
  last of them and Koukouvinos–Stylianou's Table 1 the Handbook's (1996 edition, Table 24.31). SageMath's development
  source (`hadamard_matrix.py`, fetched 2 November 2026) lists 1004 in `unknown_skew_hadamard` (the docstring list of
  the orders up to 1200 for which no skew matrix is known), and no branch of its `skew_hadamard_matrix` reaches 1004:
  Paley I, 2(q + 1) with q ≡ 5 (mod 8), complementary difference sets, Spence 1975, orthogonal designs, Williamson's
  lemma and the stored skew supplementary difference sets (no v = 251 among them) all fail, and 1004 is not a multiple
  of 8.
- The ordinary matrix of order 1004 is that of Đoković, Golubitsky and Kotsireas (arXiv:1301.3671), also a
  Goethals–Seidel matrix over ℤ₂₅₁ with the same multiplier group H = {1, 20, 113, 149, 219}. Its abstract lists
  skew-Hadamard matrices only of orders 4 · 213 and 4 · 631. It gives two supplementary difference sets
  4-{251; 125, 120, 115, 115; 224}, row sums (1, 11, 21, 21), by index sets of orbit representatives; we rebuilt both
  from the printed index sets (they satisfy λ = 224, and 0 lies in none of their blocks) and none of their eight
  blocks is of skew type, even after a translation: the blocks of 120, 115 and 115 elements have row sums 11, 21 and
  21, not ±1, and in the 125-element block (50 elements have their negative in the block too in the first, 60 in the
  second) at most 75 of the 125 pairs {c + t, c − t} have exactly one element in the block, for any centre c, where a
  skew block needs all 125. Our matrix is in the other row-sum class, (1, 9, 9, 29) against (1, 11, 21, 21) in absolute
  value, which translation, multipliers, complements, reversal and reordering of the blocks do not change, and its
  first block is skew, which none of theirs is; so it is not theirs. Egan, Kotsireas, London and Ó Catháin
  (arXiv:2610.27175, 19 October 2026) mention 1004 only as an ordinary Hadamard order found by Đoković and coauthors.
  The text of the skew row of their Table 3 (as extracted from the PDF) lists as known skew orders Paley's q + 1,
  doubling, Spence's family, "4n for odd n < 250", 276 and 292 (citing Đoković 2023), with 356 the smallest open
  order; since 356 = 4 · 89 is open, "4n for odd n < 250" cannot mean every odd n < 250 (Cati–Pasechnik's Table 5
  leaves 89, 101, 107, … , 249 open), and in any case 251 > 250. The Hadamard Matrix Explorer (math74.com/Hadamard)
  builds 1004 as `G(A(251))`, from supplementary difference sets; it tries `A` before the skew-type route `SK`, so
  this label does not show whether it holds a skew `SK(251)`, and we did not inspect its matrix file.
- The shape of this matrix, the Paley–Todd set together with a 3-{v; k₁, k₂, k₃; λ} difference family, λ = k₁ + k₂ +
  k₃ − (3v − 1)/4, is the class of Đoković and Kotsireas (Special Matrices 4 (2016) 317–325, arXiv:1511.08734), who
  list its parameter sets for v ≤ 131 (with an SDS where one is known) and construct it for v = 239, the first skew
  matrices of order 956. The other new skew order in their abstract, 1324 (v = 331), is outside this class: its
  parameter set (331; 165, 155, 155, 155; 299) has four blocks with a skew block that is not a difference set.
  (251; 121, 121, 111; 165) is a parameter set of that class (v = 251 is a prime ≡ 3 (mod 4), 125.5 > 121 ≥ 121 ≥ 111,
  λ(v − 1) = 41,250 = Σ k(k − 1)); their paper contains no family for v = 251.
- None of the known constructions (Seberry–Yamada's families SI–SXI, the census series, Cati–Pasechnik, SageMath) gives
  a skew matrix of order 1004:
  - Paley I: 1003 = 17 · 59 is not a prime power.
  - Szekeres, Whiteman, Spence and Momihara–Xiang 2(qᵉ + 1): 501 = 3 · 167 is not a prime power.
  - Wallis–Whiteman and Momihara–Xiang 4(qᵉ + 1): 250 = 2 · 5³ is not a prime power.
  - Spence 4(q² + q + 1): q² + q + 1 = 251 has no integer solution (15² + 15 + 1 = 241, 16² + 16 + 1 = 273).
  - Turyn (h − 1)ᵘ + 1 from a skew matrix of order h, u odd: 1003 is not a perfect power, so only u = 1 (the order
    itself) remains.
  - Doubling needs a skew matrix of order 502, which is not a multiple of 4.
  - Skew × amicable products need 1004 = mn with m and n Hadamard orders, and the only such factorisation is
    1 · 1004 (2 · 502 and 4 · 251 have a factor 502 or 251). Williamson's h(pᵏ + 1), with pᵏ + 1 ≡ 0 (mod 4), needs
    pᵏ + 1 = 4, that is a skew matrix of the odd order h = 251, which cannot exist (or pᵏ + 1 = 1004, and 1003 is not
    a prime power).
  - Seberry's mn(n − 1) is possible only with n = 2, which needs an OD(1004; 1, 502, 501).
  - The 3(q + 1) theorem needs 3 | 1004. Blatt–Szekeres needs complementary difference sets in a group of order
    2 · 250 + 1 = 501 = 3 · 167; the known families need 501 or 2 · 501 + 1 = 1003 = 17 · 59 to be a prime power, and
    neither is.
  - Yamada's (|t| + 1)(q + 1), with q = s² + 4t² ≡ 5 (mod 8) a prime power, has q + 1 ≡ 6 (mod 8), so it needs
    q + 1 = 502 and q = 501, not a prime power. Seberry's 2ᵗq with q = s² + 4r² ≡ 5 (mod 8) a prime power needs q = 251,
    which is ≡ 3 (mod 8). (Families SVIII and SX of the list of known skew orders in Seberry–Yamada 1992.)
  - Good matrices and skew Williamson-type matrices of order 251 are not recorded.

The sources searched are listed in [NOVELTY.md](../../NOVELTY.md). Doubling (Seberry) gives skew-Hadamard matrices of
orders 2008, 4016 and 8032, which are new as far as we know; 16064 was already known.

## Files

| File | Contents |
|---|---|
| `certificate.json` | the four sequences as orbit representatives, H, the row sums |
| `sequences.txt` | the four sequences, one per line |
| `H1004.csv.gz` | the matrix |

## Verification

```
python3 verify/verify.py skew/1004
gzip -dc skew/1004/H1004.csv.gz | ./check_matrix
```

SHA-256 of `H1004.csv`: `e4a7651a62ea59e1968a178fa5c72cff13ac8417e76d13899e5ea99130092edb`
