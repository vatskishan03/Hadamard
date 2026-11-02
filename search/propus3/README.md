# propus3: an exact three-list search for propus quadruples under a multiplier group (symmetric orders 4v)

Programs and logs behind [search/NOTES.md](../NOTES.md) §114: exhaustive searches, all with no solution, for symmetric
Hadamard matrices of order 4v from the propus array of §17, in 16 multiplier structures for v = 119, 213, 239, 247, 259,
267, 287, 343 and 397. Each of the nine values is listed as undecided for symmetric matrices by the Đoković–Seberry
census of 2025-05-28 (`search/data/symmetric_census_undecided.txt`). The engine `p3` (C, pthreads) joins three lists of
orbit-sign vectors, one for each block of the quadruple, with lossless reductions (§2) that cut the list of B by a factor
of 750 to 16,000 and the list of D by 2.6 to 19 (recomputed from the logs; for `sym343_S38` by 860 and 160); in the
largest structures the list of D still holds 2 · 10⁸ to 3 · 10⁹ members. No certificate depends on these
programs: `verify/` checks every matrix on its own, and none was found here.

Labels: CHECKED = exact computation or a line of a log; HEURISTIC = estimate, nothing depends on it for correctness.

## 1. The problem

A *propus quadruple* is (A, B, B, D), four ±1 sequences on an abelian group G of order v with A(−g) = A(g) and

    PAF_A + 2 PAF_B + PAF_D = 4v δ.

By the second array of [search/NOTES.md](../NOTES.md) §17 (the form of [symmetric/532](../../symmetric/532) and
[symmetric/876](../../symmetric/876); `verify.py`, construction `propus`) it gives a symmetric Hadamard matrix of order
4v. The row sums (r_A, r_B, r_D) are odd with r_A² + 2r_B² + r_D² = 4v: these are the *row patterns*, and a pattern is
*feasible* when each block can reach its row sum with the orbit sizes of its structure, that is, when the three pools of
sign vectors are not empty. A group H of multiplier automorphisms of G with −1 ∉ H is imposed: A is constant on the
orbits of ⟨H, −1⟩ (so it is symmetric), B and D on the orbits of H (D is free). Then each autocorrelation is constant on
the C *classes* (the orbits of ⟨H, −1⟩ on G ∖ {0}), and the problem is C equations

    F_A(c) + 2 F_B(c) + F_D(c) = 0,        F_c(y) = const_c + Σ_{i<j} T_c[i][j] y_i y_j,

in the class values F_c of three sign vectors y ∈ {±1}^O (O_A orbits of ⟨H, −1⟩ for A and O_B orbits of H for B and D,
each with the orbit {0}); the T tensors are counts of pairs of orbits at a shift. The structures below have O_A from 14
to 20 and O_B from 23 to 39.

## 2. The engine

`p3.c` (the header of the file says the same):

1. **Lists.** A: every sign vector on the O_A orbits with row sum r_A. B and D: every sign vector on the O_B orbits with
   row sum r_B and r_D, by a depth-first search over the orbit signs with incremental class values (int16 vectors) and an
   exact reachability test of the row sum, in two passes (count, then fill) so that the lists are built without
   temporary copies; the tasks are search prefixes, handed to the threads. The number of leaves of the search is the
   pool, which the logs print as `COUNT`.
2. **Reductions** at the leaves, all lossless:
   - *canonical B*: the map x ↦ tx, for t in a unit group modulo H (|Mult| of them; the multiplier group of the structure
     modulo H), takes a solution to a solution, so B can be taken least in its |Mult| images (factor |Mult|);
   - *canonical D*: reversing D alone, x(g) ↦ x(−g), leaves PAF_D unchanged (factor 2);
   - *power-spectrum bounds*: |Â(ω)|² + 2|B̂(ω)|² + |D̂(ω)|² = 4v at every character ω, so |Â|² ≤ 4v, |B̂|² ≤ 2v and
     |D̂|² ≤ 4v; the left sides are computed from the class values through the C × C matrix of Gauss periods, with a
     slack of 0.5 (the final test is exact integer arithmetic).
   An optional low-energy cut exists in `p3.c` and `p3lowe.py`; no run reported here used it.
3. **Join.** Every member carries its class-value vector F and the 64-bit hash h = Σ_c R_c F_c (R_c random odd constants
   written into the problem file; the seed is v). A solution needs h_D = −h_A − 2h_B and, on NS = 3 *split* classes S,
   F_D(S) = −F_A(S) − 2F_B(S). Since F ≡ v (mod 4) for every block, (F − v mod 4)/4 is an exact lattice coordinate, and D
   is grouped by its coordinates on S (a *σ group*). For each σ group, a table of the hashes of its members of D is built;
   for each a in A the group of B with 2F_B(S) = −σ − F_A(S) is found by index arithmetic, and every member of that group
   probes the table. The number of probes is `pairs` in the logs. Every table hit (`matches`) is re-verified on all C
   classes from the three orbit masks; `badmatches` counts hash false positives. The join is exact: no solution is lost by
   the grouping or the hashing, and the lossless reductions above are the only filters.
4. **Work.** The unit of work is a batch of σ groups: threads, `-shard i/n` for several processes or machines, and a
   checkpoint file (fingerprinted by the problem and the options) with `-resume`.

A hit would be turned by `p3cert.py` into the propus matrix of order 4v and a certificate, which `verify/verify.py`
accepts. The reductions were also checked independently of the engine: the orbits of the full symmetry group of the
canonical solutions reproduce the full solution set in all 336 + 2 validation structures (§6).

**What an empty result covers.** For each structure of §3 and each of its feasible row patterns: there is no (A, B, B, D)
with A constant on the orbits of ⟨H, −1⟩, B and D constant on the orbits of H, and PAF_A + 2PAF_B + PAF_D = 4vδ. It does
not cover other multiplier groups of the same v (the list of §3 is the list of the jobs that finished, not a claim that
every subgroup of the unit group was searched), a group H that contains −1 (then B and D are symmetric; `p3struct.py`
refuses it), the variant of the array with A and D both symmetric (§17), or solutions without a multiplier symmetry. For
v = 343 it is the group ℤ₇³ with a subgroup of the Singer cycle of GF(7³)*, not ℤ₃₄₃. It says nothing about symmetric
Hadamard matrices of order 4v in general.

## 3. Results

Every run below finished with `hits=0 matches=0 badmatches=0`. Columns: H, the generators of the multiplier group
(of ℤ_v*; for v = 343 the exponent 38 of a Singer cycle S of GF(7³)*, acting on ℤ₇³; |H| was recomputed from the
generators); C, O_A, O_B as in §1; row patterns (r_A, r_B, r_D), all the feasible ones (recomputed: the patterns run are
exactly those whose three pools are not empty, and the pools printed by `p3` equal the exact counts of a dynamic program
over the orbit sizes; for `sym343_S38`, whose log does not print the pools, the patterns and the pools were recomputed
with `p3g.py`); pairs, the probes of the join, summed over the row patterns; time, the sum over the row patterns
of `total_seconds` of `p3` (lists and join, wall clock, on a machine shared with other jobs).

| Structure | v | H | \|H\| | C | O_A | O_B | row patterns | pairs | time (s) | threads |
|---|---:|---|---:|---:|---:|---:|---|---:|---:|---:|
| `sp_v119_H4_g13` | 119 | 13 | 4 | 19 | 20 | 32 | (1, 15, 5); (3, 11, 15); (5, 9, 17); (5, 15, 1); (13, 3, 17); (15, 11, 3); (17, 3, 13); (17, 9, 5) | 314,815,197 | 202.5 | 4 |
| `sp_v213_H10_g5` | 213 | 5 | 10 | 15 | 16 | 30 | (7, 19, 9); (9, 11, 23); (9, 19, 7); (11, 1, 27); (11, 19, 3); (13, 11, 21); (27, 1, 11); (29, 1, 3) | 5,821,919 | 47.5 | 4 |
| `p239` (by hand) | 239 | 10 | 7 | 17 | 18 | 35 | (15, 1, 27) | 6,014,289,982 | 59.6 | 32 |
| `sp_v247_H9_g35` | 247 | 35 | 9 | 15 | 16 | 31 | (5, 1, 31); (5, 13, 25); (7, 7, 29); (7, 17, 19); (17, 13, 19); (17, 17, 11); (19, 1, 25); (19, 7, 23); (19, 13, 17); (19, 17, 7); (29, 7, 7); (31, 1, 5) | 69,169,107 | 160.5 | 4 |
| `sp_v247_H9_g9` | 247 | 9 | 9 | 15 | 16 | 31 | the same twelve | 60,199,500 | 166.0 | 4 |
| `sp_v247_H9_g68x87` | 247 | 68, 87 | 9 | 17 | 18 | 35 | the same twelve | 1,603,501,642 | 994.1 | 8 |
| `sp_v259_H9_g9` | 259 | 9 | 9 | 15 | 16 | 31 | (5, 5, 31); (5, 19, 17); (31, 5, 5) | 32,013,577 | 38.1 | 4 |
| `sp_v259_H9_g44` | 259 | 44 | 9 | 15 | 16 | 31 | (5, 5, 31); (5, 19, 17); (31, 5, 5) | 35,990,599 | 39.5 | 4 |
| `sp_v259_H12_g6x11` | 259 | 6, 11 | 12 | 13 | 14 | 23 | (5, 5, 31); (5, 19, 17); (17, 19, 5); (31, 5, 5) | 194,815 | 1.3 | 4 |
| `sp_v259_H12_g6x82` | 259 | 6, 82 | 12 | 13 | 14 | 23 | (5, 5, 31); (5, 19, 17); (17, 19, 5); (31, 5, 5) | 207,673 | 1.3 | 4 |
| `sp_v259_H12_g11x23` | 259 | 11, 23 | 12 | 13 | 14 | 24 | (5, 5, 31); (5, 19, 17); (17, 17, 13); (17, 19, 5); (31, 5, 5) | 109,452 | 1.7 | 4 |
| `sp_v259_H12_g36x43x60` | 259 | 36, 43, 60 | 12 | 13 | 14 | 24 | (5, 5, 31); (5, 19, 17); (17, 17, 13); (17, 19, 5); (31, 5, 5) | 104,263 | 1.7 | 4 |
| `sp_v267_H11_g4` | 267 | 4 | 11 | 13 | 14 | 27 | (1, 23, 3); (3, 23, 1) | 1,774,811 | 3.5 | 4 |
| `sp_v287_H15_g16` | 287 | 16 | 15 | 13 | 14 | 27 | (5, 9, 31); (13, 3, 31); (13, 15, 23); (25, 9, 19); (27, 13, 9); (33, 5, 3) | 43,001 | 5.2 | 4 |
| `sym343_S38` | 343 | S³⁸ | 9 | 19 | 20 | 39 | (17, 19, 19); (19, 19, 17) | 1,394,808,195,608 | 6,406.8 | 8 |
| `sp_v397_H11_g16` | 397 | 16 | 11 | 18 | 19 | 37 | (1, 23, 23) | 31,513,110,939 | 585.4 | 12 |

16 structures and 88 (structure, row pattern) cases; 1.43 · 10¹² pairs in all, 1.39 · 10¹² of them in `sym343_S38`; 8,715 s
of `p3` time summed over the structures (recomputed from the logs by `summarize_logs.py`). `results.tsv` has one line for
each case: the pools (`COUNT` lines), the list sizes A′, B′, D′ after the filters (`LISTS`), the pairs, the join and total
seconds, the threads and the counts of hits. The group name `sp_v<v>_H<|H|>_g<g₁>x<g₂>…` is the name of the scheduler job
and the file name of its log; its generators are g₁, g₂, …

- **v = 239** (`p239`). Run by hand, `THREADS=32 NS=3 sh run239.sh`, on 1 November 2026 (14:21 UTC). Pools: A 24,310, B
  2,333,606,220, D 1,855,967,520; lists after the filters: A′ 11,118, B′ 597,980, D′ 721,194,562; 14,724 σ groups;
  6,014,289,982 pairs in 2.09 s, 59.65 s in all, 1,261.6 s of user time, 14.1 GB peak resident memory (`/usr/bin/time`).
  The same problem without the canonical-form reductions (the power-spectrum filters on; `prob239_noreduce.txt`, which
  `p3problem.py --noreduce` writes) has B′ 20,331,290 and D′ 1,442,389,124 and 204,206,380,194 pairs (127.1 s, 28.0 GB):
  also no hit. It agrees with the reduced run; it cannot by itself show that a solution is not lost by a reduction (§2
  and §6 do that). The repository's estimator (`est239.py`, `search/propus_rank.estimate`, two seeds) gives 10^1.38 = 24
  solutions for this structure and pattern, 0.35 expected orbits of the symmetry group of order |Mult| · 2 = 68
  (HEURISTIC; on the hits [symmetric/532](../../symmetric/532) and [symmetric/876](../../symmetric/876) it was low by
  10^1.4 to 10^1.9, §28 and §40). The other nontrivial subgroups of ℤ₂₃₉* of odd order, ⟨6⟩ of order 17 and ⟨2⟩ of order
  119, have no feasible row pattern (`p3sweep.py 239 6` and `239 2` print 0 hits and no pattern; `logs/validation.log`).
  ℤ₂₃₉* is cyclic of order 238 and every subgroup of even order contains −1, so for v = 239 the only multiplier group with
  −1 ∉ H that is not searched is H = 1.
- **v = 397.** One pattern. The list of D has 3,298,649,676 members (pool 8,597,496,600), B′ 1,787,139, A′ 19,226;
  31,513,110,939 pairs; 585.4 s on 12 threads; peak resident memory 52.2 GB (the scheduler's periodic sample).
- **v = 343** (`sym343_S38`, `group/`). G = ℤ₇³, H = ⟨S³⁸⟩ of order 9 where S is a Singer cycle of GF(7³)* (order 342); the
  group of the problem is the group of semilinear maps x ↦ aσ(x) (order 1026), H is normal in it, the quotient has order
  114 (`|mult|=114`) and the symmetry group has order 228 (`|Gamma_amb|=228`); the only fixed point of H is 0
  (`|fix|=1`). The pools, recomputed with `p3g.py`, are A 92,378 and B and D 33,578,000,610 for both patterns, and these
  two are the only feasible patterns. `p3g_job.py --dmult` spends the quotient on D instead of B (`p3 -dperms`: D minimal
  under reversal and the 114 maps), B is reversal-canonical only. The log gives for (17, 19, 19) the lists A′ 20,691, B′ 39,068,598, D′
  206,610,930, 691,740,550,512 pairs, a join of 350.99 s and 3260.41 s in all; for (19, 19, 17) A′ 20,691, B′ 39,068,598,
  D′ 213,811,381, 703,067,645,096 pairs, a join of 365.45 s and 3146.37 s in all. The job was started twice: the first
  run (03:54 to 04:25 UTC) was ended by a restart of the scheduler of the machine before the join of its first pattern
  had begun (inferred: the hits file of that pattern was created at 05:15 UTC, in the second run), and the second run
  did everything; 8,262 s of wall time on 8 threads in all, peak resident memory 4.25 GB (the scheduler's periodic
  sample).
- **Planned and not run.** The scheduler's state records `sp_v341_H10_g2`, `sp_v341_H10_g32x35`, `sp_v341_H10_g32x47`,
  `sp_v341_H10_g32x70` as pending, and `sp_v259_H9_g100x121`, `sp_v305_H10_g9`, `sp_v1123_H33_g7x33` as cancelled; they
  produced no log and no result. Their scripts are not in `jobs/`.

## 4. Environment

The runs of §3 were made on the 32-core machine of [search/NOTES.md](../NOTES.md) §110 and §111: AMD EPYC 9R14, 32 cores,
246 GB, shared with other jobs (4, 8, 12 or 32 threads used), Ubuntu 24.04, gcc 13.3.0, Python 3.12.3, NumPy 1.26.4. The
job runs (`sp_v…`) were made on 1 November 2026 between 16:41 and 17:12 UTC, `p239` at 14:21 UTC, `sym343_S38` on 2
November between 03:54 and 06:13 UTC. The binaries that ran: `p3` (for `p239` and `sp_v…`) has md5
`8d88a11272e6ba0efcd8b5128e1b0fd0`, built on the machine on 1 November at 14:20 UTC; `group/p3` has md5
`7a8b34c3972b0783510fcee9fb9cfd1e`. The compile command that produced them was not logged; the one in the header of
`p3.c`, in `run239.sh` and in `p3run.build_binary` is `cc -O3 -march=native -ffast-math -pthread p3.c -o p3 -lm`.

## 5. Files

The files are byte for byte the ones that ran (md5 below), except 20 files in which lines that named directories of the
machine of the run, or of the Mac, were edited before publication (§7). The md5 in the table and in `MD5SUMS` are those
of the files as they ran, except `logs/validation.log`, which is listed as published (see below); `MD5SUMS.edited` has
the md5 of the 20 edited files as published.

| File | md5 | Contents |
|---|---|---|
| `p3.c` | `b3a87baa0ee8a7f489ea341a43be89b9` | the engine; usage: `./p3 PROBLEM [-t threads] [-ns k] [-shard i/n] [-ckpt file] [-resume] [-hits file] [-count] [-listonly] [-countlists] [-S c0,c1,c2] [-bs batch]` (and `-alist`, `-blist`, `-dlist`, `-slice` for tests) |
| `p3struct.py` | `4097142b945972cc5ad41db07ef7eedc` | the structure from first principles: orbits, classes, T tensors, the power-spectrum matrix, the pools, the feasible row patterns, uniform sampling, and `selfcheck` against direct FFT and autocorrelations; imports nothing of the repository |
| `p3run.py` | `ea88e81a9878cdc3d4961f7269407d8c` | the problem file writer, the runner, the exact re-check of a hit on the sequences, the orbits of the symmetry group |
| `p3problem.py` | `3c654ed8c449474483daf51cf660a8e4` | `p3problem.py V GEN --rows a,b,d --out FILE [--noreduce] [--nopsd]` |
| `p3sweep.py` | `bd90de880dc8e7b79269d945707dc93b` | all feasible row patterns of one structure, reduced and exhaustive |
| `p3plan.py` | `2d9a9eb6551cb1ba34a78a9624037cfd` | cost preview: pools, sampled power-spectrum retention, list sizes, pairs |
| `p3lowe.py` | `c60dd008d9474e5d41a52182f8408aa7` | the optional low-energy cut (heuristic, not used in §3) |
| `p3cert.py` | `e36f593738fc0f104e6414b87b7d5ca7` | a hit to the matrix and the certificate (`--matrix` writes `H<4v>.csv.gz`) |
| `p3test.py`, `p3validate.py`, `p3cert_test.py`, `p3plant.py`, `p3scale_test.py`, `p3resume_test.py`, `psd_retention_check.py`, `est239.py`, `probebench.c` | `9abc0f71e945934b51629aebf313214a`, `2e368062f3320bd695d7a627895edc4b`, `6dffa0a85a178268ba1ae05e6f747772`, `ff645919cd93a1ea2ac69fb696d7646b`, `31da8fb6ff63fa78603c88ff0d15218a`, `2a70bcb47cd4110aed4b104deaba939f`, `0f1d3ee8264a4fde6c3a720d98f29ded`, `62fb34a53fe48663920f15380f4ab139`, `5c22f7c7f7850b41553510f807ecf309` | the tests of §6, the Gaussian estimate of v = 239 and a microbenchmark of the join kernel |
| `run239.sh`, `run_all_tests.sh` | `321d034ca5b0cc8aacf3d215107d9097`, `f018833e2128ddda6a4feffdc7ab46b2` | the production run of v = 239 (resumable) and the whole test suite |
| `prob239.txt`, `prob239_noreduce.txt` | `1a035ae0cba73361ee23142953aaf0ee`, `38190bfa30ee2e0aef1815507aa0c0a9` | the problem files of the two runs of v = 239 (sha-256 of the first: `a98f371583418f15bb5e4cb2b0e7954988c31c92c6eb49c821442f9fa29b9739`); `p3problem.py` regenerates both byte for byte (CHECKED here) |
| `jobs/mkprob.py` | `b3363a84a308ed52de4093f652ce1182` | `mkprob.py V G1,G2,… --rows a,b,d [--rows …] --outdir DIR --tag TAG`: the problem files of the job scripts (any number of generators; same code path as `p3problem.py`, seed v) |
| `jobs/sp_*.sh` (14) | see `MD5SUMS` | the wrappers of the scheduler jobs: they write the problem file of each row pattern, run `./p3 … -t THREADS -ns 3 -ckpt … -resume -hits …`, record the finished pattern, and delete the problem files |
| `group/p3.c` | `5fcd023cce5e7ce776123f184e845427` | `p3.c` with two added options, `-dperms FILE` (D canonical under reversal and the listed permutations) and `-bsave FILE` (write the filtered B list) |
| `group/p3g.py`, `group/p3g_job.py` | `76b8ece4aa26218ffd3e0d6750823cc3`, `d14d52b4be9f2cec8b7026ef5b157df2` | the front end for a finite abelian group G = ∏ ℤ_p^a with the cyclic group H generated by an automorphism φ of the group ΓL(1, p^a) per prime (Singer cycle of GF(p^a)* times Frobenius): the orbits and class tensors come from `search/spence_lib.py`, the group PAF and the matrix from `search/propus_group.py`; `p3g_job.py` runs all feasible row patterns of one structure resumably and writes `HIT_<id>.json` only for a hit that passes the exact group PAF check and the matrix check |
| `group/p3g_validate.py`, `group/p3g_validate_dm.py`, `group/test_bcache.py` | `e3d09e690a1769b4282106b324853eb5`, `c3b36b6f8f700f13fbc80ec2a9cb7b73`, `61c496b9c75ee818b1f0c72abb169dcc` | the validation of the group engine against a brute force, of the `--dmult` reduction, and of the cached B list of shards |
| `group/MANIFEST.md5` | `749c8c120a1f2c378487c1a7597afd4e` | the md5 of the six files that ran for `sym343_S38` (`p3`, `p3.c`, `p3g.py`, `p3g_job.py`, `p3run.py`, `p3struct.py`; the last two are the files above) |
| `logs/` | | the logs of §3, byte for byte (`sp_v….log` is `run_01.log` of the job; `sym343_S38_run_0[12].log`; `v239_run239_0_1.log`, `v239_run239.out`, `v239_run239_noreduce.out`), and `validation.log` (see below) |
| `group/logs/` | | the validation logs of the group engine |
| `results.tsv`, `structures.tsv`, `summarize_logs.py` | | derived: the tables of §3 recomputed from `logs/` (`python3 summarize_logs.py logs > results.tsv`; `--structures`) |
| `MD5SUMS`, `MD5SUMS.edited` | | the md5 of every file in this directory except `README.md` and the two sums files, as the files ran (70 files); the md5 of the 20 edited files of §7 as published |

`logs/validation.log` is the output of `run_all_tests.sh` on an Apple M4 Mac (arm64, Python 3.11.14, NumPy 1.26.4,
Apple clang 21) on 1 November 2026, with the host name and the two paths of the directory of that machine replaced by
`<host>` and `out/` (md5 of the file as it was: `bb441aa0cf1ad009dc3055f0eb86f3e4`). It is the only logged run of the
suite; the binary of the machine of §4 was not run through the suite in a log that we found.

## 6. Validation

From `logs/validation.log` (CHECKED, arm64 Mac) and the logs of `group/logs/` (the machine of §4):

- **A reference solver** (`p3test.py`): for all 336 (v, H, row pattern) cases with v ≤ 75 and pools up to 40,000 (7,060
  solutions in all; prime and composite v, |H| from 2 to 6) the engine returns exactly the set of an independent numpy
  dictionary join, unreduced with the power-spectrum filters on and off, reduced (exactly the canonical members), with
  NS = 1, 2, 3, with 1 and 3 threads and with 1 and 3 shards; the orbits of the symmetry group of the canonical solutions
  equal the full set. No failure.
- **The two known matrices** (`p3validate.py`): v = 133, H = ⟨26⟩, rows (11, 5, 19) ([symmetric/532](../../symmetric/532)):
  72 solutions, 2 reduced hits, symmetry group of order 36, the repository's solution found, the same 72 with the filters
  off (2.1 · 10⁹ pairs); v = 219, H = ⟨4⟩, rows (3, 17, 17) ([symmetric/876](../../symmetric/876)): 288 solutions, 9 reduced
  hits, order 288. The counts are those of NOTES §23, §28, §30. `p3cert_test.py`: the certificates and matrices rebuilt from the
  engine's hits are byte for byte those of the repository, and `verify/verify.py` accepts them.
- **At the dimensions of v = 239**: a planted solution (`p3plant.py`: a random triple from the true pools defines a nonzero
  target; the join returns exactly it for NS = 1, 2, 3 and 3 shards, and nothing for a perturbed target); the filters and the
  leaf counts of the production configuration against numpy and a dynamic program (`p3scale_test.py`); checkpoint and
  resume (`p3resume_test.py`: killed after 1,464 of 4,129 batches, the resumed run completes, the union of the hits is the
  uninterrupted run's 72, and resuming with other options is refused).
- **A control from the literature**: v = 67, H = ⟨29⟩ of order 3 (the structure of Balonin and Đoković): 114 reduced hits
  in 0.3 s (`p3sweep.py 67 29`).
- **The group engine** (`group/logs/`, on the machine of §4): `p3g_validate.py` against a brute force on small structures
  over ℤ_p^a and mixed groups: 545 cases with 2,286 solutions, no failure (`validate2.log`; `validate1.log` is an earlier
  scan of 115 cases, all `OK`, without a summary line); the `--dmult` reduction: 31 cases, 2,126 solutions, all correct
  (`validate_dm.log`); the cached B list of shards: 4 cases (`bcache_test.log`). These logs were written before the files
  were frozen for the run of `sym343_S38`; the files of the validated working copy have the md5 of the frozen ones.
- **Re-run here** (CHECKED, 2 November 2026, Apple M4 Mac, Python 3.11.14, NumPy 1.26.4, Apple clang, the files as
  published after the edits of §7): `p3.c` compiles without warnings under `-Wall -Wextra`; `p3validate.py` ALL OK (72 and
  288 solutions); `p3test.py --vmax 75 --maxpool 40000`: 336 cases, 7,060 reference solutions, failures 0;
  `p3cert_test.py` ALL OK; `p3sweep.py 67 29`: 114 reduced hits; in `group/`, `test_bcache.py` ALL OK and
  `p3g_job.py --id t --v 25 --hgens 8,0` finds and verifies 9 solutions, the count of `validate_dm.log` for v = 25, H = ⟨S⁸⟩,
  rows (1, 7, 1).
- **Patterns and pools** (CHECKED): for every structure of §3 the feasible row patterns and the pools were recomputed
  (`p3struct.py`; `p3g.py` for `sym343_S38`) and equal the patterns run and, where the log has them, the `COUNT` lines.

## 7. Paths

Lines of the files that ran put directories of the machine of the run, or of the Mac, into the code. They were edited
before publication, so that the code finds the repository from the location of this directory, `search/propus3`; only
these lines changed (CHECKED: the 20 edited files have the md5 of the last column, which `MD5SUMS.edited` repeats). The
table gives the line as it ran, the edit, and the md5 of the file as it ran and as published.

| File | Line that ran | Edit | md5 as ran | md5 as published |
|---|---|---|---|---|
| `p3validate.py` | `REPO = Path(…)`, the repository directory of the machine of the run | `REPO = HERE.parent.parent` | `2e368062f3320bd695d7a627895edc4b` | `16d41dc23619928bd4b55211760b6463` |
| `p3cert_test.py` | the same | the same | `6dffa0a85a178268ba1ae05e6f747772` | `90f7512dab23b18cf99b80b0e9328dfa` |
| `est239.py` | `sys.path[:0] = […]`, the `search/` and `verify/` of that repository | the `search/` and `verify/` of the repository, from `__file__` | `62fb34a53fe48663920f15380f4ab139` | `c8e7f7456ab2220315ab2598af9a124f` |
| `run_all_tests.sh` | `PY=${PY:-…}`, the python3.11 of the Mac | `PY=${PY:-python3}` (the benchmark at its end writes `${TMPDIR:-/tmp}/p3_prob239_bench.txt`) | `f018833e2128ddda6a4feffdc7ab46b2` | `94aedbd7a22cd87a7a615f75f4670bd3` |
| `jobs/mkprob.py` | `sys.path.insert(0, …)`, the directory of `p3` on the machine of the run | the directory above `jobs/` | `b3363a84a308ed52de4093f652ce1182` | `73713b0ba2a1e1dc0005aeb1aa1ec1d1` |
| `jobs/sp_*.sh` (14) | `P3=` the directory of `p3` and `W=` the planner's work directory, both on the machine of the run; `"$W/mkprob.py"`; and the path in the comment on `p3cert.py` | `P3=$(cd "$(dirname "$0")/.." && pwd)`, `W=${W:-$PWD/work}` (the work directory of the checkpoints, hits and problem files), `"$P3/jobs/mkprob.py"`, `$P3/p3cert.py` | `MD5SUMS` | `MD5SUMS.edited` |
| `group/p3g.py` | `REPO = os.path.expanduser(…)`, the `search/` of that repository, and a default `HADAMARD_WORK` in the planner's directory | `REPO` = the `search/` above `group/`; the default of `HADAMARD_WORK` removed (`workdir.py` has its own) | `76b8ece4aa26218ffd3e0d6750823cc3` | `ec74acd040b798d3eb57b001b18fd0c9` |

The job scripts as published set `P3` and `W` relative to this directory, where those that ran named the directory of `p3`
and the planner's work directory on the machine of the run; the commands they run are the same.
`group/p3g_job.py` and the other files have no such line.

## 8. Usage

```
cc -O3 -march=native -ffast-math -pthread p3.c -o p3 -lm
python3 p3problem.py 239 10 --rows 15,1,27 --out prob239.txt        # or use the prob239.txt of this directory
THREADS=32 NS=3 sh run239.sh                                         # the run of §3; hits in work/hits239_0_1.txt
python3 p3cert.py --v 239 --gens 10 --hits work/hits239_0_1.txt --out DIR --matrix     # only if a hit exists
python3 ../../verify/verify.py DIR/symmetric/956
```

A scheduler job (any generators): `sh jobs/sp_v119_H4_g13.sh` (with `THREADS` and `W` as
environment variables; one problem file per row pattern, resumable). `python3 p3sweep.py V GEN` runs every feasible row
pattern of a cyclic structure and prints the pools, the lists and the number of reduced hits (for example `p3sweep.py 67 29`
in 0.3 s, `p3sweep.py 133 26` for the 72 solutions of [symmetric/532](../../symmetric/532)). Before a long run, `python3
p3validate.py` takes seconds; `sh run_all_tests.sh` the suite of §6 in about a minute on the Mac of §6.

The group engine (`PYTHONPATH` must contain this directory, for `p3run.py` and `p3struct.py`, and `search/`):

```
cd group && cc -O3 -march=native -ffast-math -pthread p3.c -o p3 -lm
PYTHONPATH=..:../.. THREADS=8 python3 p3g_job.py --id sym343_S38 --v 343 --hgens 38,0 --dmult --ns 3 --workdir WORK
```

`--hgens` is a `;`-separated list of generators of H, each the exponents of the Singer cycle of each prime of v and then
of the Frobenius of each prime (`38,0`: S³⁸ and no Frobenius). The job writes `HIT_<id>.json` only for a verified hit.

## 9. Not determined

- The compile command of the two binaries (§4).
- Whether the structures of §3 were already covered by the compressed propus sweep of [search/NOTES.md](../NOTES.md) §28
  and §40 (composite census values up to 435, orbits 1 to 36; §40 says its outcome is not recorded): the structures of
  119, 213, 247, 259, 267 and 287 have O_B of 23 to 35.
- How the structures were chosen: the planning files of the machine rank candidate multiplier groups by an estimate and a
  list size; we did not use them, and §3 does not claim that every multiplier group of these v was searched.
- A logged validation run on the machine of §4 (the only logged run of `run_all_tests.sh` is the Mac's).
- The first run of `sym343_S38`: `run_01.log` has only the header line, and the run left no other record than the
  scheduler's.
