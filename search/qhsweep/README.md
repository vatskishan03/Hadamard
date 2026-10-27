# qhsweep: exact sweeps of the Szöllősi arrays of order 2N with a pruned enumerator, a filtered join and exact screens

Programs behind [search/NOTES.md](../NOTES.md) §105 (QH(2N) for the eighteen Table 4 primes with N > 1000). Szöllősi's
array (`search/ssxx.py`, arXiv:2603.09572, Theorem 4) and its bordered form (`search/bsz.py`,
[complex-hadamard/BORDERED.md](../../complex-hadamard/BORDERED.md)) give a quaternary Hadamard matrix of order 2N from
four ±1 sequences constant on the orbits of a multiplier group H; for a prime n = 2N + 1 ≡ 3 (mod 4) that matrix is the
input of the complex Miyamoto construction (`verify.complex_miyamoto`,
[hadamard/COMPLEX_MIYAMOTO.md](../../hadamard/COMPLEX_MIYAMOTO.md)), a quaternary Hadamard matrix of order 4n and,
doubled, a Hadamard matrix of order 8n. `ssxx_sweep.py --exhaustive` enumerates every block by Gray code (2^O steps,
`--maxo 40`) and joins with `join4x.c`. The programs here replace both ends and add exact screens, so that a structure
with at most 64 orbits is decided whenever its lists fit in memory or the screen kills it, and a job with one list too
long for memory is joined slice by slice. Python 3 with NumPy and a C compiler with OpenMP (`enum_dfs.c` and `join4y.c`
are built on first use into `$HADAMARD_WORK/qhsweep/bin`). Lists, logs and hits are written below
`$HADAMARD_WORK/qhsweep` ([search/workdir.py](../workdir.py)), nothing inside the repository.

Labels: CHECKED = exact computation by `validate.py` in this directory.

Written in round 52 by two sweep runs (N = 1019, 1103, 1169, 1199, 1211, 1271, 1343, 1349, 1469 and N = 1125, 1143,
1185, 1335, 1341, 1439, 1443, 1451, 1499), which share the enumerator, the join, `z2fast.py` and the scheduler; the
prefix slices and `sliced.py` come from the first run, the joint quotient screen, the Williamson-type region and the
screen-only pass from the second.

## What is exact

- A **job** is a multiplier group H (at most `--maxo` orbits) and a row-sum pattern of one of three arrays: `plain`
  (Szöllősi's array over ℤ_N, −1 ∉ H, target 0), `wt` (the same with −1 ∈ H: all four blocks symmetric, Williamson type)
  and `bord` (the bordered array over ℤ_(N−1), target −4).
- The four lists of a job are *complete* (`enum_dfs.c`, equal to `bnw_enum.c`) and the join is *exact*: `TOTAL` is the
  number of solutions, nothing is sampled. A hit is verified in integers (`ssxx.complex_hadamard`, `bsz.build`: K K* =
  2N I) before it is written.
- A job is *decided* when it is joined or killed by an exact necessary condition: the Z_2 test of the bordered array
  (`z2fast.py`) or the quotient screen (`qkill.py`: `compress_check.py` joins the compressions to ℤ_d for every
  prime-power divisor d of the group order, then pairs of coprime levels jointly; an empty join proves that the job has
  no solution). A job that survives the screen and has lists too long to enumerate is *open*.
- Limits: at most 62 orbits per block (64-bit masks), at most 1000 classes, class values and normalised targets in
  `int16`. Lists up to `--maxlist` (2·10⁷ rows, 2C + 36 bytes each) are held in memory, the join input goes through
  `$HADAMARD_WORK/tmp` (a RAM disk is faster), and a join that runs out of memory is halved and rerun.

## Files

| File | Contents |
|---|---|
| `qhlib.py` | shared parts: `plan` (every job of N in a mode with its list sizes, pair work and Z_2 test), `enum_list` (cached, memory-mapped complete lists and prefix slices), `join` (exact lattice normalisation of every class, then `join4y`, halving of the larger list of the build pair on out-of-memory), `solve` and `run_job` (one job: lists, join, exact verification of a hit, one log line) |
| `enum_dfs.c` | pruned depth-first enumerator of the sign vectors of a block with a given row sum, same input and output as `search/bnw_enum.c`: orbits by decreasing size, a branch is cut when the row sum can no longer be reached, so the work is proportional to the number of vectors and not to 2^O; `fix0` fixes y₀ = +1; `nbits j` restricts the first `nbits` decisions to the bits of j (the 2^nbits slices partition the list) |
| `join4y.c` | `search/join4x.c` plus a probe-side filter: only the keys that some probe pair can complete get build entries (5 to 20 times faster on the large joins of the sweeps, measured there; the validations here use small lists); the same exact `TOTAL` |
| `z2fast.py` | the exact Z_2 compression test of the bordered array, as `bsz.z2_ok` with the even and odd slices summed separately |
| `qkill.py` | the exact quotient screen of a job (`compress_check.check` with `joint=True`); `(killed, per-level joins)` |
| `sched.py` | `plan` (the jobs, cheapest first), `run` (several N and modes at once in order of pair work; from `--screen` pair work on the quotient screen runs first; resumable through the logs; `--screen-only` screens the jobs with lists too long and records the survivors as open), `status` (decided, open and not run per N and mode) |
| `sliced.py` | the jobs with one list too long for memory: that list is enumerated in 2^k slices of at most `--maxslice` rows and each slice is joined with the three full lists |
| `finish_hit.py` | a hit to the certificate `C<2N>.json` in the format of `complex-hadamard/` (`szollosi-ssxx`, or `szollosi-bordered` with the border), rebuilt from the file by `verify/verify_complex.py` (K K* = 2N I, SHA-256), then `verify.complex_miyamoto(2N + 1, K)` and H H* = 4(2N + 1) I blockwise (float32 products of 0, ±1: exact) |
| `validate.py` | the checks of the section below |

## Usage

```
cd search/qhsweep
export HADAMARD_WORK=/path/to/work OMP_NUM_THREADS=1       # one core; the sweeps ran under ulimit -v 3000000 and nice
python3 sched.py plan 1271 --modes plain,wt,bord --show 10 # the jobs, cheapest first (x marks lists too long to enumerate)
python3 sched.py run 1271 --modes plain,bord --cap 1e11 --budget 3600   # exact joins by pair work, quotient screen from 1e8 on
python3 sched.py run 1271 --modes plain --screen-only      # exact quotient screen of the jobs with lists too long
python3 sched.py status 1271                               # joined, quotient kills, open, not run
python3 sliced.py 1271 --maxcost 1.2e11 --budget 3300      # one list too long: joined slice by slice
python3 finish_hit.py $HADAMARD_WORK/qhsweep/hits_plain_779.jsonl [--out DIR] [--skip-miyamoto]
python3 validate.py [enum] [join] [z2] [screen] [slice] [pipeline]       # all of them when none is named
```

Control: `python3 sched.py run 779 --modes plain --cap 1e7` finds a hit in 17 s (H = ⟨49⟩, rows (5, 7, 39, 39), pair
work 4.4·10⁶), and `python3 finish_hit.py $HADAMARD_WORK/qhsweep/hits_plain_779.jsonl` writes `C1558.json` and checks
`complex_miyamoto(1559)`, a matrix of order 6236 (29 s). The log of a mode and N is `<mode>_<N>.jsonl` (the job,
`status` join, quotient or open, `total`, seconds), the hits are `hits_<mode>_<N>.jsonl` (`N, mode, gens,
subgroup_order, rowsums, sequences`, and the `border` of the bordered array). A certificate is copied to
`complex-hadamard/` by hand.

## Validations (CHECKED)

`validate.py`, run from a copy of this directory in a stand-in repository root (Python 3.13.16, NumPy 2.5.3, gcc 13.3,
one thread, `ulimit -v 3000000`; 2 min 56 s):

- **Enumerator** (`enum`). `enum_dfs` and `search/bnw_enum.c` give the same set of (sign mask, class values), and the
  number of rows is that of `Struct.counts`, for 107 plain lists of 9 structures (N, H = 143 ⟨3⟩, 247 ⟨9, 35⟩, 119 ⟨8⟩,
  163 ⟨38⟩, 185 ⟨9⟩, 93 ⟨4⟩, 267 ⟨4⟩, 377 ⟨9⟩, 183 ⟨34⟩; both kinds, at most 22 orbits, 0 < r ≤ 2⌊√N⌋, at most 3·10⁶
  rows) and 146 bordered lists of 22 structures (N = 19, 23, 27, 31, 35, every H with at most 18 orbits, 0 ≤ r ≤ 6; 44
  of them with r = 0 and y₀ fixed): 253 lists.
- **Slices** (`slice`). The 2, 4 and 8 prefix slices partition the list, checked on the list of A of 779 ⟨49⟩ (r = 39,
  18480 rows), 371 ⟨13⟩ (r = 25, 20592), 163 ⟨38⟩ (r = 19, 43758) and the r = 0 list of L of the bordered N = 31, H =
  ⟨1⟩ (3432 rows, y₀ fixed): 12 partitions. The sliced join of 163 ⟨38⟩, 119 ⟨8⟩ and 143 ⟨3⟩ with 2, 4 and 8 slices of
  the largest list gives 8352, 4704 and 16 solutions, the stop mode a verified solution.
- **Joins** (`join`). 60 random instances on a lattice (47 with solutions, targets hit by construction or missed by
  one): `join4x`, `join4y` and `join4y` after the lattice normalisation give the same total. The five structures with
  known solution sets (NOTES §52 for 163, 119, 371):

  | N | H | rows | list sizes | solutions | |---:|---|---|---|---:| | 779 | ⟨49⟩ | 5, 7, 39, 39 | 120, 120, 18480, 18480
  | 96 | | 371 | ⟨13⟩ | 3, 3, 25, 29 | 120, 120, 20592, 20592 | 96 | | 143 | ⟨3⟩ | 3, 9, 11, 19 | 6, 12, 560, 448 | 16 |
  | 163 | ⟨38⟩ | 17, 19, 1, 1 | 126, 126, 48620, 48620 | 8352 | | 119 | ⟨8⟩ | 1, 5, 15, 15 | 120, 40, 105105, 105105 |
  4704 |

  `join4x`, `join4y` raw, `join4y` normalised, the normalised join with an out-of-memory failure forced on the first
  level (the halves add up) and `qhlib.solve` give these totals, and the first solution of each is a complex Hadamard
  matrix (K K* = 2N I, exact).
- **Z_2 test** (`z2`). `z2fast.z2_ok_fast` and `bsz.z2_ok` agree on 332 (structure, pattern) cases (N = 23, 31, 47, 59,
  71, 87, 91, 95, 111, 119, 143, 155, 183, 227, at most 22 orbits).
- **Quotient screen** (`screen`). A kill always agrees with the exact join. Of the jobs of pair work at most 10⁷ (N =
  119: 138 jobs, 143: 116, 163: none, 371: 3; 24 sampled per N), 49 are killed (119: 24 of 24, 143: 23 of 24, 371: 2 of
  3) and each has exactly 0 solutions; the two survivors are the 371 ⟨13⟩ job with rows (3, 3, 25, 29) (96 solutions)
  and a 143 ⟨3⟩ job with rows (3, 11, 9, 19) (none). The structures with known solutions are never killed: 779 ⟨49⟩
  survives with level joins 19: 48, 41: 2, joint 19x41: 96; 371 ⟨13⟩ with 7: 3, 53: 1088, 7x53: 1488; 143 ⟨3⟩ with 11:
  6, 13: 4, 11x13: 8; 163 (prime) has no level to screen, and the screen of 119 ⟨8⟩ runs out of memory under the 3 GB
  limit (not killed, reported as such).
- **End to end** (`pipeline`). `sched.py run 779 --modes plain --cap 1e7` finds a hit (17 s); `finish_hit.py` writes
  `C1558.json` (K K* = 1558 I exact, SHA-256 `697fecdf1d8acdfe89ef64747b8527eaba592068d3bbff0a91450e37f951df30`), and
  `complex_miyamoto(q = 1559)` gives a quaternary matrix of order 6236 with entries 1, −1, i, −i and H H* = 6236 I (22
  s). The hit is H = ⟨49⟩, rows (5, 7, 39, 39) like the certificate in the repository, but another of the 96 solutions
  (the SHA-256 differs). The same chain on N = 11 in each mode (`plain`, `wt` with H = {±1}, `bord` with H = ⟨9⟩ and its
  border) ends in a matrix of order 92 (q = 23), and `sliced.py 163 --maxslice 20000` finds a verified hit in 4 slices.

The full output is `validation.log` (not part of the repository; the command is `python3 validate.py` with
`$HADAMARD_WORK` set).

## Not ported

The round-52 runs also used the repository's `cchain4.py` through a driver (`chainsweep.py`, the same kind of sweep as
`search/cchain_sweep.py` that skips the jobs already decided; `chain2.py` patched `cchain4.final_join` to use the
normalised join), the hybrid estimates of the undecided structures (`gauss_mc`, HEURISTIC, nothing depends on them), the
character-bordered ranks and searches of `tbsz.py` (`run_rank.sh`, `tb_batch.py`), and private table and cache scripts
(`final_report.py`, `final_table*.py`, `janitor.py`; `sched.py status` and `qhlib.clear` replace them). The chain
rule-outs of NOTES §105 are therefore not reproduced here: the quotient screen is only the level-1 test of the chain. In
the original, `sliced.py` took the patterns the chain left open; here it takes every undecided job whose largest list
exceeds `--maxslice` and whose second-largest is at most `--second`, after the quotient screen.
