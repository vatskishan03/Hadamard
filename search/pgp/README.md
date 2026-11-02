# pgp: periodic Golay pairs invariant under a multiplier group, by exhaustive row enumeration (symmetric order 1556)

Programs behind [search/NOTES.md](../NOTES.md) §113 and the symmetric Hadamard matrix of order 1556 = 4 · 389
([symmetric/1556](../../symmetric/1556)), which comes from a periodic Golay pair of length 194 over ℤ₂ × ℤ₉₇ ≅ ℤ₁₉₄ by
the construction of [symmetric/PERIODIC_GOLAY.md](../../symmetric/PERIODIC_GOLAY.md). The engine is `src/pgp.c`. For a
prime w ≡ 1 (mod 4) and the subgroup H of odd order h of ℤ_w*, with m = (w − 1)/h ≤ 32 orbits of H on ℤ_w ∖ {0}, it
finds every periodic Golay pair (a, b) of ±1 sequences on ℤ₂ × ℤ_w that is constant on the orbits {(e, hx)} of H.

- **Rows.** Write a = (a₀, a₁) by layers, A = (a₀ + a₁)/2, C = (a₀ − a₁)/2, and B, D likewise from b. A, C are ternary
  sequences on ℤ_w with complementary supports, and so are B, D. Mapping ℤ[ℤ₂ × ℤ_w] onto the two characters of ℤ₂ shows that
  (a, b) is a periodic Golay pair ⇔ AA* + BB* = w δ and CC* + DD* = w δ in ℤ[ℤ_w]. A *row* is a pair (X, Y) of
  H-invariant ternary sequences with XX* + YY* = w δ, so a pair (a, b) is two rows (A, B), (C, D) with complementary
  supports. For (w, h) = (97, 3) a sequence is the value at 0 and a sign or 0 on each of the 32 orbits of size 3, and a
  pair has 132 free signs.
- **Row-sum congruences.** The entries of a row give |supp X| + |supp Y| = w (shift 0), and the sums x and y of X and Y
  over ℤ_w satisfy x² + y² = w. A sum is X(0) + h σ, with X(0) ∈ {0, ±1} and σ the sum of the signs on the orbits, so
  x ≡ X(0) (mod h), and the number α of orbits in the support has the parity of σ with |σ| ≤ α ≤ m. The program derives
  every such "type" (x, y, X(0), Y(0), σ_X, σ_Y, α_X, α_Y), with x, y > 0 after negating X or Y, which is a symmetry,
  and keeps one of each pair of types exchanged by the swap a ↔ b. For (97, 3), 97 = 4² + 9² leaves x = 4
  (X(0) = 1, σ_X = 1) and y = 9 (Y(0) = 0, σ_Y = 3), with α_X + α_Y = 32: fifteen types, α_X = 1, 3, …, 29, of which
  the one with α_X = 1 has no bucket. The sums of a and b are then 8 and 18 (8² + 18² = 2 · 194). For the lengths 106,
  178, 218 (two groups), 274 and 298 the congruences leave no type (see "Other results").
- **Buckets (stage 1).** With g the least primitive root (5 for w = 97), orbit j is {g^(j + mt)}; multiplication by g
  rotates the orbit index, and x ↦ −x is the rotation by m/2. The parity key K(S) ∈ {0, 1}^(m/2) of a support S has
  K(S)(k) = #{x : x ∈ S, x + g^k ∈ S} mod 2. It depends on S alone, and K(supp X) = K(supp Y) for every row, because
  PAF_X(k) = −PAF_Y(k) and PAF_X(k) ≡ #{x : x, x + g^k ∈ supp X} (mod 2). A bucket is a type together with a key. Three
  exact reductions cut the work: the rotation of a row by a multiplier (X and Y together) shifts the key cyclically, so
  only rows with a key that is minimal among its 16 cyclic shifts are enumerated (4,116 keys of 65,536); the reversal
  x ↦ −x of X alone or of Y alone gives a row again, so each component is kept only if X ≤ X* (the rotation by 16);
  and the swap a ↔ b is taken out in the types.
- **Rows in a bucket.** For each support of X with the bucket's key all sign patterns with the prescribed sum σ are
  enumerated depth-first, with the autocorrelations at the 16 shift classes g^k (k < 16) updated incrementally in 16
  int8 lanes and the Fourier values at the 16 character classes in float32 (Gaussian periods). A pattern is discarded
  when some |X̂|² exceeds w (by more than 0.05 for rounding), since |X̂|² + |Ŷ|² = w at every character. The side of the
  bucket with fewer candidates is stored in a hash table keyed by the 16-vector of autocorrelations, and the other side
  probes it with the negated vector; every hit is a row (X, Y). A "leaf" below is a complete sign pattern on a support,
  counted before the PSD filter.
- **Join (stage 2).** For every row (A, B), the rows whose supports are the complements of those of (A, B), taken under
  all 32 rotations, the reversal of the second component and the swap, are looked up; each match is a periodic Golay
  pair, which `tools/extract_pairs.py` rebuilds as two ±1 sequences on ℤ_2w, checks (exact integer PAF and numpy FFT),
  and sorts into equivalence classes (multiplier rotations, reversal of a or of b alone, negations, translation of a or
  of b by the element of order 2, a ↔ b, and the alternation y ↦ (−1)^y).
- **Checkpoint.** A bucket is recorded in `done.log` (type, key, rows, leaves, PSD passes, stored entries, thread CPU
  seconds) only after its rows file (`rows/t<type>_k<key>.bin`, 20 bytes a row) is written, and the same command
  resumes.

What an empty result proves: no periodic Golay pair of length 2w over ℤ₂ × ℤ_w has both sequences invariant under H. The
twisted family b(x) = (−1)ˣ a(gx) is another family ([search/NOTES.md](../NOTES.md) §66), and pairs with other or no
symmetry are not covered. The classes counted here are classes under the symmetry group above; equivalence in the
literature may be coarser, so the counts are upper bounds on the number of classes there.

Python 3 with NumPy (3.11 and 3.12 were used) and a C compiler with pthreads and the vector extensions of GCC and clang;
only clang was used (see "Files"). No certificate depends on these programs: `verify/` checks the matrix of order 1556
on its own. Binaries go to `src/`, the run directory (rows, `done.log`, certificate) is given on the command line, and
the validation scripts write `work/`, `out/` and `logs/` in `search/pgp`.

Labels: CHECKED = exact computation or direct measurement; HEURISTIC = estimate, nothing depends on it for correctness.

## Files

The files are byte for byte the ones that ran (md5 below), except two in which a line with a default path was edited
before publication (`tools/make_cert.py` and `tools/validate_all.sh`; see "Paths"): for these two the table gives the md5
of the file as it ran. "Run" marks a program the order 1556 run used; the others belong to the validation.

| File | md5 | Contents |
|---|---|---|
| `src/pgp.c` | `8445f67605aa5f26d7a119502b768af2` | run: the engine (735 lines): `info w h`, `selftest w h`, `plan w h`, `rows w h DIR`, `join w h DIR`, with `--swap-reduce`, `--threads T`, `--quiet`, and for tests `--all-keys`, `--no-rev`, `--no-psd`, `--types`, `--limit`, `--stride`, `--offset`. Built with `clang -O3 -march=native -pthread -o src/pgp src/pgp.c -lm` |
| `src/oracle.c` | `af1ebf1fb5493a57b29b10a6f61daabe` | independent brute force that shares no tables with `pgp.c` (only the orbit labelling): `oracle rows w h OUT` lists all rows and `oracle pgp w h OUT` all pairs on ℤ_2w, each with the autocorrelations from the definition. It stores them in int8, which is exact for length at most 127; above that a wrapped value could add a spurious pair but cannot hide one (autocorrelations of a ±1 sequence of length 146 are ≡ 2 mod 4, so ±128 does not occur), and the one comparison above 127, (73, 9), found none on both sides |
| `tools/driver.py` | `755751c8b3d832cd4d2e1ab2016eeff1` | run: `plan W H`, `run W H DIR [--threads T] [--no-cert]` (stage 1, stage 2, classes, certificate), `status W H DIR` |
| `tools/pgplib.py` | `35715c015f084785baa019e6383a0d83` | run: orbit tables, reading the row files, the symmetry group of a pair, exact PAF check |
| `tools/extract_pairs.py` | `08352982b36ddf090b900e6b495f570f` | run: one representative of each class of the join output, as ±1 strings on ℤ_2w, with the exact and the FFT check |
| `tools/make_cert.py` | `f8dbf27913884e26df2d5572837ead20` | run: pair → symmetric Hadamard matrix of order 4q through the repository's builder → certificate, with an independent check of H (H = Hᵀ, ±1, H Hᵀ = 4qI) and `verify/verify.py`; never writes into the repository |
| `tools/solve.py` | `89b8f11a2ae01c4a3b38e1d3b34e73ed` | stage 1, stage 2 and the exact checks in one command, optionally against the oracle's pair list |
| `tools/validate_all.sh` | `ac20a704a20022133d9a031bd3a0378a` | the validation suite, sections 1 to 9 below (about 100 CPU-seconds) |
| `tools/validate_extra.sh` | `7fb67b742815c66431d5479bb5ff0970` | section 10 below: class counts at the groups of Đoković and Kotsireas's pairs |
| `tools/check_known.py` | `3e6bd6b9fc445155605566026f5ab6c0` | is the pair of an existing certificate among the images of the engine's solutions |
| `tools/cmp_rows.py`, `tools/cmp_rows_canon.py` | `cdb30453fbeaccc052fb45b17921da2d`, `8f17af137f4b59983764a176bfb2e06c` | engine rows against the oracle's rows: all keys without the reversal reduction, and canonical keys closed under rotations and reversals |
| `tools/count_pairs_oracle.py` | `f02f75e4471a0618f9e47aed433d9ec8` | independent Python support join on the oracle's rows |
| `tools/timed.py` | `f3eac04113e36cfcc77dcf4c874f6785` | runs a command and appends its CPU time to `logs/cpu_ledger.txt` (used by `validate_all.sh`) |
| `data/pgp146_h3_classes.json` | `ca1842671c71e0c9181e90f9abed2c76` | the four classes of periodic Golay pairs of length 146 over ℤ₂ × ℤ₇₃ with H of order 3, found as the control (see "Other results") |

The other program of the run is `verify/verify.py` of the repository (md5 `970e5d9aba849fd38a7a4109a85a091d`, the file
of commits `919e73d` and `7b45d67`), which `make_cert.py` imports; on the machine of the run it was read from a copy of
the repository. The binaries that ran were built on that machine on 1 November 2026, after the last change of the
sources (`pgp.c` and `oracle.c` are older than the binaries), by `validate_all.sh`, with the commands in the table:
`src/pgp` `234818692872704b77146cda6d42a4f7`, `src/oracle` `68bed4767d989abc3c60b4d911fd18c2`. We did not rebuild them
there, so we cannot say that these sources give the same binaries; the build of the same sources on the second machine reproduces
the results below. The machine of the run: AMD EPYC 9R14 (32 cores, one thread each), 246 GB, Ubuntu 24.04.5 (Linux 7.0,
x86-64), clang 18.1.3, Python 3.12.3, NumPy 1.26.4. `driver.py`, `extract_pairs.py` and `make_cert.py` name
`python3.11`; on that machine it was a link to `python3` (3.12.3) in a directory put first on `PATH`. The second
machine: Apple M4 (arm64), macOS, Apple clang 21.0.0, Python 3.11.14, NumPy 1.26.4. Only clang was tried; gcc was not.

## Paths

Two of the files that ran put a default path into the code, the checkout of the repository on the author's Mac, where
the run used the environment variable instead:

- `tools/make_cert.py` (line 18) read `REPO = os.environ.get("HADAMARD_REPO", "<path of a checkout>")` and puts
  `<REPO>/verify` on `sys.path`. The run set `HADAMARD_REPO`, so the default was not used.
- `tools/validate_all.sh` (line 6) read `REPO=${HADAMARD_REPO:-<path of a checkout>}`. The suites reported below set it
  as well.

The two lines were edited before publication, so that the default is the repository that contains `search/pgp`:

```
tools/make_cert.py:    REPO = os.environ.get("HADAMARD_REPO") or os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", ".."))
tools/validate_all.sh: REPO=${HADAMARD_REPO:-$(cd ../.. && pwd)}
```

(the script has changed into `search/pgp` by then). The md5 of the two files as published are
`1edc419151d6e856c37496383615e8f8` (`make_cert.py`) and `6721a694d43c2b47fb038ace34786037` (`validate_all.sh`); as they ran
they were `f8dbf27913884e26df2d5572837ead20` and `ac20a704a20022133d9a031bd3a0378a`. CHECKED: with the edited
`make_cert.py`, in a tree with the layout of the repository and without `HADAMARD_REPO`, the certificate of order 1556 is
rebuilt byte for byte. `HADAMARD_REPO` still overrides the default. `validate_all.sh` writes its log into the directory
`logs`, which must exist (`mkdir -p logs`).

## Usage

```
cd search/pgp
clang -O3 -march=native -pthread -o src/pgp src/pgp.c -lm
clang -O3 -o src/oracle src/oracle.c
./src/pgp selftest 97 3
mkdir -p logs && bash tools/validate_all.sh && bash tools/validate_extra.sh
```

**The run for q = 389** (the scheduler started it from `search/pgp` with a directory holding the link `python3.11` first
on `PATH`, `HADAMARD_REPO` set, and `THREADS` = 8; `DIR` is the run directory):

```
python3 tools/driver.py plan 97 3                        # 7 s; prints the bucket and leaf counts of the table below
python3 tools/driver.py run 97 3 DIR --threads 8         # resumable: the same command resumes after an interruption
python3 tools/driver.py status 97 3 DIR                  # progress from DIR/done.log
```

`run` does the following, in this order: `src/pgp rows 97 3 DIR --swap-reduce --threads 8`; `src/pgp join 97 3 DIR
--swap-reduce > DIR/solutions.txt`; `python3.11 tools/extract_pairs.py 97 3 DIR --swap-reduce > DIR/classes.json`;
`python3.11 tools/make_cert.py --q 389 --pair DIR/classes.json --class 0 --out DIR/cert`, which writes
`DIR/cert/symmetric/1556/certificate.json` and runs `verify/verify.py` on it.

**Time and resources** (CHECKED, from the scheduler's records and the run logs; the machine was shared with other jobs
of the census, 8 of its 32 cores were this job's). The job was started and resumed four times: the scheduler that runs
the jobs was stopped three times, which terminated the engine each time, and `done.log` kept the finished buckets. The
logs of the first two stops say only "scheduler stopping"; the scheduler restarted after the second had a longer job
list (46 jobs instead of 13), and no reason is recorded for either. The third stop was a crash of the scheduler's own
script (`KeyError`, while it processed the hit of another job).

| Run | Started (UTC) | Stopped | Wall (s) | Buckets done before the run | Buckets finished in the run |
|---|---|---|---:|---:|---:|
| 1 | 2026-11-01 14:58:32 | 15:03:47 | 315 | 0 | 46 |
| 2 | 2026-11-01 15:04:58 | 16:41:01 | 5,763 | 46 | 1,615 |
| 3 | 2026-11-01 16:41:25 | 17:22:38 | 2,473 | 1,661 | 766 |
| 4 | 2026-11-01 17:23:18 | 2026-11-02 03:54:37 (finished) | 37,879 | 2,427 | 21,820 |

The scheduler's total for the job is 46,419 s (12.89 h) of wall time on 8 threads (the four times above add up to 46,430
s, 11 s more, since a stop takes a moment), 371,351 thread-seconds (103.2 h) allocated. The sum of the thread CPU times
of the 24,247 finished buckets (last column of `done.log`) is 357,760 s (99.38 h), 2.49 · 10⁷ leaves per second per
thread on average; the difference of 3.8 h is work of buckets interrupted at the stops and redone, the construction of
the support indices (about 18 s each run), and time not broken down further. The most expensive bucket took 29.6 s of
thread CPU (type 10, key 987; the buckets are started heaviest first, by their estimated leaves). Run 4 took 37,820.5 s
for stage 1 (21,820 buckets, 7.223 · 10¹² leaves, 890,184 rows) and 49 s for stage 2 and the extraction of the class;
the certificate followed in about 6 s. Peak resident memory of the job was 1.12 GB (scheduler's record); the rows files
are 12,506 files of 18.0 MB in all.

**Result** (CHECKED; the totals are the sums over `done.log`, recomputed here). All 24,247 buckets of the plan were
finished, with the plan's leaf count per type to the five digits it prints: 8,894,476,872,902 leaves, of which
15,824,365,707 (0.178 %) passed the PSD filter, 11,306,171,520 stored entries, 897,576 rows in 12,506 buckets (no
duplicates).

| Type | α_X | α_Y | Buckets | Leaves | Rows | Thread CPU (h) |
|---:|---:|---:|---:|---:|---:|---:|
| 1 | 3 | 29 | 16 | 1.0858 · 10⁹ | 0 | 0.01 |
| 2 | 5 | 27 | 1,884 | 1.0524 · 10¹¹ | 0 | 1.09 |
| 3 | 7 | 25 | 2,022 | 4.6591 · 10¹¹ | 1 | 4.95 |
| 4 | 9 | 23 | 2,055 | 1.0078 · 10¹² | 242 | 10.91 |
| 5 | 11 | 21 | 2,048 | 1.1870 · 10¹² | 9,983 | 13.27 |
| 6 | 13 | 19 | 2,066 | 8.4526 · 10¹¹ | 110,368 | 9.96 |
| 7 | 15 | 17 | 2,048 | 4.5758 · 10¹¹ | 348,693 | 5.90 |
| 8 | 17 | 15 | 2,066 | 5.2205 · 10¹¹ | 329,812 | 6.56 |
| 9 | 19 | 13 | 2,048 | 1.0168 · 10¹² | 91,116 | 11.59 |
| 10 | 21 | 11 | 2,066 | 1.4339 · 10¹² | 7,180 | 15.63 |
| 11 | 23 | 9 | 2,040 | 1.1829 · 10¹² | 181 | 12.53 |
| 12 | 25 | 7 | 2,037 | 5.4762 · 10¹¹ | 0 | 5.72 |
| 13 | 27 | 5 | 1,838 | 1.2021 · 10¹¹ | 0 | 1.23 |
| 14 | 29 | 3 | 13 | 1.0858 · 10⁹ | 0 | 0.01 |
| Total | | | 24,247 | 8.8945 · 10¹² | 897,576 | 99.38 |

Type 0 (α_X = 1) has no bucket. Stage 2 printed 2 matches (`solutions.txt`, two lines): the first rows of both lie in
type 9, in the buckets with keys 6587 and 5017, and the partner rows are of type 9 as well. They are one equivalence
class of pairs, 4,096 ordered pairs (a, b) under the symmetry group of `pgplib.gen_group_images` (CHECKED, recomputed
here); its representative is the pair of `symmetric/1556/certificate.json`, with sums 8 and 18, layer sums (13, −5, 13,
5), and exact integer PAF_a + PAF_b = 388 δ on ℤ₁₉₄. Both sequences are invariant under multiplication of ℤ₁₉₄ by 35
(and so by 61 = 35²). So for the family with |H| = 3 over ℤ₂ × ℤ₉₇, given the validation below, the pair is unique up to
the symmetries.

The cost was estimated before the run from the rate on one Apple M4 core (3.4 to 4.6 · 10⁷ leaves per second) at 54 to
73 core-hours; the run took 99.4 thread-CPU hours at 2.49 · 10⁷ leaves per second (HEURISTIC: the ratio of the two rates
depends on the machine, and the estimate was made without a measurement on x86).

## Validation

Sources: **L** is the log of `validate_all.sh` on the machine of the run (1 November 2026, 14:25 UTC); **R** is a rerun
of `validate_all.sh` and `validate_extra.sh` on the second machine on 2 November 2026, whose output is identical to L in
every count, class-set hash and certificate SHA-256 (the differences are times, paths, the git line and the ledger
total, see below); **C** was computed here, for this README, on the second machine from files copied unchanged from the
run directory.

| Check | Result | Source |
|---|---|---|
| 1. Self-test: PAF tables, parity keys, Fourier tables, rotation covariance of keys and PAF against direct computation (300 random ternary sequences), and for m ≤ 16 the support index against direct keys; (w, h) = (5,1), (13,3), (13,1), (37,3), (41,5), (61,5), (73,9), (113,7), (101,5), (97,3) | all ok, 0 problems (for (97, 3), m = 32, the index check is not run) | L, R |
| 2. Rows against the oracle, all H-invariant ternary pairs with the autocorrelations from the definition: (5,1) 200 rows, (41,5) 352, (13,1) 493,480, (37,3) 11,592; once with all keys and no reversal reduction, once with canonical keys and the reversal reduction, closed under rotations and reversals (30, 22, 23,504 and 560 canonical rows) | sets equal in all 8 comparisons | L, R |
| 3. Pairs (stage 1 + 2) against the oracle over all H-invariant ±1 pairs on ℤ_2w, with and without the swap reduction: (5,1) 6,400 ordered pairs (15 classes), (13,3) 256 (1), (41,5) 4,608 (5); (29,7), (73,9), (37,9) none | sets equal; 0 = 0 | L, R |
| 4. Stage 2 against an independent Python support join of the oracle's rows, (13,1), m = 12 (493,480 rows) | 822,016 ordered pairs of rows with complementary supports = 822,016 distinct quadruples (A, B, C, D) of the engine (all keys, no reversal reduction) | L, R |
| 5. (37,3), H = ⟨10⟩: the note of §65 gives 35,328 pairs | 35,328 ordered pairs, 66 representatives, 27 classes (class-set hash `8810877f39c9`); Đoković and Kotsireas's pair of length 74 (certificate of `symmetric/596`) is among the images (2 of its 1,536) | L, R |
| 6. Threads, checkpoint: (61,3) in one run, and in two runs with 4 threads (100 buckets, then the resume) | 370 buckets in both; the same 13 classes (hash `483a51ec4ac8`) | L, R |
| 7. Control (73,3), 4 threads: 1,456 buckets, 2,135,181,332 leaves (equal to the plan), 20,000,776 PSD passes, 12,548,350 stored, 7,452,426 probes, 43,238 rows; 8 matches | 4 classes of periodic Golay pairs of length 146; exact integer PAF and FFT both ok; stage 1 took 21.9 s wall and 87.4 CPU-s on the machine of the run, 22.8 s and 87.1 CPU-s on the second | L, R |
| 8. Rows of a sample run at m = 32 are exact rows | vacuous in L and R: 0 rows were checked (the sample run's directory holds no rows files). Replaced by C below | L, R |
| 9. Builder: q = 149 and q = 101 from the repository's own pairs; q = 149 from a pair of this engine (class 5 of 27 at (37,3)); q = 293 (order 1172) from a pair of length 146 | the first two reproduce the SHA-256 of `symmetric/596` and `symmetric/404` and the pairs; all four pass the independent check and `verify/verify.py`. The order 1172 is a test artifact: v = 293 is not undecided in the census | L, R |
| 10. Exhaustive class counts at the groups of Đoković and Kotsireas's pairs (their numbers of solutions, as recorded in `validate_extra.sh`, in brackets; not re-read in the paper): (37,3) 27 [74: 2], (41,5) 5 [82: 2], (61,5) 1 [122: 1], (101,5) 1 [202: 1], (113,7) 2 [226: 2], (61,3) 13, (73,3) 4, (109,9) 0, (181,9) 0 | the counts are the classes of this engine's group | R |
| C1. The plan's table recomputed for (97,3) | identical to `logs/plan_97_3.txt` of the run (same md5) | C |
| C2. Every row file of the run | 12,506 files, 897,576 rows, none repeated; the row count of every file equals the `done.log` entry; every row satisfies PAF_X + PAF_Y = 97 δ exactly on ℤ₉₇ (own script, integer arithmetic) | C |
| C3. Stage 2 replayed on the rows of the run (second machine, 15 s) | `solutions.txt` and `classes.json` byte for byte the files of the run (md5 `3dcbbda6a659da06379adac8ca84cfd0`, `b8a6df4fd7d548d3f96744ff3608b4fe`) | C |
| C4. The two buckets that hold the first row of each match (type 9, keys 5017 and 6587) enumerated again on the second machine | the rows files are byte for byte those of the run, with the same leaf, PSD-pass and stored counts (497,782,428 and 507,711,347 leaves; 53 and 44 rows) | C |
| C5. The certificate rebuilt by `make_cert.py --q 389` from `classes.json` and checked with `verify/verify.py` | byte for byte `certificate.json`; SHA-256 `3eea11362d17330271f31be78129dc596364c4ba99a6540f01abe7584e3c85ec` | C |

Not covered by these checks: that the PSD filter loses no row at m = 32. The filter can only discard a sequence that is
part of no row, apart from rounding; the margin of 0.05 is several hundred times a first-order estimate of the float32
error of a Fourier value (about 10⁻⁴, HEURISTIC), and no row was lost in the comparisons with the oracle (m ≤ 12), but
no oracle exists at m = 32. The log L has three more blemishes that do not touch the results: its last step, the CPU
total of the ledger, ended with a Python error because earlier lines of the ledger have labels with a space; the
repository copy on that machine is not a git checkout, so the line "working tree changes made by this suite: 0" is not a
check there (it is on the second machine, where it is 0 for the repository at commit `7b45d67`); and section 10 is not
in L, because the rerun of `validate_all.sh` on that machine overwrote the log file the earlier run had appended it to.

## Other results of the engine

- **Length 146** (w = 73, |H| = 3, m = 24), the control: 2,135,181,332 leaves, 43,238 rows, 8 matches, 4 equivalence
  classes of periodic Golay pairs of length 146 over ℤ₂ × ℤ₇₃ (CHECKED: exact integer PAF and FFT; row sums ±16 and ±6),
  in `data/pgp146_h3_classes.json`. Length 146 is one of the open lengths below 300 in Đoković and Kotsireas's list
  (arXiv:1310.5773), and we found no source that gives a pair of this length ([NOVELTY.md](../../NOVELTY.md)). It settles
  no census value: v = 293 is decided in the census.
- **No admissible type**: for (w, |H|) = (53, 13), (89, 11), (109, 9), (109, 27), (137, 17) and (149, 37), lengths 106,
  178, 218, 218, 274 and 298, no representation w = x² + y² has both x and y congruent to 0, 1 or −1 modulo |H|, so the
  plan is empty. The scheduler's jobs for these six finished within 6 s each, with no rows and no pair (CHECKED: `pgp
  info` prints no type for any of them on the second machine, and the congruences were checked by hand). Likewise (157,
  13), (229, 19), (277, 23), (29, 7), (73, 9) and (37, 9) have no type (C), and (181, 9) has 375 buckets and no pair.
- **Class counts** (the classes of this engine; R and C): (13,3) 1; (37,3) 27; (41,5) 5; (61,5) 1; (61,3) 13; (113,7) 2;
  (101,5) 1; (73,3) 4; (97,3) 1. Along |H| = 3 the counts at m = 12, 20, 24, 32 are 27, 13, 4, 1.

## Certification

`tools/make_cert.py` builds and checks the certificate and never writes into the repository: it checks that (a, b) is a
periodic Golay pair over the stated group and that t has order 2, builds H through the repository's builder
(`verify/verify.py`, construction `symmetric-involution`), re-checks H = Hᵀ, the ±1 entries and H Hᵀ = 4qI with a
float64 product, and runs `verify/verify.py` on the new directory. For this order

```
python3 tools/make_cert.py --q 389 --pair DIR/classes.json --class 0 --out OUT      # OUT/symmetric/1556/certificate.json
python3 verify/verify.py OUT/symmetric/1556
```

The same pair, as two ±1 strings of length 194 over ℤ₁₉₄, is stored in `search/sym_involution.py` (`PAIRS[389]`), and
`python3 search/sym_involution.py certify 389` writes `symmetric/1556/certificate.json` from it, byte for byte the file
in this repository.
