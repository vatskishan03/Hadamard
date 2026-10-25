# flatpool: the flat-pool k-tree for Goethals-Seidel quadruples

Programs and results of [../NOTES.md](../NOTES.md) §101: four-list joins restricted to near-flat orbit-sign vectors (the
flatness that Egan, Kotsireas, London and Ó Catháin observed in the matrices of Alpöge et al.), calibrated at v = 487, 311 and
1103, with the predicted cost for order 4412 and the other prime FRONTIER structures.

Labels: CHECKED = exact computation or direct measurement in this directory; HEURISTIC = model extrapolation; PROVED = none claimed.
All CPU times are process CPU seconds (the machine was loaded, load average 7-10 on 4 cores; wall time was 2-4x CPU time).

## 1. Verdict

* The flat-restricted idea works as a pipeline: a genuine Goethals-Seidel quadruple of order 1948 = 4 * 487 (U = <41>, row sums
  (1, 17, 17, 37), Phi = (0.129, 0.138, 0.114, 0.121)) was found from pools of near-flat sequences in about 3500 CPU-seconds
  (1219 s generation + 2283 s join; `sol487.json`), and verified exactly (CHECKED, section 3).
* For v = 1103 the model fed with measured rates predicts >= 1.5e10 CPU-seconds (optimistic bound) to 1e12 CPU-seconds
  (implemented chunking) with the 3.5 GB address-space cap, and still 1.6e7 CPU-seconds (185 CPU-days) with 2 TB of memory
  (HEURISTIC, section 5). No production job was run.
* The earlier judge's figure (10^12.9 list operations at lists of 10^9) is too optimistic by roughly 10^2 to 10^3, for three
  measured reasons: the flat fraction p(Phi) is 10-25 times smaller than its Chernoff model; the Gaussian estimate of the
  number of solutions is 10-15 times too large for the SA-selected lists (at 28 classes); and a four-list k-tree on s merged
  classes gains (sqrt 2)^s, not (6/pi)^((C-1)/2). At equal list size the measured flat-pool join (10^10.6 pair operations at
  L = 2.7e7 for v = 487) is indistinguishable from the formula of `search/triage.py` for plain lists (10^10.5).
* Generation of near-flat vectors is NOT the bottleneck (section 2): about 10^4 to 7*10^4 distinct flat necklaces per CPU-second
  at Phi <= 0.15 to 0.18 for v = 1103. Memory and the join are.

## 2. Generation of near-flat orbit-sign vectors (v = 1103, one row-sum class, one core)

Orbit sign vectors are 58 signs on the cosets rho^j U (cyclic order, rho = 5); multiplication by rho rotates them and leaves
Phi unchanged, so a "necklace" (58 rotations, 29 distinct PAF-key vectors) is the unit of generation. `flat.c`: Metropolis
on swaps (+ <-> -) with restarts every K = 10^4 proposals, tau = 0.3-0.4, every proposal with Phi <= thr is recorded and
deduplicated up to rotation (CHECKED, 2e7 proposals/s).

| Phi threshold | REM (random swap-neighbours) p(Phi <= thr) | REM rate (distinct/s) | SA distinct necklaces/s | gain |
|---|---|---|---|---|
| 0.10 | < 3e-10 (0 hits in 3.4e9) | - | 110 | - |
| 0.11 | < 9e-10 | ~0.02 | 330 | - |
| 0.12 | 3.3e-9 | 0.12 | 1.0e3 | 8e3 |
| 0.14 | 1.8e-8 | 0.67 | 5.7e3 | 8e3 |
| 0.15 | 5.5e-8 | 2 | 1.2e4 | 6e3 |
| 0.18 | 6.9e-7 | 25 | 6.8e4 | 2.7e3 |
| 0.24 | 2.7e-5 | 1e3 | 6.2e5 | 6e2 |

* p(Phi) vs the Chernoff model of `judge_flat2/ld.py`: 3.3e-9 vs 4.1e-8 (0.12), 1.8e-8 vs 4.4e-7 (0.14), 1.4e-7 vs 3.4e-6 (0.16),
  2.7e-6 vs 8.7e-5 (0.20): the model overstates the flat pool by 12-30x (CHECKED, 3.4e9 evaluations).
* Pool sizes at v = 1103, row sum 1 (necklaces = C(58,29) p / 58): 1.7e6 (0.12), 9.3e6 (0.14), 7e7 (0.16), 3.6e8 (0.18), 1.4e9 (0.20).
  Capture-recapture of the SA at v = 487 (Phi <= 0.14) gives 1.7-1.8e6 necklaces per role against 2.2-3.0e6 from REM; at
  v = 1103, Phi <= 0.11 it gives ~3e5 against <= 1.7e5 from REM: the SA covers the pool about uniformly (factor 2).
* The coordinator's data point (v = 1103, row sum 37; greedy swap-polish from random starts ends at Phi quantiles (min, 10%, 50%, 90%) =
  (0.11, 0.12, 0.15, 0.18) at ~0.1 s/sample in numpy; the least-squares phase decoding gives (0.09, 0.12, 0.14, 0.17)) agrees with the
  REM table: p(Phi <= 0.15) ~ 5e-8 per random vector and a polished median of 0.15. The C annealer reaches the same range at
  10^3-10^5 distinct necklaces/s (10^2-10^4 times the numpy rate), because each recorded sample is a swap-neighbour of an annealed state.
* `search/lowe_enum.c` cannot be the generator at 58 cosets: it visits every leaf of the sign tree (energy is tested only at the leaf,
  lines 126-136; reachability pruning only on the row sum), i.e. C(58,29) = 3.0e16 leaves for row sum 1 (CHECKED by reading the code).
  Exhaustive pools below Phi = 0.10 are tiny (<= 1.7e5 necklaces by the REM bound) and hold no solution: with 10^5 necklaces
  per role the model gives E ~ 10^-4 solution classes.

## 3. Genuine solution at v = 487 (CHECKED)

Structure: Z_487, U = <41> (order 9, 54 cosets + {0}, 27 classes), row sums (1, 17, 17, 37) (triage pattern; the other pattern
(1, 19, 19, 35) was not used). Pools: SA at Phi <= 0.14, tau = 0.4, K = 10^4: 1.0e6 necklaces for each of row sums 1, 17, 37
(generation 391 + 435 + 392 = 1219 CPU-s). Lists: all four roles rotated (27 rotations), 2.7e7 entries each; B and C from the
same pool (n = 1e6 each). Join `rj` with 9 merge classes, buckets in order of decreasing density; the first and only solution
appeared at bucket 1238 (1239 buckets, 2283 CPU-s, 1.84 s/bucket, 2.0e7 pairs/bucket).

SOL line: `675370 24 190686 1 724700 12 91673 13` (necklace index, rotation), decoded by `decode.py` into `sol487.json`:
four sequences of length 487 as orbit signs (coset index j of rho^j U, rho = 3) and as +- strings. Checks done with exact integer
arithmetic: row sums (1, 17, 17, 37), sum of squares 1948 = 4 * 487; Sum_i PAF_i(d) = 0 for all d != 0 by direct O(v^2)
summation (and by FFT in `gsverify.py`); each sequence is constant on the orbits of <41>; Goethals-Seidel array of the four
circulants (order 1948) with H H^T = 1948 I in int64 and entries +-1. Flatness of the four blocks: Phi = 0.129, 0.138, 0.114,
0.121 (published 487: 0.14, 0.14, 0.12, 0.12). 1948 was already known (Alpoge et al.); this is an independent re-derivation.
v = 311: no solution was attempted at the join level (the 311 lists that fit in 3.5 GB hold E ~ 10 tuples); 311 was used for the
model calibration only.

## 4. Calibration of the solution-count model

Gaussian-lattice estimate (full covariance of the rotation-symmetrised lists, `est.py`) vs counts of partial solutions
(`rj2` with RJ_NC = number of constrained classes; counts = hash matches verified on NC classes over the first buckets):

| v | NC (classes) | found / predicted |
|---|---|---|
| 487 | 17, 20, 22, 23 | 0.185, 0.123, 0.093, 0.033 (3 matches; 0 of 10 predicted at 24) |
| 311 | 18, 21, 24 | 0.214, 0.152, 0.159 (36 matches) |
| 1103 | 14, 16, 18, 20 | 0.320, 0.264, 0.223, 0 of 14 predicted |
| 487 full 26 classes | 26 | 1 found / 11.7 predicted over 1400 buckets = 0.085 |

The ratio decays by ~0.9 per constrained class (the SA lists are ball-like, platykurtic), so the Gaussian formula overstates the
number of solutions by ~10-15 at 26-28 classes. The model uses 0.223 * 0.9^(NC - 18). Pair counts per bucket: Gaussian/actual = 0.66
(487), 0.69 (311), 0.73 (1103). Independent-class vs full-covariance estimate: full is 10^0.5 larger (487, 311, 1103).

Join constants (`rj`, one core, 2.7e7-entry lists): 1.8e7 pairs/s (hash join of the two pair lists of a bucket), scan of the four
sorted lists 5.8e-9 s/entry/bucket (0.63 s), index build 1.6e-7 s/entry, 20 bytes/entry + 16 bytes/pair.

Model vs run at 487 (`model.py` with these constants, memory 3.3 GB): E = 1.7e3 solution tuples (63 classes), 2257 buckets and
5.6e3 CPU-s expected for the first solution, 1381 CPU-s of generation; observed 1239 buckets, 2283 CPU-s, 1219 CPU-s (Poisson-compatible:
P(N <= 1238) ~ 40%). The uncalibrated Gaussian model would have predicted ~120 buckets.

## 5. Prediction for v = 1103, U = <17> (order 19), row sums (1, 37, 39, 39)

Row-sum patterns (CHECKED, `gsflat.patterns`): r = f0 + 19 S with S even gives |r| in {1, 37, 39, 75, ...}; 4412 = r1^2+r2^2+r3^2+r4^2 has the single
solution (1, 37, 39, 39) (f0 = +1, S = 0; f0 = -1, S = 2; f0 = +1, S = 2). U of order 29 (38 cosets) or 551: no pattern.
Pools C(58,29) = 3.0e16, C(58,30) = 2.9e16, 2.9e16, 2.9e16.

Model optimum over (Phi threshold, pool fraction f <= 0.9) subject to E_classes >= 1.5 (HEURISTIC; p(Phi), SA rates and join
constants measured; E calibrated as above; "partition" = optimistic bound with sub-lists selected by exact class values, not
implemented; "chunk" = what `rj` would do with arbitrary chunks):

| memory per process | best Phi | necklaces/role | E_classes | generation CPU-s | join CPU-s | total CPU-s |
|---|---|---|---|---|---|---|
| 3.3 GB, partition bound | 0.17 | 1.2e8 | 80 | 2e4 | 1.5e10 | 1.5e10 |
| 3.3 GB, chunk | 0.16 | 5e7 | 7 | 1.8e4 | 1.0e12 | 1.0e12 |
| 64 GB | 0.16 | 4e7 | 3 | 1e4 | 1.3e9 | 1.3e9 |
| 256 GB | 0.16 | 7e7 | 20 | 3e4 | 8e7 | 8e7 |
| 2 TB | 0.18 | 3e8 | 2e3 | 5e4 | 1.6e7 | 1.6e7 |

For comparison `search/triage.py` gives 10^15.5 list operations = 10^8.2 CPU-s at 1.8e7 operations/s. Two cores for 48 hours (3.5e5 CPU-s) fall
short by a factor 4e4 (optimistic bound) to 3e6 (chunking), and by 50 even with 2 TB. Every row of the table needs lists of
>= 10^9 entries per role (E_classes >= 1 needs ~2e7 flat necklaces per role, 6e8 entries with rotations, 12 GB of keys per list).

## 6. Other FRONTIER Goethals-Seidel structures (prime v only; `table.md`, `table.json`)

Same model; best pattern = the one with the largest minimal pool. log10 CPU-seconds, memory 3.3 GB (partition / chunk), 64 GB, 256 GB, 2 TB:

| order | v | \|U\| | classes | pattern | triage log10 E | triage log10 W | 3.3 GB | 64 GB | 256 GB | 2 TB |
|---|---|---|---|---|---|---|---|---|---|
| 2956 | 739 | 9 | 41 | (1,19,35,37) | 33.5 | 22.7 | 16.9 / 20.9 | 18.6 | 18.0 | 17.0 |
| 4412 | 1103 | 19 | 29 | (1,37,39,39) | 21.1 | 15.5 | 10.2 / 12.0 | 9.1 | 7.9 | 7.2 |
| 4892 | 1223 | 13 | 47 | (25,27,27,53) | 33.7 | 32.4 | > 40 | > 40 | 30.6 | 29.5 |
| 5164 | 1291 | 15 | 43 | (1,29,29,59) | 30.2 | 29.3 | > 40 | 27.7 | 26.5 | 25.2 |
| 6668 | 1667 | 17 | 49 | (1,33,33,67) | 32.0 | 37.5 | > 40 | > 40 | > 40 | > 40 |
| 7148 | 1787 | 19 | 47 | (1,1,39,75) | 30.0 | 36.1 | > 40 | > 40 | > 40 | 33.0 |
| 11884 | 2971 | 55 | 27 | (1,1,1,109) | 14.0 | 19.2 | 14.4 / 17.4 | 14.4 | 13.0 | 11.2 |

Triage log10 W is in list operations; subtract 7.3 (1.8e7 operations/s) to compare with the CPU-s columns: 1103 -> 8.2, 2971 -> 11.9, 739 -> 15.4,
1291 -> 22.0, 1223 -> 25.1. At 2 TB the flat pool is cheaper than the plain estimate by ~1 order for 1103 (7.2 vs 8.2) and ~0.7 for 2971
(11.2 vs 11.9), and more expensive for 739 (17.0 vs 15.4); for 1223 and 1291 the flat model is 3-4.5 orders worse than the plain formula (HEURISTIC; the
triage formula assumes a selection gain (L/N)^3 that was not measured here). In no case is the gain close to the 3-4 orders that
the 10^12.9 estimate suggested. Composite v (955 = 5*191, 1257 = 3*419, 1713 = 3*571, 1893 = 3*631, 2215 = 5*443) were not modelled: their plain
log10 E is 2.9, 12.7, 5.4, 4.4, 1.9, and the flat sub-ensemble keeps a fraction of order 10^-15 to 10^-25 of the solutions, so none is expected
(HEURISTIC). For 739, 1223-1787 the REM table has fewer than 10 hits below Phi = 0.19-0.22, so p(Phi) there is a power-law extrapolation.

## 7. Programs

| file | role |
|---|---|
| `gsflat.py` | orbit/coset tables, Gaussian periods, exact PAF, row-sum options, patterns `(r1..r4)` for (v, \|U\|) |
| `mktab.py` | writes `tab_<v>.txt` (v, m, C, \|U\|, v mod 4, eta_re[m], eta_im[m]) for the C programs |
| `flat.c` -> `flat2` | `gen2` simulated-annealing generator of flat necklaces, `rem` (REM sampling of p(Phi)), `keys`, `phi` (uint64 masks, m <= 64) |
| `flatbig.c` -> `flatbig32/64` | `rem` and `gen2` with 128-bit masks (m <= 128, C <= 32 or 64) for the other FRONTIER structures |
| `genpool.py` | restart-safe pool generator: chunks `<pooldir>/c<i>.bin` + `c<i>.done`; stops at TARGET distinct necklaces |
| `mergepool.py` | merges finished chunks into `<pooldir>/pool.bin` (sorted unique canonical masks, uint64) |
| `buildrj.py` | builds the join input (int8 key vectors of the chosen necklaces) + `.map.npz` |
| `rj.c` -> `rj`, `rj2` | rotation-virtual 4-list k-tree (`rj2` = `rj` with env `RJ_NC` = number of classes constrained, for model calibration) |
| `runjoin.sh` | restart-safe driver of `rj` over bucket ranges (`.done` markers), `PAR/NPAR` to split between two processes |
| `decode.py`, `gsverify.py` | map a `SOL` line back to orbit signs, full +-1 sequences, exact PAF check, Goethals-Seidel array, `H H^T = 4v I` |
| `est.py`, `estpools.py`, `buckmodel.py` | Gaussian estimate of the number of solutions in given lists / per bucket (calibration) |
| `model.py`, `predict.py`, `tablegen.py` | cost model with measured constants; scan over (Phi threshold, pool fraction); per-structure measurements |

## Pipeline and how to resume

1. Tables: `python3 mktab.py v |U| [generator]`  -> `tab_<v>.txt` (cheap).
2. Pools (one per distinct row sum r; f0, npos from `gsflat.rowsum_options`):
   `python3 genpool.py POOLDIR tab f0 npos THR TAU K TARGET [CHUNK=2e8] [SEED0=1] [MAXCHUNKS]`
   e.g. v=487: `genpool.py p487/r1 tab_487.txt 1 27 0.14 0.4 10000 1e6 3e8 1`.  Re-running the same command resumes
   (finished chunks are those with a `.done` file; a chunk without it is regenerated).  `python3 mergepool.py POOLDIR` makes `POOLDIR/pool.bin`.
3. Join input: `python3 buildrj.py OUT.bin tab C -3 n:poolA:f0:1 n:poolB:f0:1 n:poolC:f0:1 n:poolD:f0:1 SEED`
   (all four roles rotated; T = -(v mod 4) = -3 for v = 3 mod 4).
4. Join: `./runjoin.sh $PWD/OUT.bin 9 0 NBUCKETS 200 $PWD/jdir 6e7 PAR NPAR`  (run it twice, PAR=0 and PAR=1 with NPAR=2, to use two
   cores).  Re-running the same command skips finished ranges.  A solution line is appended to `jdir/sol_all.txt`; `touch jdir/STOP` ends the drivers.
5. Verify: `python3 decode.py v |U| OUT.bin.map.npz F0A,F0B,F0C,F0D "SOL x0 k0 x1 k1 x2 k2 x3 k3" out.json`
   (exact: PAF sum zero at every shift, Goethals-Seidel array, `H H^T = 4v I`).

Memory: `rj` needs about 20 bytes per list entry (entries = C * necklaces per role) + 16 bytes per pair of one bucket; with the
3.5 GB address-space cap that means about 2.7e7 entries per list (n = 1e6 necklaces at C = 27).


Data here: `sol487.json` (the verified quadruple of order 1948), `tab_<v>.txt` (coset tables), `tg_<v>_<u>.json`
(measured p(Phi) and generation rates), `table.md` and `table.json` (the predictions of section 6). The C programs build
with `gcc -O3 -march=native -o <name> <name>.c -lm` (`flatbig.c` includes `flatbig_base.c`; see its header for the mask width).
