# flatpool: the flat-pool k-tree for Goethals–Seidel quadruples

Programs and results of [../NOTES.md](../NOTES.md) §101: four-list joins restricted to near-flat orbit-sign vectors (the
flatness that Egan, Kotsireas, London and Ó Catháin observed in the matrices of Alpöge et al.), calibrated at v = 487,
311 and 1103, with the predicted cost for order 4412 and the other prime structures of [FRONTIER.md](../../FRONTIER.md).
The flatness is Φ = (1/C) Σ_c (PSD_c/v − 1)², as defined in the header of `flat.c`.

Labels: CHECKED = exact computation or direct measurement in this directory; HEURISTIC = model extrapolation. No
statement here is a proof. All CPU times are process CPU seconds (the machine was loaded, with a load average of 7–10 on
4 cores; wall time was 2 to 4 times the CPU time).

## 1. Summary

- The flat-restricted idea works as a pipeline. A Goethals–Seidel quadruple of order 1948 = 4 · 487 (U = ⟨41⟩, row sums
  (1, 17, 17, 37), Φ = (0.129, 0.138, 0.114, 0.121)) was found from pools of near-flat sequences in about 3500
  CPU-seconds (1219 s generation + 2283 s join; `sol487.json`) and verified exactly (CHECKED, section 3). The order was
  already known.
- For v = 1103 the model, fed with measured rates, predicts at least 1.5 · 10¹⁰ CPU-seconds (optimistic bound) to
  10¹² CPU-seconds (implemented chunking) under a 3.5 GB address-space cap, and still 1.6 · 10⁷ CPU-seconds (185
  CPU-days) with 2 TB of memory (HEURISTIC, section 5). No production job was run.
- An earlier estimate of 10¹²·⁹ list operations at lists of 10⁹ is too optimistic by roughly 10² to 10³, for three
  measured reasons. The flat fraction p(Φ) is 10 to 25 times smaller than its Chernoff model. The Gaussian estimate of
  the number of solutions is 10 to 15 times too large for the annealed (SA) lists, at 28 classes. A four-list k-tree on
  s merged classes gains (√2)^s, not (6/π)^((C − 1)/2). At equal list size the measured flat-pool join (10¹⁰·⁶ pair
  operations at L = 2.7 · 10⁷ for v = 487) is indistinguishable from the formula of `search/triage.py` for plain
  lists (10¹⁰·⁵).
- Generation of near-flat vectors is not the bottleneck (section 2): about 10⁴ to 7 · 10⁴ distinct flat necklaces per
  CPU-second at Φ ≤ 0.15 to 0.18 for v = 1103. Memory and the join are.

## 2. Generation of near-flat orbit-sign vectors

Setting: v = 1103, one row-sum class, one core. An orbit sign vector has 58 signs on the cosets ρʲU (cyclic order,
ρ = 5). Multiplication by ρ rotates it and leaves Φ unchanged, so a necklace (58 rotations, 29 distinct PAF-key vectors)
is the unit of generation. `flat.c` runs a Metropolis walk on swaps (+ ↔ −) with restarts every K = 10⁴ proposals and
τ = 0.3 to 0.4. Every proposal with Φ ≤ thr is recorded and deduplicated up to rotation (CHECKED, 2 · 10⁷ proposals/s).

| Φ threshold | REM (random swap-neighbours), p(Φ ≤ thr) | REM rate (distinct/s) | SA distinct necklaces/s | gain |
|---|---|---|---|---|
| 0.10 | < 3 · 10⁻¹⁰ (0 hits in 3.4 · 10⁹) | – | 110 | – |
| 0.11 | < 9 · 10⁻¹⁰ | ≈ 0.02 | 330 | – |
| 0.12 | 3.3 · 10⁻⁹ | 0.12 | 1.0 · 10³ | 8 · 10³ |
| 0.14 | 1.8 · 10⁻⁸ | 0.67 | 5.7 · 10³ | 8 · 10³ |
| 0.15 | 5.5 · 10⁻⁸ | 2 | 1.2 · 10⁴ | 6 · 10³ |
| 0.18 | 6.9 · 10⁻⁷ | 25 | 6.8 · 10⁴ | 2.7 · 10³ |
| 0.24 | 2.7 · 10⁻⁵ | 10³ | 6.2 · 10⁵ | 6 · 10² |

- p(Φ) against the Chernoff model (code not included): 3.3 · 10⁻⁹ against 4.1 · 10⁻⁸ (Φ = 0.12), 1.8 · 10⁻⁸ against
  4.4 · 10⁻⁷ (0.14), 1.4 · 10⁻⁷ against 3.4 · 10⁻⁶ (0.16), 2.7 · 10⁻⁶ against 8.7 · 10⁻⁵ (0.20). The model overstates
  the flat pool by a factor of 12 to 30 (CHECKED, 3.4 · 10⁹ evaluations).
- Pool sizes at v = 1103, row sum 1 (necklaces = C(58, 29) · p/58): 1.7 · 10⁶ (Φ = 0.12), 9.3 · 10⁶ (0.14), 7 · 10⁷
  (0.16), 3.6 · 10⁸ (0.18), 1.4 · 10⁹ (0.20). Capture–recapture of the SA at v = 487 (Φ ≤ 0.14) gives 1.7 to 1.8 · 10⁶
  necklaces per role against 2.2 to 3.0 · 10⁶ from REM; at v = 1103 and Φ ≤ 0.11 it gives about 3 · 10⁵ against at most
  1.7 · 10⁵ from REM. The SA covers the pool about uniformly (factor 2).
- A further measurement (v = 1103, row sum 37): greedy swap-polish from random starts ends at Φ quantiles (min, 10 %,
  50 %, 90 %) = (0.11, 0.12, 0.15, 0.18) at about 0.1 s per sample in NumPy, and least-squares phase decoding gives
  (0.09, 0.12, 0.14, 0.17). This agrees with the REM table: p(Φ ≤ 0.15) ≈ 5 · 10⁻⁸ per random vector and a polished
  median of 0.15. The C annealer reaches the same range at 10³ to 10⁵ distinct necklaces/s (10² to 10⁴ times the NumPy
  rate), because each recorded sample is a swap-neighbour of an annealed state.
- `search/lowe_enum.c` cannot be the generator at 58 cosets. It visits every leaf of the sign tree (the energy is tested
  only at the leaf, in `leaf`; reachability pruning acts only on the row sum), that is C(58, 29) = 3.0 · 10¹⁶ leaves for
  row sum 1 (CHECKED by reading the code). Exhaustive pools below Φ = 0.10 are tiny (at most 1.7 · 10⁵ necklaces by the
  REM bound) and hold no solution: with 10⁵ necklaces per role the model gives E ≈ 10⁻⁴ solution classes.

## 3. A solution at v = 487 (CHECKED)

**Structure.** ℤ₄₈₇, U = ⟨41⟩ (order 9, 54 cosets and {0}, 27 classes), row sums (1, 17, 17, 37) (the triage pattern;
the other pattern, (1, 19, 19, 35), was not used).

**Pools.** SA at Φ ≤ 0.14, τ = 0.4, K = 10⁴: 1.0 · 10⁶ necklaces for each of the row sums 1, 17 and 37 (generation
391 + 435 + 392 = 1219 CPU-s). Lists: all four roles rotated (27 rotations), 2.7 · 10⁷ entries each; B and C come from
the same pool (n = 10⁶ each).

**Join.** `rj` with 9 merge classes, buckets in order of decreasing density. The first and only solution appeared at
bucket 1238 (1239 buckets, 2283 CPU-s, 1.84 s per bucket, 2.0 · 10⁷ pairs per bucket).

**Solution.** The SOL line `675370 24 190686 1 724700 12 91673 13` (necklace index, rotation) was decoded by `decode.py`
into `sol487.json`: four sequences of length 487 as orbit signs (coset index j of ρʲU, ρ = 3) and as ± strings. Checks
in exact integer arithmetic:

- row sums (1, 17, 17, 37) and sum of squares 1948 = 4 · 487;
- Σᵢ PAFᵢ(d) = 0 for all d ≠ 0, by direct O(v²) summation and by FFT in `gsverify.py`;
- each sequence is constant on the orbits of ⟨41⟩;
- the Goethals–Seidel array of the four circulants (order 1948) has H Hᵀ = 1948 I in int64 and entries ±1.

The flatness of the four blocks is Φ = 0.129, 0.138, 0.114, 0.121 (published for 487: 0.14, 0.14, 0.12, 0.12). The order
1948 was already known (Alpöge et al.); this is an independent re-derivation. For v = 311 no solution was attempted at
the join level (the 311 lists that fit in 3.5 GB hold E ≈ 10 tuples); 311 was used for the model calibration only.

## 4. Calibration of the solution-count model

The Gaussian-lattice estimate (full covariance of the rotation-symmetrised lists, `est.py`) is compared with counts of
partial solutions (`rj2` with `RJ_NC` = number of constrained classes; counts = hash matches verified on NC classes over
the first buckets):

| v | NC (classes) | found / predicted |
|---|---|---|
| 487 | 17, 20, 22, 23 | 0.185, 0.123, 0.093, 0.033 (3 matches; 0 of 10 predicted at 24) |
| 311 | 18, 21, 24 | 0.214, 0.152, 0.159 (36 matches) |
| 1103 | 14, 16, 18, 20 | 0.320, 0.264, 0.223, 0 of 14 predicted |
| 487, all 26 classes | 26 | 1 found / 11.7 predicted over 1400 buckets = 0.085 |

The ratio decays by about 0.9 per constrained class (the SA lists are ball-like, platykurtic), so the Gaussian formula
overstates the number of solutions by about 10 to 15 at 26 to 28 classes. The model uses 0.223 · 0.9^(NC − 18). Pair
counts per bucket: Gaussian/actual = 0.66 (487), 0.69 (311), 0.73 (1103). The full-covariance estimate is 10⁰·⁵ larger
than the independent-class estimate (487, 311, 1103).

Join constants (`rj`, one core, lists of 2.7 · 10⁷ entries): 1.8 · 10⁷ pairs/s (hash join of the two pair lists of a
bucket), 5.8 · 10⁻⁹ s per entry and bucket for the scan of the four sorted lists (0.63 s), 1.6 · 10⁻⁷ s per entry for
the index build, 20 bytes per entry and 16 bytes per pair.

Model against run at 487 (`model.py` with these constants, memory 3.3 GB): E = 1.7 · 10³ solution tuples (63 classes),
2257 buckets and 5.6 · 10³ CPU-s expected for the first solution, 1381 CPU-s of generation. Observed: 1239 buckets,
2283 CPU-s, 1219 CPU-s (compatible with a Poisson count: P(N ≤ 1238) ≈ 40 %). The uncalibrated Gaussian model would have
predicted about 120 buckets.

## 5. Prediction for v = 1103, U = ⟨17⟩ (order 19), row sums (1, 37, 39, 39)

Row-sum patterns (CHECKED, `gsflat.patterns`): r = f₀ + 19 S with S even gives |r| ∈ {1, 37, 39, 75, …}.
4412 = r₁² + r₂² + r₃² + r₄² has the single solution (1, 37, 39, 39) (f₀ = +1, S = 0; f₀ = −1, S = 2; f₀ = +1, S = 2).
U of order 29 (38 cosets) or 551 admits no pattern. The pools are C(58, 29) = 3.0 · 10¹⁶, C(58, 30) = 2.9 · 10¹⁶,
2.9 · 10¹⁶ and 2.9 · 10¹⁶.

Model optimum over (Φ threshold, pool fraction f ≤ 0.9) subject to E_classes ≥ 1.5 (HEURISTIC; p(Φ), SA rates and join
constants measured; E calibrated as above). "Partition" is an optimistic bound with sub-lists selected by exact class
values, not implemented; "chunk" is what `rj` would do with arbitrary chunks.

| memory per process | best Φ | necklaces/role | E_classes | generation CPU-s | join CPU-s | total CPU-s |
|---|---|---|---|---|---|---|
| 3.3 GB, partition bound | 0.17 | 1.2 · 10⁸ | 80 | 2 · 10⁴ | 1.5 · 10¹⁰ | 1.5 · 10¹⁰ |
| 3.3 GB, chunk | 0.16 | 5 · 10⁷ | 7 | 1.8 · 10⁴ | 1.0 · 10¹² | 1.0 · 10¹² |
| 64 GB | 0.16 | 4 · 10⁷ | 3 | 10⁴ | 1.3 · 10⁹ | 1.3 · 10⁹ |
| 256 GB | 0.16 | 7 · 10⁷ | 20 | 3 · 10⁴ | 8 · 10⁷ | 8 · 10⁷ |
| 2 TB | 0.18 | 3 · 10⁸ | 2 · 10³ | 5 · 10⁴ | 1.6 · 10⁷ | 1.6 · 10⁷ |

For comparison, `search/triage.py` gives 10¹⁵·⁵ list operations = 10⁸·² CPU-s at 1.8 · 10⁷ operations/s. Two cores for
48 hours (3.5 · 10⁵ CPU-s) fall short by a factor of 4 · 10⁴ (optimistic bound) to 3 · 10⁶ (chunking), and by 50 even
with 2 TB. Every row of the table needs lists of at least 10⁹ entries per role (E_classes ≥ 1 needs about 2 · 10⁷ flat
necklaces per role, 6 · 10⁸ entries with rotations, 12 GB of keys per list).

## 6. Other structures of FRONTIER.md

Prime v only; the full results are in `table.md` and `table.json`. The model is the same, with the pattern that has the
largest minimal pool. Entries are log₁₀ CPU-seconds for memory 3.3 GB (partition / chunk), 64 GB, 256 GB and 2 TB.

| order | v | \|U\| | classes | pattern | triage log₁₀ E | triage log₁₀ W | 3.3 GB | 64 GB | 256 GB | 2 TB |
|---|---|---|---|---|---|---|---|---|---|---|
| 2956 | 739 | 9 | 41 | (1, 19, 35, 37) | 33.5 | 22.7 | 16.9 / 20.9 | 18.6 | 18.0 | 17.0 |
| 4412 | 1103 | 19 | 29 | (1, 37, 39, 39) | 21.1 | 15.5 | 10.2 / 12.0 | 9.1 | 7.9 | 7.2 |
| 4892 | 1223 | 13 | 47 | (25, 27, 27, 53) | 33.7 | 32.4 | > 40 | > 40 | 30.6 | 29.5 |
| 5164 | 1291 | 15 | 43 | (1, 29, 29, 59) | 30.2 | 29.3 | > 40 | 27.7 | 26.5 | 25.2 |
| 6668 | 1667 | 17 | 49 | (1, 33, 33, 67) | 32.0 | 37.5 | > 40 | > 40 | > 40 | > 40 |
| 7148 | 1787 | 19 | 47 | (1, 1, 39, 75) | 30.0 | 36.1 | > 40 | > 40 | > 40 | 33.0 |
| 11884 | 2971 | 55 | 27 | (1, 1, 1, 109) | 14.0 | 19.2 | 14.4 / 17.4 | 14.4 | 13.0 | 11.2 |

The triage log₁₀ W is in list operations; subtracting 7.3 (1.8 · 10⁷ operations/s) gives the CPU-s scale of the other
columns: 1103 → 8.2, 2971 → 11.9, 739 → 15.4, 1291 → 22.0, 1223 → 25.1.

- At 2 TB the flat pool is cheaper than the plain estimate by about 1 order for 1103 (7.2 against 8.2) and by about 0.7
  for 2971 (11.2 against 11.9), and more expensive for 739 (17.0 against 15.4). For 1223 and 1291 the flat model is 3 to
  4.5 orders worse than the plain formula (HEURISTIC; the triage formula assumes a selection gain (L/N)³ that was not
  measured here). In no case is the gain close to the 3 to 4 orders that the earlier estimate of 10¹²·⁹ list operations
  suggested.
- Composite v (955 = 5 · 191, 1257 = 3 · 419, 1713 = 3 · 571, 1893 = 3 · 631, 2215 = 5 · 443) were not modelled. Their
  plain log₁₀ E is 2.9, 12.7, 5.4, 4.4 and 1.9, and the flat sub-ensemble keeps a fraction of order 10⁻¹⁵ to 10⁻²⁵ of
  the solutions, so none is expected (HEURISTIC).
- For 739 and for 1223 to 1787 the REM table has fewer than 10 hits below Φ = 0.19 to 0.22, so p(Φ) there is a power-law
  extrapolation.

## 7. Programs

| File | Role |
|---|---|
| `gsflat.py` | orbit and coset tables, Gaussian periods, exact PAF, row-sum options, patterns (r₁, …, r₄) for (v, \|U\|) |
| `mktab.py` | writes `tab_<v>.txt` (v, m, C, \|U\|, v mod 4, η_re[m], η_im[m]) for the C programs |
| `flat.c` → `flat2` | `gen2`: simulated-annealing generator of flat necklaces; `rem`: REM sampling of p(Φ); `keys`; `phi` (uint64 masks, m ≤ 64) |
| `flatbig.c` → `flatbig32`, `flatbig64` | `rem` and `gen2` with 128-bit masks (m ≤ 128, C ≤ 32 or 64) for the other structures of FRONTIER.md |
| `flatbig_base.c` | the 128-bit-mask version of `flat.c` with all its commands (`gen`, `keys`, `phi`, `rem`, `gen2`); `flatbig.c` is the reduced form that `tablegen.py` uses |
| `genpool.py` | restart-safe pool generator: chunks `<pooldir>/c<i>.bin` and `c<i>.done`; stops at TARGET distinct necklaces |
| `mergepool.py` | merges finished chunks into `<pooldir>/pool.bin` (sorted unique canonical masks, uint64) |
| `buildrj.py` | builds the join input of `rj` (int8 key vectors of the chosen necklaces) and `.map.npz` |
| `buildlists.py` | the same for the generic four-list k-tree `../kt4.c`: int16 key vectors, every necklace of the roles B, C, D taken through all C rotations, target −3 |
| `rj.c` → `rj` | rotation-virtual four-list k-tree |
| `rj2.c` → `rj2` | `rj` with the environment variable `RJ_NC`, the number of constrained classes, for model calibration |
| `runjoin.sh` | restart-safe driver of `rj` over bucket ranges (`.done` markers); `PAR` and `NPAR` split the work between processes |
| `test_rj.py` | writes `w/rjtest.bin`, an input of `rj` with a planted solution (3000 random necklaces per role, C = 27, T = −3) |
| `decode.py`, `gsverify.py` | map a SOL line back to orbit signs and full ±1 sequences; exact PAF check, Goethals–Seidel array, H Hᵀ = 4v I |
| `est.py`, `estpools.py`, `buckmodel.py` | Gaussian estimate of the number of solutions in given lists and per bucket (calibration) |
| `lattice.py` | rank and index (product of the pivots) of the integer lattice spanned by the differences of the first 400 key vectors of a key file on its first NC classes: `python3 lattice.py KEYS.bin C NC` |
| `model.py`, `predict.py`, `tablegen.py`, `table.py` | the cost model with measured constants; the scan over (Φ threshold, pool fraction); the per-structure measurements `tg_<v>_<u>.json`; their tabulation in `table.md` and `table.json` |

Data: `sol487.json` (the verified quadruple of order 1948), `tab_<v>.txt` (coset tables), `tg_<v>_<u>.json` (measured
p(Φ) and generation rates), `table.md` and `table.json` (the predictions of section 6).

## 8. Pipeline and how to resume

Build, in this directory:

```
gcc -O3 -march=native -o flat2 flat.c -lm
gcc -O3 -march=native -o rj rj.c -lm
gcc -O3 -march=native -o rj2 rj2.c -lm
gcc -O3 -march=native -DLW=32 -o flatbig32 flatbig.c -lm       # LW = lanes >= C, a multiple of 8 (see the header)
gcc -O3 -march=native -DLW=64 -o flatbig64 flatbig.c -lm
```

**1. Tables** (cheap).

```
python3 mktab.py v |U| [generator]                  # writes tab_<v>.txt
```

**2. Pools**, one per distinct row sum r (f0 and npos from `gsflat.rowsum_options`).

```
python3 genpool.py POOLDIR tab f0 npos THR TAU K TARGET [CHUNK=2e8] [SEED0=1] [MAXCHUNKS]
python3 genpool.py p487/r1 tab_487.txt 1 27 0.14 0.4 10000 1e6 3e8 1              # v = 487
python3 mergepool.py POOLDIR                        # makes POOLDIR/pool.bin
```

Running the same command again resumes: finished chunks are those with a `.done` file, and a chunk without one is
regenerated.

**3. Join input.** All four roles rotated; T = −(v mod 4) = −3 for v ≡ 3 (mod 4).

```
python3 buildrj.py OUT.bin tab C -3 n:poolA:f0:1 n:poolB:f0:1 n:poolC:f0:1 n:poolD:f0:1 SEED
```

**4. Join.** Run the driver twice, with PAR = 0 and PAR = 1 and NPAR = 2, to use two cores. Running the same command
again skips the finished ranges. A solution line is appended to `jdir/sol_all.txt`; `touch jdir/STOP` ends the drivers.

```
./runjoin.sh $PWD/OUT.bin 9 0 NBUCKETS 200 $PWD/jdir 6e7 PAR NPAR
```

**5. Verification.** Exact: the PAF sum is zero at every shift, then the Goethals–Seidel array and H Hᵀ = 4v I.

```
python3 decode.py v |U| OUT.bin.map.npz F0A,F0B,F0C,F0D "SOL x0 k0 x1 k1 x2 k2 x3 k3" out.json
```

**Memory.** `rj` needs about 20 bytes per list entry (entries = C · necklaces per role) and 16 bytes per pair of one
bucket. With the 3.5 GB address-space cap that means about 2.7 · 10⁷ entries per list (n = 10⁶ necklaces at C = 27).
