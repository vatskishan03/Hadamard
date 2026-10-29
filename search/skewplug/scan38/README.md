# Re-scan of the K-symmetric real census shapes

The programs behind [../../NOTES.md](../../NOTES.md) §91:

- the enumeration of every K-symmetric real skew shape of the census values (G = K × F_q, q ≡ 3 (mod 4) a prime power
  dividing v, H = H_K × squares);
- the all-slot Gaussian-estimate scan, which extends [../est_extend.py](../est_extend.py) to all shapes (partial
  coverage, resumable, with controls and tables);
- the exact all-energy row-sum scan;
- the production rotation of [../exh38.py](../exh38.py), the engine that found [skew/1540](../../../skew/1540).

Everything here runs in `search/skewplug` (`python3 scan38/NAME.py`): the scripts put `..` on the path and write
`est/`, `runs/`, `claims/` and `out/` there, after the build of that directory ([../README.md](../README.md#build)).
The terms (shape, slot, allocation, class) are explained in [../README.md](../README.md#method-and-terms).

## Programs

### Enumeration and annotation

| File | Contents |
|---|---|
| `enum38.py` | the enumeration: K every abelian group of order \|K\| (cyclic; non-cyclic by Sylow types), H_K every odd-order subgroup of Aut(K) for cyclic K (units of ℤ_m) and of GL(2, p) up to conjugacy for ℤ_p² blocks (torus, Singer, unipotent ⋊ diagonal), torus and Singer subgroups for ℤ_pᵏ, k ≥ 3, diagonal units for the other non-cyclic parts (PARTIAL, flagged), at most 25 orbits; `python3 scan38/enum38.py [NUMAX=25] [OUT]` (106 s) |
| `decided38.py`, `prep38.py`, `flag38.py` | the shapes decided EXHAUSTIVELY in the earlier sweeps (reads `data/census_shapes15.json`, `timeouts.json`, `SHAPE_TABLE.txt`, `verdict_table_cyclic.txt` and `../../negaskew/real37/data/RESULTS.txt`) → `runs/decided38.json`; the annotation of an enumeration (decided match, window allocations, pool size) → `runs/shapes38.json`; the re-annotation after `decided38.py` |
| `pshape.py` | `sk.SShape` with explicit K permutations (matrix H_K, mixed K) |
| `awin.py` | the window-aware allocation generator and lazy samplers: equal to `sk_pat.allocations` restricted to the energy window, 100 times faster at q = 3 |
| `rowfast2.py`, `test_rowfast2.py` | the vectorised exact row-sum count of a slot (equal to `sk_pat.row_patterns_skew` in feasibility; corrected, see [Notes on the code](#notes-on-the-code)); its regression test |

### Scan

| File | Contents |
|---|---|
| `scan38.py` | the scan of one shape: time-budgeted and cumulative over invocations, resumable; pools up to 60,000 canonical slots in a fixed shuffled order, larger ones sampled uniformly with the last 40% of the budget on all splits of the 30 best allocations; `WMODE=lo`: q = 3 window E ∈ [N, 1.6N]; per-slot time limits |
| `sched38.py` | the queue worker (claims in `claims/PASS/idx`; at most 3 workers) |
| `summ38.py`, `table38.py`, `mkq2.py` | the summary of the slot files (`runs/summary38.json`, `runs/SUMMARY38.txt`), the final tables (`runs/TABLE38.txt`: scanned shapes ranked by the margin min(logE, logρ + 8.5) with coverage, feasible slots, tbits and the bias-corrected margins adjMed/adjP90; `TABLE38_decided.txt`; `table38.json`), and the queue of the deeper passes |
| `rowscan38.py` | the exact all-energy row-sum scan: is there any slot (any energy, any m-split) with a row-sum pattern? Status ALL_INFEASIBLE (EXACT-dead), FEASIBLE (a slot) or PARTIAL (time cap); appends to `runs/rowscan38.jsonl` |

### Controls

| File | Contents |
|---|---|
| `ctl473.py` | 473, a solved shape with a multiplier |
| `ctl539.py` | 539 = ℤ₇ × ℤ₇ × F₁₁ with the scalars ⟨(2, 2)⟩ (the shape of [skew/2156](../../../skew/2156)) |
| `ctl309.py`, `q3probe.py` | 309 (q = 3, low window); `q3probe.py`: the row-sum-feasible share per energy at q = 3 |
| `ctl119.py` | 119 = ℤ₁₇ × F₇ (solved, H_K = 1) |
| `ctlrigid.py` | the rigid shapes 253, 285, 345, 465, 391, 275 (solved, H_K = 1) |

### Production of the explicit-list sweep

| File | Contents |
|---|---|
| `mkqueue.py`, `mkqueue2.py`, `mklow.py` | from `runs/table38.json` and the slot files: per shape the best slots by completions per class (`runs/shapes_scanned.json`); the shapes ranked by the estimated solutions in slots of cost ≤ CAP; the low-energy slot files `est/low<idx>_v<v>.jsonl` (all slots of E below the least energy of the scan's file, energy ascending, dummy estimates) |
| `rotate.py`, `rotate_launch.sh` | the time-sliced rotation of `../exh38.py` over the shapes of a queue (a pass is `[MAXCLASS, base seconds]`; from the second pass the slice follows the measured yield of new classes per second) and its launcher |
| `run_exh.sh`, `supervise.sh`, `run_job.sh` | a single detached `exh38.py` job; the supervisor that restarts `exh38.py` until ALL SLOTS FINISHED or a deadline; `prod_dedup.py` as a job |
| `summarize.py`, `fom.py`, `cmp_units.py` | the per-job table (classes, completions, hits); the figures of merit of a slot file; the check of the unit reduction of the skew block |

### Independent check

| File | Contents |
|---|---|
| `indep1540.py` | an independent check of [skew/1540](../../../skew/1540): its own group arithmetic and Goethals–Seidel array from `certificate.json`, no repository code; `python3 scan38/indep1540.py ../../skew/1540/certificate.json [H1540.csv.gz [sequences.txt]]` |

## Data

| File | Contents |
|---|---|
| `data/shapes38.json` | the annotated enumeration of the scan: 554 shapes; the shape indices of every table and queue refer to it |
| `data/TABLE38.txt`, `data/SUMMARY38.txt` | the results of the full scan, as the scan wrote them (before the correction of `rowfast2.py`) |
| `data/TABLE38_decided.txt` | the 313 decided or EXACT-dead shapes of the enumeration |
| `data/rowscan38.jsonl` | the exact all-energy row-sum scans of the full run |
| `data/queue0.json`, `queue1.json`, `queueq3.json`, `queue2.json`, `queue2b.json`, `queue1q3.json`, `queue3.json`, `queuecal.json` | queues of shape indices for `sched38.py` (see [Reproduction of search/NOTES.md §91](#reproduction-of-searchnotesmd-91)) |
| `data/queue/shapes.json` | the 36 shapes of the production rotation, with their slot files |
| `data/ctl119_rows.json`, `ctl309_rows.json`, `ctl473_rows.json`, `ctl539_rows.json`, `ctlrigid_rows.json` | the shapes of the controls |
| `data/q3probe_309.txt`, `data/dead_ids.json` | the output of `q3probe.py` for 309; an empty list |
| `data/census_shapes15.json`, `timeouts.json`, `SHAPE_TABLE.txt`, `verdict_table_cyclic.txt` | inputs of `decided38.py`: the earlier sweeps |

Source tags in `decided38.py` and in the recorded tables name the sweeps behind a verdict. `r34` to `r37` are the
sweeps documented in §88–90; `decided38.py` reads them from `verdict_table_cyclic.txt` (`r34`), `SHAPE_TABLE.txt`
(`r35`), `census_shapes15.json` with `timeouts.json` (`r36`) and `../../negaskew/real37/data/RESULTS.txt` (`r37`).
`r38` tags the scan itself (§91).

The indices of the scan (`IDX`, `data/queue*.json`) refer to `data/shapes38.json`. An enumeration made now has two
shapes fewer (385 is settled) and different indices after them, so the scan and its tables use the stored list.

Python 3 with NumPy and SymPy.

## Reproduction of search/NOTES.md §91

Run in `search/skewplug`. The result column is the output expected, trimmed; times are wall times on one core of a
loaded four-core machine.

### Enumeration

EXACT up to the stated coverage.

| Command | Time | Result |
|---|---|---|
| `python3 scan38/enum38.py 25 runs/shapes38_raw.json` | 106 s | `552 shapes (pre feasibility) Counter({3: 102, 15: 102, 9: 63, 7: 57, 5: 45, 25: 45, 21: 38, 11: 22, 17: 22, 13: 20, 19: 19, 23: 17})`; the stored list has 554 shapes: the two shapes of 385 are no longer enumerated since 385 is settled (ν = 15: 104 → 102); the order is otherwise the same |
| `python3 scan38/decided38.py` | 1 s | `642 decided shape records Counter({'r36': 452, 'r34': 80, 'r37': 64, 'r35': 46}) Counter({'cyc': 636, 'nc': 6})` (`runs/decided38.json`) |
| `python3 scan38/prep38.py runs/shapes38_raw.json runs/shapes38.json` | 505 s | `decided 296 dead identity 47 err 0`; every annotation equals the one of `data/shapes38.json` |

The scan covered 641 undecided census values, 339 pairs (v, q) with q ≡ 3 (mod 4) a prime power dividing v and a
coarse identity solution (q + 1) N_w − 4N₋ = 4|K|, and 554 shapes. Of these, 296 were decided by the earlier sweeps
(452 records of the real sweep minus its timeouts, 64 of `real37`, 46 of the shape table, 80 of the cyclic verdict
table; matched by (v, q, |K|, the unit group)), 47 of them also die by the orbit weights alone, and 258 are undecided.

### Controls of the scan

The numbers are those of [../../NOTES.md](../../NOTES.md) §90–91. The scan files are `est/S<idx>_….ext.jsonl`, written
with the `outdir` and `idx` of the control, so delete them before a real scan.

| Command | Time | Result |
|---|---|---|
| `python3 scan38/ctl473.py` | 12 s | 473 = ℤ₄₃ × F₁₁ ⟨6⟩: 1235 slots, 161 row-sum feasible, best logE +0.67 / logρ −7.04 at (0,6,7,7)/(0,3,7,7), top-40 sum 54.5; the two hit slots of [skew/1892](../../../skew/1892) rank 6 (+0.41) and 41 (−0.27) of 161 |
| `python3 scan38/ctlrigid.py` | 35 s | 253: +5.78 / −0.57, 285: +6.81 / −2.11, 345: +6.30 / −2.33, 465: +5.08 / −1.95, 391: +6.12 / −3.01, 275 = ℤ₂₅ × F₁₁ (H_K = 1): +8.74 / −8.68 (649 slots) |
| `python3 scan38/ctl309.py` | 28 s | q = 3, 309 = ℤ₁₀₃ × F₃ ⟨13⟩, `WMODE=lo`: 68,866 slots, 426 row-sum feasible, best +1.34 / −5.67 at (34,34,35,51)/(0,17,17,17), top-40 sum 219 |
| `python3 scan38/ctl539.py` | 22 s | 539 = ℤ₇ × ℤ₇ × F₁₁ with the scalars ⟨(2,2)⟩ (the shape of [skew/2156](../../../skew/2156)), top window, `est/S0000_v539_*.ext.jsonl`: 2265 slots, 266 row-sum feasible, best logE +0.94 / logρ −7.87 at allocation (6,4,6,7), the hit slot (0,6,7,7)/(0,3,4,4) third at +0.79 / −7.75 (best margin min(logE, logρ + 8.5) = +0.75, at the hit slot; the first two slots +0.63, +0.60), top-40 sum 107.2 (HEURISTIC; `est_extend.py` on this shape ran 29 slots in 240 s) |

### A scan, its summary and the tables

One shape, then the queue worker on three shapes of `data/shapes38.json`: 781 = ℤ₇₁ × F₁₁ with ⟨5⟩ is index 77, 917 is
index 101, 575 over ℤ₂₅ with H_K = 1 is index 33, and 1539 over ℤ₃⁴ is index 297.

| Command | Time | Result |
|---|---|---|
| `python3 scan38/scan38.py scan38/data/shapes38.json 77 150` | 8 s | COMPLETE, 1735 slots (pool bound 2204), 61 row-sum feasible |
| `echo '[77, 101, 33]' > runs/q_test.json; python3 scan38/sched38.py scan38/data/shapes38.json runs/q_test.json 60 TEST 5 60` | 11 s | 781 (done), 917: COMPLETE 161 slots, 575: COMPLETE 9 slots |
| `python3 scan38/summ38.py` | | 3 shapes summarised; `runs/SUMMARY38.txt` |
| `python3 scan38/table38.py` | | `runs/TABLE38.txt`: 575 (ℤ₂₅, H_K = 1): best logE +0.37, logρ −7.31, margin +0.37, tbits 24.0; 781: −1.88 / −8.68, margin −1.88, tbits 12.9, top-40 sum 0.12; 917: no feasible slot in the window |

The section "C. counts" of `runs/TABLE38.txt` reads: SCANNED 3 (here), DECIDED 249, DECIDED+identity-dead 47,
DECIDED(prose) 8, EXACT_DEAD: N_w = 1 lemma 2 (1139, 1207), N_w = 2 lemma 1 (1081), r34 big-big filter 2 (361, 529), no
row-sum pattern in any slot 1 (1539 here), NOT_SCANNED 241 (the full scan: 241 scanned, 9 EXACT-dead).

`data/TABLE38.txt` is the table of the full scan, as the scan wrote it, that is before the correction of `rowfast2.py`.
Rows with few or no feasible slots can be wrong for the shapes that shared a worker process with others: 575 over ℤ₂₅
and ℤ₅², 385 over ℤ₅₅ × F₇, 627, 957, 1539 over GF(27); the corrected values are in §91. The full scan has 241 scanned
shapes (44 COMPLETE, 197 partial). The ranking by margin has 575 over ℤ₂₅ at +0.37 with tbits 24, 245 = ℤ₃₅ × F₇ at
−0.52, 385 at −0.91, 781 at −1.88, 1645 at −2.0, 1057 at −2.25, 627 at −2.55, 2485 at −3.0 and 1267 at −3.1; the
controls give +0.67 (473) and +0.63 (539). `data/TABLE38_decided.txt` is the table of the 313 decided or EXACT-dead
shapes.

A full scan is `sched38.py` over `queue0.json` (undecided ν ≤ 15, q ≠ 3), `queue1.json` (ν ≥ 17), `WMODE=lo` over
`queueq3.json`, then the deep passes `queue2.json` and `queue2b.json` with larger cumulative budgets: 90 s per shape,
three workers, about two hours on a free machine.

### The exact all-energy row-sum scan

| Command | Time | Result |
|---|---|---|
| `python3 scan38/rowscan38.py scan38/data/shapes38.json 297 100` | 1 s | 1539 = ℤ₃⁴ × F₁₉ ⟨5⟩: ALL_INFEASIBLE (92 slots): EXACT-dead |
| `python3 scan38/rowscan38.py scan38/data/shapes38.json 101 100` | 10 s | 917 = ℤ₁₃₁ × F₇ ⟨13⟩: FEASIBLE after 7894 slots, (0,0,1,117)/(0,0,1,104) at E = 118 (alive although the top window has no feasible slot) |

`data/rowscan38.jsonl` has the scans of the full run: ALL_INFEASIBLE for 1183 (79,895 slots), 1539 over ℤ₃⁴, 1603
(1,114,897), 1859, 1919 (388), 2375, 2401 with order 57 (75,747) and 171. The status of 1793 and 2189 recorded there
('FEASIBLE at other energies') was later corrected: their orbit weights cannot realise the identity at all (§91).

### The correction of `rowfast2.py`

An EXACT regression test; see [Notes on the code](#notes-on-the-code).

| Command | Time | Result |
|---|---|---|
| `python3 scan38/test_rowfast2.py` | 4 s | 51 slots over nine shape constructions of 575, 1081, 1139, 1207 in one process: `TOTAL slots 51 mismatches 0` |
| `python3 scan38/test_rowfast2.py PATH` | 5 s | with the path of the uncorrected module (not included): mismatches 1 (1139 after 1081: the block classes of 1081 served 1139) |

### The explicit-list sweep and its rotation

`../exh38.py`; the replay of the hit slot is also in [../README.md](../README.md#order-1540).

```
cp data/repro385/est_hit_slot.jsonl .
PROD_CONTINUE=1 python3 -u exh38.py hit385 600 "[35]" 11 "[[11]]" est_hit_slot.jsonl 1 file skf
python3 cyclic.py out/skew_1540_hit385_s0_<time>_<i>.json OUT.json
python3 scan38/indep1540.py ../../skew/1540/certificate.json ../../skew/1540/H1540.csv.gz ../../skew/1540/sequences.txt
echo '[{"idx": 1, "v": 385, "q": 11, "K": [35], "hu": [[11]], "est": "est_hit_slot.jsonl"}]' > runs/rot_queue.json
(cd scan38 && ROT_QUEUE=../runs/rot_queue.json python3 -u rotate.py $(( $(date +%s) + 120 )) '[[1e30, 40]]' 1)
python3 scan38/summarize.py
scan38/supervise.sh r2 $(( $(date +%s) + 90 )) "[35]" 11 "[[11]]" est_hit_slot.jsonl 1 file skf
```

1. `exh38.py` takes 19 s and prints two HIT lines with sha256
   `9ea5b315360207c081f008440ba306014d69157810720301c9d81918c1e49dbd` and
   `1b4d57ecfdfdc1bb4e5daf2f7b2ff88a1e8ebf7ddc0605dfb439fab97a3445a8`, 100 classes, nok 2, and `ALL SLOTS FINISHED`.
2. The first HIT file, turned into the cyclic form, equals `skew/1540/certificate.json` without `sha256`.
3. `indep1540.py` takes 14 s (4 s CPU). It is an independent rebuild: subgroup of order 15, row sums [1, −7, 31, 23],
   sum of squares 1540, PAF 0 at all non-zero shifts, H Hᵀ == 1540 I and H + Hᵀ == 2I (float64 integers, exact), CSV
   sha256 `f98a3775…`; the stored matrix equals the rebuild (EXACT).
4. The rotation takes 7 s and prints `x1_385 v=385 slice 40s wall 6s: +100 classes (total 100), +2 completions
   (total 2), yield 15.47/s, finished=True`.
5. `summarize.py` prints `x1_385  slots 1  finish 1  classes 100  compl 2  hits 2`.
6. `supervise.sh` restarts `exh38.py` until `ALL SLOTS FINISHED` (6 s here).

The production runs of §91 are resumable through `runs/NAME.slots.jsonl` and `runs/NAME.seen`; their slot files are
copies of the scan's `est/S*.ext.jsonl`. They are relaunched as follows (`DEADLINE` is an epoch time).

- 245 = ℤ₃₅ × F₇ with ⟨11⟩:
  `EX_MAXCLASS=30000 scan38/run_exh.sh e245 14400 "[35]" 7 "[[11]]" est/c245_all.jsonl 4000 rho skf`.
- The rest of 385:
  `EX_MAXCLASS=1e30 scan38/run_exh.sh e385 14400 "[35]" 11 "[[11]]" est/c385_s1.jsonl 3000 rho skf`.
- The low-energy slots of 245: `scan38/mklow.py IDX`, then
  `EX_SLOT_SECS=120 scan38/supervise.sh e245L DEADLINE "[35]" 7 "[[11]]" est/low4_v245.jsonl 200000 rho skf`.
- The rotation over the other shapes:
  `scan38/rotate_launch.sh A DEADLINE '[[3000,150],[30000,300],[300000,600]]'
  '77,67,337,15,550,280,203,47,256,110,174'`.
  `data/queue/shapes.json` lists the 36 shapes with their slot files `est/q<idx>_v<v>.jsonl`, to be copied from the scan
  files. The `EX_*` variables and `runs/ex_env.json` are documented in `../exh38.py`. Stop a job with
  `touch runs/NAME.STOP` or `runs/rotate.STOP`.

## What the numbers cover

The scan is HEURISTIC: the slot weights come from Gaussian estimates that are low by 3 to 4 decades at the hits of 473
and 539, high by at least 2 decades for q = 3, and ignorant of the exact lemmas. A partially covered shape is excluded
only when its bias-corrected margin stays below 0 after the deep pass.

The exact claims are the EXACT-dead shapes (the N_w = 1 and N_w = 2 lemmas, the big–big filter at H_K = 1, the
all-energy row-sum scan, the identity by orbit weights). The classes of `exh38.py` slots cover 60 to 90% of each raw
block list, so a finished slot is HEURISTIC-complete.

## Notes on the code

- The scripts take their data from `data/` and write to `runs/` and `est/` of `search/skewplug`; they add `..` to the
  path.
- `enum38.py` takes the repository root from the environment variable `REPO` (default: three levels up).
- `decided38.py` reads its inputs from `data/`.
- The rotation scripts take `search/skewplug` as their working directory, log to `runs/` and read the queue from
  `data/queue/shapes.json` (`ROT_QUEUE`).

**A bug in `rowfast2.py`, fixed.** The module kept its caches (block classes, tau tables) in module-level dictionaries
keyed by `id(sh)`. A worker that scanned several shapes in one process reused the id of a freed shape, so slots were
wrongly marked as having no row-sum pattern: 575 over ℤ₂₅: 7 of 9 slots, 575 over ℤ₅²: 7 of 9, 385 over ℤ₅₅ × F₇: 627
of 1214, 627: 480 of 2103, 957: 1080 of 7259, 1539 over GF(27): 45 of 100. The caches now live on the shape object.
The estimates of the 575 shapes moved (ℤ₅²: −0.47 / −7.38 to +0.17 / −7.46), and no margin of the table moved above
the others (§91). `data/TABLE38.txt` and `data/SUMMARY38.txt` are the files of the scan as written then, before the
correction, with the amendments of §91.
