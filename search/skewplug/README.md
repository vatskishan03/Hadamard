# Skew Paley plug-ins in the Goethals–Seidel array

The two-level sampler and the exhaustive engines behind the skew-Hadamard matrices of orders 476, 980, 1540, 1892, 2156
and 3124 ([../NOTES.md](../NOTES.md) §87–§91 and §93). They search skew Goethals–Seidel quadruples (A, B, C, D) over
G = K × F_q, with K an abelian group and q ≡ 3 (mod 4) a prime power, whose blocks are invariant under
H = {1} × (squares of F_q\*), or under H = H_K × squares when a multiplier group H_K acts on K: the Paley products of
[../NOTES.md](../NOTES.md) §15. Each block is s ⊗ I + t ⊗ Q + e ⊗ J with Q the Paley matrix of F_q. The first block is
of skew type (s and e odd under k ↦ −k off 0, s(0) = 1, e(0) = 0, t even); B, C and D are free. §15 solved the case of
at most five w-positions (defects). This code handles many of them (nine to seventeen for order 476).

| Order | Shape | Program | Notes | Certificate |
|---|---|---|---|---|
| 476 = 4 · 119 | ℤ₁₇ × F₇, H_K = 1 | `sk_worker.py` | §87 | [skew/476](../../skew/476) |
| 980 = 4 · 245 | ℤ₃₅ × F₇, H_K = ⟨11⟩ (group ℤ₃₅ × ℤ₇, not cyclic) | `exh38.py` | §93 | [skew/980](../../skew/980) |
| 1540 = 4 · 385 | ℤ₃₅ × F₁₁, H_K = ⟨11⟩ | `exh38.py` | §91 | [skew/1540](../../skew/1540) |
| 1892 = 4 · 473 | ℤ₄₃ × F₁₁, H_K = ⟨6⟩ | `prod_ks.py`, `prod_dedup.py` | §88, §89 | [skew/1892](../../skew/1892) |
| 2156 = 4 · 539 | ℤ₇ × ℤ₇ × F₁₁, H_K = ⟨(2, 2)⟩ | `prod_dedup.py` | §90 | [skew/2156](../../skew/2156) |
| 3124 = 4 · 781 | ℤ₇₁ × F₁₁, H_K = ⟨25⟩ (group ℤ₇₈₁, cyclic) | `exh38.py` | §93 | [skew/3124](../../skew/3124) |

The pattern pipeline in `patterns/` decides rigid shapes, in which large compression entries cannot cancel. `patterns/`
and `estimates/` hold the programs behind the verdicts of [../NOTES.md](../NOTES.md) §88–89 (see
[Reproduction of search/NOTES.md §88–89](#reproduction-of-searchnotesmd-8889)), and `scan38/` the re-scan of §91.

## Contents

1. [Method and terms](#method-and-terms)
2. [Layout](#layout)
3. [Programs](#programs)
4. [Data](#data)
5. [Build](#build)
6. [Reproducing the skew certificates](#reproducing-the-skew-certificates)
7. [The pattern pipeline](#the-pattern-pipeline)
8. [Reproduction of search/NOTES.md §88–89](#reproduction-of-searchnotesmd-8889)
9. [Status labels](#status-labels)

## Method and terms

A block X is given by three functions on K: s = x₀ − e, t = (x₁ − x₂)/2 and e = (x₁ + x₂)/2, where x₀, x₁, x₂ are the
values of X on K × {0}, on K × S and on K × N (S the squares, N the non-squares of F_q\*). With A of skew type, the
quadruple is a solution exactly when the K-level system (I), (J), (Q) of `sk.py` holds. The search has two levels.

- *Level 1* works with the compressions u = s + q · e of the blocks to K. It joins the exact list of skew blocks A with
  lists of B, C, D (sampled or complete) on the equation Σ PAF_u = 4|G| δ, the equation of the trivial character of
  F_q (equation (U)).
- *Level 2* completes a level-1 instance by the signs t off the defects (an exact meet in the middle, filtered by row
  sums and spectral bounds). A solution of both levels is certified with the builders of `verify/verify.py`.

| Term | Meaning |
|---|---|
| shape | (K, q, H_K), written `'K' q 'hunits'` on the command line: the moduli of K, the field order q (a prime or 27 = 3³) and the generators of H_K acting on every block, each a unit per component of K (cyclic K: `[[g]]`) or an integer matrix (non-cyclic K) |
| defect, w-position | an H_K-orbit O of points of K on which a block is constant on F_q\* (t = 0, e = ±1), of type p if s = 0 and of type m if s = −2e; the points of O are w-positions, so O has weight \|O\|; the skew block carries defects in pairs {O, −O} |
| defect identity | (q + 1)e₊ + (q − 3)e₋ = 4\|K\|, where e₊ and e₋ count the w-positions of type p and m; with N_w the number of w-positions and N₋ those of type m it reads (q + 1)N_w − 4N₋ = 4\|K\| (§15) |
| allocation, m-split | the numbers (n_A, n_B, n_C, n_D) of w-positions per block, and the numbers (m_A, m_B, m_C, m_D) of type m among them |
| slot, energy | an (allocation, m-split) pair; its energy E is the total n_A + n_B + n_C + n_D = N_w |
| logE, logN1, logρ | base-10 logarithms of the Gaussian estimates of the number E of solutions of a slot (not its energy), of its number N1 of level-1 instances and of ρ = E/N1, the completions per level-1 instance (`est_skew.py`; HEURISTIC) |
| row-sum pattern | an exact combination of the row sums (σ_X, ε_X) of s and e over the blocks B, C, D that a slot allows (the skew block has σ = 1, ε = 0; `sk_pat.py`) |
| level-1 instance, class | a solution of level 1; a class is an instance modulo the exact symmetries of the K-level system: reflection k ↦ −k of each block, negation and translation by the H_K-fixed points of K for B, C, D, and, in `exh38.py` for cyclic K, the units of K acting on all four blocks |
| round | one step of a production loop: draw a slot, a row-sum pattern and block lists, run level 1, send the instances to level 2; rounds are seeded by (seed base, round index) and numbered from 0 in the logs |
| hit | a level-2 completion, checked at the K level and written as a certificate |

## Layout

| Path | Contents |
|---|---|
| this directory | engines, production loops, shared shape code, certificate writers and validation programs, listed under [Programs](#programs) |
| `patterns/` | the pattern pipeline and the checks of its soundness (§88–89) |
| `estimates/` | estimates and measurements around the K-symmetric loops |
| `scan38/` | the re-scan of the K-symmetric real census shapes and the production rotation of `exh38.py` (§91); [scan38/README.md](scan38/README.md) |
| `data/` | slot files, a cost file, and the states needed to replay the hits |

## Programs

All commands run in this directory. The programs write to `runs/` (logs), `out/` (certificates), `est/` (slot files)
and `tmp/` (temporary lists), none of which is tracked; the scripts in the subdirectories put `..` on the path.

### Shape and exact conditions

| File | Contents |
|---|---|
| `sk.py` | the shape (`SShape`): orbits, the skew block, the K-level equations (I), (J), (Q), exact lists of the skew block, checks |
| `sk_pat.py` | allocations allowed by the defect identity, their splits into the two defect types, exact row-sum patterns |
| `gsk.py` | the K-level system of the real Goethals–Seidel Paley plug-ins with four free blocks (`GShape`), the exact G-level check, exact lists of blocks of a given defect weight |
| `pk.py`, `idx.py` | K-symmetric orbit data (shared with the Szöllősi plug-in searches) and the index-r plug-in classes G = K × F_q |
| `sk_filter.py` | cheap exact necessary conditions on a level-1 instance for the existence of t signs, and a score for ranking |

### Level 1 and level 2

| File | Contents |
|---|---|
| `gl1.c`, `gl1.py` | level 1: the exact list of the skew block against sampled (or explicit) lists of B, C, D, class-restricted hash join on equation (U); `gl1.py` is the ctypes wrapper of `libgl1.so` |
| `gl2.c`, `sl2io.py`, `gl2io.py` | level 2 (`sl2`): exact completion of the signs t off the defects, with the skew block's own t-orbits; `sl2io.py` and `gl2io.py` write its input and read its output |
| `pl2f2.c`, `pl2f2_fast.h`, `pl2f2_main.h`, `skf_io.py` | a faster exact level 2 (`skf2`) for cyclic K with \|K\| ≤ 64, with the same input and the same set of completions as `sl2`; the headers are included by `pl2f2.c` |

### Production loops

| File | Contents |
|---|---|
| `sk_samp.py`, `sk_worker.py` | one level-1 and level-2 round; the production loop for H_K = 1 (seeded rounds, resumable, stops at the first hit) |
| `prod_ks.py` | the production loop for shapes with a multiplier on K: slots weighted by estimates or by a cost file such as `data/cost473.jsonl`, every level-1 instance to level 2, certification on a hit |
| `prod_dedup.py` | `prod_ks.py` with symmetry deduplication of the level-1 instances (canonical keys under the reflections, negations and H_K-fixed translations of the blocks), a pattern bandit and bounded level 2; found a sixth matrix of order 1892 |
| `exh38.py` | a symmetry-reduced explicit-list sweep of the slots of a K-symmetric shape: per block, defect weight and row sum the list is generated until it covers the block space, one block per orbit of the exact block symmetries, the skew block reduced by the units of K, every row-sum pattern of a slot joined with `gl1.c` (its functions `gl1_genlist`, `gl1_hashlist`, `gl1_set_hin`), new classes to the exact level 2; about 100 classes per second against 2 to 3 for `prod_dedup.py`; found skew orders 1540, 980 and 3124 |

### Estimators

| File | Contents |
|---|---|
| `est.py` | Gaussian local-limit estimator for orbit-constant plug-in systems |
| `est_skew.py` | the Gaussian local-limit estimator of the number of level-1 instances and of solutions of a slot (allocation, m-split) of a K-symmetric shape; HEURISTIC: it ranks slots, it cannot discard a shape |
| `est_extend.py`, `est_merge.py` | the Gaussian estimate over every slot of the top energy window of a shape, and its merge into a slot file `est/NAME.jsonl` (the slot weights of the run that found skew order 2156) |

### Certificates and checks

| File | Contents |
|---|---|
| `sk_cert.py` | a K-level solution to four sequences over ℤ_m × ℤ_q, checked with the builders of `verify/verify.py`, written to `out/` |
| `cyclic.py` | the same certificate over the cyclic group ℤ_mq (the format of [../../skew/820](../../skew/820)) |
| `indep_check.py` | an independent numpy check of a group-form certificate over ℤ_m × ℤ_q |
| `bf.py`, `sk_exh.py`, `validate_skew.py` | G-level and K-level brute force, the exhaustive two-level pipeline, the validation battery |
| `val_skf.py`, `val_skf.sh` | validation of the fast level 2 `skf2` against `sl2` on level-1 instances (completions, fibre sizes, sets of solutions) and its battery of 22 slots |

### Subdirectories

| Directory | Contents |
|---|---|
| `patterns/` | the pattern pipeline and the verdicts on K-symmetric and rigid shapes: see [The pattern pipeline](#the-pattern-pipeline) |
| `estimates/` | `fibre_model.py`: the fibres of the exact level 2 at H_K = 1 (fibre sizes after the single-block filters, cost of the best two-pair join); `est_slots.py`, `slot_cost.py`, `n1_est.py`, `inst_mu.py`: per-slot Gaussian estimates, level-1 yield and level-2 cost per slot (cost files), measured level-1 counts, the instance-level completion estimator; `hit_analysis.py`, `symmetry_images.py`, `trans_test.py`, `summarize_runs.py`: the round that found skew order 1892 with all completions, the 128 images of its instance, translation invariance of level 2, run summaries |
| `scan38/` | see [scan38/README.md](scan38/README.md) |

## Data

| File | Contents |
|---|---|
| `data/cost473.jsonl` | cost file of the shape of 473 (K = ℤ₄₃, q = 11, H_K = ⟨6⟩): measured level-1 yield and level-2 time per slot; the slot file of `prod_ks.py` |
| `data/est473.jsonl` | the 58 slots of the 473 shape with a finite Gaussian estimate (slot file of `prod_dedup.py`) |
| `data/ctl473d_rounds.jsonl` | the logged first 20 rounds of the `prod_dedup.py` run that found the sixth matrix of order 1892 |
| `data/repro385/` | the slot of the run that found [skew/1540](../../skew/1540) |
| `data/repro245/`, `data/repro781/` | the slots of the lowest-energy sweep that found [skew/980](../../skew/980) and [skew/3124](../../skew/3124), and in `repro245/` the certificates of the other matrices of order 980 found in the sweep (`certificate_b.json` to `certificate_j.json`) |
| `data/repro539/` | the state before the hit of the run that found [skew/2156](../../skew/2156): the log of the first seven rounds, the set of classes already run, the slot file |

## Build

In this directory:

```
gcc -O3 -march=haswell -shared -fPIC -o libgl1.so gl1.c -lm
for ps in 32 64 128; do gcc -O3 -march=haswell -DPS=$ps -o sl2_$ps gl2.c -lm; done
for ps in 32 64; do gcc -O3 -march=haswell -DPS=$ps -o skf2_$ps pl2f2.c -lm; done
mkdir -p runs tmp
```

Python 3 with NumPy; `patterns/ksym_shapes.py` also needs SymPy. The commands below that use `est/` or `out/` create
them (`mkdir -p est out`).

## Reproducing the skew certificates

All commands run in this directory after the build. A group-form certificate `out/skew_<order>_<name>.json` is over
K × F_q. `cyclic.py IN.json OUT.json` turns one over ℤ_m × ℤ_q with gcd(m, q) = 1 into the cyclic form over ℤ_mq (the
format of [skew/820](../../skew/820)). Under each order below, the result is related to the certificate in the
repository. The times are wall-clock seconds on one core; those of order 1892 were measured on a shared four-core
machine (load 6 to 10), so on a free core they are shorter.

| Order | Program | Time | Reproducible |
|---|---|---|---|
| 476 | `sk_worker.py` | about 20 s | yes: seeded |
| 1892 | `prod_ks.py`, `prod_dedup.py` | 584 s (first matrix), 61 s (sixth matrix, replay) | yes: seeded, slot weights from a fixed cost file; the sixth matrix by replay of a logged state |
| 2156 | `prod_dedup.py` | about 11 s | by replay of a logged state; a fresh run is not bit-reproducible |
| 1540 | `exh38.py` | about 14 s | yes: the sweep of the hit slot is deterministic |
| 980 | `exh38.py` | 1 to 3 minutes | yes |
| 3124 | `exh38.py` | 1 to 3 minutes | yes |

### Order 476

K = ℤ₁₇, q = 7, H_K = 1 ([../NOTES.md](../NOTES.md) §87). One core, about 20 s; the rounds are seeded, so the run is
deterministic.

```
python3 -u sk_worker.py w119a 1 900 "[17]" 7 "[[2,2,2,3]]" 2000 20000 20000 30 1.0
python3 cyclic.py out/skew_476_w119a.json out/skew_476_cyclic.json
```

The arguments are the worker name, the seed base, a time limit in seconds, K, q, the defect allocations (or `E9,10,11`
for every allocation of those total weights), the list sizes of B, C, D, the number of level-1 instances passed to
level 2 per round, and the fraction of them kept. The first command writes `out/skew_476_w119a.json` (group form over
ℤ₁₇ × ℤ₇); the second gives the certificate of [skew/476](../../skew/476) without its `sha256` field.

### Order 1892

K = ℤ₄₃ with the multiplier 6 of order 3, q = 11 ([../NOTES.md](../NOTES.md) §88, §89).

```
python3 -u prod_ks.py k473a 1 9000 "[43]" 11 "[6]" data/cost473.jsonl 40 10000 60000 60000 sl2
python3 cyclic.py out/skew_1892_k473a.json out/skew_1892_cyclic.json
```

The arguments of `prod_ks.py` are the worker name, the seed base, a time limit in seconds, K, q, the multipliers on K,
the slot file, the number of slots kept, the list sizes of B, C, D and the level-2 engine (`sl2` or `skf`). A slot file
whose name contains `cost` is read as a cost file: only its records with measured times count, weighted by estimated
completions per second. Files `engine.txt` and `lists.txt` in the working directory, if present, override the engine
and the list sizes. `PROD_CONTINUE=1` makes the loop go on after a hit (writing `runs/NAME.HIT_*` and
`out/skew_<order>_NAME_r<round>.json`). The command is deterministic: the tenth round hits, and the second command gives
the certificate of [skew/1892](../../skew/1892) without its `sha256` field. In the same class hits came about once per
2 · 10⁴ to 3 · 10⁴ level-2 runs.

Three workers found five matrices with different row sums or allocations (inequivalence not tested), in 6,058 + 11,480 +
63,828 level-2 runs. A fourth worker (`k473s3`, 51,830 level-2 runs) found nothing; together the four account for the
133,196 level-2 runs of §88. The group-form certificate of each hit carries the sha256 below; `cyclic.py` gives the
cyclic certificate, whose sha256 is that of the matrix computed with the builders of `verify/verify.py`.

| Worker | Command | Time | Hits |
|---|---|---|---|
| `k473a` | `python3 -u prod_ks.py k473a 1 9000 "[43]" 11 "[6]" data/cost473.jsonl 40 10000 60000 60000 sl2` | 584 s | tenth round (`rounds: 10, n1: 6390, l2: 6058, nok: 8`) |
| `k473b` | `python3 -u prod_ks.py k473b 2 1500 "[43]" 11 "[6]" data/cost473.jsonl 40 10000 60000 60000 skf` | 523 s | round index 16: allocation (6, 1, 3, 10), m-split (6, 1, 0, 10), pattern ((10, 1), (14, −3), (11, −4)), 8 completing instances |
| `k473s4` | `PROD_CONTINUE=1 python3 -u prod_ks.py k473s4 4 6000 "[43]" 11 "[6]" data/cost473.jsonl 40 6000 60000 60000 skf` | about 45 minutes | round indices 16, 20 and 74; after round index 101: n1 70144, level-2 runs 63828, completions 8 |

The run of `k473s4` goes to its time limit; `touch runs/k473s4.STOP` ends it. The logged level-1 counts of `k473b` in
the run that found its matrix are larger than here (other list sizes); its slots, patterns and hit are the same.

| Matrix | Worker | Group-form sha256 | Row sums | Allocation | Cyclic sha256 |
|---|---|---|---|---|---|
| 1 | `k473a` | `ae96a274…` | (1, −9, 17, −39) | (0, 4, 7, 10) | `8ab883d6…` (skew/1892) |
| 2 | `k473b` | `916b2b7a…` | (1, 21, −19, −33) | (6, 1, 3, 10) | `9f661af1…` |
| 3 | `k473s4` | `88532699…` | (1, −19, −3, −39) | (0, 3, 7, 7) | `ddc0e2fa…` |
| 4 | `k473s4` | `95c308c7…` | (1, 17, −9, −39) | (6, 1, 7, 7) | `d9fb03aa…` |
| 5 | `k473s4` | `641257b8…` | (1, 35, 15, 21) | (0, 3, 7, 7) | `6fb7c60d…` |

**The sixth matrix** (§89; `prod_dedup.py`, K = ℤ₄₃, q = 11, H_K = ⟨6⟩, slots from `data/est473.jsonl`, the best 40 by
the Gaussian estimate). The slot weights of a fresh run depend on measured times, so the logged state of the first 20
rounds of the run that found the matrix is replayed:

```
cp data/ctl473d_rounds.jsonl runs/ctl473d.jsonl
python3 -u prod_dedup.py ctl473d 4 9000 "[43]" 11 "[[6]]" data/est473.jsonl 40 10000 60000 60000 skf
python3 cyclic.py out/skew_1892_ctl473d.json out/skew_1892_ctl473d_cyclic.json
python3 indep_check.py out/skew_1892_ctl473d.json
```

The first command takes 61 s. It resumes after the 20 logged rounds and prints `HIT {order: 1892, ... verify_ok:
True}` and `n1 5246, l2 900, nok 1`. The cyclic certificate has row sums [1, −21, −9, 37], multipliers [345, 221] and
sha256 `1918413ad673c0ec5b3aaf24bf979f10d2f00a74e5a505c76b624ce1772d0206`. `indep_check.py` (12 s) prints
`{"order": 1892, "pm1": true, "HHT_nI": true, "skew_H_plus_HT_2I": true, ..., "rowsums": [1, -21, -9, 37]}`.

A fresh run (`python3 -u prod_dedup.py ctl473e 5 3000 "[43]" 11 "[[6]]" data/est473.jsonl 40 10000 60000 60000 skf`,
with `PROD_CONTINUE=1` to go on after a hit) is not bit-reproducible. Hits came about once per 4000 to 6000 distinct
classes sent to level 2; the run that found the matrix used 6,093 distinct classes and 36,834 raw level-1 instances,
in 244 s on a loaded core. `estimates/summarize_runs.py runs/ctl473d.jsonl` prints the counts of the log (raw level-1
instances, new classes, duplicate factor 6.0 for this run). The three shapes of 539 with a trivial fixed-point set,
which §89 reports as the best remaining real candidates, are run the same way (the second gave skew order 2156):
`prod_dedup.py NAME 1 3000 "[49]" 11 "[[18]]" EST 40 10000 60000 60000 skf` (`EST` builds `est/NAME.jsonl` on first
use), `"[7,7]" 11 "[[2,2]]"` and `"[[2,4]]"` with engine `sl2` (the fast level 2 is for cyclic K).

### Order 2156

K = ℤ₇ × ℤ₇ with the scalar 2 of order 3, q = 11 ([../NOTES.md](../NOTES.md) §90). The run used `prod_dedup.py` with
slot weights from `est_extend.py` over all 2265 slots of the energy window. Its state before the hit (the log of the
first seven rounds of the job with seed base 5, the set of classes already run, the slot file) is in `data/repro539`.
The replay takes about 11 s:

```
mkdir -p runs est out tmp
cp data/repro539/c539_q11_K7x7_h2-2_s5.jsonl data/repro539/c539_q11_K7x7_h2-2_s5.seen runs/
cp data/repro539/est_c539_q11_K7x7_h2-2_s5.jsonl est/c539_q11_K7x7_h2-2_s5.jsonl
PROD_EREL=5 PROD_NEST=60 PROD_MAXOUT=30000 PROD_L2CAP=3000 python3 -u prod_dedup.py c539_q11_K7x7_h2-2_s5 5 420 "[7,7]" 11 "[[2,2]]" EST 40 10000 60000 60000 sl2
```

The loop resumes at the eighth round (index 7) and writes `out/skew_2156_c539_q11_K7x7_h2-2_s5.json`, the certificate of
[skew/2156](../../skew/2156) without the fields `structure` and `subgroup_order`. The slot file itself comes from
`python3 est_extend.py NAME "[7,7]" 11 "[[2,2]]" 60 SECS 5` and `python3 est_merge.py NAME` (in `est/`). A fresh run of
`prod_dedup.py` in this mode is not bit-reproducible, since the slot weights use measured times.

### Order 1540

K = ℤ₃₅ with ⟨11⟩, q = 11 ([../NOTES.md](../NOTES.md) §91). The run swept the slots of the shape with `exh38.py` in
the order of the estimated completions per class. Its seventh slot, allocation (6, 1, 1, 4) with m-split (0, 0, 0, 1),
holds the matrix. The arguments of `exh38.py` are the job name, a time limit in seconds, K, q, H_K, the slot file, the
number of slots, their order (`rho`, `E`, `mix`, or `file` for the order of the slot file) and the level-2 engine. The
sweep of the hit slot alone is deterministic and takes about 14 s:

```
mkdir -p runs est out tmp
cp data/repro385/est_hit_slot.jsonl .
PROD_CONTINUE=1 python3 -u exh38.py hit385 600 "[35]" 11 "[[11]]" est_hit_slot.jsonl 1 file skf
python3 cyclic.py out/skew_1540_hit385_s0_<time>_<i>.json OUT.json      # the first HIT line
```

It prints two `HIT` lines (group-form SHA-256 `9ea5b315…` and `1b4d57ec…`, 100 classes of level-1 instances, 2
complete); `cyclic.py` turns the first into the certificate of [skew/1540](../../skew/1540) without its `sha256`
field.

### Orders 980 and 3124

The lowest-energy sweep of §93 found order 980 (K = ℤ₃₅, q = 7, ⟨11⟩; group ℤ₃₅ × ℤ₇, not cyclic) and order 3124
(K = ℤ₇₁, q = 11, ⟨25⟩; group ℤ₇₈₁). Each hit slot replays alone, in 1 to 3 minutes; the list caps are part of the
command:

```
cp data/repro245/est_hit245.jsonl data/repro781/est_hit781.jsonl .
EX_COV=0.99 EX_PERSIST_MB=0 PROD_CONTINUE=1 python3 -u exh38.py hit245 1200 "[35]" 7 "[[11]]" est_hit245.jsonl 1 file skf
EX_COV=0.99 EX_PERSIST_MB=0 EX_MAXDRAW=6000000 EX_JOINCAP=2e9 EX_SLOT_SECS=600 PROD_CONTINUE=1 python3 -u exh38.py hit781 1500 "[71]" 11 "[[25]]" est_hit781.jsonl 1 file sl2
python3 cyclic.py out/skew_3124_hit781_s0_<time>_<i>.json OUT.json
```

The first `exh38.py` command prints three `HIT` lines (SHA-256 `74f3a8f1…`, the certificate of
[skew/980](../../skew/980), `53465d13…` and `ed5594ef…`; 11,812 classes); the other two are `certificate_b.json` and
`certificate_c.json` in `data/repro245`. The second prints one (group form `0e0e2e2f…`, 1841 classes), which `cyclic.py`
turns into the certificate of [skew/3124](../../skew/3124) without its `sha256` field.

## The pattern pipeline

For shapes with many defects the exhaustive allocation pipeline `sk_exh.py` stops at lists of about 10⁷ blocks. The
pattern pipeline in `patterns/` first cuts a shape into classes of *defect-position patterns* (for each block the
H_K-orbits that carry a defect and the type p (s = 0) or m (s = −2e) of each; the skew block A uses pairs {O, −O}) and
decides each class by an explicit-list level 1 and level 2. A shape is **EXACTLY empty** when steps 1 and 2 already
leave no pattern, and **EXHAUSTIVELY empty** when every class has no level-2 completion (§88 and §89 use the two words
in this sense).

1. *Big–big filter* (`bbpat.py`; an exact necessary condition). With u = s + q e, |u| = 1 off the defects and q or q − 2
   on them. At a lag k ≠ 0 the products of two big entries have absolute value at least (q − 2)², and everything else
   in Σ_X PAF_{u_X}(k) is bounded by
   N(k) = Σ_X [ Σ_{y ∈ D_X} |u_y| ([y + k ∉ D_X] + [y − k ∉ D_X]) + #{x : x, x + k ∉ D_X} ].
   So |Σ_X BB_X(k)| ≤ N(k) for some choice of the orbit signs. In the rigid regime (q − 2)² > 2Eq + 4|K| a lag with
   exactly one big–big pair is impossible (count test).
2. *Row-sum filter* (`bbpat.rowfilter`): the exact row-sum patterns (`sk_pat.py`) of the (allocation, m-split) class of
   each surviving pattern.
3. *Symmetry reduction* (`patsym.py`): patterns modulo translation of B, C, D by points of K fixed by H_K, reflection of
   each block, permutation of (B, C, D) and the units of K that normalise H_K.
4. *Exhaustive class* (`pat_exh.py`): the lists of all sign configurations of the pattern, the exact join of `gl1.c`,
   the exact completion of `gl2.c`.

A shape is given as `'K' q 'hunits'` (see [Method and terms](#method-and-terms)); H_K contains no twisted generators
(`patterns/twist_check.py`). All commands run in this directory.

| File | Contents |
|---|---|
| `patterns/bbpat.py`, `patsym.py`, `pat_exh.py` | the pipeline modules (steps 1–4); `bbpat.py 'K' q 'hunits' [tlimit]` lists the alive patterns |
| `patterns/run_pat.py`, `run_bb.py` | verdict by the whole pipeline (resumable, one JSON line per class); verdict by steps 1–3 only (`run_bb.py` prints one RESULT line, `after row filter 0` = EXACTLY empty) |
| `patterns/exh_shape.py` | verdict by the allocation pipeline (explicit lists per allocation, `sk_exh.py`), with limits `MAXINST`, `PAIRS_HI` |
| `patterns/bbhk1.py` | the pattern search for H_K = 1 (defect sets up to translation and reflection, count test, signed test, row-sum filter) |
| `patterns/check_weights.py` | the orbit weights cannot realise the defect identity (an independent dynamic program) |
| `patterns/twist_check.py` | twisted generators: no group containing (h, ν), ν a non-square, admits an invariant skew block (G-level check) |
| `patterns/nw2_lemma.py`, `nw_scan.py`, `nw2_indep.py`, `lemma_check.py` | the N_w = 2 lemma: row-sum patterns that satisfy it, exact class-restricted level 1, an independent numpy join, the lemma on a repository solution |
| `patterns/rowscan.py`, `ksym_shapes.py` | exact row-sum feasibility of all slots of a shape; the list of K-symmetric census shapes (needs SymPy) |
| `patterns/test_bb_sound.py`, `test_rowfilter.py`, `val_killed.py` | soundness of the big–big filter and of the row-sum filter on actual solutions, and on killed patterns |
| `patterns/val_pat.py`, `val_n1.py`, `val_sym.py`, `pat_indep.py` | the pattern pipeline against the allocation pipeline, against an independent join, and the symmetry reduction |
| `patterns/val_new.py`, `val_q27.py` | G-level against K-level brute force against the pipeline for matrix-valued H_K, non-cyclic K and q = 27; they append to `runs/val_new.jsonl`, `runs/val_q27.jsonl` |
| `patterns/posctl_1100.py` | positive control: [skew/1100](../../skew/1100) (K = ℤ₅², H_K = Singer subgroup of order 3) is recovered from its pattern |

`val_skf.sh` also needs the `skf2_*` programs of the build.

## Reproduction of search/NOTES.md §88–89

Commands run in this directory after the build; the `patterns/` and `estimates/` scripts are called by path. The result
column is what the program prints, shortened. Times are wall-clock seconds on one core of a shared four-core machine
(load 6 to 10), so on a free core they are shorter. A shape is `'K' q 'hunits'`. The matrices of order 1892 of §88 and
§89 are reproduced under [Order 1892](#order-1892).

### §88: twisted generators are impossible

G-level check; the proof is in NOTES.

| Command | Time | Result |
|---|---|---|
| `python3 patterns/twist_check.py small` | 6 s | `TOTAL square-nu generators (tested, compatible): [756, 234]  non-square-nu generators (tested, compatible): [756, 0]  disagreements with "compatible iff odd order": 0` |
| `python3 patterns/twist_check.py` | 32 s | `TOTAL square-nu generators (tested, compatible): [2262, 551]  non-square-nu generators (tested, compatible): [2262, 0]  disagreements ... : 0` |

The first run covers the 1512 generators of seven small K: those that admit an invariant skew block are exactly those of
odd order. The second covers the 2262 generators with a non-square ν of nine small K, none of which admits one.

### §88: the fast level 2 (`skf2`) against `sl2`

`./val_skf.sh` runs the battery of 22 slots of seven shapes. Every RESULT line must end `count/status mismatches 0,
solution-set mismatches 0`; two slots print only NOINST (no level-1 instance in the sample). In the recorded run `skf2`
was 1.8 to 6.7 times faster. Two slots run here:

| Command | Time | Result |
|---|---|---|
| `python3 val_skf.py "[19]" 7 "[7]" "[0, 1, 3, 7]" "[0, 0, 3, 0]" 3000 10000 10000 400 267` | 2 s | 400 instances, 120 with completions, 288 completions, count/status mismatches 0, solution-set mismatches 0 (1.67x) |
| `python3 val_skf.py "[43]" 11 "[6]" "[6,3,4,7]" "[6,0,4,7]" 3000 30000 30000 1500 3` | 408 s | 1500 instances, 0 completions, mismatches 0 and 0, time `sl2` 296 s, `skf` 103 s (2.88x) |

### §88: anatomy of the hit of order 1892 and the estimators (HEURISTIC)

| Command | Time | Result |
|---|---|---|
| `python3 estimates/hit_analysis.py` | 141 s | level-1 instances 4032 ... instances with completions 8, completions (x2 gauge) 16 `{3634: 2, ..., 3641: 2}`, and the block hashes: A takes 2 values, B 2, C 2, D 1 over the 8 instances: one solution cluster (2 × 2 × 2) |
| `python3 estimates/n1_est.py "[43]" 11 "[[6]]" "[0,4,7,10]" "[0,4,6,10]" 10000 60000 60000 4 1` | 77 s | `mean 9.044e+06, rel.sd 0.29 ... log10 = 6.96`; the Gaussian estimate of the same slot is logN1 = 7.12 (1.3e7, in `data/est473.jsonl`): a factor 1.4 |
| `python3 estimates/inst_mu.py` | 201 s | `mean mu over non-completing sample 6.89e-06; expected completions in all 4032 instances ~ 0.0278 (observed total completions 16)` |

### §88: fibres of the exact level 2 at H_K = 1

`fibre_model.py K q NAL` takes NAL random slots of seed 1; the explicit slot is `fibre_model.py "[35]" 7 "[2,0,6,22]"
"[2,0,6,17]"`. The result is the log₂ of the best two-pair join cost. Each run takes one to two seconds. They give the
numbers of the table behind "2⁴⁴ (ℤ₂₅) to 2⁷⁷ (ℤ₄₃) pair evaluations per level-1 instance".

| v | Command | Result |
|---|---|---|
| 575 | `python3 estimates/fibre_model.py "[25]" 23 "[0,0,2,3]" "[0,0,2,3]"` | n_t per block [13, 24, 23, 22], log₂ fibre sizes [11.4, 23.4, 22.6, 21.6], join cost 44.2 |
| 667 | `python3 estimates/fibre_model.py "[29]" 23 1` | join cost 52.0 |
| | `"[31]" 19 1`, `"[33]" 19 1`, `"[33]" 23 1` | 54.0, 57.8, 58.8 |
| 805 | `python3 estimates/fibre_model.py "[35]" 23 1` | 62.8 |
| | `"[41]" 23 1`, `"[43]" 23 1` (989) | 72.5, 77.4 |

### §88: the big–big filter is sound, the pattern pipeline equals the allocation pipeline, skew/1100 is recovered

| Command | Time | Result |
|---|---|---|
| `python3 patterns/posctl_1100.py` | < 1 s | K-level check ... True, pattern alive (True, True), `pipeline on the pattern: {'n1': 64, ..., 'compl': 40, ...}`, 3 × `recovered solution K-level check: True` |
| `python3 patterns/test_bb_sound.py 0 1 2` | 99 s | E_G 176256 / 110592 / 74880000 (the exact counts of §87), each with 400 solutions checked: `pattern NOT alive: 0 direct bound holds: 400`; without arguments all 13 shapes |
| `python3 patterns/test_rowfilter.py 0 1 2` | 104 s | 300 solutions per shape: `classes rejected by the row filter: 0`; without arguments all nine shapes |
| `python3 patterns/val_pat.py "[7]" 27 "[[2]]" 6144` | 0.5 s | `patterns 1`, `ordered completions x 8 = 6144 expected 6144 MATCH` (ℤ₇ × GF(27)) |
| `python3 patterns/val_pat.py "[3,3,3]" 7 "[[[0,0,1],[0,1,1],[2,1,1]]]" 6144` | 1 s | `patterns 14` ... MATCH 6144 (ℤ₃³ with a Singer matrix of order 13 × F₇) |
| `python3 patterns/val_pat.py "[13]" 7 "[[3]]" 930816` | 17 s | `patterns 1156` ... MATCH 930816 (ℤ₁₃ with ⟨3⟩ × F₇) |
| `python3 patterns/val_n1.py "[33]" 23 "[[4]]" 124416` | 7 s | `alive patterns 306 killed(signed) 180 killed(count) 0`, `ordered level-1 instances x 8 = 124416 expected N1_G 124416 MATCH` (759 over ℤ₃₃: the filter removes 180 patterns and no level-1 instance) |
| `python3 patterns/val_sym.py "[3]" 7 "[]"` | 1 s | `patterns 38 canonical classes 6 classes with nonzero completions 2 classes with INCONSISTENT zero-ness 0`; for `"[3]" 11 "[]"` the counts are 3, 1, 1, 0; also `"[5]" 7 "[]"`, `"[13]" 7 "[[3]]"` |
| `python3 patterns/val_killed.py "[39]" 23 "[[16]]" 12` | 188 s | 897; 12 of the 7,332 killed patterns: `killed patterns with level-1 instances (must be 0): 0 total instances 0` (the recorded check used 120 patterns: about 30 minutes here) |
| `python3 patterns/val_new.py 6` | | ℤ₃³ × F₇: gcount 6144, kcount 6144, pipe 6144 MATCH; index 5 (q = 3): 0 = 0 = 0 |
| `python3 patterns/val_q27.py 0` | | ℤ₇ ⟨2⟩ × GF(27): 6144 = 6144 = 6144 |

### §88: shapes EXHAUSTIVELY empty

The last line must say `EXHAUSTIVELY empty`, or `exhaustive: true`, `E_G: 0`, `nskipped: 0` for `exh_shape.py`.

| v | Command | Time | Result |
|---|---|---|---|
| 589 | `python3 patterns/run_pat.py runs/p.jsonl "[31]" 19 "[[5]]" 600` | 1.5 s | classes 7, done 7, level-1 instances 0, completions 0: EXHAUSTIVELY empty |
| 589 | `python3 patterns/exh_shape.py runs/e.jsonl "[31]" 19 "[[5]]" 900` | 288 s | by the allocation pipeline: `"allocs": 5, "exhaustive": true, "E_G": 0, "N1_G": 0` |
| 979 | `python3 patterns/exh_shape.py runs/e.jsonl "[89]" 11 "[[2]]" 900` | 90 s | `"allocs": 31, "exhaustive": true, "E_G": 0, "N1_G": 0, "nskipped": 0` |
| 759 (ℤ₆₉, q = 11) | `python3 patterns/exh_shape.py runs/e.jsonl "[69]" 11 "[[4]]" 1500` | 367 s | `"allocs": 209, "exhaustive": true, "E_G": 0, "N1_G": 5557248` |
| 989 | `python3 patterns/run_pat.py runs/p.jsonl "[43]" 23 "[[6]]" 900` | 10 s | classes 1, done 1, level-1 instances 0, completions 0: EXHAUSTIVELY empty |
| 1209 | `python3 patterns/run_pat.py runs/p.jsonl "[39]" 31 "[[16]]" 900` | 1400 s (loaded; 17 s for the last class after a restart of the same OUT file) | classes 6, level-1 instances 1696 (544 + 1152), completions 0: EXHAUSTIVELY empty |
| 805 | `python3 patterns/run_pat.py runs/p.jsonl "[35]" 23 "[[11]]" 900` | 564 s | 8 classes, level-1 instances 2384, completions 0: EXHAUSTIVELY empty |
| 575 | `python3 patterns/run_pat.py runs/p.jsonl "[25]" 23 "[[21]]" 900` | 2 s | classes 3, level-1 instances 0 |
| 1953 (ℤ₃² × ℤ₇, order 9) | `python3 patterns/run_pat.py runs/p.jsonl "[3,3,7]" 31 "[[[1,1,0],[0,1,0],[0,0,1]],[[1,0,0],[0,1,0],[0,0,2]]]" 900` | 307 s | classes 9, level-1 instances 0: EXHAUSTIVELY empty |
| 627 (ℤ₅₇, q = 11) | `python3 patterns/exh_shape.py runs/e.jsonl "[57]" 11 "[[4]]" 1500` | 689 s | `"allocs": 209, "exhaustive": true, "E_G": 0, "N1_G": 11985408, "nskipped": 0` |
| 1539 (ℤ₈₁, q = 19) | `python3 patterns/exh_shape.py runs/e.jsonl "[81]" 19 "[[4]]" 1500` | 411 s | `"allocs": 427, "exhaustive": true, "E_G": 0, "N1_G": 0, "nskipped": 0` |
| 817 | `python3 patterns/run_pat.py runs/p.jsonl "[43]" 19 "[[6]]" 7200` | 3328 s (about 2000 s on a free core) | `alive patterns 9380, mod symmetry 317` (`run_bb.py`: alive 12726 after the big–big filter = 457 classes mod symmetry, 9380 after the row-sum filter = 317); RESULT: classes 317, done 317, level-1 instances 0: EXHAUSTIVELY empty |
| 623 | `MAXINST=2e7 PAIRS_HI=3e10 python3 patterns/exh_shape.py runs/e.jsonl "[89]" 7 "[[2]]" 9000` | about 1.5 h; not re-run here (UNVERIFIED) | 374 allocations, 1.6 × 10⁸ level-1 instances, none completes |
| 539 | `python3 patterns/val_new.py 25` | 163 s | ℤ₇² × F₁₁ with ⟨2 I, unipotent⟩ (order 21): `'gcount': 0, 'kcount': 0, 'pipe': 0, 'pipe_exhaustive': True, 'status': 'MATCH'` (G = 539) |
| 1029, 2401 | `python3 patterns/val_new.py 34,35` | 46 s | ℤ₇³ × F₃ and × F₇, Singer subgroup of order 171: `'gcount': 0, 'kcount': 0, 'pipe': 0 ... 'status': 'MATCH'` for G = 1029 and G = 2401 |
| 1323 | `python3 patterns/val_q27.py 9` | 856 s (loaded; about 340 s unloaded) | ℤ₇² × GF(27), order 21: `'gcount': 0, 'kcount': 0, 'pipe': 0, 'exh': True, 'status': 'MATCH'` (G = 1323) |

### §88: shapes EXACTLY empty

No join is run.

| v | Command | Time | Result |
|---|---|---|---|
| 851 | `python3 patterns/run_bb.py "[37]" 23 "[[10]]"` | 0.3 s | `tested 0`: the orbit weights (1 and 3) cannot realise (e_p, e_m) = (2, 5); `patterns/check_weights.py "[37]" 23 "[[10]]"` says `weights realise a target: False` |
| 897 | `python3 patterns/run_bb.py "[39]" 23 "[[16]]"` | 9 s | `tested 7332 killed(count) 0 killed(signed) 7332 alive 0` |
| 1333 | `python3 patterns/run_bb.py "[43]" 31 "[[6]]"` | 0.6 s | `tested 28 killed(count) 28 killed(signed) 0 alive 0` (every pattern has a lone big–big pair) |
| 1323, H_K of order 7 | `python3 patterns/run_bb.py "[7,7]" 27 "[[[1,1],[0,1]]]"` | 438 s | `tested 483665 ... alive 29596 after row filter 0` |
| 1323, H_K of order 3 | `python3 patterns/run_bb.py "[7,7]" 27 "[[2,2]]"` | 2 s | `alive 256 after row filter 0` |
| 1859, 2375 | `python3 patterns/check_weights.py` | 22 s | 1859 (Singer subgroup of order 21) and 2375 (Singer subgroup of order 31): `weights realise a target: False \| bbpat tests 0 patterns: agrees`; the other eight shapes: True, agrees |
| 361 | `python3 patterns/bbhk1.py "[19]" 19` | 0.7 s | `ALIVE 0` (H_K = 1; likewise `"[23]" 23` for 529) |
| 899 | `python3 patterns/bbhk1.py "[29]" 31` | 1.6 s | `ALIVE 0` |
| 961 | `python3 patterns/bbhk1.py "[31]" 31` | 1.7 s | `ALIVE 0` |
| 1763 | `python3 patterns/bbhk1.py "[41]" 43` | 2.5 s | `ALIVE 0` |
| 1269 (K = ℤ₄₇, q = 27) | `python3 patterns/bbhk1.py "[47]" 27` | 1081 s | `leaves 2276264 count-alive 713 typed tested 713 killed(signed) 713 ALIVE 0 timeout False` |

### §88: the N_w = 2 lemma

1081 = 23 · 47 and 1711 = 29 · 59 have no solution in this product for any H_K.

| Command | Time | Result |
|---|---|---|
| `python3 patterns/nw2_lemma.py` | 1 s | m = 11: 7 row-sum patterns, 2 satisfy the lemma; m = 15: 8, 2; m = 23: 4, 0; m = 29: 6, 0; m = 35: 17, 3 |
| `python3 patterns/nw_scan.py 11 23` | 0.7 s | `SUMMARY 11 23 patterns 7 total level-1 instances ... 52360` (35,640 + 16,720 in the two classes that satisfy the lemma) |
| `python3 patterns/nw_scan.py 15 31` | 3.5 s | `SUMMARY 15 31 patterns 8 ... 1571520` |
| `python3 patterns/nw_scan.py 23 47` | 117 s | 1081: `SUMMARY 23 47 patterns 4 total level-1 instances (canonical signs, classes run) 0 skipped classes 0` |
| `python3 patterns/nw2_indep.py 2 11 23` | | 35640 (independent numpy join; `3 11 23`: 16720) |
| `python3 patterns/nw2_indep.py 3` | 1189 s on a free core | 1081, class 3 (lists 2048, 1771, 26334, 170544): `sizes 2048 1771 26334 170544`, `pattern 3: level-1 key matches 0 (6405s)`; the only class of 1081 whose A + B table fits in memory |
| `python3 patterns/lemma_check.py` | | [skew/3612](../../skew/3612) (903 = 21 · 43, N_w = 2): block 1 `m` defect with rest of s = 2, block 2 `p` defect with rest −2: the lemma holds on a real solution |

### §89: census shapes of the 473 type

K-symmetric, real; the `run_pat`, `exh_shape` and `rowscan` commands of §88 with the shapes of this section.

| v | Command | Time | Result |
|---|---|---|---|
| census | `python3 patterns/ksym_shapes.py runs/shapes.json 21 200` | 363 s | `census values 737 undecided now 643`, `shapes 90` (the 90 shapes of the census, listed with K, q, H_K, nu, E) |
| 627 (ℤ₃₃ × F₁₉) | `python3 patterns/run_pat.py runs/p.jsonl "[33]" 19 "[[4]]" 900` | 4 s | classes 5, level-1 instances 0: EXHAUSTIVELY empty |
| 759 (ℤ₃₃ × F₂₃) | `python3 patterns/run_pat.py runs/p.jsonl "[33]" 23 "[[4]]" 900` | 4 s | classes 9, level-1 instances 64, completions 0: EXHAUSTIVELY empty |
| 1539 (ℤ₅₇ × GF(27)) | `python3 patterns/run_pat.py runs/p.jsonl "[57]" 27 "[[4]]" 900` | 3 s | classes 5, level-1 instances 16 |
| 1953 (ℤ₆₃ × F₃₁, ⟨4, 22⟩) | `python3 patterns/run_pat.py runs/p.jsonl "[63]" 31 "[[4],[22]]" 900` | 192 s | classes 9, level-1 instances 0: EXHAUSTIVELY empty |
| 2001 (ℤ₈₇ × F₂₃) | `python3 patterns/run_pat.py runs/p.jsonl "[87]" 23 "[[7]]" 3000` | 1581 s | classes 114, level-1 instances 0: EXHAUSTIVELY empty |
| 2349 (ℤ₈₇ × GF(27)) | `python3 patterns/run_pat.py runs/p.jsonl "[87]" 27 "[[7]]" 3000` | 951 s | classes 52, level-1 instances 0: EXHAUSTIVELY empty |
| 781 (ℤ₇₁ × F₁₁, ⟨20⟩) | `python3 patterns/run_pat.py runs/p.jsonl "[71]" 11 "[[20]]" 3000` | 179 s | classes 2047, level-1 instances 426, completions 0: EXHAUSTIVELY empty |
| 803 (ℤ₇₃ × F₁₁, ⟨2⟩) | `python3 patterns/exh_shape.py runs/e.jsonl "[73]" 11 "[[2]]" 900` | 103 s | `"allocs": 31, "exhaustive": true, "E_G": 0, "N1_G": 0, "nskipped": 0` |
| 1083 | `python3 patterns/run_bb.py "[57]" 19 "[[4]]"` | 4 s | `tested 3054 killed(count) 0 killed(signed) 3054 alive 0` (every pattern fails the big–big bound) |
| 1705 | `python3 patterns/run_bb.py "[55]" 31 "[[16]]"` | 15 s | `tested 14310 killed(count) 4600 killed(signed) 9710 alive 0` |
| 1323 (ℤ₇² × GF(27), diagonal group of order 9) | `python3 patterns/run_bb.py "[7,7]" 27 "[[1,2],[2,1]]"` | 0.1 s | `alive 16 after row filter 0`; `"[[2,2]]"` (order 3): `alive 256 after row filter 0`; `"[[2,4]]"`: `alive 16 after row filter 0` |
| 1199 | `python3 patterns/rowscan.py "[109]" 11 "[[16]]" 100000` | 84 s | `pool 3056 slots, tested 3056, with row-sum patterns 0` |
| 1397 | `python3 patterns/rowscan.py "[127]" 11 "[[22]]" 100000 99 3000` | 338 s | `pool 7870 slots, tested 7870, with row-sum patterns 0` |
| 2071 | `python3 patterns/rowscan.py "[109]" 19 "[[16]]" 100000` | 0.5 s | `pool 5 slots, tested 5, with row-sum patterns 0`; `"[73]" 11 "[[2]]"` (803): `pool 317 slots ... 0` |

### §89: symmetric copies of the 1892 instance

| Command | Time | Result |
|---|---|---|
| `python3 estimates/symmetry_images.py` | 20 s | `images 128 distinct 128`; `images with level-2 completions: skf 128 / 128, sl2 128 / 128; completion counts equal in both engines: True` |
| `python3 estimates/trans_test.py` | 2 s | `fix [0, 7, 14, 21, 28]`, `base instances 18 with differing verdict/fibres among translates: 0` |

### Validation battery

`python3 validate_skew.py` re-runs the validation battery (resumable, results in `val.jsonl`): 31 shapes, on which
G-level and K-level brute force agree. The two-level pipeline agrees with them on all but those whose counts are too
large for it (ℤ₅ × F₃, ℤ₁₃/⟨3⟩ × F₃) or that hit its time limit.

## Status labels

The labels follow [../NOTES.md](../NOTES.md).

| Label | Meaning |
|---|---|
| EXACTLY empty | steps 1 and 2 of the pattern pipeline leave no pattern for the shape, so no join is run |
| EXHAUSTIVELY empty | every class of the shape has no level-2 completion |
| HEURISTIC | an estimate or a sampled run: it ranks slots or shapes and cannot discard one |
| UNVERIFIED | the recorded result was not re-run |
