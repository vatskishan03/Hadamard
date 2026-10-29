# Origin-defect shapes of the quaternary Paley plug-in: the K-level analysis

The programs behind the theory part of [../../NOTES.md](../../NOTES.md) §90. The plug-in of §89 (G = K × F_q, q ≡ 3 mod
4, H = H_K × squares, quaternary blocks U = s ⊗ δ + t ⊗ Q + e ⊗ J, A twisted-skew) has *origin-defect* shapes: every
position off k = 0 is of type E (x₂ = −x₁, e = 0, s and t unit-valued) and the defects ε_X δ₀ of the blocks sit at the
origin of K. For these the K-level system reduces to a problem on ℤ_m with unit-valued sequences (the identities (d),
(Q), (J) of §90, the "Lock" u_B = λ u_A for the rigid pair, the row-sum identity).

The programs analyse the solution sets of `../qp.c` on such shapes, scan the family `d1=1` over K = ℤ_m, q, H_K, and
enumerate exactly the rigid shape of 517 = 11 · 47.

## Contents

| File | Contents |
|---|---|
| `qp_oo.patch` | the two-line change of `../qp.c` that adds `d1=2` (defects only at k = 0 also for singleton H_K-orbits, i.e. every non-origin orbit of type E; `d1=1` of `qp.c` restricts only the orbits with more than one point): `patch -o qp_oo.c ../qp.c qp_oo.patch` |
| `shapes_lib.py`, `oo_shapes.py` | odd-order subgroups of ℤ_m\* and their orbits; the list of shapes (m, q, H_K) whose defect costs can add up to 4m (A's half defect (q ± 1)/2 or (q − 3)/2, the other blocks' origin costs 0, (q ± 1)/2, (q − 3)/2, q ± 1, q − 3): `python3 oo_shapes.py MMIN MMAX NUMAX` |
| `run_oo.py`, `run_list.py` | `qp d1=1` over that list (`NUMIN=… NUMAX=… python3 run_oo.py MMIN MMAX NUMIN NUMAX TMO OUT.jsonl first=1 maxprint=1`) or over explicit pairs (`NUMIN=… NUMAX=… QPBIN=qp_oo python3 run_list.py TMO OUT.jsonl "m:q,m:q" d1=2 first=1`); resumable jsonl logs; `make_table.py` |
| `common.py` | loading of SOL files, the shape (orbits, sizes), PAF |
| `claims_all.py`, `verify_identities.py` | an exact pass over ALL solutions of a file: (J_k), (d_k), (Q) at every lag, no defect off the origin, Z₁(k) = Σ_X Φ_X(k) constant in k (value −z, z ∈ {0, 1}), the energies (c_s, c_t), the defect patterns; the same identities from the definitions, on every n-th solution |
| `rowid.py`, `rowsums.py`, `z1vals.py`, `anpat.py`, `zchar.py` | the row-sum identity Σ(\|R^s\|² − \|S\|²), Σ(\|R^t\|² − \|τ\|²) = (4(m − 1) + qZ, 4(m − 1) − Z); the row sums by block type; Z₁; the patterns (type and cost of the origin defects of A, B, C, D); the Legendre-affine blocks (z = 1) |
| `sym299.py`, `sym589.py`, `tA_rel.py`, `rigid_rel.py`, `overlap13.py`, `overlap5_13.py` | the symmetry group of the 299 solutions (13,056 = 136 classes of 96) and of the 589 sample; t_A against a multiplier image of s_A; block relations in the rigid solutions; the class values shared between (5, 11), (13, 19), (13, 23), (13, 27) |
| `lockcount.py`, `sstage.py`, `sstage2.py`, `sstage3.py`, `sstage_scan.py`, `sstage_rigid.py` | the local states of the Lock equation per pair ±O (library `states`), the (Lock + S) stage of the K-level problem for a defect pattern (counts of s-quadruples, against the engine's level-1 instances at 299), and its counts over shapes and for the rigid pattern |
| `rowfeas.py`, `rowfeas_scan.py` | the row-sum (trivial character) feasibility of origin-defect shapes (EXACT necessary condition): kills every shape with \|H_K\| ≥ 11 and most with \|H_K\| = 9, 13, 15 at ν ≤ 7 |
| `heur.py` | the Gaussian model of the expected number of solution classes per origin variant (HEURISTIC) with the row-sum filter: `python3 heur.py MMAX NUMAX` |
| `unif_test.py`, `unif_test2.py`, `unif_test_all.py` | the uniform cyclotomic conjecture (class values independent of the prime m) tested by extending solutions to other primes m': falsified |
| `rigid_s.py` | the exact (Lock + S) count for the rigid pattern (A, M \| E, E), K = ℤ_m, H_K = 1: all skew A, all C, D modulo translation and unit, exact linear-key meet in the middle; `python3 rigid_s.py m q` |
| `rigid_s_check.py`, `rigid_s_check2.py` | an independent reimplementation (FFT, full D, no linear keys): `rigid_s_check.py m q NSAMPLE SEED` (NSAMPLE ≤ 0: all A); `rigid_s_check2.py` covers the other origin variants |
| `rigid_j.py` | a third implementation: Gaussian-integer meet in the middle with 64-bit linear hashes and exact verification of every key match; `rigid_j.py m q [maxA] [variant]` |
| `rigid_s_test.py` | the planted-solution test of the meet in the middle |
| `rig_mod.py`, `rig_lags.py`, `rig7_struct.py`, `rig7_paf.py` | diagnostics of the rigid shapes: the 2-adic structure of the equation (none), counts per A and per lag (why the Gaussian heuristic fails there), the structure of the m = 7 solutions |
| `cpsat_plug.py` | a CP-SAT model of the K-level system (accepts the known 299 solution, does not find one at 13 × 23: a negative result; needs OR-Tools) |
| `real2156.py`, `g_check27.py` | the anatomy of the real solution 2156 (`data/c539_hit.json`: defect support and costs, the level-1 and (d) identities) and a G-level check of the (13, 27) solutions over ℤ₁₃ × F₂₇ |
| `data/` | recorded outputs: the logs `runs_*.jsonl` of the family scans (`runs_nu3`: the 109 shapes with 3 orbits; `runs_oo5`, `runs_oo7`, `runs_oo7long`, `runs_nu9`, `runs_q27`, `runs_pp`, `runs_small`, `runs_m3`, `runs_m9`, `runs_m9_d1`), `results_table.txt`, `claims_log.txt`, `heur_out.txt`, `nu3_list_all.txt`, `pp_list.txt`, and `c539_hit.json` (the solution file of skew order 2156) |

## Build

In this directory, after `gcc -O3 -march=haswell -o ../qp ../qp.c -lm` (Python 3 with NumPy and SymPy):

```
patch -o qp_oo.c ../qp.c qp_oo.patch && gcc -O3 -march=haswell -o qp_oo qp_oo.c -lm
mkdir -p shapes runs out
```

The scripts read solution files made by `../qp` (SOL lines, one digit per K-orbit for each of x₀, x₁, x₂ of the four
blocks). `QPBIN` names the engine of `run_list.py` and `run_oo.py` (default `../qp`; `qp_oo` for `d1=2`). Certificates:
`python3 ../plugcert.py m q "[H_K]" "SOL …"` (order 8mq, `../out/cert_<8mq>.json` with its sha256).

## Reproduction

Run in this directory after the build. Times are wall times on one core of a machine under load.

**The solution sets** (EXACT counts, normalised: B, C, D with x₀(0) = 1):

```
python3 ../qshape.py '[13]' 23 '[3]' shapes/p299.txt;  ../qp shapes/p299.txt l2c=1 maxprint=20000 > runs/299_all.out      # 20 s   13,056 solutions, 36,000 level-1 instances (the matrix of skew/2392)
python3 ../qshape.py '[13]' 19 '[3]' shapes/o_13_19_3.txt; ../qp shapes/o_13_19_3.txt l2c=1 d1=1 maxprint=100000 > runs/13_19_all.out   # 12 s  COUNT 1536 (v = 247)
python3 ../qshape.py '[13]' 27 '[3]' shapes/o_13_27_3.txt; ../qp shapes/o_13_27_3.txt l2c=1 d1=1 maxprint=100000 > runs/13_27_all.out   # 29 s  COUNT 2304 (v = 351, q = 27)
python3 ../plugcert.py 13 19 "[3]" "02332|32310|20132 01110|30322|32100 00303|20123|32301 00131|33103|31321"      # 2 s   skew 1976 = 8 · 247, sha256 8bddfe35…, rowsums -28, -20, 8, 48
```

Two inputs are not regenerated by these commands: the sample of 589
(`../qp shapes/p589.txt l2c=1 l2p=0.02 maxprint=3000 seed=5`, 338 solutions, about 10 minutes, as `runs/589_sample.out`)
and the full H_K = 1 solution set at (5, 11) (7.6 MB, as `runs/5_11_full.out`). `sym589.py`, `unif_test2.py` and
`overlap5_13.py` need them.

**Identities and structure of the sets** (each line is a run; the counts of §90 are in brackets):

```
python3 sym299.py                                   # 2 s    classes 136 Counter({96: 136})   (symmetry group of order 96: multiplier Z13*/H_K, x1<->x2, conjugation, reflection, S3 on B, C, D; all 'bad' counts 0)
python3 claims_all.py 13 19 "[3]" runs/13_19_all.out   # 1 s    solutions 1536, max violation of (J_k),(d_k),(Q): 0, defects off the origin 0, z = 1 for all (Z1 = -1), (c_s, c_t) = (57, 49), pattern ('M8', ('F18','F18','M8')) 1536
python3 verify_identities.py 13 23 "[3]" runs/299_all.out 13   # 19 s   checked 1005 ok 1005 bad 0 max err 4.5e-13
python3 rowid.py 13 23 "[3]" runs/299_all.out       # 1 s    (sum(|R^s|^2-|S|^2), sum(|R^t|^2-|tau|^2)) = (48, 48) for 3072 sampled solutions and (324, 36) for 192 (z = 1: Z = m - 1 = 12)
python3 tA_rel.py 13 23 "[3]" runs/299_all.out      # 2 s    13056; 9984 with t_A = rho conj u_A(mu k), shift 2; 2304 with shift 1; 768 (z = 1) with none
python3 zchar.py 13 "[3]" runs/299_all.out          # 31 s   z = 0: 256 of 12,288 per block Legendre-affine (accidental); z = 1: 0 of 768
python3 g_check27.py runs/13_27_all.out 8           # 1 s    checked 288 solutions; G-level PAF-sum failures 0; A twisted-skew failures 0
python3 overlap13.py                                # 1 s    x0x1x2 class quadruples {'299': 2176, '247': 256, '351': 336}, all pairwise overlaps 0
python3 real2156.py                                 # the defects of the real solution 2156: A none; B 2 orbits (6 positions); C origin + 2 orbits; D origin + 2 orbits; not an origin-defect solution
python3 sstage2.py                                  # 0.6 s  (Lock+S) stage of (A12 | F20, F20, E) at 299: S-solutions 448, 288 are engine level-1 instances, 288 complete = 96 + 192, 160 fail (32 by row sums)
```

**The small members** (`qp_oo`, H_K = 1, `d1=2`; the solution counts are those of `data/runs_small.jsonl`):

```
NUMIN=5 NUMAX=7 QPBIN=qp_oo python3 run_list.py 600 runs/runs_small.jsonl "5:7,5:11,7:11,7:19,7:27,7:31" d1=2 first=1 maxprint=300     # 114 s   solutions of the first chunk 428, 336, 348, 408, 324, 300 (the certificates of skew 280, 616, 1736 are made of the first lines of (5, 7), (7, 11), (7, 31))
python3 rigid_s.py 3 11; python3 rigid_s.py 5 23; python3 rigid_s.py 7 31   # TOTAL (D mod translation) 192 = 4 origin variants x 48 (x 3 translations = 144 per variant); 0; 17,472 = 4 x 4368 (x 7 = 30,576 per variant)
python3 sstage_rigid.py 7 "[1]" 31                  # 5 s   30576 for z = 0 and both signs of eps_B, 0 for z = 1
```

**The family scan** (EXHAUSTIVE per shape in the class `d1=1`; the scan of §90 covers m up to 160 with ν = 5):

```
python3 oo_shapes.py 3 130 7                        # the shape list (first lines: 3 x 7, 3 x 11, 3 x 23, 7 x 11, 7 x 19, 11 x 19, 7 x 31, 11 x 23, 11 x 31, 7 x 59 with 3 orbits)
NUMIN=5 NUMAX=5 python3 run_oo.py 13 60 5 5 300 runs/runs_oo5.jsonl first=1 maxprint=1     # 17 s   18 shapes: HITs 13 x 19 (7681 level-1 instances) and 13 x 23 (9); 13 x 103, 13 x 107, 29 x 47/59/79, 37 x 43/59/151, 49 x 59/67/79/131/199, 53 x 71/107/211 EMPTY (no verdict differs from data/runs_oo5.jsonl)
python3 rowfeas.py 13 23 3                          # 1 s    the (z = 0, 1) row-sum feasibility of the five cost combinations of 299
python3 rowfeas_scan.py 5 160                       # the row-sum test of the origin-defect shapes with 5 orbits (it kills 37 x 43, 53 x 71, … and keeps 13 x 19/23/103/107, 29 x 47/59/79, 37 x 59)
python3 heur.py 60 5                                # the heuristic counts: E_pred 234 at 299, 45.9 at 247, 98 at (7, 11) with H_K of order 3, 0 for 209, 217, 253, 341, 413, 437
python3 make_table.py                               # 62 rows -> runs/results_table.txt
```

**517 = 11 · 47: the quaternary Paley plug-in class over ℤ₁₁ × F₄₇ is empty for every H_K** (EXHAUSTIVE; the reduction
to the rigid pattern is a proof by hand, in [../../NOTES.md](../../NOTES.md) §90):

```
python3 rigid_s.py 11 47                            # 52 s   cost cases [(22, 22)]; (Lock+S) solutions with D canonical: 0 in each of the 4 origin variants; TOTAL 0
python3 rigid_s_test.py                             # 20 s   planted-solution test: found 5 of 5
python3 rigid_s_check.py 3 11 0 1; python3 rigid_s_check.py 7 31 0 1    # full-D (Lock+S) count over all A: 144; 30576
python3 rigid_s_check.py 11 47 6 1                  # 27 s   6 sampled A states: 0     (all 1024 A of origin variant 0: 0, 1831 s under load; rigid_s_check2.py for the other variants: 100 A states each, 0)
```

`rigid_s_check.py` prints a harmless `ComplexWarning` (an unused cast). The independent implementation `rigid_j.py`
agrees with `rigid_s.py` for m = 3, 5, 7 and gives 0 for all four variants at (11, 47):

```
python3 rigid_j.py 3 11     # 0.5 s  144 for each of the four (M, M) variants of rigid_s.py, plus two variants with a full defect on B, 288 each (outside that pattern); TOTAL 1152
python3 rigid_j.py 5 23     # 0
python3 rigid_j.py 7 31     # 1.2 s  30,576 for each of the four variants; TOTAL 122,304
python3 rigid_j.py 11 47    # 0 solutions in all four variants (all 1024 A states of each); 7 min 12 s CPU, 25 min wall under load; the 1.05 million C sequences are tabulated once, 7.7 s
```

## Scope of the verdicts

All verdicts are for the quaternary Paley plug-in with H = H_K × squares(F_q\*), with the stated H_K, or every H_K
for 517. They say nothing about other multiplier groups on the same v.

- EXHAUSTIVE per shape for `d1=1`: complete in the class of origin-defect solutions. For H_K of order ≥ 5 (29 × \*, 37 ×
  59, 43 × 59, 43 × 71, 19 × 79, 13 × 103/107, …) a defect off the origin costs at least |H_K|(q − 3)/2 > 4m, so the
  origin-defect run is exhaustive for the whole plug-in class with that H_K.
- HEURISTIC: the model E(m, p) and the statement that the family is finite.
- UNKNOWN: the shapes with 9 orbits (21 × 43 and 25 × 43, where 8% and 2% of the chunks ran; 25 × 67), 9 × 11, 9 × 23,
  9 × 19 with H_K = 1, 19 × 27 (58% of the chunks, 34,368 level-1 instances, no solution; the model predicted 10
  classes), 31 × 83, 19 × 23, and every shape with 11 or more orbits.

## Notes on the code

- `qp_oo.c` is not stored: it is `../qp.c` with `qp_oo.patch` applied (see Build).
- `common.py`, `cpsat_plug.py`, `sstage.py` and `sstage_rigid.py` import `qshape` from `..`.
- `run_oo.py` and `run_list.py` build shapes with `../qshape.py`, take the engine from `QPBIN` (default `../qp`) and
  create `shapes/`.
- `sym299.py`, `sym589.py`, `sstage2.py`, `overlap13.py` and `overlap5_13.py` read the solution files from `runs/`
  (`299_all.out`, `13_19_all.out`, `13_27_all.out`, `589_sample.out`, `5_11_full.out`); `make_table.py` reads the jsonl
  logs from `data/`.
- `real2156.py` reads `data/c539_hit.json` or a path given as argument.
