# The exhaustive real plug-in engine with exact symmetry reductions

The programs behind the exhaustive verdicts of [../../NOTES.md](../../NOTES.md) §90 (the real skew class of the Paley
plug-in, G = K × F_q, q ≡ 3 mod 4, H = H_K × squares, blocks ±1, A skew). `qp37r.c` is `../qp.c` (`real=1`) with exact
*existence* reductions of the K-level system and a dry run, so that shapes of 9 to 15 orbits on K that timed out in §89
become finite (EXHAUSTIVE verdicts: the join is complete and every level-1 instance goes to the exact level 2). The
quaternary class with the reductions `usym`, `tsym` is in [../quat37](../quat37); the origin-defect analysis of the
quaternary class is in [../origin](../origin).

## Contents

| File | Contents |
|---|---|
| `qp37r.c` | the engine: `../qp.c` with `real=1` and the reductions and options below |
| `qshape37.py` | `../qshape.py` plus, at the tail of the shape file, the permutations of the K-orbits by the diagonal units of K (`uA`) and by the translations and reflections through the H_K-fixed points (`tF`); `python3 qshape37.py '[m]' q '[H_K generators]' shapes/NAME.txt`, for non-cyclic K the moduli `'[7,7]'` and unit tuples |
| `qshape37m.py` | the same for K = ℤ_p² with a **matrix** group H_K ≤ GL(2, p) of odd order (Singer subgroups of F_{p²}\*, tori): `python3 qshape37m.py p q '[[a,b],[c,d]]' ... shapes/NAME.txt`; `singer(p, order)` gives the generator. Only the diagonal units that normalise H_K are written (see Matrix H_K) |
| `symtest.py`, `symtest2.py` | validation of the reductions by **orbit coverage**: the full solution set equals the union of the group orbits of the reduced one (`symtest.py`: `uA rB sym s3`; `symtest2.py`: `uA tF sym s3`) |
| `unitcheck.py` | every unit permutation of a shape file is a bijection of the K-orbits preserving their sizes, and the set is closed under composition |
| `enum37.py` | shapes the sweeps of `../sweep_real_census.py` do not cover: K cyclic with q \| m, non-cyclic K with diagonal-unit H_K, prime-power q (27, 343, …), for the census values not settled in `skew/` |
| `runq.py`, `verdicts.py` | the resumable queue (one job = shape + time limit; status EMPTY / HIT / INCOMPLETE / TIMEOUT / ERROR in a jsonl log) and the table of its records; `dry2.py`, `probe2.py`: exact level-1 probe totals and level-2 timing probes of a job list; `estl.py`: weight-DP bound of the list sizes |
| `dc_check.py` | the elementary row-sum proof that five q = 3 shapes with 11 orbits (and others whose level 1 is empty) have no solution, independent of any engine |
| `matrix_shapes.sh` | the three matrix-H_K shapes of §90 (575, 1183, 1083) and the control 275, reduced and unreduced |
| `data/RESULTS.txt` | the 66 jobs of the recorded queues (65 EMPTY, 1 TIMEOUT) |
| `data/q*.json` | the job lists that ran: `q1` (ℤ_m with H_K and q = 3, 7, 11), `q2` (ν = 15 with q ≥ 19), `q3` and `q3b_keep` (the new families), `q4_p0` and `q4_p1` (the two halves of 539 = ℤ₄₉ × F₁₁), `q5` (q = 3 with 11 orbits) |
| `data/smoke37.json` | nine quick jobs of these lists |

Options of `qp37r` in addition to those of `../qp.c`:

- `uA=1`: the skew block A canonical under the units of K, which act on all four blocks at once.
- `tF=1`: B, C, D canonical under negation, reflection k ↦ −k and translation by the H_K-fixed points; replaces `rB=1`
  and the normalisation x₀(0) = 1.
- `rB=1`: reflection only.
- `s3=1`: class(B) ≤ class(C) ≤ class(D).
- `sym=1`: C ≤ D.
- `dry=1`: exact level-1 probe totals, no join.
- `prog=SECONDS`: progress lines.
- `part=i nparts=n`: split level 1.
- `resume=IC`: continue at the C-class row IC.

## Build

In this directory (Python 3 with NumPy and SymPy). Every command runs here; shapes go to `shapes/`, logs to `runs/`,
`logs/` and `hits/`:

```
gcc -O3 -march=haswell -o qp37r qp37r.c -lm
mkdir -p shapes runs logs hits
```

`QP=path` overrides the binary name `./qp37r` of the Python drivers. A SOL line of `qp37r real=1` has the format of
`../qp.c`; `python3 ../cert_real.py m q "[H_K generators]" "SOL ..."` turns it into a certificate (writes
`out/cert_<4v>.json` and builds the matrix with the builders of `../../../verify/verify.py`).

## Exhaustive run of one shape

With every reduction (`first=1` stops at the first solution). A finished run with `l2_toobig=0` and no SOL line is an
EXHAUSTIVE verdict for the class: K, q, H_K, real alphabet, solutions up to the symmetries of the K-level system.

```
python3 qshape37.py '[57]' 11 '[4]' shapes/c627.txt
./qp37r shapes/c627.txt real=1 l2c=1 first=1 maxprint=1 sym=1 uA=1 tF=1 s3=1 prog=60
```

The reductions are *existence* reductions: a solution has an image with A minimal in its unit orbit, B, C, D minimal
under negation, reflection and the translations by Fix(H_K), C ≤ D and class(B) ≤ class(C) ≤ class(D). `symtest*.py`
checks them against the unreduced engine on shapes with solutions. The counts of a reduced run are numbers of reduced
solutions, not the normalised counts of `../qp.c`.

## Reproduction

Run in this directory after the build. The comment of a line is the expected output, trimmed. Times are wall times on
one core of a loaded four-core machine (shorter on a free core).

**The reductions, by orbit coverage** (EXACT: the closure of the reduced solutions is the full set):

```
python3 qshape37.py '[3]' 7 '[]' shapes/z3_f7.txt
./qp37r shapes/z3_f7.txt real=1 l2c=1 maxprint=0          # COUNT 22032  (the real count of NOTES §87)
python3 symtest.py shapes/z3_f7.txt                       # 0.3 s   full 22032, reduced 436, closure 22032 equal: True      (also '[3]' 11: 13,824 / 384)
python3 qshape37.py '[5]' 23 '[]' shapes/t_5_23_.txt; python3 symtest.py shapes/t_5_23_.txt    # 5 s   240,000 / 2,800 equal: True
python3 qshape37.py '[5]' 19 '[]' shapes/t_5_19_.txt; python3 symtest.py shapes/t_5_19_.txt    # 8 s   384,000 / 4,480 equal: True
python3 symtest2.py shapes/t_5_23_.txt                    # 7 s   (tF, Fix = 5)   240,000 / 80 equal: True
python3 symtest2.py shapes/t_5_19_.txt                    # 10 s  384,000 / 128 equal: True
python3 qshape37.py '[7]' 31 '[2]' shapes/t_7_31_h2.txt; python3 symtest2.py shapes/t_7_31_h2.txt     # 0.3 s  576 / 12 equal: True
python3 qshape37.py '[13]' 7 '[3]' shapes/t_13_7_h3.txt; python3 symtest2.py shapes/t_13_7_h3.txt     # 4 s    116,352 / 706 equal: True
```

The unreduced and reduced runs of the quaternary engine of `../quat37` give the same weighted counts. `../qskew` (the
orbit join over ℤ_v, no Paley decomposition) agrees with the unreduced counts 22,032 (ℤ₂₁, H = ⟨4⟩) and 13,824 (ℤ₃₃).

**Matrix H_K** (`./matrix_shapes.sh`; K = ℤ_p², H_K a Singer subgroup, H = H_K × squares). diag(u₁, u₂) permutes the
H_K-orbits only if it normalises H_K, which fails for u₁ ≠ u₂ and a Singer group. `qshape37m.py` therefore writes only
the normalising units as the unit action of `uA=1` (275 and 575: 7 of 15 pairs; 1183: 7 of 47; 1083: 3 of 35 pairs), and
`unitcheck.py` checks them.

```
./matrix_shapes.sh
  # 575 (Z5^2 x F23, order 3):    reduced |A| 16 -> 2, no level-1 instance; unreduced |A| = 16, |L| = 768, no level-1 instance    0.02 s
  # 1183 (Z13^2 x F7, order 21):  reduced |A| 144 -> 18 (7 unit perms): 5806 level-1 instances, all fail the level-2 DC test, COUNT 0;  unreduced 1,999,872 level-1 instances, all dcfail, COUNT 0 (24 s)
  # 1083 (Z19^2 x F3, order 45):  reduced |A| 1116 -> 264 (3 unit perms), |L| = 347,034, no compatible class pair, COUNT 0;  unreduced |A| = 1116 (5 classes), |L| = 686,754, COUNT 0
  # ctl275 (Z5^2 x F11, order 3, the class of skew/1100):  COUNT 9216 unreduced, 24 with sym uA tF s3 (7 unit perms); python3 symtest2.py shapes/ctl275.txt: 9216 / 24 equal: True   (44 s)
python3 unitcheck.py 'shapes/m*.txt'                       # checked 3 shapes with unit perms; invalid: 0
```

The unreduced runs of 575, 1183 and 1083 give the same verdicts as the reduced ones; the verdict at 1083 is "no
compatible class pair", not "empty list A". The three verdicts are EXHAUSTIVE: there is no real skew matrix of orders
2300, 4732 and 3252 with this symmetry.

**Quick verdicts of the queues** (nine jobs of `data/smoke37.json`, 19 s):

```
python3 runq.py data/smoke37.json runs/q_smoke.jsonl; python3 verdicts.py
  # 9 jobs; 9 EMPTY; 0 HIT:  539 (Z77 x F7, <4>: 44,578 level-1 instances), 623 (Z89 x F7, <2>: 53,034 / 21,851 level-2), 627 (Z57 x F11, <4>), 759 (Z69 x F11, <4>), 781 (Z71 x F11, <20>: 184 level-1), 783, 959 (Z137 x F7, <16>: 22,464), 1269, 1631 (Z233 x F7, <2>: 4688)
```

The rows of `data/RESULTS.txt` for 623, 781, 959, 1631 and 539/ℤ₇₇ have the same counts. Those of 627 (ℤ₅₇ × F₁₁) and
759 (ℤ₆₉ × F₁₁) were recorded without `tF` (3816 and 8136 level-1 instances, against 156 and 326 with it; the same
verdict).

**539 = ℤ₄₉ × F₁₁ with H_K = ⟨8⟩ (order 7, 13 orbits, Fix = 7): EXHAUSTIVELY empty** (two halves of 25 CPU minutes under
load; 4.4 minutes each on a free core). The dry run and one of 16 slices:

```
python3 qshape37.py '[49]' 11 '[8]' shapes/z49h8_q11.txt
./qp37r shapes/z49h8_q11.txt real=1 l2c=1 dry=1 sym=1 uA=1 tF=1 s3=1     # 49 s   lists |A| = 276, |L| = 507,636;  DRY chunks=5 total_probes_this_part=3.37026e+10  (= 1.6866e10 + 1.6836e10 of the two halves)
./qp37r shapes/z49h8_q11.txt real=1 l2c=1 first=1 maxprint=1 sym=1 uA=1 tF=1 s3=1 part=3 nparts=16     # 48 s CPU: l1_instances=3598 l2_instances=3596 probes=2018142833 l2_solutions=0 COUNT 0
# the full verdict: part=0/1 of nparts=2 (data/q4_p0.json, q4_p1.json, run by runq.py): l1_instances 14,286 + 13,800 = 28,086, all through level 2, 0 solutions, probes 16,866,291,593 + 16,836,277,711
```

**Positive controls** (the whole real pipeline, with every reduction, finds a matrix and certifies it; the control
solutions are not the repository matrices):

```
python3 qshape37.py '[11]' 23 '[]' shapes/ctl_253.txt; ./qp37r shapes/ctl_253.txt real=1 l2c=1 first=1 maxprint=1 sym=1 uA=1 tF=1 s3=1      # 0.06 s  1 level-1 instance, 72 level-2 solutions: skew 1012
python3 qshape37.py '[15]' 19 '[]' shapes/ctl_285.txt; ./qp37r shapes/ctl_285.txt real=1 l2c=1 first=1 maxprint=1 sym=1 uA=1 tF=1 s3=1      # 13 s    1 level-1 instance, 138 level-2 solutions: skew 1140 (also 345: Z15 x F23 -> 1380, 465: Z15 x F31 -> 1860)
python3 qshape37.py '[103]' 3 '[13]' shapes/ctl_309.txt; ./qp37r shapes/ctl_309.txt real=1 l2c=1 first=1 maxprint=1 sym=1 uA=1 tF=1 s3=1     # 7.6 s   q = 3, |H_K| = 17: 526,344 level-1 instances, 28,426 level-2 instances, 4 solutions: skew 1236
python3 ../cert_real.py 11 23 '[]' "SOL 02022200020|00222222220|22000000002 22222022200|02022002202|20200220020 02020222220|20222220202|22000002020 22002000020|20220002202|22002220020"
   # 1.6 s  PAF peak 1012 max off-peak 0, A skew True; order 1012 sha256 d9b617ec...; own check +-1 True HH^T=nI True skew True   (the builders of ../../../verify/verify.py accept it)
python3 ../cert_real.py 103 3 '[13]' "SOL 0220200|0220200|2220200 2220200|2022022|0200022 2200022|0220020|2002202 2220002|2020202|2202020"     # order 1236 sha256 14eb697f...
```

**Row sums alone** (`python3 dc_check.py`, 1 s): `DC-infeasible (explained by row sums alone): [361, 529, 573, 927, 961,
1473, 1713, 1719, 1849, 1893, 2151, 2209, 2433]`. The five q = 3 shapes with 11 orbits (573, 1473, 1713, 1893, 2433) and
eight more have no level-1 instance for a reason visible in the row sums of the blocks (1 + R_B² + R_C² + R_D² = 4v has
no solution with the possible block sums). The other zero-level-1 verdicts (817, 897, 1083, 1705, 1953, 2001, 2277,
2365, 2413, 1407, 1947) rest on the engine.

**The census shape families** (`python3 enum37.py`, writes `shapes37.json` in the current directory; `NUMAX=11` by
default): the shapes that the §89 sweep did not cover, for the census values not in `skew/`:
`28 new shapes Counter({3: 10, 9: 8, 7: 5, 5: 3, 11: 2}) Counter({'A': 17, 'C': 11})`, `21 distinct v` (1.4 s;
`NUMAX=13`: 30 shapes; `NUMAX=15`: 51 shapes for 33 values, including three of the non-cyclic family B). The recorded
queue listed 29 shapes with at most 11 orbits for 22 values (the extra one is 539 = ℤ₇₇ × F₇, settled since by
skew/2156), all EMPTY, seconds each; `data/q3.json` and `data/q3b_keep.json` also hold the shapes with 13 to 15 orbits.

## Scope of the verdicts

- EXHAUSTIVE: the verdicts of `data/RESULTS.txt`, for the real class, K, q, H_K as written, H = H_K × squares, with the
  reductions above.
- Validated by orbit coverage and by unreduced re-runs: 627, 759, 1631, 623, 959, 1183, 539/ℤ₇₇, 575 and 1083; with a
  later build of the engine, also the larger shapes 627/ℤ₂₀₉, 1045, 1407, 1947, 2185/ℤ₁₁₅, 2277, 1659 and 759/ℤ₂₅₃ (in
  full).
- 201 = ℤ₆₇ × F₃: the original run was killed at class 19 of 54; a later complete run is EXHAUSTIVELY empty
  (2.28 · 10⁹ level-1 and 2.5 · 10⁸ level-2 instances).
- UNKNOWN (time, not memory): the shapes with 11 orbits and q = 7 (917, 1057, 1477, 2317), 2185 = ℤ₉₅ × F₂₃, and the
  rigid shapes with 15 orbits whose level 2 needs 2¹⁴ to 2¹⁵ sign fibres per block (2485, 1645, 1209, 2109, 2093, 805).

## Notes on the code

- The Python drivers call `./qp37r` (override with `QP`). `qshape` is imported from the parent directory.
- `enum37.py` (census list, `skew/`) and `dc_check.py` (`data/RESULTS.txt`) locate their inputs relative to the
  repository.
- `runq.py`, `dry2.py` and `probe2.py` create their directories.
