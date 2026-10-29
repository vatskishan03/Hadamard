# The quaternary plug-in engine with exact symmetry reductions

The programs behind the quaternary verdicts of [../../NOTES.md](../../NOTES.md) §90: the Paley plug-in of §89 (G = K ×
F_q, q ≡ 3 mod 4, H = H_K × squares, quaternary blocks, A twisted-skew, skew-Hadamard matrices of order 8v through the
negacyclic ring) for the census values that the sweeps of §88 left open, with K = ℤ_m, H_K any subgroup of odd order and
up to nine orbits on K. `qp37q.c` is `../qp.c` with two exact reductions and exact weights, so that the shapes of 7 and
9 orbits, whose lists have 10⁶ to 10⁷ blocks, become finite. The real class with its own reductions is in
[../real37](../real37).

## Contents

| File | Contents |
|---|---|
| `qp37q.c` | `../qp.c` with exact symmetry reductions and exact weights (options below). Every list entry carries its orbit size as a weight, so `WCOUNT` and `l1_weighted` equal the unreduced normalised counts |
| `qshape37q.py` | `../qshape.py` plus, at the tail of the shape file, the permutations of the K-orbits by the units and by the translations through the fixed points of H_K; `python3 qshape37q.py '[m]' q '[H_K generators]' shapes/NAME.txt`. **Not compatible with the shape files of `../real37`** (another tail) |
| `enum37q.py`, `enum37qb.py` | the shapes: K = ℤ_m, q a prime ≡ 3 mod 4 with q ∤ m and q ≤ 8m + 3 (A's half defect (q − 3)/2 ≤ 4m), H_K every odd-order subgroup of ℤ_m\*, 3 ≤ ν ≤ NUMAX orbits, whose defect identity is realisable by the orbit weights, and the origin-defect combinations (`oof`); for the census-undecided v in a range not in `skew/`, or the Table 5 targets (m ≥ 4); `enum37qb.py`: non-cyclic K (ℤ_a × ℤ_b, diagonal units) and prime-power q = pᵉ, e odd |
| `est37.py`, `lister.py`, `run37.py` | the weight-DP bound of the list size (before the spectral cut), the list sizes by `listonly=1`, and the resumable driver: shapes sorted by \|A\|\|L\|, status EMPTY (no solution, `l2_toobig = 0`) / HIT / INCOMPLETE_TOOBIG / TIMEOUT / ERROR in a jsonl log |
| `classical.py` | which classical constructions (Paley I at 8v, 2ᵘ(q + 1), Spence's 4(q² + q + 1), Paley I at 4v) could give a skew order 8v: the v that already have one are not new orders |
| `go.sh`, `chain1.sh`, `d1run.sh`, `t_tsym.sh` | one job; the chain of long runs of §90 (ν = 7 and 9 shapes, and 943); the origin-defect family (`d1=1`: defects only on H_K-fixed orbits); the weighted-count validation of `usym` and `tsym` |
| `data/r36_have.txt` | the v with a skew matrix of order 4v or 8v in the repository when the sweeps of §88–89 ran (the enumerations also read `skew/`) |
| `data/runs/`, `data/p943_b.out`, `data/p943_b.err` | the jsonl logs of the runs, and the log of the 943 run |

Options of `qp37q` in addition to those of `../qp.c`:

- `usym=1`: A canonical under the unit group ℤ_m\*/H_K and the conjugation A ↦ −i·conj A, B, C, D ↦ conj.
- `tsym=1`: B, C, D canonical under translation by the H_K-fixed points, renormalised to x₀(0) = 1.
- `sym=1`: C ≤ D.
- `lcap=N`: capacity of the left table (default 1.6 · 10⁷).
- `est=1`: class-level probe count of the join, without the join.
- `resume=IC`: continue at the C-class row IC (rows are processed in increasing order).

## Build

In this directory (Python 3 with NumPy and SymPy):

```
gcc -O3 -march=haswell -o qp37q qp37q.c -lm
mkdir -p shapes runs
```

`QP=path` overrides the binary name `./qp37q`. A SOL line has the format of `../qp.c`;
`python3 ../plugcert.py m q "[H_K generators]" "SOL ..."` certifies it (order 8mq, `../out/cert_<8mq>.json`).

## Reproduction

Run in this directory after the build. Times are wall times on one core of a machine under load.

**Weighted counts, with and without the reductions** (EXACT: `WCOUNT` of the reduced run = the normalised count of the
unreduced run):

```
python3 qshape37q.py '[13]' 23 '[3]' shapes/p299.txt
./qp37q shapes/p299.txt l2c=1 maxprint=0 usym=1                # 3 s   |A| 8 (64 unreduced), |L| 2303;  COUNT 1632  WCOUNT 13056 l1_weighted=36000   (the 13,056 solutions of skew/2392; unreduced 36,000 level-1 instances)
./qp37q shapes/p299.txt l2c=1 maxprint=0 usym=1 tsym=1 sym=1   # 1.6 s  COUNT 816  WCOUNT 6528 (exactly half: C <-> D)
python3 qshape37q.py '[3]' 7 '[]' shapes/z3_f7.txt
./qp37q shapes/z3_f7.txt l2c=1 maxprint=0 usym=1               # 27 s   WCOUNT 79164864  (the unreduced count of §89)
./qp37q shapes/z3_f7.txt l2c=1 maxprint=0 usym=1 tsym=1        # 1.6 s  |L| 1681 -> 561 (Fix = 3)  COUNT 751728  WCOUNT 79164864
python3 qshape37q.py '[7]' 11 '[2]' shapes/z7_f11.txt
./qp37q shapes/z7_f11.txt l2c=1 maxprint=0 usym=1              # 0.3 s  WCOUNT 3840
```

`./t_tsym.sh` repeats this on ℤ₉ with ⟨4⟩ and q = 23, 19, 11, with and without the reductions. The runs without the
reductions exceed the limit of 300 s; for q = 23 only the reduced count, 549,504, was recorded.

**The shapes: 452 values, 204 shapes, 130 of them realisable**
(`python3 enum37q.py runs/shapes37.json 999 2499 9 census`, 7 s): 46 shapes with 3 orbits, 25 with 5, 36 with 7 and 23
with 9 (K = ℤ_m, q ∤ m). No census value of this list is settled in `skew/`. The weight-DP bound:
`python3 est37.py runs/shapes37.json`.

**The sweep of the shapes with 5 orbits** (25 identity-feasible; EXHAUSTIVE: 23 run, 2 have an empty A list):

```
python3 -c "import json; s=json.load(open('runs/shapes37.json')); json.dump([x for x in s if x['nu']==5], open('runs/shapes37_nu5.json','w'))"
python3 -u lister.py runs/shapes37_nu5.json 40 runs/lists37_nu5.jsonl                                   # 14 s
INCLUDE_PALEY=1 python3 -u run37.py runs/lists37_nu5.jsonl runs/run_nu5.jsonl 900 5 usym=1 tsym=1 sym=1   # 14 s   23 jobs, all EMPTY
```

- Level-1 instances exist for 1363 = ℤ₂₉ × F₄₇ (696), 1711 = ℤ₂₉ × F₅₉ (540), 2291 = ℤ₂₉ × F₇₉ (576 weighted = 72
  reduced without `sym`, 36 with it), 1903 = ℤ₁₇₃ × F₁₁ (86 reduced, 1344 in the unreduced log) and 1603 = ℤ₂₂₉ × F₇
  (10,130). None completes at level 2. With `sym=1` the weighted counts differ from the unreduced ones when instances
  with C = D exist; without it they agree. `data/runs/run_nu5.jsonl` holds the same verdicts with the unreduced counts,
  recorded before the reductions (A = 64 for 2291).
- EMPTY in `data/runs/` (`run_nu3`, `run_nu7`, `run_cheap`, `run_1333`): the shapes with 3 orbits (46: 42 have an empty
  A list, 4 are empty by the join), 7 orbits (11 of 36: 1273, 1333 = ℤ₄₃ × F₃₁, 1397, 1529, 1661, 1793, 1817, 1917,
  2077, 2369, 2449), 9 orbits (1271, 1763 and **943**, a Table 5 target) and the prime-power q = 27 (1269, 1323, 1593,
  1647, 1917).
- Not finished: 1333 = ℤ₃₁ × F₄₃ (36,864 list entries, 120 s), 1133 = ℤ₁₀₃ × F₁₁ (stopped at class 768 of 1911: 640
  level-2 instances, none solves) and 869 (ℤ₇₉ × F₁₁, 3000 s).
- Not run (UNKNOWN: lists of 10⁷ and more entries, pair space 10¹⁴ to 10¹⁵): 2189, 22 shapes with q = 7, 3 (ν = 7), 21
  of the census shapes with 9 orbits, and 779, 759 and 627 among the Table 5 targets.

**943 = 41 · 23, H_K = ⟨10⟩ of order 5, 9 orbits: EXHAUSTIVELY empty** (Table 5: 943(4); one hour of CPU; the recorded
run `data/p943_b.out`: `probes=60046760078 l1_instances=0 COUNT 0 time=3577.5`). The class-level estimate gives the
probe total and its distribution over the rows in 45 s:

```
python3 qshape37q.py '[41]' 23 '[10]' shapes/p943.txt
./qp37q shapes/p943.txt est=1 usym=1 sym=1 lcap=50000000     # lists |A| = 64 (10 classes), |L| = 1,724,685 (1731 classes);  EST rows<=255: 18.9%, <=511: 62.3%, <=767: 99.5%;  EST total probes 6.0047e+10
./qp37q shapes/p943.txt l2c=1 first=1 maxprint=2 usym=1 tsym=1 sym=1 lcap=50000000 resume=1700    # 47 s: the last rows; STATS l1_instances=0 ... COUNT 0 WCOUNT 0 (usym=1 tsym=1 nperm=8 ntrans=1)
./go.sh p943_b shapes/p943.txt l2c=1 first=1 maxprint=2 usym=1 tsym=1 sym=1 lcap=50000000                # the full run (3577 s; detached: setsid nohup ./go.sh ... &); resume=IC continues a stopped run
```

The estimate equals the measured probe count of the run (6.0047 · 10¹⁰ against 60,046,760,078) and its row profile the
progress lines (62.4% after row 512). The 943 verdict is for the quaternary class K = ℤ₄₁, q = 23, H = ⟨10⟩ × squares: a
skew matrix of order 7544 = 8 · 943 with this symmetry does not exist.

**The origin-defect family** (`./d1run.sh 13 19 3 t247`, 0.6 s): the shape of a skew order 1976 = 8 · 247 (247 is
census-decided and 4 · 247 is known, so this is not a new order). `d1=1` restricts the defects to H_K-fixed orbits (the
origin). Result: 961 level-1 instances, 4 reduced solutions (WCOUNT 32), SOL line
`02233|30231|22013 00111|32032|30210 03030|22103|30321 01013|33310|31132`. The family is analysed in
[../origin](../origin).

**Classical routes** (`python3 classical.py 1081 1101 943 2271`): 1081 and 1101 have 8v − 1 a prime power (Paley I), 943
and 2271 do not. The 17 census-undecided v of the shape lists with 8v − 1 a prime power are 1081, 1101, 1141, 1329,
1333, 1363, 1441, 1473, 1561, 1603, 1671, 1711, 2101, 2149, 2433, 2443, 2469: their EMPTY verdicts are not about new
orders.

## Scope of the verdicts

- EXHAUSTIVE for the class K = ℤ_m (or the listed K), q, H_K, H = H_K × squares, quaternary alphabet, when the
  symmetry-reduced join is complete (A canonical under units and conjugation, B, C, D under the fixed-point
  translations, C ≤ D) and every level-1 instance goes to the exact level 2 (`l2_toobig = 0`).
- Validated: the reductions, by the weighted counts above and by independent recomputation of the weighted counts (ℤ₃ ×
  F₇ 79,164,864, ℤ₇ with ⟨2⟩ × F₁₁ 3840, 299 at 13,056); and 35 of the 36 quick EMPTY shapes, re-run unreduced with the
  engine of §89, agree (1661 timed out unreduced).
- UNKNOWN: everything listed above as not finished or not run.

A shape of this class says nothing about other multiplier groups on the same v.

## Notes on the code

- The default binary is `./qp37q` (override with `QP`). `qshape` is imported from this directory.
- The census list, `skew/`, `t5.json` (the copy in `..`) and `data/r36_have.txt` are located relative to the repository.
  `go.sh`, `chain1.sh`, `d1run.sh` and `t_tsym.sh` take their directory from their own location.
