# Skew Goethals–Seidel arrays over the negacyclic ring

The search that found the skew-Hadamard matrices of orders 1528 = 8 · 191 and 2216 = 8 · 277
([../../skew/NEGACYCLIC.md](../../skew/NEGACYCLIC.md); [../NOTES.md](../NOTES.md) §88). It looks for four sequences
u_A, u_B, u_C, u_D on ℤ_v (v odd) with values in {1, i, −1, −i}, constant on the orbits of a multiplier group H of odd
order, with Σ u_X u_X* = 4v δ and u_A(−g) = i · conj(u_A(g)). Each such quadruple is a negacyclic Goethals–Seidel
quadruple of length 2v with a palindromic first sequence, so a skew-Hadamard matrix of order 8v.

| File | Contents |
|---|---|
| `qskew.c` | exact orbit-level join of the four lists (mode 1 quaternary, mode 0 real skew Goethals–Seidel), row sums bucketed, hash join on the exact autocorrelation vectors |
| `qdrv.py` | the shape file of (v, H) for `qskew`, and `run(v, gens, mode, maxprint, firstonly)` |
| `certify.py` | a solution (orbit exponents) to the ±1 negacyclic sequences, the checks, and `out/cert_<8v>.json` |
| `indep_nega.py` | an independent check of a negacyclic certificate (no complex arithmetic) |
| `enum_shapes.py`, `sweep.py`, `t5.json` | the multiplier groups of odd order with at most n orbits for every target v, and the sweep over them (resumable) |
| `qshape.py`, `qp.c` | the Paley plug-in in the same ring (G = K × F_q, H = H_K × squares of F_q*): K-level shape files, and the two-level exact engine (level 1: the compressions to K, class-bucketed hash join; level 2: the t parts, meet in the middle with exact Gaussian integers); `real=1` gives the real skew class of search/skewplug |
| `plugcert.py` | a solution of `qp` over K × F_q to ℤ_v (Chinese remainder theorem), then `certify.py` and the SHA-256: `out/cert_<8v>.json` |
| `gfshape.py` | multiplier shapes over GF(pᵏ) (G = the additive group of F_q, H ≤ F_q* of odd order), run through `qskew2` |
| `dcnp.py` | the row-sum (DC) condition of a shape: `plain_dc`, and `spin_dc` for the spin structure; an exact necessary condition |
| `est.py`, `costest.py` | Gaussian estimate of the number of solutions of a shape (HEURISTIC); exact cost (table entries and probes) of the four-list join of `qskew` |
| `rank_shapes.py`, `dc_rank.py` | the shapes with 15 to 21 orbits ranked by the estimate; the row-sum test of the best |
| `qsa.c`, `cpsat.py` | simulated annealing and a CP-SAT model of the orbit-level system (neither recovers the solution at v = 191) |
| `coverage.py` | the Table 5 entries m ≥ 4 sorted by the least number of orbits of a multiplier shape (77 / 18 / 45 of §88) |
| `spin/` | the spin structure C = B ∘ μ, D = B ∘ μ²: engines `qspin.c`, `qspin_fast.c` (n ≤ 21 orbits), shape lists, row-sum filter, drivers; its [README](spin/README.md) |
| `paley_targets.py`, `paley_dc.py`, `idfeas.py`, `alloc.py`, `paley_shapes.py` | the defect identity of Paley shapes Z_K × F_q: the least number of defect positions, the row-sum condition of the rigid shapes, feasibility of a shape, the defect patterns |
| `pal_lists.py`, `qjoin.c`, `validate_sym.py` | the Paley-pattern pipeline for K with at most 5 orbits: per defect-pattern class explicit block lists and an exact join (row-sum aware cost, symmetry reduction) |
| `g_count.py`, `indep_l1.py`, `indep_l2.py` | checks of `qp`: the G-level count of `qskew` at a plug-in shape; independent numpy reimplementations of level 1 and of both levels |
| `enum_plug.py`, `enum_oof.py`, `census_shapes.py` | the plug-in shapes (K = ℤ_m, q, H_K) of the Table 5 targets and of the census-undecided v; the identity-feasible ones; the origin-defect family |
| `sweep_real.py`, `sweep_real_census.py`, `cert_real.py` | the real-alphabet sweeps of `qp real=1` (Table 5 targets, census values) and the certificate of a real hit (a skew order 4v) |
| `plug_verdicts.sh`, `plug_sampled.sh`, `anat.py` | the exhaustive verdicts of the plug-in shapes of §89, the sampled runs (UNKNOWN verdicts), the anatomy (defect types and costs) of a solution |
| `ncenum.py`, `gshape.py` | shapes over non-cyclic abelian groups G (homocyclic parts, k ≤ 3) with odd-order automorphism groups |
| `octic/` | the octic ring ℤ[x]/(x^{4v}+1) = ℤ[ζ₈][ℤ_v] (order 16v): arithmetic, checks, exhaustive counts, the orbit-level engine `oskew.c`, row-sum feasibility; its [README](octic/README.md) |

Build (in this directory):

```
gcc -O3 -march=haswell -o qskew qskew.c
cp qskew qskew2
gcc -O3 -march=haswell -o qp qp.c -lm
gcc -O3 -march=haswell -o qjoin qjoin.c
gcc -O3 -march=haswell -o qsa qsa.c -lm
mkdir -p shapes out runs
```

The programs of `spin/` and `octic/` are built in their directories (see their READMEs). Python 3 with NumPy and SymPy; `cpsat.py`
needs OR-Tools. All commands are run in this directory; shape files go to `shapes/`, certificates to `out/`, logs and lists to `runs/`
(none of them is committed).

The runs that found the two matrices (exhaustive joins, stopping at the first solution):

```
python3 -c "from qdrv import *; print(run(191, [180], 1, 1, 1)[1])"     # about 4 s
python3 certify.py 191 "[180]" "32032213013 03001230212 02232120000 03303111220"
python3 -c "from qdrv import *; print(run(277, [273], 1, 1, 1)[1])"     # about 3 minutes
python3 certify.py 277 "[273]" "3313220001321 0301322112300 0023212321000 0033103132201"
python3 indep_nega.py out/cert_1528.json
```

`certify.py` writes `out/cert_1528.json` and `out/cert_2216.json`: the certificates of
[../../skew/1528](../../skew/1528) and [../../skew/2216](../../skew/2216) without their `sha256` field (`plugcert.py`
adds it).

The plug-in runs that found [../../skew/2392](../../skew/2392) and [../../skew/4712](../../skew/4712)
([../NOTES.md](../NOTES.md) §89):

```
python3 qshape.py '[13]' 23 '[3]' shapes/p299.txt
./qp shapes/p299.txt l2c=1 maxprint=20000           # exhaustive: 36,000 level-1 instances, 13,056 solutions, 11 s
python3 plugcert.py 13 23 "[3]" "03120|31320|23102 00313|10322|32100 01233|20130|22312 01022|30132|32310"
python3 qshape.py '[19]' 31 '[7]' shapes/p589b.txt
./qp shapes/p589b.txt l2c=1 l2p=0.01 first=1 maxprint=1 seed=3   # level 2 on a 1% sample, first solution after 2.4 s
python3 plugcert.py 19 31 "[7]" "0031210|1123120|0301302 0010322|0031232|0213010 0223120|2100213|0322031 0311200|2013232|2231010"
```

The arguments of `qshape.py` are the moduli of K, q and the generators of H_K; a SOL line of `qp` gives, for each
block, the exponents of i of x₀ (value at (k, 0)), x₁ (on the squares) and x₂ (on the non-squares), one digit per
H_K-orbit of K. Options of `qp`: `l2c=1` (level 2 bucketed by row-sum classes), `l2=0` (level 1 only), `l2p=p`
(level 2 on a random fraction p of the level-1 instances, with `seed=s`), `first=1` (stop at the first solution),
`maxprint=n`, `pA=p pL=p` (thinned lists, sampling), `part=i nparts=n` (split level 1), `real=1` (real alphabet:
skew order 4v).

The sweep: `python3 enum_shapes.py 13 2499` lists the shapes in `enum_n13_v2499.json`, and
`python3 sweep.py 13 2499 1500 runs/sweep.jsonl 1` runs them (1500 s per shape; mode 1 quaternary, 0 real). The list
leaves out the v that the repository settles, so it shrinks as `skew/` grows: 373 shapes when the sweep ran, 371 now,
without 473 and 933, which [../../skew/1892](../../skew/1892) and [../../hadamard/3732](../../hadamard/3732) settle.

## Reproduction of search/NOTES.md §88–89

The programs behind §88 (Negacyclic Goethals–Seidel arrays) and §89 of [../NOTES.md](../NOTES.md) are in this directory, in `spin/` and in
`octic/`. The commands below are run here after the builds above (`spin/` and `octic/`: in their directories); the comment of a line is the output
expected (trimmed). Times are wall times on a loaded core. A sweep line that depends on a time limit (TIMEOUT) is UNKNOWN, not a result. The
K-symmetric shapes and the "census shapes of the 473 type" of §88 and §89 are in [../skewplug](../skewplug). A line marked "expected from the notes" takes hours of CPU; its
command is the one that was run, the number is the one of the notes, and it was not repeated when this block was written.

### §88 Negacyclic Goethals–Seidel arrays

*Validation.* The real mode of `qskew` gives the exact counts of §87; quaternary solutions over ℤ₃ and ℤ₅ give skew matrices of orders 24 and 40.

```
python3 -c "from qdrv import *; print(run(21, [16], 0, 0, 0)[1])"   # COUNT 22032 (x 8 = 176,256)     Z3 x F7
python3 -c "from qdrv import *; print(run(33, [4], 0, 0, 0)[1])"    # COUNT 13824 (x 8 = 110,592)     Z3 x F11
python3 -c "from qdrv import *; print(run(35, [16], 0, 0, 0)[1])"   # COUNT 9360000 (x 8 = 74,880,000), 9 s    Z5 x F7
python3 qdrv.py 3 1 "[]"                                           # real order 24 (True, True, True)
python3 qdrv.py 5 1 "[]"                                           # real order 40 (True, True, True)
```

*Sweep* (373 shapes when it ran, 371 now; the real class on 510, now 508): the shape lists are `enum_shapes.py`, the sweeps `sweep.py`
(commands above).

```
python3 enum_shapes.py 13 2499     # first line: 278 v with shapes 371 shapes     (2 to 4 minutes; the listing after it is long)
python3 enum_shapes.py 17 2499     # first line: 364 v with shapes 508 shapes
```

*Row sums.* The primes v = 1 + k(k + 1)/2 with k ≡ 3 (mod 4):

```
python3 -c "from sympy import isprime; print([v for k in range(3,60,4) for v in [1+k*(k+1)//2] if isprime(v)])"   # [7, 29, 67, 191, 277, 379, 631, 947, 1129, 1327]
```

*GF(pᵏ) shapes* (361, 841, 1369, 1849, 2197 empty): `gfshape.py` with no argument runs the six shapes (the last one, 13³ with |H| = 183 and 13
orbits, takes 80 s; each of the others 2 s).

```
python3 gfshape.py 19 2 45         # 361 H 45 n 9 mode 1 ['INCOMPLETE_TUPLES 0', 'COUNT 0 (normalized: ...']
python3 gfshape.py                 # the same COUNT 0 for 1369 (H 171), 1849 (231), 841 (105), 2197 (549, n 5) and 2197 (183, n 13)
```

*Paley products in this ring.* The least number N of positions at which the defect weights can add up to 4|K| (with A's forced half defect): single
defect N = 1 (413 = 7 · 59: only A's own half defect fits, dead), the rigid shapes with N = 2, q = 4|K| − 1 or 4|K| + 3: 473, 517, 1139, 1207.
`qp` itself finds the A list of 413 empty (below).

```
python3 paley_targets.py 413 473 517 1139 1207
# (413, 4, 59, 7, None, 1)  (473, 4, 43, 11, 1, 2)  (517, 4, 47, 11, 1, 2)  (1139, None, 67, 17, 1, 2)  (1207, None, 71, 17, 1, 2)    (v, Table 5 m, q, K, N real, N quaternary)
python3 paley_dc.py      # 473 11 43 n 33 DC feasible True; likewise 517, 1139, 1207, 299, 589, 413, 575, 35 (40 s): the row sums do not kill the rigid shapes
```

*Spin shapes* (C = B ∘ μ, D = B ∘ μ²; at most 19 orbits they fail the row-sum condition or are empty (817); 589 with ⟨574⟩, μ = 5 at 21 orbits is empty).
In `spin/` after its build:

```
(cd .. && python3 enum_shapes.py 21 2499)       # the multiplier shapes with at most 21 orbits (2 minutes)
python3 enum_spin.py 17                          # 114 Counter({7: 42, 17: 38, 13: 34})
python3 spin_est.py && python3 spin_filter.py                     # n = 17: "2 of 38" pass the row sums: 817 <814> and 817 <769> (mu 4)
python3 spin_sweep.py 7 13 900 runs/spin_sweep.jsonl runs/spin_shapes_n17.json     # 76 shapes (42 + 34): all status 'empty', 5 s
python3 spin_sweep.py 14 17 3000 runs/spin_sweep17.jsonl runs/spin_rank_dc.json   # both 817: 'COUNT 0 (part 0/1 passDC 6830208 ...)', 3 s each
python3 spin_drv.py 589 "[574]" 5     # COUNT 0  (part 0/1 passDC 61471872 pass1 3663072); 2 minutes (setup 70 s)
```

(`passDC` and `pass1` are the numbers of B that pass the row-sum test and the first class filter; `qspin_fast` equals `qspin` on v = 7, 13.) At n = 19
all 13 shapes fail the row sums, and at n = 21 18 shapes pass (20 of the 128 with n = 17 to 21) (`enum_spin.py 21; spin_est.py 21; spin_filter.py` below).

*Coverage of Table 5* (172 entries m ≥ 4; 463 not census-undecided; 29 with a matrix of order 4v; 142 others; 191 and 277 the hits; of the remaining 140 :
77 with a shape of at most 13 orbits, 18 with only 15 to 21, 45 with none):

```
python3 coverage.py
# entries with m >= 4: 172 / not census-undecided: [463] / with a matrix of order 4v in skew/ (at the time of the sweep): 29 / others: 142  hits among them: [191, 277]  remaining: 140
# a multiplier shape with at most 13 orbits: 77 / only shapes with 15 to 21 orbits: 18 [295, 299, 335, 337, 401, 515, 539, 543, 611, 655, 667, 695, 707, 851, 881, 899, 927, 995]
# no shape with at most 21 orbits: 45 [225, 245, 257, 265, 305, ..., 979]
python3 enum_shapes.py 21 2499 && python3 rank_shapes.py && python3 dc_rank.py    # the shapes with 15..21 orbits: only 101 is rich (+37 bits) and has a row-sum pattern; 337, 461, 401, 487 die by the row sums
python3 costest.py 191:[180] 277:[273] 299:[289]   # probes 8.07e9 (n 11), 7.91e12 (n 13), 6.84e15 (n 15)
```

*Simulated annealing and CP-SAT do not recover the solution at v = 191* (11 orbits):

```
python3 -c "from qdrv import *; shape(191, [180], 'shapes/t_191.txt')"
./qsa shapes/t_191.txt 1 60 300 3 2000000                          # runs 9 bestf 72: no SOL line (seed, seconds, T0, T1, moves per run)
python3 cpsat.py 191 "[180]" 240                                   # status UNKNOWN 240.2 s: no solution
```

### §89 Paley plug-ins in the negacyclic ring; the octic ring

*Engine validation* (normalised counts of `qp` against `qskew` over ℤ_v; three empty shapes; the real alphabet):

```
python3 qshape.py '[3]' 7 '[]' shapes/z3f7.txt; ./qp shapes/z3f7.txt maxprint=0       # COUNT 79164864   (160 s; qskew2: python3 g_count.py 3 7 "[]", 13 min)
python3 qshape.py '[3]' 11 '[]' shapes/z3f11.txt; ./qp shapes/z3f11.txt maxprint=0    # COUNT 3753216    (9 s)
python3 qshape.py '[7]' 11 '[2]' shapes/p77.txt; ./qp shapes/p77.txt maxprint=0       # COUNT 3840       (1 s; python3 g_count.py 7 11 "[2]": COUNT 3840, 50 s)
./qp shapes/z3f7.txt real=1 maxprint=0                                                # COUNT 22032
./qp shapes/z3f11.txt real=1 maxprint=0                                               # COUNT 13824
python3 qshape.py '[5]' 7 '[]' shapes/z5f7.txt; ./qp shapes/z5f7.txt real=1 maxprint=0   # COUNT 9360000   (45 s)
# the three empty shapes: COUNT 0 for  qshape.py '[3]' 19 '[]',  '[3]' 23 '[]',  '[7]' 19 '[2]'
python3 validate_sym.py 77 7 11 "[36,23]" v77s 3840           # weighted total 3840 expected 3840 OK  (the symmetry reductions of pal_lists.py)
PAL_NOSYM=1 python3 pal_lists.py 21 3 7 "[16]" 3 t21         # TOTAL ordered-pattern count ...: 79164864   (the pattern pipeline without symmetry halving, 1.5 to 2.5 minutes;
                                                               # with the halving the total is 33882408, the weighted count of validate_sym.py is the full one)
```

The independent reimplementations (numpy, no code shared with `qp.c`): level 1 and both levels at ℤ₁₃ with ⟨3⟩ × F₂₃, level 1 at ℤ₁₉ with ⟨7⟩ × F₃₁.

```
python3 indep_l1.py 13 23 "[3]" --psd                          # LEVEL1_INSTANCES (ordered, B,C,D normalised x0(0)=1): 36000  (hash false positives: 0), 5 s
python3 indep_l2.py 13 23 "[3]" --out runs/p299_indep.txt      # level-1 instances 36000; LEVEL2: ... completing instances 4800 solutions 13056   (6 minutes)
./qp shapes/p299.txt l2c=1 maxprint=20000 > runs/p299_qp.out   # after  python3 qshape.py '[13]' 23 '[3]' shapes/p299.txt ; the same 13,056 SOL lines as a set:
sort runs/p299_indep.txt | grep ^SOL | cmp - <(grep ^SOL runs/p299_qp.out | sort) && echo equal
python3 indep_l1.py 19 31 "[7]" --psd                          # after psd: |A|= 256 |L|= 36863; LEVEL1_INSTANCES ... 261984 (hash false positives: 0)   (13 minutes unloaded)
```

*Hits* (the commands of the top of this file): `qp` gives for 299 `l1_instances=36000 ... l2_solutions=13056` (11 s) and for 589 `l1_instances=261984`
with 84 of 2,540 sampled instances completing (168 solutions) when it is run on the 1% sample `l2p=0.01 seed=3` without `first=1`. The pattern pipeline
finds the 299 class independently: `PAL_CLASSES=2 python3 pal_lists.py 299 13 23 "[185, 209]" 5 p299` (class 2 of 9, the pattern of the certified
solution: `"count": 80`, about 5 minutes; `PAL_MAXPRINT=100` lists them all, and the certified solution is among them up to units of ℤ₂₉₉, reflections of
B, C, D and their order; class 0 is another pattern, with 32 solutions).
The anatomy (every defect at k = 0):

```
python3 anat.py 13 23 "[3]" "03120|31320|23102 00313|10322|32100 01233|20130|22312 01022|30132|32310"   # A W= 10, B W= 0, C W= 20, D W= 22; sum W 52 4|K| 52
python3 anat.py 19 31 "[7]" "0031210|1123120|0301302 0010322|0031232|0213010 0223120|2100213|0322031 0311200|2013232|2231010"   # A 16, B 32, C 0, D 28; sum W 76
```

*Other plug-in shapes, EXHAUSTIVELY empty* (K = ℤ_m, q, H_K; `plug_verdicts.sh` runs `qshape.py` and `qp l2c=1` for each name and prints the `STATS` line,
whose `l1_instances` is the number of level-1 solutions and `l2_solutions` the number of solutions):

```
./plug_verdicts.sh 667 851 899 589z31 817a 817b 989 413a 413b   # l1_instances=0 and COUNT 0 for each (667, 851, 899, 817b, 413a, 413b: at once; 817a: 8 s; 989: 2 s;
                                                               # 589z31: 94 s CPU, probes=2336393706)
./plug_verdicts.sh 707      # l1_instances=290880 l2_solutions=0 COUNT 0         (2 minutes)
python3 qshape.py '[49]' 11 '[2]' shapes/v539.txt; ./qp shapes/v539.txt l2=0 maxprint=0   # level 1 only: l1_instances=10856160 (probes 1.8e9, 60 s CPU)
./plug_verdicts.sh 539      # both levels: l1_instances=10856160 ... l2_solutions=0 (level 2 takes hours; not rerun here)
```

The names are 667 (ℤ₂₉, 23, ⟨7⟩), 851 (ℤ₃₇, 23, ⟨7⟩), 899 (ℤ₂₉, 31, ⟨7⟩), 589z31 (ℤ₃₁, 19, ⟨2⟩), 817a (ℤ₄₃, 19, ⟨4⟩), 817b (ℤ₁₉, 43, ⟨7⟩), 989 (ℤ₄₃, 23, ⟨4⟩),
413a and 413b (ℤ₇, 59, H_K = 1 and ⟨2⟩), 707 (ℤ₁₀₁, 7, ⟨5⟩), 539 (ℤ₄₉, 11, ⟨2⟩); `299` and `589` run the two hits. The row-sum verdict of 667 and 851
(every defect pattern dies by the row sums), and of 1403, in the pattern pipeline (`PAL_DRY=1`: the cost 0 of a class means no row-sum combination):

```
PAL_DRY=1 python3 pal_lists.py 667 29 23 "[645, 117]" 5 d667      # ordered patterns 444 classes 21; every class "cost3b+p 0"; TOTAL ... 0
PAL_DRY=1 python3 pal_lists.py 851 37 23 "[599, 186]" 5 d851      # ordered patterns 336 classes 16; every cost 0
PAL_DRY=1 python3 pal_lists.py 1403 61 23 "[622, 1221]" 5 d1403   # ordered patterns 288 classes 13; every cost 0
```

*EXACTLY empty.* No defect pattern of weight 4|K| (K = ℤ₅ with q = 59, 67, 103, 131, 139, 199 for 295, 335, 515, 655, 695, 995; ℤ₉ with q = 103 for 927;
ℤ₁₃ with q = 47 for 611, also 899 = ℤ₂₉ × F₃₁ with H_K = ⟨7⟩):

```
python3 idfeas.py 5 59          # K 5 q 59 |HK| 1 kappa 5 sizes [1] infeasible          (likewise 5 67, 5 103, 5 131, 5 139, 5 199)
python3 idfeas.py 9 103 4       # K 9 q 103 |HK| 3 kappa 5 sizes [1, 3] infeasible
python3 idfeas.py 13 47 3       # K 13 q 47 |HK| 3 kappa 5 sizes [1, 3] infeasible
python3 idfeas.py 29 31 7       # K 29 q 31 |HK| 7 kappa 5 sizes [1, 7] infeasible
python3 idfeas.py 13 23 3       # IDENTITY-FEASIBLE (the hit 299);  python3 idfeas.py 7 59: IDENTITY-FEASIBLE (413: only A's own half defect)
python3 enum_shapes.py 21 2499 && python3 paley_shapes.py && python3 idfeas.py   # all kappa <= 7 shapes: the kappa 5 lines of 299, 539, 667, 707, 851, 1363, 1403, 1711, 2291 are IDENTITY-FEASIBLE
python3 alloc.py 295 5 59 "[$(python3 -c 'from qshape import cyc_gens; print(*cyc_gens(5,59,[]))')]" 5   # defect-config quadruples with weight sum = 4K: 0
```

The census beyond Table 5: 1363 = ℤ₂₉ × F₄₇ and 2291 = ℤ₂₉ × F₇₉ (both with H_K = ⟨7⟩) are EXHAUSTIVELY empty (all 8 and 4 classes searched, 0 solutions; about
30 CPU minutes each), and 1403 dies by the row sums (above). The pattern pipeline, one class of 1363 as a check (class 6: 1 to 4 minutes, count 0):

```
PAL_CLASSES=6 python3 pal_lists.py 1363 29 47 "[471, 871]" 5 t1363     # {"cls": 6, ... "count": 0, ...}
python3 pal_lists.py 1363 29 47 "[471, 871]" 5 t1363                    # all 8 classes, count 0 in each:  TOTAL ordered-pattern count ...: 0   (expected from the notes, 30 CPU minutes)
python3 pal_lists.py 2291 29 79 "[1660, 88]" 5 t2291                    # all 4 classes, count 0                                                  (expected from the notes)
```

(1711 = ℤ₂₉ × F₅₉: three of six classes were run, UNKNOWN.) *Sampled without a hit* (UNKNOWN; thinned lists, every job at most 15 or 25 minutes, seeded):

```
./plug_sampled.sh 1      # 869, 737, 943, 779, 759b, 627b: pL = 0.05 .. 0.01; no SOL line      (expected from the notes: 6 jobs of up to 15 minutes)
./plug_sampled.sh 2      # the same with block B all-E (bE=1), and 721, 973 (q = 7); no SOL line
./plug_sampled.sh 3      # larger samples of 721 and 973; no SOL line
./plug_sampled.sh 299    # control, 1 s: the shape of the hit with bE=1 (no thinning): runs/sampled299.log has l1_instances=2784 ... l2_solutions=2560 and *** HIT 299
```

*Origin-defect shapes* (4|K| = A's half defect + at most three origin costs):

```
python3 enum_oof.py | tail -1                       # 9 shapes; distinct v: 4     (299, 413, 517, 589)
PLUG_SETTLED=all python3 enum_oof.py | tail -1      # 21 shapes; distinct v: 11  (also 209, 253, 341, 345, 391, 473, 903, settled in the repository)
```

*The real class* (`qp real=1`, skew order 4v; the lists of targets and census values are the repository without skew/2392 and skew/4712, which `qp` found later):

```
python3 enum_plug.py                                 # 45 shapes 28 identity-feasible  Counter({9: 13, 7: 10, 5: 5})   (shapes with 5 to 9 K-orbits)
NUMAX=5 TMO=60 python3 sweep_real.py                 # the 54 shapes with 3 to 5 orbits: all empty, 3 s (the same STATS lines as the full sweep)
NUMAX=10 TMO=120 python3 sweep_real.py               # 85 shapes with 3..10 orbits: 81 empty, 4 TIMEOUT (623, 627, 759, 959), no HIT      (expected from the notes, about 1 hour)
NUMIN=11 NUMAX=15 TMO=30 LOGF=runs/real_hi.jsonl python3 sweep_real.py    # 21 shapes: empty 517 649 737 869 913 (Z_11, trivial H_K), 299 611 (Z_13), 589 (Z_31, <5>), 851 (Z_37, <10>); the rest TIMEOUT   (expected from the notes)
NUMAX=10 python3 census_shapes.py; NUMAX=15 python3 census_shapes.py   # 643 census values not in skew/; 398 shapes (nu 3: 182, 5: 77, 7: 88, 9: 51) and 549 (+ 11: 41, 13: 27, 15: 83)
NUMIN=7 NUMAX=10 TMO=40 python3 sweep_real_census.py    # 139 shapes: 128 empty, 11 TIMEOUT, no HIT       (expected from the notes; its first shapes: 329, 413, 489, 553, 581 empty)
NUMIN=11 NUMAX=15 TMO=20 LOGF=runs/real_census_hi.jsonl python3 sweep_real_census.py    # 151 shapes: 65 empty, 86 TIMEOUT, no HIT   (expected from the notes)
python3 qshape.py '[3]' 7 '[]' shapes/z3f7.txt && ./qp shapes/z3f7.txt real=1 maxprint=1     # a real solution (SOL line, 4v = 84)
python3 cert_real.py 3 7 "[]" "<the SOL line>"       # a control: order 84, own check: +-1 True HH^T=nI True skew True
```

*Spin shapes with 21 orbits* (in `spin/`; each is EXHAUSTIVELY empty, `COUNT 0`; the numbers are `passDC` and `pass1`; setup 45 to 80 s, then about 0.2 µs per B that passes the
row sums). The shape lists are `python3 enum_spin.py 21; python3 spin_est.py 21; python3 spin_filter.py; python3 spin_cost.py` (the DC-feasible shapes at n = 21 and their `passDC`: 589 6.15e7, 329 1.73e10,
737 7.68e8, 589 <577,514> 1.31e9, 869 2.05e9, 817 <815> 5.5e8, 413 5.19e9, 817 <809,767> 4.78e8, ...; `spin_cost.py` takes about half an hour for all):

```
python3 spin_drv.py 589 "[574]" 5          # passDC 61471872 pass1 3663072        (v = 589, H of order 45; 2 minutes)
python3 spin_drv.py 589 "[577,514]" 4      # passDC 1306277280 pass1 28.8M        (expected from the notes; about 10 minutes)
python3 spin_drv.py 737 "[729]" 4          # passDC 768398400 pass1 111218724     (expected from the notes)
python3 spin_drv.py 817 "[815]" 6          # passDC 549831744 pass1 268690264     (expected from the notes)
python3 spin_drv.py 817 "[809,767]" 4      # passDC 478114560 pass1 583548        (7 minutes under load)
python3 spin_drv.py 869 "[852]" 4          # passDC 2049062400 pass1 97.4M        (expected from the notes)
python3 spin_drv.py 1273 "[1261,1255]" 4   # passDC 51226560                      (expected from the notes)
python3 spin_drv.py 1273 "[1270]" 4        # passDC 215151552                     (expected from the notes)
python3 spin_drv.py 1397 "[1357]" 9        # passDC 768398400 pass1 14252136      (expected from the notes)
python3 spin_drv.py 1661 "[1620,1076]" 4   # passDC 768398400 pass1 6137172       (expected from the notes)
```

*The octic ring* (order 16v; in `octic/` after its build):

```
python3 verify_all.py                          # 23 PASS lines, ALL PASS (17 s): ring isomorphism, the 16-letter alphabet, moduli 4 +- 2 sqrt2, the peak equation, skew <=> palindromic
./exsmall 4 0 0                                # TOTAL_ORDERED_FULL 6144           (v = 1)
./exsmall 12 0 0                               # TOTAL_ORDERED_FULL 2057674752     (v = 3; COUNT 148848)
./exsmall 20 0 0                               # TOTAL_ORDERED_FULL 75652141056000 (v = 5; COUNT 1182064704; 7 minutes, 1.1 GB)
python3 cert_exh.py 4; python3 cert_exh.py 12; python3 cert_exh.py 20    # skew matrices of orders 16, 48, 80, built and checked by verify/verify.py (n = 20: 7 minutes)
./exsmall 12 0 1                               # DOUBLING_TEST tested 148848 doubled 10368   (7.0%: the solutions made of two pairs (X, m theta X))
./exsmall 20 0 1                               # DOUBLING_TEST tested 1182064704 doubled 521856   (0.044%; 24 minutes under load, 15 unloaded)
python3 odrv.py 61 "[12]"                      # COUNT 4352 (class-ordered ...; ORDERED total of (A,B,C,D) with B,C,D normalised by w^k = 15360 ...), 2 s
python3 -c "import sys; sys.path.insert(0,'..'); from qdrv import shape; shape(61,[12],'shapes/f1_61_12.txt')"; ../qskew shapes/f1_61_12.txt 1 0   # COUNT 0: the order-8v family has no solution there
python3 cert_oct.py 61 "[12]" "09c27 0c325 2ed30 259e2"     # ... all checks True, order 976: out/o976 (a control: the order is known)
python3 targets.py; python3 shapes_oct.py 65; python3 dcfeas.py 65    # 61 targets; 177 shapes with n <= 65; 83 dc_infeasible, 63 dead_selfneg, 22 dc_feasible, 9 dense_not_evaluated
python3 dcpat.py 337 '[2]'                     # DC class patterns: 44, every large block has a coordinate +-41 = 2*21 - 1
python3 oest.py 61 "[12]"                      # logN 3.95 (HEURISTIC; the exact log2 is 13.91)
```

The last three lines give the row-sum statement of §89: with |H| = h a block has R = U(0) + hS with S of even coordinates, so Σ|R_X|² = 16v needs (2h − 1)² + 15 ≤ 16v;
all shapes with at most 14 orbits of the 61 targets are infeasible (`python3 -c "import json; r = json.load(open('runs/dcfeas_n65.json')); print(min(x['n'] for x in r if x['status'] in ('dc_feasible', 'dense_not_evaluated')))"`
prints 15, the dense shapes of 667 and 851; the least dc_feasible n is 17 at 337), and the feasible ones (17 orbits at 337, 21 at 101 and 821, 23 at 947, ...) need lists of 2·16^(n−1) ≥ 2^57 blocks.
