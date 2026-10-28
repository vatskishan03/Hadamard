# tbplain: plain character-bordered cores of orders 358 and 498

Programs behind the quaternary Hadamard matrices of order 358
([../../complex-hadamard/C358.json](../../complex-hadamard/C358.json)) and 498
([../../complex-hadamard/C498.json](../../complex-hadamard/C498.json)). Both come from a core of the
character-bordered Szöllősi array of [../../complex-hadamard/TBORDERED.md](../../complex-hadamard/TBORDERED.md) with t
= 3, in the plain form: the four sequences A, B, C, D take the values ±1 everywhere, including 0 (contrast
[../tbcore/README.md](../tbcore/README.md), the zero-diagonal case). The core of order 358 is found by a tower of
exact levels with sampled fibres of C and D and an exact leaf join. The core of order 498 is found by an exhaustive
search of a band of sparse C and D. Python 3 with NumPy, and gcc; `n498/band.c` needs a CPU with AVX-512F.

Labels: PROVED = argument in the cited document; CHECKED = exact computation or direct measurement, with the command
that does it in this directory; HEURISTIC = estimate or sampling statistics, nothing depends on it for correctness;
UNKNOWN = not determined.

## 1. What was found

A core is four ±1 sequences on an abelian group G of order N − 3, A and B symmetric, all four constant on the orbits
of a multiplier group H, with Σ_k PAF_k(g) = −Γ(φ(g)) for g ≠ 0 and 4|G| at g = 0, where φ: G → Q is onto and Γ = 4
Σ_{χ∈T} χ (TBORDERED.md). With a border (Xb, Rh, Ch) of the compressions it gives K of order 2N with K K\* = 2N I.

| | QH(358) | QH(498) |
|---|---|---|
| N | 179 | 249 |
| G | Z2⁴ × Z11 (176) | Z2 × Z3 × Z41 (246) |
| H | order 3: ω on both F4 factors of Z2⁴ = F4², trivial on Z11 | order 5: z ↦ 10 z on Z41, trivial on Z2 × Z3 |
| Q, φ | Z2² = the second F4 factor | Z3 = the Z3 coordinate |
| T | the three nontrivial characters of Q | all three characters of Q |
| Σ PAF at g ≠ 0 | −12 on ker φ, +4 off ker φ | −12 on ker φ, 0 off ker φ |
| row sums (A, B, C, D) | 18, 0, 14, 14 (Σ r² = 716 = 4N) | 2, 0, 2, 2 (Σ r² = 12 = 4t) |
| PAF columns | 35 (36 classes of ⟨H, −1⟩, minus the class of 0) | 27 (28 classes minus the class of 0) |
| search | tower, C and D sampled at one level | exhaustive band, no sampling |
| time to the hit | 4 minutes | 9 to 10 minutes (the bucket takes 15 minutes) |
| SHA-256 of the canonical text of K | `de730aa0…3f48` | `1164a9fc…5200` |

Both regenerated certificates are byte-identical to the ones in `complex-hadamard/`, and `verify/verify_complex.py`
accepts them (section 5). The matrix of order 358 comes from the first of four hits of the production run (section 3).
The other three hits were regenerated the same way (section 5) and give three more matrices of order 358 (`eaff0478…`,
`931384ac…`, `ac800775…`), which are not in the repository.

## 2. The border and the assembly (`border.py`, `assemble.py`)

The border of TBORDERED.md depends on the compressions only (the fibre sums y_k(a) of the four sequences over the
cosets of ker φ). `border.py` decides it exactly at the level of Q: the rows of Rh (columns of Ch) are the quaternary
vectors of length 2q whose two halves have Fourier support in T, normalised to first entry 1 (324 candidates for Q =
Z2², 1024 for Q = Z3), and a depth-first search in lexicographic order over rows, with bitsets of compatible columns,
finds the first Rh, Ch with Rh Rh\* = Ch\* Ch = 2q I and Xb = −Rh Mc\* Ch / 2q quaternary. It takes 0.1 s for QH(358)
and 0.4 s for QH(498). `assemble.py` checks a core by the definition of the PAF over G (A and B symmetric, Σ PAF =
target at every g), computes the compressions, calls the border search, builds K = [[Xb, R], [R′, Kc]] with Kc = [[X,
Y], [Y\*, −X\*]], checks that all entries are in {1, −1, i, −i} and that K K\* = K\* K = 2N I in integer arithmetic,
and writes `complex-hadamard/C<2N>.json` in the format of C322.json and C358.json, with the SHA-256 of the canonical
text of K (n lines of n digits, digit e for the entry i^e). `verify/verify_complex.py` rebuilds K from the certificate
by its own code.

## 3. QH(358): the tower (`n358/`)

The core is searched through a chain of quotients of G, every level exact except the sampling of C and D at the Z11
level. Block sizes: A and B have 36 cells (the ⟨H, −1⟩-orbits of G), C and D have 66 (the H-orbits); a column is a PAF
class, 35 in all.

**Level 0, Q = Z2² (`qlevel.py`).** The fibre sums of a block over the cosets of ker φ have the form y = (y0, y1, y1,
y1), since H permutes the three nonzero elements of Q; the row sum is r = y0 + 3 y1. A dynamic program over the cells
counts exactly how many sequences of each kind have a given y (2025 values of y per kind, 2³⁶ A-type and 2⁶⁶ C-type
sequences in all), and the two Q-level PAF equations are joined exactly over the four blocks. There are 1440 ordered
Q-tuples in 120 ordered patterns |r_A|, |r_B|, |r_C|, |r_D| (32 up to A ↔ B, C ↔ D); a pattern has one Q-tuple up to
block signs. The pattern used is (18, 0, 14, 14), with pools of 2.83·10⁸, 3.18·10⁸, 5.43·10¹⁷ and 5.43·10¹⁷ sequences
(A, B, C, D). The border exists for it (assemble.py finds it).

**Level 1, E = G/Z11 = Z2⁴ (`rlevel.py`, `rtuples.py`).** The compression r(e) = Σ_z X(e, z) of a block is constant on
the six classes of E (sizes 1, 3, 3, 3, 3, 3), and the number of sequences with a given r is a product over the
classes. For the pattern there are 4444, 4624, 4444 and 4444 vectors r with the Q-compression of the tuple; their PAF
vectors over E (six columns) are joined exactly (dictionary join). There are 56,916 R-tuples, which keep 10^−6.36 of
the product of the pools, and 196 distinct weights (the weight is the number of sequence quadruples of an R-tuple,
10^42.3 for the heaviest). The R-tuples are ranked by weight, one representative (the first in a NumPy sort) per
distinct weight; the production run used ranks 0 to 10. The hit is in rank 1 (index 19593, weight 10^42.20): the
fibres have 10⁶ and 2.5·10⁵ members for A and B and about 10^15.5 and 10^15.3 for C and D.

**Level 2, the Z11 level (`engine.py`, `pjoin.c`).** For z in Z11 let w(z) = Σ_c |P_c| x_c(z), the weighted sum of the
six class values; the PAF of w over Z11 (shifts 0 to 5) gives six more columns, whose target for shift s is the sum of
the target of G over the elements with Z11 coordinate s. The A and B fibres are listed completely and grouped by w
(7350 and 6500 classes for rank 1). For C and D, `ns` = 10⁵ members of each fibre are drawn uniformly (exact sampler
`RSampler`) and grouped by w (99,977 and 99,968 classes). `pjoin` joins the classes exactly on the six columns: 38,323
(R,w)-tuples for rank 1. This is the only sampled step: C and D classes that are not drawn are not searched.
`tuple_priority` orders the (R,w)-tuples by an estimate of their fibre product (A and B exact, C and D by a Gaussian
proxy); the proxy decides the order and nothing else. The hit is in the (R,w)-tuple of priority 1291, with |A_t| = 934
and |B_t| = 467 members.

**Levels 3 and 4, line 1 and line 2 (`tower.c`, `lines.py`).** The sixteen characters of E form six orbits under H;
the transform of a sequence along an orbit is an integer sequence on Z11, a line. Line 0 is w. The 36 class equations
are equivalent to the six PAF equations (shifts 0 to 5) of each of the six lines (a Fourier transform over E).
`tower.c` takes one (R,w)-tuple, lists for C and D the possible values of line 1 for every z (from the letters of
weight w(z); at most 4 values per z), splits the A_t and B_t members by their line-1 values, and joins the classes
exactly on the six line-1 columns. Each match is a child. Inside a child it does the same for line 2 (grandchildren).
The tuple of the hit has 16 children, 114 consistent grandchildren and 104 non-empty leaves.

**Leaf (`tower.c`).** In a leaf the A and B members are listed, and the C and D members are enumerated completely by
depth-first search over the 11 letters with the row sums r and the line values fixed; the PAF class vector of each
member is accumulated from a table K of the 64 × 64 letter pairs. The four lists (here 2·10⁵ C and D members over all
104 leaves, 2·10⁶ pair operations) are joined exactly on all 35 columns: the pair sums (A, C) are hashed, (B, D) are
looked up and every key match is compared on the 35 columns. A match is a core; `produce.py` decodes it, verifies it
by the definition of the PAF over G, and writes it to `hits/`. One (R,w)-tuple takes 0.05 to 0.6 s in the production
log.

**Production run (statistics from the log of the search session; the four hits were re-run one by one, section 5).**
Two processes (seeds 1 and 2, `ns` = 10⁵, the first 3000 (R,w)-tuples of each R-tuple) were stopped after about 35
minutes of wall time on two cores. They had resolved 23,332 (R,w)-tuples of the R-tuple ranks 0 to 10 (2539 s of tower
time; a stage 1 takes about 2 minutes per R-tuple and the member preparation 15 to 48 s) and found four hits, in the
ranks 1 (two), 5 and 6. The four cores give four different matrices of order 358, all with K K\* = 358 I; the first
one is the matrix of `C358.json`. That is one hit per 17 core-minutes. HEURISTIC: a Gaussian model of the leaf density
(10^−13.1 per covered quadruple, not packaged) predicted 11.7 hits; the ratio found to predicted is 0.34, from four
events.

## 4. QH(498): the band search (`n498/`)

**Structure.** G = Z2 × Z3 × Z41, element (a, b, z). H acts by z ↦ 10 z, so each row (a, b) has nine H-orbits ({0} and
eight of size 5), a sequence is six 9-bit words, and ⟨H, −1⟩ has 28 classes. The row sums satisfy Σ r_k² = 4t = 12
with even r_k, so they are ±2, ±2, ±2, 0 in some order; the search uses 0 for B (the other placement of the zero was
not searched), with fibre sums y_A = y_C = y_D = (2, 0, 0) and y_B = 0 after translating C and D by elements of Z3 (a
translation of one block does not change its PAF).

**P and M parts (`plevel.py`, `mcomp.py`).** Put p = (x(0, ·) + x(1, ·))/2 and m = (x(0, ·) − x(1, ·))/2 on Z123 = Z3
× Z41. Then p takes the values −1, 0, 1 and m is a sign exactly where p = 0. The characters of G that are trivial on
Z2 see only p and the others only m, so the 27 equations split into 14 P columns (target 243, −6 four times, 0 nine
times; the 243 is the total support of the four p) and 13 M columns (target 0). The fibre sums force every slice of p
(one per b, nine orbit values) to be a balanced ternary vector: the eight values on the orbits of size 5 sum to 0
(1107 vectors), and the value at 0 is 1 for b = 0 of A, C, D and 0 elsewhere. A C- or D-type object is a triple of
slices (1107³ objects), an A- or B-type object a symmetric slice for b = 0 (19 per centre) and a free slice y whose
reflection is the third slice (21,033 objects). An object with its p fixed has n free signs, the number of orbit
positions where p vanishes (reflected pairs counted once for A and B); it has 2ⁿ completions.

**P-tuples (`band.c`).** A P-tuple is four objects whose P class vectors add up to the P target. The program lists all
of them for one band, C and D with exactly 16 free signs each (`CBAND=16:16`): 18,169,200 canonical C objects and the
same for D, A canonical under the multipliers (5277 of 21033 objects) and all 21,033 B objects. The join is bucketed:
pair sums (A, C) and (B, D) are keyed by the P columns 6, 7, 8, so a bucket is a value s0 of the first key (the sum
over A and C), then (s1, s2); inside a bucket the pairs (A, C) are hashed and (B, D) looked up, with every match
compared on all 14 columns. The buckets are taken in the order s0 = −2, −3, −1, −4, 0, …. The join is exact; the
symmetry reduction is exact (PROVED: the group acting on P-tuples has order 32, generated by the multipliers modulo
⟨H, −1⟩, the independent reflections of C and D, and the swap of C and D; CHECKED by `check_sym.py` on a band: 672
tuples without reduction, 21 with, every orbit met once).

**M completion (`band.c`, `mlib.c`).** A P-tuple with at least 41 free signs in all (the gate, `NMIN`) is completed.
For each of the four objects the 2ⁿ sign vectors are listed in Gray order with their 13 M columns (quadratic forms in
the signs). The lists are pruned by an exact necessary condition: the sum over the four blocks of |m̂_k(χ)|² is
exactly N = 249 for every character χ of Z123 (the characters of G that are nontrivial on Z2 are never lifts of T), so
|m̂_k(χ)|² ≤ 249 − Σ_{j≠k} min_j |m̂_j(χ)|²; members that violate it are dropped, up to six rounds, with AVX-512. The
surviving lists are joined exactly on the 13 columns with target 0, with the pairing of the lists that has the fewest
pair operations. A tuple with a block of more than 22 free signs is written to `.skip` and not completed (none
occurred in the buckets run). Every solution is printed as `SOL` with the object ids and the Gray indices of the
signs; `verify_sol.py` rebuilds the four sequences over G and checks them by the group convolution, independently of
`band.c`.

**What was run.** Bucket s0 = −2 took 15 minutes (28,454 P-tuples, no solution); the hit is in bucket s0 = −3 (26,693
P-tuples, all above the gate, 64 `SOL` lines = 4 P-tuples × 16 sign copies; the first line is the core of the
certificate). The production run was stopped at the first solution, so the rest of the band is UNKNOWN. Restricted to
that bucket, the program is deterministic: the first 32 `SOL` lines are identical to the ones of the original run.

## 5. Files and reproduction

| File | Contents |
|---|---|
| `assemble.py` | core → exact core check → border → K → exact K K\* = 2N I → certificate (structures 358, 498 and the small analogues 102, 18) |
| `border.py` | exhaustive Q-level border search (any small Q and character set T) |
| `reproduce.sh` | build, controls, both hits, comparisons with `complex-hadamard/` and `verify_complex.py` (stages `tests`, `358`, `358alt`, `498`, `sym`, `all`) |
| `n358/structure.py` | G = Z2⁴ × Z11 (and Z2⁴ × Z3), H, orbits, PAF classes, target, exact PAF (FFT and by definition) |
| `n358/qlevel.py`, `rlevel.py`, `rtuples.py` | levels 0 and 1: Q-level pools and tuples, R-level vectors and join, explicit R-tuples, exact sampler |
| `n358/fibres.py`, `lines.py` | letters, cell maps, complete A and B fibres, C and D sampler; lines of E and their targets |
| `n358/engine.py`, `produce.py` | stage 1, tuple priority, tower specs, hit decoding; restartable production driver |
| `n358/pjoin.c`, `pio.py` | exact bucketed four-list join with prefix histogram (stage 1 and controls); its input and output |
| `n358/tower.c` | one (R,w)-tuple resolved completely: children, grandchildren, leaf enumeration and join |
| `n358/test_pjoin.py`, `test_tower.py`, `small_analogue.py` | controls |
| `n498/structure.py`, `plevel.py`, `mcomp.py` | tables of G = Z2 × Z3 × Z41, P-level slices, numpy M part |
| `n498/band.c`, `mlib.c`, `run_band.sh` | the band search, the M part (also a shared object for the test), restartable launcher |
| `n498/verify_sol.py`, `check_sym.py` | exact check of `SOL` lines over G; check of the symmetry reduction |
| `n498/test_mlib.py`, `small_analogue.py` | controls |

Machine of the timings: one Intel Xeon core at 2.1 GHz, `nice -n 12`, `ulimit -v 3000000`, a second job on another
core; NumPy 2.5.3, Python 3.13, gcc 13.3. Everything is run with `OMP_NUM_THREADS=1`. Build:

```
cd search/tbplain/n358
gcc -O2 -march=native -o pjoin pjoin.c
gcc -O2 -march=native -o tower tower.c
cd ../n498
gcc -O2 -march=native -o band band.c -lm
gcc -O2 -march=native -shared -fPIC -o mlib.so mlib.c -lm            # only for test_mlib.py
```

QH(358) (the output directory OUT is anywhere; the verifier is a copy of `verify/verify_complex.py` next to the
certificate):

```
cd search/tbplain/n358
python3 produce.py 18,0,14,14 1 1e5 3000 OUT/p358 2 --only 1291     # 4 min: R-tuple rank 1 (stage 1 2 min, members 40 s), tuple 1291 (0.14 s)
python3 ../assemble.py 358 OUT/p358/hits/hit_r1_t1291_0.json --out OUT/c358 --text      # 1 s
mkdir -p OUT/c358/verify && cp ../../../verify/verify_complex.py OUT/c358/verify/
python3 OUT/c358/verify/verify_complex.py 358       # 358 ok K K* = 358 I  sha256 de730aa01ed6c4387914916d60f39553aad9fd07b72153378e2f76e3747a3f48
```

`produce.py` prints `*** HIT *** r1 t1291 quad [3889, 3410, 63496, 63837] verified=True`; the hit file holds the four
sequences of `C358.json`. The R-tuple rank is a position in a NumPy sort and the (R,w)-tuple rank one in another, so
on a machine with a different sort order they may name other tuples; the same hit is selected by index with `python3
produce.py 18,0,14,14 1 1e5 3000 OUT/p358alt 2 --ridx 19593 --quad 3889,3410,63496,63837`. The other three hits of the
production run come from `python3 produce.py 18,0,14,14 R 1e5 3000 OUT/rR SEED --only T` with (R, SEED, T) = (1, 2,
1292) (the same R-tuple as above), (5, 2, 214) and (6, 1, 1464) (note the seed 1), 3 to 4 minutes each; `assemble.py`
turns each hit into a certificate. The production run was two processes, `python3 produce.py 18,0,14,14
0,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30 1e5 3000 OUT/prodA 1` and the same with ranks `1,3,5,…,31`, seed 2 and
OUT/prodB, stopped after about 35 minutes (ranks 0 to 10 reached).

QH(498):

```
cd search/tbplain/n498
./run_band.sh OUT/p498 b16 16:16 16:16 -3 -3                       # 15 min, 1.1 GB; the SOL lines appear after 9 to 12 min
python3 verify_sol.py OUT/p498/b16.sol --out OUT/p498/verified_solutions.json     # 64 lines, 1 s
python3 ../assemble.py 498 OUT/p498/verified_solutions.json --index 0 --out OUT/c498 --text     # 2 s
mkdir -p OUT/c498/verify && cp ../../../verify/verify_complex.py OUT/c498/verify/
python3 OUT/c498/verify/verify_complex.py 498       # 498 ok K K* = 498 I  sha256 1164a9fcc7fafbcf5ba2994f5b8c9e40ea1b687af927d4909443a47141be5200
```

`run_band.sh OUT NAME CLO:CHI DLO:DHI [S0LO S0HI]` is restartable (it skips the buckets listed in `OUT/NAME.prog`);
without the last two arguments it runs all buckets, from the centre outwards, at about 15 minutes each.

`bash reproduce.sh OUT all` runs the builds, the controls of section 7 and both reproductions with the comparisons (a
certificate in `OUT/c<n>/` equal to `complex-hadamard/C<n>.json`, `verify_complex.py` output); `bash reproduce.sh OUT
358` and `bash reproduce.sh OUT 498` run one hit each.

## 6. What is exact, what is sampled

| Step | Status |
|---|---|
| QH(358) levels 0 and 1: pools, Q-tuples, R-vectors, R-tuples | exact (integer dynamic programs and dictionary joins) |
| A and B fibres of an R-tuple | complete lists (10⁶ and 2.5·10⁵ members for the hit) |
| C and D at the Z11 level | **sampled**: 10⁵ uniform members per fibre, classes by w; a hit is a solution, a miss proves nothing |
| (R,w)-tuple order and the cut at 3000 per R-tuple, ranks 0 to 10 | HEURISTIC (Gaussian proxy); only a restriction of the search |
| inside an (R,w)-tuple: lines, children, leaf enumeration and join | exact (no sampling) |
| each hit | verified by the definition of the PAF over G (`produce.py`, `assemble.py`) |
| QH(498) P-tuples of a bucket | exact join of the complete band lists, symmetry reduction by 32 (checked) |
| M completion | exact: complete sign lists, pruning by a necessary condition, exact join; the gate of 41 free signs is a restriction |
| band and buckets | restriction: C and D with exactly 16 free signs, buckets s0 = −2 and −3 only; the rest of the band and other bands UNKNOWN |
| each solution | verified over G by the group convolution (`verify_sol.py`, `assemble.py`) |
| K | K K\* = K\* K = n I exactly (`assemble.py`; `verify_complex.py` independently) |

## 7. Controls (CHECKED)

| Command | Compares | Result |
|---|---|---|
| `python3 n358/structure.py` | FFT PAF with the PAF by definition, class constancy, on random sequences of both structures | OK |
| `python3 n358/test_pjoin.py` | `pjoin` with a numpy brute force over all quadruples: 12 random instances, 1, 8, 16 and 64 buckets, planted and noisy targets | ALL OK (2 s) |
| `python3 n358/test_tower.py 3e4 0` | `tower.c` at real size: a quadruple planted in a (R,w)-tuple of rank 0 (targets set to its values) | found in the right child, equal to the planted one, exact PAF check (2 min 39 s) |
| `python3 n358/small_analogue.py --out OUT/c102` | the whole chain on G = Z2⁴ × Z3 (N = 51): exhaustive pools, `pjoin` on 11 columns, every core checked over G, border, QH(102) | 15,552 cores in the pattern (0, 2, 10, 10), K K\* = 102 I, `verify_complex.py` ok (13 s) |
| `python3 n498/mcomp.py` | the quadratic forms of the M part with the exact PAF over G, on random objects | OK |
| `python3 n498/test_mlib.py` | `mlib.c` (Gray lists, partial counts, `join4`) with numpy, planted quadruple | OK (1 s) |
| `PLANTM=1 CBAND=20:99 DBAND=20:99 ./band 41 OUT/plant -12 12` | the whole M path of `band.c` on 21 real P-tuples with a random planted completion each | `PLANT ok=21 fail=0` (113 s) |
| `python3 n498/small_analogue.py --out OUT/c18` | border and assembly for Q = Z3 on G = Z2 × Z3: brute-force cores, QH(18), `verify_complex.py` ok | 6912 cores, K K\* = 18 I for the first five (4 s) |
| `bash reproduce.sh OUT sym` | `band.c` without and with the symmetry reductions on band 20:99, s0 in [−12, 12]; `check_sym.py` | 672 tuples without, 21 with, all 21 orbits met once (9 min) |
| `bash reproduce.sh OUT 358alt` | the hit selected by `--ridx 19593 --quad …` with the hit selected by `--only 1291` | same core (4 min) |
| the three other hits (commands in section 5) | regenerated cores, `assemble.py` certificates against the ones recorded in the search session | all three equal, K K\* = 358 I each (3 to 4 min each) |

## 8. Not packaged

The pilots and calibrations that led to the tower and to the band are not part of this directory: the Gaussian
estimates of the number of solutions per tuple and their comparison with exact partial-agreement counts, the pilot
joins at the Q, R and Z11 levels, the cost estimates of a plain four-list join at the Q and R level (10^18 to 10^24
pair operations per hit; HEURISTIC, which is why the tower exists), the first (python-only) version of the tower and
its planted control, the scan of the Q-level border over all 120 patterns of the structure of QH(358) and the
alternative structure Z4² × Z11 for N = 179, and the pilot runs of the band search (other bands, sampled variants).
The census that selected the two structures is not part of this directory either; nothing here depends on it. The
statements of this section and the production statistics of section 3 (tuple counts, times, predicted hits) are from
the log of the search session and were not re-run.
