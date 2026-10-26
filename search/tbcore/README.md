# tbcore: the character-bordered zero-diagonal Szöllősi core of order 478 and the Hadamard matrix of order 3820 = 4 · 955

Programs behind [hadamard/3820](../../hadamard/3820) (n = 955 = 5 · 191, open in Table 4 of Cati–Pasechnik and in [FRONTIER.md](../../FRONTIER.md)) and behind the
quaternary Hadamard matrix of order 1910 on the way (an open complex order, `complex_status.py`). The core is the zero-diagonal case of the character-bordered Szöllősi
array of [complex-hadamard/TBORDERED.md](../../complex-hadamard/TBORDERED.md); the levels (quotient classes, exact lists, exact join) are the idea of [bgschain](../bgschain/README.md).
Python 3 with NumPy and SciPy; a C compiler for `join4.c` (built on first use into the output directory).

Labels: PROVED = argument in the cited document; CHECKED = exact computation or direct measurement in this directory; HEURISTIC = sampling, ranking or estimate only, nothing
depends on it for correctness; UNKNOWN = not determined.

## 1. What it does

**Chain.** Four ±1 sequences A, B, C, D on G = ℤ₂³ × ℤ₂₉ (232 elements) → a quaternary conference matrix K of order 478 = (955 + 1)/2 (entries 1, −1, i, −i and exactly one 0 in every
row and column, K K* = 477 I) → the quaternary Hadamard matrix of order 1910 = 2(q + 2) of the complex Yamada–Kiyasu construction with q = 953, r = 3, ξ = (1, 2)
(`verify.build_complex_yamada_kiyasu`) → the Hadamard matrix of order 3820 by Turyn's doubling H = [[Kr + Ki, Kr − Ki], [−(Kr − Ki), Kr + Ki]]. Only the first arrow is a search.

**Structure** (CHECKED by `check_core.py`, independent of the search code). Index of g = (e₀, e₁, e₂, z): ((e₀·2 + e₁)·2 + e₂)·29 + z. K₀ = 0 × ℤ₂₉, φ: G → Q = ℤ₂³ the projection onto e, T = the 7 nontrivial characters of Q (t = 7),
N = 232 + 7 = 239, core order 2N = 478. H = ⟨h⟩, h = (τ, ×16): (e, z) ↦ (τe, 16z), τ the companion matrix of x³ + x + 1 over GF(2) (order 7, the Singer cycle of GF(8)), 16 of order 7 modulo 29.

* A and B symmetric (x(−g) = x(g)) with A(0) = B(0) = 0, ±1 elsewhere; C, D ±1; all four invariant under H. A, B are constant on the 17 free orbits of ⟨H, −1⟩ (2¹⁷ sequences, enumerated), C, D on the 34 orbits of H (2³⁴).
* Σ_k PAF_k(g) = 926 = 4·232 − 2 at g = 0, −28 on K₀ ∖ {0}, +4 off K₀ (the 18 classes of ⟨H, −1⟩; the class of 0 is automatic, 17 columns are joined). Equivalently Σ_k |x̂_k(ψ)|² = 954 = 4N − 2 for the characters ψ of G that are not lifts of a nontrivial character of Q (the trivial one included) and 26 = 4t − 2 for the 7 lifts.
* Row sums (21, 21, 6, 6) (441 + 441 + 36 + 36 = 954), coset sums over Q: y_A = y_B = (0, 3), y_C = y_D = (−1, 1), written (y(0), y(a) for a ≠ 0).
* The border (14 × 16 and 16 × 14 matrices of quaternary units and a 14 × 14 quaternary conference matrix Xb) is built from the compressions alone and is the same for every solution (section 6). K = [[Xb, R], [R′, Kc]] with Kc = [[X, Y], [Y*, −X*]], X = αX_A + ᾱX_B, Y = αX_C + ᾱX_D R, α = (1 + i)/2, is TBORDERED.md with A(0) = B(0) = 0; the constants 4N and 4t of its proof become 4N − 2 and 4t − 2, which gives K K* = (2N − 1) I. The identity is CHECKED for the matrix built (`assemble_K.py`); the proof is the one of TBORDERED.md.

**Levels.** All exact until the sampled fibres of C and D (section 5).

1. *Q = ℤ₂³* (`structure955.py`, `zq.feasible`). The dynamic program over the orbits allows exactly two row-sum patterns with Σ r² = 954, (7, 21, 8, 20) and (21, 21, 6, 6), and each has a single feasible tuple of Q-compressions. Pattern 1 = (7, 21, 8, 20) has no border of the τ-seed form of section 6 (`border8.py`: 0 seed pairs), so pattern 2 = (21, 21, 6, 6) is the one searched. CHECKED.
2. *ℤ₂₉* (`s1_tuples.py`). Compress also over ℤ₂₉: each block gives a Q-compression and five class values on ℤ₂₉ (the value at 0 and one per orbit of ⟨16⟩ on ℤ₂₉*). The A-type and C-type joint classes are listed with their multiplicities
   (414 and 107,892 classes), the PAF of the ℤ₂₉-compression is reduced to three numbers, and the four lists are joined exactly: **3316 tuples** (mass 10²²·¹³), sorted by weight; the index of a tuple is its rank. The same 3316 come out of the independent dense join `zq.feasible` over the quotient ℤ₂₉ (`structure955.py`). CHECKED.
3. *Fibres* (`s1_tstats.py`). |A| and |B| are 882 or 1960 for every tuple (A-type sequences with the Q-compression (0, 3) and the class of the tuple, by exhaustive enumeration of the 2¹⁷), the fibres of C and D (class and Q-compression (−1, 1)) are counted exactly: 2.6 · 10⁴ to 9 · 10⁶. 196 tuples have |A| = |B| = 882.
4. *Sampler and join* (`s1_prod.py`, `join4.c`). A tuple with |A| = |B| = 882 and both C fibres ≥ 3 L is drawn; A and B are the exhaustive lists, C and D are L = 400,000 uniform samples of the fibres (`s1_prepare.sample_C`); the class vectors (17 columns) of the four lists are joined exactly: pair sums (A, C) hashed by a random linear key, pairs (B, D) looked up, every key match re-checked on all 17 columns, in 2⁵ key chunks (memory). Iteration 0 of the production run: tuple 60, |A| = |B| = 882, fibres 5,882,450 and 3,001,250 (6.8 % and 13.3 % sampled), 1.24 · 10¹⁷ quadruples in 215 s. **The only non-exhaustive step is the choice of the samples of C and D**: a hit is a solution, a miss proves nothing.
5. *Assembly* (`assemble_K.py`, `cyk955.py`). The solution gives K (checked exactly), the Yamada circulants of q = 953 give the quaternary Hadamard matrix of order 1910 (K K* = 1910 I, exact), the doubling gives H (H Hᵀ = 3820 I, exact, and the repository's `verify/check_matrix`).

## 2. Files

| File | Contents |
|---|---|
| `zs.py`, `zq.py` | the model: `ZS` (group ℤ_m₁ × …, multiplier group, orbits, PAF classes by FFT, target, dynamic program for the sign counts and the exact sampler of sign vectors with a prescribed row sum, row-sum patterns, Gaussian estimates), `Quot` (compression to a quotient, exact join `feasible`, conditioned sampler and estimate) |
| `s1_lists.py` | the n = 955 structure (`st`), the exhaustive A-type enumeration (`enum_A`), `build` (orbit values → sequence of G), `zq_compress` |
| `s1_tuples.py`, `s1_tstats.py` | step 2 (3316 tuples) and step 3 (list sizes) of section 1 |
| `s1_prepare.py`, `s1_export.py` | `sample_C`, `classvec`; driver for the lists of one tuple (`.npz`) and its conversion to the binary input of `join4` |
| `s1_prod.py`, `join4.c` | the production sampler and the exact four-list join (C, compiled by `tbio.build_join4` into the output directory with `gcc -O2`) |
| `structure955.py` | row-sum patterns and the exact compression tuples over ℤ₂³ and ℤ₂₉ (the levels above `s1_*`) |
| `border8.py` | the k = 8 border search (all stages in one program) |
| `assemble_K.py`, `cyk955.py` | solution + border → K (order 478) → quaternary Hadamard matrix of order 1910 → H(3820), with the exact checks and the SHA-256 of the matrix and of its canonical CSV |
| `check_core.py` | independent exact check of a solution (group table, no FFT, no use of `zs.py`) |
| `test_join4.py`, `test_planted.py`, `control850.py`, `small_e2e.py` | controls (section 7) |
| `reproduce.sh` | the whole reproduction of section 4 with the comparisons |
| `tbio.py` | output directory, tuple file, on-demand build of `join4` |
| `data/` | the two solutions found (`SOLUTION_0_329_743_226805_104054.json`, `SOLUTION_0_746_137_184467_29180.json`: record `{"tuple", "A", "B", "C", "D", "ok"}`, sequences indexed as above) and `border_k8_pattern2.json` (3 KB each) |

Every program takes `--out DIR` (default: the current directory) and writes only there: `tuples.json`, `trows.npy`, `join4`, the solutions, `K478.npy`, `H3820.npy`, the border files. Nothing is written into the source tree
(Python caches aside). Lists and matrices are regenerated; only the small JSON files above are shipped.

## 3. Usage

```
cd search/tbcore; export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1; ulimit -v 3500000
python3 s1_tuples.py --out OUT                              # 3316 tuples -> OUT/tuples.json
python3 s1_tstats.py --out OUT                              # list sizes -> OUT/trows.npy
python3 s1_prod.py LC MAXIT SEED LCH WALL --out OUT         # sampler + join; a hit -> OUT/SOLUTION_<it>_<a>_<b>_<c>_<d>.json
python3 s1_prod.py LC MAXIT SEED LCH WALL --out OUT --replay IT,A,B,C,D     # lists of iterations 0..IT drawn, join skipped, the given indices taken as the hit
python3 s1_prepare.py TUPLE LC SEED LISTS.npz --out OUT     # the four lists of one tuple (class vectors and sequences)
python3 s1_export.py OUT/LISTS.npz OUT/LISTS.bin            # binary input of join4: ./join4 OUT/LISTS.bin LOG2CHUNKS
python3 check_core.py OUT/SOLUTION_....json                 # independent exact check
python3 assemble_K.py OUT/SOLUTION_....json --out OUT       # K478.npy (border: data/border_k8_pattern2.json, or --border FILE)
python3 cyk955.py --out OUT [--csv-gz FILE]                 # H3820.npy, sha256 of the matrix and of its canonical CSV
python3 border8.py --out OUT                                # the border search (data/border_k8_pattern2.json)
```

`LC` is the size of the sampled C and D lists, `MAXIT` the number of iterations (tuple draws), `SEED` the seed of `numpy.random.default_rng`, `LCH` the log₂ of the number of key chunks of the join (5: about 0.45 GB for `join4`), `WALL` a wall-clock budget in seconds
(checked between iterations). The run stops after the first iteration with a hit. `bash reproduce.sh [OUT] [--border]` does everything of section 4.

## 4. Reproduction (CHECKED; bit-exact)

The production run was `python3 s1_prod.py 400000 4 101 5 900`, so it was seeded: all randomness (the tuple draw, then the samples of C and D) comes from `numpy.random.default_rng(101)`, and a run is determined by (LC, SEED) up to the NumPy version (it depends on the PCG64 stream of `Generator.choice` and `Generator.integers`; NumPy 2.5.3 here; `join4` hashes with a fixed key, which changes nothing but the order of the printed hits). The solution records are shipped in `data/` as well, so nothing below depends on the stream for correctness. From a fresh output directory, one core, `nice`, a loaded machine (load average 7 on 4 cores):

```
python3 s1_tuples.py --out OUT                              #   4.5 s: 3316 tuples, 10^22.13
python3 s1_tstats.py --out OUT                              #   1.5 s
python3 s1_prod.py 400000 4 101 5 900 --out OUT             # 252 s (join 215 s); 161 s (join 136 s) on a less loaded machine
```

```
usable tuples: 882-type 196 1960-type 120
it 0 tuple 60 |A|=882 |B|=882 L=400000 combos 1.24e+17 total 1.24e+17 join 215s sols 2 elapsed 250s
verify: PAF sum classes equal target: True
verify: PAF sum classes equal target: True
```

* `s1_tuples.py` prints the digest of the canonical text of the tuple list, `ec465c5c36d3f45abcec27d596c4aeeb4dddabb1d6305dd59cd95a6bf37f067d` (3316 tuples with their weights); `trows.npy` has sha256 `c1fd29cf45915d9ef3ccb2ad0c5e40c71f25354f449b809ef751d737d64f53d4`.
* Iteration 0 draws tuple 60 = ((−7, −2, 4, −2, 4), (−7, 4, −2, 4, −2), (6, 2, 0, 0, −2), (−8, 0, −2, 0, 4)) (weight 1.37 · 10¹⁹) and the join returns **two** quadruples, indices (a, b, c, d) = (329, 743, 226805, 104054) and (746, 137, 184467, 29180) in the lists A, B (882 each), C, D (400,000 each).
  Both regenerated files are byte-identical to `data/SOLUTION_0_329_743_226805_104054.json` (sha256 `49fac1f1…`) and `data/SOLUTION_0_746_137_184467_29180.json` (`8d68fbfc…`). `check_core.py` accepts both.
* `--replay 0,329,743,226805,104054` (or `0,746,137,184467,29180`) draws the same lists without the join and writes the same file, byte for byte (2.6 s).

```
python3 assemble_K.py OUT/SOLUTION_0_746_137_184467_29180.json --out OUT     # 0.3 s: K K* = 477 I, one zero per row and column, entries 1, i, -1, -i
python3 cyk955.py --out OUT --csv-gz OUT/H3820.csv.gz                        # 12 s
```

The two hits give two different Hadamard matrices of order 3820:

| hit | `K478.npy` sha256 | int8 matrix (row-major) sha256 | canonical CSV sha256 |
|---|---|---|---|
| B = 746, 137, 184467, 29180 | `fdf17efe…8bfa0` | `66f69fad…50c75` | `145a4984…87f49` |
| A = 329, 743, 226805, 104054 | `dae04363…b9374` | `c9b9cc08…e470b` | `79c29912…1482e` |

(full digests in `reproduce.sh`). The first matrix of order 3820 that was verified (`H3820.csv.gz`, CSV digest `145a4984…`) was built from hit B. A gzip file is not a reproducible object (level, header); the decoded CSV is: `--csv-gz` writes it with mtime 0, level 9.
`cyk955.py` prints: `core order 478 C C* = (c-1) I exactly: True`, `Yamada circulants: M^2+N^2 = qI OK`, `K K* = 1910 I exactly: True`, `H H^T = 3820 I exactly: True`;
`gzip -dc H3820.csv.gz | ./check_matrix` (`cc -O2 verify/check_matrix.c`) prints `n = 3820, 7298110 row pairs, 0 violations: H H^T = nI` for both (15 s each).
The canonical CSV is the one of `verify.canonical_csv` (the bytes of `numpy.savetxt(H, fmt="%d", delimiter=",")`).

`bash reproduce.sh OUT --border` runs all of it and the comparisons (hashes above, byte-identity of the two solutions and of the border file) and ends with `reproduce.sh: ALL OK`.

## 5. What is exact, what is sampled

| Step | Status |
|---|---|
| orbit tables, PAF classes, target, row-sum patterns, Q- and ℤ₂₉-level tuples | exact (integer arithmetic; PAF by FFT rounded to integers, cross-checked in `check_core.py` by the group table) |
| A and B lists | exhaustive: all 2¹⁷ A-type sequences, filtered by compression and class |
| fibre sizes of C and D | exact counts |
| C and D lists | **sampled**: 400,000 uniform draws (with repetition) from fibres of 2.4 · 10⁶ to 5.9 · 10⁶ sequences for the 196 tuples used (`sample_C`: the cells of the fibre are drawn with their exact weights, then the orbit signs uniformly within the cell) |
| join | exact (hash of pair sums, every key match compared on all 17 class values; a lost solution is impossible, a false one is rejected) |
| solution | exact verification of the four sequences (`s1_prod.py`: sum of PAF classes = target; `check_core.py`: invariance, symmetry, row and coset sums, PAF over G by the group table, Fourier form) |

The search is RESTRICTED in the sense of bgschain: one iteration joins 882 · 882 · 4 · 10⁵ · 4 · 10⁵ = 1.24 · 10¹⁷ quadruples out of 882 · 882 · 5.9 · 10⁶ · 3.0 · 10⁶ = 1.4 · 10¹⁹ in the fibres of its tuple (0.9 %), and only 196 of the 3316 tuples are used.
The yield is HEURISTIC, from four joins of the shipped code (each 1.24 · 10¹⁷ quadruples, 2.5 to 3.5 minutes on one loaded core; every hit is verified): SEED 101, iteration 0 (tuple 60): 2 hits; SEED 102, iteration 0 (tuple 86): 1 hit, (a, b, c, d) = (774, 328, 25766, 117949);
SEED 103 (`s1_prod.py 400000 4 103 5 1100`): iteration 0 (tuple 108) no hit, iteration 1 (tuple 182) 2 hits, (286, 487, 307357, 379309) and (752, 340, 266335, 224668). Five hits in four joins, about one per join; each of the five solutions passes `check_core.py` and completes with the same border to a Hadamard matrix of order 3820
(five different matrices, H Hᵀ = 3820 I exactly). The sample is too small for a yield statistic. The number of solutions of the whole problem, the yield of the tuples with |A| = 1960 and the yield of larger C and D lists are UNKNOWN.

## 6. The border (CHECKED, `border8.py`)

The border of TBORDERED.md for t = 7, q = 8 depends on the compressions only: a 14 × 16 matrix R̂ and a 16 × 14 matrix Ĉ of quaternary entries, whose rows (columns) have two zero-sum halves on Q = ℤ₂³ (Fourier transform vanishing off T) and are mutually orthogonal, with
X_b = −R̂ Mc* Ĉ/16 quaternary or zero (a quaternary conference matrix of order 14, zeros a permutation matrix). Ansatz: the 7 rows of one orbit are the τ-translates of a seed vector (two seeds for the rows, two for the columns), τ the Singer cycle acting on Q.
`border8.py` (1.5 minutes, deterministic):

```
zero-sum quaternary vectors of length 8: 4900          signatures 105, seeds (h0 | h1) with total autocorrelation 0 and h0[0] = 1: 17640
p1 seed pairs good at all k: 0                         p2 seed pairs good at all k: 8820   (pairs (r, c) with -r Mc* tau^k c / 16 zero or quaternary, k = 0..6)
row seeds 1064, orthogonal row pairs 6272; col seeds 1064, orthogonal col pairs 6272
borders (ansatz tau-invariant, seeds): 28812
R R* = 16 I: True   C* C = 16 I: True   Xb Xb* = 13 I: True   zero pattern is a permutation matrix: True
```

The first border in the order of the search is stored in `data/border_k8_pattern2.json` (sha256 `76e1a3d0…237ab`; the exponents of i of R̂, Ĉ, X_b, and the number 28,812 of seed borders); the regenerated file is byte-identical. `border8.py` also writes `good_p1.npy`, `good_p2.npy` and the 28,812 borders (`borders_p2.npy`) into `--out`. `assemble_K.py` lifts R̂, Ĉ through φ and checks K K* = 477 I, which contains the identities of the proof of TBORDERED.md.
Pattern 1 (row sums 7, 21, 8, 20) has no border in this ansatz, which is why the search is run for pattern 2.

## 7. Controls (CHECKED)

| Command | Compares | Result |
|---|---|---|
| `python3 test_join4.py` | `join4` with a numpy brute force: 60 random lattice-structured instances (1–8 columns, steps, offsets, constant columns, negative values, planted and noisy targets), 1, 4 and 16 key chunks | ALL OK, 1,262,103 solutions compared (17 to 28 s) |
| `python3 test_planted.py` | the production code path (exhaustive A, B of the tuple, C and D lists of 20,000 samples with the recorded C, D appended, class vectors, binary file, `join4`, verification) with the recorded solution of tuple 60 planted | planted (329, 743, 20000, 20000) found at 1 and 32 chunks, 1 solution, verified (9 s) |
| `python3 control850.py` | `zs.py` against an independently found solution: the zero-diagonal core of order 850 of [hadamard/6796](../../hadamard/6796) over ℤ₅² × ℤ₁₇ (H of order 24): orbit structure, target = PAF classes of the four sequences, row sums among `patterns()`, DP counts, sampler | ALL OK (1 s) |
| `python3 check_core.py data/SOLUTION_*.json` | independent exact check of both recorded solutions | ALL CHECKS PASSED |
| `python3 structure955.py` | the exact join over ℤ₂₉ (`zq.feasible`) with the hand-written ℤ₂₉ join of `s1_tuples.py` | 3316 = 3316 tuples (2 s) |
| `python3 small_e2e.py P [U]` | negative control on small analogues G = ℤ₂³ × ℤ_P, N = 8P + 7, H = ⟨(τ, U)⟩: all A- and C-type sequences listed (up to 3.3 · 10⁴ and 6.6 · 10⁴), every row-sum pattern joined by a dictionary join with the same classes and targets (`zs.ZS`) | no solution in any of the 19 cases (P, U) = (3, 1), (3, 2), (5, 1), (5, 4), (7, 1), (7, 6), (9, 1), (9, 8), (11, 10), (13, 3), (13, 4), (13, 9), (13, 10), (13, 12), (15, 2), (15, 7), (15, 14), (17, 4), (17, 13): 163 pattern joins; no positive small case was found, so this exercises the code paths and the exhaustive lists only |
| `python3 border8.py`, `assemble_K.py`, `cyk955.py` | exact identities of the border, K and H listed in sections 4 and 6 | all True |

## 8. Not packaged

The census that selected the group, the multiplier group and the pattern (scans over structures and Gaussian estimates `ZS.estimate`, `Quot`/`cond_estimate`: HEURISTIC rankings, nothing depends on them) is not part of this directory; the estimator code is in `zs.py` and `zq.py` but is not used by the chain.
Diagnostic variants (a join that histograms the length of the agreeing prefix of near misses; a k = 4 border search over ℤ₂²) were not needed for the result and are not included.
