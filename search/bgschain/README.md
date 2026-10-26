# bgschain: the exact compression chain for bordered Goethals–Seidel quadruples of any index

Programs behind [hadamard/3356](../../hadamard/3356) (n = 839, index 8, G does not split over the subgroup) and the regression on
[hadamard/3964](../../hadamard/3964) (n = 991, index 4); see [../NOTES.md](../NOTES.md) §18–19 and §102. Python 3 with NumPy; a C compiler
for `hcore.c` (built on first use into `libhcore.so`) and for `e2e/z7search.c` (built on demand into a temporary directory).

Labels: CHECKED = exact computation or direct measurement in this directory; HEURISTIC = ranking or estimate only, nothing depends on it for correctness;
RESTRICTED = exact only inside a stated restriction of the search (the search is then not exhaustive).

## 1. What it does

**Problem.** G = ℤ_{m₁} × … × ℤ_{m_r} (mixed radix, first factor slowest, as `verify/verify.py`), S ≤ G of index k (`S_gens`), Q = G/S with the cosets of S
numbered in order of their smallest element (the numbering of `verify.build_bordered_gs`). Four ±1 sequences x₀, …, x₃ on G, constant on the orbits of
H = ⟨gens⟩ (integer matrices acting on coordinate columns, row i reduced mod m_i; at most 64 orbits), with prescribed coset sums y_i (k numbers each) and

```
Σᵢ PAF_xᵢ(d) = t_in on S ∖ {0},   t_out off S     (t_in = −4(k − 1), t_out = 4: −28 / +4 for k = 8, −12 / +4 for k = 4)
```

are the input of the bordered Goethals–Seidel array of Egan, Kotsireas, London and Ó Catháin (Theorem 7) of order 4(|G| + k − 1) once a border exists.

**Levels.** A level is a quotient G → G/N by an H-invariant subgroup N given by generators (`{"ker": [...], "name": ...}`; N need not be a coordinate
subgroup, so non-split G works: N = 2ℤ₄³ in ℤ₄³ × ℤ₁₃). The compression Y(a) = Σ_{π(g) = a} x(g) of an H-invariant sequence is constant on the orbits of
H in G/N, and Σᵢ PAF_{Yᵢ}(b) = Σ_{s ∈ b + N} Σᵢ PAF_xᵢ(s) on G/N: the same kind of problem on a smaller group, with explicit targets (CHECKED against independent
tuple arithmetic, `tests/test_nonsplit.py` A). Cosets of N are numbered by smallest element; translation and negation of G/N come from the coset table.

**Tree.** The root lists, for each block, every compressed vector of the first level with its multiplicity (a dynamic program over the H-orbits with
the coset sums fixed; exact) and joins the four lists exactly on the PAF classes of the quotient. Each solution tuple is a node. Inside a node the next level is
listed by the same DP restricted to the fibre of all earlier compressions and joined again. A **branch** is a node at depth `bdepth`, the restart unit
(one `.done` file per branch). A **leaf** (all levels fixed) enumerates the four fibres (meet in the middle over the quotient groups), optionally groups them through the
`split` levels, joins the four lists exactly on all PAF classes of G and verifies every hit over G. Blocks with equal coset sums are interchangeable
(`sym_swap`); an explicit normalizer of H (matrices that normalize H and preserve the target and the coset sums, checked at start-up) keeps one node per orbit,
with multiplicities. The ingredients are exact: DP (orbit and grouped strategies), four-list join (random linear key, hashed matches re-checked on the
full vector), bytewise deduplication, verification over G; a branch is `complete` unless a restriction of section 5 was active, a fibre exceeded `maxfibre`, or `--max-leaves` / `--branch-budget` stopped it.

## 2. Files

| File | Contents |
|---|---|
| `hcore.c`, `hlib.py` | C kernels (DP over orbits, join4); the model: `Group`, `QGroup`, `Level` (kernel or kept coordinates), `Problem` (S, k, targets, DP, join, fibre enumeration, exact `verify` over G). `libhcore.so` is built from `hcore.c` on first import |
| `hchain.py` | driver: `prepare`, `run`, `status`, `verify` |
| `modlattice.py` | Howell span over ℤ/p^k and the exact obstruction screen (section 6) |
| `configs/` | `n839_P1.json`, `n839_P2.json`, `n839_P3.json` (the three live row-sum patterns of 839), `normalizer_839.json`, `n3964_b.json` (the k = 4 control), `n1039_prod.json`, `n1039.json`, `n3964.json` (legacy k = 4 structures used by the tests) |
| `mkcfg2.py`, `mk_prod839.py` | builders of the 839 configs (`mk_prod839.py DIR` rewrites `n839_P*.json`, byte for byte as shipped) |
| `e2e/cert839.py` | an 839 solution record → exact check over G → border → certificate → `verify.build_checked` of the repository, in place |
| `e2e/cert_lib.py`, `e2e/border8.py`, `e2e/z7search.c` | certificate writer; k = 8 border from the exhaustive ℤ₇-invariant search `z7search` (4 s, deterministic; first border found) with the border identities checked exactly |
| `e2e/border4.py`, `e2e/e2e_k4_nonsplit.py` | k = 4 border search and the end-to-end control on the non-split G = ℤ₄² (n = 19, order 76) |
| `e2e/struct_k8.py`, `e2e/independent_check_3356.py` | structural k = 8 check with random sequences (border blocks of H Hᵀ), independent re-assembly of the matrix of order 3356 with an exact Gram check |
| `control_3964.py` | the k = 4 control: engine solution for 3964 → certificate format of `hadamard/3964` → H Hᵀ = 3964 I |
| `tests/` | validation (section 8) |

Output goes to `state/<name>/` (branch list `branches.npz`, `br/<id>.done`, `solutions.jsonl`, optional `leaves/`) and `results/` under the current directory, or to `--state`.
Nothing is written into the source tree except `libhcore.so` (`.gitignore`d) and Python caches.

## 3. Usage

```
cd search/bgschain; export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1; ulimit -v 3500000
python3 hchain.py prepare configs/CFG.json [--state DIR] [--force]       # tree down to bdepth -> branch list
python3 hchain.py run     configs/CFG.json [--state DIR] [--worker i/n | --dynamic] [--shuffle SEED] [--ids i,j] [--branches a:b]
                                           [--stop-first] [--budget SEC] [--max-branches N] [--leaf-mode count|search]
                                           [--max-leaves N] [--branch-budget SEC] [--sample N --seed S]
python3 hchain.py status  configs/CFG.json [--state DIR]
python3 hchain.py verify  configs/CFG.json state/NAME/solutions.jsonl
```

Restart: run the same command again (finished branches are skipped). After a crash remove stale claims of `--dynamic` workers:
`for c in state/NAME/br/*.claim; do [ -e "${c%.claim}.done" ] || rm "$c"; done`. `--leaf-mode count` replaces the leaf search by the exact count of leaves and fibre sizes.
Config keys: `mods`, `S_gens` (legacy `nq`), `gens`, `ys`, `t_in`, `t_out` (or an explicit `tgt`), `levels`, `split`, `bdepth`, `sym_swap`, `normalizer` (matrices or `"units"`), `maxfibre`,
`maxstates`, `dp_cache_states`, and the restrictions of section 5.

## 4. Order 3356 = 4 · 839 (CHECKED)

**Structure.** G = ℤ₄ × ℤ₄ × ℤ₄ × ℤ₁₃ (832 = 839 + 1 − 8), S = 2ℤ₄³ × ℤ₁₃ of index 8 (every involution of G lies in S), H = ⟨(B, 3)⟩ of order 21,
B = [[2,2,1],[1,0,0],[2,3,3]] mod 4. ℤ₄³ has 10 orbits of H and ℤ₁₃ has 5: 50 orbits per block, 25 PAF classes. The H-invariant subgroups are {0, 2ℤ₄³, ℤ₄³} × {0, ℤ₁₃}.
Pattern P1: coset sums a + dδ₀ with a = (2, 4, 4, 4), d = (4, −2, −2, 2), row sums (20, 30, 30, 34). `n839_P2.json` (a = (2, 4, 4, 4), d = (−4, −2, 2, 2)) and `n839_P3.json` (a = (0, 4, 4, 4), d = (2, −2, 2, 4)) are the other two patterns of the family that were kept.
Tower of `configs/n839_P1.json`, `bdepth` 2:

| Level | N | Measured |
|---|---|---|
| 1 | ℤ₄³ (quotient ℤ₁₃) | lists of 13,930 / 7,852 / 7,852 / 13,741 compressed vectors, exact join: 80 tuples, 20 up to the normalizer (order 96 modulo H) |
| 2 | 2ℤ₄³ (quotient ℤ₂³ × ℤ₁₃) | 350 nodes = branches |
| 3 (leaves) | ℤ₁₃ (quotient ℤ₄³) | exhaustive lists would have 1.07–2.1 · 10⁶ vectors per block: restricted, section 5 |

**Reproduction** of the first solution (branch 72 of 350, from a fresh state directory; one core, `nice`, loaded machine):

```
cd search/bgschain; export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1; ulimit -v 3500000
python3 hchain.py prepare configs/n839_P1.json --state state/n839_P1                        # 18 s: 350 branches
python3 hchain.py run configs/n839_P1.json --state state/n839_P1 --ids 72 --stop-first      # 35 s
python3 e2e/cert839.py state/n839_P1/solutions.jsonl --out results/certificate_3356.json    # 7 s (incl. z7search, 3.5 s)
```

Output (the `branch` line is the engine's; the solution record is deterministic: the tuples of a node are taken in a pseudo-random order seeded by the node data):

```
tree to depth 2: 350 branches  ... nodes_L1=1 sols_L1=20 nodes_L2=20 sols_L2=350 [18.3s]
branch 72/350: weight 10^35.6 wall 34.9s complete=False found=1 leaves=17 leaves_w=17 nodes_L3=1 sols_L3=2.11e+06
verified over G: H-invariance, coset sums [6, 2, 2, 2, 2, 2, 2, 2] ..., sum PAF = -28 on S\{0}, +4 off S;  pattern a [2, 4, 4, 4] d [4, -2, -2, 2] rowsums [20, 30, 30, 34]
REPO VERIFY OK (verify.build_checked): order 3356, H H^T = 3356 I, sha256 235a0f685b2f29daf207ac0012298376ba12b6c170431915c0d0483d942d45b5
```

* The branch list is the same in every run (350 branches ordered by weight). Branch 72 has 2,110,808 group tuples (level 3, 2000 vectors per block), of which 40 are taken; the search stops
  at the 17th of them: 17 fibre enumerations (984,362 sign vectors in all), 13 tuples rejected by the mod-8 screen, 4 joins of 4 · 10⁸ pair operations each, the hit in the fourth.
  Solution record (`masks` = the four sequences as orbit-sign bit masks): `[934370981655782, 452118242878979, 239244054975593, 305026013998834]`.
* `e2e/cert839.py` re-checks the record over G with independent code (H-invariance, coset sums, Σ PAF = −28 on S ∖ {0} and +4 elsewhere, all 832 shifts), takes the border for
  d = (4, −2, −2, 2) from `z7search` (the first border of its deterministic exhaustive search), writes the certificate and has `verify.build_checked` rebuild the matrix, check the entries, H Hᵀ = 3356 I and the
  SHA-256 of the canonical CSV. The certificate is the object of [hadamard/3356/certificate.json](../../hadamard/3356/certificate.json). `python3 e2e/independent_check_3356.py results/certificate_3356.json`
  assembles the matrix once more with separate code (float32 Gram, exact) and gives the same SHA-256.
* The border of 3356 is the one found by `z7search` for these d, so the SHA-256 is reproduced by the search alone, with no stored border.
* Branch 72 gave the first solution found (runs with `--dynamic --shuffle 7 --stop-first`, one worker per core, section 5); `--ids 72` runs just that branch.

## 5. Restrictions: the search is not exhaustive

The last level of 839 cannot be done exactly. Two options restrict the tree; a branch in which either one cut something is written `complete = false` in its `.done` file, and a run with them
is RESTRICTED: exact only inside the restriction, and a negative result proves nothing ("no solution in branch X" is not a statement the engine can make).

| Option | Effect |
|---|---|
| `"cap": [0, 0, 2000]` | per tree level, keep only the `cap` compressed vectors with the largest fibre in each block's list (deterministic) |
| `"child_cap": [0, 0, 40]` | per tree level, process only `child_cap` solution tuples of a node |
| `"child_order": "random"` | which tuples `child_cap` keeps: the first of a pseudo-random permutation seeded by the node data (restart-safe, the same in every run); default: the heaviest first |

`n839_P*.json` use all of them: 2,000 vectors per block and 40 random tuples per node on level 3 (about 2.1 · 10⁶ group tuples per node). The search is a lottery over branches and tuples:
`--dynamic --shuffle 7 --stop-first` with one worker per core claims the branches in a random order. Original production runs (3 workers; not repeated here): the first
solution after about 12 core-minutes of P1 search, five distinct solutions in about 90 core-minutes, 4 solutions in 4,561 core-seconds over 116 branches (one per 19 core-minutes).
For P2 and P3 the Gaussian estimate of the yield per join is 25–50 times lower (HEURISTIC; 80 joins, 0 solutions).

## 6. The exact mod-8 screen

`"screen_mods": [[2, 3]]` tests every leaf tuple before the expensive join. With V_b[i] the 25-class PAF vectors of the fibre of block b and T the target, a solution needs
T = Σ_b V_b[i_b], hence T − Σ_b V_b[0] lies in the ℤ/8-span of all differences V_b[i] − V_b[0]; `modlattice.py` computes that span (Howell echelon form over ℤ/p^k) and tests
membership. A tuple that fails has no solution: the screen is exact and cannot lose a solution (`tests/test_modlattice.py`: 125,560 membership tests against brute force,
planted sums never rejected; `tests/test_screen.py`: engine with the screen = brute force on planted instances, including 800 obstructed leaves). The obstruction is 2-adic:
for a ±1 sequence with w entries −1, PAF_x(g) = |G| − 4w + 4|Y ∩ (Y − g)| (Y = the −1 positions), so the condition constrains the parities of the autocorrelation counts.
It costs 0.15 s per tuple against 3.7 s for the join; in P1 about 91 % of 5,777 screened tuples were obstructed (production measurement), 13 of the 17 tuples of branch 72 (reproduction above).
`"score": true` adds a Gaussian local-limit log-score of the tuple (HEURISTIC ranking, `score_min` filters) and `"leaflog": true` a per-leaf log `state/<name>/leaves/<branch>.jsonl`.

## 7. Control: order 3964 = 4 · 991, k = 4 (CHECKED)

Structure of [hadamard/3964](../../hadamard/3964): G = ℤ₂ × ℤ₂ × ℤ₁₃ × ℤ₁₉ (ℤ₂₄₇ by the Chinese remainder theorem), S = 0 × 0 × ℤ₁₃ × ℤ₁₉ of index 4, H of order 27 with 46 orbits, levels ℤ₁₃ → ℤ₂² × ℤ₁₃ → ℤ₂² × ℤ₁₉, `split` ℤ₁₃ × ℤ₁₉
(`configs/n3964_b.json`, legacy coordinate form of the levels). No restriction: the count is exhaustive.

```
python3 hchain.py prepare configs/n3964_b.json --state state/n3964_b                          # 220 s (nice, loaded machine): 32,520 branches
python3 hchain.py run     configs/n3964_b.json --state state/n3964_b --leaf-mode count         # 114 s
python3 hchain.py status  configs/n3964_b.json --state state/n3964_b                          # 32520/32520 branches done, incomplete 0
#   sols_L2: 304   nodes_L3: 304   sols_L3 (= leaves): 124416
mkdir -p state/ctl_b12227/br; cp state/n3964_b/branches.npz state/ctl_b12227/
python3 hchain.py run configs/n3964_b.json --state state/ctl_b12227 --ids 12227 --stop-first   # 44 s: 5 leaves, 1 verified quadruple
python3 control_3964.py state/ctl_b12227/solutions.jsonl                                       # order 3964 H H^T == 3964 I: True
```

The 304 solutions at the second level (N = ℤ₁₉, quotient ℤ₂² × ℤ₁₃) are the 304 of step 1 in the README of [hadamard/3964](../../hadamard/3964), here reached through the extra first level ℤ₁₃
(32,520 branches); the 124,416 combinations at the third level (N = ℤ₁₃) are its step 2 summed over the 304 (160 for branch 12227). The quadruple found in branch 12227 is not the one of the
certificate (`same as certificate sequences: [False, False, False, False]`): it is a different solution, assembled with the border of `hadamard/3964/certificate.json` through
`verify.build_bordered_gs` (the border depends only on the coset sums, which the configuration fixes), and H Hᵀ = 3964 I is checked exactly.

## 8. Validation tests

Run from `search/bgschain` (`OMP_NUM_THREADS=1`, `ulimit -v 3500000`, `nice`, two cores shared with other jobs; every test prints `ALL OK` / `-> OK` lines, and `MISMATCH` or `FAILURES` on a failure). CHECKED:
the table gives the wall time of the last run and what was compared. "Brute force" is an independent enumeration of all sign vectors of the orbits and all quadruples with the right PAF sum
(`small_lib.py`, `nonsplit_lib.py`); the engine must return the same set of solutions (expanded under the normalizer and block exchanges) and the same orbit-weighted count, with every branch `complete`.

| Command | Compares | Result |
|---|---|---|
| `python3 tests/test_modlattice.py` | `ModSpan` (Howell echelon over ℤ/p^k) with brute-force span membership for all targets in (ℤ/q)^n, 7 (p, k, n); planted sums never rejected | 125,560 membership tests OK; 0 planted sums rejected (2 s) |
| `python3 tests/test_join.py` | `join4` with a numpy brute force on 60 random lattice-structured instances (negative values, constant columns, planted solutions), at two chunk budgets (50 and 100,000) | ALL OK (2 s) |
| `python3 tests/test_dp.py` | DP (orbit and grouped strategies), fibre enumeration, fibre counts and class values, on the k = 4 structures of 1039 and 3964 with planted vectors | ALL OK (4 s) |
| `python3 tests/test_screen.py` | engine with the exact screens (moduli 8, 32, 9, 5, 7 together) and the score switched on = brute force, 16 planted instances: ℤ₄² × ℤ₃ with a non-split normalizer element (5 seeds, with and without `split`), ℤ₄² × ℤ₃ (3 seeds), ℤ₄³ with k = 8 (3 seeds); the largest instance screens 2,096 leaf tuples and obstructs 800 of them | ALL OK, no solution lost (49 s) |
| `python3 tests/test_nonsplit.py A B` | A: level structure (coset labels, compressions, quotient classes, targets, orbit weights) and the identity Σ PAF_Y = Σ_{s ∈ b + N} Σ PAF on every quotient element, for ℤ₄², ℤ₄² × ℤ₃, ℤ₄³ × ℤ₃, ℤ₄³ × ℤ₅ and the real ℤ₄³ × ℤ₁₃ (839); B: DP lists = multiset of compressed vectors of all sign vectors of the fibre | 21 checks OK (4 s) |
| `python3 tests/test_nonsplit.py C1 C2 C3 C6` | engine = brute force. C1: G = ℤ₄², k = 4, n = 19, all 14 consistent patterns (8 with 352 solutions, 6 with none); C2, C3: ℤ₄² × ℤ₃, 28 patterns × several towers, all with 0 solutions (C3: 10 skipped, brute-force list too long); C6: non-split normalizer elements (1 + 2T, −1 on ℤ₃, both) with exchanged blocks, towers with kernel levels, splits, bdepth 1–2: 36 checks, up to 2,112 solutions | 126 checks OK, 10 skipped, 0 mismatches (305 s) |
| `python3 tests/test_nonsplit.py C4` | k = 8 planted instances: ℤ₄³ (3 seeds: 4,720 / 1,064 / 1,040 solutions) and ℤ₄³ × ℤ₃ (3 seeds with 290,040 / 738,576 / 5,192 solutions, up to four towers each; the towers whose join does not fit in 3.5 GB are skipped) | 11 checks OK, 4 skipped, 0 mismatches (2,637 s, nearly all of it the two large ℤ₄³ × ℤ₃ instances) |
| `python3 tests/test_small.py`, `test_small3.py`, `test_small4.py` | the k = 4 suite of the original engine on planted instances over ℤ₂² × ℤ₃ × ℤ₅ (towers, splits, bdepth 1–2, equal coset sums, `sym_swap` off, normalizer) | 48 / 32 / 63 checks OK (1,058 s / 12 s / 372 s) |
| `python3 tests/test_small2.py` | the same with the normalizer "units" of ℤ₃ × ℤ₅; only the instances planted invariant under it (the others are refused at start-up: "normalizer does not preserve the target") | 18 checks OK (225 s) |
| `python3 tests/test_planted1039.py` | a quadruple planted in the real structure of order 4156 = 4 · 1039 (k = 4, G = ℤ₂² × ℤ₇ × ℤ₃₇, `configs/n1039_prod.json`): the compressions on the first three levels are fixed, the engine runs the fourth level and the leaves | the planted quadruple is found, branch complete (253 s) |
| `python3 e2e/struct_k8.py P1` (P2, P3) | k = 8 chain on the real group of 839 with the repository's `build_bordered_gs` (autocorrelation test off, random sequences with the pattern's coset sums): the border blocks of H Hᵀ are exactly 3356 I and 0, and M Mᵀ = I₄ ⊗ (Σ PAF) | OK for all three patterns (7 s each) |
| `python3 e2e/e2e_k4_nonsplit.py` | k = 4 non-split G = ℤ₄² (n = 19): engine solution → border (`border4.py`) → certificate → `verify.build_checked` | order 76, H Hᵀ = 76 I, sha256 `543afdaa…` (0.4 s) |

The k = 8 border/certificate chain is thus checked structurally for all three patterns and end to end on the real solution (section 4); no small k = 8 non-split case with a solution is known in this family.

The tests of the original k = 4 engine (`test_small*.py`) were not changed except for `test_small2.py` (above); `tests/test_nonsplit.py` and `test_screen.py` use `tests/nonsplit_lib.py`.
