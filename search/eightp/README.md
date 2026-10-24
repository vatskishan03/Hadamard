# The orders 8p: Spence's construction over odd groups, the plug-in census, plain and block-symmetric structures

The programs behind the paragraph "The orders 8p" of [../NOTES.md](../NOTES.md) §92. Egan, Kotsireas, London and Ó Catháin list six open orders 8p up to 10 000 ([../../FRONTIER.md](../../FRONTIER.md), "Orders 8p"):
p = 599, 719, 739, 839, 1031, 1087, that is n = 2p = 1198, 1438, 1478, 1678, 2062, 2174 and the orders 4n = 4792, 5752, 5912, 6712, 8248, 8696. A Hadamard matrix of order 4n is a Goethals–Seidel array over a group of order n,
or Spence's bordered array over a group of order n − 1 = 2p − 1 (odd for every p), or a bordered array over a group of order n − t, or a recursive construction. All of them were searched or decided
as far as the tools below reach; **no matrix of any of the six orders was found, and none is claimed.** The labels are those of the repository: EXHAUSTIVE (a complete search of a stated class, validated against an
independent brute force), EXACT (a proof or an exact computation), VALIDATED, HEURISTIC, UNKNOWN.

## What was searched

| | structure | what decides it | label |
|---|---|---|---|
| (a) | real Goethals–Seidel over ℤ₂ₚ = ℤ₂ × ℤ_p with a multiplier group H ≤ ℤ_p\* | the ℤ₂-compression identity Σ(s₀² + s₁²) = 8p, Σ s₀s₁ = 0 with s_a ≡ ±1 (mod 2\|H\|) (`census/c8p_basic.py`, `check/dc_all.py`): every H except those listed below is dead; the live ones have at least 2¹⁶⁶ sign vectors per block | EXACT kills; unreachability HEURISTIC |
| (b) | negacyclic Goethals–Seidel of order 8p = quaternary sequences on ℤ_p with a multiplier group | Σ\|u₀ + \|H\|S\|² = 4p (S ∈ ℤ[i], \|a\| + \|b\| ≤ (p − 1)/\|H\|), the same live groups, at least 2⁹⁴ sign vectors per block | EXACT kills; unreachability HEURISTIC |
| (c) | Spence's bordered construction over an abelian group G of order v = 2p − 1 (odd), blocks invariant under a multiplier group | the odd-\|G\| construction ([below](#spences-construction-over-an-odd-group-and-its-certificate-path)); four kinds of structure: the Paley plug-ins G = K × F_r (`gs/`, `census/`), plain cyclic ℤ_v with every subgroup H ≤ ℤ_v\* (`plain/`), non-cyclic ℤ₃² × ℤ₁₃₃, and block-symmetric structures (A, B, μB, μ²B), (X, μX, μ²X, μ³X), (X, μX, Y, μY) (`spin/`) | EXACT filters and EXHAUSTIVE verdicts on the shapes listed in the reproduction block; the rest UNKNOWN |
| (d) | bordered Goethals–Seidel over \|G\| = n − t | the border t = 1 is (c); for t = 2, 4, 5, 6, 7 (Q = ℤ_t, all characters) the exact filters leave 54 shapes, all with t = 2, of which 52 have Gaussian estimates, all below 10⁻²⁵ except one (10⁷·⁵, p = 1031, ℤ₄ × ℤ₁₀₃, r = 5, 28 orbits, 76 sign bits); the other two (p = 1031, ℤ₄ × ℤ₁₀₃ and ℤ₂² × ℤ₁₀₃, r = 5, 8 orbits) have none (singular covariance) (`data/census/censC_*`, `estC8p.jsonl`) | EXACT filters, HEURISTIC estimates; not searched (UNKNOWN) |
| (e) | complex routes | the Craigen–Holzmann–Kharaghani decomposition 8(x + y), x + y = p, with complex Golay lengths: none under the survey's closure of the lengths; with the 'generous' closure under all products mn it would exist at five of the six p, but that closure is false (3 · 3 = 9: no complex Golay pair of length 9, `check/golay9.py`, 4⁹ sequences); the plain Szöllősi array over ℤ₂ₚ is a subset of (a) | EXACT (length 9: EXHAUSTIVE) |
| (f) | the recursive and direct theorems of Tables 3 and 4 of the survey (Paley I, II, blocks; Miyamoto I–III; Seberry–Yamada; Cooper–Wallis; Craigen–Holzmann–Kharaghani; T-sequences; Spence, Mathon, Yamada, Agaian) | `census/recursive8p.py`: no hypothesis is met at any of the six p | EXACT from the stated hypotheses |

Live groups in (a) and (b), from the judge's census of every subgroup order (`check/dc_all.py`; `census/c8p_basic.py` covers the odd orders): |H| = 2 at every p (300 to 544 quaternary orbits); |H| = 13 at p = 599 (quaternary only: 47 orbits, 4⁴⁷ = 2⁹⁴ sequences per block);
|H| = 3, 6, 9 at p = 739 (247, 124, 83 quaternary orbits); |H| = 5 at p = 1031 (207) and |H| = 10 (quaternary only: 104); |H| = 3 and 6 at p = 1087 (363, 182). Every other subgroup is dead (EXACT). The smallest live pool is 2⁹⁴ sequences per block (2¹⁶⁶ in the real array), far beyond any search (HEURISTIC).

## Spence's construction over an odd group, and its certificate path

Spence's bordered construction ([../NOTES.md](../NOTES.md) §3.1) is stated for even |G| (row sums (2, 0, 0, 0)). It also works for **odd |G| = v** (VALIDATED, `oddspence/oddspence_test.py`; rebuilt by the judge, `check/odd_spence_build.py`):
four ±1 sequences a, b, c, d on an abelian group G of order v with Σ PAF = −4 at every g ≠ 0 have row sums (r_a, r_b, r_c, r_d) with Σ r² = 4v − 4(v − 1) = 4, each odd (r ≡ v mod 2), so (±1, ±1, ±1, ±1). With
S = [[a, b, c, d], [−b, a, d, −c], [−c, −d, a, b], [−d, c, −b, a]] the 4 × 4 Hadamard matrix of the row sums, P = Q = S, K = −S and M the group-developed Goethals–Seidel array (X[g, h] = x(h − g), R the inversion),

    H = [ K        Q ⊗ 1ᵀ ]
        [ P ⊗ 1    M       ]

is a Hadamard matrix of order 4(v + 1) = 4n, n = 2p. H Hᵀ = nI was checked for cyclic G of order 3, 5, 7, 9, 11, 13 (orders 16 to 56), the solutions found by an exhaustive join. **Consequence** (EXACT, derivation): the row sums
(±1)⁴ allow balanced blocks, so for odd |G| there is no row-sum bound on the multiplier group, unlike (a) and (b); there is no character of order 2 to compress to, and the Parseval sum of a character ψ ≠ 1 is only Σ_X \|X̂(ψ)\|² = 4v + 4.

Certificate path. A hit of a sampler (`gs/gs_samp.py`) is verified at G level (exact integer PAF = −4 off 0, all row sums ±1) and passed to `gs/gs_cert.py`, which for odd \|G\| writes the certificate
`{order: 4v + 4, construction: spence, group, sequences}` (P and Q are placeholders). `verify/verify.py` of the repository accepts only the row sums (2, 0, 0, 0); the patch
[`oddspence/verify_odd_spence.patch`](oddspence/verify_odd_spence.patch) (apply it to a scratch copy of `verify/`, `patch -p1` from the copy's root) lets the `spence` builder accept (±1, ±1, ±1, ±1) and build P = Q = S, K = −S. It is not applied to the repository.
A hit would be checked by that patched `verify.py`, by the independent construction in `check/odd_spence_build.py` (`build`), and against Cati–Pasechnik's Table 4, [../../FRONTIER.md](../../FRONTIER.md) and the survey's Table 5. The three certificates of
`oddspence/` (orders 88, 112, 144 = 4 · 22, 4 · 28, 4 · 36: G = ℤ₇ × F₃, ℤ₉ × F₃, ℤ₅ × F₇, found by the sampler in seconds) pass the patched `verify.py` and are refused by the unpatched one.

## Programs

Subdirectories (each program runs in its directory; Python 3 with NumPy, SymPy, SciPy; gcc):

### `gs/`: the real Goethals–Seidel-type Paley plug-ins (the engine of §86) and the exact filters of §92

The shape is G = K × F_r, K abelian, r an odd prime power, a block is X = s ⊗ I + t ⊗ Q + e ⊗ J (Q the Paley matrix of F_r) with s, t, e functions on K constant on the orbits of a multiplier group H_K (twisted generators allowed); the target is
Σ PAF = 4|G|δ − (4 or 0 or Γ) off 0: `gs` (φ = 0, order 4|G|), `spence` (φ = −4, order 4|G| + 4, **odd or even \|G\|**), `bord` (φ = −Γ, order 4(|G| + t)). The K-level system, its defect identity (r + 1)e₊ + (r − 3)e₋ = 4|K| + φ(0) and the row-sum identities are those of §81 and §86.

| File | Contents |
|---|---|
| `gsk.py` | `GShape` (shape, targets, orbit structure from `pk.Kind`), the exact K-level check and the exact G-level `verify_g`, `enum_weight` (exact lists of the blocks of a defect weight), `row_patterns` (exact (σ, ε) row-sum patterns of an allocation, τ only by parity and range); **§92**: `block_classes3`, `row_patterns3` (the same with the exact t-sums τ) and `enum_pattern` (the exact list of one row-sum pattern) |
| `gl1.c`, `gl1.py` | level 1: the four lists (explicit or regenerated from seeds), class-restricted hash join (A, C) × (B, D) on the compressed autocorrelations, per-block exact constraints (defect weight, m-weight, row sum, (Σs, Σe)), canonical block sign; `libgl1.so` |
| `gl2.c`, `gl2io.py` | level 2: the t signs off the defects by a meet in the middle with the row-sum and PSD filters, any abelian K up to 640 points, four free blocks; builds `gl2_32`, `gl2_64`, `gl2_128` by the number of PSD channels |
| `gs_exh.py` | the exhaustive two-level solver: every allocation, exact lists, exact level 2; the count E_G = 16 × (ordered K-level count) is the number of ordered quadruples of ±1 sequences on G invariant under the plug-in group; allocations it cannot do are listed (`skipped`: join too big, too many level-1 instances, time) and make the verdict non-exhaustive |
| `gbf.py`, `gbrute.py` | the **independent** G-level brute force (the orbits of the multiplier group computed numerically, PAF over G in integer arithmetic) against which `gs_exh.py` is validated |
| `gs_samp.py`, `gs_worker.sh`, `alloc_scan.py` | the sampling production (one round per process, restart-safe loop; `GS_PATTERNS=1`: one exact row-sum pattern per round; **§92**: `GS_EXACT_CAP` (pattern-restricted exact lists, verdict EXHAUSTED per allocation and pattern) and `GS_PCACHE` (pattern cache)), automatic G-level verification and certification; the allocations that pass the identity and the row-sum patterns |
| `gs_cert.py` | certificates for `verify/verify.py` in a scratch copy of the repository (`repo_copy/`, see Build): `gs`, `spence` (even \|G\|: row sums normalised to (2, 0, 0, 0); **odd \|G\|**: P = Q = S, see above), `bord` |
| `pk.py`, `idx.py`, `est.py`, `est2.py`, `k4.py`, `est4.py` | shape code and the Gaussian/lattice estimator (HEURISTIC): `pk`, `idx`, `est` are identical to `../../plugh/{pk,idx,est}.py`; `k4.py` is the K-level system of the four-free-block plug-in for modes a, b, c (`Shape4`); `est4.py` the estimator per allocation (log₁₀ E, level-1 count, t bits) |
| `tau_filter.py`, `tau_patterns.py`, `tau_sweep.py`, `tau_check.py` | the **exact τ filter** (§92): the exact t-sum sets of one block (c, σ, ρ, τ) joined over four blocks on Σc = 4\|K\| + φ(0), Σ(σ² + rτ²) = 4\|G\| − φ(0), Σστ = 0 (r ≡ 1 mod 4), Σρ² = R², pairs bucketed by the exact c-sum; a shape with no 4-tuple has no solution (EXACT). `tau_patterns.py` lists the feasible row-pattern classes; `tau_sweep.py` applies the filter to the live shapes of the census (resumable); `tau_check.py` re-derives the recorded verdicts |
| `exh_sweep.py`, `exh_table.py` | the exhaustive verdicts of the live shapes (resumable jsonl), and the table of the recorded verdicts |
| `prod_tau.py`, `cells.py`, `small_alive.sh` | the sampling production for one τ-alive shape (allocation tokens from `tau_patterns`, exact-pattern lists when small), the size of its exact search per (allocation, pattern) cell, and the run on the eight smallest alive shapes (360 s each) |
| `val_odd.py`, `calib_odd.py`, `validate.py`, `run_validate.sh`, `run_validate2.sh` | validation of the engine on odd-\|G\| Spence shapes (`val_odd.py`: K-level engine against the G-level brute force; **new in §92**), calibration of `est4` on the same exact counts, and the battery of §86 (gs, spence, bord) |
| `t_eq.py`, `t_l1.py`, `t_l1b.py`, `t_l2.py`, `t_identity.py`, `t_pp.py`, `test_samp.py`, `test_samp2.py`, `run_tests.sh` | unit tests: K-level equations against G-level PAF, level 1 against a numpy pair join, level 2 against brute force over all t, the defect and row-sum identities on 300 solutions of each of six shapes, prime-power r, the sampler and its certificates; `run_tests.sh` is the quick regression of the stack |
| `t_bigN.py` … `t_bigN5.py` | large K (66 to 330 points, 3 to 5 orbits): the exhaustive engine against the G-level brute force, positive controls, and at \|K\| = 283 a pure-numpy level 2 (no `gl2.c`, no filters) on sampled level-1 instances; not part of the quick tests |

### `census/`: the plain multiplier census, the recursive theorems, the plug-in shape census and its estimates

| File | Contents |
|---|---|
| `c8p_basic.py` | the exact census of (a) and (b): multiplier groups of odd order h dividing (p − 1)/2, orbits 2 + 2(p − 1)/h (real) and 1 + (p − 1)/h (quaternary), the row-sum conditions; prints the DC-live and DC-dead groups (`../data/census/basic_census.txt`) |
| `recursive8p.py` | (f): every theorem of Tables 3 and 4 of the survey decided from its stated hypotheses, for the six p (including the complex Golay closure assumptions) |
| `census_a.py`, `shapes4.py`, `shapes4b.py` | the enumeration of the plug-in shapes (r an odd prime power dividing \|G\|, K abelian of order \|G\|/r, cyclic or product multiplier groups H_K with twists, 2 to 40 orbits) for modes a (unbordered), b (Spence), c (bordered with Q = ℤ_t or ℤ₂²) and the exact filters of `dp_filter` (defect identity and row sums in orbit units, by a DP on (c, ρ)); `run_census.sh` runs the six n |
| `run_est.py`, `check_est.py` | the Gaussian/lattice estimate (`gs/est4.py`, HEURISTIC) of every live shape (resumable), and the check of a re-run against the recorded estimates |
| `merge_cand.py` | the merged census files `cand8p`, `cand3only`, `cand3Ponly`, `candC8p`, `est_all8p` of `../data/census` from the per-n outputs (check against the recorded ones) |

### `plain/`: plain cyclic Spence over ℤ_v, v = 2p − 1, every multiplier subgroup H ≤ ℤ_v\*

| File | Contents |
|---|---|
| `plain_spence3.py` | four ±1 sequences on ℤ_v constant on the orbits of H, all row sums +1 (the gauge), Σ PAF = −4 off 0; the orbits of H, the classes of non-zero orbits under ±H, a quadratic form per class, the list L of sign vectors with row sum +1, and the exact count of ordered quadruples by a bucketed hash join of pair sums (every candidate verified exactly). `python3 plain_spence3.py v Omax [Omin [maxL]]` does every subgroup with Omin ≤ orbits ≤ Omax; statuses `noL` (no sign vector with row sum +1: EMPTY), `EMPTY`, `SOLUTION`, `bigL`. This is the **corrected** solver of §92 (see Deviations): the first version lost matches when it split a large bucket |
| `plain_nc3.py` | G = ℤ₃² × ℤ_m or ℤ₉ × ℤ_m (m prime to 3), automorphism group GL(2, 3) × ℤ_m\* (or ℤ₉\* × ℤ_m\*): random 2-generator subgroups, deduplicated by the orbit partition, each solved exactly (the sample, not all subgroups: HEURISTIC coverage) |
| `plain_census_L.py`, `est_plain.py` | the shapes with a +1 row-sum vector by a subset-sum DP; the Gaussian estimate (HEURISTIC) of a plain shape |
| `test_plain.py`, `cmp_logs.py` | validation of the solver against an independent dictionary join, with and without forced bucket splitting and forced chunking; comparison of logs ignoring the time field |
| `run_plain.sh`, `run_plain_hi.sh`, `run_plain_hi2.sh`, `tally.py` | the queues (O ≤ 22 for the six v; O = 23 and 24 for v = 1197, 1677, 2061, 2173 and the re-verification of 1477 and of the ℤ₃² × ℤ₁₃₃ sample; 2061 and 2173 with the final solver), and the table of all logs of `../data` |

### `spin/`: block-symmetric Spence structures over ℤ_v

A is invariant under K = ⟨H, μ⟩, B under H, [K : H] = d. Blocks (A, B, μB, μ²B) (d = 3: PAF_A(g) + PAF_B(g) + PAF_B(μg) + PAF_B(μ²g) = −4, `spin3.c`, `spin_drv.py`), (X, μX, μ²X, μ³X) (d = 4) and (X, μX, Y, μY) (d = 2) (`spinD.c`, `spinD_drv.py`),
all row sums +1; Gray-code enumeration of B (or X) in C, a hash lookup of A (d = 3), the self-join of the +1 list (d = 2). The **COUNT** of `spin3` is the number of pairs (B, distinct PAF vector of an admissible A): admissible A with equal PAF vectors count once, so it is not the number of ordered (A, B) (zero and non-zero are exact);
the d = 2 counts are 64-bit-hash matches, every printed solution is verified exactly. `run_spin.sh`, `run_spin_hi.sh` (d = 3), `run_spinD.sh` (d = 2, 4) run all pairs (H, μ) with O_H in a range for the v where the subgroup of index d exists.

### `oddspence/`, `check/`, `data/`

| File | Contents |
|---|---|
| `oddspence/oddspence_test.py` | the odd-\|G\| construction for ℤ_v, v = 3, 5, 7, 9, 11, 13: sequences by an exhaustive join, H Hᵀ = nI |
| `oddspence/verify_odd_spence.patch`, `cert_88.json`, `cert_112.json`, `cert_144.json` | the patch of `verify/verify.py` and three certificates (sampler hits at \|G\| = 21, 27, 35) |
| `check/dc_all.py`, `golay9.py`, `odd_spence_build.py`, `gcount.py`, `plain_shapes.py` | independent checkers (code of the judge of §92, sharing nothing with the engines): the DC census for all subgroup orders; the complex Golay pairs of length 2 to 9 by exhaustive search (none of length 7 and 9); the odd-\|G\| construction for ℤ_v, v = 3 to 11; a chunked G-level count of Spence quadruples over ∏ ℤ_d with a multiplier group (re-counts the shapes of p = 719 and 739 and two controls); the census of plain cyclic shapes (subgroups, orbit counts, +1 row-sum vectors) |
| `data/census/` | `cens*_<n>.jsonl` (outputs of `census_a.py`, 30 files), the merged `cand*.jsonl`, the estimates `est*.jsonl` (`est8p`, `est3only`, `est3Ponly`, `estC8p`, `est_all8p`), `basic_census.txt` |
| `data/gs/` | `val_odd.jsonl` (166 odd-\|G\| Spence shapes: 159 MATCH, 7 partial), `VALIDATION.txt` (battery of §86), `tau_sweep2.jsonl` (the sweep of the exact τ filter: 213 shapes), `tau_sweep_v1.jsonl` (an earlier, weaker run of the filter, kept as the origin of `cand_alive.jsonl`; superseded), `cand_alive.jsonl` (the 139 shapes alive in that earlier run; the input of the `exh_alive` run), `exh_small.jsonl`, `exh_r3.jsonl`, `exh_alive.jsonl` (+ `.log`: the records of `exh_sweep.py`), `alive_ranked.json` (the 99 τ-alive shapes with the modelled cost of the sampler; the script that ranked them was not kept), `calib_odd.json`, `alloc_1478_o11*.txt` (`alloc_scan.py` on ℤ₂₁₁, r = 7, H_K = ⟨34⟩, 11 orbits), `runs/sa_1198_*` (the small-alive production), `runs/p739a`, `runs/pilot*` (sampler pilots at p = 739) |
| `data/plain/`, `data/spin/` | the logs of every verdict: `v<v>.O22.log` (all subgroups with at most 22 orbits), `r3_*.log` (corrected solver: 1477 up to 22, 1197 and 1677 at 23 to 24 orbits, 2061, 2173, ℤ₃² × ℤ₁₃₃), `cal*.log` (positive controls v = 105 to 255), `D*.log` and `v*.log` (spin d = 2, 4 and d = 3: `O26` = O_H ≤ 26, `O27-30` = 27 to 30) |

## Build

```
cd search/eightp/gs && gcc -O3 -march=haswell -shared -fPIC -o libgl1.so gl1.c -lm && for ps in 32 64 128; do gcc -O3 -march=haswell -DPS=$ps -o gl2_$ps gl2.c -lm; done && mkdir -p runs tmp out work
cd search/eightp/spin && gcc -O3 -march=haswell -o spin3 spin3.c && gcc -O3 -march=haswell -o spinD spinD.c && mkdir -p runs
cd search/eightp/gs && mkdir -p repo_copy && cp -r ../../../verify repo_copy/verify && ln -sfn ../../.. repo_copy/search && (cd repo_copy && patch -p1 < ../../oddspence/verify_odd_spence.patch)     # the scratch repository of gs_cert.py
```

Rules for every job: `nice -n 5`, `OMP_NUM_THREADS=1`, `ulimit -v 3500000`; detach long jobs with `setsid nohup`, record the PID, kill only recorded PIDs.

## Reproduction of search/NOTES.md §92 (the orders 8p)

Run in the directory named, after the builds. The comment of a line is the output expected, trimmed. Times are wall seconds on one core of a machine whose four cores were all busy (CPU times are about 0.7 of them; shorter on a free core).

**(a), (b), (f), (e): the plain census, the theorems, the complex Golay pairs**

```
cd census
python3 c8p_basic.py | diff - ../data/census/basic_census.txt    # 6 s  identical: for each p the DC-live and DC-dead groups of odd order (p = 599: h = 13 live as quaternary, 23 and 299 dead; 739: h = 3, 9 live; 1031: h = 5; 1087: h = 3)
python3 recursive8p.py                                            # 0.4 s  every theorem 'fails' at the six p; the only 'AVAILABLE' lines (p = 599, 719, 739, 839, 1087) are 'same with the generous closure under mn as well': Craigen-Holzmann-Kharaghani 8(x + y) with the complex Golay lengths closed under plain products, a closure that is false (9 = 3 * 3)
cd ../check
python3 dc_all.py                                                 # 8 s  the same census for all subgroup orders, odd and even (|H| = 2 live at every p; |H| = 6 at 739 and 1087, 10 at 1031 quaternary live)
python3 golay9.py                                                 # 2 s  n count: 2 64, 3 128, 4 512, 5 512, 6 2048, 7 0, 8 6656, 9 0 (ordered pairs of quaternary sequences): no complex Golay pair of length 7 or 9
```

**The census of the plug-in shapes (c), (d)** (EXACT filters; 30 output files, byte for byte the recorded ones; HEURISTIC estimates):

```
cd census && bash run_census.sh          # 102 s (65 s CPU)  cens_<n>, censP_<n> (r >= 5, cyclic and product H_K), cens3_<n>, cens3P_<n> (r >= 3), censC_<n> (bordered) for n = 1198 ... 2174;  'census: 30 files identical to ../data/census, 0 different'
python3 merge_cand.py                    # 0.1 s  cand8p 489 rows (229 live), cand3only 111, cand3Ponly 152, candC8p 164, est_all8p 420: all identical to the recorded files
python3 check_est.py 1478                # 15 s   the Gaussian/lattice estimates of the 7 live shapes of p = 739 (ℤ₂₁₁ × F₇): '7 of 7 identical (logE of the best shapes: [6.76, 9.97, 22.79])'
```

Live shapes after the filters (Spence mode, r ≥ 3, cyclic and product H_K, 2 to 40 orbits; from `data/census`): p = 599: 337 of 464, 719: 1 of 3, 739: 7 of 25, 839: 49 of 169, 1031: 26 of 69, 1087: 0 of 22 (every shape dead by the defect identity: no allocation of orbit-unit defects exists).
Bordered shapes (mode c): 21, 0, 21, 0, 12, 0 live, all with Q = ℤ₂.

**The exact τ filter and the exhaustive plug-in verdicts** (EXACT / EXHAUSTIVE; the engine is validated first):

```
cd gs
bash run_tests.sh                                       # 99 s on a less loaded run, 141 s in the smoke test  t_eq: 'mismatches 0' for six shapes; t_l1 ALLOK; t_l2 'mismatch 0'; the exhaustive pipeline equals the G-level count (gs Z3 x F3: 286,685,568; spence Z2 x F3: 124,416); the sampler certificates verify (orders 84, 132, 52)
python3 t_pp.py                                         # 1 s   prime-power r: 7 shapes MATCH (248,832; 1,920; 768 nonzero)
python3 val_odd.py runs/val_odd_test.jsonl 60           # the odd-|G| battery: K in Z3 ... Z15, Z3^2, Z5 x Z3; r = 3, 5, 7, 9; twisted and untwisted H_K: ../data/gs/val_odd.jsonl is the full record (159 MATCH, 35 nonzero, largest 741,520,000, 0 MISMATCH, 7 partial); the first 19 rows in 150 s equal it (18 MATCH identical, one time-limited 'partial' row differs in its K count)
python3 tau_check.py 12                                 # 66 s  97 shapes with 3 to 12 orbits re-derived: 27 alive, 70 dead, 0 mismatches with ../data/gs/tau_sweep2.jsonl
python3 tau_filter.py '[211]' 7 '[[[4],1]]'             # 0 s   (False, {'blockstates': 8, 'features': 2}): p = 739, K = 211, H = <4>, dead
python3 gs_exh.py '[7,19]' 9 spence '[[[6,1],1],[[1,2],-1]]' | tail -1      # 4 s   'allocs': 5 ... 'E_G': 0 ... 'exhaustive': True   (p = 599, 8 orbits)
python3 gs_exh.py '[3,3,7]' 19 spence '[[[2,1,1],1],[[1,2,1],1],[[1,1,3],1]]' | tail -1    # 5 s  'allocs': 140, E_G 0, exhaustive   (p = 599, 8 orbits; the one of the 14 that is tau-alive)
python3 gs_exh.py '[7,19]' 9 spence '[[[2,1],1],[[1,2],1]]' | tail -1       # 4 s   'allocs': 833, E_G 0, exhaustive   (p = 599, 6 orbits)
python3 gs_exh.py '[7,19]' 9 spence '[[[6,1],-1],[[1,2],1]]' | tail -1      # 4 s   E_G 0, exhaustive   (p = 599, 8 orbits)
python3 gs_exh.py '[3,3,19]' 7 spence '[[[2,1,1],1],[[1,2,1],1],[[1,1,2],-1]]' | tail -1   # 53 s  E_G 0, exhaustive   (p = 599, 8 orbits, sign twist)
python3 exh_sweep.py runs/x_1478.jsonl 2 8 300 1478     # 21 s  K = 211, r = 7: 3 orbits (0.0 s) and 7 orbits (20.6 s; recorded 12.4 s): E_G 0, exh True  (p = 739)
python3 exh_sweep.py runs/x_2062.jsonl 2 8 300 2062     # 1.4 s K = 229, r = 9: 5 and 7 orbits: E_G 0, exh True  (p = 1031)
CAND=../data/census/cand3only.jsonl EST=../data/census/est3only.jsonl python3 exh_sweep.py runs/x_1438.jsonl 2 8 300 1438    # 6.5 s  K = 479, r = 3, H = <2>, 3 orbits: E_G 0, exh True  (p = 719; the independent G-level count of check/gcount.py agrees)
python3 exh_sweep.py runs/x_1198.jsonl 2 6 1500 1198    # 179 s (the fourth argument is a wall-clock limit per shape; with 300 s a loaded machine can stop the 126 s shape early: exh False)  p = 599, the four shapes with at most 6 orbits (K = 7 x 19, r = 9: three shapes; K = 9 x 19, r = 7: 126 s): E_G 0, exh True
python3 exh_table.py                                    # the table of the recorded verdicts: 14 distinct shapes EXHAUSTIVE-empty (p = 599: 9, 719: 1, 739: 2, 1031: 2), 12 not finished
```

The 14 EXHAUSTIVE-empty shapes (`../data/gs/exh_*.jsonl`; K, r, H_K, orbits): p = 599: ℤ₇ × ℤ₁₉, r = 9 (**six** shapes: three with 6 orbits, three with 8), ℤ₉ × ℤ₁₉, r = 7 (6 orbits), ℤ₃ × ℤ₃ × ℤ₇, r = 19 (8), ℤ₃ × ℤ₃ × ℤ₁₉, r = 7 with a sign twist (8);
p = 719: ℤ₄₇₉, r = 3, H_K = ⟨2⟩ (3 orbits); p = 739: ℤ₂₁₁, r = 7, H_K = ⟨4⟩, ⟨5⟩ (3, 7); p = 1031: ℤ₂₂₉, r = 9, H_K = ⟨3⟩, ⟨4⟩ (5, 7). Thirteen of the fourteen are also killed by the τ filter; the fourteenth is τ-alive and empty.
Not finished (UNKNOWN): at p = 599 three shapes (ℤ₇ × ℤ₁₉ with 8 and 9 orbits, ℤ₃ × ℤ₃ × ℤ₁₉, r = 7, 8 orbits: 300 s limit, E_G = 0 so far); at p = 839 three shapes (ℤ₁₃ × ℤ₄₃, r = 3, 9 and 12 orbits: MemoryError of the list generator); at p = 1031 six shapes (ℤ₃ × ℤ₂₂₉, r = 3, `N <= 640` of `gl1.py`).
The τ sweep (`../data/gs/tau_sweep2.jsonl`): every live shape with at most 20 orbits and one of 21 orbits, 213 in all: p = 599: 165 shapes, 87 alive (8 to 21 orbits), 75 dead, 3 undecided (DP memory, 18 orbits); p = 719: 1, dead; 739: 4, 1 alive (ℤ₂₁₁, r = 7, |H| = 15, 15 orbits); 839: 23, 11 alive (9 to 20 orbits); 1031: 20, all dead;
1087: none live. The sweep stopped at 21 orbits: 207 of the 208 live shapes with 21 or more orbits are not swept (UNKNOWN).

**Plain cyclic Spence over ℤ_v** (EXHAUSTIVE per shape, corrected solver; the six v are 1197, 1437, 1477, 1677, 2061, 2173):

```
cd plain
python3 test_plain.py                                   # 4 s   104 shapes (v = 9..75, <= 11 orbits), 52 with a nonzero count, 103 with the independent dictionary join: 0 mismatches between the solver, the forced-splitting solver, the forced-chunking solver and the join
mkdir -p runs; for v in 105 143 189 205 231 255; do python3 plain_spence3.py $v 18 2 1e5 > runs/cal$v.log; done; python3 cmp_logs.py ../data/plain runs/cal*.log    # 6 of 6 identical: the positive controls (v = 143: counts 16, 24, 16, 400, 40)
mkdir -p runs && OMAX=22 bash run_plain.sh              # 140 s  all subgroups with 2 to 22 orbits for the six v
python3 cmp_logs.py ../data/plain runs/v1197.O22.log runs/v1437.O22.log runs/v1477.O22.log runs/v1677.O22.log runs/v2061.O22.log runs/v2173.O22.log    # 6 of 6 identical (statuses, counts; time field ignored)
python3 plain_spence3.py 1197 24 23 1e5 > runs/r3_1197_23_24.log    # 362 s wall (218 s of solver time in the recorded run): 20 shapes, 12 EMPTY, 8 noL
python3 plain_spence3.py 1677 24 23 1e5 > runs/r3_1677_23_24.log    # 471 s wall (412 s): 7 shapes, 3 EMPTY, 4 noL
python3 cmp_logs.py ../data/plain runs/r3_1197_23_24.log runs/r3_1677_23_24.log     # 2 of 2 identical to the recorded logs of the corrected solver
python3 plain_nc3.py z3z3 133 3000 16 1e5 > runs/nc_small.log; python3 cmp_logs.py --subset ../data/plain/r3_nc_z3z3_133.log runs/nc_small.log   # 7 s  88 partitions, 28 solved, all EMPTY, every record in the recorded log
python3 tally.py                                        # the table of every recorded log
python3 plain_census_L.py 2061 23 30                    # 1 s  the shapes with a +1 row-sum vector beyond 22 orbits (the same for the other v): 23 (three, L = 291,720), 25 (two), 27 orbits
```

Shapes with at most 22 orbits: 55, 5, 33, 45, 26, 22 for v = 1197, 1437, 1477, 1677, 2061, 2173, of which 12, 2, 6, 2, 2, 4 have a +1 row-sum vector (the others are `noL`, EMPTY exactly); every one has count 0 (the judge re-counted all of them with a different join and its own census of the subgroups: `check/plain_shapes.py`).
O = 23 and 24: v = 1197: 20 shapes (12 EMPTY, 8 noL; 218 s of solver time), v = 1677: 7 shapes (3 EMPTY, 4 noL; 412 s), v = 1477: 3 shapes (all noL). v = 1437 is complete: no subgroup of ℤ₁₄₃₇\* has 10 to 479 orbits (the orbit counts are 4, 5, 6, 9, 480, 719, 720, 958, 1437).
Open (UNKNOWN): v = 2061 has 5 shapes with 23 or 24 orbits: two are `noL`, the three with |H| = 114 and 23 orbits have a +1 list of 291,720 vectors each, and their exact join (stopped after 3 h 13 min of CPU and more than 2.2 GB with no shape finished; the first attempts at these v ended in MemoryError at the 3.5 GB limit) did not finish; `solve` regenerates the pair list for every residue part, so its time grows like the square of the number of pairs, and the 240 s of the 41,760-vector shape of v = 1677 extrapolate to days here (HEURISTIC); v = 2173 has one (23 orbits, 154,000 vectors; MemoryError of the first solver, not run with the final one).
Beyond those, `python3 plain_census_L.py v 23 30` lists the shapes with a +1 vector that nobody has solved: v = 1197: 26 shapes with 25 to 30 orbits; 1477: one with 28 orbits (4.5 · 10⁶ vectors); 1677: four with 27 to 30 orbits (up to 3.5 · 10⁶); 2061: two with 25 and one with 27 orbits (up to 6.1 · 10⁶); 2173: shapes with 25, 27 and 30 orbits (up to 3.7 · 10⁶).
Non-cyclic ℤ₃² × ℤ₁₃₃ (v = 1197): 40,000 random 2-generator subgroups of GL(2, 3) × ℤ₁₃₃\*, 583 distinct partitions with 3 to 24 orbits, 283 with a +1 row-sum vector, all EMPTY (3076 s of solver time): HEURISTIC coverage (a sample), EXHAUSTIVE per shape. ℤ₃² × ℤ₂₂₉ (v = 2061) was not run.

**Block-symmetric structures** (EXHAUSTIVE for the stated range, VALIDATED by the self-test):

```
cd spin
python3 spinD_drv.py selftest                           # 3 s   'selftest: pairs 126 mismatches 0 with nonzero count 46' (d = 2, 4 against a brute force on v = 9..63; the printed solutions verify)
python3 spin_drv.py 1477 26 2                           # 6 s   v = 1477: 31 pairs (H, mu) with 2 <= O_H <= 26, 4 live, count 0
python3 spinD_drv.py 1477 2 26                          # 6 s   v = 1477, d = 2: 36 pairs, 16 live, count 0
python3 spinD_drv.py 1677 4 26                          # 6 s   v = 1677, d = 4: 14 pairs, 2 live, count 0
python3 spin_drv.py 1197 26 2                           # 73 s  v = 1197, d = 3: 73 pairs, 20 live (max O_H = 26), count 0
python3 spinD_drv.py 1197 2 26                          # 113 s v = 1197, d = 2: 79 pairs, 34 live (max O_H = 26), count 0
cd ../plain && python3 tally.py ../data                 # the table of all recorded logs
```

Recorded (from the logs of `../data/spin`, not re-run beyond the commands above): d = 3, O_H ≤ 26 and 27 to 30: v = 1197 (73 + 46 pairs, 20 + 14 live, max 30, 76 s + 2970 s), 1477 (31 + 10, 4 + 1, max 28), 1677 (34 + 15, none live), 2061 (28 + 1, 2 + 1, max 27); v = 1437 and 2173 have no subgroup of index 3 in ℤ_v\* (3 does not divide 956 and 2080).
d = 2, O_H ≤ 30: v = 1197 (79 + 27 pairs, 34 + 24 live, max 30), 1437 (6, 4 live), 1477 (36 + 3, 16 + 1), 1677 (78 + 48, 9 + 12), 2061 (35 + 3, 15 + 3), 2173 (34 + 7, 18 + 4). d = 4: v = 1197, 1437, 1477 no admissible (H, μ) with O_H ≤ 30; 1677 (14 + 7, 2 + 3 live), 2061 (10, 4), 2173 (24 + 6, 14 + 3). Count 0 in every case.

**The odd-|G| construction and the certificate path** (EXACT / VALIDATED):

```
cd oddspence && python3 oddspence_test.py               # 2 s  v 3 order 16 ... v 13 order 56: HHT=nI True, rowsums [-1, -1, -1, -1]
cd ../check && python3 odd_spence_build.py 2>/dev/null  # 0.1 s  v = 3, 5, 7, 9, 11: orders 16, 24, 32, 40, 48, entries +-1 True, H H^T = n I True, rowsums [1, 1, 1, 1]  (the judge's code; RuntimeWarnings of the 64-bit hash are harmless)
cd ../gs && python3 gs_samp.py spence '[7]' 3 '[]' null "1.0:1.0:2.0:2.0,0.0:2.0:2.0:2.0,1.0:1.0:1.0:3.0" 300 3000 60 1 runs/oddtest_7_3
   # 10 s  a hit in the first round (G = Z_7 x F_3, |G| = 21): '*** certificate {order: 88 ... verify_ok: True ... hadamard/88 ok H H^T = 88 I sha256 ...}' written to repo_copy/hadamard/88 and out/spence_88_oddtest_7_3.json
cd repo_copy && for o in 88 112 144; do mkdir -p hadamard/rec$o; cp ../../oddspence/cert_$o.json hadamard/rec$o/certificate.json; python3 verify/verify.py hadamard/rec$o; done    # 'ok  H H^T = 88 I  sha256 dd090020...', 112 (d0e3198b...), 144 (31a5a620...)
```

The unpatched `verify/verify.py` answers `FAIL  row sums must be (2, 0, 0, 0)` on the same certificates.

**Judge's independent checkers**

```
cd check
python3 gcount.py                                       # 0.6 s  A  p=719 K=[479] r=3 <2>: (9, 9, 0, 239); B  p=739 K=[211] r=7 <4>: (9, 9, 0, 315): count 0; controls Z3 x F3 trivial (9, 10, 6259194, 1) and 16 * 6259194 = 100,147,104, Z3 x F7 (9, 9, 62208, 3)
python3 plain_shapes.py 22 1197 1437 1477 1677 2061 2173    # 33 s  the shape censuses 55, 5, 33, 45, 26, 22 and the shapes with a +1 row-sum vector 12, 2, 6, 2, 2, 4
```

## Recorded runs that were not re-run (labelled 'from the recorded log')

* The full validation battery of §86 (`VALIDATION.txt`: gs 30 MATCH shapes (20 nonzero), spence 47 (9), bord 30 (23); battery 2 with the final binaries 22 + 32 + 22 MATCH; 0 MISMATCH), `run_validate.sh`, `run_validate2.sh`; the full `val_odd.py` (166 shapes, hours), and `calib_odd.py` (`est4` on the 35 nonzero and 124 empty odd-\|G\| shapes: estimate − exact between −315 and +2.8 decades on the nonzero ones, 43 of the 124 empty shapes get log₁₀ E > 0, at most 7.2: HEURISTIC with false alarms up to 10⁷).
* `t_identity.py` (the defect and row-sum identities on 300 solutions of each of six spence and bordered shapes: 0 violations; 251 s, in the long smoke group), the 8-orbit exhaustive verdict ℤ₇ × ℤ₁₉, r = 9, H = ⟨(3, 1), (1, 8)⟩ (246 s in the recorded run), the full τ sweep (588 CPU s in all; `tau_check.py 12` re-derives the shapes up to 12 orbits).
* The sampling production (HEURISTIC; `prod_tau.py`, `small_alive.sh`, `data/gs/runs/`): on the eight smallest alive shapes, 360 s each (seven of p = 599 completed; the eighth, p = 839, ℤ₁₃ × ℤ₄₃, r = 3, produced no output in its 420 s): sa_1198_0 (ℤ₃ × ℤ₃ × ℤ₇, r = 19) 372 of 374 allocation tokens EXHAUSTED, 0 instances; the first production attempt at p = 739 (`p739a`) stalled on huge pattern lists and was killed. No round of any sampler produced a complete instance at p = 599, 739, 839: the sampled coverage of the alive space is about 0 (UNKNOWN, not evidence of absence).
  The modelled cost of the best alive shapes (`alive_ranked.json`, HEURISTIC): p = 599, ℤ₉ × ℤ₁₉, r = 7, |H| = 18, 13 orbits: 10⁵·⁸ s on one core; p = 739, ℤ₂₁₁, r = 7, |H| = 15, 15 orbits: 10⁶·¹ s.
* The ℤ₃² × ℤ₁₃₃ sample (3076 s of solver time: the recorded log `data/plain/r3_nc_z3z3_133.log`; a 3000-sample subset is re-run above), and the spin runs of v = 1197 with 27 ≤ O_H ≤ 30 (2970 s and 4822 s) and the other long spin logs of `data/spin`.


## Relaunch of the long jobs (every queue is resumable; when this was written only `plain/run_plain_hi2.sh` of the original runs was still going)

```
cd census
mkdir -p runs; python3 run_est.py ../data/census/cand8p.jsonl runs/est8p.jsonl 700 150 6          # the Gaussian/lattice estimates of the 229 live shapes (the same for cand3only, cand3Ponly, candC8p; resumable by key)
cd ../gs
python3 tau_sweep.py runs/tau_sweep2.jsonl ../data/census/cand8p.jsonl ../data/census/cand3only.jsonl ../data/census/cand3Ponly.jsonl     # the exact tau filter, ascending orbit counts, 150 s per shape (588 CPU s up to 20 orbits); kill the recorded PID beyond 21 orbits
(CAND=../data/gs/cand_alive.jsonl EST=../data/census/est_all8p.jsonl setsid nohup python3 exh_sweep.py runs/exh_alive.jsonl 2 13 300 > runs/exh_alive.log 2>&1 & echo $! > runs/exh_alive.pid)   # exhaustive verdicts of the shapes with 2 to 13 orbits, 300 s each
setsid nohup bash small_alive.sh > /dev/null 2>&1 &                                  # the exact-pattern production on the eight smallest alive shapes (about 6 minutes each; sa_1198_* in ../data/gs/runs; the eighth, at p = 839, gave no output)
setsid nohup python3 prod_tau.py spence '[9,19]' 7 '[[[2,4],1]]' 3600 SEED runs/p599a 2000 200000 > runs/p599a.log 2>&1 &    # sampling production for one tau-alive shape: one JSON line per round in runs/p599a.jsonl; a verified hit writes runs/p599a.HIT_<time>.json and certifies (use a new SEED per relaunch)
cd ../plain
OMAX=22 setsid nohup bash run_plain.sh > /dev/null 2>&1 &                            # all subgroups up to 22 orbits for the six v (143 s)
setsid nohup bash run_plain_hi.sh > /dev/null 2>&1 &                                 # O = 23, 24 for v = 1197, 1677 (and 2061, 2173), the 1477 re-verification, the Z_3^2 x Z_133 sample: hours; .done files in runs/
setsid nohup bash run_plain_hi2.sh > /dev/null 2>&1 &                                # v = 2061 and 2173 at O = 23, 24 with the final solver: the queue that was still running (see below)
python3 plain_nc3.py z3z3 229 40000 24 1e5 > runs/nc_z3z3_229.log                    # the sample of Z_3^2 x Z_229 (v = 2061): never run
cd ../spin
OMAX=26 bash run_spin.sh; OMIN=2 OMAX=26 bash run_spinD.sh; OMIN=27 OMAX=30 bash run_spin_hi.sh; OMIN=27 OMAX=30 bash run_spinD.sh     # all block-symmetric structures up to O_H = 30 (about 3 hours in all, dominated by v = 1197)
```

Hit handling: the independent numpy check H Hᵀ = nI of the final matrix (G-level exact integer PAF = −4 off 0, row sums ±1), the patched `verify.py` (above), and the order against Cati–Pasechnik's Table 4, [../../FRONTIER.md](../../FRONTIER.md) and the survey's Table 5.

## Deviations from the sources

The code of the original searches is carried over unchanged except: docstrings and headers (NOTES sections instead of rounds and lenses), paths relative to the directory (`census/` and `plain/` scripts import the shared modules from `../gs` and `plain_spence3`; the census data in `../data/census`), the spin drivers import `plain_spence3` instead of
`plain_spence2` (they use only the subgroup and orbit helpers, which are identical), the run scripts write to `runs/`. New: `gs/{exh_table,tau_check}.py`, `census/{merge_cand,check_est}.py`, `census/run_census.sh`, `plain/{test_plain,cmp_logs}.py`, `oddspence/verify_odd_spence.patch`, this README. Not carried over: the buggy `plain_spence2.py` and `plain_nc.py`, the scratch variants, `gs_exh2.py`
(`gs_exh.py` with an unused argument), the r31 benchmarks and throughput scripts, `rank4.py`, `charfilter.py`, `exhaust4.py` (the round-30 level 2), the 700-point build of the engine (used for no verdict).
