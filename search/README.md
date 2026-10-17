# Search programs

The programs that found the matrices in this repository, and the experiments recorded in [NOTES.md](NOTES.md) and
[../NOTES.md](../NOTES.md). This is research code: each program was written for one family of searches and is
documented by its docstring and by the notes. No certificate depends on it; `verify/` checks every matrix on its own.

**Requirements.** Python 3 with NumPy; JAX for `generate.py`; a C compiler with OpenMP. A few experiments use
OR-Tools or PySAT. The C helpers are compiled on first use by `workdir.py`. Hit logs, temporary list files and
binaries go to `$HADAMARD_WORK` (default `<tmp>/hadamard-work`); only the `*_certify.py` scripts write into the
repository.

**Commit hashes.** The order READMEs give the commit at which each search was run. Until commit `b7bfde3` this
directory was named `attack/`, so in a checkout of an earlier commit read `attack/` for `search/`.

## The Goethals–Seidel search over ℤ_p

The pipeline of [METHOD.md](../METHOD.md), which found the orders 2092, 3532, 4492 and 5932.

| Program | Role |
|---|---|
| `triage.py` | ranks the multiplier structures of a target n by predicted matching work |
| `generate.py` | generates candidate sequences by batched local search (JAX on TPU, GPU or CPU) |
| `dedupe.c` | removes duplicate candidates |
| `expand.c` | adds the multiplier images of each candidate |
| `match.c` | extended k-tree matcher; prints complementary quadruples |
| `certify.py` | writes `certificate.json`, `sequences.txt` and the matrix for a match |

**Example: order 4492.**

```
cc -O2 search/dedupe.c -o dedupe
cc -O3 -march=native -fopenmp search/match.c -o match -lm
python3 search/triage.py 1123
```

The best row is H = ⟨60⟩ (order 33, 35 orbits, 17 classes) with row sums (1, 1, 1, 67). Generate candidates for each
row sum, remove duplicates, and match. Roles A and D form the left pairs; B supplies both middle sequences, so the
repeated row sum goes in the second file. The target is −(p mod 4).

```
python3 search/generate.py 1123 60 --rowsum 1  --count 160000000 --out r1.bin  --greedy 3 --polish 1 --width 4 --rounds 0
python3 search/generate.py 1123 60 --rowsum 67 --count 160000000 --out r67.bin --greedy 3 --polish 1 --width 4 --rounds 0
./dedupe 17 r1.bin u1.bin
./dedupe 17 r67.bin u67.bin
./match -s 6 -t -3 u1.bin u1.bin u67.bin 17 > solutions.txt
python3 search/certify.py 1123 60 "$(head -1 solutions.txt)" 4492
```

Small structures (pools up to about 10¹⁰ sequences) need only light descent, as here. For orders like 2092 use the
default settings of `generate.py`, which add iterated local search; they run at about 1.3·10⁵ candidates per second on
one TPU v6e chip for p = 523. `match` holds about 20 bytes per list element plus the pair lists of one value of τ, so
lists of 10⁹ elements fit in 172 GB. Choose the clamp `-s` so that the expected number of pairs per τ, printed at
start-up, is close to the list size.

## Index

**Exact joins and local search** (shared by many searches)

| Program | Role |
|---|---|
| `find4.c`, `find4b.c`, `join4x.c` | exact four-list joins of class-value vectors |
| `count4.c`, `count_exact.py`, `exact4.py` | exact solution counts |
| `kt4.c`, `match4.c` | four-list k-tree on sampled lists |
| `join3.c` | three-level exact join for cyclotomic-product structures |
| `lattice.py` | integer lattices: Hermite normal form and membership |
| `orbit_ls.py`, `orbit_ls.c`, `gs_ls.c`, `ls/` | tabu search, simulated annealing and parallel tempering for sequence families |
| `dm/` | relaxed-reflect-reflect (Elser) searches for Goethals–Seidel quadruples |
| `psd_pool.c` | orbit vectors whose power spectrum stays below a bound |
| `workdir.py` | the work directory and the compilation of the C helpers |
| `gsgroup.py`, `gdfs.c`, `osa.c` | orbit tables of H-invariant blocks over any abelian group, their exhaustive lists (PAF and cross vectors), and orbit annealing |

**Goethals–Seidel structures for the unresolved orders**

| Program | Role |
|---|---|
| `gs_lib.py` | Goethals–Seidel structures over abelian groups with a group of automorphisms |
| `gs_sweep.py`, `gs_sweep_subgroups.py`, `gs_exhaustive.py`, `had_sweep.py` | exhaustive sweeps of small multiplier structures |
| `lottery_rank.py` | ranks small structures by Gaussian estimate and exact pair work |
| `had_chain_sweep.py` | the compression chain of `chain.py` for plain Goethals–Seidel structures |
| `cyclo_product.py`, `cp_lowe.py`, `cp_enum.c`, `cp_field.py`, `cp_sweep.py`, `cp_hunt.py`, `cp_breadth.py` | cyclotomic products ℤ_m × GF(q) and the low-energy k-tree (3732, 5020, skew orders) |
| `spin_walk.c` | the walk with neighbour lookups of [METHOD.md](../METHOD.md) (2284) |
| `d3_model.py`, `d3/` | the d = 3 structure (X₀, X₁, μX₁, μ²X₁): Gaussian model, and the walk for n = 1893 |
| `gs_1719.py`, `gs_anti.py`, `gs_array_gf2.py` | individual structures and array tests |

**Bordered arrays** (2148, 3964)

| Program | Role |
|---|---|
| `bordered.py`, `border3.c` | bordered Goethals–Seidel arrays with Q = ℤ₃ and a border of size 12, and their border solver |
| `bgs.py`, `bgs_gauss.py`, `bgs_uw.py`, `bcp.py` | bordered cyclotomic products, their Gaussian estimates and the ℤ₂-split of order 2148 |
| `border_probe.py`, `z4triage.py` | border problems over ℤ_q, and the families over ℤ₄ |

**Spence's construction** (2292)

| Program | Role |
|---|---|
| `spence_spin/` | the search for order 2292; see [spence_spin/README.md](spence_spin/README.md) |
| `spence_lib.py`, `spence_singer.py`, `spp.py` | Spence structures with symmetry, over ℤ₂³ × ℤ_m and over ℤ_m × GF(q) |
| `spence_sweep.py`, `spence_rank.py`, `spence_slice_check.py`, `calibrate.py` | sweeps, ranking, admissibility, and exact counts against the Gaussian estimate |

**Skew-Hadamard matrices**

| Program | Role |
|---|---|
| `skew_triage.py`, `skew_rank.py`, `skew_rank2.py` | triage and ranking of skew multiplier structures |
| `skew_exhaustive.py`, `skew_sweep.py`, `skew_ktree.py`, `skew_lowE.py`, `lowE_plan.py`, `lowE_sweep.py`, `skew_gen.c` | exact and k-tree searches, and low-energy candidate lists |
| `chain.py`, `chain_sweep.py`, `census_screen.py`, `screen_chain.py`, `skew_compress.py`, `comp_ls.c`, `comp_margins.py` | the exact compression chain and the census screen |
| `paley_product.py`, `paley_config.py`, `pcg.py` | Paley-product structures |
| `skew_kmult_sweep.py` | cyclotomic products with an extra multiplier group on ℤ_m |
| `skew_group.py`, `skew_group_sweep.py`, `skew_field.py`, `skew_noncyclic.py` | non-cyclic groups and fields |
| `spin_skew.py`, `spin_compress.py`, `spin_ls.c`, `spin_ls_gls.c`, `spin_pt.c`, `spin_sol.py`, `twin_skew.py`, `pairspin.py` | spin structures: blocks permuted by a multiplier (3676) |
| `spin_ls2.c`, `spin_rrr.py`, `spin_check.py` | the vectorised spin tabu search (`SPIN_LS=spin_ls2`), the difference map for spin structures, and a checker that rebuilds a spin solution |
| `spin_fibre.py` | exact fibre test for spin structures over composite ℤ_v (subgroup classes) |
| `spin_group.py` | spin structures over ℤ_p² × ℤ_q (automorphisms from GL(2, p) × ℤ_q*): ranking, tabu search, and the fibre test on subgroup classes (2107) |
| `tri_skew.py`, `tri_ls.py` | three skew blocks in one multiplier orbit (A, A∘μ, A∘μ², D): exact joins, estimates and tabu search |
| `szekeres.py`, `szek_pair.c`, `szek_group.c` | Szekeres difference sets and Szekeres pairs |
| `szek2.py` | two-block skew (and symmetric) matrices of order 2(w + 1) with a general second block |
| `spence_plane.py` | the plane pair of Spence's skew 4(1 + q + q²) construction, from the Singer set of PG(2, q²) |
| `make_yamada.py` | the Yamada–Kiyasu certificates (24124, 24412, 29404, 51580) |

**Symmetric Hadamard matrices**

| Program | Role |
|---|---|
| `propus.py`, `propus_group.py`, `propus_chain.py`, `propus_gchain.py`, `propus_bgroup.py`, `propus_paley.py` | propus quadruples (A, B, B, D): cyclic, non-cyclic, with compressions, bordered |
| `propus_rank.py`, `propus_sat.py`, `propus_ls.c` | estimates, a CP-SAT model, local search |
| `prow.c`, `p65.c`, `pz5.c`, `hit2rec.py` | dedicated row searches (732) |
| `qborder_clique.c` | the border problem of the ℤ_q-bordered propus array |
| `sym_mm1.py`, `spropus.c`, `symod.c` | symmetric Hadamard matrices of order m(m − 1) and m(m + 3) (4828, 5700) |
| `symod_g.py`, `symod_engine.py`, `symod_join.py`, `symod_reach.py`, `symtpl.c` | symmetric ODs in Balonin form over abelian groups for the weighted plug-in |

**Orthogonal designs for Seberry's theorem**

| Program | Role |
|---|---|
| `odm_engine.py`, `odmx.c`, `odm_join.py`, `odm_join2.py`, `odmc.c` | exhaustive searches for skew-type OD(4t; 1, m, 4t − m − 1) |
| `odm_g.py`, `odm_estimate.py` | the same over any abelian group (`gsgroup.py`), and a Gaussian estimate per template |
| `od13.c`, `od13b.c`, `od13x.c`, `odm_sa.c` | annealing and exhaustive searches for small m |
| `thm7_census.py` | census values reachable by Seberry's theorem with m ≥ 3 |
| `od_anti.py` | involutive anti-automorphisms of the theorem's designs and amicable pairs (the symmetric 7956 test) |

**Complex Hadamard matrices**

| Program | Role |
|---|---|
| `ssxx.py`, `ssxx_group.py`, `ssxx_sample.py`, `ssxx_sweep.py`, `ssxx_table.py` | Szöllősi's array: cyclic, non-cyclic and sampled searches, and the sweep |
| `bsz.py`, `bsz_group.py`, `bsz_sweep.py` | the bordered Szöllősi array |
| `tsz.py`, `tsz.c` | twisted Szöllősi quadruples (tsz.c optionally restricts the streamed side to one compression) |
| `cchain4.py`, `cchain_sweep.py` | exact compression chains for four-block problems (310, 518) |
| `ssxx_paley.py` | exact search on Paley-product structures ℤ_m × GF(q): skeletons, then joins over ℤ_m |
| `ssxx_paleyk.py` | the same with an extra multiplier group K on ℤ_m (exhaustive over K-orbit skeletons) |
| `kstream.py`, `kstream.c` | the streaming k-tree (382) |
| `lowe_sz.py`, `lowe_enum.c`, `lowe_scan.py`, `sls_sz.c` | low-energy k-tree and annealing for Szöllősi structures |
| `comp_est.py`, `comp_triage.py` | estimates conditioned on the block compressions to ℤ_d (exact feasibility, lattice Gaussian per tuple), and the triage of the open orders with them (plain and bordered) |
| `lowe_comp.py` | the low-energy k-tree on one compression class at a time (lowe_enum with fibre-class constraints) |
| `comp_scan.py` | fast conditioned scans of the product families: bordered ℤ₂ × ℤ_q, plain d · q, propus d · q for the symmetric census, and `plainall`/`bordered` over every multiplier group and every listable quotient |
| `tsz_comp.py` | compression-conditioned triage of the twisted Szöllősi structures, and their exact search restricted to the feasible compressions (`sweep`: triage and search every open N) |
| `gcomp.py` | the same conditioning for the bordered array over non-cyclic groups (bsz_group's spin structures on ℤ₂² × ℤ_m, quotient ℤ₂²) |
| `tbsz.py` | the character-bordered Szöllősi array (complex-hadamard/TBORDERED.md): structures on any G with a quotient Q and character set T, conditioned estimates (`rank`), the low-energy search on one compression class, the border solver, certificates (`certify`), and `test` |
| `border_exh.c` | the exhaustive border search of the character-bordered array (compatibility and orthogonality bit matrices, rows depth-first, then a column clique); called by `tbsz.solve_border_exh` (322) |
| `qlp_bin.py` | quaternary Legendre pairs with one binary sequence |
| `paley_complex.py` | quaternary Hadamard matrices of order (q + 1)/2 from Paley's skew Hadamard matrix, q ≡ 3 (mod 4) a prime power (complex-hadamard/PALEY.md): `build` (exact check), `certify`, `list` |
| `miyamoto_wt.py` | Williamson-type matrices of order q and complex Hadamard matrices of order 2q from a symmetric commuting pair of order (q − 1)/2 (Seberry–Yamada Lemma 8.6), the pair from a doubly symmetric propus quadruple, and the completion: a complex Hadamard matrix of order 2q from any quaternary Hadamard matrix of order (q − 1)/2 ([../NOTES.md](../NOTES.md) §19), and Cooper–Wallis with the Williamson-type matrices of Corollary 8.8(ii) from a symmetric Hadamard matrix of order (q − 1)/2 (`cw`, order 21036): `selftest`, `vpatterns`, `complete`, `cw` |
| `quaternion_type.py` | Paley's matrix in the basis of a semiregular Q₈ (quaternion type): block pattern, base-point changes, plugging into Cooper–Wallis arrays ([../NOTES.md](../NOTES.md) §21): `pattern`, `rebase`, `plug` |
| `quaternary_check.py` | exact check of a quaternary Hadamard matrix in the text format of Zaitsev's page, and of its real double (Turyn): `--double` |
| `complex_status.py`, `data/` | known, constructed and open complex orders; tables parsed from Seberry–Yamada (1992) |

**T-matrices, Williamson-type matrices and Legendre pairs**

| Program | Role |
|---|---|
| `tmatrix.py`, `tmat_mult.py`, `tmat_scan.py`, `tmat_composite.py`, `tmat_t3_model.py` | T-matrices by support splits and multiplier structures |
| `bt/` | bordered T-matrices |
| `williamson.py` | Williamson matrices with multiplier symmetry |
| `bnw.py`, `bnw_enum.c` | bordered Williamson-type matrices ([williamson-type](../williamson-type)) |
| `wt_twist.py`, `twist_join.c`, `twist_scan.py` | the order-3 twist of near-Williamson matrices (111) |
| `nw_prep.py`, `nw_exh.c`, `nw_certify.py` | near-Williamson quadruples on ℤ_v with a multiplier group, searched exhaustively: every power-spectrum-feasible triple of the symmetric blocks in a hash table, every member of the arbitrary block streamed against it (split into symmetric and antisymmetric parts, whose autocorrelations add); Cooper–Wallis certificates from a quadruple |
| `lp_search.py`, `lp_sweep.py`, `lp_hunt.py`, `lp_lowe.py`, `lp_mult.py`, `lp_open.py`, `lp_rank_big.py` | Legendre pairs for the unresolved orders |
| `tlp_model.py`, `tlp_rank.py`, `tlp_search.py` | twisted Legendre pairs |
| `two_circ.py`, `two_chain.py` | two-block searches: D-optimal designs and Legendre pairs |

**Paley's matrix in sixteen circulant blocks** (3156, 5260, …, 11804)

| Program | Role |
|---|---|
| `paley_blocks.py` | Paley's matrix of order q + 1 (q ≡ 3 mod 8, a prime or a prime power) as a 4 × 4 array of circulants along the orbits of a torus element, times T-matrices of order t from base sequences or from a file (Seberry–Yamada 1990, Theorem 4), over ℤ_t × ℤ_s when gcd(s, t) > 1; writes the certificate ([../hadamard/PALEY_BLOCKS.md](../hadamard/PALEY_BLOCKS.md)): `check`, `quad`, `certify` |

**Certificates**

| Program | Writes |
|---|---|
| `certify.py`, `gs_certify.py` | Goethals–Seidel orders in `hadamard/` |
| `skew_certify.py`, `skew_certify_group.py`, `skew_readme.py` | skew orders in `skew/` |
| `sym_certify.py`, `border_certify.py` | symmetric orders in `symmetric/` |
| `bnw_certify.py` | Cooper–Wallis orders with bordered Williamson-type matrices |
| `ssxx_certify.py` | complex orders in `complex-hadamard/` |
