# The triangular family v = 1 + k(k + 1)/2: estimate, partitioned join, structure tests

The tools of [../../NOTES.md](../../NOTES.md) §92 for the skew-Hadamard matrices of order 8v from the negacyclic array of [../README.md](../README.md) in the family of §88: v a prime of the form 1 + k(k + 1)/2 with k ≡ 3 (mod 4)
(7, 29, 67, 191, 277, 379, 631, 947, 1129, 1327, 2017, …), H the subgroup of ℤ_v\* of order k, e = (k + 1)/2 cosets g^j H (the orbits of the non-zero points) and one orbit {0}. Both hits of §88, 191 (k = 19) and 277 (k = 23), are of this form: three row sums are units
and the fourth is u(0) + k(±1 ± i), so |R|² = 4v − 3 = k² + (k + 1)² for one block; A (the skew block) is that block iff k ≡ 3 (mod 8). The four quaternary blocks are constant on the cosets (one digit per orbit), A twisted-skew with u(0) ∈ {1, −i}.

| File | Contents |
|---|---|
| `tri_est.py` | the Gaussian and lattice estimate (HEURISTIC) of the number of normalised solutions of the pattern (A big with coset sum m(1 + i), B, C, D small with coset sum 0 and u(0) = 1): sampled PAF vectors of big and small blocks, their covariance and the covolume of the lattice they span; `python3 tri_est.py [v k ...]` (always 67 11 and 191 19) |
| `tri_join.c` | the partitioned four-list join: the A list (all skew blocks with coset sum m(1 + i), m = ±1 aligned with u(0): u(0) = 1 with m = +1, or u(0) = −i with m = −1; pairs of cosets j, j + e/2 carry c and 1 − c) and the list L of all B, C, D with u(0) = 1 and coset sum 0 (as many digits 0 as 2 and as many 1 as 3 among the e cosets); each element carries its PAF vector (per class of coset shifts the differences of the digits counted over the pairs of orbits: two integers), a 64-bit random linear key of it and a residue π mod M, M a power of two; the left pairs (A, B) with π = r against the right pairs (C, D) with π = −r by a hash join on the key, every match verified exactly; `tri_join v k g M r0 r1 maxsol` runs the residues r0 … r1 − 1 (g a primitive root mod v; one residue = one partition, so the work splits and resumes) |
| `tri_cert.py` | **new in this port**: certificate of a SOL line of `tri_join` (digits in coset order, converted to the orbit order of `../qdrv.py`), `../certify.py`, the SHA-256 of the matrix: `python3 tri_cert.py v k g "SOL r=R A B C D"` writes `../out/cert_<8v>.json` |
| `ana.py`, `spec.py`, `rel.py`, `flat.py`, `psd.py`, `leg.py`, `coset.py` | the structure tests: the row-sum patterns of the solutions (`ana.py`), the coset-DFT supports and magnitudes of the blocks (`spec.py`), relations between the small blocks (equal up to unit and the tilde B(g) ↦ conj B(−g), self-dual blocks: `rel.py`), the number of distinct spectral values per block on the cosets and whether a block is flat (`flat.py`), the spectra themselves (`psd.py`), whether any skew block has PAF = 3 off the origin (a Legendre-type A, `leg.py`) and whether 0 lies in the coset of the lattice of the PAF vectors (the PAF equation is not excluded modulo the lattice: `coset.py`) |
| `PROOFS.md`, `proofs_check.py` | **new in this port**: the proofs P1 to P7 of §92 (parity of the PAF, the row sums, the skew block as two sign vectors, the line on which its spectrum lies with the corrected P4(c), the symmetry group of order 96e and the clusters of solutions, the period and Gauss-sum forms of the equations, the lattice of the PAF vectors); `proofs_check.py` tests every statement on the exact data (`python3 proofs_check.py [P1 … P7]`, 15 s, needs `tri2`) |
| `tri2.c`, `run_tri2.sh` | the cache-resident partitioned four-list join of the second analysis of §92: the partition is by a residue modulo a PRIME M of a linear form of the PAF vector, the left table of one partition stays in the cache, every hit is verified in integers inside the program; mode 3 (k ≡ 3 mod 8: A big, B, C, D small) and mode 7 (k ≡ 7 mod 8: A small, B the big free block with S = 1 + i, C, D small), optional reduction `red=1` (A over the shift orbits, C, D over the tilde classes, exact weights); `tri2 v k g mode red M r0 step cnt idx0 maxsol seed` (header of `tri2.c`; e ≤ 15); `run_tri2.sh TAG v k g mode red M r0 step cnt maxsol seed` is the resumable driver (appends to `runs/TAG.out`, resumes at the number of its STAT lines). Output: `SOL i=I r=R w=W A B C D` (digits u(0), c_0 … c_{e−1} in coset order), one `STAT` line per partition, `TOTAL` |
| `est_fast.py`, `est_all.py`, `legest.py`, `hyb.py`, `cop.py`, `rho_r.py`, `est_logs.py`, `tricore.py` | the estimates (HEURISTIC): `est_fast.py KMAX` the Gaussian estimate with the covolume 2^(e−1) of P7 for every member up to k = KMAX (`python3 est_fast.py 100`); `est_all.py v,k …` with the lattice step; `legest.py` the Gaussian expectation of solutions with a Legendre-type flat free block; `hyb.py v k [FILE]` the exact list of A with the Gaussian density of the free triple, per A against an exact set, and the share of the mass of the A with a spectral peak above x·v; `cop.py v k [FILE]` its Gamma-marginal Gaussian-copula correction (scipy); `rho_r.py v k N1 N2 seed [rmax]` the collision sampler of ρ(r) = (empirical probability that the first r PAF coordinates of a random quadruple vanish)/(Gaussian); `est_logs.py M FILE …` the estimate Σw·M/(partitions) from the STAT lines of `tri2` logs (and `--combine`); `tricore.py` the sampler class they share |
| `trilib.py`, `mk_shape.py` | `trilib.py`: coset structure of ℤ_v, `load_sols` (SOL lines of `../qskew` and of `tri2`), spectra and PAF of a block; `mk_shape.py v k` writes `../shapes/tri_<v>.txt` for `../qskew` |
| `check379.py`, `resid379.py`, `test_trijoin379.py`, `find277.py`, `classes_red.py` | `check379.py FILE [NMAX [v k]]` the floating-point (FFT) check of every SOL line of a `tri2` log and its row-sum pattern; `test_trijoin379.py` the residue 170258 of the cluster under the hash of `tri_join.c`; `resid379.py` the distribution of the hit residues mod 2, …, 64; `find277.py [FILE]` the published solution of 2216 reduced as `tri2 red=1` does, and its line in a reduced 277 log; `classes_red.py v k FILE` the classes of a reduced log under the group of order 96e and the sums of their weights |
| `symm.py`, `symm2.py`, `stab.py` | the symmetry group of P5 on a solution set: closure under the shifts, the tilde of each free block and the swaps (`symm.py`), under A → −A with the classes of the full group (`symm2.py`), the stabilisers of a sample (`stab.py`) |
| `stats1.py`, `stats2.py`, `gauss_dom.py`, `spinrel.py`, `astat.py`, `astruct.py`, `charblk.py` | the structure tests of the second analysis: real cosets and positions of the big block (`stats1.py`), flat free blocks and the number of distinct spectral values per kind of block (`stats2.py`), sparsity of the ℤ_e-DFT, spectral twins and the reality of PAF sums (`gauss_dom.py`), free blocks related by a shift, the tilde or the conjugation (`spinrel.py`), which A occur: score, spectral peak, per decile (`astat.py`), invariants of the sign vectors of A (`astruct.py`), blocks with c_j = i^(aj+b) (`charblk.py`, new in this port) |
| `data/run379_partial.out`, `data/run379_residue170258.out` | the partial run of the main session at v = 379 (see below) and the 32 solutions of the residue that holds the matrix of order 3032 |
| `data/c277_red.out`, `data/c277_red_M1013.out`, `data/c277_unred.out` (each with `.err`) | the recorded logs of the three complete runs of `tri2` at v = 277 (k = 23, e = 12, g = 5, mode 7): reduced, M = 1009, seed 1 (all 646 reduced solutions printed); reduced, M = 1013, seed 2; unreduced, M = 8191, seed 3 (7.29 · 10¹¹ probes); each `TOTAL … wsum=59904` |
| `data/s379.out`, `data/c379_red.out` (each with `.err`, `.log`) | the recorded logs of `tri2` at v = 379 on the same sequence of partitions (r0 = 12345, step 100003, M = 524287, seed 1): unreduced, 1484 partitions (14 hit partitions, 448 solutions printed); reduced (`red=1`), 7454 partitions (5 hits of weight 224); the run was stopped there |
| `data/sol379_qskew.txt`, `data/cand277.txt`, `data/r_sol.txt` | the solution of the matrix of order 3032 in the orbit order of `../qskew` (`python3 ../certify.py 379 "[87]" "$(cat data/sol379_qskew.txt)"` from `..`), the candidate printed by `find277.py`, the residue printed by `test_trijoin379.py` |
| `data/rho_191_a.out`, `data/rho_379_a.out`, `data/est_fast_100.out` | recorded output of `rho_r.py 191 19 2e7 2e7 …` (165 s) and of the 379 run (713 s), and the output of `est_fast.py 100` of this port |

Build (in this directory; the structure scripts read SOL files of `../qskew`, built in `..`, or of `tri2`; Python 3 with NumPy and SymPy, SciPy for `cop.py` only):

```
gcc -O3 -march=haswell -o tri_join tri_join.c
gcc -O3 -march=haswell -o tri2 tri2.c
(cd .. && gcc -O3 -march=haswell -o qskew qskew.c && mkdir -p shapes out)
mkdir -p runs
```

The structure scripts read SOL lines of `../qskew` (orbit order of `../qdrv.shape(v, gens)`) or of `tri2` (coset order, `SOL i=I r=R w=W A B C D`; `trilib.load_sols` recognises the format); the SOL lines of `tri_join` (coset order, field `r=R`) are for `tri_cert.py`, which converts them.
Every program of this directory takes its paths from its own location (it can be run from any directory); the logs of fresh runs go to `runs/`, the recorded ones are in `data/`. One core, `nice -n 5`, `ulimit -v 3500000`, `OMP_NUM_THREADS=1` for every command below; the times are wall times on a machine shared with another job (CPU time in parentheses when it differs a lot).

## Reproduction of search/NOTES.md §92 (triangular family)

Run here after the builds. The tools were written for v = 379 (k = 27), the first member after the two hits that is in reach; the smaller members are the controls. Every command below was re-run in the staging tree and printed the output shown, except those marked "not re-run": the long runs, which are given with their recorded logs (`data/`), and the independent numpy construction of the certificate (80 s), which keeps the output documented when it was written.
The proofs of the statements marked PROVED in §92 are in [PROOFS.md](PROOFS.md); `python3 proofs_check.py` (15 s) checks them on the exact data.

**Estimates** (HEURISTIC). The exact counts are 5952 (29), 6192 (67), 99,360 (191), 359,424 (277, below) and about 10⁵ at 379 (below); the Gaussian estimate with the covolume 2^(e−1) (P7) gives 10^3.891 = 7780, 10^3.807, 10^5.304, 10^5.946 and 10^5.524, so exact/Gaussian = 0.77, 0.97, 0.49, 0.41 and 0.24 to 0.47 at 379 (from the two samples of the logs, below; 0.28 combined):

```
python3 tri_est.py                      # 3 s    67 11:  log10E 3.807 (exact 3.792), rank 5, covol 55.4;   191 19: log10E 5.293 (exact 4.997), rank 9, covol 1144.9
python3 tri_est.py 379 27 631 35        # 9 s    379: log10E 5.562 (rank 13, covol 21674);   631: 5.029 (rank 17, covol 393216)
python3 coset.py                        # 1.5 s  67, 191, 379, 631, 947: rank 5, 9, 13, 17, 21, '0 in coset: True' for all (no congruence obstruction)
python3 est_all.py 67,11 191,19         # 2 s    67: log10E 3.806, rank 5, coset_ok True;   191: 5.304, rank 9, coset_ok True   (the lattice step; est_fast.py omits it)
python3 est_fast.py 100                 # 11 s   log10E for every v = 1 + k(k+1)/2 prime, k = 3 mod 4, k <= 99 (data/est_fast_100.out)
      #  v:   7    29     67     191    277    379    631    947    1129   1327   2017   2557   2851   4561    4951
      #  E: 1.92  3.891  3.807  5.304  5.946  5.524  5.027  4.700  5.074  3.893  3.128  1.858  0.243  -3.250  -4.891
python3 legest.py                       # 3 s    expected number of solutions with a Legendre-type flat free block: 67: 179, 191: 60.9, 379: 0.757, 631: 0.0018, 947: 7.1e-6
python3 hyb.py 191 19 runs/s191.out     # 3 s    exact A list, Gaussian free triple: E_hybrid 1.64 x exact (see 'Exact counts' for runs/s191.out); max peak/v among the used A 3.033, among the unused >= 3.797; A with peak > 4v: 18.8% of the list, 5.4% of the hybrid mass
python3 hyb.py 67 11 ; python3 hyb.py 379 27   # 1 s, 3 s   A with a peak > 4v (dead by positivity, PROOFS.md P4): 25.0% of the list at 67 (5.7% of the mass), 51.6% at 379 (16.6%); E_hybrid 6525 (67), 3.0e5 (379)
python3 cop.py 191 19 runs/s191.out     # 8 s    Gamma-marginal copula: E_cop = 0.939 x exact at 191;   python3 cop.py 67 11 runs/s67.out: 0.737;   python3 cop.py 379 27: E_cop = 1.00e5 (E_gauss_hyb 3.02e5), 4296 of the 8960 A dead by positivity   (11 s)
python3 rho_r.py 191 19 2e6 2e6 1 9     # 45 s   rho(r), r = 1 … 8: 0.979 0.962 0.922 0.893 0.823 0.752 0.688 0.603  (the recorded run with 2 · 10⁷ pairs, 165 s: data/rho_191_a.out: 0.980 0.965 0.924 0.896 0.825 0.751 0.700 0.606 0.586; at 379, data/rho_379_a.out, 713 s: 0.992 1.018 1.006 1.003 0.937 0.902 0.855 0.729 0.62 (± 25%))
```

The first four equations behave like independent Gaussians (ρ ≈ 1 at 379), the loss accumulates over the last ones; there is no exactly dependent equation beyond the weighted sum rule (P6). The expected number of inequivalent solutions E ρ/(96 e) with ρ = 0.4: 24.6 (631), 9.5 (947), 20.6 (1129), 1.25 (1327), 0.17 (2017), 0.008 (2557) (arithmetic on the table above; with ρ = 0.28, the value that the 379 samples give, 0.7 times these).

**Exact counts, two ways** (the exact orbit-level join `qskew` and the partitioned join; the counts are normalised: B, C, D with u(0) = 1, A with u(0) ∈ {1, −i}, ordered):

```
(cd .. && python3 -c "from qdrv import *; d,o,e=run(67,[64],1,100000,0); open('tri/runs/s67.out','w').write(o)")   # 0.3 s  H = <64> of order 11: 'COUNT 6192 (normalized: B,C,D have w0=1 ...)', 6192 SOL lines
./tri_join 67 11 2 32 0 32 100000 > runs/j67.out        # 0.04 s  |A| = 48, |L| = 400; 6192 SOL lines (grep -c SOL), the qskew count
./tri_join 191 19 19 1024 0 1024 1000000 > runs/j191.out  # 133 s  |A| = 640, |L| = 63,504; found = 99,360 (4.03e9 probes), the count of `tri2` below (qskew at 191 needs more than four minutes here)
```

The same counts with `tri2` (`tri2 v k g mode red M r0 step cnt idx0 maxsol seed`, header of `tri2.c`; M = 7, 31, 1009 are primes, the partitions are independent):

```
./tri2 29 7 2 7 0 7 0 1 7 0 0 1                       # 0.01 s  mode 7: TOTAL hits=992 (x 6 positions and signs of the big block = 5952, the qskew count)
./tri2 67 11 2 3 0 31 0 1 31 0 0 1                    # 0.03 s  mode 3: TOTAL hits=6192
./tri2 191 19 19 3 0 1009 0 1 1009 0 100000 1 > runs/s191.out   # 13 s (7 s CPU)  TOTAL parts=1009 probes=4032758016 hits=99360: the count of tri_join, 1.8 ns per probe against about 70 ns for tri_join (its run of the residue 170258 below: 1.06e9 probes in 74 s CPU); the SOL lines of runs/s191.out are the exact set used below
(cd .. && python3 -c "from qdrv import *; d,o,e=run(29,[16],1,100000,0); open('tri/runs/s29.out','w').write(o)")   # 0.3 s  COUNT 5952 (normalized: B,C,D have w0=1 ...)
```

**v = 277, k = 23 (e = 12, g = 5, H = ⟨273⟩ of order 23, 13 orbits): the exact count 359,424 = 6 · 59,904** [EXHAUSTIVE, three complete runs of `tri2` in mode 7, from the recorded logs, not re-run; mode 7 counts B = the big free block with S = 1 + i only, so the normalised count is 3 positions of the big block x 2 (S = 1 ± i) x the weight]. The three runs use different moduli and hash seeds, two reduce by symmetry and one does not:

```
./run_tri2.sh c277_red 277 23 5 7 1 1009 0 1 1009 100000 1     # data/c277_red.out:       TOTAL parts=1009 probes=91314224925 hits=646 wsum=59904   |L| = 853,776, |A| = 2560 (216 shift orbits), |Big| = 731,808
./tri2 277 23 5 7 1 1013 0 1 1013 0 0 2                        # data/c277_red_M1013.out: TOTAL parts=1013 probes=91314224925 hits=646 wsum=59904
./tri2 277 23 5 7 0 8191 0 1 8191 0 0 3                        # data/c277_unred.out:     TOTAL parts=8191 probes=728933458176 hits=59904 wsum=59904   (no symmetry reduction, 7.29e11 probes; the hits per partition are multiples of 8, 16 in 2322 of the 2996 hit partitions: P5, group of order 16 in mode 7)
python3 est_logs.py 1009 data/c277_red.out ; python3 est_logs.py 8191 data/c277_unred.out    # 'COMPLETE scan, all 1009 (8191) partitions: ... exact normalised count (sum of the weights) = 59904'
```

Re-run here: the first residues (4 s, the output is the first 5 partitions of the recorded log, byte for byte) and the whole reduced run (6 min 40 s wall, 4.6 CPU-minutes; its `runs/c277_red.out` is byte-identical to `data/c277_red.out`):

```
./tri2 277 23 5 7 1 1009 0 1 5 0 100000 1                      # 4 s   STAT i=0 r=0 left=156646 probes=90505020 hits=0;  i=1 hits=0;  SOL i=2 r=2 w=96 3330311221200 0231020012020 0210203313210 0033012100222 (and one more), hits=2 wsum=192;  i=3 hits=1 wsum=96;  i=4 hits=1 wsum=96;  TOTAL parts=5 hits=4 wsum=384
./run_tri2.sh c277_red 277 23 5 7 1 1009 0 1 1009 100000 1      # 6 min 40 s (4.6 CPU-min): TOTAL parts=1009 probes=91314224925 hits=646 wsum=59904; diff runs/c277_red.out data/c277_red.out: identical
python3 check379.py data/c277_red.out 2 277 23                  # 1 s    checked 646 all ok True ; row-sum patterns {(1, 1105, 1, 1): 646}   (FFT check of every reduced solution; 4v - 3 = 1105)
python3 classes_red.py 277 23 data/c277_red.out                 # 1 s    distinct classes under the full group (Z_e x tilde^3 x S_3 x negA): 313 ; reduced hits per class [(2, 308), (6, 5)] ; class weight sums [(96.0, 2), (192.0, 311)] total 59904:   311 classes of 1152 solutions and 2 of 576 make 359,424 (the 5 classes with 6 reduced hits are exactly those whose A has a shift symmetry, stabiliser {4, 8}: six hits of weight 32 instead of two of weight 96)
python3 find277.py data/c277_red.out                            # 1 s    candidate (1021359, 59776364, (8548308, 14886472)); found in data/c277_red.out : SOL i=672 r=672 w=96 3322111233000 0321131000123 0113332120020 0201212030230 ; matches in the log: 1  (the published solution of 2216 is among the 646)
```

**v = 379, k = 27 (e = 14, H = ⟨87⟩ of order 27, 15 orbits): the matrix of order 3032 = 8 · 379** [EXACT existence]. One residue of M = 2¹⁸ holds a cluster of 32 solutions (the tilde of B, C, D and the swap of C, D give 2 · 2 · 2 · 2 images, A ↦ −A with u(0) 1 ↔ −i the factor 2):

```
./tri_join 379 27 2 262144 170258 170259 1000 > runs/j379.out 2> runs/j379.err     # 2 min 16 s wall (1 min 14 s CPU): |A| = 8960, |L| = 11,778,624, left = 806,598 pairs, probes 1,058,559,234, found 32
python3 tri_cert.py 379 27 2 "SOL r=170258 021120313003120 013021333112200 020023311133131 013110010232332"
   # orbit exponents (qdrv order): 021110220330013 013201021330132 020103321311313 013310312002312;  G-level sum PAF: peak 1516, max|off-peak| 0.0;  A twisted-skew True;  row sums [24, -18, -18, -42]
   # sha256 a142deedaab955e19ba601251dc09d814db24141e6965a6f6201cd32dbeebddb -> ../out/cert_3032.json
```

The same sequences from the orbit-order solution of the data (from `..`: `python3 certify.py 379 "[87]" "$(cat tri/data/sol379_qskew.txt)"`, 3 s): G-level sum PAF peak 1516, max|off-peak| 0.0, A twisted-skew True, ring-realisation matrix of order 3032 (True, True, True), row sums [24, −18, −18, −42].

Checked: the certificate in a scratch copy of the repository (`skew/3032/certificate.json`): `python3 verify/verify.py skew/3032` → `ok  H H^T = 3032 I  sha256 a142deed…  (skew: H + H^T = 2I)` (3 s); and an independent numpy construction from the sequences alone
(negacyclic matrices N[i, j] = x_{j−i} for j ≥ i, −x_{n+j−i} for j < i, the Goethals–Seidel array of [../../../skew/NEGACYCLIC.md](../../../skew/NEGACYCLIC.md), int64): entries ±1, H Hᵀ = 3032 I exactly, H + Hᵀ = 2I, a₀ = 1, a palindromic (80 s).
The same matrix (same sha256) was found independently with a different join (`tri2`, below: the 14 hit partitions of a random sample of partitions hold 448 solutions, the cluster of 32 at i = 37 is the one above). The order 3032 is not new: 379 is absent from Tables 4 and 5 of Cati–Pasechnik (an absent n means m = 2: Hadamard and skew-Hadamard matrices of order 4 · 379 = 1516 are known, so 3032 follows by doubling) and from the 95 open orders ≤ 10 000 of Egan–Kotsireas–London–Ó Catháin.


**The sample of partitions at v = 379 and the independent clusters** (`tri2` in mode 3; the logs are `data/s379.out` (unreduced, 27 minutes of the original run by the timestamps of its files) and `data/c379_red.out` (reduced, 73 minutes), from the recorded logs, not re-run in full; the first partitions were re-run):

```
./run_tri2.sh s379 379 27 2 3 0 524287 12345 100003 100000 1000000 1      # data/s379.out: 1484 partitions, 448 SOL lines in 14 hit partitions, every one with hits=32
./run_tri2.sh c379_red 379 27 2 3 1 524287 12345 100003 524287 1000000 1  # data/c379_red.out: 7454 partitions (stopped), 5 hit partitions each with hits=1 wsum=224
./tri2 379 27 2 3 0 524287 12345 100003 38 0 1000 1                       # 68 s wall (37 s CPU): |L| = 11,778,624, |A| = 8960; the output is the first 38 partitions of data/s379.out byte for byte: STAT i=37 r=42447 left=202380 probes=264602746 hits=32 wsum=32 and its 32 SOL lines;  TOTAL parts=38 probes=10055645186 hits=32
./tri2 379 27 2 3 1 524287 12345 100003 98 0 1000 1                       # 54 s wall (28 s CPU): |L'| = 5,891,028 tilde classes, 640 shift orbits of A; the first 98 partitions of data/c379_red.out byte for byte: SOL i=97 r=275470 w=224 …, hits=1 wsum=224 (the same partition i = 97, r = 275470 as the unreduced hit)
python3 check379.py data/s379.out                         # 1 s    checked 448 all ok True ; row-sum patterns {(1513, 1, 1, 1): 448}   (FFT; the integer check is inside tri2)
python3 test_trijoin379.py                                # 1 s    sum of PAF vectors [0 0 ... 0]; M 262144 left residue r = 170258  right residue 91886 = -r   (the residue of tri_join 379 27 2 262144 170258 170259: the cluster of the matrix of order 3032 is the cluster of i = 37 of tri2)
python3 resid379.py                                       # 1 s    the residues of the first 300 hit solutions under the hash of tri_join: r mod 2 all 0 (an even residue), mod 4: 140/160, mod 64: 10 values of 32 or 12 solutions: spread almost at random
python3 est_logs.py 524287 data/s379.out data/c379_red.out   # data/s379.out: 1484 partitions, 14 hit partitions, weight 448: E = 1.58e+05 +- 27%;   data/c379_red.out: 7454 partitions, 5 hit partitions, weight 1120: E = 7.88e+04 +- 45%
python3 est_logs.py --combine 524287 data/s379.out data/c379_red.out   # 9.45e+04 +- 35%: the unreduced log for the partitions 0 to 1483 (all solutions), the reduced log for 1484 to 7453 (weight 896)
python3 symm2.py 379 27 data/s379.out                     # 1 s    closed under negA: True; classes under the full group: 12 (class sizes [(32, 11), (96, 1)]: the sample of 14 clusters, not the set of classes)
```

The normalised count at 379 is HEURISTIC (a sample of partitions, not a count): the two logs give 1.6 · 10⁵ ± 27% (14 clusters in 1484 partitions) and 7.9 · 10⁴ ± 45% (5 reduced hits in 7454 partitions). The two samples overlap (the same first 1484 partitions; the only reduced hit among them, i = 97, lies in the cluster of the unreduced hit at i = 97, which also holds the images under A → −A), so they are not pooled by adding; using the unreduced log for the partitions 0 to 1483, where it holds every solution, and the reduced log for the others gives 9.5 · 10⁴ ± 35%, that is exact/Gaussian = 0.28 (0.24 to 0.47 from the two samples; the pooled value 1.25 · 10⁵ ± 25% quoted earlier, ratio 0.38, is inside the combined interval but is not reproduced by these three ways of pooling). About 70 inequivalent solutions (9.5 · 10⁴ / 1344).
The exact count at 379 would be the sum of the weights of the reduced scan over all 524,287 partitions: `./run_tri2.sh c379_red 379 27 2 3 1 524287 12345 100003 524287 1000000 1` continues the recorded log (about 11 core-hours for the remaining partitions, not run).

**The partial run `data/run379_partial.out` is not a negative result.** It is `tri_join 379 27 2 262144 0 28 …` of the main session: the residues 0 to 27 of M = 2¹⁸ (|A| = 8960, |L| = 11,778,624, about 10⁹ probes per even residue: 14 even residues, about 807,000 left pairs each), found = 0 in each.
It is 14 of the 2¹⁷ even residues; every solution without a stabiliser comes with 31 others in its residue (PROOFS.md P5; the 14 hit partitions of the sample of `tri2` hold exactly 32 each), so the number of independent clusters expected there is E/(32 · 2¹⁷) · 14 = 0.3 to 0.5 (the exact-count estimates 9.5 · 10⁴ to 1.6 · 10⁵) or 1.1 (Gaussian estimate), and P(no hit) is 0.3 to 0.7. The count at the hit residue 170258 above is the same as that of `tri2` (found = 32: `test_trijoin379.py`).
The partition is by a linear residue of the PAF vector, so the solutions are spread over the residues almost at random (in the random sample of partitions of `tri2` all 14 hit partitions contain exactly 32 solutions).

**Structure of the solutions at v = 67** (the first lines of `ana.py`, `leg.py`, `flat.py`, `rel.py`, `spec.py`, `psd.py` on `runs/s67.out`; `python3 ana.py 67 "[64]" runs/s67.out` etc.):

```
ana.py 67 [64] runs/s67.out        # 0.5 s  v 67 k 11 n 7 sols 6192;  row-sum pattern |R|^2 = (265, 1, 1, 1) for all 6192  (4v - 3 = 265)
leg.py 67 [64]                     # 0.3 s  the Legendre PAF off 0 is -1 (not 3); 'A with PAF=3 off 0: 0 best max deviation 5.657'
flat.py 67 [64] runs/s67.out       # 2 s    blocks with 6 distinct spectral values on the 6 cosets: 3600 of 6192 (all four), 576 + 576 + 576 with one block with 3, 192 x 3 with one with the single value 68, ...
rel.py 67 [64] runs/s67.out        # 0.8 s  the three small blocks are never equal up to unit and tilde (equalpairs 0); self-dual small blocks: 5184 solutions none, 864 one, 144 two
spec.py 67 [64] runs/s67.out 1     # 3 s    support sizes of the coset DFT per block (4, 5, 5, 5) 2592, (6, 5, 5, 5) 1728, …
psd.py 67 [64] runs/s67.out 1      # 3 s    distinct spectral values per block: 1 value 576, 2 values 288, 3 values 1728, 6 values 22,176 (over all 24,768 block spectra)
```


**Classes and the symmetry group** (P5; the exact sets of 29 and 67 from `qskew`, of 191 from `tri2`; `python3 symm2.py v k FILE`):

```
python3 symm2.py 29 7 runs/s29.out      # 1 s    5952 solutions closed under negA: True;  P_A = P_negA(A) on 200 random solutions: True;  24 classes under Z_e x tilde^3 x S_3 x negA, sizes [(96, 2), (192, 14), (384, 8)]
python3 symm2.py 67 11 runs/s67.out     # 2 s    13 classes, sizes [(144, 1), (288, 3), (576, 9)]   (symm.py: closed under shifts, tilde of each free block, swaps; 26 classes without negA)
python3 symm2.py 191 19 runs/s191.out   # 33 s   99,360 solutions: 105 classes, sizes [(480, 3), (960, 102)]   (symm.py: 210 classes without negA, 6 + 204; 420 A used of 640, 42 up to shifts; 5830 free blocks used of 63,504)
python3 stab.py 67 11 runs/s67.out      # 16 s   stabilisers: a single tilde (351 solutions), two tildes (50) among the first 401 stabilised ones;  stab.py 191 19 runs/s191.out (24 s): single tildes only (29 solutions)
```
The sizes are 96e (384, 576, 960) except for the solutions with a self-dual free block and, at e ≤ 6, equal free blocks. At 379 every hit partition of the sample holds one orbit of 32 (`proofs_check.py P5`).

**Structure of the solutions at v = 29, 191, 277, 379** (the tests of the second analysis; the 277 and 379 files are the recorded logs, so they hold the reduced list (277) and 14 clusters (379), not the whole set; HEURISTIC statistics of EXACT sets):

```
python3 stats2.py 29 7 runs/s29.out       # 1.6 s  flat small blocks: 960 of 5952 solutions have one
python3 stats2.py 67 11 runs/s67.out      # 2 s    flat free blocks: 576 solutions have one, 5616 none
python3 stats2.py 191 19 runs/s191.out    # 27 s   num flat small blocks per solution {0: 99360};  python3 stats2.py 277 23 data/c277_red.out: {0: 646};  python3 stats2.py 379 27 data/s379.out: {0: 448}   (a flat free block only at v = 29 and 67, P4(c) excludes three of them)
python3 spinrel.py 29 7 runs/s29.out      # 1.3 s  free blocks related by shift/tilde: 1056 of 5952 solutions;  67: 576 of 6192;  191: 0 of 99,360 (33 s);  277: 0 of 646;  379: 0 of 448   (no spin or shift relation beyond e = 6)
python3 gauss_dom.py 379 27 data/s379.out # 2 s    support of the Z_e DFT: A 14 (448 of 448), small 13 (1344 of 1344), random small blocks 13 (1994 of 2000); no spectral twins (0 of 1344 pairs);  191 (2000 solutions): support e − 1 = 9 for 97.0% of the small blocks (98.75% of random small blocks), A 10 or 6;  277: 2 of 1938 pairs of free blocks with equal |c-hat|
python3 charblk.py 191 19 runs/s191.out   # 4 s    blocks with c_j = i^(aj+b): none at 191, 277 (646), 379 (448); at 67 only a = 2 (the Legendre type: 288 blocks in each of the positions B, C, D), at 29 also a = 1
python3 astat.py 191 19 runs/s191.out     # 2 s    the 420 used A of 640: max_m |F_A(m)|^2 <= 579.25 = 3.033 v, the unused >= 725.2 = 3.797 v;  at 67: used <= 177.5 = 2.65 v, unused >= 233.4 = 3.48 v   (positivity: > 4v is dead)
python3 astruct.py 191 19 runs/s191.out   # 2 s    used/unused A by (u0, <sigma, sigma'>) and by the parity of the minus signs of sigma': no separation
python3 stats1.py 191 19 runs/s191.out    # 4 s    big block at position A in all 99,360; at 277: position B in all 646 (mode 7)
```

**The proofs** (`python3 proofs_check.py`, 15 s: P1 parity of 70,200 (sequence, shift) pairs; P2 the pattern (4v − 3, 1, 1, 1) on the exact sets of 29 (992) and 67 (6192), on 448 solutions at 379 and 646 at 277, and the two patterns (25, 1, 1, 1) x 80 and (13, 13, 1, 1) x 168 at v = 7; P3 the list sizes 16, 48, 640, 2560, 8960, 129,024 of A at 29 … 631; P4 `conj U = (1+i) − iU` to 10⁻¹², Σρ and Σρ², and the arithmetic of P4(c) for every k up to 10⁶; P5 the 15,360 and 23,040 images of 40 solutions at 29 and 67 and the 14 orbits of 32 at 379; P6 the period form and the Gauss-sum form at 29, 67, 379; P7 rank e − 1 and index 2^(e−1) at 29, 67, 191, 277, 379, 631, 947, 1129): `ALL CHECKS PASSED`.

The second analysis shows with these tools (EXACT or VALIDATED on the sets above): no flat or Legendre-type free block beyond v = 67, no spin or shift relation between the free blocks beyond e = 6, no sparsity of the coset DFT, no character-type block beyond e = 6; a free block's spectrum on the cosets is not flat (the expected number of flat blocks is 179, 61, 0.76 at 67, 191, 379).
The family is expected to be finite: the expected number of inequivalent solutions (E ρ/(96 e), ρ ≈ 0.4) is 25 (631), 9.5 (947), 21 (1129), 1.3 (1327), 0.17 (2017), 0.008 (2557) (HEURISTIC, table above); whether 631, 947, 1129, 1327 carry solutions is UNKNOWN. The probes needed to meet one solution with random sub-lists of the free blocks (fraction f of each list with E f³ = 1, so f²|L|² probes with |L| = C(e, e/2)²): 4.6 · 10¹⁵ at 631 (|L| = 2.36 · 10⁹, E ρ = 4.3 · 10⁴), 3 · 10²⁰ at 947 (|L| = 4.98 · 10¹¹, E ρ = 2 · 10⁴); the plain scan of all partitions costs |L|² = 5.6 · 10¹⁸ and 2.5 · 10²³.

### Deviations from the sources

The code of the main session (`r39/tri`) is carried over unchanged except: `ana.py`, `spec.py`, `rel.py`, `flat.py`, `psd.py` and `leg.py` import `qdrv` from `..` (they put `.` on the path before) and have a docstring; `coset.py` and `tri_join.c` have a docstring/header; no change of any computation.
`tri_cert.py` is new. The tools of the second analysis (`tri2.c`, `run_tri2.sh`, the estimators, the structure and symmetry tests, the logs of the 277 and 379 runs) are carried over from its working directory with these changes only: docstrings that name the sections of NOTES.md instead of rounds; the paths of the program's own directory (`sys.path` and `cd` instead of the working directory and the scratch directory); `lib.py` → `trilib.py` and `core.py` → `tricore.py` (generic names); the driver `run_tri2.sh` calls `./tri2` (the binary `tri2b` of the original was the same source), without `taskset` and with `nice -n 5`; `tri2.c`: a guard that refuses e > 15 (the digits of a block are packed into 32 bits; the original would overflow silently), header comments; `trilib.load_sols` also reads the SOL lines of `tri2`; `check379.py` takes optional `v k` (the source was fixed to 379 27); `find277.py` writes `runs/cand277.txt` and takes an optional log; `resid379.py` and `test_trijoin379.py` read `data/s379.out` by default; `mk_shape.py` writes `../shapes/`; `est_fast.py`: `float()` around two printed numbers (numpy 2 printed `np.float64(...)`); `classes_red.py`: one more output line (the class weight sums). No computation changed. New: `PROOFS.md`, `proofs_check.py`, `est_logs.py`, `charblk.py`.
Not carried over: the binaries (`tri2b` is `tri2`), `__pycache__`, the repository copy in the source directory of the second analysis, its copies of `ana.py … tri_est.py` and `tri_join.c` (identical to the files here), its notes (`STATUS.md`, `PRODUCTION.md`: the commands are in the reproduction block above), `run379.pid`, the logs `c277_unred_M1013_partial.out`, `trijoin191_bigM.out`, `test379_red.*` (partial or test runs), `trijoin379_rsol.*` (equal to `data/run379_residue170258.*` but for an exit line), `verify3032.out`, `*.pid`. The statement of the README that "every solution comes with 31 others in its residue" is made precise (P5: 31 others unless the solution has a stabiliser; 15 in mode 7).
