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
| `data/run379_partial.out`, `data/run379_residue170258.out` | the partial run of the main session at v = 379 (see below) and the 32 solutions of the residue that holds the matrix of order 3032 |

Build (in this directory; the structure scripts read SOL files of `../qskew`, built in `..`; Python 3 with NumPy and SymPy):

```
gcc -O3 -march=haswell -o tri_join tri_join.c
(cd .. && gcc -O3 -march=haswell -o qskew qskew.c && mkdir -p shapes out)
mkdir -p runs
```

The structure scripts read SOL lines of `../qskew` (orbit order of `../qdrv.shape(v, gens)`), not the SOL lines of `tri_join` (coset order, with a field `r=R`); `tri_cert.py` converts the latter.

## Reproduction of search/NOTES.md §92 (triangular family; to be completed with the main session's section)

Run here after the builds. The tools were written for v = 379 (k = 27), the first member after the two hits that is in reach; the smaller members are the controls.

**Estimates** (HEURISTIC; the theory analysis of round 39 measured exact/Gaussian = 0.97 at 67, 0.49 at 191, 0.41 at 277 and 0.38 at 379 with the counts 6192, 99,360, 359,424 (277, EXHAUSTIVE there, not re-run here) and 1.25 · 10⁵ ± 25% (379)):

```
python3 tri_est.py                      # 3 s    67 11:  log10E 3.807 (exact 3.792), rank 5, covol 55.4;   191 19: log10E 5.293 (exact 4.997), rank 9, covol 1144.9
python3 tri_est.py 379 27 631 35        # 9 s    379: log10E 5.562 (rank 13, covol 21674);   631: 5.029 (rank 17, covol 393216)
python3 coset.py                        # 1.5 s  67, 191, 379, 631, 947: rank 5, 9, 13, 17, 21, '0 in coset: True' for all (no congruence obstruction)
```

**Exact counts, two ways** (the exact orbit-level join `qskew` and the partitioned join; the counts are normalised: B, C, D with u(0) = 1, A with u(0) ∈ {1, −i}, ordered):

```
(cd .. && python3 -c "from qdrv import *; d,o,e=run(67,[64],1,100000,0); open('tri/runs/s67.out','w').write(o)")   # 0.3 s  H = <64> of order 11: 'COUNT 6192 (normalized: B,C,D have w0=1 ...)', 6192 SOL lines
./tri_join 67 11 2 32 0 32 100000 > runs/j67.out        # 0.04 s  |A| = 48, |L| = 400; 6192 SOL lines (grep -c SOL), the qskew count
./tri_join 191 19 19 1024 0 1024 1000000 > runs/j191.out  # 133 s  |A| = 640, |L| = 63,504; found = 99,360 (4.03e9 probes), the count of the theory analysis of round 39 (qskew at 191 needs more than four minutes here)
```

**v = 379, k = 27 (e = 14, H = ⟨87⟩ of order 27, 15 orbits): the matrix of order 3032 = 8 · 379** [EXACT existence]. One residue of M = 2¹⁸ holds a cluster of 32 solutions (the tilde of B, C, D and the swap of C, D give 2 · 2 · 2 · 2 images, A ↦ −A with u(0) 1 ↔ −i the factor 2):

```
./tri_join 379 27 2 262144 170258 170259 1000 > runs/j379.out 2> runs/j379.err     # 2 min 16 s wall (1 min 14 s CPU): |A| = 8960, |L| = 11,778,624, left = 806,598 pairs, probes 1,058,559,234, found 32
python3 tri_cert.py 379 27 2 "SOL r=170258 021120313003120 013021333112200 020023311133131 013110010232332"
   # orbit exponents (qdrv order): 021110220330013 013201021330132 020103321311313 013310312002312;  G-level sum PAF: peak 1516, max|off-peak| 0.0;  A twisted-skew True;  row sums [24, -18, -18, -42]
   # sha256 a142deedaab955e19ba601251dc09d814db24141e6965a6f6201cd32dbeebddb -> ../out/cert_3032.json
```

Checked: the certificate in a scratch copy of the repository (`skew/3032/certificate.json`): `python3 verify/verify.py skew/3032` → `ok  H H^T = 3032 I  sha256 a142deed…  (skew: H + H^T = 2I)` (3 s); and an independent numpy construction from the sequences alone
(negacyclic matrices N[i, j] = x_{j−i} for j ≥ i, −x_{n+j−i} for j < i, the Goethals–Seidel array of [../../../skew/NEGACYCLIC.md](../../../skew/NEGACYCLIC.md), int64): entries ±1, H Hᵀ = 3032 I exactly, H + Hᵀ = 2I, a₀ = 1, a palindromic (80 s).
The same matrix (same sha256) was found independently by the theory analysis of round 39 with a different join (`tri2.c`: 448 solutions in all in its random sample of partitions). The order 3032 is not new: 379 is absent from Tables 4 and 5 of Cati–Pasechnik (an absent n means m = 2: Hadamard and skew-Hadamard matrices of order 4 · 379 = 1516 are known, so 3032 follows by doubling) and from the 95 open orders ≤ 10 000 of Egan–Kotsireas–London–Ó Catháin.

**The partial run `data/run379_partial.out` is not a negative result.** It is `tri_join 379 27 2 262144 0 28 …` of the main session: the residues 0 to 27 of M = 2¹⁸ (|A| = 8960, |L| = 11,778,624, about 10⁹ probes per even residue: 14 even residues, about 807,000 left pairs each), found = 0 in each.
It is 14 of the 2¹⁷ even residues; every solution comes with 31 others in its residue, so the number of independent clusters expected there is 0.4 (exact count) to 1.1 (Gaussian estimate) and P(no hit) is 0.33 to 0.66. The count at the hit residue 170258 above is the same as the one of the theory analysis (found = 32).
The partition is by a linear residue of the PAF vector, so the solutions are spread over the residues almost at random (in that analysis's random sample of partitions all 14 hit partitions contain exactly 32 solutions).

**Structure of the solutions at v = 67** (the first lines of `ana.py`, `leg.py`, `flat.py`, `rel.py`, `spec.py`, `psd.py` on `runs/s67.out`; `python3 ana.py 67 "[64]" runs/s67.out` etc.):

```
ana.py 67 [64] runs/s67.out        # 0.5 s  v 67 k 11 n 7 sols 6192;  row-sum pattern |R|^2 = (265, 1, 1, 1) for all 6192  (4v - 3 = 265)
leg.py 67 [64]                     # 0.3 s  the Legendre PAF off 0 is -1 (not 3); 'A with PAF=3 off 0: 0 best max deviation 5.657'
flat.py 67 [64] runs/s67.out       # 2 s    blocks with 6 distinct spectral values on the 6 cosets: 3600 of 6192 (all four), 576 + 576 + 576 with one block with 3, 192 x 3 with one with the single value 68, ...
rel.py 67 [64] runs/s67.out        # 0.8 s  the three small blocks are never equal up to unit and tilde (equalpairs 0); self-dual small blocks: 5184 solutions none, 864 one, 144 two
spec.py 67 [64] runs/s67.out 1     # 3 s    support sizes of the coset DFT per block (4, 5, 5, 5) 2592, (6, 5, 5, 5) 1728, …
psd.py 67 [64] runs/s67.out 1      # 3 s    distinct spectral values per block: 1 value 576, 2 values 288, 3 values 1728, 6 values 22,176 (over all 24,768 block spectra)
```

The theory analysis of round 39 shows with these and its own tools (EXACT or VALIDATED there): no flat or Legendre-type free block beyond v = 67, no spin or shift relation between the free blocks beyond e = 6, no sparsity of the coset DFT, no character-type blocks; a free block's spectrum on the cosets is not flat (the expected number of flat blocks is 179, 61, 0.76 at 67, 191, 379).
The family is expected to be finite: the expected number of inequivalent solutions (E ρ/(96 e), ρ ≈ 0.4) is 25 (631), 9.5 (947), 21 (1129), 1.3 (1327), 0.17 (2017), 0.008 (2557) (HEURISTIC); whether 631, 947, 1129, 1327 carry solutions is UNKNOWN (the join needs 5 · 10¹⁵ probes at 631).

### Deviations from the round-39 sources

The code of the main session (`r39/tri`) is carried over unchanged except: `ana.py`, `spec.py`, `rel.py`, `flat.py`, `psd.py` and `leg.py` import `qdrv` from `..` (they put `.` on the path before) and have a docstring; `coset.py` and `tri_join.c` have a docstring/header; no change of any computation.
`tri_cert.py` is new. Not carried over: the binary `tri_join`, `__pycache__`, the sparse `search/`, `skew/` copies of the repository in the source directory, `run379.pid`. The tools of the theory analysis of round 39 (`tri2.c`, the estimators `est_fast.py`, `cop.py`, `hyb.py`, the proofs and checks `check379.py`, …, with the 277 and 379 runs) are not part of this port.
