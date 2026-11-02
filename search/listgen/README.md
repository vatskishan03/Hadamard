# listgen: generated and twin-reduced lists for the low-energy k-tree (skew orders 1004 and 3044)

Programs behind [search/NOTES.md](../NOTES.md) §111 and the skew-Hadamard matrices of orders 1004 = 4 · 251
([skew/1004](../../skew/1004)) and 3044 = 4 · 761 ([skew/3044](../../skew/3044)). Both are Goethals–Seidel families over
a prime field ℤ_p with a multiplier group H of odd order and a skew first block, searched by the low-energy k-tree of
[NOTES.md](../NOTES.md) §20 and §110 (`search/cp_lowe.py`, `search/cp_enum.c`, `search/kt4.c`): a member of a block has
a key, its vector of autocorrelations on the orbits of ⟨H, −1⟩; the members whose key is closest to the conditional mean
given that the four keys sum to zero (whitened energy) are kept, and `kt4` joins the four lists. The programs here
change how the lists are made:

- **v = 251** (`lg_lowe.py`, `lg_lib.py`, `gen_run.py`, `lg_pt.c`). The blocks B, C and D have C(50, 24) = 1.2 · 10¹⁴
  and C(50, 22) = 8.9 · 10¹³ members, too many for `cp_enum.c`, so the low-energy members are generated: Metropolis
  chains over the orbit signs (`lg_pt.c`) collect distinct members below an energy threshold, and the Gray-code
  enumeration of `lg_pt.c` lists the 2²⁴ skew first blocks up to reversal. The join is `kt4`, as in `cp_lowe.py`.
- **v = 761** (`cp_lowe_dd.py`, `dd_run.sh`). The lists can be enumerated by `cp_enum.c`, but reversing a block
  (x(t) → x(−t)) leaves its autocorrelations unchanged, so every list holds the twins of each member. `cp_lowe_dd.py` is
  `cp_lowe.py` with the option `--dedup`, which keeps one member of each twin: it halves every list and cuts the pair
  work of the join. `dd_run.sh` is the resumable wrapper that ran it.

Python 3 with NumPy and SciPy, and a C compiler with OpenMP. No certificate depends on these programs: `verify/` checks
both matrices on their own. Hit logs, lists and binaries go to `$HADAMARD_WORK` or to `--workdir`, and `lg_lowe.py`
writes `HIT_<v>_<e>_<cfg>.json` into the current directory.

Labels: CHECKED = exact computation or direct measurement; HEURISTIC = estimate or sampling statistics, nothing depends
on it for correctness.

## Files

The files are byte for byte the ones that ran (md5 below), except `lg_lib.py`, `gen_run.py` and `dd_run.sh`: in the
files that ran, a few lines name the working directory on the machine of the run, and here those lines are replaced (md5
of the published file, then the md5 of the file that ran; see "Paths").

| File | md5 | Contents |
|---|---|---|
| `lg_pt.c` | `8e9fab9877257775e78d0dd42e455d18` | the list generator: for one block, written as a quadratic form in the ±1 class signs (`key_c(y) = const_c + Σ lin y + Σ_{i<j} S y_i y_j`), Metropolis chains at one or several temperatures (swap moves of a +1 with a −1 keep the class sum, single flips for the skew block) that collect distinct members of whitened energy at most τ, one of each reversal twin (class mask and its rotation by e/2), as records (class mask, int16 key); with mode 3 the exhaustive Gray-code enumeration of the 2^(e/2) skew blocks. Built with `cc -O3 -march=native -ffast-math -fopenmp lg_pt.c -o lg_pt -lm` |
| `lg_lib.py` | `2662b4539b3cb78670bc6d9cae013633` (as run `92721f7d60d12fe8b874b0479d632402`) | the helpers of the generator: the quadratic form of a block, the keys of sign vectors, the mean and covariance of the keys of uniform members (`moments`), the conditional mean and the whitening (`energy_setup`), the export file read by `lg_pt`, the reversal-canonical class masks |
| `gen_run.py` | `d0f14e225bf572f7669c0101f4159335` (as run `72aa3b2851863e1bb05115fcc175e7bc`) | `tau_for`, the energy threshold at which the Gaussian model puts a given number of members of a block inside the ball (noncentral χ² quantile), and a stand-alone run of `lg_pt` on one block for measuring its rate: `gen_run.py v e cfg k Lball seconds R Tmin Tmax nprop threads [seed] [maxdistinct]` |
| `lg_lowe.py` | `415f7f57091a6e61ad4ba423c3f1f2a7` | the driver: stage 1 and the moments, the four lists (block A by Gray code, B, C and D by chains; `--share-bc` uses the list of B for C when the two row sums agree; `--kB`, `--kC`, `--kD` keep the members of lowest energy), the Gaussian estimate inside the lists, the ranking of the split values, the join by `kt4` (resumable through `kt4_cells.err`; `--join-core-h` caps the core-hours, `--K` the number of solutions), and the decoding of a solution to four ±1 sequences, checked (autocorrelations, skew first block, invariance under H) and appended to `$HADAMARD_WORK/skew_hits.jsonl` |
| `cp_lowe_dd.py` | `2b9f668d4a99969bebc8ea4e2b3bc317` | `search/cp_lowe.py` (md5 `4e2f4c952a56d1d0942b5c1c1a7d4d70`, 801 lines) with 816 lines: `--dedup` keeps one member of each reversal twin in every list (the skew block with sign + on its first class, B, C and D the members whose class mask is at most its rotation by e/2; needs `--stream`), and `--only-sigma` joins one given split value (a validation option) |
| `dd_run.sh` | `be9c5fafeb317108cde6f8ca2f5a58af` (as run `f9a321939f884403ede551e9a0355432`) | the resumable wrapper of `cp_lowe_dd.py`: `dd_run.sh WORKROOT <arguments of cp_lowe_dd.py>` sets `HADAMARD_WORK` to WORKROOT, puts the list files on a RAM disk, writes `round_NN.log` and `cmd.txt`, and passes the split values finished in earlier rounds as `--skip-sigmas` |

The other programs of the two runs are those of commit `919e73d`, byte for byte: `search/cp_lowe.py`
(`4e2f4c952a56d1d0942b5c1c1a7d4d70`), `search/cp_enum.c` (`4c2094d92d03faff5b970b136c837958`), `search/kt4.c`
(`d22c49d986b9b3ba9ee6d3e62ded6fea`), `search/cyclo_product.py` (`98be16555ffe42590a0cd2185add9f3e`),
`search/workdir.py` (`b845fde4e2f09a184b42ec79d7776f4b`), `verify/verify.py` (`970e5d9aba849fd38a7a4109a85a091d`), and,
for the certificates, `search/skew_certify.py` (`813b89691dad146f8d68b217b292659a`) with `search/certify.py`
(`f27dd2834f970c8424a1d6f45a94e896`). `cp_enum` and `kt4` are compiled on first use by `workdir.py`; `lg_pt` is built by
hand and `lg_lowe.py` calls the binary next to it. The machine of both runs: AMD EPYC 9R14, 32 cores, 246 GB, gcc
13.3.0, Python 3.12.3, NumPy 1.26.4. The binaries that ran have the md5 `84a3b7bed9e13368d64e6eca4fbd5435` (`lg_pt`),
`ffe1b35a47b8bc1fb8d2cbb5fda2e986` (`kt4`) and `fd442e1ee5f538eda5e8dd8be1ae524e` (`cp_enum`, run for 3044 only);
rebuilt from these sources on that machine, `lg_pt` and `kt4` have the same md5.

## Paths

The three files named above differ from the files that ran only in the lines below, which put absolute paths of the
machine of the run into the code:

- `lg_lib.py` put two directories of its working copy (`search/` and `verify/`, which hold `cp_lowe.py`,
  `cyclo_product.py`, `workdir.py` and `verify.py`) at the head of `sys.path`; that line is now a comment, and
  `PYTHONPATH=search` supplies the modules (`lg_lowe.py` puts its own directory first, and `cyclo_product.py` adds the
  repository's `verify/` itself).
- `gen_run.py` put the directory that held `lg_lib.py` on `sys.path`; it now puts its own directory.
- `dd_run.sh` took `HERE`, the directory of `cp_lowe_dd.py`, from a fixed path; it is now the directory of the script,
  `HERE=$(cd "$(dirname "$0")" && pwd)`, set before the wrapper changes into WORKROOT.
- `cp_lowe_dd.py` is as it ran. It puts its own directory and the directory `verify` beside it on `sys.path`, as
  `cp_lowe.py` does; in `search/listgen/` neither holds `workdir.py` or `cyclo_product.py`, so it needs
  `PYTHONPATH=search`, as above.

`dd_run.sh` needs Linux (it puts the list files on `/dev/shm`) and absolute paths for WORKROOT and for `PYTHONPATH`,
since it changes into WORKROOT. When it starts it deletes `$WORKROOT/tmp` and `/dev/shm/hadamard-dd-<basename of
WORKROOT>` (and puts a link to the second at the place of the first), and it deletes the second again when it exits. On
another system, or when WORKROOT has a `tmp` of its own, run `cp_lowe_dd.py` directly: the same program without the
round logs, the resumption and the RAM disk.

CHECKED here, with the files as published: with `PYTHONPATH=search`, `lg_lowe.py --help` and `cp_lowe_dd.py --help`
start, and the two controls below run to a verified solution (the second also through `dd_run.sh`, called by a relative
path, on a system with `/dev/shm` replaced by another directory in a copy of the script).

## Usage

```
cc -O3 -march=native -ffast-math -fopenmp search/listgen/lg_pt.c -o search/listgen/lg_pt -lm
export HADAMARD_WORK=$HOME/hadamard-work PYTHONPATH=$PWD/search OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
```

`lg_lowe.py` looks for `lg_pt` next to itself, so the binary is built into `search/listgen/`; it is not part of the
repository.

**Controls** (v = 79, H = ⟨23⟩ of order 3, rows 1, 1, 5, 17; each finds a skew Goethals–Seidel family in a few seconds
on two threads). The first is the generator, the second the twin-reduced k-tree.

```
python3 search/listgen/lg_lowe.py 79 26 0 --LA 1e9 --LB 2e5 --LC 2e5 --LD 2e5 --nB 3e4 --nC 3e4 --nD 3e4 \
    --nsig 20 --samples 100000 --seed 3 --threads 2 --mcpairs 500000 --workdir $HADAMARD_WORK/w79
python3 search/listgen/cp_lowe_dd.py 79 --q 79 --e 26 --config 0 --L 300,300,300 --stream --groups 5 \
    --samples 100000 --nsig 5 --mcpairs 500000 --maxpairs 3e8 --K 1 --seed 1 --threads 2 --dedup
```

**Skew 1004, v = 251** (`--workdir` is the directory for the lists and the join log; the hit is appended to
`$HADAMARD_WORK/skew_hits.jsonl`):

```
python3 search/listgen/lg_lowe.py 251 50 0 --LA 3e5 --LB 7e7 --LC 7e7 --LD 4e8 --nB 1.4e7 --nD 1e7 --kC 4e6 --kD 3e6 \
    --share-bc --resume --keep-lists --secs-gen 1800 --threads 12 --T 0.42 --nsig 6000 --mcpairs 4000000 \
    --maxpairs 3e8 --join-core-h 120 --K 1 --seed 2 --workdir $HADAMARD_WORK/w251
```

H = ⟨219⟩, the 50th powers modulo 251 (order 5, e = 50), config 0 = row sums (1, 9, 9, 29). `--LA`, `--LB`, `--LC`,
`--LD` give the thresholds through the nominal number of members of the Gaussian ball (3 · 10⁵, 7 · 10⁷, 7 · 10⁷ and
4 · 10⁸); `--nB` and `--nD` the number of distinct members the chains collect (14,000,000 and 10,000,000); `--kC` and
`--kD` the members of lowest energy kept (4,000,000 of the list of B for C, 3,000,000 of D); `--T 0.42` the temperature
of the chains (one chain per thread); `--secs-gen` a time cap per block, not reached; `--nsig 6000` split values ranked
by expected solutions per unit of work; `--join-core-h 120` stops the join after 120 core-hours; `--K 1` at the first
solution; `--seed 2` seeds the 4 · 10⁵ moment samples (`--samples`, the default), the chains (10 · seed + k for block k,
and 1000003 times that plus t for thread t) and the 4 · 10⁶ pairs per side that rank the split values. The job ran in a
directory of its own, with `lg_pt` beside `lg_lowe.py`.

Result (CHECKED, from the logs): lists A 343,764 (energy at most 10.69), B 14,000,000 (at most 3.677), C the 4,000,000
lowest of B (cut 3.289), D 3,000,000 of 10,000,008 (cut 3.884), in 557 s (5,617 CPU seconds); the join over five merge
classes (0, 7, 11, 14, 22) searched 52 split values without a solution and found one in the 53rd, 1918 s after the start
on 12 threads, about 4.2 GB at peak. The Gaussian estimate of key-quadruples inside the lists, 10^1.18 (HEURISTIC), is
the one printed by the program. The lists are not reproducible bit for bit (the threads race into one hash table and the
chains stop when the target is reached): a rerun gives lists of the same sizes and thresholds with other members, and
the ranking of the split values follows from the lists.

**Skew 3044, v = 761** (one round; `W` is the absolute path of the run directory, see "Paths"):

```
bash search/listgen/dd_run.sh $W 761 --q 761 --e 40 --config 0 --L 6e8,3e7,3e7 --dedup --stream --groups 5 \
    --samples 2000000 --mcpairs 10000000 --maxpairs 300000000 --nsig 100000 --K 1 --seed 1 --threads 12
```

H = ⟨498⟩, the 40th powers modulo 761 (order 19, e = 40), config 0 = row sums (1, 1, 39, 39). The thresholds keep
590,647,863, 33,937,598 and 30,532,180 members of B, C and D, and the twin reduction 295,323,935, 16,968,798 and
15,266,090, with the 524,288 skew blocks of sign + on the first class. Result (CHECKED, from the logs): the lists in
1076 s, the ranking of the split values at 1083 s, the first 2702 split values of five merge classes (2, 7, 8, 12, 14)
empty, two solutions in the 2703rd, one printed (its limit is 1) and verified, 9955 s (2.77 h) on 12 threads, 30.7 GB at
peak. Without the reduction the same lists needed 29.5 s per split value instead of 3.2 s (377 split values in 11,116 s,
no solution). The Gaussian estimate inside the reduced lists is 10^1.21 per twin class (HEURISTIC).

## Certification

`lg_lowe.py` and `cp_lowe_dd.py` rebuild a solution as four ±1 sequences, check the autocorrelations, the skew first
block and the invariance under H, and append it to `$HADAMARD_WORK/skew_hits.jsonl`. The certificates are written by the
repository's `search/skew_certify.py`, which reads that file, and checked by `verify/verify.py`:

```
HADAMARD_WORK=<work directory of the run> python3 search/skew_certify.py 251      # skew/1004
HADAMARD_WORK=<work directory of the run> python3 search/skew_certify.py 761      # skew/3044
python3 verify/verify.py skew/1004 && python3 verify/verify.py skew/3044
```

`skew_certify.py` takes the first hit for v, asserts the skew first block, the autocorrelation identity (4v at shift 0,
zero elsewhere) and the invariance, writes `certificate.json`, `sequences.txt` and `H<order>.csv.gz`, and rebuilds the
matrix with `verify.py` (H Hᵀ = nI, H + Hᵀ = 2I, SHA-256). The two certificates here were written this way from the hit
records of the runs. The matrices also pass `verify/check_matrix.c` (504,510 and 4,634,490 row pairs, no violation).
