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

Build (in this directory):

```
gcc -O3 -march=haswell -o qskew qskew.c
cp qskew qskew2
gcc -O3 -march=haswell -o qp qp.c -lm
mkdir -p shapes out runs
```

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
