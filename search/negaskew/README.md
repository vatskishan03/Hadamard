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

Build (in this directory):

```
gcc -O3 -march=haswell -o qskew qskew.c
cp qskew qskew2
mkdir -p shapes out
```

The runs that found the two matrices (exhaustive joins, stopping at the first solution):

```
python3 -c "from qdrv import *; print(run(191, [180], 1, 1, 1)[1])"     # about 4 s
python3 certify.py 191 "[180]" "32032213013 03001230212 02232120000 03303111220"
python3 -c "from qdrv import *; print(run(277, [273], 1, 1, 1)[1])"     # about 3 minutes
python3 certify.py 277 "[273]" "3313220001321 0301322112300 0023212321000 0033103132201"
python3 indep_nega.py out/cert_1528.json
```

`certify.py` writes `out/cert_1528.json` and `out/cert_2216.json`, whose sequences are those of
[../../skew/1528](../../skew/1528) and [../../skew/2216](../../skew/2216).

The sweep: `python3 enum_shapes.py 13 2499` lists the shapes in `enum_n13_v2499.json`, and
`python3 sweep.py 13 2499 1500 runs/sweep.jsonl 1` runs them (1500 s per shape; mode 1 quaternary, 0 real).
