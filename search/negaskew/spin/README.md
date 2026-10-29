# Spin shapes in the negacyclic Goethals–Seidel problem

The quaternary skew problem of [../README.md](../README.md) with the spin structure C = B ∘ μ, D = B ∘ μ², where μ ∈
ℤ_v*, μ ∉ H and μ³ ∈ H. The problem reduces to a skew block A and one free block B: a two-list problem of 2^(n−1) ·
4^(n−1) quadruples instead of four lists ([../../NOTES.md](../../NOTES.md) §88, §89). The programs read the shape files
and tables of the parent directory.

## Contents

| File | Contents |
|---|---|
| `qspin.c` | the plain exact search: the list of all skew blocks A hashed by the autocorrelation vector against every B (4^(n−1), B(0) = 1), with the row-sum filter \|R_A\|² + 3 \|R_B\|² = 4v; for n ≤ 13 |
| `qspin_fast.c` | the same search for up to n = 21 orbits: the row-sum condition met in the middle over the low and high digits of B, a two-class presence filter, a bilinear fast path for the class sums; equal counts and statistics |
| `spin_drv.py` | the shape file of (v, H, μ) for the engines, `run_spin`, and a driver that checks a solution at the group level and as a real skew matrix (v ≤ 60) |
| `enum_spin.py`, `spin_est.py`, `spin_filter.py`, `spin_cost.py` | the spin shapes (H, μ) of the census-undecided v from `../enum_n21_v2499.json`, their ranking by the spin estimate, the row-sum filter (`../dcnp.py`), and the number of B that pass the row sums at n = 21 |
| `spin_sweep.py` | runs the engine on a shape list (resumable, one json line per shape) |
| `spin_cert.py` | a solution (A, B) to the four orbit-exponent strings and `../certify.py` |

## Build and controls

In this directory:

```
gcc -O3 -march=haswell -o qspin qspin.c
gcc -O3 -march=haswell -o qspin_fast qspin_fast.c
mkdir -p shapes runs
python3 spin_drv.py 7 "[]" 2        # COUNT 1820, real skew matrix of order 56
python3 spin_drv.py 13 "[]" 3       # COUNT 110448, order 104 (3 s; the plain qspin takes 30 s)
```

`SPIN_ENGINE=./qspin` selects the plain engine in `spin_drv.py`. Shape files go to `shapes/`, shape lists and logs to
`runs/`; neither is stored in the repository.

## Reproduction

The commands that reproduce the spin results of [../../NOTES.md](../../NOTES.md) §88 and §89 are in the section
[Reproduction](../README.md#reproduction) of [../README.md](../README.md).
