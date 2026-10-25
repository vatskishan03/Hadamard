# Order 36w from a split of an Hadamard matrix of order 4w

The programs behind [../NOTES.md](../NOTES.md) §97. The AG(2, 3) hyperplane blocks of
[../conference118](../conference118/README.md) and [../flag9p](../flag9p/README.md), applied to an Hadamard matrix K of
order N = 4w instead of a conference core, give H_ij = K_ij(J − 2M_ij). In the ℤ₃²-invariant case this is a known
plug-in. The matrices Y_a = J − 2L_a (L_a the subgroup of ℤ₃² orthogonal to the direction a) satisfy:

- Y_a Y_bᵀ = J for a ≠ b;
- Σ Y_a Y_aᵀ = 36 I.

They form Seberry and Whiteman's regular 4-set of order 9 (Seberry–Zhang, Graphs Combin. 9 (1993)). If K = K₀ + K₁ + K₂
+ K₃, with the K_a disjoint weighing matrices W(4w, w), then H = Σ K_a ⊗ Y_a is an Hadamard matrix of order 36w. H is
symmetric if every K_a is.

The open orders it would reach are n = 1719 = 9 · 191 (from an H(764)) and the census values v = 9w with a known
symmetric H(4w). No split of an Hadamard matrix of order 764 was found, and none of a symmetric matrix at all.

| File | What it does |
|---|---|
| `split4.py` | the blocks Y_a, the plug-in `plug(parts)` and the weighing check |
| `csplit.py` | CP-SAT for circulant splits: Goethals–Seidel quadruples (one layer label per position; each layer needs Σ PAF = wδ) and propus quadruples (symmetric) |
| `demo_small.py` | the small cases: splits of orders 3 and 11 with the H(108), H(396) checks; the infeasible propus cases |
| `gen_small.py`, `gsgen.c` | small Goethals–Seidel and propus quadruples (meet in the middle; simulated annealing) |
| `sa_split.c` | simulated annealing for circulant splits, with label swaps that keep every layer at weight w; slots can be H-orbits |
| `symw.py`, `symsplit_any.py` | symmetric weighing matrices W(n, k); a symmetric Hadamard matrix with a symmetric split (CP-SAT) |
| `dk764.py` | Đoković's two Goethals–Seidel quadruples of order 191 (arXiv:math/0703312), rebuilt and checked |
| `split1719.py`, `mkinput.py`, `verify1719.py` | H-invariant splits of them (CP-SAT), input files for `sa_split`, and the check of a split and of the H(6876) it would give |

Reproduction of the checks of §97:

```
python3 demo_small.py                      # orders 3 and 11 split (H(108), H(396) checked); propus w = 5..11 infeasible
python3 dk764.py                           # both quadruples: sum of PAFs = 764 delta
gcc -O2 -o sa_split sa_split.c -lm && python3 mkinput.py 0 orb && ./sa_split in_0_orb.txt 1 100000000 30 0.5 4
```
