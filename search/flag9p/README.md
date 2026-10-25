# The flag equation for symmetric conference matrices of order 9p + 1

The programs behind [../NOTES.md](../NOTES.md) §96. They extend the ansatz of [../conference118](../conference118/README.md)
(points ℤ_p × AG(2, 3) and ∞, blocks −χ(j − i)(J − 2M) with M a hyperplane of AG(2, 3)², diagonal blocks J − I − 2R)
from p = 13 to every prime p ≡ 1 (mod 4). A symmetric conference matrix of order 9p + 1 is a conference graph on 9p
vertices; Brouwer's table (read 25 October 2026) lists those on 117 (J. Y. Z. Ho) and 153 (Chen Wang) vertices as
existing and those on 261, 333, 369, 477, 549, 657 and 801 vertices (p = 29, 37, 41, 53, 61, 73, 89) as open.

## The equation

On the nontrivial characters of ℤ₃², C² = (9p)I becomes Ψ² = p I for the Hermitian matrix Ψ = D + 2Φ on the 4p flags
(point, direction of AG(2, 3)): D = +1 on the directions (0, 1), (1, 0) and −1 on (1, 1), (1, 2); Φ has exactly one
entry χ(j − i)ω^c in each 4 × 4 block (i, j), i ≠ j, joining the two flags whose characters the hyperplane block maps
onto each other. Every connected component of the flag graph must satisfy Ψ_C² = p I. Consequences:

- a component has an even number of flags, half with D = +1 (its trace is a multiple of √p);
- a bipartite component (edges only between D = +1 and D = −1 flags) is X X* = d I, d = (p − 1)/4, for its biadjacency
  matrix X; a complete bipartite component K_{a,a} has D constant on each side and X is a Butson matrix BH(a, 6) whose
  sign pattern is the Paley pattern;
- a structure invariant under ℤ_p is an 8 × 8 matrix over ℤ[ω][ℤ_p] (fs3's formalism), bipartite ones a 4 × 4 matrix X
  with a real structure;
- if a multiplier group K of odd order k (k | p − 1) is an automorphism, every orbit of differences carries one edge type,
  with a cube-root phase twist when 3 | k (without a twist, X(1) has entries in kℤ[ω], which needs k² | d).

## Programs

| File | What it does |
|---|---|
| `bipcyc.c` | all bipartite ℤ_p-invariant structures (stage 1 rows with a pruned DFS and phase backtracking, stage 2 join on complementary masks); `bipcyc p flip [1]`, the last argument restricting to decoupled rows (one orientation) |
| `build_bip.py` | the conference matrix of order 9p + 1 from a `bipcyc` solution, with the exact check C = Cᵀ, C² = 9p I |
| `alam.py` | the Alamouti sub-family [[P₁, P₂], [−λP₂*tᵏ, λP₁*tᵏ]]: supports in 2k · {1, …, (p − 1)/2} |
| `mulcyc.c` | all ℤ_p-invariant structures (bipartite or not) with a multiplier group of odd order k, up to 4 orbits of ±pairs |
| `kn.c`, `kg.c` | component templates: Hermitian D + 2Φ with Φ of sixth roots on a given support (complete graphs, or any graph), Ψ² = (1 + 4·degree) I |
| `bh.c`, `bhclass.py` | normalised Butson matrices BH(n, 6) and their sign classes (switching and permutations; pynauty) |
| `k8classes.py`, `k8place.py`, `k8full.py` | the switching classes of the K₈ templates; the 8-subsets of ℤ₂₉ with the template's Seidel spectrum; those whose Paley two-graph and D-labelling match the class (`./kn 8 4 100000 > k8_all.txt; python3 k8classes.py k8_all.txt; python3 k8place.py 29; python3 k8full.py`) |
| `kddplace.py` | biclique placements (X₊, X₋) whose Paley sign matrix is in the sign class J − 2P (BH(7, 6)) or rank one (`--rank1`, BH(3, 6)) |

## Validation

- p = 5 (`mulcyc 5 1`): 576 structures, Mathon's order 46.
- p = 13: no ℤ₁₃-invariant structure (`bipcyc`, `mulcyc 13 3`, and the general search of round 44); the matrix of
  [../conference118](../conference118/README.md) is not cyclic.
- p = 17: 55,296 bipartite cyclic structures (27,648 with the first row canonical); every one built and checked with
  `build_bip.py` is a symmetric conference matrix of order 154.

## Results

See [../NOTES.md](../NOTES.md) §96. Bipartite cyclic: p = 29 none; Alamouti: p = 29, 37, 41 none; odd multiplier groups:
none for (p, k) = (13, 3), (29, 7), (37, 9), (53, 13), (61, 15), (109, 27), (73, 9). K₈ templates: one class, 4466
placements on ℤ₂₉, at most 4 of them pairwise compatible (exact maximum clique, §97); BH(7, 6) sign patterns found: one class, J − 2P, with no placement on ℤ₂₉.
