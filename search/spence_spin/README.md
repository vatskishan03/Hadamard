# Spence spin search (order 2292)

The programs that found the Hadamard matrix of order 2292 = 4 · 573 ([2292](../../2292)).
Structure: Spence's construction over G = ℤ₂ × ℤ₂ × ℤ₁₄₃, with blocks
(A, B, B∘μ, B∘μ²), H = ⟨27⟩ and μ = (α, ×3). See [NOTES.md](../../NOTES.md) §8.

Requirements: Python 3 with NumPy, and a C compiler.

| File | Role |
|---|---|
| `ss_core.py` | groups, orbits, classes and class keys for Spence spin structures |
| `struct573.py` | the structure for 573 |
| `levels.py`, `exact_small.py`, `lattice.py` | quotient levels G/N, exact gates (product value sets) and the lattice test |
| `hashjoin.py`, `stage1.py`, `fillB.c` | solutions at the level ℤ₂ × ℤ₂ × ℤ₁₃ (energy-matched search) |
| `liftlab.py`, `lift_poly.py`, `mklift.py`, `liftsolve.c` | the lift to G: Legendre fibres, the polynomial form of the 27 class keys, the problem file, the exhaustive solver |
| `verify_lift.py` | rebuilds the four sequences and the Hadamard matrix from a lift solution |
| `wit_s11.npy`, `wit_s12.npy` | the two level solutions that lifted |

Reproduce the matrix in the repository from the saved level solution (about a minute):

```
cc -O2 -o liftsolve liftsolve.c
python3 mklift.py wit_s11.npy p11.txt
./liftsolve p11.txt 12              # prints: SOL 41286 126925355812639
python3 verify_lift.py wit_s11.npy 41286 126925355812639 s11
```

The last command checks the Spence condition at every shift and H Hᵀ = 2292 I. It
writes `certificate_s11.json` and prints the SHA-256 recorded in
[2292/README.md](../../2292/README.md). `wit_s12.npy` gives a second, different
matrix the same way.

New level solutions: `cc -O2 -o fillB fillB.c`, then `python3 stage1.py <seed> 12 120 2 1e6`.
Each run writes `wit_u3_s<seed>_<k>.npy`.

Completeness check of the lift solver: `python3 mklift.py wit_s11.npy plant.txt 7` replaces
the target by the class keys of random signs and prints them as `PLANTED`. `./liftsolve plant.txt 12`
must return them.
