# The octic ring: Goethals–Seidel arrays of order 16v

ℤ[x]/(x^{4v} + 1) ≅ ℤ[ζ₈][ℤ_v] for v odd: a ±1 sequence of length 4v becomes a sequence on ℤ_v over 16 letters c₀ +
c₁ζ + c₂ζ² + c₃ζ³ (c_j = ±1) of squared moduli 4 ± 2√2, and the negacyclic Goethals–Seidel condition is Σ U_X U_X* = 16v
δ, with a peak equation that has no analogue in the order-8v family ([../../NOTES.md](../../NOTES.md) §89). A
palindromic first block is a skew block: U(0) = (1, c, d, c) and U(−g) = −conj U(g).

## Contents

| File | Contents |
|---|---|
| `octic.py` | exact arithmetic of the ring: ℤ[w], the isomorphism φ (sequence → U) and its inverse ψ, the adjoint, the alphabet |
| `verify_all.py` | every numerical check of the derivations, one PASS line each (ring isomorphism, alphabet, peak equation, skew ⇔ palindromic, row sums, multiplier classification, θ-pairs) |
| `mult_classify.py` | the symmetries a(x) = ε x^s a(x^t) that admit invariant skew sequences: only t ≡ 1 (mod 8), unshifted, −1 ∉ H |
| `exsmall.c` | EXHAUSTIVE counts of the ordered skew quadruples at v = 1, 3, 5 (no multiplier group; meet in the middle over the monomial orbits) and the census of the θ-paired ("doubled") solutions |
| `dbl.py`, `cert_exh.py` | the θ-pairing test of a solution; the matrices of orders 16, 48, 80 from the exhaustive solutions |
| `oskew.c`, `odrv.py` | the orbit-level exhaustive join for a multiplier group H (row-sum classes in ℤ[√2], hash join); the driver |
| `cert_oct.py`, `mkcert.py` | an orbit-level solution to ±1 sequences by ψ, the checks, and the certificate built with `../../../verify/verify.py` (`out/o<16v>/`) |
| `targets.py`, `shapes_oct.py`, `dcfeas.py`, `dcpat.py` | the 61 Table 5 targets (m ≥ 5, not improved), their multiplier shapes, the row-sum feasibility of each shape (\|H\| = h needs (2h − 1)² + 15 ≤ 16v), and the row-sum patterns of one shape |
| `oest.py` | a Gaussian and lattice estimate of the number of solutions of a shape (HEURISTIC) |

## Build

In this directory:

```
gcc -O2 -o oskew oskew.c
gcc -O2 -o exsmall exsmall.c
mkdir -p shapes out runs
```

`oskew` needs the parent directory only for `../qskew` (the order-8v count at the same shape).

## Results

No order 16v was found. The matrices built from exhaustive or orbit-level solutions have known orders and serve as
controls: orders 16, 48 and 80 at v = 1, 3, 5, and orders 112, 208, 464, 592 and 976 at v = 7, 13, 29, 37 and 61. The
last one comes from |H| = 15 at v = 61, a shape without solutions in the order-8v family.

## Reproduction

The commands that reproduce the octic results of [../../NOTES.md](../../NOTES.md) §89 are in the section
[Reproduction](../README.md#reproduction) of [../README.md](../README.md).
