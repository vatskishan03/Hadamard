# Hadamard matrix of order 34196

34196 = 4 · 8549, with 8549 = 83 · 103.

## Construction

Cooper–Wallis: T-matrices of order 103 ([`t-matrices/T103.txt`](../../t-matrices/T103.txt)) give an orthogonal design
OD(412; 103, 103, 103, 103) through the Goethals–Seidel array; substituting four Williamson-type matrices of order 83
for its variables gives a Hadamard matrix of order 4 · 103 · 83 = 34196 (J. Cooper and J. Wallis, Bull. Austral. Math.
Soc. 7 (1972)).

The Williamson-type matrices are bordered ([`williamson-type`](../../williamson-type)): four ±1 sequences x₁, …, x₄ on
ℤ₈₂, the last three symmetric, with PAF₁ + PAF₂ + PAF₃ + PAF₄ = −4 off 0 and row sums 0, 2, 0, 0, all constant on the
orbits of the multiplier group ⟨37⟩ ≤ ℤ₈₂* of order 5, give

    Wₖ = [[1, bₖ eᵀ], [cₖ e, Mₖ]],

with Mₖ the circulant matrices (M₁R for the first) and the border signs of the construction. `certificate.json` holds
the T-matrices and the four sequences; `build_cooper_wallis_bordered` in [`verify/verify.py`](../../verify/verify.py)
rebuilds W₁, …, W₄, checks WₖWₗᵀ = WₗWₖᵀ and Σ WₖWₖᵀ = 332 I, assembles the matrix in int8 and checks H Hᵀ = 34196 I.

## Search

`python3 search/bnw.py search 83 37` lists the four blocks completely (12,870, 70, 70 and 70 members) and joins them
exactly: 320 solutions, the first of which is in the certificate. Under a second.

## Status

**New, beyond the tables.** The order is listed as a candidate in [NOTES.md](../../NOTES.md) §7
("4 · 8549 = 4 · 83 · 103, which needs only T103 and the classical WT(83)"); this matrix is that construction.

No construction we know of reaches 34196:

- Paley I and II need 34195 = 5 · 7 · 977 and 17097 = 3 · 41 · 139 to be prime powers.
- Miyamoto does not apply: 8549 is not a prime power. Neither is 8547 = 3 · 7 · 11 · 37 (Yamada, Kiyasu).
- 8549 is not of Spence's form q² + q + 1.
- A product route needs an object of order 103 or of order 83 · 103. T-matrices of order 103 are from 2026
  (JD-Jones-ASES/Hadamard-T, decoded from Schneider's matrix of order 2060, and found independently in
  [2060](../2060)). The alternative T(83) × WT(103) would need T-sequences or T-matrices of order 83 and
  Williamson-type matrices of order 103, the latter known only from the near-Williamson matrices of 2026.

We found no construction or statement of order 34196. The sources searched are listed in
[NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | T-matrices of order 103, the four sequences on ℤ₈₂, the multiplier, SHA-256 |

The matrix is not stored (2.3 GB as CSV).

## Verification

```
python3 verify/verify.py 34196
```

This rebuilds and checks the matrix in 4 minutes 40 seconds on 4 cores, using about 2.5 GB of memory.

SHA-256 of `H34196.csv`: `f4c7b260ccde989a5e591c7978c4d3431d229e6b8a83290b60572754c5a9ca5a`
