# Hadamard matrix of order 2292

2292 = 4 · 573, with 573 = 3 · 191.

## Construction

Spence's bordered construction over the group

```
G = ℤ₂ × ℤ₂ × ℤ₁₄₃,   |G| = 572 = 573 − 1   (ℤ₁₄₃ ≅ ℤ₁₁ × ℤ₁₃).
```

Four ±1 sequences a₀, …, a₃ on G with row sums 2, 0, 0, 0 satisfy

```
Σᵢ PAF_aᵢ(d) = −4   for every d ≠ 0.
```

The Goethals–Seidel array developed over G has order 4 · 572 = 2288, and a border of four rows and columns built from
two 4 × 4 Hadamard matrices P and Q completes it to a Hadamard matrix of order 2292 (`build_spence` in
`verify/verify.py`).

**Spin structure.** Write H for the multiplier subgroup ⟨27⟩ = {1, 14, 27, 53, 92} of ℤ₁₄₃* (it acts as the quadratic
residues on ℤ₁₁ and trivially on ℤ₁₃). Let μ = (α, ×3), where α(x₁, x₂) = (x₂, x₁ + x₂) is the automorphism of order 3
of ℤ₂ × ℤ₂ and ×3 is multiplication on ℤ₁₄₃. Then μ has order 15, μ³ = ×27 lies in H, and K = ⟨H, μ⟩ = ⟨μ⟩. The four
sequences are

```
a₀ = A,   a₁ = B,   a₂ = B ∘ μ,   a₃ = B ∘ μ²,
```

with A constant on the 54 orbits of K and B constant on the 156 orbits of H. The condition then holds on the orbits of
⟨K, −1⟩, which gives 27 equations.

## Search

1. **Exact quotient gates.** For each subgroup N invariant under H and μ, the compressions to G/N (sums over the cosets
   of N) satisfy the same kind of equations. Each orbit lies over one orbit of G/N and meets every fibre in the same
   number of points, so the attainable compressed vectors are products of per-orbit value sets. The structure passes all
   six quotients: ℤ₂², ℤ₁₁, ℤ₁₃, ℤ₂² × ℤ₁₁, ℤ₁₄₃ and ℤ₂² × ℤ₁₃. It also passes a lattice test on the class values.
2. **Level ℤ₂² × ℤ₁₃** (52 points, fibre ℤ₁₁). Solutions at this level come from a hash join at level ℤ₁₃, followed by
   an energy-matched search. The A candidates are drawn at a fixed energy, and the B columns are restricted to the
   complementary energies. This takes about one verified solution per 20–60 CPU seconds.
3. **Lift to G.** Over a compressed value ±1, the fibre ℤ₁₁ carries a "Legendre" pattern (b₀ at 0, s·χ(ρ) elsewhere, χ
   the quadratic character) with one free sign s. Compressed values ±9 and ±11 force the fibre completely. That leaves
   61–63 free signs. The 27 equations then split into 8 linear equations and quadratic ones.
   - An exact four-list meet-in-the-middle search lists every sign vector that satisfies the linear equations, about 10⁷
     per lift.
   - Each of these candidates is then checked against all 27 equations.
   - The solver was checked first on planted solutions, and recovered every one.
4. Of the first five solutions at level ℤ₂² × ℤ₁₃, two lifted to full solutions. This matrix comes from the first. The
   second is a different matrix of the same order.

## Status

**New.** The entry 573(3) is unresolved in Cati–Pasechnik Table 4.

- Paley's constructions do not apply: 4 · 573 − 1 = 2291 = 29 · 79 and 2 · 573 − 1 = 1145 = 5 · 229.
- Miyamoto's theorem does not apply: a symmetric conference matrix of order 574 cannot exist (573 = 3 · 191 is not a sum
  of two squares).
- Yamada's theorem does not apply: 573 − 2 = 571 ≡ 3 (mod 8).
- The orders announced in August 2026 stop below 2000.
- After 2060, 2092, 2148 and 2284 were settled, 2292 was the smallest order we knew to be unresolved.

The sources searched are listed in [NOVELTY.md](../../NOVELTY.md). The matrix also settles 9172 = 4 · 2293 by Miyamoto's
theorem ([9172](../9172)).

## Files

| File | Contents |
|---|---|
| `certificate.json` | group, the four sequences (`+` = 1, `-` = −1, elements in mixed-radix order of the group), P, Q, SHA-256 |
| `sequences.txt` | the four sequences, one per line |
| `H2292.csv.gz` | the matrix |

## Verification

```
python3 verify/verify.py 2292
gzip -dc hadamard/2292/H2292.csv.gz | ./check_matrix
```

SHA-256 of `H2292.csv`: `99acfee8855745a3e696ac1cc3dd8dd17ee479e17c8adbe68bfb2671978d9142`