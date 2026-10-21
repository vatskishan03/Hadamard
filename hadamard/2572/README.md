# Hadamard matrix of order 2572

2572 = 4 · 643, with 643 prime.

## Construction

The same three steps as [6796](../6796), with the notation defined there, each checked exactly by `verify/verify.py`.

**1. A quaternary conference matrix of order 322.** The zero-diagonal Szöllősi array
([search/NOTES.md](../../search/NOTES.md) §74–76) over the cyclic group ℤ₁₆₁ = F₇ × ℤ₂₃. A and B map ℤ₁₆₁ to
{0, ±1}, are symmetric and vanish exactly at 0; C and D map it to {±1}; and

```
PAF_A(g) + PAF_B(g) + PAF_C(g) + PAF_D(g) = 0   for g ≠ 0.
```

With ω = (1 + i)/2, X_f[g, h] = f(h − g) and the reversal R[g, −g] = 1, the matrices X = ω X_A + ω̄ X_B and
Y = ω X_C + ω̄ X_D R give Q = [[X, Y], [Y*, −X*]] of order 322 with zero diagonal, entries ±1, ±i elsewhere, and
Q Q* = 321 I.

The sequences are constant on the orbits of H = ⟨93⟩, of order 3: 93 ≡ 2 (mod 7) and 93 ≡ 1 (mod 23). So H is the
group of squares {1, 2, 4} of F₇*, acting trivially on ℤ₂₃. C and D are constant on the 69 orbits of H, and A and B
on the 35 orbits of ⟨H, −1⟩. The row sums are −20, −4, 15 and −1: 400 + 16 + 225 + 1 = 642 = 4 · 161 − 2.

This is a *Paley plug-in* ([search/NOTES.md](../../search/NOTES.md) §81). On ℤ₂₃ × {0} each block is a function x₀,
and on ℤ₂₃ × F₇* it is x₁ on the squares and x₂ on the non-squares. Off 14 *defect* positions, x₂ = −x₁, so the block
there is x₀ ⊗ δ + x₁ ⊗ χ with χ the Legendre symbol of F₇. The defects are at k = 0, 10, 13 (A), 0, 4, 19 (B),
9, 10, 15, 19 (C) and 0, 12, 19, 20 (D), where k is the ℤ₂₃ coordinate.

**2. A quaternary Hadamard matrix of order 1286.** q = 2 · 322 − 3 = 641 is a prime ≡ 1 (mod 8), and the complex
Yamada–Kiyasu template ([NOTES.md](../../NOTES.md) §31, §34) gives K of order 2(q + 2) = 1286 with K K* = 1286 I
([complex-hadamard/C1286.json](../../complex-hadamard/C1286.json)).

**3. Turyn's doubling** H = [[K_r + K_i, K_r − K_i], [−(K_r − K_i), K_r + K_i]] gives order 2572.

## Search

The core was found by a two-level sampler on the Paley plug-in structure, in the manner of the FLC classification
method ([17] in the [README](../../README.md#references)). A generic two-level sampler on the same group and
multiplier group had found nothing: about 10⁻⁷ completions per instance
([search/NOTES.md](../../search/NOTES.md) §78).

The class came from a calibration census, in which the plug-in class ℤ₃² × ℤ₁₇ of order 153 has tens of thousands of
solutions against 5 predicted, all with the defect structure above ([search/NOTES.md](../../search/NOTES.md) §81). For
ℤ₂₃ × F₇ the identity of the trivial character forces 2n_p + n_m = 20, where n_p and n_m count the defects with
x₀ = x₁ and with x₀ = −x₁ (this core: n_p = 8, n_m = 4); two further defects sit at the zeros of A and B. So a
solution has 12 to 22 defects. The class was estimated to hold about 10¹² ordered solutions, the only target class
with a large expected count ([search/NOTES.md](../../search/NOTES.md) §81–82).
- **Level 1** (`search/plug7_gen.c`) samples x₀ quadruples with Σ PAF_{x₀} ≡ 0 (mod 3). It places the defects in the
  allocation (3, 3, 4, 4) and joins the four blocks on the trivial-F₇-character equation by a class-restricted hash
  join; two row-sum identities are exact filters.
- **Level 2** (`search/plug7_l2.c`) chooses the signs t = (x₁ − x₂)/2 off the defects by a meet-in-the-middle join on
  the remaining equations. `search/plug7_sol.py` turns a solution into the four sequences.

Calibration at the analogues ℤ₁₃ × F₇ and ℤ₁₇ × F₇ predicted about 12 CPU-hours per core in the row-sum class
τ_C² + τ_D² = 10. The core was found after 42451 instances of that class, about 13.5 CPU-hours on four cores, plus 2
CPU-hours of earlier test runs ([search/NOTES.md](../../search/NOTES.md) §82–83).

**Reproduction.**

| | |
|---|---|
| Code | `search/plug7_gen.c` (level 1, `plug7_gen p a b c d nquads seed logWthr outfile`), `search/plug7_l2.c` (level 2, `plug7_l2 file [maxsol_print] [tlimit_per_instance]`; x86-64 with AVX2), `search/plug7_sol.py` (turns a solution into the four sequences) |
| Instances | the level-1 output feeds level 2 in the instance format of the sampler of [search/NOTES.md](../../search/NOTES.md) §78; the converter between them is not included |
| Cost | about 12 CPU-hours per core in the row-sum class τ_C² + τ_D² = 10 (calibrated; range 4 to 30); the search was not rerun |
| Check | `python3 verify/verify.py 2572` gives `ok` and SHA-256 `ebf232c3f7a8b99e47c25c004cd3a6fd759710c7ce88b856ee4931ac5ed6f8b6` |

The certificate alone proves the result.

## Status

**New.** The entry 643(3) is unresolved in Cati–Pasechnik Table 4 (v1 and v2), and Table 5 of Egan, Kotsireas, London
and Ó Catháin (arXiv:2610.27175, 19 October 2026) lists 643 as open. No later construction was found. Table 4 already
knows 8 · 643, so doubling settles nothing further. Đoković's set Δ and his conversions of 2010 do not contain 643,
and SageMath's development version (21 October 2026) has no branch that reaches 2572. None of the standard recursive
constructions applies:

- Paley: 2571 = 3 · 857 and 1285 = 5 · 257 are not prime powers, and no symmetric conference matrix of order 1286 is
  known (Brouwer's table lists srg(1285, 642, 320, 321) as undecided).
- Miyamoto: 643 ≡ 3 (mod 4).
- Yamada: it needs a symmetric conference matrix of order 322, which cannot exist, because 321 = 3 · 107 is not a sum
  of two squares.
- Kiyasu: q = 641 ≡ 1 (mod 8), not 5.
- Product theorems: 2572 ≡ 4 (mod 8).
- Spence, Mathon, Seberry–Whiteman: 1 + q + q² = 643, q²(q + 2) = 1285 and 5 · 9^k = 1285 have no solutions.
- Cooper–Wallis: 643 is prime, so it needs T-matrices or Williamson-type matrices of order 643. Neither is known.
- The block routes of [PALEY_BLOCKS.md](../PALEY_BLOCKS.md) and [MIYAMOTO_BLOCKS.md](../MIYAMOTO_BLOCKS.md): they need
  2571 a prime power or 643 ≡ 1 (mod 4).
- Legendre pairs and cocyclic matrices: no binary Legendre pair of length 1285 is known, and a cocyclic H(4p) with
  p ≡ 3 (mod 4) is of Williamson or Ito type, neither of which exists for 643.

**Credits.** As for [6796](../6796/README.md#status). Miyamoto's Corollary 4 (1991, as transcribed in the
JD-Jones-ASES Hadamard-M note) reaches the same orders 4(q + 2) from a C₂-matrix of order q + 3. In his displayed
form, such a matrix is the same thing as a quaternary conference matrix Q of order (q + 3)/2 with Qᵀ = Q. The template
needs no such symmetry, and the core here has Qᵀ ≠ Q as built; it is not Hermitian either. The sampler adapts the FLC
classification method (fibre-lattice completion, [17] in the [README](../../README.md#references)), an exhaustive
classifier, to sampling.

**Ingredients.**
- The quaternary conference matrix of order 322 is new to our knowledge. No real one exists: a symmetric one would
  need 321 to be a sum of two squares, and a skew one needs order ≡ 0 (mod 4). The families of complex ones
  (Paley-type, Berman, BGW and GC, Seberry–Whiteman, Et-Taoui) need 321 to be a prime power. Fallon and Iverson's
  ancillary library (arXiv:2610.09169) has a conference signal over ℤ₃₂₁, which gives a complex conference matrix of
  order 322. That matrix is not equivalent to a quaternary one: some of its invariants Q_ij Q_kl Q̄_il Q̄_kj lie 0.77
  from {±1, ±i}, the largest distance possible.
- The quaternary Hadamard matrix of order 1286 settles an order that is blank in Seberry–Yamada's Table 11.2 (1992)
  and on the unknown list of Zaitsev et al.'s page (last updated 14 October 2026, unchanged on 21 October). Its other
  routes fail: 1285 and 2571 are not prime powers (Turyn, Ito–Egan), and no Williamson or Williamson-type matrices of
  order 643 are known.

**Risk.** Zaitsev et al. use the same template and do not date their matrices; their page lists 1286 as unknown.
Middletownbooks/Hadamard lists 643 as open (last commit 20 October, README only) and records a screen of
Goethals–Seidel and Legendre-quadruple structures over ℤ₆₄₃ (6 October 2026) without a hit. The Seberry–Yamada book
of 2020, the Handbook's Table 1.53 row for 643 and Miyamoto's lists of 1991 were not read; his Corollary 4 at m = 321
would need a symmetric quaternary conference matrix of order 322. The sources searched are listed in
[NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | the group, the reversal, the symmetry group, the four sequences, and q, r, ξ for Yamada's circulants |

The matrix is not stored. `python3 verify/verify.py --write 2572` writes `H2572.csv`.

## Verification

```
python3 verify/verify.py 2572
python3 verify/verify_complex.py 1286
```

SHA-256 of `H2572.csv`: `ebf232c3f7a8b99e47c25c004cd3a6fd759710c7ce88b856ee4931ac5ed6f8b6`

An independent rebuild from the certificate alone (its own code for the core, Yamada's circulants, K and H; program
not included) writes a byte-identical `H2572.csv` and checks H Hᵀ = 2572 I exactly.
