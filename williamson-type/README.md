# Bordered Williamson-type matrices

Williamson-type matrices of order w are four ±1 matrices W₁, …, W₄ of order w with WₖWₗᵀ = WₗWₖᵀ (k ≠ l) and Σ WₖWₖᵀ =
4w I. With T-matrices of order t they give a Hadamard matrix of order 4tw (Cooper–Wallis), and with t = 1 the
Williamson array gives order 4w. This directory holds Williamson-type matrices built from sequences on the group
ℤ_(w−1), with one border row and column. The point is the group: when w is prime, the multipliers of ℤ_w may be
useless (ℤ₁₇₉* and ℤ₂₆₃* have no small odd subgroups), while ℤ_(w−1) has many.

## Construction

Let m = w − 1 and let x₁, x₂, x₃, x₄ be ±1 sequences on ℤ_m with x₂, x₃, x₄ symmetric (x(−g) = x(g)) and

    PAF₁(g) + PAF₂(g) + PAF₃(g) + PAF₄(g) = −4      for every g ≠ 0.

Summing over g gives Σ rₖ² = 4 for the row sums, so rⱼ = ±2 for one index j and rₖ = 0 otherwise. Let Mₖ be the
circulant matrix of xₖ, except M₁ R in place of M₁ (R the back-diagonal permutation i ↦ m − 1 − i, so M₁R is
symmetric; the permutation g ↦ −g works as well and gives another matrix, but the certificates use the
back-diagonal), let e be the all-ones vector of length m and s = rⱼ/2. Put

    Wₖ = [ 1     bₖ eᵀ ]        cₖ = 1, bₖ = s (k ≠ j),      cⱼ = −1, bⱼ = −s.
         [ cₖ e   Mₖ   ]

**Proof.** WₖWₗᵀ has top-left entry 1 + m bₖbₗ, top-right block (cₗ + bₖrₗ) eᵀ, bottom-left block (cₖ + rₖbₗ) e and
bottom-right block cₖcₗ J + MₖMₗᵀ.
- Amicability: M₁R Mₗᵀ = M₁MₗR = MₗM₁R = Mₗ(M₁R)ᵀ for every circulant Mₗ (R Mₗᵀ R = Mₗ), and MₖMₗᵀ = MₗMₖᵀ for
  symmetric circulants. The border blocks need cₗ + bₖrₗ = cₖ + bₗrₖ, which holds for the choice above (both sides are
  1 when j ∉ {k, l}; when k = j both sides are 1 = −1 + 2s²).
- Sum: Σ (1 + m bₖ²) = 4w; Σ (cₖ + bₖrₖ) = 3 + (−1 − 2s²) = 0; and Σ (cₖ² J + MₖMₖᵀ) = 4J + (4m I − 4(J − I)) = 4w I,
  since Σ MₖMₖᵀ has 4m on the diagonal and −4 elsewhere.

Over ℤ_m = ℤ₂ × ℤ_(m/2) (m/2 odd) this is the "Spence-W" structure of [NOTES.md](../NOTES.md) §7: Spence's condition
for a Hadamard matrix of order 4w with three symmetric sequences.

**Search** (`search/bnw.py`). All four sequences are constant on the orbits of a multiplier group H ≤ ℤ_m*, the
symmetric ones on the orbits of ⟨H, −1⟩; the condition then holds on the ⟨H, −1⟩-orbits of ℤ_m \ {0}. Each block is
listed completely in Gray-code order with incremental class values (`search/bnw_enum.c`), the lists are joined exactly
(`search/join4x.c`) or by the k-tree (`search/kt4.c`), and every hit is rebuilt as W₁, …, W₄, checked, and checked
again through the Williamson array. `verify/verify.py` builds Hadamard matrices from these quadruples with the
construction `cooper-wallis-bordered`.

## Matrices

| w | group | H | row sums | file | used for | Status |
|---:|---|---|---|---|---|---|
| 47 | ℤ₄₆ | ⟨45⟩ = {±1}, all four symmetric | 2, 0, 0, 0 | [WT47.json](WT47.json) | [12596](../hadamard/12596) = 4 · 67 · 47, [19364](../hadamard/19364) = 4 · 103 · 47 | Known |
| 83 | ℤ₈₂ | ⟨37⟩, order 5 | 0, 2, 0, 0 | [WT83.json](WT83.json) | [34196](../hadamard/34196) = 4 · 103 · 83 | Known |

Williamson-type matrices of orders 47 and 83 were known before (near-Williamson matrices of 2026 for 47;
Seberry–Yamada 1992, Table A.1 for 83), so both are Known; the matrices here are instances of the bordered form. Order
47 came from uniform samples of 4 · 10⁵ of the 1.35–1.41 · 10⁶ members of each block, matched by the k-tree (4
solutions in 288 s). The structure holds 320 solutions (exact join), while the Gaussian model of `bnw.py triage`
predicts 10^−2.1.

## Calibration and limits

Exact counts against the Gaussian model on all structures with at most 17 and 15 orbits for 27 ≤ w ≤ 133: the model
was within one to two orders of magnitude where solutions exist (w = 27: 10^2.2 predicted, 10^3.4 found), and every
structure it rated below 10^−2 was empty, with the exception above (w = 83, ⟨37⟩).

For the primes that would settle open Hadamard orders through T-matrices of small order (Cooper–Wallis with t = 3, 5,
7, …), the counts are against this route:

| w | settles (with T(t)) | log₁₀ E, no symmetry | best with \|H\| ≥ 2 |
|---:|---|---:|---:|
| 167 | 1169 (7), 2171 (13) | 10.7 | −4.4 |
| 179 | 1969 (11), 2327 (13) | 10.4 | −5.1 |
| 191 | 573 (3), 955 (5), 1719 (9), 2865 (15) | 10.1 | −4.1 |
| 223 | 669 (3), 1115 (5), 2899 (13) | 8.3 | −4.0 |
| 263 | 789 (3), 1315 (5), 1841 (7), 2893 (11) | 5.7 | −4.2 |

The four orders of 263 no longer need Williamson-type matrices: Paley's matrix of order 1052 in sixteen circulant
blocks gives them with the same T-matrices ([../hadamard/PALEY_BLOCKS.md](../hadamard/PALEY_BLOCKS.md)).

Solutions exist in quantity without symmetry, in a space of about 2⁴⁷⁸ sign choices for w = 191, but every multiplier
group removes them. For w = 191 the structures with |H| = 24, 12 (two) were searched exhaustively and are empty.
