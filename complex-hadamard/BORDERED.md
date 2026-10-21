# The bordered Szöllősi array

A construction of complex (quaternary) Hadamard matrices of order 2N from four ±1 sequences on an abelian group of
order N − 1. To our knowledge it has not been stated before. It relates to Szöllősi's array
(arXiv:2603.09572, Theorem 4) the way Spence's bordered construction relates to the Goethals–Seidel array.

## Statement

Let G be an abelian group of order m = N − 1 and let A, B, C, D : G → {±1} with:

1. A and B symmetric: A(−g) = A(g) and B(−g) = B(g).
2. Σ_k PAF_k(g) = −4 for every g ≠ 0, where PAF_X(g) = Σ_t X(t) X(t + g).

Summing condition 2 over g gives r_A² + r_B² + r_C² + r_D² = 4m − 4(m − 1) = 4 for the row sums. So one row sum is ±2
and the other three are 0.

Put α = (1 + i)/2 and write X_A for the group-developed matrix X_A[g, h] = A(h − g). Let R be the permutation matrix
of g ↦ −g, and e the all-ones vector. Define

    X′ = α X_A + ᾱ X_B,            Y′ = α X_C + ᾱ X_D R,
    X = [[x, s eᵀ], [t e, X′]],     Y = [[y, u eᵀ], [v e, Y′]],     K = [[X, Y], [Y*, −X*]],

with x, s, t, y, u, v ∈ {1, i, −1, −i}. Let ρ_X = α r_A + ᾱ r_B and ρ_Y = α r_C + ᾱ r_D. Then K has entries in
{±1, ±i}, and K K* = 2N I if and only if

    s v = u t,
    x u + s ρ_Y = y s + u ρ_X,
    t y + v ρ_X = v x + t ρ_Y,
    x t̄ + s ρ̄_X + y v̄ + u ρ̄_Y = 0,
    x̄ s + t̄ ρ_X + ȳ u + v̄ ρ_Y = 0.

Each of the eight row-sum patterns (±2 in one of the four places) has exactly 32 solutions (x, s, t, y, u, v). For
example, when r_A = 2:

    x = 1,  s = −t̄,  u = −i y t̄,  v = i t y      (t, y arbitrary).

## Proof

All four entries of X′ are of the form αa + ᾱb with a, b = ±1, and αa + ᾱb ∈ {±1, ±i}. Likewise for Y′. So K is
quaternary.

**Commutation and norms of the cores.**
- Group-developed matrices commute with each other, and R X_Z R = X_Zᵀ for every Z.
- A and B are symmetric, so X_A and X_B are symmetric, X′ is symmetric, and X′R = RX′. Hence X′Y′ = Y′X′.
- X′X′* = (X_A X_Aᵀ + X_B X_Bᵀ)/2, because α² + ᾱ² = 0.
- Y′Y′* = (X_C X_Cᵀ + X_D X_Dᵀ)/2 + (α² X_C X_D R + ᾱ² X_D R X_Cᵀ). The bracket vanishes, since
  X_D R X_Cᵀ = X_D X_C R = X_C X_D R.
- So X′X′* + Y′Y′* = (Σ_k X_k X_kᵀ)/2 = 2N I − 2J by condition 2, where J is the all-ones matrix.
- All of these matrices are normal, so the same identities hold for X′*X′ and Y′*Y′.

**The border.** e is an eigenvector of X′, Y′ and of their adjoints, with eigenvalues ρ_X, ρ_Y (and their conjugates):
- e is fixed by R;
- the row and column sums of X_Z are r_Z.

Multiplying out the blocks of K K* = [[XX* + YY*, XY − YX], [Y*X* − X*Y*, Y*Y + X*X]]:
- The corner entries are 2(1 + m) = 2N.
- The core blocks are J + J + X′X′* + Y′Y′* = 2N I.
- The border rows and columns give the five scalar equations above: two from XX* + YY* and X*X + Y*Y, and three
  from XY = YX (the J-coefficient tu = vs of XY − YX repeats the first equation).

∎

## Verification

`python3 search/bsz.py test` runs two checks:
- It enumerates every solution of conditions 1–2 on ℤ₆ (N = 7). All 2304 quadruples give complex Hadamard matrices
  of order 14 (exact integer arithmetic).
- It lists the 32 border solutions of each row-sum pattern.

`search/bsz_group.py` builds the array over ℤ₂² × ℤ_m. For N = 13, 21 and 25 it finds and verifies complex
Hadamard matrices of orders 26, 42 and 50. At N = 13 and 21 the solution counts, 12672 and 4.6 · 10⁷, agree with
the Gaussian estimates 10^4.25 and 10^7.02. `verify/verify_complex.py` accepts certificates with constructions
`szollosi-bordered` (G = ℤ_m) and `szollosi-bordered-group`.

## Results

**Order 454 (N = 227).** 227 = 2 · 113 + 1 is a safe prime: besides {1}, the only odd-order subgroup of ℤ₂₂₇* is the
group of squares, with three orbits, so Szöllősi's array over ℤ₂₂₇ has no usable structure. Over ℤ₂₂₆ = ℤ₂ × ℤ₁₁₃
the multiplier 49 has order 7 (trivial on ℤ₂, order 7 on ℤ₁₁₃) and 34 orbits.
- Rows (0, 0, 2, 0): A and B share one complete list of 25740 members; C (row sum 2) and D (row sum 0) have pools of
  6.0 · 10⁸ and 1.2 · 10⁹, of which the low-energy k-tree (`search/lowe_sz.py`) kept 3 · 10⁷ and 2.77 · 10⁷.
- `kt4` found 512 solutions at the 851st split value of the solution law. The per-class Gaussian estimate of the
  structure was 10^0.63 (10^−1.06 inside the kept lists); the hybrid estimate of `search/gauss_mc.py` was 10^4.4.
- Border: x = s = t = y = 1, u = v = −i. Certificate [C454.json](C454.json) (construction `szollosi-bordered`).

**Order 366 (N = 183).** 183 = 3 · 61. Over ℤ₁₈₂ = ℤ₂ × ℤ₇ × ℤ₁₃ the multiplier 3 has order 6 (trivial on ℤ₂, order
6 on ℤ₇, order 3 on ℤ₁₃) and 36 orbits.
- Rows (0, 0, 2, 0). The pools of A and B have 55968 members each; those of C (row sum 2) and D have 1.28 · 10⁹ and
  2.57 · 10⁹.
- Compressions to ℤ₂, ℤ₇ and ℤ₁₃ that occur in no level-1 solution rule out 93% of A and B, 97.9% of C and 96.4% of
  D (fibre filters, `search/lowe_sz.py`, [search/NOTES.md](../search/NOTES.md) §55). That left 3880 members in A and
  B, all 2.68 · 10⁷ allowed members of C, and 3 · 10⁷ of the 9.35 · 10⁷ allowed members of D.
- `kt4` found 1536 solutions at the 146th split value. Without the filters, with lists kept to 2.3% and 1.2% of C and
  D, it processed about 4700 split values without a solution.
- The per-class Gaussian estimate of the structure was 10^−0.1; the hybrid estimate was 10^3.2.
- Border: x = s = t = y = 1, u = v = −i.
- **Certificate.** [C366.json](C366.json) is the cyclic quadruple over ℤ₁₈₃ (⟨34⟩, row sums 1, 23; 9, 11) of
  [search/NOTES.md](../search/NOTES.md) §42, checked by `verify/verify_complex.py`. The bordered solutions found here
  are not stored.

## Context

Constructions from four ±1 sequences, in two families:

| symmetric blocks | on ℤ_N, ΣPAF = 0 | on a group of order N − 1, ΣPAF = −4 |
|---|---|---|
| none | Goethals–Seidel: H(4N) | Spence: H(4N) |
| two | Szöllősi: complex Hadamard of order 2N | **bordered Szöllősi: complex Hadamard of order 2N** |
| three | near-Williamson (Kharaghani, Mohammadian, Tayfeh-Rezaie, Zaitsev 2026): Williamson-type of order N | bordered Williamson-type ([williamson-type](../williamson-type)): Williamson-type of order N |

The bordered versions move the search to a group of order N − 1. That helps when ℤ_N has no useful multiplier group,
for example the safe primes N = 2q + 1. Order 454 is of this kind: N = 227 = 2 · 113 + 1, on ℤ₂₂₆ with H = ⟨49⟩ of
order 7, the complex Hadamard matrix of order 454 ([C454.json](C454.json), rows 0, 0; 2, 0; border
(1, 1, 1, 1, −i, −i)). The order itself is known (4 · 227 − 1 = 907 is a prime; [README.md](README.md#status)).

**Negative results** (NOTES.md §9). A 3-point border with border rows and columns constant on the cosets of an
index-4 subgroup is impossible for Williamson-type matrices: algebraically for character borders, and by SAT at
p = 11. For the complex array, the character ansatz is also contradictory: the off-diagonal border equations force
ρ_i = −ρ_j for all three pairs i ≠ j. The general coset-constant case at N = 7 was left undecided by a 15-minute SAT
run.
