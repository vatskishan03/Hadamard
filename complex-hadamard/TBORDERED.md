# The character-bordered Szöllősi array

A construction of complex (quaternary) Hadamard matrices of order 2N from four ±1 sequences on an abelian group of
order N − t and a border of size 2t. For t = 1 it is the bordered array of [BORDERED.md](BORDERED.md). It is the
complex analogue of the bordered Goethals–Seidel arrays with a character set T
([search/NOTES.md](../search/NOTES.md) §18), the family of the Hadamard matrix of order 3964 and of Alpöge et al.
To our knowledge it has not been stated before.

## Statement

Let G be an abelian group of order m, φ: G → Q a homomorphism onto an abelian group of order q with kernel K, and T
a set of t characters of Q that is a union of Galois orbits (so that Γ(a) = 4 Σ_{χ∈T} χ(a) is an integer function on
Q; in particular T is closed under complex conjugation). Put N = m + t. Let A, B, C, D : G → {±1} with

1. A and B symmetric: A(−g) = A(g), B(−g) = B(g);
2. Σ_k PAF_k(g) = −Γ(φ(g)) for every g ≠ 0.

Write y_X(a) = Σ_{φ(g)=a} X(g) for the compressions, X_Z for group-developed matrices (X_Z[g, h] = Z(h − g)), R for
the matrix of g ↦ −g, α = (1 + i)/2, and

    Kc = [[X, Y], [Y*, −X*]],   X = α X_A + ᾱ X_B,   Y = α X_C + ᾱ X_D R          (Szöllősi's array over G),

and Mc for the same array over Q built from the compressions. Let V_T be the space of pairs (f₁, f₂) of functions on
Q whose Fourier transforms vanish off T (dimension 2t). Suppose there are quaternary matrices (entries ±1, ±i)

- R̂ of size 2t × 2q whose rows lie in V_T and are mutually orthogonal,
- Ĉ of size 2q × 2t whose columns lie in V_T and are mutually orthogonal,
- Xb = −R̂ Mc* Ĉ / (2q) of size 2t × 2t.

Put R[a, (i, g)] = R̂[a, (i, φ(g))] and R′[(i, g), b] = Ĉ[(i, φ(g)), b]. Then

    K = [[Xb, R], [R′, Kc]]

is a complex Hadamard matrix of order 2N.

## Proof

In Fourier terms condition 2 says Σ_k |X̂_k(ψ)|² = 4N for every character ψ of G that is not a lift of a member of
T, and 4t for the lifts of T (condition 2 at g = 0 holds automatically: Σ_k PAF_k(0) = 4m = 4N − Γ(0)). On Q this
gives Σ_k |ŷ_k(χ)|² = 4t for χ ∈ T.

- **Core rows.** Group-developed matrices commute, R X_Z R = X_Zᵀ, and X is symmetric, so X and Y commute and are
  normal, and Kc Kc* = I₂ ⊗ M with M = X X* + Y Y* = (Σ_k X_k X_kᵀ)/2 (as in [BORDERED.md](BORDERED.md)). M is
  group-developed with eigenvalue 2N off the lifts of T and 2t on them, so 2N I − M = 2m P, P the orthogonal
  projection onto the lifts of T: P[g, h] = (1/m) Σ_{χ∈T} χ(φ(g) − φ(h)).
- **Border columns.** The columns of Ĉ are an orthogonal basis of V_T, each of norm 2q, so Ĉ Ĉ* = 2q P_T with P_T
  the projection onto V_T, and R′ R′*[(i, g), (j, h)] = (Ĉ Ĉ*)[(i, φg), (j, φh)] = 2 δ_ij Σ_{χ∈T} χ(φ(g) − φ(h)).
  Hence R′ R′* = I₂ ⊗ 2m P and Kc Kc* + R′ R′* = 2N I.
- **Border rows.** R R* = |K| R̂ R̂* = (m/q) · 2q I = 2m I. Mc commutes with P_T (P_T is I₂ ⊗ a symmetric
  Q-developed matrix) and Mc* Mc = I₂ ⊗ (Σ_k Y_k Y_kᵀ)/2, which is 2t on V_T. So Xb Xb* = R̂ Mc* Ĉ Ĉ* Mc R̂*/(4q²) =
  R̂ Mc* Mc P_T R̂*/(2q) = (t/q) R̂ R̂* = 2t I (the rows of R̂ lie in V_T), and Xb Xb* + R R* = 2N I.
- **Cross terms.** Summing a group-developed matrix over a fibre of φ gives the Q-developed matrix of the
  compression (also for X_D R and the conjugate transposes), so R Kc* = (R̂ Mc*) ∘ φ and Xb R′* = (Xb Ĉ*) ∘ φ. With
  Xb = −R̂ Mc* Ĉ/(2q): Xb Ĉ* = −R̂ Mc* P_T = −R̂ P_T Mc* = −R̂ Mc*. So Xb R′* + R Kc* = 0.

So K K* = 2N I; K is square, hence also K* K = 2N I. ∎

## Instances

- t = 1, Q trivial, T = {1}: the bordered array, ΣPAF = −4 off 0 ([BORDERED.md](BORDERED.md)).
- Q = ℤ_q, T = all characters (t = q): ΣPAF = −4q on K \ 0 and 0 off K. Needs q | m, i.e. q | N.
- **Q = ℤ₄, T = {1, χ, χ̄} (t = 3)**: Γ = 4 · (3, 1, −1, 1), so ΣPAF = −12 on K \ 0, −4 on the two cosets
  φ = ±1 and +4 on the coset φ = 2. Needs 4 | N − 3: every N ≡ 3 (mod 4). In particular it applies to the primes N for
  which neither ℤ_N nor ℤ_(N−1) has a usable multiplier group (107, 167, 179: orders 214, 334, 358).
- Q = ℤ₂², T = the three nontrivial characters, or {1, χ₁, χ₂}: also t = 3 and N ≡ 3 (mod 4), on groups with a
  ℤ₂² quotient.

For the ℤ₄ family the compressions satisfy Σ_k r_k² = 12 (r_k the row sums), Σ_k |ŷ_k(χ)|² = 12, and
Σ_k ŷ_k(χ²)² = 4N: three ±2 row sums and one 0 when |K| is even. V_T has 36 quaternary vectors per half, so the border
problem has 1296 candidate rows and columns (324 after normalising the first entry) and is solved by a small
depth-first search (`search/tbsz.py`).

**Relation to earlier negative results.** [BORDERED.md](BORDERED.md) and NOTES.md §9 record that a 3-point border
constant on the cosets of an index-4 subgroup is impossible for bordered Williamson-type matrices, and that for the
complex array a border made of single characters is contradictory. The border here is not a character ansatz: the
rows of R̂ and the columns of Ĉ are arbitrary quaternary vectors in V_T, and Xb is determined by them.

## Verification

`python3 search/tbsz.py test` lists every core on ℤ₈ and ℤ₂ × ℤ₄ (N = 11) and ℤ₁₂ (N = 15) with Q = ℤ₄ and
T = {1, χ, χ̄} (32768, 49152 and 1,419,264 cores), completes samples of them with the border solver, and checks
K K* = 2N I exactly: complex Hadamard matrices of orders 22 and 30. `verify/verify_complex.py` accepts certificates
with construction `szollosi-tbordered` (group, fibre labels of φ, sequences and the three border matrices).

## Order 214

N = 107 is a safe prime (107 = 2 · 53 + 1): ℤ₁₀₇ has only the Paley multiplier group, and ℤ₁₀₆ = ℤ₂ × ℤ₅₃ only
H of order 13, which the ℤ₂ compressions rule out. N − 3 = 104 = 8 · 13 does have room:

- G = ℤ₁₀₄ = ℤ₈ × ℤ₁₃, φ(x, y) = x mod 4 (|K| = 26), Q = ℤ₄, T = {1, χ, χ̄};
- H = ⟨(1, 3)⟩ of order 3 (trivial on ℤ₈, the fourth powers {1, 3, 9} on ℤ₁₃): 40 orbits for C and D, 21 for A
  and B, 20 classes;
- row sums r_A, r_B; r_C, r_D = 2, 0; 2, 2.

On ℤ₄ the compressions have 128 feasible tuples with a positive conditioned estimate (logE 5.7 in all). The first,
y_A = (0, 2, −2, 2), y_B = (2, −2, 2, −2), y_C = (−4, 4, −4, 6), y_D = (0, 0, 2, 0), has pools of 5040, 2520,
1.97 · 10⁷ and 1.92 · 10⁸ sequences and a conditioned estimate of 10^3.8 solutions. A, B and C were listed in full and
D cut to its 2.9 · 10⁷ members of lowest energy (retention 0.5); the 52,762 split values on five merge classes are
exact (sumset intersection), and kt4 found 5 solutions in the first 256 of them, 338 seconds after the start. Each
completes to a complex Hadamard matrix of order 214 with a border found by the depth-first search in under a second.
The certificate is [C214.json](C214.json) (`verify/verify_complex.py 214`).

    python3 search/tbsz.py rank 107 --family z4cyc                 # the structures of N = 107 and their estimates
    python3 search/tbsz.py search 107 --family z4cyc --gen 1,3 --rows 2,0,2,2 --tuple 0 --log hits.jsonl
    python3 search/tbsz.py certify --hits hits.jsonl               # writes complex-hadamard/C214.json
