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
  the projection onto V_T, and R′ R′\*[(i, g), (j, h)] = (Ĉ Ĉ*)[(i, φg), (j, φh)] = 2 δ_ij Σ_{χ∈T} χ(φ(g) − φ(h)).
  Hence R′ R′* = I₂ ⊗ 2m P and Kc Kc* + R′ R′* = 2N I.
- **Border rows.** R R* = |K| R̂ R̂* = (m/q) · 2q I = 2m I. Mc commutes with P_T (P_T is I₂ ⊗ a symmetric
  Q-developed matrix) and Mc* Mc = I₂ ⊗ (Σ_k Y_k Y_kᵀ)/2, which is 2t on V_T. So Xb Xb* = R̂ Mc* Ĉ Ĉ* Mc R̂*/(4q²) =
  R̂ Mc* Mc P_T R̂*/(2q) = (t/q) R̂ R̂* = 2t I (the rows of R̂ lie in V_T), and Xb Xb* + R R* = 2N I.
- **Cross terms.** Summing a group-developed matrix over a fibre of φ gives the Q-developed matrix of the
  compression (also for X_D R and the conjugate transposes), so R Kc* = (R̂ Mc*) ∘ φ and Xb R′* = (Xb Ĉ*) ∘ φ. With
  Xb = −R̂ Mc* Ĉ/(2q): Xb Ĉ* = −R̂ Mc* P_T = −R̂ P_T Mc* = −R̂ Mc*. So Xb R′* + R Kc* = 0.

So K K* = 2N I; K is square, hence also K* K = 2N I. ∎

## Equivalent character sets, and t is odd

**T and λT.** Let λ be a real character of Q (λ² = 1). Multiplying the four sequences by λ ∘ φ multiplies
PAF_k(g) by λ(φ(g)) and keeps A and B symmetric, so it turns a core for T into a core for λT = {λχ : χ ∈ T}. The
compressions become λ y_k, Mc becomes (I₂ ⊗ Λ) Mc (I₂ ⊗ Λ) with Λ = diag(λ), and R̂ (I₂ ⊗ Λ), (I₂ ⊗ Λ) Ĉ and the same
Xb are a border for it: the new K is D K D with D = diag(I, I₂ ⊗ diag(λ ∘ φ)), an equivalent matrix. So {1, χ, χ̄} and
{χ², χ, χ̄} on ℤ₄ are one family, as are the two choices of T on ℤ₂² below and the two on ℤ₆.

**t is odd.** If 1 ∈ T and t is even, then m = N − t is odd, every row sum r_k is odd and Σ r_k² ≡ 4 (mod 8), while
Σ r_k² = 4t ≡ 0 (mod 8). If 1 ∉ T, each half of a row of R̂ is a quaternary vector with zero sum, so q is even, hence
m and t = N − m are odd. Since T is closed under conjugation and t is odd, T contains a real character λ, and λT
contains 1. So we may always assume 1 ∈ T; then Σ_k r_k² = 4t.

## Instances

- t = 1, Q trivial, T = {1}: the bordered array, ΣPAF = −4 off 0 ([BORDERED.md](BORDERED.md)).
- Q = ℤ_q, T = all characters (t = q): ΣPAF = −4q on K \ 0 and 0 off K. Needs q | m, i.e. q | N.
- **Q = ℤ₄, T = {1, χ, χ̄} (t = 3)**: Γ = 4 · (3, 1, −1, 1), so ΣPAF = −12 on K \ 0, −4 on the two cosets
  φ = ±1 and +4 on the coset φ = 2. Needs 4 | N − 3: every N ≡ 3 (mod 4). In particular it applies to the primes N for
  which neither ℤ_N nor ℤ_(N−1) has a usable multiplier group (107, 167, 179: orders 214, 334, 358).
- Q = ℤ₂², T = {1, χ₁, χ₂} (equivalently the three nontrivial characters): also t = 3 and N ≡ 3 (mod 4), on groups
  with a ℤ₂² quotient. Up to equivalence the t = 3 families are T = {1, χ₁, χ₂} on ℤ₂² and T = {1, χ, χ̄} with χ of
  order 3, 4 or 6.
- **Q = ℤ₆, T = all characters but the sign (t = 5)**: Γ = 4 · (5, 1, −1, 1, −1, 1), so ΣPAF = −20 on K \ 0, −4 on
  the cosets with φ odd and +4 on those with φ even. Needs 6 | N − 5: every N ≡ 5 (mod 6), among them 161, 167 and
  179 (orders 322, 334, 358). Here Σ_k r_k² = 20: row sums 4, 2, 0, 0 in some order.

For the ℤ₄ family the compressions satisfy Σ_k r_k² = 12 (r_k the row sums), Σ_k |ŷ_k(χ)|² = 12, and
Σ_k ŷ_k(χ²)² = 4N: three ±2 row sums and one 0 when |K| is even. V_T has 36 quaternary vectors per half, so the border
problem has 1296 candidate rows and columns (324 after normalising the first entry) and is solved by a small
depth-first search (`search/tbsz.py`).

**The border search.** The border depends only on the compressions. `search/border_exh.c` decides it exactly: the
candidates are the normalised quaternary vectors of V_T ⊕ V_T, a bit matrix records which row candidate u and column
candidate c have u Mc* c / 2q quaternary and another which candidates are orthogonal, and a depth-first search over
rows keeps the set of columns compatible with all chosen rows and finally looks for 2t orthogonal ones among them.
"No solution" proves that a compression tuple has no border. For Q = ℤ₆ there are 400 quaternary vectors per half
of V_T and 40,000 candidates; building the bit matrices takes about a minute and the search a fraction of a second.
A randomised search over bases of the form (v_i, ±c w_i) finds no ℤ₆ border at all; the exhaustive search finds one
for every ℤ₆ compression tuple we gave it.

**Relation to the negative results for the bordered array.** [BORDERED.md](BORDERED.md) and NOTES.md §9 record that a
3-point border constant on the cosets of an index-4 subgroup is impossible for bordered Williamson-type matrices, and
that for the complex array a border made of single characters is contradictory. The border here is not a character
ansatz: the rows of R̂ and the columns of Ĉ are arbitrary quaternary vectors in V_T, and Xb is determined by them.

## Zero-diagonal form: quaternary conference matrices

Keep the statement, but let A(0) = B(0) = 0 (A, B symmetric, ±1 elsewhere) and require Xb to be a quaternary
conference matrix of order 2t (one zero in each row and column). Condition 2 is unchanged; at g = 0 it now reads
Σ_k PAF_k(0) = 4m − 2. Then K = [[Xb, R], [R′, Kc]] is a quaternary conference matrix of order 2N: zero diagonal
(after a permutation), units elsewhere, K K* = (2N − 1) I. The proof is the one above with 4N and 4t replaced by
4N − 2 and 4t − 2: M = (2N − 1) I − 2m P, Mc* Mc = 2t − 1 on V_T, so Xb Xb* = (2t − 1) I, and the cross terms cancel as
before. For t = 1 this is the bordered zero-diagonal array of [hadamard/4204](../hadamard/4204). With Q = ℤ₂³, T its
seven nontrivial characters and G = ℤ₂³ × ℤ₂₉ it gives the conference matrix of order 478 behind
[hadamard/3820](../hadamard/3820) and [C1910.json](C1910.json); `verify/verify.py` checks all the conditions
(`_szollosi_zero_diagonal_core`, field `character_border`).

## Verification

`python3 search/tbsz.py test` lists every core on ℤ₈ and ℤ₂ × ℤ₄ (N = 11) and ℤ₁₂ (N = 15) with Q = ℤ₄ and
T = {1, χ, χ̄} (32768, 49152 and 1,419,264 cores), and on ℤ₆ (N = 11) and ℤ₁₂ (N = 17) with Q = ℤ₆ and T all characters
but the sign (9216 and 474,624 cores). It completes samples of them with the border search and checks K K* = 2N I
exactly: complex Hadamard matrices of orders 22, 30 and 34. `verify/verify_complex.py` accepts certificates
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

## Order 322

N = 161 = 7 · 23. Over ℤ₁₆₁ the plain array's multiplier groups with at most 64 orbits are empty or rate below
10^−17, and the bordered array over ℤ₁₆₀ is empty on ℤ₂ for every multiplier group with 33 to 50 orbits.
N − 5 = 156 = 12 · 13 has a ℤ₆ quotient:

- G = ℤ₂ × ℤ₂ × ℤ₃ × ℤ₁₃ (≅ ℤ₂ × ℤ₇₈), φ(x) = (x₀, x₂ mod 3) ∈ ℤ₂ × ℤ₃ = ℤ₆ (|K| = 26), Q = ℤ₆, T = all characters but
  the sign;
- H = ⟨M⟩ of order 6, M = (x₀, x₁, x₂, x₃) ↦ (x₀, x₀ + x₁, 2x₂, 3x₃): a shear on ℤ₂² that fixes K, −1 on ℤ₃ and an
  element of order 3 on ℤ₁₃; −1 ∉ H. 35 orbits for C and D, 20 for A and B, 19 classes;
- row sums r_A, r_B; r_C, r_D = 2, 0; 4, 0.

On ℤ₁₅₆ the best multiplier group, ⟨29⟩ of order 6, rates 10^0.3 conditioned on ℤ₆. On this group the pattern has 8
feasible ℤ₆ compression tuples and a conditioned estimate of 10^2.5 (plain Gaussian 10^−3.7). The first tuple,
y_A = (0, 2, −2, 2, −2, 2), y_B = (−2, 2, −2, 2, −2, 2), y_C = (2, 0, 0, 2, 0, 0), y_D = (2, −4, 4, −2, 4, −4), has
complete lists of 864, 432, 8.2 · 10⁶ and 1.3 · 10⁶ sequences and a conditioned estimate of 10^1.8. The 1367 split
values on five merge classes are exact. kt4 found 1024 solutions at the 246th of them, in order of the solution law,
80 seconds into the join. Every solution checked completes, with a border found by the exhaustive search, to a
complex Hadamard matrix of order 322. The certificate is [C322.json](C322.json) (`verify/verify_complex.py 322`).
In code, with `tbsz.py`, the structure is fixed as

    G = [2, 2, 3, 13]; phi(x) = (3 (x0 mod 2) + 4 (x2 mod 3)) mod 6; T = [(0,), (1,), (2,), (4,), (5,)]
    M = [[1, 0, 0, 0], [1, 1, 0, 0], [0, 0, 2, 0], [0, 0, 0, 3]]
    tbsz.search(tbsz.TStruct(G, [M], [6], phi, T), (2, 0, 4, 0), [0], log="hits.jsonl")
    python3 search/tbsz.py certify --hits hits.jsonl               # writes complex-hadamard/C322.json
