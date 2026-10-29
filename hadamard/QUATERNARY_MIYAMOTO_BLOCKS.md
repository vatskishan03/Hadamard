# Miyamoto's quaternary array in blocks, times T-matrices

[COMPLEX_MIYAMOTO.md](COMPLEX_MIYAMOTO.md) turns a quaternary Hadamard matrix K of order q − 1 (q ≡ 3 mod 4) into a
quaternary Hadamard matrix of order 4q made of sixteen bordered blocks. [MIYAMOTO_BLOCKS.md](MIYAMOTO_BLOCKS.md) shows
that for q ≡ 1 (mod 4) and a real K made of circulant or negacyclic blocks, Miyamoto's sixteen blocks commute and are
transposed by one permutation, so that T-matrices multiply them. Here the quaternary case is brought to the same point.
The sixteen quaternary blocks can be made to commute, but T-matrices cannot multiply them linearly. Their Turyn doubles
can be multiplied, and the product is a real Hadamard matrix of order 8qt.

**Theorem.** Let q ≡ 3 (mod 4) be a prime power, m = (q − 1)/2, and let K be a quaternary Hadamard matrix of order
q − 1 whose four blocks of order m are circulant (or all negacyclic). If T-matrices of order t exist, a Hadamard matrix
of order 8qt exists.

(`verify/verify.py` and `search/cmblocks` implement prime q; an independent check covered q = 27 = 3³ as well.)

The ingredients for K:

- *Paley's quaternary matrix of order (p + 1)/2*, p = 2q − 3 a prime
  ([complex-hadamard/PALEY.md](../complex-hadamard/PALEY.md)). It has exactly this form: the non-split torus τ of SL(2,
  p) has order p + 1; g = τ^m acts on Paley's skew matrix as the complex structure i, and ρ = τ⁴ of order m acts
  semiregularly on the 2m pairs {f, gf} with two orbits. In the basis P_ρ^j e_f of the two orbits, K is circulant in
  each block (`search/cmblocks/paley_circ.py`).
- *K = C + iI*, C Paley's symmetric conference matrix of order p + 1, p = q − 2 a prime ≡ 1 (mod 4), in the same way
  with ρ = τ² (`search/cmblocks/conf_circ.py`).

The quaternary matrices of order q − 1 found by search in this repository (Szöllősi's arrays) are not of this form.

## Construction

1. *Points and core.* List GF(q)* along a primitive root g: the squares g^(2j), then the non-squares g^(2j+1),
   j ∈ ℤ_m. The blocks of W = i [χ(y − x)] are circulant in j (χ(−1) = −1 makes the diagonal blocks skew circulant),
   U = diag(W, W), V = [[I, K], [−K*, I]], both cut into 4 × 4 blocks U_rc, V_rc of order m.
2. *Bordered blocks.* N_rc = [[A_rc, B_rc e], [B_rc eᵀ, [[U_rc + V_rc, U_rc − V_rc], [U_rc − V_rc, U_rc + V_rc]]]]
   (order q), with the same matrix B on the border row and the border column, and A, B 4 × 4 matrices of units with

       A A* = B B* = 4I,   A B* + 2 B R* = 0,   A − 2R = B,
       R = [[0, i, 0, 0], [−i, 0, 0, 0], [0, 0, 0, i], [0, 0, −i, 0]],

   R being the matrix of block row sums of U. There are 128 such pairs, and 128 more with A − 2R = −B, which work as
   well (`search/cmblocks/border_search.py`); the certificates use A = [[1, i, 1, 1], [−i, −1, i, −i], [−i, 1, i, i],
   [−i, −1, −i, i]], B = [[1, −i, 1, 1], [i, −1, i, −i], [−i, 1, i, −i], [−i, −1, i, i]]. The borders of
   [COMPLEX_MIYAMOTO.md](COMPLEX_MIYAMOTO.md) (B′ = R P F ≠ B) do not work here: 102 of the 120 pairs of blocks fail to
   commute.
3. *Turyn doubles.* X_rc = Re N_rc + Im N_rc, Y_rc = Re N_rc − Im N_rc, B_rc = [[X_rc, Y_rc], [−Y_rc, X_rc]]
   (order 2q, entries ±1).
4. *T-matrices.* With circulant T-matrices T_1, …, T_4 of order t, A_j = Σ_k T_k ⊗ B_kj (order 2qt), R″ = R_t ⊗ R′
   with R_t: s ↦ −s on ℤ_t, R′ = diag(1, −1) ⊗ R_B and R_B: j ↦ −j on each half of the core (border fixed), the
   Goethals–Seidel array

       [[A1, A2R″, A3R″, A4R″], [−A2R″, A1, A4ᵀR″, −A3ᵀR″], [−A3R″, −A4ᵀR″, A1, A2ᵀR″], [−A4R″, A3ᵀR″, −A2ᵀR″, A1]]

   is a Hadamard matrix of order 8qt.

## Proof

- *(N_rc) is a quaternary Hadamard matrix of order 4q.* This is the proof of [COMPLEX_MIYAMOTO.md](COMPLEX_MIYAMOTO.md)
  with B′ = B: the border identities there are A A* + 2m B B* = 4qI and A B′* + 2 B R* = 0, and the core rows close with
  any border B′ of unit entries with B′B′* = 4I.
- *The sixteen blocks commute.* Write N = [[a, b eᵀ], [b e, C]] with C = J₂ ⊗ U_rc + D₂ ⊗ V_rc. All the circulants U_rc,
  V_rc commute, so C C′ = C′ C, and C has equal row and column sums 2ρ (ρ = R_rc). Then N N′ − N′ N is zero except on
  the border, where it is (a b′ + 2 b ρ′ − a′ b − 2 b′ ρ) eᵀ: the blocks commute exactly when (A_rc − 2R_rc)/B_rc is one
  constant λ (necessarily ±1). The certificates take λ = 1.
- *One permutation transposes them.* R_B reverses every circulant and fixes the border, and the border row and
  column carry the same B_rc, so R_B N_rc R_B = N_rcᵀ (a linear transposition, not a conjugate one).
- *The linear plug-in fails.* With quaternary blocks the Goethals–Seidel array of A_j = Σ_k T_k ⊗ N_kj leaves cross
  terms (−A_1 conj(A_2) + A_2 conj(A_1))R + … in G G*; at q = 7, t = 3 the largest entry of G G* − 84I is 42
  (`search/cmblocks/cx_counter.py`), and an exhaustive search over Goethals–Seidel-type arrays with entries
  ±X, ±Xᵀ, ±conj X, ±X* times R finds none for generic quaternary blocks (`search/cmblocks/gs_search2.py`).
- *The Turyn doubles work.* Re N_rc and Im N_rc are each transposed by R_B, and diag(1, −1) turns the transpose of
  [[X, Y], [−Y, X]] into the double of N_rc*, so the sixteen real blocks B_rc commute and R′ B_rc R′ = B_rcᵀ.
  Since the double of a quaternary Hadamard matrix is Hadamard, Σ_c B_rc B_r′cᵀ = 8q δ_rr′. These are exactly the
  hypotheses of the Goethals–Seidel proof: commuting A_j, R″ A_j R″ = A_jᵀ, R″ a symmetric signed involution, and
  Σ_j A_j A_jᵀ = Σ_k T_k T_kᵀ ⊗ 8qI = 8qtI.

## Checks

`verify/verify.py` builds the certificates (construction `quaternary-miyamoto-blocks`): it checks K K* = (q − 1)I, the
border identities and the T-matrices, assembles the matrix and checks H Hᵀ = 8qt I. `search/cmblocks` checks the
intermediate claims (unit entries and H H* = 4qI of the quaternary array, commutation and R′-transposition of the
real blocks) exactly for q = 7, 11, 19, 23, 31 with both kinds of K and t = 3, 5, and for negacyclic K at q = 7, 11.
An independent check re-derived the proof and tested about 1000 further cases (q = 7 to 43 and 27, six sources of K,
every primitive root, all 256 borders at q = 7 and 11, t = 1 to 13). Every deliberate violation of a hypothesis (points
out of order, non-circulant K, borders with non-constant λ, R_B without the sign) breaks the construction.

## Consequences

Agaian's multiplication gives 8qt from Hadamard matrices of orders 4q and 4t; the theorem needs only a quaternary
one of order q − 1, so it matters when 4q is open. The q ≡ 3 (mod 4) with 4q open ([FRONTIER.md](../FRONTIER.md)) and
K of the right form include 751, 863, 907, 1103, 1187, 1291, 1367, 1423, 1427, 1471, 1487, …
(`search/cmblocks/scan_t4.py`). The matrix built here:

| Order | n | q | t | K | Table 4 |
|------:|---|---:|---:|---|---|
| [20712](20712) | 2589 = 3 · 863 | 863 | 3 | Paley's quaternary matrix of order 862 (p = 1723) | 2589(4) |

Within Table 4 the theorem reaches twelve entries n(m ≥ 4) = q t. For eleven of them (933, 1169, 1793, 1981, 2171, 2215,
2315, 2429, 2855, 2899, 2973) H(4q) is known, so Agaian's multiplication gives them already
([FRONTIER.md](../FRONTIER.md)); 2589 is the only one with 4q open. Beyond Table 4 it gives 75 orders 8n with
3001 ≤ n ≤ 9999 (for example 8 · 3 · 1103, 8 · 5 · 751, 8 · 5 · 863), and infinitely many above. Without Handbook
Table 1.53's values of t we cannot tell which of them were open, so they are not built.
