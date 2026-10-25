# Search notes: structures, obstructions and searches

Notes on the searches for the unresolved orders of [FRONTIER.md](../FRONTIER.md) and for the skew, symmetric and
complex orders: structures, obstructions, estimates and negative results. §1–§36, §38–§40, §48–§49 and §53:
Goethals–Seidel, Spence, skew and symmetric structures; §37 (in part), §41–§47, §50–§52 and §54–§59: complex orders
(Szöllősi arrays); §60–§83: decoupling constructions, conference-matrix cores and their obstructions (§74–§83 for the
complex Yamada–Kiyasu targets). Statements are proved here or measured with a program named in the text or listed in
[README.md](README.md). Results marked (program not included), and computations with no program named, which include
most of §60–§83, rest on code that is not in this repository. §k alone refers to this file; [NOTES.md](../NOTES.md) §k
to the structural notebook.

Theorems named by a letter: A is [symmetric/PROPUS_WEIGHING.md](../symmetric/PROPUS_WEIGHING.md), B is
[symmetric/YAMADA_KIYASU.md](../symmetric/YAMADA_KIYASU.md), F is §65, G is §67 and K is §68 (Spence's Corollary 2 at
r = 1). E and W are stated in §70.

## 1. The bottleneck is the coarse characters

For a group-developed construction over an abelian group G, write X_i(χ) for the
character sums of the four sequences. The Goethals–Seidel condition is

    Σ_i |X_i(χ)|² = 4|G|   for every nontrivial character χ,

and Σ_i X_i(1)² = 4|G| (the row sums). If the sequences are constant on the orbits
of H ≤ Aut(G), then X_i(χ) = Σ_O y_O η_O(χ) with periods η_O(χ) = Σ_{g∈O} χ(g).

Call a character *coarse* when its periods are rational integers with large
common factors, as with the trivial character (η_O = |O|) and characters of order
2. For a coarse character the condition is a four-squares equation in integers
that are fixed modulo the orbit sizes. That forces orbit sizes ≲ √|G|, and so at
least about √|G| orbits.

Over ℤ_p this is the √p bound of METHOD.md. Section 3 shows that the same
mechanism reappears in every variant tried here.

## 2. Cost scaling

With k ≈ 2C orbits, C classes, per-class standard deviation σ_c of the class
values, and pools N ≈ 2^k/σ_r, the model in search/triage.py gives
E ≈ N⁴ Π p_c with p_c ≈ 0.8/σ_c. The k-tree work per solution with lists of
size L is

    W ≈ N³/(E L) = 1/(N Π_c p_c L) ≈ (σ_r / L) · Π_c (σ_c / 3.2).

For σ_c ≈ √|G| this is (√|G|/3.2)^C · σ_r/L, exponential in C. The minimum
possible C is about √|G|/2. An order is cheap exactly when its group admits
symmetry with orbits of size about √|G| and an admissible row-sum pattern, as in
the m² + m + 1 family. For prime p this is a property of the divisors of p − 1.

## 3. Variants that do not escape the bound

### 3.1 Spence's construction and the slice character

Take G = ℤ₂ × K with |K| odd, and H acting on K. The character χ(z, k) = (−1)^z
has X_i(χ) = a_i − b_i, the difference of the two slice sums. The Spence
condition Σ_i |X_i(χ)|² = 4n with a_i + b_i = r_i = (2, 0, 0, 0) becomes

    (a₁ − 1)² + a₂² + a₃² + a₄² = n.

Each a_i is a signed sum of K-orbit sizes, so a_i ≡ ±1 (mod 2h) when all
nontrivial orbits have size h. This is again a √n bound, on K, which has half the
order. Spence structures therefore need about 2√n orbits, twice as many as
Goethals–Seidel over ℤ_n. The triage behind FRONTIER.md omits this condition.

Measured (`calibrate.py`, `count4.c`):

| structure | Gaussian model | exact count |
|---|---:|---:|
| ℤ₂ × ℤ₅², \|H\| = 3 | 10^9.5 | 10^10.5 |
| ℤ₂ × ℤ₁₁², \|H\| = 15 | 10^4.1 | 0 |
| ℤ₂ × ℤ₁₃², \|H\| = 21 | 10^2.9 | 0 |

In both zero cases the slice equation has no solution. `spence_slice_check.py`
checks every multiplier subgroup for each FRONTIER Spence entry:

- 3356 (ℤ₂ × ℤ₄₁₉, listed at log W 15.0) has no admissible structure with ≤ 400 orbits.
- 8828 and 11420 have no admissible structure with ≤ 400 orbits.
- 4748 (listed at 15.2) is admissible only with |H| = 8 and 150 orbits.

The cheapest feasible entry of the table is therefore 4412 (Goethals–Seidel, 10^15.5).

### 3.2 Full automorphism groups (Singer cycles) on elementary-abelian parts

In characteristic 2, orbits of Singer subgroups that lie in a subfield are
subspaces. For example, with H of order 3 on ℤ₂⁴ the orbits are the lines of a
spread, and every Walsh coefficient collapses to Ŵ(a) = c₀ + 4u_{O(a)}. On
ℤ₂⁴ × ℤ₃ the model predicts 10^7.15 solutions. Exact enumeration (`count_exact.py`,
68 s) finds 0. In odd characteristic the slice constraint of 3.1 applies.

### 3.3 Non-cyclic groups for Goethals–Seidel

The only open order with a non-cyclic group of order n is n = 1719 = 3² · 191
(`gs_1719.py`). On ℤ₃² × ℤ₁₉₁ the shifts inside ℤ₃² see whole ℤ₁₉₁-orbits, so
their class values move in steps of the orbit size, and the row sums fix the
residues. Every structure fails this lattice test except two, and those have
log E ≈ −30.

### 3.4 An 8-block Goethals–Seidel array

A ratio-16 construction over an odd group with a border would remove every real
character. `gs_array_gf2.py` reduces the existence of a Goethals–Seidel-type
realization (entries ±X, ±Xᵀ, ±XR, ±XᵀR) of an orthogonal design to a linear
system over GF(2):

- OD(4; 1,1,1,1) is solvable. This is the Goethals–Seidel array.
- OD(8; 1⁸) is inconsistent: 336 equations, rank 105.

### 3.5 Central elements of non-abelian groups

Class functions commute, so the Goethals–Seidel array works for them. For
G = ℤ_q ⋊ ℤ_r the conjugacy classes outside ℤ_q have size q, and the linear
characters then force (a_i + 2q ε_i)² and (a_i − q ε_i)² to both sum to 4n. This
is impossible since |a_i| ≤ q.

## 4. What is left: small structures and Legendre pairs

- The FRONTIER table lists only structures with E ≥ 10, and the exhaustive pass
  in METHOD.md stopped at 18 cosets. Re-running triage.py over all open orders
  with 8–70 orbits turns up small-pool Goethals–Seidel structures with
  E ≈ 10^−1 to 10^−5. Each is unlikely to hold a solution, but each can be searched exhaustively
  (`find4.c`, `gs_exhaustive.py`, validated by constructing H(148) and H(172)).
  The best is n = 933 (order 3732), |H| = 62, 22 orbits, E ≈ 10^−0.9.
- Legendre pairs (ℓ = 2n − 1) have no coarse character except the trivial one,
  and the border absorbs it. They are 2-list problems, linear in the pool. The
  repository never searched them. `lp_search.py` recovers the classical pairs and
  many others for prime ℓ. For composite ℓ, the only case open orders produce,
  small-structure pairs are rare.

## 5. Exhaustive searches of small structures

All searches below are exhaustive (`gs_exhaustive.py`, `find4b.c`). Each verified
solution would have been rebuilt and checked as a Hadamard matrix.

| n | \|H\| | orbits | classes | pattern | pools | pairs | log₁₀ E (Gaussian) | solutions |
|---:|---:|---:|---:|---|---|---:|---:|---:|
| 933 | 62 | 22 | 11 | (1,1,3,61) | 1.16·10⁵ ×3, 1.12·10⁵ | 2.7·10¹⁰ | −0.9 | 0 |
| 1981 | 141 | 21 | 10 | (1,1,1,89) | 4.0·10⁴ ×3, 1.9·10⁴ | 2.4·10⁹ | −3.2 | 0 |

`lottery_rank.py 4 40` ranks every Goethals–Seidel structure of every open order
with 4 to 40 orbits (−1 ∉ H). Eleven structures are admissible. Besides the two
above, the ones with log₁₀ E > −6 are:

| n | \|H\| | orbits | log₁₀ E | exhaustive pairs |
|---:|---:|---:|---:|---:|
| 933 | 31 | 33 | 0.3 | 2.8·10¹⁷ |
| 1893 | 90 | 30 | −3.1 | 1.3·10¹⁵ |
| 955 | 38 | 33 | −3.1 | 3.1·10¹⁶ |
| 2973 | 110 | 38 | −4.6 | 6.5·10¹⁹ |
| 1255 | 50 | 33 | −5.0 | 6.3·10¹⁶ |

None of these is affordable here.

`spence_rank.py 44` does the same for Spence structures ℤ₂ × ℤ_m (n ≡ 3 mod 4)
that pass the slice condition of §3.1. There are 502 with at most 44 orbits.

| n | \|H\| | orbits | pairs | log₁₀ E |
|---:|---:|---:|---:|---:|
| 599 | 33 | 30 | 1.9·10¹⁴ | −6.8 |
| 947 | 35 | 42 | 1.7·10²¹ | −8.4 |
| 647 | 24 | 36 | 7.7·10¹⁷ | −10.5 |

Every other structure is below 10⁻¹¹. METHOD.md reports that the Gaussian model
for Spence can be off by a factor of 100, which does not change the conclusion. So
neither small Goethals–Seidel nor small Spence structures leave a structure worth
searching. The only exception is n = 933 with |H| = 31, which has about two expected
solutions but needs about 10¹⁷ pair operations.

## 6. Cooper–Wallis routes: T-matrices for the open orders

FRONTIER.md lists Goethals–Seidel and Spence structures only. Every composite
open n = t·w with w a Williamson order (w ≤ 63) has t ≥ 167. So these routes need
T-matrices of a prime order t ≥ 167:

- T(179) would give 2148, 7876 and 9308.
- T(191) would give 2292, 3820, 6876 and 11460. 9172 then follows by Miyamoto from 2292.
- T(223) would give 2676 (and through it 10708), 4460 and 11596.
- T(263) would give 3156, 5260, 7364 and 11572.

**Reformulation.** Put (A, B, C, D) = H₄ (T₁, T₂, T₃, T₄). T-matrices of order t
are exactly the Goethals–Seidel quadruples of length t with A·B·C·D ≡ 1
pointwise; the repository's T(97), T(103) and T(197) all satisfy this. With
E = AB = CD the condition becomes f_E(A) + f_E(C) = 0, where
f_E(X) = P_X + P_{XE}. For fixed E this is a 2-list self-join.

**Symmetric T-matrices.** With H ≤ ℤ_t* of order h, the block containing 0 has
row sum ±1 + hj and the other three have row sums hj, with Σ rᵢ² = t. For a d = 3
twist (T₁ invariant under K = ⟨H, μ⟩, and T₃, T₄ the images of T₂ under μ, μ²),
r₁² + 3r² = t. A K-orbit can either lie in T₁, or have T₂ on exactly one of its
three H-cosets, so each K-orbit has 8 states.

`tmat_scan.py` lists the admissible structures with at most 40 classes:

| route | n | t | structure | K-orbits | classes |
|---|---:|---:|---|---:|---:|
| 7572 | 1893 = 3·631 | 631 | d = 3, \|H\| = 7, \|K\| = 21, rows (22; 7, 7, 7) | 30 | 15 |
| 9260 | 2315 = 5·463 | 463 | plain, h = 11, rows (10, 11, 11, 11) | 42 | 21 |
| 5020 | 1255 = 5·251 | 251 | plain, h = 5 | 50 | 25 |
| 3732 | 933 = 3·311 | 311 | plain, h = 5 | 62 | 31 |
| 2676, 4460, 11596 | 669, 1115, 2899 | 223 | plain, h = 3 | 74 | 37 |

Structures that are dead:

- T(179) and T(263): (t − 1)/2 is prime, so ℤ_t* has no useful subgroup.
- T(191): every subgroup, plain or twisted, fails the row-sum equation.

For T(631) with the d = 3 structure, `tmat_t3_model.py` predicts 10^4.7
solutions, most with 5 to 9 K-orbits in T₁. This is the only route among the open
orders with 15 classes and a large E. But it couples the blocks:

- Once the 30 K-orbit types are fixed, only two independent sign lists remain, T₁
  and the μ-orbit of T₂.
- In Goethals–Seidel form (A, B, μB, μ²B), the condition A = B·B^μ·B^{μ²} makes it
  a 1-list problem.

The best split (k₁ = 15, fixed types with parity, a 2-list join over the signs)
costs about 10^13.6 operations per solution. This is no better than the d = 3
structure over ℤ₁₈₉₃ in METHOD.md (10^14.7).

**A four-list split** (`tmatrix.py`). Fix E = supp T₁ ∪ supp T₂ and put X₁ = T₁ + T₂, X₃ = T₁ − T₂
(±1 on E) and X₂ = T₃ + T₄, X₄ = T₃ − T₄ (±1 off E). Then ΣPAF(Xᵢ) = 2·ΣPAF(Tᵢ), and conversely any four
such ±1 sequences give T-matrices. So for fixed E the four lists are independent, and the
Goethals–Seidel matchers apply (find4b exactly, match4 by k-tree). Base sequences are the case where E is
an interval. Validation: for each of the three ways to pair up the blocks of T(113), the exact join on
that E finds the known solution.

The split does not change the count. A plain structure has 3 bits per orbit and C ≈ O/2 classes. With P
the chance that one class sums to zero, it holds (64P)^C T-matrices, but only (16P)^C of them on a fixed
E. Here P ≈ 1/(2.5·√(t/4)), which falls below 1/16 near t ≈ 160. Per E, the Gaussian model predicts:

| t | h | classes | solutions per E |
|---:|---:|---:|---:|
| 113 | 7 | 8 | 10^0.6 |
| 251 | 5 | 25 | 10^−2.4 |
| 463 | 11 | 21 | 10^−4.8 |

With full lists of 2^{O/2}, the cost is 10^17 or more per solution for every open-order route in the table
above. This matches where the T-matrices of [t-matrices](../t-matrices) stop (t ≤ 197, at most 36 orbits).

## 7. Other routes checked

- **Recursive constructions.** Yamada (J. Austral. Math. Soc. 1989) gives
  H(4(q + 2)) for a prime power q = n − 2:
  - from a skew-Hadamard matrix of order (q + 3)/2 when q ≡ 5 (mod 8), a result
    credited to Kiyasu;
  - from a symmetric conference matrix of order (q + 3)/2 when q ≡ 1 (mod 8).
  
  For the open orders:
  - The skew orders needed are 940, 1108 and 1192 (through 596). All three are
    unknown in Cati–Pasechnik (n = 235, 277, 149). ℤ₁₄₉ and ℤ₂₇₇ have no
    multiplier structure for a skew first block.
  - The conference route needs (n − 1)/2 to be a sum of two squares. This holds
    only for n = 955, which needs a symmetric conference matrix of order 478, that
    is a conference graph on 477 = 3²·53 vertices. Abelian Paley partial difference
    sets need a prime power order, or 9m⁴ or m⁴, so no Cayley construction exists.
- **Legendre pairs** (ℓ = 2n − 1, H(2ℓ + 2) = H(4n)). Two ±1 sequences give 2ℓ
  signs against (ℓ − 1)/2 equations, so 4 signs per equation instead of the 8 of
  Goethals–Seidel. Each equation P_A(s) + P_B(s) = −2 has spread about √(2ℓ) on a
  lattice of spacing 4. The Gaussian count is then about
  [4·(1.13/√ℓ)^{1/2}]^{(ℓ−1)/h} for a structure with orbits of size h, and the
  bracket is below 1 for ℓ > 326. Legendre pairs of length ℓ ≥ 1073, the smallest
  length an open order needs, can only come from algebraic structures (Paley, twin
  primes) or from a few orbits by luck. `lp_sweep.py 3 24` searches every
  multiplier structure with at most 24 orbits for every open n, exhaustively:
  3,703 structures across all 151 orders, with no Legendre pair.
- **Projection methods.** The Elser RRR / difference-map iteration on
  A = {±1}^{4v}, B = {Σᵢ|X̂ᵢ(k)|² = 4v ∀k} (`dm/rrr_gs.py`) solves unstructured
  quadruples for v ≤ 21 in about 10⁴ iterations, but none for v = 23 to 31 in
  5·10⁵. The symmetric version stays in the H-invariant subspace exactly, since
  X(hk) = X(k) (`dm/rrr_orb.py`). It fails in 2·10⁵ iterations on the order-4492
  structure (1123, |H| = 33), which the k-tree solves in seconds. Not competitive.
- **Limits on the number of lists.** Plug-in constructions (Cooper–Wallis,
  Baumert–Hall, orthogonal designs with circulant variables) give independent
  k-tree lists only through the variables of an orthogonal design of order 4m with
  m odd. By Radon–Hurwitz there are at most ρ(4m) = 4. Splitting a sequence into
  disjoint parts leaves bilinear cross terms, and T-matrix-type cancellation needs
  4 × m sign matrices M with MᵀM = 4I, so m ≤ 4. Four lists is the ceiling for
  orders 4·odd.
- **Schur rings.** Every Schur ring over ℤ_p is an orbit ring of a subgroup of ℤ_p*
  (Pöschel). So for prime n the only sequence families whose autocorrelations stay
  in a fixed proper subring are the multiplier structures, and the √p bound of
  METHOD.md covers all of them.

## 8. Twisted Legendre pairs, local search, base sequences

- **Twisted Legendre pairs** (`tlp_model.py`, `tlp_search.py`). This is the pair
  (A, A∘μ) with μ² ∈ ⟨H, −1⟩ and μ ∉ ⟨H, −1⟩. It is the form of the cyclotomic
  pairs for primes ℓ ≡ 1 (mod 4), and composite ℓ = 2n − 1 has many such μ. The
  condition P_A(s) + P_A(μs) = −2 holds on the classes of ⟨H, μ, −1⟩. The search
  is a 1-list problem, linear in the pool. It recovers verified pairs for
  ℓ = 37, 41, 61, 73 and 89.
  - For 2148 (ℓ = 1073) the best structures are |H| = 63, 25 orbits, 6 classes,
    pool 4.6·10⁵, with Gaussian E ≈ 10⁻¹·⁴ for each of two choices of μ.
  - All twisted structures with ≤ 28 orbits for n = 537 and 573 were searched
    exhaustively, with no pair.
  - Every other open order has E < 10⁻³·⁴ in its best twisted structure.
- **Open Legendre-pair lengths.** The lengths below 200 open as of the RICAM
  report on lengths 85 and 87 are 115, 145, 159, 161, 169, 175, 177, 185, 187 and
  195. Plain and twisted multiplier structures with ≤ 26 orbits were searched
  exhaustively for all of them (`lp_open.py`), with no pair.
- **Local search.** These are steepest-descent tabu searches with all single flips,
  in C (`ls/`).
  - Unstructured Goethals–Seidel: v = 41 is solved in seconds, v ≥ 51 is not solved
    in 10⁸ steps.
  - The steps needed grow like (solution density)^−0.42. That predicts about 10¹⁹
    steps for the 4412 structure and about 10⁴⁰ for an unstructured length of 167.
  - In orbit coordinates the search fails even on the 4492 structure (17 classes),
    which the k-tree solves in seconds. Large orbits make every move a large jump in
    the class values.
  - Base sequences BS(m + 1, m): m = 15 is solved in 4 s, and m ≥ 20 stalls at
    energy 4.
- **Base sequences as a universal route.**
  - Turyn: BS(m + 1, m) gives T-sequences of length 2m + 1, which give T-matrices.
    So BS(90, 89) would give T(179) and, with Williamson matrices of order 3, 11
    and 13, the orders 2148, 7876 and 9308.
  - Gaussian count for BS(90, 89): 358 signs against 89 aperiodic equations whose
    spread shrinks with the shift, about 2³⁸ solutions. They should exist.
  - Known methods (end-in backtracking, compression) stop near m = 40.
- **Composite-t Cooper–Wallis** (`tmat_composite.py`). Combining the heuristic
  E_T ≈ E_GS·2^{−(k−1)} for the condition A·B·C·D = 1 with the Goethals–Seidel
  estimates, every structure over a composite t = n/w (w ≤ 17) has E_T < 10⁻⁷.
- **Throughput here.** The repository's JAX generator runs at 2.3·10⁴ candidates/s
  on this 4-core CPU. Calibration on the 4492 structure gave 7.6·10⁵ distinct
  candidates per row sum and no match in 600 s; the original runs used 1.6·10⁸.
  The 4412 search (pilot: 10⁻⁹ expected solutions per τ) is about four orders of
  magnitude beyond the computing used here.
- **The d = 3 structure for 7572, run here** (`d3/prep1893.py`, `d3/walk.c`). This
  is the METHOD.md structure over ℤ₁₈₉₃: H = ⟨8⟩ of order 30, K = ⟨8, 2⟩ of
  order 90, 30 and 86 orbits, 15 classes, row sums (87; 1, 1, 1).
  - The complete first-block pool has 24,608,948 sequences, all stored in an exact
    hash table.
  - The model gives log₁₀ E = 2.15, and log₁₀(success per exactly centred lookup)
    = −14.7, matching METHOD.md.
  - The tabu walk over the second block checks every row-sum-preserving swap
    neighbour: 6.5·10⁶ lookups/s on 4 threads. That is at least 900 days per
    expected solution on 4 threads.

## 9. Structural ideas checked: the construction of the August 2026 matrices, and obstructions

- **How the August 2026 matrices are built.** E. Piza (Zenodo 10.5281/zenodo.23244254)
  verified all twelve and identifies their structure:
  - 668, 716, 1676 and 1772 are bordered (Wallis–Whiteman/Spence). The 668 seed is
    four cyclic sequences over ℤ₁₆₆ with sums (2, 0, 0, 0), a
    (166; 82, 83, 83, 83; 164) supplementary difference set.
  - 892, 1132, 1244, 1948 and 1964 are Goethals–Seidel.
  - The other three are not identified.
  
  So the 716 matrix carries no T-matrices of order 179, and the 2148 Cooper–Wallis
  route does not open.
  
  Their Goethals–Seidel orders have small odd multiplier subgroups. The Gaussian
  model gives E ≈ 10^39.3 and k-tree work ≈ 10^10.5 for p = 223 with h = 3. The
  triage in FRONTIER.md stops at 110 orbits, but for p ≳ 400 small h means
  hundreds of orbits and astronomically large W, so nothing comparable exists for
  the open orders.
- **Bordered T-matrices: a parity obstruction.** Suppose T₁..T₄ are 0/±1 sequences
  over a group of order t − 1 with partitioning supports, Σ P = −1 and row sums
  (±1, 0, 0, 0). The Cooper–Wallis array with border −Q, Q⊗1ᵀ, GS(c)⊗1 would then
  be an OD(4t; t, t, t, t), hence H(4tw) for every Williamson order w; the algebra
  checks.
  - They cannot exist. The trivial character forces Σ r_k² = 1, so A = ΣT_k is a
    ±1 sequence of even length with odd row sum.
  - Exhaustive checks over all groups of order 4 and 8 find none.
  - (The order-2 character would separately force t ≡ 1 mod 4.)
- **Recursive chains with the 2026 matrices.**
  - Miyamoto: H(q − 1) ⇒ H(4q) for a prime power q ≡ 1 (mod 4).
  - Yamada: H((q − 1)/2) ⇒ H(4q) for q ≡ 1 (mod 8).
  
  Applied to all newly known orders, these settle no FRONTIER entry:
  - 7796 (n = 1949) is already credited (Miyamoto from 1948).
  - 5732 and 11108 were never open.
  - 17972 and 15332 are Paley I orders.
  - 15716, 14132, 14708, 21748 and 31412 lie beyond n = 2999, where no authoritative
    open list exists.
- **Singer cycles on the 2-part of Spence groups.** Example: 2148, G = ℤ₂³ × ℤ₆₇, with
  τ ∈ GL(3,2) of order 7 combined with a multiplier.
  - The seven order-2 characters collapse into one condition,
    (1 − 4α₁)² + 16Σα² = 537, which is solvable.
  - But the Gaussian estimates for τ-twisted structures are below 10⁻¹⁷⁰, the
    characteristic-2 degeneracy of §3.2. Plain multipliers on ℤ₆₇ reach only
    E ≈ 10⁻⁰·⁵ with 31 classes.
- **Algebraic-curve sequences.** Take A(x) = χ(x² − a), the Legendre symbol of a
  quadratic.
  - P_A(s) = −1 − a_p(E_{a,s}) + O(1) boundary terms, where E_{a,s} is the elliptic
    curve y² = (x² − a)((x + s)² − a).
  - Quadratic twists (λ ↦ 1 − λ, with χ(−1) = −1) cancel the traces in pairs. Such
    quadruples land naturally on Σ P = −4, the Spence target for a group of order
    n − 1 = p. They never land on 0.
  - For the bordered Spence group ℤ₂ × ℤ_m with Legendre slices, the cross-slice
    correlations have spikes of size m that cannot be cancelled.

## 10. Skew-Hadamard orders

Skew-Hadamard matrices (H + Hᵀ = 2I) have their own open problem: Table 5 of Cati–Pasechnik
lists 234 odd v ≤ 999 for which no skew matrix of order 4v is known. The FRONTIER obstruction
(large coarse-character classes) is milder here, so small multiplier structures still work.

- **Structure.** Goethals–Seidel over ℤᵥ, all four blocks invariant under H ≤ ℤᵥ* with −1 ∉ H, and the first block of
  skew type (A(0) = 1, A(−x) = −A(x)). Negation pairs the H-orbits, so there are 2^{(O−1)/2} skew first blocks, and
  the row sums satisfy 1 + b² + c² + d² = 4v.
- **Triage** (`skew_triage.py`, `skew_rank.py`) uses the Gaussian model of §2 with the skew A pool. **For composite v
  it is not to be trusted at all**: it is off by 10²–10⁵ in both directions, which is the same lesson as order 3964.
  Composite v are searched exactly (`chain.py`, below); the model only orders work for prime v.
  - Under: 209 with ⟨26⟩ was predicted at 3 solutions and has 288.
  - Over: 201 with ⟨22⟩ (|H| = 11, rows 1, 1, 19, 21) was predicted at 10³ and has none (exhaustive `find4b`).
  - Under: 377 with ⟨16⟩ was predicted at 0.1 and gave 64 in the first τ; 469 with ⟨9⟩ was predicted at 0.2 and has
    192; 889 with ⟨2, 25⟩ came from the exact sweep, not from the model.
- **Matching.** `match4.c` is `search/match.c` with four distinct lists: left pairs (A, largest pool), right pairs
  from the other two. With lists of very different sizes, the yield per τ grows with the smaller pair list, so the
  clamp S is set as loose as memory allows (pair lists ≤ 2.5·10⁷). Pools beyond 4·10⁶ are replaced by uniform samples
  (`--maxpool`), at a cost of L/N per list.
- **New orders**: 836, 1036, 1100, 1148, 1236, 1364, 1412, 1508, 1612, 1852, 1876, 2404, 3556, 3676 and 1372, i.e. v =
  209, 259, 275, 287, 309, 341, 343, 353, 377, 403, 463, 469, 601, 889, 919. 889 came from the exact sweep, on H = ⟨2,
  25⟩ of order 63 with 21 orbits. Each took seconds to minutes on 4 cores.
  - 259 used sampled lists: 3.9·10⁶ of 3–4·10⁷ per pool.
  - 353 (⟨22⟩) also used sampled lists: 1.2·10⁷ of 5.7·10⁸ per pool, 2 solutions in 33 values of τ.
  - 403 (⟨100⟩) used sampled lists of 1.05·10⁷ of 4–6·10⁷: 4 solutions in 154 values of τ on rows (1, 5, 25, 31). A
    119 run with low-energy lists (8·10⁶ per role) found nothing.
  - 469 (⟨9⟩, predicted 0.2) gave 192 solutions in the exact matcher. The k-tree had missed them because 36 values of
    τ with the top-6 values per class do not cover the solutions. For pools of ≤ 10⁵, run the exact matcher.
- **Spin structure** (`spin_skew.py`, `spin_ls.c`), the type of Đoković's skew 4·397.
  - A is invariant under K = ⟨H, μ⟩ with μ of order 3, and C = B∘μ, D = B∘μ².
  - The condition P_A + P_B + P_B∘μ + P_B∘μ² = 0 lives on the ⟨K, −1⟩-classes, a third as many as the ⟨H, −1⟩-classes.
    It is a 2-list problem: B ranges over 10¹³–10¹⁵ choices, and a tabu search on (A, B) in C does about 5·10⁴
    moves/s.
  - The row sums are forced: 1 + 3b² = 4v. For open v ≤ 999 that leaves v = 469, 721, 817 and 919.
  - 919 (|H| = 17, μ = 52, 9 classes, predicted 2·10⁵ solutions) was solved in 37 s.
  - 721 (|H| = 17, μ = 46, predicted 200): no solution after 6·10⁶ moves (best residual: two classes off by 4); the
    structure is infeasible by the fibre test of §36 (μ = 46 is equivalent to μ = 18).
  - 817: μ³ = 1 gives no structure with E > 10⁻³, but μ³ ∈ H suffices: the blocks B, B∘μ, B∘μ² close up since B is
    H-invariant. With μ³ ∈ H, H = ⟨11⟩ (order 21) and μ ∈ {6, 9, 17, 24} give spin structures with 8 classes and model
    E ≈ 10³·². Tabu searches on them found no solution; for μ = 9 and 6 the fibre test of §36 shows the structure
    infeasible.
  - Since a skew block has row sum 1, spin needs 4v = 1 + 3b², i.e. v = 3m² + 3m + 1, a centred hexagonal number.
    Đoković's 397 and the DGK order 631 are of this form, as are the four values above.
  - `spin_skew.py exhaust` joins all A against all B exactly when the pools are small. For 817 (|H| = 63, μ = 49 or
    552) and 721 (|H| = 51) it finds no solutions; most of these structures have no B with row sum b at all.
- **General local search.** The same C engine handles plain skew GS: B, C, D are stacked into one vector with a
  block-diagonal form, and swaps stay inside a block (`spin_skew.py gsearch`). It solves 209/⟨26⟩ in under a second.
- **Compression obstruction** (`skew_compress.py`). This explains the Gaussian model's misses.
  - For d | v, the compressions to ℤ_d are complementary: Σ_X PAF_x(k) = 4v·[k = 0] for k ∈ ℤ_d. This is an identity
    of integers, equivalent to Σ_X |X̂(χ)|² = 4v at every character. The compressions are also constant on the
    H-orbits of ℤ_d, and the skew block has x₀ = 1 and x₋ⱼ = −xⱼ.
  - Enumerating each role's achievable compressed vectors and meeting the integer equations is exact and cheap for
    small d. These are the coarse characters of §1 again.
  - **The test must be integral.** For d ≥ 5 the norms |X̂(χ)|² are not rational integers (they lie in ℚ(ζ_d)⁺), so
    comparing them after rounding to integers is not exact; only its d = 3 verdicts are sound. The verdicts below use
    the integer test.
  - It proves there are no solutions for:
    - 201/⟨22⟩, rows (1, 1, 19, 21): d = 3, predicted 10³;
    - 301/⟨4⟩, rows (1, 1, 19, 29): d = 7, predicted 4;
    - 265/⟨16⟩, rows (1, 1, 23, 23): d = 5.
  - The solved structures pass. 209/⟨26⟩ has 11 and 21 feasible combinations at d = 11 and 19; 469/⟨9⟩ has 1 at d = 7.
  - Any triage of composite v should apply this test before the Gaussian estimate.
- **Which tool when.** The k-tree enumerates every left pair (A, D) and every right pair (B, C) once over all values
  of τ. Its work per solution is therefore (N_A·N_D + N_B·N_C)/E′, with E′ the yield of the (sampled) lists. At about
  10⁷ pairs/s:
  - 353 (E′ = 10⁴, lists 1.2·10⁷): about 10 minutes;
  - 763 (E′ = 10²·⁶): about 10 hours, and it gave no hit in 400 τ;
  - 703 (E′ = 10¹·⁸): days.
  
  Low-energy lists (`skew_gen.c`) narrow each list's class spread by 2–2.4×. With A uniform
  this adds roughly 0.21 per class to log₁₀ E′, which matters most at 15–25 classes.
  
  Joint tabu search on (A, B, C, D) works up to about 10 classes:
  - 209 (10 classes) and the spin structure of 919 (9 classes): seconds;
  - 353 (16 classes, E = 10⁹, solved by k-tree): no solution in 20 minutes;
  - 205 (22 classes, E = 10¹¹): best residual about 6 classes off by 4 after 8·10⁶ moves.
  
  For the spin structure, a walk that looks up all 512 skew A's at every swap neighbour (the
  METHOD.md technique for 2284, `spin_walk.c`) did 4·10⁹ lookups in 300 s without finding 919,
  which the joint tabu finds in under a minute.
- **Second-tier sweep** (pair work 1.5·10⁹–5·10¹⁰, ≤ 37 orbits). Its structures are covered exactly by the
  compression-chain sweep below; verdicts of the rounded-norm test are not used.
- **Compression-chain engine** (`chain.py`, `chain_sweep.py`), after the search for order 3964. For v = p·q
  (or a prime power):
  1. A DP over the H-orbits lists every block's pairs of compressions (ℤ_d1, ℤ_d2) with multiplicities.
  2. The ℤ_d1 compressions are joined exactly. Inside each branch, the ℤ_d2 compressions are joined.
  3. Each consistent combination is expanded by meet in the middle, and the four lists are joined on all
     classes with find4b.
  
  No step samples, so with no truncation a structure that yields nothing is ruled out completely.
  **Validation:** on 209/⟨26⟩, rows (9, 15, 23), the exhaustive chain counts 288 solutions, exactly the
  count of the direct four-list join. The sweep runs over the 140 composite unsolved v of Table 5,
  for every H with −1 ∉ H and up to 46 orbits, ignoring the Gaussian model.
- **No hits:**
  - 201 with ⟨22⟩ and 301 with ⟨4⟩ (both exhaustive; both obstructed, see above);
  - 265 with ⟨16⟩ (400 values of τ);
  - 623 with ⟨2⟩ (400 values of τ on both patterns).
  - `skew_sweep.py` runs the exact matcher on every structure with ≤ 33 orbits and pair work
    ≤ 1.5·10⁹.
- **Hadamard orders through Yamada's theorem.** If q ≡ 5 (mod 8) is a prime power and a skew matrix
  of order (q + 3)/2 exists, then H(4(q + 2)) exists. Among FRONTIER orders this needs exactly the
  skew orders 940, 1108 and 1192, for n = 1879, 2215 and 2383. All three are out of reach of skew
  multiplier structures:
  - **940 = 4 · 235.** −1 ∉ H forces the projections of H to ℤ₅* and ℤ₄₇* to avoid −1, so
    H ≤ {1} × QR₄₇. There are then 15 orbits (five fixed points and ten of size 23), and
    b² + c² + d² = 939 has no solution with b, c, d ∈ ±{1, 3, 5} + 46ℤ.
  - **1108 = 4 · 277.** The orders of H with −1 ∉ H are 1, 3, 23 and 69. Order 23 leaves row sums
    ±1 + 46ℤ, and 1107 is not a sum of three such squares. Order 3 has 46 classes and pools near
    10²⁷, beyond the k-tree.
  - **1192 = 8 · 149.** A skew Goethals–Seidel block needs a group of odd order, so this order must
    come from skew 596 by doubling. For v = 149, H has order 37 (row sums ±1 + 74ℤ: no pattern) or
    1 (unstructured).

## 11. Hadamard orders: the next tier of small structures

FRONTIER.md's exact sweep stopped at 22 orbits. `had_sweep.py` extends it.
- **Scope.** All 151 unresolved n, every subgroup H ≤ ℤ_n* with 19–40 orbits, every row-sum
  pattern a² + b² + c² + d² = 4n with exact pair work ≤ 2·10¹¹.
- **Filter.** Each pattern passes through the compression test of §10, now for plain (non-skew)
  blocks. On n = 63 it agrees with exhaustive search: patterns marked feasible have solutions
  (49,104 for H = ⟨2⟩, rows (1, 1, 5, 15)), and patterns marked dead have none.
- **Result.** Only two structure–pattern pairs are within the work bound, and both are dead. Both
  verdicts were re-checked with the integer test of §10:
  - n = 1981, H = ⟨4⟩ (|H| = 141, 21 orbits), rows (1, 1, 1, 89): the ℤ₇ compressions cannot
    be complementary;
  - n = 933, H = ⟨20⟩ (|H| = 62, 22 orbits), rows (1, 1, 3, 61): dead at d = 3.

So no Goethals–Seidel multiplier structure with at most 40 orbits that can be searched exactly
holds an unresolved Hadamard order. Larger structures are the FRONTIER regime (≥ 10¹⁴·⁵ operations).

## 12. Skew orders over non-cyclic groups

For p² | v, ℤ_v can be replaced by ℤ_p × ℤ_p × ℤ_{v/p²}. Its automorphism group contains GL(2, p),
whose odd-order subgroups (Singer, split torus, unipotent, Borel) give orbit structures that ℤ_v
does not have (`skew_group.py`). A subgroup of odd order always admits a skew first block. The
Goethals–Seidel array works over any abelian group, and `verify.py` checks it under the construction
`goethals-seidel-group`.

- **Survey** of the 26 open v ≤ 999 with a square factor (10–60 orbits):
  - The only structures with E near 1 are 275 (Singer-3 × ⟨3⟩, E ≈ 10⁰·⁴, 13 classes) and
    423 (E ≈ 10⁻¹·³).
  - Structures containing the unipotent group get model estimates of 10⁻⁴⁰ or below.
- **275 is solved:** skew-Hadamard order 1100, 32 solutions in four minutes of k-tree matching on
  the complete pools.
- **Cyclotomy over GF(pᵏ)** (`skew_field.py`). For v = q = pᵏ, the additive group of GF(q) with H an
  odd-order subgroup of GF(q)* is the field analogue of a multiplier structure.
  - **343 = 7³:** H of order 19 has 19 orbits, 9 classes and pair work 2·10⁹. The exact matcher finds
    4,211,424 solutions, giving skew-Hadamard order 1372 (Table 5: 343(4)).
  - **625 = 5⁴ and 729 = 3⁶:** the small cyclotomic structures (orders 39 and 13, and 91 and 13)
    have no admissible row sums.

## 13. Hadamard orders from the skew orders (Yamada's Theorem 4, after Kiyasu)

Yamada, J. Austral. Math. Soc. (A) 46 (1989) 371–383, Theorem 4: if q ≡ 5 (mod 8) is a prime power and a
skew-Hadamard matrix of order (q + 3)/2 exists, then a Hadamard matrix of order 4(q + 2) exists.
`verify/verify.py` has the `yamada` builder, and `search/make_yamada.py` nests one of our skew certificates in it,
after k doublings. For a skew order s (≡ 0 mod 4), q = 2s − 3 is automatically 5 (mod 8), so only its
primality matters.

| skew order s | q = 2s − 3 | Hadamard order 4(q + 2) | status of that order (October 2026) |
|---|---|---|---|
| 836 | 1669 | 6684 (n = 1671) | known; n ≤ 2999 not in the open list; built as a pipeline test |
| 1036, 1100, 1148, 1372 | 2069, 2197 = 13³, 2293, 2741 | 8284, 8796, 9180, 10972 | known (n ≤ 2999, not open) |
| 1564 | 3125 = 5⁵ | 12508 (n = 3127 = 53 · 59) | known: T-sequences of length 59 × Williamson-type of order 53 (Handbook "."); built as a test of the prime-power case |
| 1612 | 3221 | 12892 | known: T-matrices of order 11 × Williamson-type of order 293 |
| 1852 (Đoković 1994) | 3701 | 14812 | known: T(7) × Williamson-type 529; also Kiyasu from Đoković's 1852 |
| 2200 = 2·1100 | 4397 | 17596 | known: T-sequences of length 53 × Williamson-type 83 |
| 3016 = 2·1508 | 6029 | 24124 | implied since 2026-08-30 by T-matrices of order 163 (JD-Jones-ASES/Hadamard-T) × Williamson 37, but stated nowhere we found |
| 3556 | 7109 | 28444 | known: Paley II (q = 14221) |
| 3676 (Schumer 2026) | 7349 | **29404** (n = 7351) | implied since 2026-09-17 by Schumer's skew 3676, stated nowhere we found; Handbook Table 1.53 marks 7351 ".", apparently in error (every twin-prime n ≡ 3 mod 4 without a Paley route between 4000 and 9999 is marked "."); built from our skew 3676: [29404](../hadamard/29404) |
| 4808 = 2·2404 | 9613 | 38460 | known: Paley I (q = 38459) |
| **6448 = 4·1612** | **12893** | **51580 (n = 12895 = 5·2579)** | **probably open**: Paley I/II fail (51579, 25789 composite); no Williamson-type family reaches 2579 or 12895; Seberry–Yamada 1992 Table A.2 gives only 8·2579; n > 10000 is beyond Handbook Table 1.53; no 2026 deposit mentions it |
| 7112, 7408, 14816 | 14221, 14813, 29629 | 56892, 59260, 118524 | known (Paley I, Paley II, Williamson-type routes) |

So the new skew order 1612 (4·403), doubled twice, gives the Hadamard order 51580, which is probably new, and
Schumer's (and our) skew order 3676 gives 29404, whose Handbook entry we could not trace.

**Prime powers.** `build_yamada` and `make_yamada.py` also accept q = p^k: GF(q) is given by a primitive polynomial
over GF(p) in the certificate (fields `p`, `poly`), and r and ξ are encoded by base-p digits. Checked on q = 13, 29
and 5³ (orders 60, 124, 508) and on q = 5⁵ with skew 1564 (order 12508, H Hᵀ = 12508 I, 42 s).

`verify/verify.py 51580` was run separately (2026-10-12, 40 minutes, about 3.5 GB): `51580 ok`, with the recorded
digest.

## 14. The census of undecided skew orders: leads and calibration

The Đoković–Seberry census lists 737 odd v < 2500 with no skew-Hadamard matrix of order 4v (May 2025;
`search/data/skew_census_undecided.txt`). Leaving out the 13 orders of this repository and Schumer's 919,
723 remain: 233 primes and 490 composites.

**Spin structures on the centred hexagonal census primes.** A skew block has row sum 1, so the spin
structure needs 4v = 1 + 3b². Then |H| is forced by b ≡ ±1 (mod |H|). The census primes of this form above
1000 all have one spin structure with about 10⁵ predicted solutions (`spin_skew.py rank`):

| v | 4v | \|H\| | μ | H-orbits | K-classes | log₁₀ E |
|---:|---:|---:|---:|---:|---:|---:|
| 1657 | 6628 | 23 | 27 | 73 | 12 | 5.3 |
| 1801 | 7204 | 25 | 21 | 73 | 12 | 5.1 |
| 1951 | 7804 | 25 | 13 | 79 | 13 | 5.2 |
| 2269 | 9076 | 27 | 59 | 85 | 14 | 5.1 |
| 2437 | 9748 | 29 | 31 | 85 | 14 | 4.8 |

The composite ones (1141, 1261, 1387, 2107) have no structure with E > 10⁻¹. A literature check on
2026-10-12 found none of these five orders: Đoković (arXiv:2404.14375v3) lists spin families only up to
v = 631, and Schumer's 919 is the largest published.

The obstacle is the search, not the count. Restarts (of 3000 steps) that the joint tabu search `spin_ls.c` needs
per solution, measured over several seeds:

| v | K-classes | restarts per solution |
|---:|---:|---|
| 397 | 6 | 1, 2, 1 |
| 631 | 7 | 52, 30, 38, 46, 89, 45 |
| 919 | 9 | 4049, 5885, 6628, 2244 (and 261, 399 from csearch) |

That is roughly a factor 8 per class, or about 2·10⁶ restarts for 12 classes. At about 3 restarts/s per core
this is a week per order. The runs on 1657 and 1801 (16,000 restarts each) and 1951 (10,000) found nothing.
Their best residual had six (1657) and eight (1801) classes off by ±4. The runs were stopped in favour of
the four-list searches below.

Restarting half the time from the best state so far, perturbed by six swaps, cut the restarts on 631 to
10, 8, 29, 15, 2 and 103. That is a factor 2–4, which does not change the conclusion.
Two alternatives did worse on 919:
- a B-only tabu search scored against the nearest of all 2048 skew A: no solution in 2 × 400 s;
- a CP-SAT model with the products of signs as XNOR literals (OR-Tools 9.12): unknown after 600 s on 2
  workers.

**Census primes ranked for the four-list search** (`skew_rank.py`, lists of 1.2·10⁷): from 1009 to 2473
only 1123 (|H| = 33, rows 1, 1, 1, 67) and 1483 (|H| = 39, rows 1, 1, 1, 77) have admissible structures.
Below 1000 the best is 613 (|H| = 17, rows 1, 1, 35, 35). The rank's low-energy column (0.21 per class)
is optimistic for these.

**Calibration of `skew_lowE.py`.** The script prints the Gaussian estimate of the number of solutions
inside its actual lists, using their class means and variances.
- On 353 (low-energy lists of 4.7·10⁵ and full A) it estimated 10^1.7, and the k-tree found 16 solutions in
  80 s of matching. So the estimate is good to a factor of a few when v is prime.
- With lists of about 10⁶, 10⁶ and 2·10⁶, the estimates are 10^−1.3 for 613 and 10^−2.7 for 1123. Neither run
  is worth finishing at that size.
- `skew_gen` saturates: of 4·10⁶ requested minima only 2.1·10⁶ are distinct. Larger lists need shallower
  descents, which trade spread for diversity.

**Szekeres difference sets** (`szekeres.py`). A Legendre pair (a, b) on a group of order w, with a of skew
type and b symmetric (both with row sum 1), gives a skew-Hadamard matrix of order 2(w + 1) through a bordered
two-block array. The script checks H Hᵀ = (2w + 2)I and H + Hᵀ = 2I.
- Validation: with H the fourth powers it finds Szekeres' cyclotomic sets for w = 13, 29 and 53. It finds
  none for 37 and 61, where q = s² + 4t² has t = ±3.
- The census through this route: w = 2u − 1 for every undecided u, cyclic groups with odd multiplier
  subgroups, and the split GL(2, p) structures. All 988 structures with at most 18 sign pairs give nothing.
  767 have no symmetric b with row sum 1, and in the other 221 the exact join is empty.
- The three Kiyasu routes are covered with up to 22 pairs. They need skew 596 (Z₂₉₇, for 9532 after a
  doubling), 940 (Z₄₆₉, for 7516) and 1108 (Z₅₅₃, for 8860).
- The generic count is against this family: two free signs per class equation leave (4P)^C ≪ 1. Only
  algebraic coincidences could appear, and none did in these structures.

**Composite census orders: exact screen, then the chain.**
- `census_screen.py` lists, for each composite census v, every multiplier structure with 8–50 orbits and
  every row-sum pattern. It applies the exact compression test at each prime divisor. 903 of 1643
  structure-patterns survive, spread over 131 values of v.
- `screen_chain.py` runs `chain.py` on the survivors, ordered by the model estimate. That model is
  unreliable for composite v, so it only sets the order.
- **763 = 7 · 109 is solved:** H = ⟨9⟩ of order 27, rows (1, 1, 5, 55). Of the 40 level-1 branches, the
  first held 31 second-level combinations. Their final joins counted 3 solutions, and two were verified as
  skew-Hadamard matrices: skew order 3052, see [skew/3052](../skew/3052). Table 5 lists 763(4), and the
  census lists 763 as undecided.
- The model predicted 10^5.5 solutions for this structure, so it is off by orders of magnitude here as
  well. The sampled k-tree (400 values of τ, §10) had missed the order.
- Through Yamada's theorem, skew 3052 gives q = 6101, which is prime and ≡ 5 (mod 8), so a Hadamard
  matrix of order 4 · 6103 = 24412 exists. It is **probably new**, see [24412](../hadamard/24412):
  - Handbook Table 1.53 lists 6103 with t = 3. Paley I gives 8 · 6103.
  - No later update or 2026 deposit mentions it.
  - Cooper–Wallis would need Williamson-type or T-matrices of order 359.

  Its doubled skew input, 6104, gives q = 12205, which is not a prime power. The next one, 12208, gives
  q = 24413, prime, and so H(97660).
- **703 = 19 · 37 is solved**, in the breadth-first pass, which tries at most two level-1 branches per
  structure.
  - Structure: H = ⟨16, 26⟩ of order 27, a non-cyclic subgroup of ℤ₇₀₃* ≅ ℤ₁₈ × ℤ₃₆, with rows (1, 1, 1, 53).
  - Branches: 704 level-1 branches. The heaviest held 49 combinations, which counted 3 solutions in
    47 minutes.
  - This gives skew order 2812; see [skew/2812](../skew/2812). Table 5 lists 703(3).
  - No Yamada consequence: 2s − 3 is composite for s = 2812, 5624, 11248, 22496 and 44992.
- In both hits, the solutions sat in the heaviest level-1 branch. That is why the pass tries branches
  breadth first, and within a branch the largest second-level combinations first.
- 253 = 11 · 23 with ⟨12⟩ (33 orbits): the heaviest branches are truncated at the fibre cap
  (`--maxfibre 3e6`), so the chain cannot search it. The low-energy k-tree solved it (§15).
- **613 (prime) is solved** by the unequal-list low-energy k-tree: skew order 2452, see [skew/2452](../skew/2452).
  - Structure: H = ⟨37⟩ of order 17, rows (1, 1, 35, 35).
  - Lists: all 262,144 skew A; 1.29·10⁶ low-energy B (row sum 1); 7.1·10⁶ low-energy C = D (row sum 35, one
    shared list). The estimate inside the lists was 10^1.65.
  - The matcher found 4 solutions in 3.6 hours (match4 stops at 4).
  - The first 40,000 values of τ hold only about 4% of the τ-weight, which predicts about 2. In this
    regime the k-tree needs about 4 hours of work per solution.
  - An earlier run with smaller and deeper lists (estimate 10^−1.3) found nothing in 3000 values of τ. List
    size and diversity matter more than descent depth.
  - Yamada: the skew order 8 · 2452 = 19616 gives q = 39229, prime, and so H(156924). That order lies far
    beyond every table and was not pursued.
- **Single-level structures** were a blind spot of the chain. For 623 = 7 · 89 and 451 = 11 · 41, the joint
  DP over both primes does not fit in memory (the 451 joint DP was killed at 15 GB). The chain then used one
  level, and each branch ended in a single find4b over fibres of up to 3·10⁶, which never finishes.
  - `split_join` now groups the enumerated fibres of each branch by their compressions to the primes the
    joint DP left out. It joins those compressed autocorrelations exactly and runs find4b only on the
    consistent sub-lists.
  - It is exact: on 209/⟨26⟩, rows (9, 15, 23), the one-level chain plus the on-the-fly level-19 split counts
    the same 288 solutions as the two-level chain.
  - The per-prime screen also uses the full state cap now. Intermediate DP states can exceed the final count
    a hundredfold; for 451 at level 41 there are 2·10⁵ final states.
  - `split_join` turned out not to be needed: the exhaustive joiner `search/join4x.c` makes the single-level
    final joins cheap (next item), and the extra level is now off by default (`--split`).
- **The exhaustive joiner `join4x.c`.** It partitions both pair sides exactly on the three classes with the
  largest spread and probes small open-addressing tables keyed by a 64-bit hash of all classes. Partition keys
  are processed in lexicographic chunks that bound the build side; since lexicographic order is translation
  invariant, the probe cells that complete a chunk form an interval, found by binary search. The build phase
  inserts in parallel with compare-and-swap. It agrees with `find4b` on random instances and reproduces the
  288 solutions of 209/⟨26⟩. On the heaviest 623 branch (fibres 1,120, 3.4·10⁵, 2.7·10⁵, 2.2·10⁵) it
  finishes in 103 s; `find4b` had timed out at 900 s.
- **Exhaustive results with the new joiner** (`python3 search/chain.py v gens rows --budget 30000`):
  - 623 = 7 · 89, ⟨2⟩, rows (1, 5, 21, 45): **ruled out** (16 level-1 branches, 16 combinations, 0 solutions,
    28 minutes). The model had predicted 10^1.6 solutions.
  - 451 = 11 · 41, ⟨16, 42⟩, rows (1, 1, 29, 31): **ruled out** (20 branches, 12 minutes).
  - 451, ⟨16, 42⟩, rows (1, 1, 11, 41): **found** in the fifth of 60 level-1 branches, 16 solutions, 4 minutes:
    skew order 1804, see [skew/1804](../skew/1804).

## 15. Paley-product structures, and where the Gaussian model fails

**253 = 11 · 23 fell at once.** `skew_lowE.py 253 12 5,19,25` (seed 41) built all 65,536 skew A and 406,765,
506,562 and 872,916 low-energy candidates with row sums 5, 19 and 25 (`skew_gen` saturated far below the
requested sizes). The Gaussian estimate inside these lists was 10^−1.8. `match4` found 96 solutions at the
first value of τ, after 1.7 s, and all were verified: skew order 1012, see [skew/1012](../skew/1012).
Table 5 lists 253(4), SageMath lists 1012 as unknown, and the census lists 253 as undecided.

**Calibration with exact counts** (`join4x` counts every solution inside the lists):

| lists (A, B, C, D) | exact solutions | diagonal model | full covariance |
|---|---:|---:|---:|
| 253/⟨12⟩ (5, 19, 25): 65,536, 33,230, 42,858, 42,429 | 19,320 (10^4.29) | 10^−5.0 | 10^−2.4 |
| same, distinct class vectors only: 32,736, 14,362, 20,045, 20,467 | 245 | 10^−6.4 | 10^−4.1 |
| 253/⟨12⟩ (1, 13, 29), class vectors 32,736, 22,113, 3,251, 8,032 | 0 | 10^−7.6 | 10^−5.3 |

- The models of `skew_rank.py` and `skew_lowE.py` multiply per-class densities. For 253 the 16 classes are
  strongly correlated: the covariance of a block's class values has eigenvalues from 31 to 4,858, against 494
  to 633 for the prime 613. A full-covariance Gaussian gains two to three orders and is still off by six.
- Solutions also come in clumps. H fixes ℤ₁₁ × {0} pointwise, so B, C and D can each be translated by it
  without changing their autocorrelations: a solution has up to 11³ translates.
- For composite v the model fails in both directions: it gave 10^1.6 for 623/⟨2⟩, which has no solution.
  It should only order a list of structures; every structure deserves a cheap exact test.

**The mechanism: Paley products.** Let v = m q with q ≡ 3 (mod 4) prime, q ∤ m, and H = {h : h ≡ 1 (mod m),
h a square mod q}, of odd order (q − 1)/2. Then H has 3m orbits {(a, 0)}, {a} × QR, {a} × NQR. With χ the
Legendre sequence of ℤ_q and J the all-ones sequence, every block is

    X = u ⊗ 1 + w ⊗ J + z ⊗ χ,    u, w, z sequences on ℤ_m,

where at each a either w(a) = 0 and u(a), z(a) = ±1 (a z-position), or z(a) = 0, w(a) = ±1 and
u(a) ∈ {0, −2w(a)} (a w-position). Using J² = qJ, Jχ = 0 and χχ* = q·1 − J, the condition ΣX_kX_k* = 4v
splits into three equations over ℤ[ℤ_m]:

    (I)   Σ u_k u_k* + q Σ z_k z_k*                   = 4mq · 1
    (II)  Σ (u_k w_k* + w_k u_k*) + q Σ w_k w_k*       = Σ z_k z_k*
    (III) Σ (z_k u_k* − u_k z_k*)                       = 0

- **A counting identity.** Let N_w be the number of w-positions in all four blocks and N₋ the number with
  u = −2w. The coefficient of the identity in (II) gives (q + 1)N_w − 4N₋ = 4m with 0 ≤ N₋ ≤ N_w. So there
  is **no solution** unless some integer N_w satisfies it. In particular there is none when q > 4m + 3. This
  excludes every census value v = 3q, 5q, 7q, 9q with a prime q ≡ 3 (mod 4).
- **Census values that pass** (`paley_product.structures`, all factorisations): 119 (N_w = 9–17), 253 (2),
  285 (3), 345 (3), 391 (3), 465 (2), 473 (1), 475 (5), 517 (1), 575 (5), 621 (5), 667 (5), 899 (4), 903 (2),
  1081 (2), 1139 (1), 1207 (1), 1209 (5), 1333 (3), 1501 (1), 1645 (3), 1711 (2), 1763 (4), 2193 (5),
  2279 (5), 2365 (5), 2485 (2).
- **N_w = 1** (q = 4m − 1 or 4m + 3). A skew A has its w-positions in pairs ±a, so the single one lies in
  another block j, translated to a = 0. Since |Σ PAF_u(s)| ≤ 4m < 2q, u_j must be skew-type away from 0,
  and then (u_k) and (z_k) are complementary quadruples and (III) holds. The row-sum identity of the
  compressions gives r_B² + r_C² = 4m − 2, so **2m − 1 must be a sum of two squares**. For m = 11 (21) and
  m = 17 (33) it is not. Hence **473, 517, 1139 and 1207 have no solution in this structure**. Only 1501
  (m = 19, 37 = 1 + 36) remains among the N_w = 1 cases.
- **N_w = 2 in two different blocks** forces q = 2m + 1 and N₋ = 1. The same bound gives complementary
  (u_k) and (z_k), with u₁(0) = 0, u₂(0) = −2 and u₁(s) + u₁(−s) = −(u₂(s) + u₂(−s)).
- In both cases q enters only through the counting identity: the problem is a search over ℤ_m.

- **N_w = 3.** Two w-positions in one block (or a pair ±a in A) put a term of size about q² into the autocorrelation
  of the compressions y_k = u_k + q w_k, which only another such pair could cancel. So for N_w ≤ 3 each w-position
  sits alone in its own non-A block. Without any decoupling assumption the system is then: stage 1, the compressions
  y_k (±1 with one entry q or q − 2 at 0 in each w-block) have Σ PAF_y = 0; stage 2, the z_k have Σ PAF_z(s) = R(s) =
  Σ_{w-blocks} (u_k(s) + u_k(−s)) and zero cross terms. This general search reproduces the u-quadruple counts of the
  decoupled N_w = 2 search exactly (9,520 and 209,536).
- **Precedent.** The counting identity is classical in another setting. Seberry Wallis (J. Combin. Theory A, 1976),
  Lemma 8, fills an orthogonal design of order 2ᵗ with J, J − 2I and a Paley-core block under a(v + 1) + b(v − 3) =
  2ᵗ, which is our identity with a = N_w − N₋ and b = N₋. Seberry (Ars Combin., 1978), Theorems 6 and 7, give skew
  orders 2q(q − 1) and mn(n − 1) from orthogonal designs, and Wallis (1969, 1970) uses the same J, J − 2I, Paley-core
  blocks. In every version we found, the outer structure is an orthogonal design or a Kronecker product. Here it is a
  Goethals–Seidel quadruple over ℤ_m, and no located theorem produces the orders below (literature search of
  2026-10-12).
- **Prime powers.** Nothing above uses that q is prime: with χ the quadratic character of GF(q), q = pᵏ ≡ 3 (mod 4),
  the group becomes ℤ_m × GF(q)⁺ (non-cyclic for k > 1). The census has one such case, 513 = 19 · 27 (N_w = 3, N₋ =
  2); `--field 3,3` lifts to ℤ₁₉ × ℤ₃³ and verifies with the `goethals-seidel-group` builder.

**`search/paley_config.py`** solves the reduced problem exactly. Stage 1 lists the u-quadruples (or, in the general
search, the y-quadruples) by an exact join. Stage 2 joins, for each of them, all z_A (symmetric), the z_k with a zero
at each w-position and the free z_k on the PAF_z and cross classes with `join4x`. Each hit is lifted to ℤ_v and
verified.

| v = m · q | N_w | N₋ | stage-1 quadruples | z-solutions of the first | time | result |
|---|---:|---:|---:|---:|---:|---|
| 253 = 11 · 23 | 2 | 1 | 9,520 | 100 | 1 s | skew 1012 again |
| 465 = 15 · 31 | 2 | 1 | 209,536 | 28 | 16 s | **skew 1860**, see [skew/1860](../skew/1860) |
| 285 = 15 · 19 | 3 | 0 | 2,046,720 | 28 | 20 s | **skew 1140**, see [skew/1140](../skew/1140) |
| 345 = 15 · 23 | 3 | 3 | 2,040,192 | 22 | 21 s | **skew 1380**, see [skew/1380](../skew/1380) |
| 391 = 17 · 23 | 3 | 1 | 1,620,096 | 18 | 252 s | **skew 1564**, see [skew/1564](../skew/1564) |
| 513 = 19 · 27 (ℤ₁₉ × GF(27)) | 3 | 2 | 12,124,320 | 6 | about 60 min | **skew 2052**, see [skew/2052](../skew/2052) |
| 1501 = 19 · 79 | 1 | 1 | 95,552,064 | 24 | about 40 min | **skew 6004**, see [skew/6004](../skew/6004) |

- **1140 was within classical reach.** Seberry's theorem (Ars Combin. 6 (1978), Thm 7: an OD(mn; 1, m, mn − m − 1) and
  amicable designs of order n give a skew order mn(n − 1)) produces it with m = 3, n = 20 from an OD(60; 1, 3, 56),
  which a Goethals–Seidel search over ℤ₁₅ finds in 0.1 s; we rebuilt that matrix and it verifies. The design was not
  recorded anywhere we looked, and every list had 1140 as unknown. In our terms this is the case where (q + 1)/4
  divides m (285 = 15 · 19 with (19 + 1)/4 = 5). For 1012, 1380, 1564, 1860 and 6004 it does not divide m, and none of
  Seberry's theorems applies.
- **Open lead: 475 = 25 · 19** has the same classical route through an OD(100; 1, 5, 94). The Goethals–Seidel layout
  "a at A(0); b at {0, 5} in B and C and at 0 in D; c elsewhere" has no solution invariant under the multipliers of
  order 5; without that symmetry the exact join (about 7·10¹³ probe pairs) is out of reach.

Over 40 random u-quadruples of 253 every one extends (mean 82 z-solutions, minimum above 0, maximum 124); over 12 of
465 the mean is 20 (minimum 4, maximum 42). The counts are stable, so the first quadruple is representative.

**Extra multiplier symmetry K ≤ ℤ_m*** (`--kgen`, u and z constant on K-orbits) shrinks the search but has found
nothing: no stage-1 quadruple for 903/⟨4⟩, 1081/order 11, 1711/order 7 or 1501/order 9; four for 1333/⟨5⟩ and four
for 1333/⟨2⟩, none of which extends; and no extension among the first 20 of 4,992 (1645/⟨11⟩) or the first 50 of
2,560 (2485/⟨11⟩).

The measured means (82 at m = 11, 20 at m = 15) fall by about 0.7 per unit of m. That predicts about 2.5
z-solutions per u-quadruple for m = 21 (903), 1.2 for m = 23 (1081), 0.15 for m = 29 (1711) and 0.02 for m = 35
(2485). Stage 2 costs about 2^{2m−2} probes per u-quadruple, so 903 and 1081 are in reach and the larger two need a
better stage 2.

**Paired w-positions (N_w = 4, 5).** Two w-positions in one block put y(a)y(b) ≈ ±q² into PAF_y at their distance,
so pairs must cancel: {0, d} in two blocks with opposite products, or {a, −a} in A against {0, 2a} in another
block; a fifth w-position sits alone. `paley_config.py` enumerates these layouts. With an extra multiplier group
fixing the pair positions (order 5 on ℤ₂₅ for 475 and 575, order 3 on ℤ₂₇ for 621) no layout has a stage-1
quadruple. For 1209 (m = 39, ⟨16⟩ of order 3, pairs at distance 13) six layouts have 544 or 1,152 stage-1
quadruples, none of the first 50 of which extends.

**Quartic analogue (lead, not searched).** For q ≡ 5 (mod 8), H = {h ≡ 1 (mod m), h a fourth power mod q} also
avoids −1. The same compression identity Σ_k Σ_a y_k(a)² = 4mq gives, with e = (q − 1)/4,
a₂(q + 3)/4 + b₂(q − 5)/4 + a₄(q + 1) + b₄(q − 3) = 4m for the numbers of positions where the four class values
sum to ±2 or ±4. For q = 53 it rules out m = 5 and 11 (265, 583) and allows m = 13 and 15 (689, 795); 112 census
values pass it. Single-big-position cases need r_B² + r_C² = 4m − 2 again: 1903 = 11 · 173 fails (42), 1341 = 9 · 149
passes (34 = 3² + 5²) with the single row pattern (1, 3, 5, 73).
  - **The equations.** With C₀..C₃ the classes, χ = C₀ − C₁ + C₂ − C₃ (symmetric), κ = C₀ − C₂, λ = C₁ − C₃
    (antisymmetric) and ψ = κ + iλ the quartic character, a block is X = aδ + wG + bχ + Dκ + Eλ with a = f − s/4,
    w = s/4, b = (g₀ − g₁ + g₂ − g₃)/4, D = (g₀ − g₂)/2, E = (g₁ − g₃)/2. The products close with Jacobi sums
    J = J(ψ, ψ), J′ = J(χ, ψ) (Gaussian integers of norm q): ψψ = Jχ, ψψ̄ = G − qδ, χψ = J′ψ̄, χ² = qδ − G. The
    condition becomes five equations over ℤ[ℤ_m]:
    δ: Σ aa* + q bb* + (q/2)(DD* + EE*) = 4mq; G: Σ aw* + wa* + q ww* − bb* − (DD* + EE*)/2 = 0;
    χ: Σ ab* + ba* − (Re J/2)(DD* − EE*) − (Im J/2)(DE* + ED*) = 0;
    κ: Σ (Da* − aD*) + Re J′(Db* − bD*) + Im J′(Eb* − bE*) = 0;
    λ: Σ (Ea* − aE*) + Im J′(Db* − bD*) − Re J′(Eb* − bE*) = 0.
    δ + q·G is again the compression condition. A position with s = 0 is either χ-type (b = ±1, D = E = 0) or
    ψ-type (b = 0, D, E = ±1); both add 2 + 2q to the δ identity coefficient. Skewness forces A to be ψ-type at
    a = 0. If that were the only ψ-type position in all four blocks, the κ-equation would need
    |±Re J′ ± Im J′| = 1, which is impossible, so ψ-type positions occur elsewhere too.
  - Stage 2 has six channel patterns per position against two in the Paley case. A direct low-energy test on 1341
    (O = 45) needed 8 GB in the exact join and was stopped; no dedicated search has been written.


## 16. Cyclotomic products of higher index, and the order 3732

**Structure.** Let v = m q with q a prime not dividing m, e | q − 1 and F = (q − 1)/e odd, so −1 is not an e-th power.
H = {h : h ≡ 1 (mod m), h mod q an e-th power} has order F and avoids −1. Over ℤ_m × GF(q) an H-invariant block is
given at each a ∈ ℤ_m by f(a) (the value at (a, 0)) and x₀(a), …, x_{e−1}(a) (the values on (a, C_j), C_j the
cyclotomic classes of index e). e = 2 is the Paley product of §15. The compression to ℤ_m is y(a) = f(a) + F s(a),
s(a) = Σ_j x_j(a), and Σ_k Σ_a y_k(a)² = 4mq becomes

  Σ over the 4m positions of (F s² + 2 f s) = 4 m e.

A position with s ≠ 0 (a spike) costs at least 4F − 4. In the character basis a block is fδ + Σ_i ĉ_i Ψ_i
(Ψ_i the multiplicative characters of order dividing e, extended by 0), and the products Ψ_i Ψ_j are Jacobi sums
J(χ^i, χ^j) Ψ_{i+j} of absolute value √q. Unlike the Paley case there is no separation of the equations by size, so
the search uses the two compressions as exact intermediate levels instead (`search/cyclo_product.py`, `join3.c`):

1. stage 1: the compressions to ℤ_m, up to negation, translation and reversal of the blocks;
2. level B: the compression to GF(q) (summing over ℤ_m), whose autocorrelations on the classes C_t, t < e/2, are
   the group key;
3. level C: the autocorrelations at the shifts (b, C_c), 1 ≤ b ≤ (m − 1)/2; stage 1 + B + C is the full condition.

**Single-spike family.** One spike with s = ±2 and f = −sign(s) costs 4F − 4, which equals 4me iff F = m e + 1, i.e.
q = m e² + e + 1. For m = 3: q = 53 (e = 4, v = 159), 311 (e = 10, v = 933), 991 (e = 18, v = 2973), 1753 (e = 24,
v = 5259). For m = 9: 149 (e = 4, v = 1341), 331 (e = 6, v = 2979).

**933 = 3 · 311, e = 10 (order 3732).** The single spike leaves one compression quadruple up to symmetry, rows (1, 61,
3, 1). The lists have 8,064 (skew block), 13,335,840 (spike block) and 16,003,008 members; level B groups them into
512, 55,915, 57,892 and 57,917 groups. join3 found 24 solutions among the first 6,978 of about 6 · 10⁷ matching group
quadruples (5 minutes on 3 threads): of the order of 10⁵ skew solutions in all. A Gaussian model over the 14
independent equations conditioned on the compressions predicts 10^4.2 to 10^4.6 (skew) and 10^8.1 (plain). The
estimate in §5 for the same subgroup was 10^0.3 with 10¹⁷ pair operations. This settles the order 3732; see
[3732](../hadamard/3732).
  - Since 3733 is prime and ≡ 1 (mod 4), Miyamoto's theorem gives H(4 · 3733) = H(14932) from H(3732). The
    repository's Miyamoto builder rebuilds and checks it (H Hᵀ = 14932 I) in about 4 minutes.

**Other skew census values.** Exact stage-1 results (no compression quadruple at all, so no solution in the
structure): 201 = 3 · 67 (e = 6), 265 = 5 · 53 (e = 4), 305 = 5 · 61 (e = 4), 335 = 5 · 67 (e = 6), 515 = 5 · 103
(e = 6), 549 = 9 · 61 (e = 4), 573 = 3 · 191 (e = 10), 695 = 5 · 139 (e = 6), 707 = 7 · 101 (e = 4) and 909 = 9 · 101
(e = 4). 1341 = 9 · 149 (e = 4) has 36 stage-1 classes with a Gaussian estimate of about 1 solution each, but its
GF(149) level has only two dimensions, so the final joins would cost about 10¹⁴ each; not run. 261 = 9 · 29 (e = 4)
has 3.3 · 10⁶ compression options per block and was not finished.

**Plain orders.** Among the composite open n of FRONTIER.md, 933 is the only one whose cyclotomic product has lists
small enough to enumerate: log₂ of the per-block list is 24 for 933, 40 for 1255 = 5 · 251 (e = 10, two spikes) and 47
for 2973 = 3 · 991 (e = 18, a single spike). For Spence's bordered construction the compressed target off zero is
−4q while the row sums satisfy Σ r_k² = 4, so the spikes must cancel inside a block, and their products, of size
about (2F)², cannot give −4q at every shift. This rules out, for example, 947 (ℤ₂₂ × GF(43), e = 2) and 537
(ℤ₈ × GF(67), e = 6).

**A multiplier on ℤ_m: 1255 = 5 · 251, e = 10 (order 5020).** `--kgen k` adds the multipliers ⟨k⟩ of ℤ_m (acting on
the ℤ_m coordinate) to H. For plain structures K_m may contain −1. For 1255, H = {h ≡ ±1 (mod 5), h a 10th power mod
251} has order 50 and 33 orbits per block; this is the |H| = 50 structure that the ranking of §5 rated at log₁₀ E =
−5.0. With F = 25 a spike costs 96 or 104 at position 0 and twice that on the orbits {±1}, {±2}, so the identity (Σ =
200) forces one spike of each kind at position 0 of two blocks, and stage 1 has a single class (rows 1, −1, −53, −47).
The lists have 1.3–1.6 · 10⁷ members, but level B splits them into about 10⁶ groups, so a complete join is out of
reach. Keeping the 8,000 largest groups of every list (about 8 % of the members) gave 32 solutions in 17,273 group
quadruples (6 minutes); see [5020](../hadamard/5020). An FFT-based Gaussian estimate with the K_m symmetry gives log₁₀
E ≈ 1.5 for the whole structure, which the run shows to be far too low.

**Sweep of the plain structures for the open orders.** For every composite n of FRONTIER.md (47 orders) and every
product Z_m × GF(q), even index e with (q − 1)/e ≥ 3 and at most about eight spikes, and every cyclic multiplier
group K_m ≤ ℤ_m* (including those containing −1) whose lists have at most 2 · 10⁷ members per block, stage 1 was
run exactly: 239 structures. Only 1255 (e = 10, K₅ = {±1}) has a compression quadruple; 933 had been solved
already. Every other structure has no solution. Structures with larger lists (for instance 2973 = 3 · 991 with
e = 18, a single spike and 4.9 · 10⁴ patterns per position) are not covered.

**1081 = 23 · 47 (Paley, N_w = 2).** `paley_config.py 1081 --maxu 3` went through all 6¹¹ joint (u₁, u₂) options in
2,759 s: stage 1 is empty, so 1081 has no solution in the Paley product with N_w = 2 (the only admissible N_w).

**Skew sweep with a multiplier on ℤ_m (`search/skew_kmult_sweep.py`, 2026-10-13).** Every census value v = m q with q
prime, e | q − 1 even and (q − 1)/e odd, and every cyclic K_m ≤ ℤ_m* without −1 such that the structure has at most 22
position orbits and at most 2·10⁷ members per block: stage 1 run for 210 structures, 8 with stage-1 classes, all 8
joined (`cyclo_product.py`, `--sub 8000`, K = 1 per class). One hit: **v = 1417 = 13 · 109**, e = 4, K_m = ⟨3⟩, H =
⟨1405, 328⟩ of order 81: skew-Hadamard matrix of order 5668 ([skew/5668](../skew/5668)). The structure is tiny (lists
of 144, 7,776, 5,184 and 5,184 members, 12 stage-1 classes); an exhaustive join of all 12 classes gives 8 solutions in
each of the classes 2, 3, 6 and 11 and none in the others, 1 s per class. No hit for 1313 = 13 · 101 and 481 = 13 · 37
(e = 4, K_m = ⟨3⟩), 759, 261 and 1333.

**1341 = 9 · 149 (skew, e = 4), partial.**
`cyclo_product.py 1341 --q 149 --e 4 --config 0 --K 3 --threads 3 --bkey 3 --sub 8000` (level-B key extended by the
shifts (3, C_c)): lists of 5,184, 10,077,696, 10,077,696 and 6,718,464 members, of which the 8,000 largest groups kept
100 %, 13 %, 22 % and 3.3 %; 2,604,080 group quadruples, 2.0 · 10¹⁰ member operations, 1,961 s, no solution. This
covers about 10⁻³ of the configuration's product space and rules nothing out.

## 17. Symmetric Hadamard matrices from (A, B, B, D)

The Đoković–Seberry census also lists 595 odd v < 2500 for which no *symmetric* Hadamard matrix of order 4v is
known; the smallest are 65, 89, 101, 107, 119, 133, 149, 153 and 163, and only 12 values (29, 43, 47, 59, 67, 73,
93, 103, 109, 127, 151, 191) were settled by computer searches.

**Arrays.** A search over 4 × 4 arrays with blocks ±X, ±XR, ±Xᵀ, ±XᵀR (X ∈ {A, B, D}, C = B, R the back identity)
for arrays that are symmetric and orthogonal for every admissible quadruple gives two useful forms. With
AAᵀ + 2BBᵀ + DDᵀ = 4v I for circulants A, B, D:

- A and D symmetric: [[A, D, BR, BR], [D, −A, BR, −BR], [BR, BR, −D, −A], [BR, −BR, −A, D]];
- only A symmetric: [[A, B, B, D], [Bᵀ, DR, −A, −BR], [Bᵀ, −A, −DR, BR], [Dᵀ, −BR, BR, −A]].

Both are symmetric Hadamard matrices of order 4v. In sequence terms: PAF_A + 2 PAF_B + PAF_D = 4v δ with A symmetric
(and D symmetric for the first form). `search/propus.py` searches multiplier structures exhaustively: A on the orbits
of ⟨H, −1⟩, B and D on the orbits of H (or of ⟨H, −1⟩), every row-sum pattern r_A² + 2r_B² + r_D² = 4v, exact
three-list joins with join4x, and every hit is rebuilt as a matrix and checked (H = Hᵀ, H Hᵀ = 4v I).

**Checks.** v = 13: thousands of solutions in both forms (orders 52, all verified). v = 29: with A and D symmetric and
B arbitrary there is no cyclic solution at all (both row-sum patterns (3, 7, 3) and (7, 3, 7) exhausted, pools of 2–3
· 10³ symmetric and 3.5–6.8 · 10⁷ arbitrary sequences), so the known symmetric matrix of order 116 is not a cyclic
instance of the first form. It is an instance of the second: the family of Di Matteo, Đoković and Kotsireas (Spec.
Matrices 3 (2015) 227–234, arXiv:1503.04226) with −1 entries A = {0, 4, 5, 6, 7, 9, 13, 16, 20, 22, 23, 24, 25}, B = C
= {0, 1, 2, 5, 6, 8, 11, 12, 14, 16, 22}, D = {0, 1, 2, 3, 5, 8, 10, 13, 14, 15, 18, 22, 25} gives a symmetric
Hadamard matrix of order 116 in `symmetric_matrix_v3` (checked).

**Literature (survey of October 2026).** Every computer-found symmetric matrix of the census is of the second kind:
one non-repeated block symmetric, the repeated block B = C neither symmetric nor skew (except v = 43, where it is the
skew Paley set), the fourth block arbitrary (Seberry–Balonin's propus construction, Australas. J. Combin. 69 (2017),
and Balonin–Đoković's "propus difference families", arXiv:1803.08787). Both blocks A and D symmetric is known only for
v = 13 and Turyn's series (arXiv:1708.05098). Cyclic propus searches are exhaustive for odd v ≤ 63 only. Since the
census, Đoković (arXiv:2508.17141, August 2025) proved that a symmetric Hadamard matrix of order q(q + 1) exists for
every prime power q ≡ 3 (mod 8), which settles 4 · 189 (q = 27) and 4 · 473 (q = 43); the paper also builds q = 67,
which gives 4 · 1139, an order on the census list. Nothing settles 65, 89, 101, 107, 119, 133, 149, 153, 163 or 167;
the symmetric orders 668, 716, 892, … remain open (the matrices of Alpöge et al. are not symmetric).

**Sweep.** `propus.py v --sweep --maxo 24 --maxpool 2e7` runs, for every subgroup of ℤ_v* (cyclic subgroups and
products of two) whose structures have at most 24 orbits, all three variants (A, D symmetric with B symmetric or
not; only A symmetric). For v = 65, 89, 101, 107 and 119 (24, 6, 6, 2 and 17 subgroups) every row-sum pattern either
has an empty pool for some block (the orbit sizes cannot reach the row sum) or was joined completely: no symmetric
Hadamard matrix of order 4v in these structures. All pools fit in 2 · 10⁷, so nothing was skipped for size.

## 18. Bordered Goethals–Seidel arrays with a character set T

**Construction.** Let G be abelian, K ≤ G of index q, Q = G/K, and T a set of t characters of Q that is a union of
Galois orbits; put Γ = 4 Σ_{χ∈T} χ (an integer function on Q, Γ(0) = 4t). If four ±1 sequences on G satisfy
Σᵢ PAFᵢ(g) = −Γ(g + K) for g ≠ 0, a border of 4t rows and columns that are constant on the cosets of K completes
the Goethals–Seidel array developed over G to a Hadamard matrix of order 4n, n = |G| + t, provided the *border
problem* has a solution: X (Hadamard of order 4t), R̂ (4t × 4q, orthogonal rows), Ĉ (4q × 4t, Ĉ Ĉᵀ = I₄ ⊗ Γ) with
X Ĉᵀ + R̂ Mcᵀ = 0, where Mc is the Goethals–Seidel array over Q of the compressions yᵢ (sums over the cosets of K).
Equivalently the matrix [[X, √|K| R̂], [√|K| Ĉ, Mc]] of order 4t + 4q is orthogonal with Gram 4n I; the border
problem does not depend on |K| beyond |K| > 18 (which forces X and R̂ to be Hadamard-like). Consequences:
Σᵢ |ŷᵢ(χ)|² = 4t for χ ∈ T and 4n for χ ∉ T; rows of R̂ and columns of Ĉ have Fourier support in T, so for odd q
the trivial character must lie in T. Spence's construction is Q = 1, T = {1}; the matrices of Alpöge et al. and
the matrix of order 3964 use Q = ℤ₂ × ℤ₂ with T the three nontrivial characters (n ≡ 3 mod 4).

**The family Q = ℤ_q, T = all characters.** Here Γ = 4q δ_K: the targets are −4q on K \ 0 and 0 off K, the border has
size 4q and X, R̂, Ĉ are Hadamard of order 4q. For q = 3 it applies to every n divisible by 3 (|G| = n − 3), in
particular to n ≡ 1 (mod 4), where Alpöge's family does not exist. The compressions to ℤ₃ satisfy Σᵢ PAF_{yᵢ} = 12 δ;
for even |K| they are three single ±2 entries in three different blocks (4 · 27 · 8 = 864 patterns).
`search/bordered.py patterns` (C solver `search/border3.c`) finds a border for all 864 patterns (4 min).
`search/bordered.py validate` assembles complete examples: orders 36 (G = ℤ₃ × ℤ₂) and 60 (G = ℤ₃ × ℤ₄ and ℤ₃ × ℤ₂ ×
ℤ₂), each with H Hᵀ = nI. Negation on the ℤ₃ factor seems unusable as a symmetry: the structures tried have no
solution (n = 15, 21, 27, exhaustive) or no row pattern whose lattice classes contain the targets (n = 573, all single
generators).

**2148 = 4 · 537 in this family.** G = ℤ₅₃₄ = ℤ₆ × GF(89) is forced (t = 3). As a cyclotomic product (multiplier
classes of index 8, F = 11) the compression identity becomes Σ over the 24 positions of (11 s² + 2 f s) = 96: exactly
two spikes, ±23 at (0, 0) and ∓23 at (0, 1) of the block whose ℤ₃-compression is zero (`search/bgs.py`, one stage-1
class). No further multiplier is compatible: −1 on ℤ₃ forces equal compression values at ±1 and breaks the delta
condition, and twisted multipliers force s ∈ 4ℤ at the spikes. The lists then have 1.2 · 10¹¹ members per block, too
many for the joins of §16. (The diagonal Gaussian estimate of 10⁷ given here earlier ignores the strong correlations
between the classes; the full-covariance estimate is 10^13.2.) §19 solves this structure: order 2148. The other orders
537 − t with a usable quotient (t = 9: |G| = 528, t = 15: 522) give much smaller estimates. `bgs.py scan` (cyclic G =
ℤ_m × GF(r), m ≤ 30, t ≤ 30, at most four spikes) found no other order with a feasible stage 1 among the FRONTIER
orders up to n = 1423.

**Q = ℤ₂² × ℤ₃, T = (nontrivial on ℤ₂²) × (all on ℤ₃), t = 9.** This product of Alpöge's set and the ℤ₃ set applies to
every n ≡ 9 (mod 12), i.e. to all open n ≡ 1 (mod 4) divisible by 3, and allows Alpöge's order-3 automorphism α of ℤ₂²
coupled with multipliers of K. Gaussian estimates (`search/bgs_gauss.py`, summed over all row patterns): 789 (K = ℤ₅ ×
ℤ₁₃, σ = (α, ×2, ×4), 84 orbits, 47 classes) 10^8.1; 789 with (α, ×4, ×2) 10^3.4; 537 at best 10^−3.1 (G = ℤ₂⁴ × ℤ₃₃);
669, 1437, 2661 below 10^−10. The border problem becomes: two orthogonal bases of the 36-dimensional space V_T made of
"fibre-character" vectors (one signed character of ℤ₂² on each of the 12 fibres (block, a)), with all mutual inner
products ±2 (in units where the norm is 12) for the single-spike compressions, i.e. a pair of unbiased weighing
matrices W(36, 12) with this fibre structure. The twisted-product ansatz (rows h_j(f) e_{π_f(u)} from a Hadamard
matrix of order 12) cannot work: it would need four columns of a Hadamard matrix of order 12 with constant product. A
meet-in-the-middle search over all 6¹² candidate columns found none for random bases of this kind. Open.

## 19. The ℤ₂-split of a bordered cyclotomic product: order 2148

**Structure** (§18): G = ℤ₆ × GF(89) ≅ ℤ₅₃₄, K = {a ≡ 0 mod 3} × GF(89), Q = ℤ₃ with T = all characters (border of
size 12), sequences invariant under the eighth powers of GF(89)* (e = 8, F = 11), and the single stage-1 class
y = (1, 1, −1, 1, −1, 1), (1, −1, −1, −1, 1, −1), (1, 1, 1, 1, −1, −1), (−23, 1, 1, 23, −1, −1).

**The split.** ℤ₆ ≅ ℤ₂ × ℤ₃ pairs the positions a and a + 3. Put U = (x(a) + x(a + 3))/2 and W = (x(a) − x(a + 3))/2
on ℤ₃ × GF(89) (α = a mod 3). Both are 0/±1 with complementary supports, x is recovered as U ± W, and the bordered
condition is equivalent to

- Σ_k PAF_{U_k}(β, d) = −6 [β = 0] for (β, d) ≠ 0 (the compression of x to G/⟨3⟩), and
- Σ_k PAF_{W_k}(β, d) = 0 for (β, d) ≠ 0 (the part odd under a ↦ a + 3).

At (0, 0) the first reads Σ_k |supp U_k| = 531 = 3 + 11 · 48: exactly 48 of the 96 class cells are nonzero in the U's.
The other U-classes are (0, C_t) for t < 4, (1, 0) and (1, C_t) for t < 8, with two linear relations from stage 1 (12
independent equations); the W-condition has 11 independent classes, and its values lie on a parity lattice of index
2¹¹ once the supports are fixed. Given the U's, the W of each slice is a sign vector on the zero cells of U with a sum
fixed by stage 1 (C(4, 2) = 6 choices for four zero cells), so each U-solution leaves a small exact four-list problem
of about 200 vectors per block. This is the four-list split of T-matrices (§6) applied inside one block.

**Calibration on n = 177** (G = ℤ₆ × GF(29), e = 4, F = 7: the same stage-1 pattern with spikes ±13):
- the structure holds exactly 26,095,424 solutions in its stage-1 class (`join3.c` on the full lists of 46,656 and
  20,736 members, 46 minutes); the full-covariance Gaussian estimate is 10^6.87, the diagonal one 10^4.17;
- 157,482,016 U-solutions (`join4x.c`), with 0.166 W-completions each on average;
- reconstructed solutions give Hadamard matrices of order 708 (H Hᵀ = 708 I, border from `border3.c`).

**The search for 537** (`search/bgs_uw.py`, `search/kt4.c`):
1. Each block's U-vectors with exactly 12 nonzero class cells: 1.55 · 10⁸ (all U-vectors: 1.36 · 10⁹, 1.27 · 10⁹ for
   the spike block). Lists of 10⁷ random members per block, keyed by the 14 U-classes.
2. `kt4.c`: a four-list k-tree merging on five of the classes (1, C_t), t = 0..4, at the split values nearest the
   centre. Merging on four classes made 7 · 10⁸ pairs per split value and ran out of memory; five give 9 · 10⁶. One
   split value yields about 142,000 U-solutions in 8 s, sixty times the Gaussian estimate (10^14.7 U-solutions with 12
   nonzero class cells per block).
3. W-completion: an exact join of the four W-lists per U-solution (`UW.wjoin_fast`, about 9 ms; it reproduces
   `UW.wjoin` exactly on 60 U-solutions each for n = 177 and 537), four processes.
4. The Gaussian estimate of the W-completion probability is 10^−8.7 per U-solution, or 10^−5.4 with the parity
   lattice; on n = 177 the same estimate gives 10^−2.2 against an observed 10^−0.68. The first 1,000 U-solutions
   (seed 1) gave no completion.
5. `python3 search/bgs_uw.py search 89 8 --L 10000000 --sigmas 2 --maxu 300000 --threads 4 --seed 2`: 284,292
   U-solutions, **one completion**, verified as a Hadamard matrix of order 2148 (197 s on 4 cores); see
   [2148](../hadamard/2148). The observed rate 3.5 · 10⁻⁶ matches the estimate with the parity lattice.

So the structure that §18 judged out of reach takes minutes once the ℤ₂-split separates the conditions: the U-problem
has 12 equations and very many solutions, and each U-solution leaves a W-problem small enough to solve exactly.

**Generalization** (`UWG` in `search/bgs_uw.py`). For any even m = 2m′ the pairs (α, α + m′) give
PAF_x(β + j m′, d) = 2 (PAF_U(β, d) + (−1)^j PAF_W(β, d)), with PAF_U cyclic and PAF_W negacyclic over ℤ_m′, so the
targets of U and W are (t(β) ± t(β + m′))/4. A multiplier group K_m = ⟨k⟩ on ℤ_m maps pairs to pairs (k is odd), with
W changing sign when a slice lands on the second position of its image pair; an orbit that maps a slice to itself with
a flip forces W = 0 there. Validation: on small structures with K_m (n = 24 with ℤ₄ × GF(5), n = 30 with ℤ₄ × GF(7))
the U/W counts equal exhaustive counts over all x-quadruples (36 and 128), and for n = 177 the class `UWG` reproduces
the completion rates of `UW`.

**Other open orders in this shape.** `python3 search/bgs_uw.py scan <the 148 n of FRONTIER.md>` looks at cyclic
ℤ_m × GF(r) with m even, K_m trivial, t ≤ 40, F ≥ 3 and at most six spikes. It finds nonempty stage 1 only for 537,
1367 and 2027 (m = 2, Spence-type, 44 to 62 classes per fibre) and 1689 = 6 · 281 + 3 (the 537 pattern with e = 14 or
20). It does not cover about 4,800 structures (145 orders) with more than 22 position orbits or more than 4 · 10⁵
stage-1 options per block. For 1689, F = (281 − 1)/e is even,
so −1 lies in C₀, the classes are no longer paired by negation, and the condition has 55 (e = 14) or 79 (e = 20)
independent classes instead of 23: the full-covariance estimates (computed in the same level-B/level-C coordinates
that give 10^6.88 for 177, against an exact count of 10^7.42, and 10^13.18 for 537) are 10^−9.8 and 10^−11.4. The
537 structure is special: an odd F, few classes and exactly two spikes.


## 20. Low-energy k-trees on cyclotomic products: skew orders 820, 1044 and 1204

**Why the earlier sweeps stopped short.** The cyclotomic-product engine of §16 lists every member of each block and
groups the lists by their level-B keys, so it needs blocks of at most about 2 · 10⁷ members. The sweeps of §16 and §19
therefore skipped every structure with many spikes or wide fibres. Yet a structure's expected number of solutions
grows with its size. In the Gaussian model each condition is met with probability ~4/(2.5 · 2√v), so small v with big
blocks are the most promising. For example, the census survivor 205 with H = ⟨16⟩ (45 orbits; `census_screen.py`) has
log₁₀ E = 11 for the whole structure.

**The tool** (`search/cp_lowe.py`, `search/cp_enum.c`).

- **Additive keys.** In ℤ_m × GF(q) a member of a block is one slice pattern per K_m-orbit representative. Its key,
  the autocorrelations of level B (the compression to GF(q), classes C_0 .. C_{nb−1}) and level C (the shifts
  (b, C_c), 1 ≤ b ≤ m/2), is a sum of pair terms between representatives. So `cp_enum.c` enumerates every member of a
  block, up to about 10¹⁰ of them at about 2 · 10⁷ per second on three threads, and updates the key incrementally.
- **Low-energy filter.** With μ_k and Σ_k the mean and covariance of block k's keys (from 2 · 10⁵ uniform samples),
  a member is kept when |Σ_k^{−1/2}(key − c_k)|² ≤ τ_k. Here c_k = μ_k + Σ_k Σ⁺(t − Σ_j μ_j) is the conditional
  mean given that the four keys sum to the target t, and τ_k is the quantile that keeps about L members.
- **Join.** The kept lists go to `kt4.c`, which now also takes explicit split values (`@file`). The split values σ of
  the five merge classes are ranked by AB(σ) CD(σ) ρ(σ) / (AB + CD + |A| + |C|):
  - AB and CD are estimated from 4 · 10⁶ random pairs per side, which respects sublattices;
  - ρ is the ratio of the Gaussian law of the split of a solution (conditional on all classes, since the level-B
    classes are sums of level-C ones) to the product of the pair laws.
- **Twists.** Representatives carry a class twist when K_m acts as (a, y) ↦ (k a, g^t y). For skew blocks, t must
  be a multiple of 2^{v₂(e)}, or position 0 has no pattern that is both skew and twist-invariant.
- **Sampled stage 1.** When stage 1 is too large to list, `stage1_sample` hashes block 3's compressions by (cost,
  autocorrelations) and completes random triples from the other blocks. In 205 every class has nearly the same
  estimate (10^8.4–10^8.6 over the sampled classes), so a few sampled classes suffice.
- **Verification.** Every reported solution is rebuilt as four ±1 sequences and checked: the PAF sums, skewness, and
  invariance under the generators of H.

**Validation.**
- 933: a known solution planted in the lists is found at its split value.
- 1417 (K₁₃ = ⟨3⟩, config 2): all 8 solutions of the class are found.
- Pitfall, now fixed: ranking split values by AB · CD alone holds a vanishing share of the solutions. In 933, the
  conditional spread of a solution's split is about 49 per merge class (lattice step 4), so solutions are spread over
  about 10⁷ split values, centred where ρ puts them.

**Results.**

| v | structure | lists (A; B, C, D) | Gaussian E (class; in lists) | result |
|---|---|---|---|---|
| 205 = 5 · 41 | e = 8, F = 5, K_m = 1, H = ⟨51⟩ of order 5; 1084 stage-1 classes, class 0 | 78,400; 5.2, 4.8, 4.9 · 10⁶ of 1.7 · 10⁹, 1.1 · 10⁹, 4.3 · 10⁸ | 10^8.47; 10^3.84 | 6 solutions at the first split value, 20 in 268 s: **skew order 820**, [skew/820](../skew/820) |
| 261 = 9 · 29 | e = 4, F = 7, K_m = 1, H = ⟨190⟩ of order 7; stage 1 sampled (3.3 · 10⁶ options per block), class 0 of 5 | 5,184; all 6.7, 4.5, 4.5 · 10⁶ | 10^3.62 | 2 solutions within the first 65 split values: **skew order 1044**, [skew/1044](../skew/1044) |
| 301 = 7 · 43 | e = 6, F = 7, K_m = 1, H = ⟨127⟩ of order 7; stage 1 sampled (seed 1), class 0 of 3; chosen by `cp_hunt.py` | 64,000; 9.9, 10.1, 9.9 · 10⁶ of 9.6, 5.4, 7.2 · 10⁸ | 10^6.03; 10^2.52 | first solution at the 24th split value, about 6 minutes: **skew order 1204**, [skew/1204](../skew/1204) |

**Negative and cost results.**
- **553 = 7 · 79** (e = 6, F = 13; 36 stage-1 classes, all ≈ 10^4.5). Over the list allocations tried, the Gaussian
  estimate inside the lists is at most 10^0.16. The pair work per solution is at least about 1.4 · 10¹³, roughly
  three days at 5 · 10⁷ pairs per second. Not run.
- **589 = 19 · 31 and 817 = 19 · 43.** Untwisted K₁₉ = ⟨7⟩ and the twisted variants (order-3 twists of e = 6) give
  10^−3.5 per class.
- **301 = 7 · 43** with K₇ = ⟨2⟩ and e = 14 needs a correlation table over 15,873 slice patterns (7 GB). Not run
  (301 itself was then settled with e = 6 and K₇ trivial, above).
- **Exact chain** on the census-screen survivors with at most 33 orbits, smallest first
  (`screen_chain.py 3600 33 0 --order orbits`):
  - ruled out: 539/⟨4⟩ with rows (15, 29, 33) and (3, 11, 45), 1333/⟨66⟩ (11, 13, 71), 1269/⟨4⟩ (3, 29, 65) and
    621/⟨4⟩ (23, 27, 35);
  - 1313/⟨16⟩ (7, 51, 51) was not finished in 3600 s.

**Sweeps.** `search/cp_hunt.py` runs this estimate-then-search loop over every census composite and, in plain mode,
every composite n of FRONTIER.md (`search/cp_sweep.py` records the estimates only).

**Prime v (m = 1).** For v = q prime there is no ℤ_m to enumerate over, so `PrimeBlock` cuts the e classes of a
block into G groups (default 4) that play the part of positions. A member is one sign pattern per group; its key,
PAF(d_c) for c < nb, is a quadratic form in the class signs whose coefficients are the cyclotomic numbers of order e,
so it splits into pair terms between groups. The class sum s is fixed through tags in `cp_enum.c`: the search
prunes on the reachable tag range and buckets the last group's patterns by tag. Stage 1 is exact
(`stage1_prime`: row sums y = f + F s with Σ y² = 4v).
- Stage 1 is empty, i.e. no four admissible row sums have Σ y² = 4v, for every structure that `cp_hunt.py` admits
  (blocks of 10⁵ to 3 · 10¹⁰ members) at these census primes, with their indices e: 101 (20), 191 (38), 409 (24),
  419 (22, 38), 443 (26, 34), 457 (24), 541 (36), 571 (30, 38), 599 (26), 647 (34, 38), 673 (32), 751 (30),
  829 (36), 929 (32), 1051 (42), 1103 (38), 1117 (36), 1303 (42), 1321 (40) and 1471 (42). 883 (e = 42) ran out of
  memory building the pair tables of 4 groups.
- 761 (e = 40, F = 19): the best class has 10^9.5 solutions, but only 10^−0.2 inside lists of 10⁷ members, at a
  predicted 10^14.2 pairs per solution. Not run. 1483 (e = 38): 10^6.35 and 10^−1.7. Not run.
- 1123 (e = 34, F = 33, rows 1, 1, 1, 67, one class): 10^6.63 for the class and 10^2.24 inside lists of 6 · 10⁷, 10⁷
  and 10⁷ members (all 131,072 members of the skew block). Predicted 10^11.8 pairs per solution. A run with lists of 3
  · 10⁷, 10⁷ and 10⁷ found no solution in 6000 split values (§22); the run with 6 · 10⁷ ran out of memory after 72
  split values.

**The Yamada–Kiyasu routes of FRONTIER.md have no usable multiplier structure.** H(7516), H(8860) and H(9532) need
skew orders 940 = 4 · 235, 1108 = 4 · 277 and 1192 = 2 · 596 = 8 · 149. A skew Goethals–Seidel family over ℤ_v with
multipliers H needs every h ∈ H to move each nonzero t off −t, which leaves H of odd order:
- v = 149 (148 = 4 · 37): H = 1 or the fourth powers (F = 37). With F = 37 a spike has |y| ≥ 73 and 73² > 596, so
  every class sum is 0 and Σ y² = 4 ≠ 596.
- v = 235 = 5 · 47 (ℤ₂₃₅* ≅ ℤ₄ × ℤ₄₆): H = 1 or the squares of GF(47) (h ≡ 1 mod 5), with F = 23 and 45² > 940.
  Even H ∋ −1 with a sign character on the skew block fails: the character must be trivial on the stabilisers of
  (a, 0) and (0, y), which forces χ(−1) = 1, unless |H| ≤ 2.
- v = 277 (276 = 4 · 3 · 23): the orders 69 and 23 (F = 69, 23) fail in the same way; order 3 leaves 92 classes, so
  2⁹² members per block, far beyond exhaustive enumeration.
So these routes need families with (almost) no symmetry. ℤ₂₃₅ is the only group of order 235 and ℤ₁₄₉ the only one
of order 149.

## 21. Legendre pairs for the open orders (negative), prime-power fibres, faster enumeration and joins

**Legendre pairs over the same structures do not reach the open orders.** A Legendre pair of length l = 2n − 1
gives a Hadamard matrix of order 4n, and a two-list join keeps far more solutions in low-energy lists than a
four-list one: for equal covariances the key of A in a solution has covariance Σ_A/2, and the whitened offsets of A
and B coincide up to sign, so a list of 10⁻³ of a block with rank 40 still holds about a third of the solutions.
`search/lp_lowe.py` implements it (stage 1 on ℤ_m with spike costs summing to 2(m − 1)e, PAF sums −2q; level B target
−2m, level C −2; exact hash join) and `search/lp_hunt.py` runs it over the FRONTIER orders. It finds Legendre pairs
of length 143 (ℤ₁₁ × GF(13), e = 4, K₁₁ = ⟨3⟩: the twin-prime class). But the count is against it:
- An H-invariant pair has two blocks of O orbits each, so about 2 · 0.3 · O decimal digits of freedom, against about
  O/2 conditions that each cost log₁₀(√(4πl)/4). That gives log₁₀ E ≈ O (0.6 − 0.5 log₁₀(√(4πl)/4)), which is
  negative for l ≳ 320. Goethals–Seidel structures have twice the freedom per condition (four blocks of length n
  against two of length 2n − 1), which is why their FRONTIER estimates are 10²⁰ and more.
- Measured: l = 1501 (n = 751), ℤ₇₉ × GF(19), e = 6, K₇₉ of order 13, the only stage-1 class: blocks of
  7.2 · 10⁸ and 2.2 · 10⁸ members, Gaussian estimate 10^−24.8 (rank 20).
So the route is closed for these lengths except through algebraic families.

**Prime-power fibres** (`search/cp_field.py`, `cp_lowe.py --field p,k`). The cyclotomic product ℤ_m × GF(p^k) lives
on the group ℤ_m × ℤ_p^k: field elements are coordinate vectors in the basis 1, w, …, w^{k−1} (w a root of a
primitive polynomial), the correlation table comes from a k-dimensional FFT, and the generators of H are integer
matrices on the coordinates, as `verify.build_goethals_seidel_group` checks them. Hits are checked by the group PAF,
the group negation and the automorphisms, and `skew_certify_group.py` writes the certificate.
- Validation: 205 = 5 · 41 as ℤ₅ × ℤ₄₁ (k = 1) gives 2 verified solutions in the class that gave skew 820.
- 75 census values have a prime-power fibre with F ≥ 3 odd and at most 2 · 10⁵ patterns per position, e.g.
  GF(25) with e = 8 (F = 3) for 325, 475, 575, 725, …, GF(49) with e = 16 for 245, 441, 539, 637, …, GF(121) for
  605, 1089, 1573, GF(169) for 845, 1183, 1859.
- 325 = 13 · 25 with K₁₃ = ⟨3⟩: class 0 of the sampled stage 1 has blocks of 0.6–8.6 · 10⁸ members and a Gaussian
  estimate of 10^−0.44 (rank 19). Not promising; the hunt ranks the other classes and structures.

**Speed.**
- `cp_enum.c` scans the last position in a transposed structure-of-arrays layout, with the patterns of each class sum
  contiguous. The kept members are bit for bit those of the generic recursion (checked on blocks of 261, 301, 409 and
  571); prime-field blocks run 2.3 times faster, composite ones about 1.2 times.
- `kt4.c` sorts all four lists by their merge tuples and finds the matching groups by one merge pass (by binary
  search when one list has few groups). The pairs of two groups then read contiguous memory instead of a binary
  search and a random hash load per member. On the 261 lists the same 12 solutions come in 79 s instead of 395 s.

**The 1123 run.** With lists of 5.8 · 10⁷, 10⁷ and 10⁷ (10^2.23 solutions expected inside), the old `kt4` processed 72
of 6000 split values in 830 s without a solution before it ran out of memory. The run with lists of 3 · 10⁷, 10⁷, 10⁷
(10^1.98) and the new `kt4` is in §22: no solution.

**Bordered structures through the same machinery** (`search/bcp.py`: `bgs.Structure` seen as a cyclotomic product,
level-B/level-C targets −(m/q) Σ_α Γ(α) and −Γ(b mod q)). On the 2148 structure (ℤ₆ × GF(89), e = 8, Q = ℤ₃) the class
estimate is 10^12.6 and the estimate inside low-energy lists of 10⁷ members per block is 10^3.98, so the direct
four-list search would find it without the U/W split. (In-list estimates need about 10⁶ samples per block when the
lists keep less than 10⁻⁴ of a block; with 6 · 10⁴ samples the 2148 figure came out as 10^−56.)
- Total Gaussian estimates (all row patterns and stage-1 classes; random H-invariant sequences, PAF at the orbits of
  ⟨H, −1⟩ on nonzero shifts) for 573: with Q = ℤ₃ (border 12) every structure with at most 64 orbits per block is
  below 10^−0.6; Spence-type (t = 1) ℤ₅₂ × GF(11), e = 2, K₅₂ of order 3: 10^10.9 over 60 orbits, spread over very
  many stage-1 classes.
- Larger totals need character sets whose border problem is not known to be solvable: 573 with ℤ₈ × GF(71), e = 10, Q
  = ℤ₈, T = characters of order 1 and 8 (t = 5): 10^16.6 (55 orbits); 643 with ℤ₁₂ × GF(53), e = 4, Q = ℤ₁₂, T of
  orders 1, 3, 12 (t = 7): 10^16.3 (60 orbits); 599 with ℤ₁₆ × GF(37), Q = ℤ₈, orders 1, 4, 8: 10^7.6. For Q = ℤ₈ and
  T = {1} ∪ {order 8}, a ±1 row of R̂ is on each block either constant or antiperiodic (f(α + 4) = −f(α)), since its
  compression to ℤ₄ must be constant. Whether such borders exist is open here.

**Reach of the low-energy four-list search.** In the tail regime a list of L members drawn from a block of N members
keeps a share of about g·L/N of the solutions (g ≈ (4/3)^{r/2} for rank r and equal covariances; the factor is
(1/2)^{−r/2} for the two-list case). So the expected number inside the lists is E_class · Π_k g L_k/N_k, while the
join costs about |A| L_B + L_C L_D pairs. Enlarging a structure (bigger N, larger r) at a fixed list size therefore
lowers the yield; the useful structures are those whose blocks can be kept nearly whole, N ≲ 10¹⁰, with E_class ≳ 10³.
The census map (`census_reach.py`, program not included; necessary norm and row-sum condition only, blocks of 10⁵–10¹¹
per class) finds such structures for 44 of the undecided composite v up to 1089 — 119, 225, 319, 325, 369, 429, 473,
475, 481, 539, 553, 575, 589, 609, 623, 627, 637, 657, 689, 693, 707, 725, 737, 777, 781, 783, 793, 817, 845, 851,
891, 897, 899, 925, 927, 933, 943, 949, 957, 989, 1073, 1075, 1081, 1089 — and for none of the others. Many of these
fail the autocorrelation part of stage 1 or have class estimates below 1 (e.g. 481: 10^1.5; 589: at most 10^−2.6), so
the method has a natural frontier in the census; most undecided v admit no multiplier structure of searchable size at
all.
- Non-cyclic K_m (two generators, no twist, m ≤ 120) add a single census structure that passes the norm and row-sum
  condition with blocks of 10⁵–10^11.5 per class: 2223 = 117 · 19, e = 2, |K| = 9, 25 position orbits.


## 22. Bordered arrays over ℤ₄, pair-spin structures, and calibrations

**The skew 4492 run is negative.**
`cp_lowe.py 1123 --q 1123 --e 34 --config 0 --L 3e7,1e7,1e7 --nsig 6000 --threads 3 --K 5` (seed 1): lists of 131,072
(the whole skew block), 2.9 · 10⁷, 1.0 · 10⁷ and 1.0 · 10⁷ members, 6000 split values holding a predicted share of
0.0746 of the 10^1.98 solutions inside the lists, 7.4 · 10⁸ AB and 2.3 · 10¹¹ CD pairs in 4932 s: no solution. With
about 7 expected, the prime-field estimate inside low-energy lists is optimistic by a factor of at least 2.4 (95 %).
For comparison, the four blocks of the known plain solution of 4492 (same H and rows 1, 1, 1, 67) lie at energy
quantiles 2 · 10⁻⁴ to 1 · 10⁻³ of their plain blocks (ranks 5 · 10⁵ to 2 · 10⁶ of 2.3 · 10⁹), well inside lists of
10⁷. That solution was found by a search that favours low energy, so this does not test the estimate.

**Faster joins and resumable runs.** On the lists of the 319 run (AB tables of 3–5 · 10⁶ pairs per split value), the
radix-partitioned join of `kt4.c` takes 27 s per split value against 39 s (same four values, same solutions).
`cp_lowe.py --skip-sigmas FILE` drops split values already searched (kt4 log lines) when a run is restarted; with
`--S` fixed the merge classes stay the same. The 319 run was restarted this way (131 of its first 267 values recur in
the new ranking).

**Spin local search on composite v.** `spin_ls.c` swaps signs only between H-orbits of equal size, and each run
started from one random B: for composite v the number of + signs per orbit size then stays that of the start for the
whole run. For 817 (orbit sizes 1, 3, 7, 21) four such splits are possible, with 67, 27, 5 and 1.4 % of the B pool, so
a run explored one of them. `spin_skew.py csearch` now passes 4096 uniform starts and each restart draws one. (721 has
a single split: sizes 1 and 17.) Runs: 817, |H| = 21, μ = 9 (8 K-classes, estimate 10^3.2), and 721, |H| = 17, μ = 58
(10^2.3), 3000 steps per restart, about 24 restarts per second.

**Pair-spin structures** (`search/pairspin.py`). With μ ∉ H, μ² ∈ H and K = ⟨H, μ⟩, a pair of blocks (X, X ∘ μ)
contributes T_X(s) = P_X(s) + P_X(μs), which is K-invariant.
- A skew block invariant under K needs −1 outside the image of K in ℤ_p^* for every prime p | v (otherwise an orbit
  of elements of order p is closed under negation), so every such image has odd order, and so has K (the kernel of
  ℤ_v^* → Π ℤ_p^* has odd order). Hence no layout with a K-invariant skew block and [K : H] = 2 exists; the spin
  structure of §10 (μ of order 3) is consistent with this.
- The remaining layouts are (A, A ∘ μ, C, D) with C, D K-invariant (2 + c² + d² = 4v) and (A, A ∘ μ, C, C ∘ μ)
  (2v − 1 = c²). Both are exact two-list joins.
- Validation: v = 63, |H| = 3, μ = 2, rows (1, 1, 9, 13): skew-Hadamard matrices of order 252 (the Gaussian estimate
  was 10^−2.5). Among the composite v < 100 with 2v − 1 a sum of two squares, solutions exist for 15, 21 and 63; for
  25, 27, 33, 45, 55, 57, 75, 77, 87, 91 and 93 every structure within the caps (2²² skew members, 3 · 10⁸ pairs) is
  empty (`pairspin.py sweep ... --maxA 22 --maxpairs 3e8 --maxlist 2e7`). A census sweep (largest H first) found
  nothing for 201, 253, 261 and 265 and was paused: the family thins out quickly.

**Border problems over ℤ_q** (`search/border_probe.py`). The rows of R̂ and the columns of Ĉ are ±1 vectors of
V_T (four blocks, Fourier support in T) and form orthogonal bases of it, and Mc Mcᵀ = 4t on V_T. For a fixed basis
R̂, the border condition becomes Ĉ = −(1/4t) Mc R̂ᵀ X: the columns of X are the sign vectors ε with
|Mc R̂ᵀ ε| = 4t entrywise (a set S_R found by enumerating 2^{4t−1} vectors), and a border exists for R̂ iff S_R holds
4t mutually orthogonal vectors. Only the components of the compressions at the characters in T enter.
- Calibration: Q = ℤ₃, T = all characters (solved by `border3.c`): random bases R̂ give a border at once.
- **Q = ℤ₄ with T = {1, χ, χ³} (the trivial and the two characters of order 4) and with T = {χ, χ², χ³}:** both
  border problems are solvable. `border_probe.py validate 4 3 2 1,4` and `... 2,4` search all sequences on
  G = ℤ₁₂ = ℤ₄ × ℤ₃ (targets −Γ(g mod 4): −12, −4, +4, −4 and −12, +4, +4, +4), find borders for their compressions
  and assemble Hadamard matrices of order 60 (H Hᵀ = 60 I for all six tested). The second set has the targets of
  Alpöge et al.'s ℤ₂ × ℤ₂ family. So every n ≡ 3 (mod 4) has bordered structures with border 12 over the cyclic group
  ℤ_{n−3} (and over ℤ₄ × K), where the ℤ₂ × ℤ₂ family needs a non-cyclic group. The census of these structures for the
  FRONTIER orders 599, 643 and 647 is in §23.
- Q = ℤ₈, T = {1} ∪ {characters of order 8} (t = 5, the 573 lead of §21): the compressions have row sums (±4, ±2, 0,
  0) and an antiperiodic part with five ±1 cells and Σ NPAF = 0 off 0 (48,640 patterns). Random bases give S_R = ∅
  (each R̂ is compatible with about 12.5 % of all ±1 vectors of V_T per row). Unresolved.

**The large bordered totals of §21 were artifacts.** The totals of `btotal_scan.py` (program not included) average
over random sequences and ignore the exact compression budget Σ (F s² + 2 f s) = 4(m − t)e of `bgs.stage1`:
- 573, ℤ₈ × GF(71), e = 10, T = {1} ∪ {order 8}: stage 1 is empty for K_m = ⟨3⟩ and ⟨7⟩ (the 10^16.6 structures); with
  K_m trivial (88 orbits per block) it is too large for the exact join.
- 643, ℤ₁₂ × GF(53), e = 4, T of orders {1, 3, 12} or {1, 4, 12}: stage 1 empty (10^16.3).
- 599, ℤ₄ × GF(149), e = 4 (T of orders {1, 4}): the budget is 16 while a spike costs at least 4F − 4 = 144, so stage
  1 is empty (the total was 10^3.8).

**Prime-power fibres.** `cp_hunt.py --fibre-only --estimate-only` over the 76 census values with a fibre GF(p^k):
GF(25) and GF(49) structures (F = 3) come out between 10^−67 and 10^−0.04 per class; the best lead is 845 = 5 · 13²,
ℤ₅ × GF(169), e = 8, F = 21, K₅ trivial, class 1 of the sampled stage 1, rows (1, 43, −3, −39): 10^4.54 for the class
and 10^1.48 inside lists of 5 · 10⁷, 2.5 · 10⁷ and 2.5 · 10⁷ members (blocks of 1.3–1.7 · 10⁹). Searched after 319:
skew order 3380 (§23).

## 23. Skew orders 1276 and 3380, a symmetric order 532, and calibrations

**Skew-Hadamard matrix of order 1276** ([skew/1276](../skew/1276)). The 319 run of §21–22 (ℤ₁₁ × GF(29), e = 4, F = 7,
H = ⟨45⟩ of order 7, stage-1 class 0 of the sampled three, rows 1, 25, 19, 17; lists of 31,104 (skew block),
60,466,176 (the whole second block), 9,966,434 and 10,012,823 members; estimate 10^1.64 inside the lists, share 0.09
in 6000 split values) found a solution at the 321st split value of the resumed run (the first run had searched 267
values; 131 recurred), about 450 distinct split values in all. With about 4 solutions expected over 6000 values and
the densest values first, this is in line with the estimate for composite v (unlike the prime-field 1123 run of §22).
Census: undecided; Cati–Pasechnik Table 5: 319(3); no classical construction gives skew 1276 (its README).

**Skew-Hadamard matrix of order 3380** ([skew/3380](../skew/3380)), the lead of §22: ℤ₅ × GF(169), e = 8, F = 21,
H = ⟨w⁸⟩ of order 21, stage-1 class 1, rows 1, 43, −3, −39. Lists of 78,400 (skew block), 50,002,409, 25,105,989 and
25,135,175 members (estimate 10^1.48 inside them); four solutions at the 26th split value, 732 s in all on three
cores. The group is non-cyclic; the certificate uses the `goethals-seidel-group` builder with the automorphism
w ↦ w⁸ written as a matrix on ℤ₅ × ℤ₁₃². Census: undecided; Table 5: 845(4), so 6760 is new as well.

**Layer separation for F > 4m.** In a cyclotomic product ℤ_m × GF(q) of index e (class size F), the autocorrelation of
a block at the shift (b, 0) is P(b, 0) = PAF_f(b) + F · Σ_i PAF_{ε_i}(b), where f(a) is the value at (a, 0) and
ε_i(a) the value on (a, C_i). Since |Σ_k PAF_{f_k}(b)| ≤ 4m, the condition Σ_k P_k(b, 0) = 0 splits when F > 4m:
- the f-layer (f₁, …, f₄) is a Goethals–Seidel quadruple over ℤ_m, with f₁ of skew type in the skew case (a
  skew-Hadamard matrix of order 4m of Goethals–Seidel type);
- the 4e class sequences ε_{k,i} of length m are periodic complementary (Σ PAF = 4em δ).
Both hold in the 3380 solution (m = 5, F = 21): its f-layer has rows 1, 1, −3, 3, and the 32 class sequences of
length 5 have PAF sum 160 δ. Stage 1 already implies the first statement (y = f + F s fixes f), so this is a
description, not a new filter; it says which ℤ_m-level objects a large-F solution contains.

**A symmetric Hadamard matrix of order 532** ([symmetric/532](../symmetric/532)). The exhaustive propus sweep of §17,
with the orbit cap raised from 24 to 30 (`propus.py 133 --sweep --mino 1 --maxo 30 --maxpool 5e7`), found 72 solutions
in the structure H = ⟨26⟩ of order 6 (A on the 14 orbits of ⟨H, −1⟩; B and D on the 26 orbits of H), only A
symmetric, rows (11, 5, 19). The census lists 133 among the undecided symmetric values; ordinary and skew matrices of
order 532 are known. `verify.py` gained the construction `propus` (it checks H = Hᵀ).
- **Gaussian estimate (§29).** Without the lattice factor of the class values (about C log₁₀ 4), `propus_rank.py`
  gives log₁₀ E = −7.33 for exactly this structure and row pattern (pools 10^3.0, 10^6.5, 10^6.4; 12 independent
  classes), against 72 solutions; with the factor the estimate is 10^0.5, and 72 is a factor of about 10^1.4 above it.
  The other row patterns of the structure have none.
- Sweeps (`--maxo 30`) over the symmetric census values 89 < v ≤ 400 and 400 < v ≤ 1000
  (`$HADAMARD_WORK/logs/sym_queueA.log`, `sym_queueB.log`) were not completed; the sweeps of small structures for
  v > 275 were stopped (§30). v = 89 has no structure with a non-empty pool pattern within the cap.

**Symmetric Hadamard matrices: propus estimates without the lattice factor.** `search/propus_rank.py` (Gaussian
estimate of the propus structures of `propus.py`) gives, without the lattice factor (§29), for the structures in which
Balonin and Đoković found propus difference families: v = 67, H of order 3 (A symmetric, D free) 10^−2.9; v = 109, H
of order 3, 10^−6.6; v = 151, H of order 5, 10^−6.5. Solutions exist in all three (six, three and one families in
their Table 1). These numbers are too low by roughly C log₁₀ 4 for C classes.
- v = 65 (the smallest open symmetric order, 260): every cyclic propus structure scores at most 10^−4.7 (H of order 3,
  25 orbits); the exhaustive sweep `propus.py 65 --sweep --maxo 30 --maxpool 1e8` (all structures up to 30 orbits, 207
  row patterns, including the 25- and 26-orbit structures beyond the earlier cap of 24) finds no solution.
- Larger structures (v = 163 with H of order 3: 28/55/55 orbits, 10^−12.1; v = 193: 10^−15.5; both without the lattice
  factor of §29; for 163 with it see §40) need a low-energy three-list search, which has not been written.

**Parallel tempering for the spin structures** (`search/spin_pt.c`, same input as `spin_ls.c`; replicas at geometric
temperatures, swap or pair-flip moves in O(C), replica exchange after every sweep). On 631 it needs 5–23 s per
solution (tabu: about 1 s); on 919 none in 10 CPU-minutes with 24 replicas between T = 20 and 20,000. It does not
beat the tabu search with restarts.

**ℤ₄ bordered families (§22) for 599, 643 and 647** (`search/z4triage.py`, both character sets T, every multiplier
group with 14–44 orbits): 30, 437 and 192 structures, stage 1 empty in every one.

**Census hunt** (`cp_hunt.py skew`, composite census values up to 627 in the order of the list, `--mine 0.5`): beyond
301, 319 and 553 (already run) no structure reaches the threshold; 481 (ℤ₁₃ × GF(37), e = 4) has 10^1.46 for its best
class but only 10^−3.95 inside lists of 10⁷; most other values have empty stage-1 samples or blocks above the
enumeration cap (627: 1.15 · 10¹¹ members). The sweep was stopped at 627.

## 24. Symmetric and non-cyclic sweeps

**Tools.**
- `search/propus_group.py`: the propus search over non-cyclic groups G = ℤ_p² × … × ℤ_m (group-developed blocks, R[g,
  h] = [g + h = 0]; the arrays of §17 hold over any abelian group). Structures: subgroups of GL(2, p) of any order
  (Singer, diagonal, torus, unipotent, Borel) on each split prime times any subgroup of ℤ_m*. Check: symmetric
  matrices of order 100 over ℤ₅² (Singer subgroup of order 3).
- `search/propus.py --plain`: Goethals–Seidel quadruples (A, B, B, D) on H-orbits for ordinary Hadamard matrices (hits
  in `gs_propus_hits.jsonl`); `--maxwork` skips row patterns whose pair work exceeds a cap, the two smallest lists now
  form the pairs, and `--only-freeD` sweeps only the variant with A symmetric and B, D free (the form of every
  computer-found symmetric matrix of the census, §17).
- `search/skew_group_sweep.py`: every skew structure over a non-cyclic group (`skew_group.py`, odd-order subgroups)
  with 8–30 orbits and every row pattern, exactly (find4b) when the pair work is at most 3 · 10¹⁰ and otherwise by the
  k-tree of `skew_ktree.run` with 2000 split values (it finds the 1100 solutions at the 810th).
- Speed: `Structure.class_values` is one matrix product per block of rows (7 times faster, identical values) and
  propus pools are enumerated by meet in the middle (`pool_exact`; 400 times faster at 26 orbits). The 133 pattern
  that holds the 532 solutions takes 26 s instead of 175 s.

**Symmetric Paley products** (`search/propus_paley.py`). For v = m q with q ≡ 1 (mod 4) a prime power and
H = {h ≡ 1 (mod m), h a square in GF(q)}, the quadratic character χ is symmetric, and every block is
X = u ⊗ δ + w ⊗ J + z ⊗ χ as in §15. The propus condition with weights c = (1, 2, 1) becomes
(I) Σ c_k [u_k u_k* + q z_k z_k*] = 4mq δ, (II) Σ c_k [u_k w_k* + w_k u_k* + q w_k w_k*] = Σ c_k z_k z_k*,
(III) Σ c_k [u_k z_k* + z_k u_k*] = 0, and the identity coefficient of (II) gives (q + 1) N_w − 4 N₋ = 4m with
c-weighted counts, so N_w is even. N_w = 2 needs q ∈ {2m − 1, 2m + 1, 2m + 3}; with one w-position in B the size bound
forces u_B to be skew away from 0 and R = 0, so the z-part is itself a propus-type system with z_B(0) = 0.
- Equations checked on the symmetric matrices of order 60 (v = 15, which use N_w = 4).
- N_w = 2, exact (both layouts): **no solution** for 65, 119, 153, 377, 435 and 629 (all stage-1 triples tried; 703
  and 779 were stopped). For m = 13 the u-systems (1,664) and z-systems (74,776) exist separately but no pair has
  zero cross terms (III); the cross values are ≡ 0 (mod 4) with no further congruence, so this is a count, not an
  obstruction. Unlike the skew products of §15, the symmetric ones look thin.

**Ordinary Hadamard: (A, B, B, D) for the open orders** (`propus.py n --sweep --plain --maxo 30 --maxwork 5e9`). For
the first 47 open orders of FRONTIER.md (n = 573 to 1553; the run was stopped there, and n = 1169 was interrupted),
every multiplier structure with at most 30 orbits has an empty pool for some block in every row pattern (the orbit
sizes cannot reach the row sums), so nothing was joined. Small structures do not exist for these n, as for the
four-block searches of §5.

## 25. Compressed propus searches, D-optimal designs and Legendre pairs

**`search/propus_chain.py`: compressions make large symmetric structures cheap.** For d | v the compressions of the
three blocks to ℤ_d satisfy PAF_a′ + 2 PAF_b′ + PAF_d′ = 0 off 0 and are constant on the H-orbits of ℤ_d. Stage 1
lists every compression each block attains (a dynamic program over the orbits, vectorised over mixed-radix codes) and
joins the three lists; stage 2 lists the members of each compatible triple by meet in the middle (the two halves of
the orbits are built once per structure and keyed by row sum and compression) and joins them with join4x. Nothing is
sampled. The modulus is the largest prime factor whose H-orbits on ℤ_d number at most 10.
- Check: 133/⟨26⟩, rows (11, 5, 19): with d = 19, 15,096 compatible triples and pair work 1.1 · 10⁸ against
  2.3 · 10⁹ for the direct join, and the 532 solution is found; with d = 7, two triples hold all 72 solutions.
- 119 = 7 · 17, structures with 31–40 orbits: pair work 6 · 10⁷ – 5.5 · 10⁸ per row pattern against 3 · 10¹² direct
  (four to five orders of magnitude). 65 = 5 · 13 with 31–44 orbits: 50 compatible triples, 3.1 · 10⁸ against
  1.6 · 10¹⁴.
- `--plain` runs the same search for Goethals–Seidel quadruples (A, B, B, D) with A on the H-orbits (ordinary
  Hadamard matrices); checked on n = 21 and 33.
- The sweep over all composite symmetric census values v ≤ 1000 (structures with at most 36 orbits, freeD form; a
  37–44-orbit pass is not recorded) found order 876 (§28) and covered the values up to 435 before it was stopped
  (§40). No result is recorded for the sweep of the composite open orders of FRONTIER.md with `--plain` and 20–44
  orbits.

**D-optimal designs and Legendre pairs** (`search/two_circ.py`, `search/two_chain.py`). Two ±1 sequences with
PAF_A + PAF_B = 2 off 0 give a D-optimal design of order 2v ([[A, B], [−Bᵀ, Aᵀ]], over any abelian group); with −2
they form a Legendre pair. Open according to Đoković (arXiv:2210.03931): D-optimal designs of order 2v for
v = 99, 115, 123, 135, 141, 147, 153, 159, 163, 167, 169, 175, 177, 185, 187, 195, 199 (his Conjecture 1 allows any
abelian group); Legendre pairs of lengths 115, 145, 159, 161, 169, 175, 177, 185, 187, 195 (arXiv:2111.02105v3).
`two_circ.py` joins complete pools exactly (sorted 64-bit keys) and sends larger pools through the power-spectrum test
|X̂(χ)|² ≤ 2v − T of `psd_pool.c` (meet in the middle in C); `two_chain.py` is the compressed version.
- Checks: D-optimal designs for v = 13, 25, 41 and Legendre pairs for 13 and 21 are found and verified.
- **No solution** in the cyclic structures with at most 48 orbits whose pools (or PSD survivors) fit in 2 · 10⁷:
  D-optimal for all 17 values above; Legendre pairs for 115, 145, 159, 161, 169, 175, 177, 185 and 187 (195 was
  stopped). Non-cyclic groups (ℤ₃² × ℤ₁₁ for 99): no D-optimal design with at most 48 orbits. Structures whose pools
  exceed the caps were skipped and are listed in the logs. The power-spectrum test alone empties several pools
  (for instance Z145⟨28, 86⟩ and the non-cyclic 99 structures singer2 and unip with ℤ₁₁* of order 10: no sequence
  passes), which rules those structures out at once.

**Williamson matrices of orders 65 and 67** (`search/williamson.py`; the only odd orders below 71 for which
Kharaghani, Mohammadian, Tayfeh-Rezaie and Zaitsev, arXiv:2605.08661, list Williamson matrices as not yet found):
four symmetric ±1 sequences on the orbits of ⟨H, −1⟩, all row patterns, exact four-list joins. Checked on n = 13, 25
and 37. No result is recorded.

## 26. Local search without symmetry, and the skew hunt beyond 627

**Why unstructured searches are tempting.** A crude Gaussian count for cyclic families without multipliers (each
equation P(s) = 0 met with probability about 4/(σ√(2π)), σ² the variance of P(s); about 10 bits for the row-sum
pattern) gives log₂ E ≈ v (2.5 − ½ log₂(0.63 √(7v))) − 10 for propus triples (A symmetric, B and D free), that is
about 2³³ triples for v = 65 and 2³⁵ for v = 89, positive up to v ≈ 370; for skew Goethals–Seidel quadruples
log₂ E ≈ v (3.5 − ½ log₂(0.63 √(5v))) − 13, about 2¹³⁴ for v = 89. Under the same count a multiplier group of order 2
leaves about 2⁻⁴ propus triples for v = 65: symmetry discards almost all of the expected solutions, and for the
smallest open symmetric order (260) and skew order (356) only the unstructured space is large. The question is
whether a solution can be *found* there.

**It cannot, by local search** (tabu search, every move evaluated exactly, restart after 20,000 steps without a new
best; `search/propus_ls.c`, the general `search/gs_ls.c` for the types skew / symmetric / free with weights, and
`search/orbit_ls.c` with `orbit_ls.py` for multiplier-orbit structures). Steps per solution:

| v | 25 | 29 | 33 | 35 | 37 | 41 | 45 |
|---|---:|---:|---:|---:|---:|---:|---:|
| propus (A symmetric, B, D free) | 3.0 · 10⁴ | 3.8 · 10⁵ | 9.4 · 10⁵ | 3 · 10⁶ | 5.1 · 10⁶ | 1.0 · 10⁸ | none in 9 CPU-min |
| skew Goethals–Seidel | | 5 · 10³ | | | 1.0 · 10⁶ | | none in 7 CPU-min |

Every solution is checked (the propus ones as matrices: v = 29 gives a symmetric Hadamard matrix of order 116 with
H = Hᵀ and H Hᵀ = 116 I). The cost grows by about a factor 1.7 per unit of v for propus and 2 for skew quadruples,
which puts v = 65 near 10¹³ steps (about three CPU-years at the rate reached for v = 41) and v = 89 far beyond.
Two refinements do not help:
- **Fixed compressions** (`search/comp_ls.c`, `comp_margins.py`): for v = p q every block is kept on the arrays with
  given compressions to ℤ_p and ℤ_q (margin-preserving 2 × 2 switches, mirrored for symmetric and skew blocks), so the
  compressed equations hold throughout. With the margins of a known v = 35 solution it finds nothing in 40 CPU-s
  (6.2 · 10⁵ steps, best energy 320) where the unconstrained search finds two solutions.
- **Multiplier orbits** (`orbit_ls`): on the 1276 structure (ℤ₃₁₉, ⟨45⟩ of order 7, 55 orbits, 27 classes, where the
  k-tree of §23 found a solution) no solution in 3.1 · 10⁶ steps (300 s; best energy 2176), and on the 532 structure
  (ℤ₁₃₃, ⟨26⟩, 26 orbits) none in 3.8 · 10⁶ steps (best 288), against 26 s for the exact join.
So the exact joins and low-energy k-trees remain the tools; local search is only competitive below about 40 free
positions per block.

**Compression DP.** `Comp.dp` (propus_chain.py, used by two_chain.py) drops partial row sums that can no longer reach
the target and gives up above 3 · 10⁷ states (numpy codes) or 3.75 · 10⁶ states (dictionary); the structure is then
reported as skipped instead of exhausting memory. The 133 check gives the same 15,096 compatible triples and the
solution. The v = 99 and 123 D-optimal runs that died in the DP were not completed.

**Skew census hunt beyond 627** (`cp_hunt.py --trivial-k`, then all K_m). Of the 410 remaining composite census values
in 629..2499, 213 have cyclotomic-product structures within the member caps (477 structures, 89 with K_m trivial). The
trivial-K structures have large class sizes (F = 11 … 45) and their sampled stage 1 is empty except for 1565
(ℤ₅ × GF(313), e = 8, F = 39: 10^1.93 for the class but 10^−2.45 inside lists of 10⁷, blocks of 10^9.2), 1313
(10^−1.86) and 927 (a block of 3.8 · 10¹¹ members, above the enumeration cap). With F > me a spike alone exceeds the
budget Σ y² = 4v, so large-F structures are infeasible as soon as F ≳ me. No result is recorded for the structures
with K_m nontrivial.

## 27. Bordered propus: symmetric Hadamard matrices of order 4(m + 1) from groups of order m

**Construction.** Let G be an abelian group of order m, A, B, D ±1 functions on G with A(−g) = A(g), and
PAF_A + 2 PAF_B + PAF_D = −4 at every g ≠ 0. Let P be the propus array of §17 on the group-developed blocks
(P = Pᵀ and P² = I₄ ⊗ (AAᵀ + 2BBᵀ + DDᵀ) = I₄ ⊗ (4(m + 1) I − 4J)), S the same 4 × 4 array on the row sums
(a, b, d), K a Hadamard matrix of order 4 and L = −K S Kᵀ / 4. Then
    H = [[L, K ⊗ eᵀ], [Kᵀ ⊗ e, P]]
is symmetric, and H² = 4(m + 1) I exactly when Kᵀ K = 4I (the J-blocks of P² are cancelled), L K + K S = 0 (this
defines L) and L² = 4I, which holds iff a² + 2b² + d² = 4. For even m the row sums are even, so (a, b, d) is
(±2, 0, 0) or (0, 0, ±2); a K with L ∈ {±1}^{4×4} exists in every case we met. A census value v is thus reached
through any abelian group of order v − 1, in place of ℤ_v.
- Checks: every triple of ℤ_m, m = 2, 4, …, 10, satisfying the conditions (16, 96, 576, 4096 and 1600 triples)
  gives a symmetric Hadamard matrix of order 4(m + 1); the searches below find and verify orders 36 … 84 on ℤ_m
  and 52 on ℤ₂ × ℤ₆.
- Real characters: at a character χ of order 2 the transforms are even integers with Â² + 2B̂² + D̂² = 4v. For
  v = 65 this forces B̂(χ) = 0 (65 − 2β² is a sum of two squares only for β = 0), so B sums to 0 on every coset of
  2G; on ℤ₂⁶ every character is real, B would be constant, and the group is excluded.
- Literature: the propus papers we checked (Seberry–Balonin 2017, Balonin–Đoković arXiv:1803.08787, Đoković
  arXiv:2101.05429, whose propus conjecture is verified for v ≤ 63) use ℤ_v only; we have not seen this bordered
  form, but have not made an exhaustive search for it.
- Character-based blocks χ(αᵗ + c), χ(αᵗ + α⁻ᵗ + c), χ(α²ᵗ + c) and χ((αᵗ + c)(α⁻ᵗ + c)) on GF(q)* ≅ ℤ_{q−1}
  (q prime, 13 ≤ q ≤ 113, both values at the zero of the argument, both signs) give no triple: no cheap algebraic
  family of this shape.

**Tools.** `propus.py --border` (direct joins over ℤ_m), `propus_chain.py --border` (compressions to any divisor
of m), `propus_bgroup.py` (non-cyclic groups: all automorphisms by enumeration of the matrices with
M_ij ∈ (n_i / gcd(n_i, n_j))ℤ, cyclic and 2-generated subgroups up to an invariant of the action), verify.py
construction `propus-bordered`, `border_certify.py`.

**Order 260 (v = 65, the smallest open symmetric order).** Settled by a quadruple with no multiplier symmetry, decoded
from Verma's quaternary Legendre pair of length 64 (arXiv:2609.25573; [NOTES.md](../NOTES.md) §32); the structured
searches below could not find it.
- ℤ₆₄, direct joins: every multiplier structure with at most 25 orbits is empty (both row patterns).
- ℤ₆₄, compressions: ⟨31, 33⟩ (25 orbits) has no compatible ℤ₁₆-compression triple for either pattern (ruled out);
  ⟨17⟩ (17/32 orbits): pattern (0, 0, 2) has no compatible ℤ₈ triple, pattern (2, 0, 0) has 640 (pair work 5.2 · 10¹⁰;
  no result is recorded); no result is recorded for ⟨31⟩ (25/33) either. ⟨33⟩ (25/48 orbits) is out of reach: 55,616
  compatible ℤ₈ triples with pair work 3 · 10¹⁸, and H acts trivially on ℤ₁₆ and ℤ₃₂, so no finer compression stays
  small. For the whole ℤ₆₄ family §30 gives at most 10⁻⁵ expected solution orbits.
- Non-cyclic groups of order 64, cyclic automorphism subgroups, at most 30 orbits: ℤ₂ × ℤ₃₂ (6 structures),
  ℤ₄ × ℤ₁₆ (12), ℤ₂² × ℤ₁₆ (14), ℤ₂ × ℤ₄ × ℤ₈ (21), ℤ₈², ℤ₄³ (17), ℤ₂³ × ℤ₈, ℤ₂² × ℤ₄² — no solution. No
  result is recorded for the 2-generated subgroups (≤ 30 orbits) and the quotient chain (`propus_gchain.py`, 26–46
  orbits).

**All census values v ≤ 1000** (197 values): ℤ_{v−1} with at most 30 orbits (direct) and 31–44 orbits (compressions),
and the 319 non-cyclic groups of order v − 1 whose endomorphisms number at most 2 · 10⁶, were queued
(`border_queue.sh`, `bgroup_census.sh`; programs not included). For the census primes these queues searched the
structures of at most 18 orbits exhaustively (§30).

## 28. A symmetric Hadamard matrix of order 876, and small product structures

**Order 876 = 4 · 219** ([symmetric/876](../symmetric/876)). The compressed propus sweep of the composite census
values (`propus_chain.py v --sweep --mino 1 --maxo 36`, smallest structures first) found it in ℤ₂₁₉ with H = ⟨4⟩ of
order 9: H fixes the ℤ₃ coordinate of ℤ₃ × ℤ₇₃ and acts on ℤ₇₃* as the eighth powers, so the structure is the
cyclotomic product ℤ₃ × GF(73) of index 8 (27 orbits for B and D, 14 for A, 13 classes). Row pattern (3, 17, 17):
compressions to ℤ₇₃ give 768 compatible triples (pair work 3.4 · 10⁷ against 6.9 · 10⁹ for the direct join); the
573rd holds 18 solutions, and all 768 together hold 288, in 16 triples (the 16 cosets of H in ℤ₂₁₉* permute them).
Verified by `verify.py`, by `check_matrix` and by an independent rebuild; reproduced from scratch with the same
sequences. Literature status (its README): 219 is undecided in the Đoković–Seberry census of 2025-05-28 and in no
series of it; no paper, database or 2026 deposit we could read has a symmetric matrix of order 876, every classical
route fails by arithmetic, and the 2018 Balonin–Seberry table (Math. Scientist 43(2), Table 2) lists 219 with a dash.

**The propus estimate without the lattice factor.** Without the lattice factor of §29, `propus_rank.py` gives 10^−8.1
solutions for this structure and pattern (pools 10^3.0, 10^6.9, 10^6.9; 13 classes, full covariance), against 288
found: a gap of 10^10.6. With the factor the estimate is 10^0.6, a factor 10^1.9 below the 288 found. For 532 the gap
without the factor is 10^9.2 (§23), and with it the estimate is 10^0.5 against 72. Both structures are products
ℤ_m × GF(q) in which H meets GF(q)* in a cyclotomic subgroup: for 876, H is the eighth powers and fixes ℤ₃; for 532
(ℤ₇ × GF(19)), H = {(h, φ(h))} is the graph of a homomorphism φ: ℤ₇* → ⟨7⟩ ≤ GF(19)* of order 3, so its orbits on
ℤ₇* × GF(19)* are indexed by GF(19)* (an index-18 cyclotomic product with a twist). There the class values are linear
in a few quadratic invariants of short ℤ_m-sequences (the f-layer autocorrelations and the cross-correlations of the
class sequences, weighted by cyclotomic numbers), so they sit on a sparse lattice that a Gaussian with the right
covariance underweights near 0 unless the lattice factor is included. With the factor the excess over the estimate is
10^1.4–10^1.9, and the estimator matches exact counts within ±0.3 on small v (§29).

The sweeps run alongside were the compressed propus sweep of the composite census values (it covered the values up to
435 before it was stopped, §40), the bordered sweeps of §27 (outcome in §27 and §30) and the exact skew chain
(`screen_chain.py 3600 36 0 --order orbits`) over the census-screen survivors (results for at most 33 orbits in §20).

## 29. ℤ_q-bordered propus: symmetric Hadamard matrices of order 4(m + q)

**Construction.** Let G be abelian of order m, φ: G → ℤ_q a surjection with kernel K, and A (symmetric), B, D ±1
functions on G with

    PAF_A + 2 PAF_B + PAF_D = 4m at 0,  −4q on K \ 0,  0 off K.

Let P be the propus array of §17 over G, y_X(c) = Σ_{φ(g)=c} X(g) the compressions to ℤ_q, Mc the propus array over
ℤ_q of (y_A, y_B, y_D) (Mc² = 4q I), and R̂ a Hadamard matrix of order 4q with R̂ Mc R̂ᵀ = −4q X for a ±1 matrix X.
With R[:, (i, g)] = R̂[:, (i, φ(g))],

    H = [[X, R], [Rᵀ, P]]

is a symmetric Hadamard matrix of order 4(m + q): P² = I₄ ⊗ (4(m + q) I − 4q [g − h ∈ K]) and RᵀR = I₄ ⊗ 4q [φ(g) =
φ(h)] cancel; R Rᵀ = (m/q) R̂ R̂ᵀ = 4m I; compressing a group-developed block (also its transpose and its product with
the reflection g ↦ −g) over the fibres of φ gives the ℤ_q-developed block of the compression, so R P = (R̂ Mc)∘φ and X
R + R P = 0 from X R̂ = −R̂ Mc; X² = 4q I. q = 1 is §27. The row sums satisfy a² + 2b² + d² = 4q, all even when v = m
+ q is odd. A census value v is thus reached through every abelian group of order v − q with a ℤ_q quotient, for each
odd q dividing v.
- **Compressions.** z = y/2 is integral with PAF_zA + 2 PAF_zB + PAF_zD = q δ on ℤ_q (z_A symmetric). q = 3: 48
  patterns, single spikes (z_A = ±δ₀ with one B spike, or one B spike and one D spike). q = 5: 120 patterns, and in
  all of them z_B = 0, so B̂ vanishes on the characters of G that are trivial on K.
- **The border problem** is a clique problem: 4q mutually orthogonal ±1 rows r_i (first entry +1) with r_i Mc r_jᵀ =
  ±4q for all i, j. q = 3: 2048 candidate rows, all 48 patterns solved at once (`propus.qborder_solve`;
  Goethals–Seidel arrays of ℤ₃-circulants also work). q = 5: 91,648 candidate rows, mean degree ≈ 3,200; a plain
  depth-first search stalls, a colouring-bounded clique search on each vertex neighbourhood
  (`search/qborder_clique.c`, Tomita's MCQ) solves the patterns met so far in about 150 s (the first after 34 start
  vertices). Goethals–Seidel borders do not exist for the single-spike ℤ₅ pattern (z_A = δ₀, z_D = 2δ₀): none among
  1536 block conjugates × 120,000 GS quadruples. A second version avoids the global adjacency (random adjacent pairs
  (s, t), MCQ on N(s) ∩ N(t), |N(s) ∩ N(t)| ≈ 250): every ℤ₅ pattern met so far is solved in about a second. Simulated
  annealing over the signed column permutations of a fixed Paley matrix of order 20 does not work (cost stuck far from
  0). q = 7: 17.4 million candidate rows, |N(s)| ≈ 370,000, |N(s) ∩ N(t)| ≈ 15,000; 329 start pairs on the pattern
  (z_A, z_B, z_D) = (δ₀, δ₀, (1, 1, 0, 1, 0, 0, −1)) found no clique (open). Hits without a
  border are kept in `qborder_pending.jsonl`.
- **Checks.** ℤ₃: every one of the 576 triples on ℤ₆ gives a symmetric Hadamard matrix of order 36; the chain finds
  and verifies orders 108 (ℤ₂₄), 132 (ℤ₃₀) and 36 on ℤ₂ × ℤ₃; ℤ₁₂ and ℤ₁₈ have no triple (exhaustive, ℤ₁₂ also by
  brute force). ℤ₅: ℤ₁₀ and ℤ₂₀ have none (exhaustive), ℤ₃₀ has several (six compression patterns), and each pattern
  gives a verified symmetric Hadamard matrix of order 140. ℤ₇: sequences exist on ℤ₁₄ and ℤ₂₈ (55 compression
  patterns), but no border is known. Existence depends on m, as for §27.
- **Tools.** `propus.py --qborder q`, `propus_chain.py --qborder q` (per-class targets in stage 1 and in join4x),
  `propus_bgroup.py --qborder q [--phi c₁,…]` (automorphisms that preserve ker φ), `verify.py` construction
  `propus-qbordered` (fields q, quotient = coefficients of φ, R), `border_certify.py` (re-solves the border in the
  verifier's convention: the reflection g ↦ −g, where `propus.py` uses the back diagonal).
- **Literature.** We have not seen this bordering in the propus papers we read (Seberry–Balonin 2017, Balonin–Đoković
  arXiv:1803.08787, Đoković arXiv:2101.05429); it is the symmetric analogue of the bordered Goethals–Seidel arrays of
  §18 with T = all characters of ℤ_q, with the extra requirement that the border be R̂ and R̂ᵀ.
- **Runs.** The 117 census values divisible by 3 (`qborder_pool.sh symQ3.txt … 3`: direct joins up to 30 orbits, the
  chain for 31–44) and the 98 divisible by 5; v = 65 (order 260, the smallest open symmetric order) on ℤ₆₀ (chain,
  1–44 orbits) and on ℤ₂ × ℤ₃₀ (`propus_bgroup.py 2,30 --qborder 5 --phi 0,1 --pairs`). On ℤ₆₀ every structure up to
  24 orbits is empty (direct joins), and the chain found no compatible compression triple in any structure except one
  of 30 orbits (60 triples, no solution). The ℤ₃- and ℤ₅-bordered sweeps of small structures for v > 275 were stopped
  (§30).

**Symmetric Legendre pairs (checked, not pursued).** In the standard Legendre-pair array of order 2ℓ + 2,
[[X, E], [Eᵀ, M]] with X = [[1, 1], [1, −1]], E = [[eᵀ, eᵀ], [eᵀ, −eᵀ]], M = [[A, B], [Bᵀ, −A]] and row sums −1, H is
symmetric exactly when A is (all LPs with A symmetric for ℓ = 3, 5, 9, …, 15 give symmetric Hadamard matrices; this is
known, Đoković–Kotsireas arXiv:1801.07627). For a census value it needs ℓ = 2v − 1 (composite: prime powers are Paley
II). The two published Legendre pairs of length 129 (Kotsireas–Koutschan arXiv:2101.03116, H = {1, 49, 79}) give no
symmetric matrix of order 260: neither sequence satisfies X(r − i) = X(i) for any r. A count of free bits against
constraints (A on ⟨H, −1⟩: ≈ ℓ/2|H| bits, B ≈ ℓ/|H|, about ℓ/2|H| classes) leaves the symmetric LP with a negative
margin where propus has a positive one, so the CPU stays with the propus families.

**Speed.** The stage-1 compressed autocorrelations of `propus_chain.py` are now vectorized (`cpaf_batch`): lists of
757,000 compressions took 15 minutes per pattern in a Python loop (the ℤ₃ chain pass of v = 153 timed out on them);
stage 1 of v = 219 dropped from 18.8 s to 1.5 s.

**The lattice factor in the propus estimate.** The sampled covariance of the class values of a propus structure is
always singular: the row sums fix Σ_c |c| · val_c. `propus_rank.estimate` therefore takes the density on the
hyperplane w^⊥, w = (|c|), times the covolume ∏ g_c · |w| / gcd(g_c w_c) of the lattice points there. Without this
factor the estimates are too low by about Σ_c log₁₀ g_c ≈ C log₁₀ 4. With it the estimator matches exhaustive counts
on ℤ_v with trivial H within ±0.3 in log₁₀ (v = 7, 9, 11, 13, every row pattern; one arithmetic zero at v = 7, rows
(1, 1, 5), predicted 10^1.8). Figures without and with the factor:

| structure | without the factor | with the factor | found |
|---|---|---|---|
| 876: ℤ₂₁₉, ⟨4⟩, rows (3, 17, 17) | 10^−8.1 | 10^0.6 | 288 |
| 532: ℤ₁₃₃, ⟨26⟩, rows (1, 1, 23) / (1, 11, 17) | 10^−7.4 | 10^0.5 | 72 (rows (11, 5, 19)) |
| v = 65, H trivial, rows (11, 3, 11) | 10^−8.3 | 10^10.9 | (not searchable) |

So the two finds sit in structures where about one solution was expected, with a structural excess of 10^1.4–10^1.9,
and the small structures (≤ 30 orbits) we swept mostly have expectations below one. The large ones are where solutions
are plentiful: the propus ranking can be trusted to choose them, and the problem becomes search technology for
structures of 45–80 orbits (multi-level compression for composite v, low-energy lists for prime v), as it was for the
skew orders. The skew estimators (`cp_lowe.estimate_in_lists`, census screen) are separate code that does not use
`propus_rank.estimate`; they are checked separately.

**Calibration of the estimator on the swept structures.** For v = 133 and 219 (all structures up to 30 orbits searched
exhaustively) the estimates stay below 10^0.9, and the two structures that hold solutions are among the highest
(133/⟨26⟩: 10^0.77, 72 found; 219/⟨4⟩: 10^0.57, 288 found; the best empty one, 133/⟨4⟩ with H of order 9, 10^0.81). So
in the swept range each structure holds at most about one expected solution, and a sweep succeeds or fails by chance.
v = 65 is the exception:

| v = 65 structure | orbits A, B/D | estimate (best rows) | searched | found |
|---|---|---|---|---|
| ⟨16⟩ (ℤ₅ × GF(13), index 4) | 13, 25 | 10^2.9 | exhaustive (direct, chain, and an independent numpy hash join over the full pools of rows (11, 3, 11): 420 · 1.76 · 10⁶ · 1.43 · 10⁶) | 0 |
| ⟨51⟩ | 21, 35 | 10^2.0 | exhaustive (chain, d = 13) | 0 |
| ⟨21⟩ | 12, 20 | 10^0.4 | exhaustive | 0: the target is off the coset of the lattice generated by the attainable class-value vectors (rows (3, 11, 3)), an exact obstruction |
| ⟨14⟩ (−1 on ℤ₅, trivial on ℤ₁₃) | 21, 39 | 10^4.1 | no: 6 compatible ℤ₅ triples with pair work 1.05 · 10¹³ | — |

In the cyclotomic products the solution counts are not Gaussian: ℤ₃ × GF(73) (876) holds 70 times its estimate,
ℤ₅ × GF(13) (65/⟨16⟩) none of about 800. The joint-lattice test (`target − s₀` in the ℤ-span of sampled
differences) explains ⟨21⟩ but not ⟨16⟩ or ⟨51⟩, whose zeros must come from a non-linear (norm or cyclotomic-number)
constraint we have not identified. The open v = 65 structure with the highest estimate, ⟨14⟩, is beyond the chain (its
ℤ₁₃ compressions number about 10¹⁰).

## 30. Solution orbits, and dedicated row searches: order 260 over ℤ₆₅ is out

**Count orbits, not solutions.** In a propus structure (ℤ_v, H) the solutions fall into orbits of the group Γ
generated by the multipliers ℤ_v*/H (acting on A, B and D together), the translations of B and of D separately by
the fixed points Fix(H) = {x : hx = x for all h ∈ H} (they keep H-invariance and the PAF; A cannot move, it must stay
symmetric), and the reversals of B and of D separately (the reversal of all three is the multiplier −1). A generic
orbit has |Γ| = |ℤ_v*/H| · 2 · |Fix(H)|² members. The estimator of §29 counts solutions; what decides whether a
structure is worth a search is the expected number of orbits E/|Γ|. A cyclotomic product ℤ_r × GF(p) whose H fixes
ℤ_r has |Fix(H)| = r and loses a factor r².
- 876: the 288 solutions in ℤ₂₁₉/⟨4⟩ are one orbit (|Γ| = 16 · 2 · 9 = 288), where 10^0.57/288 ≈ 0.01 orbits were
  expected; 532: the 72 solutions are two orbits (|Γ| = 18 · 2 = 36). The "structural excess" of §28/§29 was these
  orbit sizes plus selection: several hundred structures with ≈ 0.01–0.1 expected orbits each were swept, and two
  of them paid.
- 65/⟨16⟩: |Γ| = 16 · 2 · 25 = 800, so its 10^3.4 expected solutions (all row patterns) are 2.8 orbits, not 10³;
  ⟨51⟩: 0.29 orbits. The v = 65 "anomaly" of §29 is mostly this miscount.
- Re-ranked (rank_v3, census composites ≤ 391, 10–46 orbits, summed over row patterns): 65/⟨14⟩ 5.3, 65/⟨16⟩ 2.8,
  183/⟨34⟩ 0.69, 213/⟨25⟩ 0.44, 65/⟨51⟩ 0.29, 213/⟨37⟩ 0.08, then 247 (0.13 over all its structures) and 259 (0.09);
  the other 575 structures together hold 0.34 orbits, and none with v > 341 holds more than 10⁻³. The
  sweeps of small structures for v > 275 (plain, ℤ₃- and ℤ₅-bordered) were stopped.
- 2-adic check (`local_density.py`, program not included): for 65/⟨16⟩, rows (11, 3, 11), the fraction of pool triples
  that solve the class equations modulo 2^k, against the uniform 2^−(C−1) per level, is 0.53, 0.35, 0.90 for k = 1, 2,
  3 (a factor of about 1/6) and 0 modulo 16. But modulo 16 the class-value maps are nearly injective (170,948 distinct
  vectors against 170,960 over ℤ), so "none modulo 16" is the integer non-existence itself, not a separate
  obstruction. The 20,800 triples that solve the equations modulo 8 are 26 orbits of size 800 (25 from the ℤ₅
  translations of B and D, 32 from multipliers and reversal) and hit the integer targets in 1–10 of the 12 classes.

**Dedicated row searches.** When v = r · m with gcd(r, m) = 1 and H acts on one factor only, an H-invariant sequence
is a few rows over the other factor, and the class values are row autocorrelations plus cross-correlations of row
pairs; with tables of the pair correlations a class vector is a dozen lookups. The B and D lists are cut by canonical
forms under their own symmetries and by PSD filters (|B̂(χ)|² ≤ 2v, |D̂(χ)|² ≤ 4v, Â(χ)² ≤ 4v), and the join is a hash
join. Each program has a `test` mode that plants random rows, makes their class values the target (PSD filters off,
and multiplier canonicalization off since a planted target is not multiplier-invariant) and must recover the planted
triple; all do. Hits are rebuilt on ℤ_v and checked before they are printed; `search/hit2rec.py` turns them into
`sym_hits.jsonl` records (checked on the 532 solution).
- `search/p65.c` — ℤ₆₅ = ℤ₅ × ℤ₁₃, H = ⟨14⟩ (−1 on ℤ₅, trivial on ℤ₁₃): rows x₀, x₁, x₂ over ℤ₁₃ (X(a, y) =
  x_{|a|}(y)); B and D canonical under the dihedral group of ℤ₁₃ (translations by Fix(H) and reversal), A under the
  multipliers y ↦ μy. The 32 ℤ₅-compression triples (the chain's 6 + 6 + 10 + 10; 16 up to x ↦ 2x) have B lists of
  6·10⁴–1.4·10⁵ and D lists of 10⁶–6·10⁷ after the filters, and joins of at most 4·10⁹ probes, against pair work 10¹³
  per row pattern for the chain. **No solution: ℤ₆₅ with H = ⟨14⟩ holds no propus quadruple.** With ⟨16⟩, ⟨51⟩ and
  ⟨21⟩ empty as well, every multiplier structure on ℤ₆₅ with at most 39 orbits is ruled out, against 8.4 nominal
  expected orbits (about 1.4 if each loses the 2-adic factor measured for ⟨16⟩). What remains on ℤ₆₅ is H = ⟨−1⟩ (B
  symmetric) and H = 1, both far beyond any exhaustive method.
- `search/pz5.c` — the same with ℤ₅ × ℤ_N for N = 12 or 13, optionally with the ℤ₅-border targets of §29. Order 260 as
  a ℤ₅-bordered propus over ℤ₆₀ with H = ⟨49⟩ (−1 on ℤ₅, trivial on ℤ₁₂; 0.37 expected orbits by
  `propus_rank.py --qborder 5`): the compression triples are A = (2, 0, 0), D = (4, 0, 0) and the reverse, with every
  row of B balanced. **No solution** (both triples, 1.8·10¹⁰ probes in all). ⟨11⟩ on ℤ₆₀ (−1 on ℤ₁₂, trivial on ℤ₅;
  0.19 expected orbits) by the chain (`propus_chain.py 60 11 --qborder 5`, 80 compatible ℤ₁₅ triples): **no
  solution**.
- `search/prow.c` — ℤ₃ × GF(p) with H the e-th powers on GF(p) and trivial on ℤ₃ (183/⟨34⟩, 213/⟨37⟩, 213/⟨25⟩): three
  rows of 1 + e bits; the row sums of B and of D are taken sorted (their translations and reversal permute the rows),
  A canonical under the e coset rotations. D is hashed in parts chosen by an additive key of its s = 0 classes (every
  row PAF is 1 mod 4), so all parts together cost one enumeration of D, and each part is probed by every pair (A, B).
  - 213/⟨37⟩ (e = 10): two ℤ₃-compatible triples, exhaustive in two minutes: **no solution**. The row pattern (11, 1,
    27), on which the chain had spent hours (12,140 ℤ₇₁-compatible triples), has no ℤ₃-compatible compression.
  - 183/⟨34⟩ (e = 12): 8 triples, B ≈ 4·10⁶ and D ≈ 4–6·10⁸ per triple, about 2·10⁹ probes each; this gave order 732
    (§31).
  - 213/⟨25⟩ (e = 14): about 3,000 rows per row sum, B ≈ 3·10¹⁰ before the filters per triple: beyond this join.

**Order 260 by the bordered families**, in expected orbits: ℤ₅-bordered on ℤ₆₀, ⟨49⟩ 0.37 and ⟨11⟩ 0.19 (both empty,
above), the rest ≤ 0.02; ℤ₁₃-bordered on ℤ₅₂, ⟨9⟩ 0.22 (20 orbits, but the border is a clique problem with rows of
length 52 that the solver of §29 cannot handle); bordered (q = 1) on ℤ₆₄, at most 10⁻⁵ (the character obstruction
of §27). For the census primes (89, 101, 107, 149, …) the bordered family over ℤ_{v−1} holds at most 0.3 expected
orbits, all in structures of ≤ 18 orbits that the border queues of §27 searched exhaustively.

**Where the solutions are: large structures.** Ranking 47–90 orbits (`propus_rank.py … --mino 47 --maxo 90`) moves the
expected orbits by orders of magnitude, as §25 predicted (below v ≈ 430 the propus count grows with the number of
classes):

| structure | H acts on | orbits A, B | classes | expected solutions | orbit size \|Γ\| | expected orbits |
|---|---|---|---|---|---|---|
| 247/⟨68⟩, 247/⟨87⟩ | order 3 on ℤ₁₃ and on ℤ₁₉ | 42, 83 | 41 | 9.2·10³ | 144 | 64 each |
| 183/⟨13⟩ | trivial on ℤ₃, order 3 on GF(61) | 32, 63 | 31 | 3.4·10⁴ | 720 | 47 |
| 259/⟨121⟩, 259/⟨100⟩ | order 3 on ℤ₇ and on ℤ₃₇ | 44, 87 | 43 | 2.5·10³ | 144 | 18 each |
| 119/⟨69⟩ | −1 on ℤ₇, trivial on ℤ₁₇ | 36, 68 | 35 | 5.6·10⁴ | 27,744 | 2.0 |
| 153/⟨35⟩ | −1 on ℤ₉, trivial on ℤ₁₇ | 45, 85 | 44 | 1.7·10⁴ | 27,744 | 0.6 |
| 189/⟨37⟩, ⟨100⟩ | order 3 on ℤ₂₇ and ℤ₇ | 35, 69 | 34 | 2.9·10³ | 5,832 | 0.5 each |

None of these is within reach. The row form of 183/⟨13⟩ has rows of 21 bits (about 10¹⁵ members of B per row-sum
class), its GF(61) compression leaves about 10¹³ compatible stage-1 triples, and the diagonal actions of 247 and 259
have no row form at all; 119/⟨69⟩ is p65's shape with four rows of 17 bits (10¹⁶ members of B per class). A CP-SAT
model (`search/propus_sat.py`: one Boolean per orbit, every class equation a weighted sum of pairwise XORs) does not
even recover the known 532 solution (2⁶⁶ points, 2 orbits) in 10 minutes, where the join takes seconds; local search
was already negative for spaces of this size (§26). These are the targets for a better search method.

**Symmetric orthogonal designs as plug-in arrays (order 260).** Đoković's symmetric OD(1 + q; 1, q) (arXiv:2508.17141,
Balonin array with circulant blocks, blown up by J_q and the Paley core) suggests building order 260 = 20 · 13 from a
full symmetric OD of order 20 in Balonin form over ℤ₅ with sequences of length 13 plugged into its variables (a propus
solution over ℤ₆₅ with tensor structure, outside every multiplier structure). `search/symod.c` enumerates such designs
(blocks A symmetric, B = C, D; every entry ±x_i; conditions c_ii(s) = 0 and c_ij(s) + c_ji(s) = 0 off 0, with the
pairs where c_ij(s) itself vanishes marked). It returns the propus design (1, 2, 1) for n = 1 and 274 (two variables)
and 76 (three variables) designs for n = 3. For n = 5 there are 304 designs, of types (1, 19) and (2, 18) only, and
none with three variables. Neither admits circulant plug-ins of length 13: (1, 19) needs PAF_X = −19 PAF_Y with
|PAF_X| ≤ 13, and (2, 18) needs PAF_X = −9 PAF_Y, i.e. PAF_Y = +1 at every shift (PAF ≡ 13 mod 4), forcing the row sum
of Y to be ±5 against r_X² + 9 r_Y² = 130. So this route to 260 is closed as well.

**Unsearched small structures.** 247/⟨35⟩, ⟨9⟩, ⟨49⟩, ⟨30⟩, ⟨68, 87⟩, ⟨107⟩, ⟨126⟩ and 259/⟨44⟩, ⟨9⟩ (0.02–0.04
expected orbits each) never appeared in any sweep log; the run with the chain then current had found no solution when
these notes were written; its outcome is not recorded.

## 31. Symmetric Hadamard matrices of order m(m − 1), and orders 5700 and 732

**The construction.** A symmetric orthogonal design OD(m; 1, m − 1) is Y = xP + yN with P a symmetric signed
permutation matrix (P² = I), N symmetric with zeros exactly where P is non-zero, N² = (m − 1)I and PN + NP = 0.
Equivalently (C = NP): S = I + C is a skew-Hadamard matrix and P an involutive signed permutation with PCP = −C, an
anti-automorphism of S of order 2 (then Y = (xI + yC)P is symmetric). Taking x = y = 1, M = P + N is a symmetric
Hadamard matrix of order m (trace 0, since P and N have trace 0 — in each 2-dimensional invariant subspace of the pair
they act as two reflections). Normalising M symmetrically so that row and column 0 are (−1, 1, …, 1) leaves a
symmetric D of order m − 1 with DJ = J and D² = mI − J, and

    H = P ⊗ J_{m−1} + N ⊗ D,   H Hᵀ = I ⊗ (m − 1)J + (PN + NP) ⊗ J + N² ⊗ D² = m(m − 1) I,

is a symmetric Hadamard matrix of order m(m − 1) (`search/sym_mm1.py`, verify.py construction `skew-propus-plug`).
This is Đoković's proof of arXiv:2508.17141, Theorem 2, with two changes: D may come from any symmetric Hadamard
matrix of order m (there it is (I + Q)R with Q the Paley core), and the design need not come from the
Xia–Xia–Seberry–Wu families, so m − 1 need not be a prime power. The plug-in only needs J and D symmetric and
commuting with J² = (m − 1)J, D² = mI − J; for a symmetric OD(kh; k, k(h − 1)) the same argument with J_{h−1} and the
core of a symmetric Hadamard matrix of order h gives order kh(h − 1) (unused: see below).

**Sources of the design.**
- *Paley*: q ≡ 3 (mod 4) a prime, S Paley's skew matrix of order q + 1, P = (−1) ⊕ (x ↦ −x) (the sign on ∞ is
  needed). Checked for q = 3, 7, 11, 19, 23 (orders 12, 56, 132, 380, 552). For q ≡ 7 (mod 8) the order q(q + 1) is
  8 · odd, nothing new for the census.
- *Skew-propus quadruples* (Đoković's Theorem 1): A3 symmetric, A1 used twice, A0 of skew type on ℤ_n (n odd),
  PAF_A3 + 2 PAF_A1 + PAF_A0 = 4nδ, in the Balonin array with D = xI + y(A0 − I). The row sum of A0 is 1, so
  a3² + 2a1² = 4n − 1 is necessary. `sym_mm1.py propus n` searches exhaustively (A1 hashed on its autocorrelations,
  up to shifts; (A3, A0) pairs probed); solutions exist for n = 3, 5, 7, 11, 13, 15 (all verified up to order 3540),
  none for n = 9 (35 is not a² + 2b²). These existence results were published before: Abuzin, Balonin, Đoković and
  Kotsireas (Inform. Control Syst. 2019, no. 5, 2–9, Table 1, "X1 skew, X2 = X3, and X4 symmetric") list such families
  for every odd v < 50 with 4v − 1 = a² + 2b² except 37, 43 and 47, for which their Table 3 has only families with X4
  neither symmetric nor skew; for 37 there is none with X4 symmetric (`sp_mitm.c`, exhaustive, §49). Their families of
  lengths 19, 25 and 31 are equivalent to ours (a multiplier, a shift of A1, and for 31 a reversal and sign change of
  A1).
- *Szekeres-type pairs*: a Legendre pair (M, N) of odd length v' with M of skew type and N symmetric gives the
  bordered two-block skew-Hadamard matrix [[1, 1, e, −e], [−1, 1, −e, −e], [−e, e, A, B], [e, e, −B, Aᵀ]] of
  order 2v' + 2 (orthogonality forces the border signs up to the choices a, α, β with γ = aβ, δ = −aα), and
  P = (border swap) ⊕ diag(R, −R) is an involutive anti-automorphism. `search/szek_pair.c` (cyclic) and
  `search/szek_group.c` (any abelian group) search exhaustively: pairs exist for v' = 5, 9, 13, 21, 29, 33, 41 and
  not for 17, 25 (ℤ₂₅ and ℤ₅²), 37, 45 (ℤ₄₅ and ℤ₃ × ℤ₁₅). For q ≡ 5 (mod 8) prime, M = C₀ ∪ C₁ (quartic classes) and
  N = C₀ ∪ C₂ work for q = 5, 13, 29 (q = s² + 4) but not for 37, so the classical Szekeres difference sets do not
  have a symmetric N in general.

**The census.** v = n(4n − 1) with n odd: 189, 473 and 1139 are in the 2025 census and were settled by
arXiv:2508.17141 (q = 27, 43, 67); 315, 663, 885, 1743, 2475 are not in the census. **1425 = 19 · 75**: the exhaustive
search of length 19 finds two skew-propus quadruples up to equivalence (row sums 5, 5, 1; 7 s), so a symmetric OD(76;
1, 75) and a symmetric Hadamard matrix of order 5700 = 4 · 1425 ([symmetric/5700](../symmetric/5700)). **2093 = 23 ·
91** stays open: 91 = 7 · 13 is not a² + 2b², so there is no skew-propus quadruple of length 23. More generally, for a
Goethals–Seidel skew matrix over ℤ₂₃ whose automorphism group has a normal Sylow 23-subgroup, an involutive
anti-automorphism normalises the ℤ₂₃ action (it normalises the automorphism group, and the Sylow subgroup is
characteristic): it permutes the four blocks with a constant sign on each and acts on each by g ↦ ug + t with the same
u = ±1 on all blocks; the diagonal blocks force u = −1 (a reflection transposes a circulant, so A − I goes to −(A −
I)), and the block equations then force two of B, C, D to have equal autocorrelations (one is a shift of the other or
of its transpose: the transpositions (12), (13), (14), (23), (24), (34)) or are contradictory (the double
transpositions give ε₁ε₃ = 1 and ε₁ε₄ = −1 with ε₃ = ε₄). With b² + c² + d² = 91 only {1, 3, 9} is possible, so this
fails. The Szekeres-type structures over the two groups of order 45 need a symmetric partner (an involutive multiplier
u other than −1 would need u ± 1 invertible, which no involution of ℤ₄₅ or ℤ₃ × ℤ₁₅ satisfies), and there is none.
Order 92 needs another source.

**Beyond the census.** `search/spropus.c` (the same search on the orbits of a multiplier group H of odd order, −1 ∉ H,
hashed on the orbit-level autocorrelations) over every odd n ≤ 267 with 4n − 1 not a prime power and a3² + 2a1² = 4n −
1 solvable, every H with at most 29 orbits: the only solution is n = 31, H = ⟨5⟩ (order 3, 11 orbits), with row sums
5, 7, 1 (q = 123 = 3 · 41 = 5² + 2 · 7²), e.g. A3 = −++−−++−+++−+−+−−+−+−+++−++−−++, A1 =
+++−++++−−++−−−−+−+++−−−+++++−+, A0 = ++++++−+−−+−−+−+−+−++−++−+−−−−− (`sym_mm1.py mpropus 31 5`), hence a symmetric
Hadamard matrix of order 124 · 123 = 15252 = 4 · 3813 (built and checked by `sym_mm1.py build`). The quadruple is
Abuzin–Balonin–Đoković–Kotsireas's (31; 15, 12, 12, 13; 21) over ⟨5⟩ (2019, Table 1) up to reversing and negating
A1, so 15252 is obtainable from published ingredients. Paley II (q = 7625),
XXSW (q = 15251 = 101 · 151) and Đoković's q(q + 1) (q = 123) do not reach 15252, and the census stops at v < 2500.
The Szekeres-type pairs on multiplier orbits (`sym_mm1.py mszek`, every odd subgroup with at most 34 orbits) for v' =
45, 57, 73, 77, 85, 93, 97, 101, 109, 117, 129, 133, 137, 145, 149, 153, 157 found nothing.

**What the plug-in can reach.** If every variable of a symmetric OD of type (s_1, …, s_k) is replaced by a symmetric
±1 matrix X_i of order r with X_i² = α_i I + β_i J (all commuting with J), the result is Hadamard iff Σ s_i β_i = 0.
The off-diagonal Gram value β_i of a ±1 matrix whose Gram matrix is α I + β J satisfies r + β(r − 1) = (row sum)² ≥ 0,
so β_i ≥ −1, with equality exactly for the core of a symmetric Hadamard matrix of order r + 1; β = 0 needs a
symmetric Hadamard matrix of order r (impossible for r odd); X = J has β = r, and 1 ≤ β < r means X = ±(J − 2N)
for a symmetric design N with a polarity. With r odd the cores must balance the positive β: Σ_{β_i > 0} s_i β_i =
Σ_{cores} s_i. With J alone (weight k) this is a symmetric OD(k(r + 1); k, kr) and order kr(r + 1), so for the census
the route gives v = k · n(4n − 1) (r + 1 = 4n); designs in place of J would need symmetric designs whose ±1
matrices commute with a core, which we have not looked for.

**Order m(m + 3): the plug-in J − 2I (order 4828 = 4 · 1207).** (The non-symmetric version is classical: Wallis, J.
Austral. Math. Soc. 11 (1970), Corollary 10 (iv), and Linear and Multilinear Algebra 3 (1975), Corollary 4, which
lists 1207, 1827, 2185, 2575 as new Williamson-type orders; our symmetric matrix is (P ⊗ I) times Wallis's. Found
during the literature search; the same holds for X = J and Williamson's k(k − 1).) In the balance above, X = J − 2I
has β = r − 4, so a symmetric OD(m; 1, m − 1) with x ↦ J − 2I and y ↦ a core of a symmetric Hadamard matrix of order r
+ 1 works when r − 4 = m − 1: H = P ⊗ (J − 2I) + N ⊗ D is a symmetric Hadamard matrix of order m(m + 3) whenever a
symmetric OD(m; 1, m − 1) and a symmetric Hadamard matrix of order m + 4 exist ((J − 2I)² = 4I + (m − 1)J; (m − 1)D² =
(m − 1)(m + 4)I − (m − 1)J). Checked for Paley q = 7, 11, 19, 23 (orders 88, 180, 460, 648; `sym_mm1.py test`). For m
= 4n, n odd, the order is 4 · n(4n + 3); in the census: **1207 = 17 · 71** from Paley q = 67 (m = 68) and the
symmetric Hadamard matrix of order 72 (Paley II for q = 17, doubled) — [symmetric/4828](../symmetric/4828), verify.py
construction `sym-od-plug` — and 2185 = 23 · 95, which again needs a symmetric OD(92; 1, 91). The other balances with
J, J − 2I and cores give r = (m − a + 3b)/(a + b) for weights a on J and b on J − 2I; with a + b = 2 the core has
order 2n or 2n + 4 (≡ 2 mod 4) except for (a, b) = (1, 1), which needs a symmetric OD(4n; 1, 1, 4n − 2): none in
Balonin form (`symod.c` finds only the 3-variable types (2, 8, 2) and (9, 2, 1) at order 12 and none at order 20), and
the Paley designs for q = 7, 11, 19, 23 have no second involution P₂ with P₂N + NP₂ = 2I anticommuting with P
(exhaustive backtracking). (a, b) = (2, 1) would give 391, 1645 and 2255 but needs a symmetric OD(4n; 2, 1, 4n − 3);
in Balonin form that is A3 symmetric with A1 and A0 both of skew type, which forces a3² = 4n − 3, incompatible with n
≡ 2 (mod 3). Designs in place of J − 2I (X = J − 2N, N a symmetric design commuting with the core) would need the
parameters (r, κ, λ) with κ − λ = (r − m + 1)/4, e.g. (71, 21, 6) for 213 and (103, 18, 3) for 1133; not pursued.

**Skew 5700 from the same design (amicable pairs).** A symmetric OD xP + yN is also an amicable pair: S = I + PN is
skew-Hadamard, M = P + N symmetric, and S Mᵀ = M Sᵀ = mP. (It is the symmetric half of a strong amicable orthogonal
design in Seberry's sense, Australas. J. Combin. 55 (2013), whose known orders 2^t, p^r + 1, (n − 1)^r + 1 and
products exclude 76 and 92.)
Normalised at a fixed point of P with sign +1, S and M give amicable cores W (skew, W1 = 0) and D (symmetric,
D1 = −1) with WD + DW = 0, and Seberry's 1978 theorem with m = 1 — a skew-type OD(m; 1, 1, m − 2), here from
G-matrices of order 19 in the Goethals–Seidel array — gives the skew-Hadamard matrix
I ⊗ (I + W) + B ⊗ J + C ⊗ D of order m(m − 1) = 5700 ([skew/5700](../skew/5700), verify.py `skew-amicable-plug`).
1425 is undecided in the skew census as well. The J − 2I trick has no skew analogue here: the plug-in of the weight-1
diagonal variable is I + W, and J − 2I would have to anticommute with W.

The same route works for every skew-propus length n for which G-matrices of order n exist (their row sums need
4n − 2 = s² + t²): the skew census values n(4n − 1) are 473 (n = 11), 663 (13), 1139 (17), 1425 (19), 2093 (23) and
2475 (25). For 473 and 1139 the classical amicable pairs exist (Paley, q = 43, 67) but 42 and 66 are not sums of two
squares, so there are no G-matrices; 2093 lacks the design. **663: skew-propus quadruples of length 13 exist, so
amicable Hadamard matrices of order 52 (not a classical amicable order) and, with G-matrices of order 13, a
skew-Hadamard matrix of order 52 · 51 = 2652** ([skew/2652](../skew/2652)). **2475: the exhaustive skew-propus search
of length 25 finds a quadruple (row sums 9, 3, 1; 18 minutes), so amicable Hadamard matrices of order 100 and, with
G-matrices of order 25, a skew-Hadamard matrix of order 100 · 99 = 9900** ([skew/9900](../skew/9900); not stored,
verify.py rebuilds it). The amicable orders 52, 76, 100 (and 124 from n = 31) are ≡ 4 (mod 8) but not of the
classical form q + 1 with q ≡ 3 (mod 8) a prime power. No Szekeres-type pair of length 49 exists in ℤ₄₉ or ℤ₇²
(`szek_group`), so 100 comes only from the skew-propus quadruple.

**Weighted variant (not realised).** v = k n (4n − 1)-type census values (475 = 5 · 95 with h = 20, 1989 = 3 · 663
with h = 52, 693, 1235, 2185, 2375, and 33k for k = 27, 33, 41, …) would need symmetric OD(kh; k, k(h − 1)) with
k ≥ 3 odd. In Balonin form over ℤ₅ (order 20) only the types (1, 19) and (2, 18) exist (`symod.c`, §30), and a
circulant weighing block CW(n, k) for the x-part needs k to be a square; no such design is known to us.

**Order 732 = 4 · 183** ([symmetric/732](../symmetric/732)). The row search `prow 61 12` (ℤ₃ × GF(61), H = ⟨34⟩, the
twelfth powers on GF(61), trivial on ℤ₃; 20 orbits for A, 39 for B and D, 19 classes; §30) found two solutions in the
fifth ℤ₃-compatible triple (rows (21, 1, 17); A (−1, 11, 11), B (−9, 1, 9), D (−1, 9, 9) by rows; |A| = 464,
B 3.8 · 10⁶, D 4.3 · 10⁸, 1.75 · 10⁹ pairs), in part 3 of 29, after about 1.5 hours with two processes on eight
triples. The estimate of `propus_rank.py` (§29) for the structure is 0.69 expected solution orbits (|Γ| = 432), 0.06
for this row pattern. The remaining triples were not finished; a separate run of the fifth triple alone (commit
`913f7b5`) completes it with exactly the two solutions.

**prow.c part keys (fix).** The D parts of `join3`/`join4` were chosen by the additive key
(F₀ − 3)/4 + 7 (F₁ − 3)/4 mod np of the first two x = 0 classes. Those class values span only a few hundred residues,
so with np in the thousands almost all D landed in a few parts: the 213/⟨25⟩ run (np = 1811, table sized for the
average part) stopped with "part 0 overflows its table" after its B passes. The key is now Σ_k μ_k (F_k − 3)/4 over
all e/2 classes with fixed pseudo-random μ_k (μ₀ = 1), which is still additive over the rows of D and splits into an
A part and a B part for a needed vector. The planted tests (join3 and join4) recover their solutions; 213/⟨25⟩ was
restarted with it.

**Status of 5700, 2652 and 9900** (statuses in the READMEs). Seberry's tables of amicable orders 2ᵗq (1978, 1991,
1992) give t = 3 for q = 13, 19, 25, so they list no amicable pair of order 52, 76 or 100. These pairs follow from the
skew-propus families of Abuzin, Balonin, Đoković and Kotsireas (2019) by the step (I + PN, P + N), which we found
nowhere in print, so 5700 and 9900 were obtainable classically; the first explicit amicable pairs of orders 52, 76 and
100 are ours.
The 1992 summary line "AV 2(q + 1), q ≡ 5 (mod 8); J. Wallis" is the |t| = 1 case
(q = s² + 4) of the line above it; q = 37 = 1 + 4 · 3² gives (|t| + 1)(q + 1) = 152, consistent with the table entry
"19 3". 2652 is the exception: the census lists 663 as undecided, but Seberry–Yamada 1992 (Table 7.1, a blank
meaning t = 2) and Cati–Pasechnik (Table 5, m = 2 not shown) both imply that 4 · 663 is known, without a construction
that we could find (every classical family fails by arithmetic). Both tables also disagree with the census on other
values (783, 879; 451, 921), so we read them as errors, but say so in the README. G-matrices of orders 13, 19, 25
(and 31) are in Georgiou–Koukouvinos, JCMCC 40 (2002) (Đoković, arXiv:0903.5053, Table 1).

### Seberry's Theorem 7 with m ≥ 3: skew 7956 = 4 · 1989

The theorem holds for any m: an OD(mn; 1, m, mn − m − 1) aI + bB + cC and amicable cores (W, D) of order n − 1 give
H = I ⊗ (I + W) + B ⊗ J + C ⊗ D of order mn(n − 1), since (I + W)(I + W)ᵀ + m(n − 1)J + (mn − m − 1)(nI − J)
= mn(n − 1)I (the J terms cancel exactly because of the weight m). `skew_seberry` and the verify.py builder
(`"od"` form) take any m. With n = 4k (k odd) the order is 4 · m k(4k − 1), so the census values in reach are
v = m k(4k − 1), 4k an amicable order, m odd. Among the undecided ones (`search/thm7_census.py`; n = 4 gives every
multiple of 3, but needs an OD of order 4v/3, as hard as the original):

| n (amicable) | v (m) |
|---|---|
| 52 (first explicit) | **1989 (3): done, [skew/7956](../skew/7956)** |
| 20 (q = 19) | 475 (5), 1045 (11), 1235 (13), 2185 (23), 2375 (25); 285 (3) is [skew/1140](../skew/1140), 1425 (15) is [skew/5700](../skew/5700) |
| 12 (q = 11) | 429 (13), 627 (19), 693 (21), 759 (23), 891 (27), 957 (29), 1089 (33), 1221 (37), 1353 (41), 1617 (49), 1683 (51), 1749 (53), 1815 (55), 1947 (59), 2013 (61), 2145 (65), 2277 (69), 2409 (73); 2475 (75) is [skew/9900](../skew/9900) (n = 100) |
| 28 (q = 27) | 1323 (7), 2457 (13) |
| 44 (q = 43) | 2365 (5) |

**m = 3 in Goethals–Seidel form.** The b part must be three single entries in A₂, A₃, A₄. A block with two or more b
entries has off-peak b² autocorrelation that, at total weight 3, nothing can cancel; three entries in one block never
cancel (on a coset of order 3 the off-peak sum is ±1 or 3); and A₁ can only hold paired, skew b entries. Shift each to
0. Write A₁ = aI + cW₀ (W₀ skew) and Aᵢ = bI + cWᵢ. The bc terms then say that on every pair {s, −s} the symmetric
parts of W₂, W₃, W₄ are (0, 0, 0) or a permutation of (1, −1, 0): 20 patterns per pair, an exact invariant of the
moves in `od13.c`. The cross terms between the symmetric and skew parts of one sequence cancel, so PAF(Wᵢ) = PAF(Sᵢ) +
PAF(Kᵢ). Row sums: r₂ + r₃ + r₄ = 0 and Σ rᵢ² = 4t − 4, so with rᵢ = 2ρᵢ, ρ₂² + ρ₂ρ₃ + ρ₃² = (t − 1)/2 must be a value
of the Eisenstein norm form. Annealing finds nothing at t = 5, 11 and 13 (impossible by this condition) and solves t ≤
27 at once. `od13b.c` holds the row sums fixed, with moves that keep them; it still does not reach t = 39 in 20
minutes. `od13x.c` is exhaustive over H-invariant designs: hash the 2^(orbits) skew W₀ on their PAF at orbit
representatives, then enumerate the pattern assignments with the right row sums (one target suffices: permutations of
W₂, W₃, W₄ and a global sign act transitively on them), with the PAF updated from tables of orbit-pair contributions.
Results:

- t = 19, H = ⟨7⟩: solutions.
- t = 27, H = ⟨10⟩: none (1.0 · 10⁷ assignments).
- t = 33, H = ⟨4⟩ (order 5): none.
- **t = 39, H = ⟨16⟩ = {1, 16, 22}**: 2688 solutions out of 4.3 · 10⁶ assignments, in 6 s. This is the largest
  admissible H: it may not contain −1 mod 3 or mod 13.

With the amicable pair of order 52 this gives the skew-Hadamard matrix of order 7956 = 4 · 1989. 1989 is undecided in
both census lists. The symmetric 7956 would need a symmetric OD(156; 3, 153) (the weighted variant above), which is
still open.

**m ≥ 5: skew 1900 = 4 · 475.** The b part is no longer forced. Since the b part of A₁ is skew (row sum 0), m must
be a sum of three squares: m ≡ 7 (mod 8) is impossible in this form (1323, 759, 1815, 2185 are out). For m = 5, the
b parts δ₀ + δ_d, δ₀ − δ_d and δ₀ in A₂, A₃, A₄ have autocorrelation sum 4δ + δ. The bc condition
X(s) = Σᵢ Σₐ uᵢ(a)[Wᵢ(a + s) + Wᵢ(a − s)] = 0 couples W at s and d ± s. Summing it over s gives 2r₂ + r₄ = 0, so the
row sums are (3, 7, −6) at t = 25 and (3, 13, −6) at t = 55. `search/odm_sa.c` (annealing with the b part fixed, the
row sums held, energy Σ cc² + 2 Σ X²) finds **OD(100; 1, 5, 94)** for d = 10 at restart 21 (about 20 s). With the
classical amicable pair of order 20 (Paley, q = 19) this gives the skew-Hadamard matrix of order 5 · 20 · 19 = 1900
([skew/1900](../skew/1900)). 475 is undecided in both census lists, and Table 5 has 475(4).

**The rest of the list is harder.**
- **2365 (t = 55, m = 5).** Annealing stalls at energy about 300. `search/odmx.c` (exhaustive over designs invariant
  under H = ⟨36⟩ of order 5, the only admissible group: lists of c parts per block with fixed row sums, vectors of
  PAF and X at the pair-orbit representatives, 4-list meet in the middle) finds nothing for four b templates (pairs
  at d = 11 and 22; a skew pair in A₁ cancelled by a pair in A₂).
- **429 (t = 39, m = 13; b = a ×3-invariant CW(13, 9) on 3ℤ₃₉ plus δ₀ ± δ₁₃).** The 16 invariant CW(13, 9) exist,
  but H = ⟨16⟩ is incompatible. The orthogonality of row sums forces 3ρ_CW + 2ρ_pair = 0, i.e. ρ_pair ≡ 0 (mod 3),
  while the pair block has one free fixed point and so ρ_pair ≡ ±1 (mod 3).
- The others (1045, 1235, 627, 693, …) have t ≥ 55.

**Multiplier-invariant designs for every b part: skew 9828 = 4 · 2457.** `search/odm_engine.py` automates the search
for a given t, H and m. It enumerates all H-invariant b parts of weight m: ternary combinations of H-orbits, with
block 1 skew, found by a meet in the middle on the PAF vectors at the pair-orbit representatives; the enumeration is
weight-bounded and uses tables of orbit-pair counts. It reduces them modulo multipliers, permutations of blocks 2–4
and signs, screens the row sums (each block's c-part row sum must be a signed sum of its free orbit sizes, and
(0, r₂, r₃, r₄) must satisfy Σ rᵢ² = 4t − m − 1 and Σ rowsum(uᵢ) rᵢ = 0), and runs `odmx` on each survivor.

**t = 91, H = ⟨79, 29⟩ (order 9), m = 13 (v = 2457, n = 28).** There are 18 templates; four of them carry 2
solutions each, in 12 seconds. The b parts live on 7ℤ₉₁ ≅ ℤ₁₃: u₂ = −δ₀ + 1_S, u₃ = −1_{7·{1,3,9}} + 1_{7·{4,10,12}},
u₄ = −1_S with S = 7 · {7, 8, 11}, weights 4 + 6 + 3. With the Paley amicable pair of order 28 (over GF(27);
verify.py now builds Paley pairs over GF(pᵏ)) this gives the skew-Hadamard matrix of order 13 · 28 · 27 = 9828
([skew/9828](../skew/9828)). 2457 is undecided in the skew census. Its Szekeres-type candidate 2(q + 1) has
q = 17³ ≡ 1 (mod 8), outside the series.

**No H-invariant b part at all** (0 templates) for:
- t = 57, |H| = 9, m = 19 (627);
- t = 59, 67, 83, 103 with H the squares (n = 4, m = t: 177, 201, 249, 309);
- t = 67 and 103 with |H| = 11 and 17;
- t = 125, |H| = 25, m = 25 (2375);
- t = 81, |H| = 27, m = 27 (891).

Theorem 7 with n = 4 (m = t = v/3, OD(4t; 1, t, 3t − 1)) reaches every undecided v divisible by 3. With prime t the
squares give only three orbits and no b part of weight t.

**b part on the subgroup of order m (`--bsub`): skew 1716 = 4 · 429.** In both hits the b part lives on the subgroup
of order m (t = km): 7ℤ₉₁ ≅ ℤ₁₃ for 9828 and 3ℤ₃₉ ≅ ℤ₁₃ for 1716. Restricting it there shrinks the template space
(110 763 candidate b parts per block at t = 39 without the restriction, 243 with it). With `--bsub 13`,
`odm_engine.py 39 16 13` gives 10 templates and 18 solutions in 3 s (the b parts u₂ of weight 10 and u₃ = −1 on the
×3-orbit {1, 3, 9}, in units of 3). With the Paley pair of order 12 this gives the skew-Hadamard matrix of order
13 · 12 · 11 = 1716 ([skew/1716](../skew/1716)). 429 is undecided in the skew census, and Table 5 and Table 7.1 both
give t = 3. The run at t = 91 with `--bsub 13` also finds solutions (10 templates).

With `--bsub` there are no templates (0) for 1045 (t = 55, |H| = 5, m = 11) and 957 (t = 87, |H| = 7, m = 29).

**Eades' condition prunes the list.** A 4-type OD(4t; 1, m, c) (four commuting blocks, as in Goethals–Seidel) needs
an integer matrix E with EEᵀ = diag(1, m, c) (Handbook of Combinatorial Designs, Thm 2.96). With the unit row e₁ this
means orthogonal r, ρ ∈ ℤ³ with |r|² = m and |ρ|² = c = 4t − m − 1: the row sums of the b and c parts, which is the
engine's row-sum screen. Example: for 627 (OD(228; 1, 19, 208)), 19 forces r = (3, 3, 1) up to signs and order, 208
forces ρ = (12, 8, 0), and no pair is orthogonal. That is why all 97 templates at t = 57 failed the screen.

The condition fails for 627, 759, 1045, 1221, 1235, 1323, 1425, 1683, 1815, 2145, 2185 and 2375: no 4-type design
exists for these. (1425 is [skew/5700](../skew/5700) by another route.) It holds for 429, 475, 693, 891, 957, 1089,
1353, 1617, 1749, 1947, 1989, 2013, 2277, 2365, 2409 and 2457. Of these, 429, 475, 1989 and 2457 are done; the rest
are open.

**Harvest table (Theorem 7, Goethals–Seidel form over ℤ_t, H-invariant).** "templates" are the H-invariant b parts of
weight m with autocorrelation sum mδ, up to equivalence. The c parts are exhaustive over H-invariant blocks
(`odm_engine.py`; `odm_join.py` with join4x for large lists; `--bjoin` builds the templates with join4x from numpy
candidate arrays).

| v | t | m | n | H | b part | outcome |
|---:|---:|---:|---:|---|---|---|
| 1989 | 39 | 3 | 52 (first explicit) | {1, 16, 22} | three singletons | **2688 solutions → skew 7956** |
| 475 | 25 | 5 | 20 | none (annealing) | δ₀ ± δ₁₀, δ₀ | **found → skew 1900** |
| 429 | 39 | 13 | 12 | {1, 16, 22} | on 3ℤ₃₉ | **18 solutions → skew 1716** |
| 2457 | 91 | 13 | 28 | ⟨79, 29⟩, order 9 | on 7ℤ₉₁ (found also unrestricted) | **8 solutions → skew 9828** |
| 2365 | 55 | 5 | 44 | ⟨36⟩, order 5 (the only one) | all 3 templates | none |
| 627 | 57 | 19 | 12 | order 9: no template; order 3: no template passes | — | Eades fails |
| 957 | 87 | 29 | 12 | ⟨7⟩, order 7 (the only one) | no template (110 763 candidates per block) | none |
| 1749 | 159 | 53 | 12 | ⟨10⟩, order 13 (the only one) | no template | none |
| 2013 | 183 | 61 | 12 | ⟨16⟩, order 15 | no template | none (order 5: 39 orbits, too large) |
| 1947 | 177 | 59 | 12 | ⟨4⟩, order 29 (the only one) | no template | none |
| 1089 | 99 | 33 | 12 | ⟨4⟩, order 15 | 278 templates (b on 33ℤ) | no c part |
| 177, 201, 249, 309 | 59, 67, 83, 103 | t | 4 | the squares (and orders 11, 17) | no template | none |
| 2375 | 125 | 25 | 20 | ⟨6⟩, order 25 | no template | Eades fails anyway |
| 891 | 81 | 27 | 12 | ⟨4⟩, order 27 | no template | none (order 9: no design, capped search, §33) |
| 693 | 63 | 21 | 12 | ⟨4⟩ = {1, 4, 16} | on 3ℤ₆₃; 1848 templates from 200 000 of 1 054 944 b quadruples (cap) | **4 solutions → skew 2772** |
| 1617 | 147 | 49 | 12 | ⟨4⟩, order 21 | 3169 templates from 200 000 of 81 450 144 (cap) | none passes the row-sum screen |

**Order 92: no Goethals–Seidel skew matrix with a block-monomial involutive anti-automorphism.** A symmetric
OD(92; 1, 91) would give three census values: skew 4 · 2093 (amicable 92 with G-matrices of order 23), symmetric
4 · 2093 (P ⊗ J + N ⊗ D) and symmetric 4 · 2185 (X = J − 2I). It is equivalent to a skew-Hadamard matrix S of order
92 with an involution P (a symmetric signed permutation) such that PSP = Sᵀ.

`search/gs_anti.py` takes S of Goethals–Seidel form over ℤₙ and P block monomial. On the diagonal, Q A Qᵀ = Aᵀ for a
skew-type A forces every block of P to be a reflection R_c: g ↦ c − g. That leaves an involution π of the blocks,
signs ε_{π(i)} = ε_i and centres c_{π(i)} = c_i; block translations with zero sum normalise all centres but one.
Every entry of PSP = Sᵀ is a signed equation between sequence entries, and a signed union-find gives the induced
structure or a contradiction. π = identity is always contradictory (it needs ε_iε_j = −1 for every pair).

| n | configurations | free classes | GS solutions |
|---:|---:|---|---|
| 3, 5, 7 | 36, 60, 84 | (n + 3)/2 | 288, 960, 3360; each checked to give a symmetric OD(4n; 1, 4n − 1) |
| 9 | 108 | 6, 9, 18 | 0 |
| 11 | 132 | 7 and 22 | 10 560 (small configurations) |
| 13–21 | — | (n + 3)/2, and 2n for 12 configurations | 0 in the small ones |
| 23 | 276 | 13 (264 configurations), 46 (12) | 0 in all 264 × 2¹³ |

The 12 large configurations always have the skew-propus layout: A of skew type, B symmetric, and D tied to C by a
reversal, so PAF(A) + PAF(B) + 2 PAF(C) = 4nδ. That problem was exhausted for n = 23 earlier (§31, no skew-propus of
length 23). So the order-92 question needs S outside the Goethals–Seidel form or P outside block-monomial
reflections.

**Tooling for the larger cases.** `search/odm_join2.py` replaces odmx's in-memory join, which rebuilt every c list
from scratch per template: about 25 minutes per template at t = 63 with |H| = 3.
- `search/odmc.c` builds each list of H-invariant c parts once per (block type, b part, row sum). A depth-first
  search over the orbit variables updates the PAF and bc vectors incrementally from orbit-pair tables, and visits only
  branches that can still reach the row sum.
- The lists are cached across templates, and each template's row-sum vectors are joined by `join4x.c`.
- Validation: at t = 39 it reproduces odm_join's counts exactly (120 per row-sum vector for template 0), and a decoded
  solution checks as an OD.

## 32. The symmetric weighted plug-in in Balonin form; generic group engines

**The plug-in.** A symmetric OD(4t; μ, 4t − μ) xP + yN (P, N symmetric with complementary supports, P² = μI,
N² = (4t − μ)I, PN + NP = 0) and a symmetric Hadamard matrix of order h with 4t = μh give the symmetric Hadamard
matrix P ⊗ J_{h−1} + N ⊗ D of order μh(h − 1), D the core of the symmetric matrix (§31: the J terms cancel exactly
when the weight of P is μ = 4t/h). No amicable pair is needed, only a symmetric matrix of order h = 4k, so on paper
the census values v = μk(4k − 1) are in reach: 125 symmetric-undecided values in 139 splits, 122 in 134 without the
values already settled here (k = 1, h = 4 alone covers every multiple of 3), e.g. 475 (t = 25, μ = 5, h = 20), 891
(t = 81, μ = 27, h = 12), 1989 (t = 39, μ = 3, h = 52).

**Balonin form.** With A symmetric and B, D any circulants (or group-developed matrices),

    [ −A   BR   BR   DR ]
    [ BR   RD   A   −RB ]
    [ BR   A   −RD   RB ]
    [ DR  −RB   RB   A  ]

is symmetric and orthogonal whenever AAᵀ + 2BBᵀ + DDᵀ is scalar (the off-diagonal blocks cancel identically, also with
variables, provided both parts of A are symmetric). Splitting each block into an x part and a y part, the design
conditions are PAF(A_x) + 2PAF(B_x) + PAF(D_x) = μδ, the same for y with 4t − μ, and X_A + 2X_B + X_D = 0 for the
cross terms X(s) = Σ_g x(g)[y(g + s) + y(g − s)]. Summing the cross condition over s gives the row-sum screen r_Aρ_A +
2r_Bρ_B + r_Dρ_D = 0 with r_A² + 2r_B² + r_D² = μ and ρ_A² + 2ρ_B² + ρ_D² = 4t − μ (`search/symod_reach.py`). It
removes 1989 (μ = 3), 693 (μ = 21), 1235, 1353, 2013, 2093 (μ = 1, as in §31) and most of the n = 4 splits; it keeps
475, 891, 1089, 2185, 2375, 1617, 1749, 1815, 1947, 2277, 513, 657, 669 and 723.

Small x parts are rigid. For μ = 3 the cross condition forces B and D to be of skew type around their single x
entries, so PAF(A) + 2PAF(B) + PAF(D) = 4tδ with A symmetric and a² = 4t − 3 for the row sum of A.

**Exhaustive results (`search/symod_g.py`, H trivial, every Balonin design over the group).**

| t | μ with a design | μ without one |
|---:|---|---|
| 3 | 1, 3 | — |
| 5 | 1 | 3, 5 |
| 7 | 1, 3 | 5, 7 |
| 9 | — | 1, 3, 5, 7, 9 |
| 11 | 1 | 3, 5, 7, 9, 11 |
| 13, 15 | 1 | 3, 5, 7 |
| 17 | 1 | 3, 5 |
| 19 | 1 | 3 |
| 21 | 1 | 3, 7 |
| 25 (ℤ₂₅ and ℤ₅²) | — | 5 (symmetric 475) |

No design with μ ≥ 5 exists in any of these cases, and μ = 3 only for t = 3 and 7. So the weighted variant in Balonin
form looks barren; 475 is out over both groups of order 25 (2 and 1 x templates up to equivalence, one row-sum class
each, joins empty). The H-invariant run for 891 (ℤ₈₁, H = ⟨10⟩, x part on 3ℤ₈₁) found 151 templates in a capped sample
of 200 000 x triples and no design.

**Generic engines.** `search/gsgroup.py` (groups ∏ ℤ_n^k, automorphisms as matrices, orbit and cross tables),
`search/gdfs.c` (depth-first lists with incremental PAF and cross vectors), `search/symod_g.py` (symmetric Balonin)
and `search/odm_g.py` (skew Goethals–Seidel, Theorem 7) work over any such group with H of odd order:
- templates by join4x;
- canonical forms by row-wise lexicographic minima over the normaliser, block signs, and translations by H-fixed
  elements and negation of the non-skew blocks;
- `--hb` puts the b (x) part under a larger group than the c part;
- `--b4w` bounds the b weight of one of blocks 2–4, which turns the template join into a three-list join;
- `--anneal` replaces the exhaustive c lists by `search/osa.c`, annealing over the orbit variables with the row sums
  held.

Checks: the engines reproduce the 1716 designs (12 solutions over ℤ₃₉ with H = ⟨16⟩, b on 3ℤ₃₉; the annealer finds one
too), μ = 1 skew-propus designs, and the ℤ₂₅ negative. In the engines the signed row-sum size of a skew pair is 0, not
its weight. verify.py's Theorem-7 "od" form takes an optional "group" (blocks over ∏ ℤ_n^k, R the negation).

**Theorem 7 with n = 4.** A skew OD(4t; 1, t, 3t − 1) gives a skew-Hadamard matrix of order 12t, so t = v/3 reaches
every undecided skew v divisible by 3 (172 of them). Exhaustively over ℤ_t (H trivial, Goethals–Seidel form):
- t = 3, 5, 9: designs exist;
- t = 7: impossible, since 7 is not a sum of three squares (Eades), as for every t ≡ 7 (mod 8).
- t = 11: a design exists (template from a capped sample of 2 · 10⁶ b quadruples).
- t = 13: no result is recorded.
In every design found the b parts sit in a short window around 0, where periodic and aperiodic autocorrelations agree.
They are four short ternary complementary sequences with u₁ antisymmetric. For t = 11:
- u₁ = (1, 0, −1) on positions (−1, 0, 1);
- u₂ = −(δ₀ + δ₁ + δ₂);
- u₃ = (1, −1, −1) on (−1, 0, 1);
- u₄ = (−1, 1, −1) on (−2, −1, 0).

Their aperiodic autocorrelations at shifts 1 and 2 sum to 0, and the weights sum to 2 + 3 + 3 + 3 = 11. So b parts
of any weight are plentiful (from complementary ternary sequences); the c part is the hard half.

The census sweep (`n4sweep.py`; program not included) covers every cyclic H of odd order with at most 15 orbits for
the 120 targets with t a sum of three squares: 3-orbit groups (squares) and most 5–9-orbit groups have no b template
at all, and no c part was found in any run.

**Theorem-7 targets: what the H-invariant searches can and cannot reach.**
- Non-cyclic groups, exhaustive. Over ℤ₃² × ℤ₁₁ (H = unipotent × ⟨3⟩, order 15, 15 orbits) 1089 has 5 templates
  (88 128 b quadruples, complete) and no c part. Over ℤ₃² × ℤ₂₃ (order 33, 15 orbits) 2277 has no b template. Over
  ℤ₃ × ℤ₇² (order 21, 15 orbits) 1617 has 42 templates in a capped sample of 2 · 10⁵ quadruples, and none passes the
  row-sum screen. Over ℤ₃⁴ with the order-13 automorphism (9 orbits) 891 has no b template.
- Cyclic groups with 27 orbits per block (1353/⟨10⟩, 1089/⟨91⟩, 2409/⟨4⟩, 2277/⟨163⟩): the c lists reach about
  2 · 10⁷ per block. A single template's join ran for more than 40 minutes (1353, 233 templates), so these are beyond
  exhaustive joins. With one b block of weight ≤ 3 there is no template at all for 2409/⟨4⟩ and 2277/⟨163⟩.
- `search/odm_estimate.py` (Gaussian model: exact list sizes times the per-coordinate probability of a zero sum)
  under-counts badly. The known hits sit 3 (1716) to 5 (2772) orders of magnitude above the model, so it cannot rule
  a structure out; for 1353/⟨10⟩ it gives 10⁻¹¹ per template and for 1089/⟨91⟩ 10⁻¹⁵.
- Orbit annealing (`osa.c`) finds the 1716 design. For 1089 with H of order 3 (55 H-orbits, about 52 variables per
  block) the energy never drops below about 10⁵ (random level), at T₀ = 40 or 3000. This agrees with §26: local search
  is not competitive above about 40 free positions per block.

**A second plug for Theorem 7 (negative).** In H = I ⊗ (I + W) + B ⊗ X + C ⊗ D, any symmetric ±1 X with WX = −XW,
XD = DX and X² = αI + βJ works, and the weights must satisfy c = m(r − 4n') − 1, N = m(r + 1 − 4n') for
X² = 4n'I + (r − 4n')J (r = n − 1). X = J (n' = 0) is Theorem 7. X = J − 2P (n' = 1, P an involutive permutation)
would give skew orders m r(r − 3), e.g. skew 4 · 575 from the OD(100; 1, 5, 94) of 1900 with the Paley cores of
order 23, and skew orders 28t from the n = 4 designs. But no permutation P with PWP = −W and PDP = D exists for the
Paley cores, q = 7, …, 47 (exhaustive backtracking). Among the back-circulant ±1 matrices on GF(q), which all
anticommute with the Jacobsthal W, the only ones commuting with D = (I + Q)U with a two-level square are ±J and the
core itself (q = 7, 11, 19, 23, all 2^q sets). So the plug needs amicable cores with a different symmetry.

**Two-block skew matrices with a general second block (`search/szek2.py`).** The bordered array
[[1, 1, e, e], [−1, 1, e, −e], [−eᵀ, −eᵀ, A, B], [−eᵀ, eᵀ, −Bᵀ, Aᵀ]] is skew-Hadamard of order 2(w + 1) whenever A is
of skew type, both row sums are 1 and PAF_a + PAF_b = −2 off 0. H + Hᵀ = 2I only uses A + Aᵀ = 2I, because the
off-diagonal blocks pair B with −Bᵀ; every row product reduces to the PAF condition and the row sums. So b need not
be symmetric. `szekeres.py` (§15) required it (the classical Szekeres difference sets), and 767 of its 988 census
structures had no symmetric b at all.

With b only H-invariant the free signs per class equation rise from 2 to 3. The symmetric analogue (A symmetric, b
any, row sums −1, the Legendre-pair array of §29) has the same count. Exhaustively with H trivial:
- skew: w = 5, 9, 13, 21 give matrices, 17 and 25 none;
- symmetric: w = 5, 9, 13, 17, 21 give matrices, 25 none.

The census sweeps (w = 2v − 1 for every undecided v, cyclic ℤ_w, every odd-order cyclic H with both lists ≤ 3 · 10⁷)
are complete and empty: 631 skew structures and 504 symmetric ones, no hit (most symmetric ones have no symmetric a
with row sum −1 at all, nA = 0). The Kiyasu routes are empty in this form for the largest groups: 235 (ℤ₄₆₉, |H| =
33), 277 (ℤ₅₅₃, |H| = 39) and 149 (ℤ₂₉₇, |H| = 45). verify.py has a construction "two-block" (kinds skew and
symmetric, any abelian group) ready for certificates.

## 33. Plug-in balances, two weight-1 symmetric variables, and cyclotomic plugs

**Theorem 7 outside Goethals–Seidel form (negative).** Eades' condition (§31) only rules out 4-type designs. A
skew-type OD(n; 1, m, c) of order n ≡ 4 (mod 8) in any form gives a rational family, which exists only if
⟨1, m, c, mc⟩ ≅ 4⟨1⟩, i.e. the Hilbert symbols (−m, −c)_p are 1 at every odd prime p (the transpose is the canonical
involution of the quaternion algebra (−m, −c) generated by the two skew variable matrices, so the trace form is a
multiple of its norm form). For the twelve Theorem-7 census targets where Eades fails this test fails as well:
627 (p = 13), 759 (7), 1045 (13), 1221 (7, 29), 1235 (3, 41), 1323 (47), 1683 (7), 1815 (151), 2145 (7, 17), 2185
(109), 2375 (3, 79). So no array of any kind (two-block, non-abelian groups, bordered) can give these designs. The
test passes for 429, 475, 693, 891, 957, 1089, 1353, 1617, 1749, 1947, 1989, 2013, 2277, 2365, 2409 and 2457, exactly
Eades' list. (program not included; the same test for symmetric designs is not valid: their transpose is an orthogonal
involution, and the symmetric OD(76; 1, 75) of §31 would fail it.)

**No symmetric OD(4n; 1, 1, 4n − 2) for n > 1.** In the symmetric plug-in (§31) the balance s_J r + s_P (r − 4) =
s_core with weight 1 on J, weight 1 on J − 2I and the rest on cores gives r = 2n + 1, i.e. symmetric Hadamard matrices
of order 4n(2n + 1) (census: 903, 1953) from a symmetric OD(4n; 1, 1, 4n − 2) and any symmetric Hadamard matrix of
order 2n + 2. No such design exists. Let P₁, P₂ be its weight-1 matrices (symmetric signed permutations with
underlying involutions τ₁, τ₂, which commute because P₁P₂ = −P₂P₁) and N the rest, N anticommuting with both. Evaluate
NP₁ + P₁N at (i, τ₂(i)) and NP₂ + P₂N at (i, τ₁(i)): with x = N[i, τ₁τ₂(i)] and y = N[τ₁(i), τ₂(i)] they give x
s₁(τ₂i) = −s₁(i) y and x s₂(τ₁i) = −s₂(i) y, while P₁P₂ = −P₂P₁ gives s₁(i)s₁(τ₂i) = −s₂(i)s₂(τ₁i). So x = 0. Since N
vanishes only on the supports of P₁ and P₂, every i is fixed by exactly one of τ₁, τ₂. On X₁ = Fix(τ₂), P₂ is a
diagonal ±1 matrix D₂ with D₂(τ₁i) = −D₂(i), and NP₂ = −P₂N forces N[i, j] = 0 whenever D₂(i) = D₂(j); but N[i, j] ≠ 0
for all j ∉ {i, τ₁(i)} in X₁, so |X₁| ≤ 2, likewise |X₂| ≤ 2, and the order is at most 4 (where the design exists).
This also explains the empty Balonin-form searches for this type (§31).

**Weighted cyclotomic plugs (negative).** Generalising Theorem 7's plug, a skew-type OD(4t; 1, s₂, …, s_k) with
I ⊗ A + Σ X_i ⊗ Y_iR, A group-developed of skew type a and Y_i of symmetric y_i over GF(p), is skew-Hadamard of order
4tp iff PAF_a + Σ s_i PAF_{y_i} = 0 off 0 (R-type symmetric blocks commute with each other and are amicable with any
group-developed A, so nothing else is needed). With a, y_i constant on cyclotomic classes of order e (−1 ∈ C_{e/2})
the condition is e/2 linear equations in the class PAF vectors; a meet in the middle over all sign patterns
(program not included; primes p ≤ 700, e = 4, 6, 3- and 4-variable types of orders 4–28 with square products)
finds only p = 5, 7, 13, 19 (orders up to 4 · 91). The class PAFs are q/e² (Σ f)² plus a few Jacobi-sum
parameters times √q, so the equations force relations like q = s² + 4 (the classical families) rather than behaving
like independent Gaussian constraints. The 3-sequence version with free weights (m, c solved from the equations)
gives only q = 5 and 7.

**Theorem 7 for 891 with b on 3ℤ₈₁ (negative, capped).** `odm_g.py 81 10 27 --bsub 27` with b-block weights
(9, 3) and (12, 5) in blocks 3 and 4 (template joins capped at 2 · 10⁵ b quadruples): 5 and 9 templates up to
equivalence, 0 and 8 pass the row-sum screen, and the exhaustive c joins (39 distinct c lists) give no design.

**Difference map for spin structures (`search/spin_rrr.py`, `search/spin_check.py`).** In Fourier form the spin
condition reads |Â(k)|² + |B̂(k)|² + |B̂(μk)|² + |B̂(μ²k)|² = 4v, so per K-class the three H-orbit values of B̂ form a
vector in ℂ³ of prescribed norm 4v − |Â|² (A from the precomputed pool of K-invariant skew blocks). RRR alternates
this projection (rescale each triple, best A per step) with the ±1 projection (top-k rounding at the forced row sum).
It finds Đoković's spin structure for 397 (|H| = 11, μ = 34, 6 classes) in 18 200 iterations (seconds; rebuilt and
checked as a skew-Hadamard matrix of order 1588 by `spin_check.py`), but not 631 (7 classes) in 10 minutes or 919
(9 classes) in about 45 restarts of 3 · 10⁴ iterations, which the tabu search of `spin_ls.c` solves in about a
minute and twenty minutes. Treating A as continuous variables too (one 7-vector per class) fails even for 397. The
general skew Goethals–Seidel version (four independent blocks, 7-vector per character class) does not find 353/⟨22⟩
(16 classes, 10⁹ solutions) in 20 restarts. So the projection method is not a substitute for the tabu search here.

**Faster spin tabu (`spin_ls2.c`, `SPIN_LS=spin_ls2 spin_skew.py csearch ...`).** The same moves and acceptance rule
as `spin_ls.c` (identical trajectories for a given seed: the best residual after 60 restarts on 1657 is 26 496 in
both), but the swap scan evaluates, for each + orbit i, all partner orbits j at once in single precision with the
j-dependent terms (2g_cj + 2M_cjj and the row i of M_c) laid out contiguously, so the compiler vectorises it; the
exact integer objective is recomputed after every move. 2.75 times faster: 60 restarts of 3000 steps on 1657 take 4 s
instead of 11 s. Measured restarts per solution on 919 (|H| = 17, μ = 8, 9 classes): 1361, 3753, 8848 (steps 3000),
328 and 1695 (steps 8000), 14 543 (steps 1000), i.e. one to five minutes per solution on a loaded machine, against
twenty minutes before; on 631 (|H| = 7, μ = 2, 7 classes) 7–50 restarts (under a second). Two variants did not help on
631: jumping A to the best of the 2^NP pool every 10 or 50 steps, and restarting half the time from the best state so
far perturbed by six swaps. The census queue run with it (817 and 721 on their 8-class structures, then 1657, 1801 and
1951 in two-hour slices) found no solution (817 and 721 are infeasible, §36; earlier `spin_ls` runs on 1657, 1801 and
1951 found nothing, §14).

## 34. Spence's skew 4(1 + q + q²) construction made explicit: it needs a cyclic CDS

Spence (Canad. J. Math. 27 (1975) 555–560) builds a skew-Hadamard matrix of Goethals–Seidel type of order
4(1 + q + q²) from a cyclic projective plane of order q and complementary difference sets (CDS) in ℤ_v,
v = 1 + q + q². The plane part, as implemented in `search/spence_plane.py`:
- Q = q², N = Q² + Q + 1 = v w with w = q² − q + 1 (coprime to v). D = {i mod N : Tr(αⁱ) = 0} is the Singer
  (N, Q + 1, 1) difference set of PG(2, Q); it is fixed by the multiplier q.
- f(r) = #{d ∈ D : d ≡ r mod v} takes values 0, 1, 2 (the fibres are orbits of the subgroup of order w, arcs).
  Every difference of D occurs once, so Σ_r f(r) f(r + s) = w for s ≢ 0, and u = f − 1 has PAF_u = w − 2(Q + 1) + v =
  0 off 0: a perfect ternary sequence.
- With D₁ = {f = 1} and D₂ = {f = 2}, R = +1 on D₁ ∪ D₂ and S = +1 on D₂ are R = u + 1_{D₁}, S = u − 1_{D₁}, so
  PAF_R + PAF_S = 2 λ_{D₁}(s). D₁ is a planar (v, q + 1, 1) difference set, so the sum is the constant +2. Row sums 1
  and −(2q + 1).
- Checked for q = 2, 3, 4, 5, 7, 8, 9 and 23. With the Paley CDS (a = b = Legendre, a(0) = 1) the Goethals–Seidel
  matrix GS(a, b, R, S) is skew-Hadamard for q = 2 (order 28) and q = 5 (order 124).

(R, S) is the classical D-optimal pair of order 2v. So Spence's theorem is "CDS + D-optimal pair over the same group".
The census-undecided Spence values are v = 553 (q = 23), 1057 (32), 1407 (37), 1893 (43) and 2257 (47). For none of
them is there a classical CDS: v is not a prime ≡ 3, 5, 7 (mod 8), and 2v + 1 = 1107, 2115, 2815, 3787, 4515 are not
prime powers. All five v are squarefree, so ℤ_v is the only abelian group of order v.

**CDS search (negative).** A CDS in ℤ_v is a two-block skew pair (a of skew type, b any, PAF_a + PAF_b = −2), which
is exactly what `szek2.py` searches. It was run on every odd-order multiplier subgroup with at most 34 H-orbits:
- the structures: 553 with |H| = 117 and 39 (four), 1057 with |H| = 225 and 75 (four), 1407 with |H| = 99, 1893 with
  |H| = 315 and 105, 2257 with |H| = 135;
- this includes the non-cyclic H of orders 117 and 225 that the census sweep (§33) skipped;
- all are empty;
- 1057/⟨2, 85⟩ (|H| = 45) and 1893/⟨91⟩ (|H| = 63) have b lists of 1.2 · 10⁸ and 4.7 · 10⁸ and were skipped.

The heuristic count says why. Per ⟨H, −1⟩-class a CDS has 3 free signs (one for a, two for b). The class equation
costs log₂(√(πv)/2) ≈ 4.4 bits at v = 553: a sum of two PAFs has variance about 2v on the lattice 4ℤ + 2. So the
expected count falls like 2^(−1.4 · classes) for every multiplier structure (about 10⁻³ at four classes), and only
algebra can supply these pairs. The same holds for any twist of the plane pair: signs y on D₁ give
PAF_R + PAF_S = 2ε(s) with ε(d − d′) = y(d)y(d′), which only changes the target of the same two-block problem. The
route stays closed for these v unless a CDS family for v = pq (or 3pq) appears.

## 35. Three skew blocks in one multiplier orbit (negative for the census)

**The structure.** Take A of skew type and H-invariant, μ with μ³ ∈ H, K = ⟨H, μ⟩ with −1 ∉ K, and D K-invariant.
The blocks A, A∘μ, A∘μ², D have autocorrelation sum T_A + P_D, where T_A(s) = P_A(s) + P_A(μs) + P_A(μ²s). Both
terms are K-invariant, so the condition is an exact join of two lists over the ⟨K, −1⟩-classes: the 2^pairs skew
blocks A, and the K-invariant D with row sum d. The row sums are (1, 1, 1, d), so d² = 4v − 3: v = k² + k + 1 and
d = 2k + 1. This is the pattern of 4492 and 5932, but with the three skew blocks tied. The census values of this form
are 553, 1057, 1123, 1261, 1333, 1407, 1483, 1561, 1893, 1981, 2071, 2163, 2257 and 2353.

**Tools.**
- `tri_skew.py`: the exact join over any abelian group (gsgroup tables, gdfs lists, join4x), and `--sweep v` over
  every odd H ≤ ℤ_v* and every admissible μ.
- `tri_ls.py`: a Gaussian estimate (`rank`), and the tabu search of `spin_ls2.c` with the roles exchanged (its
  swap block is D, its pair block is A with the quadratic form of T_A).

Both are validated on small v:
- The joins with H trivial give skew matrices of orders 28, 52, 124 and 172 (v = 7, 13, 31, 43).
- There are none at v = 21, nor at 57, 73 and 91 with |H| = 3.
- The tabu search solves v = 43 at once.

**Census (negative).**
- 553:
  - |H| = 13 with μ = 18 (10 classes): nA = 2²⁴, nD = 17 424, exact join empty, Gaussian estimate 10⁻⁶·⁰.
  - The other order-13 and order-39 structures have no D with row sum 47.
  - |H| = 3 (40 classes): 10⁻¹⁴·⁸ to 10⁻⁷·⁹.
- 1123 (|H| = 11, 17 classes): 10⁻²·⁹.
- 1483 (|H| = 13, 19 classes): 10⁻⁴·⁵.
- 1057: no joinable structure (16 without D, 24 above 6 · 10⁷ members).
- 1261: only |H| = 3, with 210 pairs.
- 1333: 16 tiny structures, all empty.

The estimates explain the failure. Per ⟨K, −1⟩-class the tri structure has 3 free signs in A and 2 in D. The spin
structure (A K-invariant, B, B∘μ, B∘μ²) has 1 + 6. At these v a class equation costs about 5.3–5.5 bits. So spin
gains about 1.7 bits per class, while the tri structure loses a few tenths and works only for small v:
- with H trivial, 10^5.3 at v = 133 (22 classes) and 10^4.2 at v = 241 (40 classes);
- with |H| = 5, 10^1.9 at v = 211.

## 36. Spin structures: an exact fibre test, and the landscape of the tabu search

**Fibre test (`spin_fibre.py`).** Over a composite ℤ_v, take a shift s in a proper subgroup dℤ_v. Its
autocorrelation splits over the d cosets a + dℤ_v (the fibres): P_X(s) = Σ_a P_{X_a}(s), where
X_a(y) = X(a + dy). The fibre orbits under a group G of multipliers are independent, and the fibre g·a contributes
the term of fibre a at the pulled-back shift g⁻¹s. So a DP over fibre orbits gives exactly:
- the attainable values of P_A(s) for skew K-invariant A (fibres a and −a tied by the skew condition);
- the attainable values of T_B(s) = Σ_j P_B(μʲs) for H-invariant B, together with the row sum of B.

The spin condition P_A(s) + T_B(s) = 0 needs the two sets to meet on every subgroup class. The DP agrees with brute
force (all 256 A of 721) and contains every sampled T_B.

The two composite hexagonal structures of the tabu queue (§33) are **infeasible**:
- 721, H = ⟨8⟩, μ = 2 or 18: on the class of 7, P_A ∈ {−7, 57, 65, 257, 273, 593, 617}, but T_B ≥ 11.
- 817, H = ⟨11⟩, μ = 9 or 6: on the class of 19, P_A ∈ {−19, 53, 77, …}. The value needed, T_B = 19, is not
  attainable (attainable values run −41, −5, 7, 43, …), and T_B ≥ −41 rules out the larger P_A.

This explains why their low-energy tabu states always kept the residual −12 on the class of 19 (817), and why
2 × 40 minutes each found nothing, although the Gaussian model gave them more solutions per unit of search space
than 919. The test passes for 919 (prime: no subgroup classes).

**Landscape of the tabu search** (per-restart best objective, 3000 steps, levels in units of the smallest
nonzero objective):
- 919 (E = 10^5.4, 55 + 9 variables, 9 classes): level 1 in 1.5% of restarts, level 2 in 14%, and a solution
  about every 3000 restarts.
- 1657 (E = 10^5.3, 73 + 12 variables, 12 classes): over 7160 restarts the minimum was level 2 once and level 3
  once.

With the step-down ratios of 919 near zero, a solution at 1657 needs about 3 · 10⁶ restarts, roughly 60 core-hours
at 15 restarts per second. The four other census hexagonal primes (1801, 1951, 2269, 2437; 12–14 classes) should
be similar or harder.

Per-restart step count barely matters on 919: 3000 versus 10 000 steps give 8.7 · 10⁶ versus 7.9 · 10⁶ total steps
per solution (five seeds each).

**Deep neighbourhoods do not help.** Take the tabu states of 919 at levels 1–2 (230 states from 1500 restarts) and
test every B within two swaps against all 512 skew A at once (exact lookup of −T_B). This found no solution. The
checker recovers solutions from known ones perturbed by up to three swaps. So the low-energy states are not near
solutions in swap distance, and the restarts stay.

## 37. Near-Williamson quadruples by exhaustive search: order 191 and the complex orders

**The route.** Cooper–Wallis needs Williamson-type matrices of order w and T-matrices of order t to give H(4tw).
Near-Williamson quadruples (A, B, C, D) — circulant, B, C, D symmetric, ΣPAF = 4wδ — give the Williamson-type
matrices (AR, B, C, D) ([near-williamson](../near-williamson)), and since three blocks are symmetric they are also
Szöllősi quadruples, i.e. complex Hadamard matrices of order 2w. For w = 191, T-matrices of orders 5, 9 and 15
(from base sequences BS(3, 2), BS(5, 4), BS(8, 7)) would settle three entries of FRONTIER.md at once: 3820
(n = 955), 6876 (n = 1719) and 11460 (n = 2865). §6 and williamson-type/README.md rule out T-matrices of order 191 and
the bordered form on ℤ₁₉₀; the plain form on ℤ₁₉₁ is ranked here.

`bnw.py triage 191 --plain`: one admissible structure, H = ⟨39⟩ (order 5), A on the 39 H-orbits, B, C, D on the 20
L-orbits (L = ⟨39, −1⟩), 19 classes, row sums 21 | 11, 11, 9, log₁₀ E = 1.38. Pools: A 2.9 · 10¹⁰, each symmetric
block 9.2 · 10⁴. (The same group with two symmetric blocks holds the complex matrix of order 382 in
[complex-hadamard](../complex-hadamard), at about 10⁷ predicted solutions; the third symmetric block costs about 2¹⁹.)

**Splitting the arbitrary block.** On each L-orbit O ∪ −O the block A is either symmetric (equal values on O and −O)
or antisymmetric (opposite values). With A_s, A_a the two ternary parts (disjoint supports), Â_s is real and Â_a
imaginary, so PAF_A = PAF_{A_s} + PAF_{A_a} and |Â|² = Â_s² + |Â_a|² exactly. For a fixed split E (the symmetric
L-orbits) the parts range independently over about 2^|E| and 2^(19 − |E|) sign vectors, so every member of A comes
out of a join of two short lists, with its PAF and power-spectrum vectors, at one vector addition per member.

**The power-spectrum filter on triples.** Every quadruple satisfies |Â(t)|² + B̂² + Ĉ² + D̂² = 4w at every
character, so every partial sum is at most 4w. For w = 191, 32% of each symmetric list passes, 5.7% of the pairs
(B, C), and only 3.2 · 10⁷ of the 1.3 · 10¹³ triples (B ≤ C, D): acceptance 2.6 · 10⁻⁶. So all feasible triples fit
in a hash table keyed by their PAF vector, and the structure is searched completely by streaming every member of A
against it (`search/nw_exh.c`, tables from `search/nw_prep.py`). The D candidates of a pair are the AND of two
bitsets (D compatible with B, and with C): 73 s for the triples, then 2.5 · 10¹⁰ members of A (those whose
symmetric and antisymmetric parts pass the spectrum bound separately) and 2.1 · 10¹⁰ table probes in 17 minutes on
three threads (with the table slots prefetched in blocks of 16).

Check: on ℤ₆₁ with H = ⟨13⟩, rows 1 | 11, 11, 1, the run finds 320 quadruples (160 up to A ↔ reversed A, which have
the same PAF), against the Gaussian prediction 10^2.83 = 676 for ordered (B, C) and both A, i.e. 169; each was
checked against direct autocorrelations, and the T-matrices of orders 5 and 9 turn one of them into Hadamard
matrices of orders 1220 and 2196 (`search/nw_certify.py`, construction `cooper-wallis-near-williamson` of verify.py).

**Results (all empty).**

| w | H | rows | log₁₀ E | feasible triples | members of A | hits |
|---:|---|---|---:|---:|---:|---:|
| 191 | ⟨39⟩, order 5 | 21 \| 11, 11, 9 | 1.38 | 3.2 · 10⁷ | 2.5 · 10¹⁰ | 0 |
| 183 | ⟨34⟩, order 5 | 1 \| 19, 19, 3 | −1.12 | 1.5 · 10⁶ | 1.3 · 10¹⁰ | 0 |
| 183 | ⟨34⟩, order 5 | 1 \| 17, 21, 1 | −1.14 | 4.4 · 10⁶ | 1.3 · 10¹⁰ | 0 |
| 183 | ⟨34⟩, order 5 | 9 \| 17, 19, 1 | −1.15 | 3.0 · 10⁶ | 1.5 · 10¹⁰ | 0 |
| 213 | ⟨37⟩, order 7 | 15 \| 13, 13, 17 | −2.81 | 2.4 · 10⁴ | 2.0 · 10⁸ | 0 |

So no near-Williamson quadruple of order 191 is constant on the orbits of ⟨39⟩. The other nontrivial subgroups of
ℤ₁₉₁* fail the row-sum equation (symmetric rows ±1 + 38ℤ or ±1 + 190ℤ), so the plain near-Williamson route to 3820,
6876 and 11460 needs a quadruple without multiplier symmetry. For 183 (complex order 366) and 213 (426) the best
near-Williamson structures are empty too; both orders are constructed from other structures (§41 and §42).
Single-precision sums of spectra are compared with 4w + 1, so rounding cannot reject a quadruple; the 191 run used
this margin (a run with 4w + 10⁻⁶ was also empty).

**Where the route stands.** The same triage over the other FRONTIER composites n = t·w with known T-matrices of
order t: w = 463 (|H| = 11, rows 23 | 21, 21, 21, log₁₀ E = −2.3, would give n = 2315), 631 (|H| = 15,
1 | 29, 29, 29, −3.5, n = 1893), 991 (−31); the other w (167, 179, 223, 227, 239, 263, 269, 283, 347, 359, 383, 419,
431, 443, 479, 491, 503, 563, 571, 599, 643, 863, 887, 911, 941, 947, 971) have no multiplier structure with few
classes or fail the row-sum equation. Over the 214 open complex orders the best near-Williamson structures are the
ones above (183, 213) and 463 (926); the rest are below 10⁻⁴. The Gaussian model is accurate here (ℤ₆₁), so these
structures were long shots; the method itself (lists of 10⁴–10⁵ symmetric blocks, an arbitrary block of 10¹⁰–10¹¹)
is cheap, and any near-Williamson structure with a predicted 10 or more solutions can be settled in hours.

The search for near-Williamson matrices of Kharaghani, Mohammadian, Tayfeh-Rezaie and Zaitsev (arXiv:2605.08661v2,
§4) works without symmetry: it enumerates A and B and then solves for the symmetric pair (C, D) by linear algebra
modulo 4. For a structure like ours, where the arbitrary block dominates, streaming A against all feasible triples is
cheaper than enumerating pairs (A, B).

## 38. Anti-automorphisms of the Theorem 7 designs: no symmetric 7956 this way

**The idea.** The skew matrix of order 7956 = 4 · 1989 ([skew/7956](../skew/7956)) is H = I ⊗ (I + W) + B ⊗ J + C ⊗ D,
with aI + bB + cC a skew-type OD(156; 1, 3, 152) in Goethals–Seidel form over ℤ₃₉ and (I + W, D) the amicable cores
of order 51. 1989 is also undecided in the symmetric census. If a symmetric signed permutation P₁ of order 156 (P₁² =
I) satisfies P₁(aI + bB + cC)P₁ = aI − bB − cC, then P₁, P₁B and P₁C are symmetric (P₁B = −BP₁ gives (P₁B)ᵀ = −BP₁ =
P₁B) with disjoint supports, they anticommute pairwise, and X = P₁B (weight 3) with Y = P₁ + P₁C (weight 153) is a
symmetric OD(156; 3, 153). The weighted plug-in of §32 then gives the symmetric Hadamard matrix X ⊗ J₅₁ + Y ⊗ D′ of
order 7956 for D′ the core of any symmetric Hadamard matrix of order 52 (the J terms cancel because 3 · 52 = 156).
The Balonin form of §32 is ruled out for 1989 by the row-sum screen, but P₁(aI + bB + cC) need not be in that form.
The tensor version P = P₁ ⊗ Q on H itself would also need Q W Q = −W and Q D Q = ±D on the cores.

**Results (`search/od_anti.py`).**
- *The amicable pair of order 52.* No signed permutation σ satisfies σSσᵀ = Sᵀ and σMσᵀ = ±M for the pair (S, M)
  of the certificate (backtracking over all 52 points: none). The joint automorphisms (σSσᵀ = S, σMσᵀ = M) form a
  group of order 3 (up to sign). So the tensor version fails for this pair; the plug-in version needs only P₁.
- *Block-monomial P₁.* Let P₁ map block i to block π(i) by ε_i times a signed permutation g ↦ c_i + k_i g of ℤ₃₉. The
  diagonal blocks (all equal to the skew-type first block) need k_i ∈ −Stab(A), the off-diagonal blocks force all k_i
  equal, and P₁² = I needs k_i² = 1. When Stab(A) has odd order (here {1, 16, 22}) the only such k is −1: reflections,
  as in `gs_anti.py`. The b part (one entry at 0 in each of blocks 2–4, so ±R in every off-diagonal block of B) forces
  c_i + c_j to be the same for all i ≠ j, hence one common centre, and P₁ is block-monomial in ±R: an involution π and
  signs ε, 80 candidates. For every sign pattern of the b entries, the signed union-find of P₁GP₁ = Gᵀ (G the design
  in the normal form of `od13.c`) leaves 48 consistent configurations (computed for t = 5, 7 and 39; the form does not
  depend on t). In every one of them one of W₂, W₃, W₄ is symmetric and the other two are tied, W_l(g) = ±W_k(±g). The
  bc condition of §31 asks that the symmetric parts of W₂, W₃, W₄ add to 0 on every pair {s, −s}; the symmetric block
  contributes ±1 and the tied pair S_k + S_l with S_l = ±S_k, which is 0 or ±2. So **no Goethals–Seidel OD(4t; 1, 3,
  4t − 4) has a block-monomial involutive anti-automorphism** when the stabiliser of the first block has odd order,
  for any t. For the design of the certificate, the direct test over all reflection configurations with arbitrary
  centres also finds none (`od_anti.py cert`).

The same question for the other Theorem 7 designs of the repository, whose orders are symmetric-census values too:
- [skew/1900](../skew/1900) (OD(100; 1, 5, 94), b parts δ₀ + δ_d, δ₀ − δ_d, δ₀; plug-in order 5 · 20 · 19 = 4 · 475):
  the certificate's design has no block-monomial reflection P₁, and for every d, every placement of the singleton
  and every sign pattern of the b entries (t = 25) no P₁ is consistent with the a, b and c parts at all;
- [skew/2772](../skew/2772) (OD(252; 1, 21, 230), plug-in order 21 · 12 · 11 = 4 · 693): the certificate's design
  has none (the only k with k² = 1 and −k ∈ Stab(A₁) is −1 for all five designs, so reflections are the whole
  block-monomial family here).

So symmetric matrices of orders 7956 and 1900 along these lines need a P₁ that mixes blocks within a row of blocks,
a design outside Goethals–Seidel form, or a symmetric OD built some other way.

**The skew matrices themselves.** A skew-Hadamard matrix S of order 4v with an involutive signed anti-automorphism P
gives the symmetric Hadamard matrix SP of the same order ((SP)ᵀ = PSᵀ = SP), and 24 of the skew orders of this
repository have v in the symmetric census. For the 29 cyclic Goethals–Seidel certificates, `od_anti.py gs` tests
every block-monomial P with affine maps r ↦ d_i + kr (k² = 1, every d_i; each block of S is circulant or
back-circulant, so each pair condition depends only on d_j − d_i or d_i + d_j). None has one. (The test finds the
anti-automorphisms of the n = 5 and n = 7 solutions of `gs_anti.py`.) This is what §31 predicts: such a P either
leaves (n + 3)/2 free classes or forces the skew-propus layout, a special propus quadruple, which the symmetric
sweeps search directly.

## 39. Spin structures over ℤ₇² × ℤ₄₃: the order 8428 = 4 · 2107 is out of reach this way

2107 = 7² · 43 is in the skew census, and 4 · 2107 = 1 + 3 · 53², so spin structures (row sums 1, 53, 53, 53) are
admissible. Over the cyclic group they are all below 10⁻¹ (§14). Over ℤ₇² × ℤ₄₃, with automorphisms
(x, y) ↦ (Mx, uy) from GL(2, 7) × ℤ₄₃*, `spin_group.py rank 7 43` lists 401 pairs (H, μ) and gives 38 a finite
estimate, up to 10^18.5 (7 K-classes). The exact fibre test of §36, generalised to this group (`spin_group.py check
7 43 H μ`), decides every class whose K-orbit generates a proper subgroup U from the attainable values of P_A and of
T_B on the cosets of U.

The 19 structures with the highest estimates (all with log₁₀ E > −0.5, the other 19 lower) are **infeasible**, each on
the class of the subgroup {0} × ℤ₄₃ (|U| = 43): the skew block A only reaches values of P_A that −T_B cannot match.
Typically also the class of a subgroup of order 7 is dead. Example, H = ⟨(0,1;6,6), 4⟩, μ = ((0,1;5,3), 6), estimate
10^18.5: on the class of U = ℤ₄₃, P_A ∈ {−49, 191, 215, 1631, 1799, 1871, 1895, 2039} while T_B ≥ −111 never meets
−P_A; on the class of order 7, P_A ∈ {1763, 1931, 2099} against T_B ≥ −891.

| H | μ | K-classes | log₁₀ E | dead classes |
|---|---|---:|---:|---|
| (0,1;6,6), 4 | (0,1;5,3), 6 | 7 | 18.52 | U = ℤ₄₃, \|U\| = 7 |
| (0,1;3,5), 10 | (0,2;3,6), 1 | 7 | 18.32 | the same |
| (0,2;5,5), 10 | (0,1;3,5), 6 | 7 | 18.32 | the same |
| four more with 7 classes | | 7 | 15.9–17.7 | the same |
| six with 11 classes | | 11 | 12.0–12.5 | the same |
| (0,2;3,6), 9 | (0,1;6,6), 6 | 15 | 4.02 | the same |
| (0,1;6,6), 10 | (0,2;3,6), 9 | 15 | 4.00 | the same |
| (0,1;3,5), 10 | (0,4;3,3), 1 | 13 | 2.76 | the same |
| (0,2;3,2), 15 | (0,1;3,5), 1 | 13 | 1.04 | the same |
| two with 19 classes | | 19 | −0.4 | U = ℤ₄₃ |

(M written a,b;c,d for the matrix [[a, b], [c, d]] over GF(7), then the unit of ℤ₄₃.) As for 721 and 817 (§36), the
Gaussian model ignores the subgroup classes, where the values of the blocks live on a coarse lattice. So the spin
route to 4 · 2107 needs a structure without these subgroup classes, i.e. the cyclic group, where every structure is
below 10⁻¹.

## 40. What the rankings say about the remaining symmetric and FRONTIER leads

**The symmetric chain sweep beyond v ≈ 440.** The compressed propus sweep of §28 had covered the composite census
values up to 435 (a coverage map rebuilt from the chain logs) before it was stopped. It was restarted on the 85
composite values from 445 to 995 (orbits 1–36); its outcome is not recorded. The Gaussian estimate (§29) puts every
ranked structure with at most 45 orbits for v ≥ 440 below 10^−1.5 (best: 589/⟨28⟩ 10^−1.5, 451/⟨5⟩ 10^−2.3, 445/⟨16⟩
10^−3.1), so the sweep mostly produces exact negatives. The estimates were low by 10^1.4–10^1.9 for the two hits of
this family (532, 876), so a hit is not excluded.

**163 (order 652).** The estimates disagree with the exact count where the count is available. With H = ⟨38⟩ (order
9; A on 10 orbits, B and D on 19, 9 classes) the estimate is 10^1.47 for the row patterns (17, 1, 19) and
(19, 1, 17), but the exhaustive joins (pools 126, 48 620, 43 758) are empty. This is not a congruence obstruction:
the class values of A + D and of 2B reach compatible residue vectors modulo 2, 3, 4 and 8. With H = ⟨58⟩ (order 3)
the estimate is 10^4.2 (rows 19, 11, 7), but the pools are 10^7.2, 10^15.2 and 10^15.3, beyond any join, and 163 is
prime, so there is no compression.

**FRONTIER 1169 by Spence spin structures** (NOTES.md §13 lists it as an untested survivor). The structure is
G = GF(16) × ℤ₇₃ with μ(x, y) = (σx, hy), σ a primitive element of GF(16)* and h of order 9 in ℤ₇₃*,
K = ⟨μ⟩ of order 45 and H = ⟨μ³⟩ of order 15 (`ss_core.Spin`): 34 K-orbits, 100 H-orbits, 17 classes,
log₁₀ E = 4.94 (diagonal) and 6.78 (full covariance), pools 10^8.2 (A) and 10^27.9 (B). Its quotients are only ℤ₇₃
(compressing GF(16); 4 classes) and GF(16) (1 class): σ is irreducible on GF(16), so there is no intermediate
subgroup to play the part of ℤ₂² × ℤ₁₃ for 573. Over a ℤ₇₃-level solution each fibre GF(16) × {y} splits into
{0} and the three cosets of the subgroup of order 5, so a compressed value ±16 or ±14 fixes the fibre and ±6 or ±4
leaves the choice of the coset with the minority sign: a 3-colouring of about 33 fibre orbits, with the 12
remaining class equations quadratic in it. Nothing like the rigid ℤ₁₁ fibre of 573 exists here, so the
level-and-lift method does not transfer; a direct local search (joint space 10^36, about 10^−30 solutions per point)
is far beyond what the tabu searches of §36 reach.

## 41. Paley products for the complex orders: order 330

**Structure.** For N = m q with q a prime or prime power, gcd(m, q) = 1, and H = {h : h ≡ 1 (mod m), h a square in
GF(q)}, every block of Szöllősi's array (A, B symmetric, C, D arbitrary) is X = u ⊗ δ + w ⊗ J + z ⊗ χ as in §15, and
the condition splits on δ, J and χ (ε = χ(−1), (f ⋆ g)(s) = Σ_a f(a + s) g(a)):

    (I)   Σ_k PAF_{u_k} + q PAF_{z_k}                 = 4N δ
    (II)  Σ_k (u_k ⋆ w_k + w_k ⋆ u_k) + q PAF_{w_k}   = Σ_k PAF_{z_k}
    (III) Σ_k ε (u_k ⋆ z_k) + (z_k ⋆ u_k)              = 0.

The symmetric blocks have u, w even and z(−a) = ε z(a), so for q ≡ 3 (mod 4) position 0 of A and of B is a
w-position. The compressions y_k = u_k + q w_k to ℤ_m satisfy Σ PAF_y = 4N δ, which is (I) + q (II); given it, (I)
and (II) are equivalent, and its value at 0 is the counting identity (q + 1) N_w − 4 N₋ = 4m of §15.

**Solver** (`search/ssxx_paley.py`, exact). The skeletons (w-positions and their values ±q, ±(q − 2)) are listed up
to multipliers of ℤ_m, translations and reversal of C and D, block signs and the swaps A ↔ B, C ↔ D. Stage 1 joins,
per skeleton, the ±1 entries of y on Σ PAF_y = 0 (join4x); its solutions are reduced to one per orbit of the same
group; stage 2 joins the z_k on (II) and (III). Every hit is lifted to ℤ_N and checked as a complex Hadamard matrix.
- Validation against the exact joins of `ssxx.py` on the same structures, without symmetry reduction: N = 77 = 7 · 11
  gives 451,584 = 16 · (2 · 2352 + 4 · 1176 + 4 · 2352 + 4 · 2352) solutions (its four row patterns, times block signs
  and orderings), and N = 65 = 5 · 13 (q ≡ 1 mod 4) gives 38,400 = 16 · 2400. With the reduction, N = 77 has the same
  6 stage-1 orbits, of which the same 4 extend.

**Order 330.** N = 165 = 15 · 11, H = ⟨136⟩ of order 5. The counting identity allows (N_w, N₋) = (5, 0), (6, 3) and
(7, 6). N_w = 5 (34 skeletons) has no solution. Among the 2455 skeletons of N_w = 6 the run found three verified
solutions when it ran out of memory (it kept every stage-1 orbit as a tuple; they are now packed bytes). No exhaustive
count is recorded. This gives the complex Hadamard matrix of order 330
([C330.json](../complex-hadamard/C330.json), rows 7, 9; 1, 23). The order is blank in Seberry–Yamada's Table 11.2,
329 = 7 · 47 is not a prime power, and no quaternary Legendre pair of length 164 is known (neither family of
Jedwab–Pender covers it).

**Empty structures** (every admissible (N_w, N₋), all skeletons, all stage-1 orbits):

| N | m × q | (N_w, N₋) | skeletons | stage-1 solutions (orbits) | stage 2 |
|---:|---|---|---:|---|---:|
| 171 | 9 × 19 | (2, 1) | 1 | 486 (1) | 0 |
| 253 | 11 × 23 | (2, 1) | 1 | 0 | – |
| 275 | 11 × 25 | (2, 2) | 6 | 1,120 (14) | 0 |
| 323 | 17 × 19 | (4, 3) | 10 | 480 (14) | 0 |
| 345 | 15 × 23 | (3, 3) | 1 | 240 (1) | 0 |
| 391 | 17 × 23 | (3, 1) | 2 | 0 | – |
| 513 | 19 × 27 | (3, 2) | 2 | 1,368 (1) | 0 |
| 629 | 17 × 37 | (2, 2) | 6 | 24,832 (194) | 0 |

So the orders 342, 506, 550, 646, 690, 782, 1026 and 1258 have no Szöllősi matrix invariant under the Paley-product
group (nor under any larger multiplier group containing it). The counting identity alone empties 207 = 9 · 23 (no
admissible N_w ≥ 2). Small q relative to m (q = 11 for m = 15) means more w-positions and many more skeletons; that is
where the solutions were. The q ≈ 2m cases, with two or three w-positions, have few stage-1 solutions and none
extends.

**With an extra multiplier group K on ℤ_m** (`search/ssxx_paleyk.py`: every block constant on the K-orbits, skeletons
and both joins exact over them). Empty, exhaustively: N = 171 (q = 19, K = ⟨4⟩), 209 (q = 11, K = ⟨4⟩ and ⟨7⟩; the
latter 1380 skeletons, 120 stage-1 solutions, none extends), 253 (q = 11, ⟨2⟩; q = 23, ⟨3⟩), 299 (q = 13, ⟨2⟩), 319
(q = 11, ⟨7⟩), 341 (q = 11, ⟨2⟩: 1944 skeletons, 48 stage-1 solutions; ⟨7⟩), 391 (q = 17, ⟨2⟩), 403 (q = 13, ⟨2⟩:
2124 skeletons; ⟨7⟩), 407 (q = 11, ⟨7⟩) and 437 (q = 19, ⟨2⟩). Larger cases (341/⟨5⟩, 357/⟨4⟩, …) exceeded ten
minutes each and the sweep was stopped.

**The bordered array** (`--bordered`). On ℤ_m × GF(q) of order N − 1 (m even, bsz.py's array: Σ PAF = −4 off 0)
the same splitting gives Σ PAF_y = −4q off 0, Σ PAF_z = R + 4 and (III); the counting identity becomes
(q + 1) N_w − 4 N₋ = 4(m − 1), and for ε = −1 both fixed points 0 and m/2 of a ↦ −a are w-positions of A and B.
Every one of the 50,688 bsz solutions for N = 31 (ℤ₆ × GF(5), counted by bsz's exact joins in all four row
patterns) satisfies the bordered equations, and the solver's hits pass bsz.build. Empty, exhaustively: N = 171
(ℤ₁₀ × GF(17)) and 407 (ℤ₁₄ × GF(29)), 16 skeletons each.
- **N = 599, the order 2396 of FRONTIER.md through a complex matrix of order 1198** (ℤ₂₆ × GF(23)): the counting
  identity allows only N_w = N₋ = 5, so A and B carry ±21 at 0 and 13 and the fifth w-position lies in C or D. On the
  compression x to ℤ₁₃ the shifts 0 and 13 give Σ x² = Σ y² + Σ PAF_y(13) = 2304 − 92 = 2212. If A or B has equal
  signs at 0 and 13 its x(0) = ±42 costs 1764, so at most one of them does; the exact ℤ₁₃ join is empty for every
  sign pattern in which one of them does not, and with both equal the cost 2 · 1764 + 400 exceeds 2212. **The bordered
  Paley structure of 599 is empty.**

**Order 426.** The split values of `lowe_sz.py` are ordered for kt4's pairing (A, C), (B, D) (§50, "Split values").
The first structure ranked by the mid-size triage, N = 213 = 3 · 71 with H = ⟨37⟩ of order 7 (33 orbits for C and D,
16 classes; rows 13, 17; 13, 15; estimate 10^2.2 in the structure and 10^0.95 inside lists of 6435, 6435, 3 · 10⁷ and
3 · 10⁷ kept from pools of 4.4 · 10⁸), gave 72 solutions within the first split values (six minutes in all), all
checked. This is the complex Hadamard matrix of order 426 ([C426.json](../complex-hadamard/C426.json)); 425 = 5² · 17
is not a prime power and no quaternary Legendre pair of length 212 is known. The 183 ⟨34⟩ join (rows 17, 19; 1, 9)
found no solution in the first 6833 split values, 17% of the pair-product mass (§51).

## 42. Compressions decide which row patterns live: orders 366 and 454

**The defect.** `ssxx.gauss` (and `bsz.gauss`, and the low-energy runs built on them) treats the class sums of a
structure as independent continuous variables. For a divisor d of N the fibre sums x(a) = Σ_{g ≡ a (mod d)} X(g) of
every block obey, exactly,

    Σ_X PAF_{x_X}(s) = 4N δ_s          on ℤ_d               (bordered: 4m + t(m/d − 1) at 0, t m/d elsewhere, t = −4)

and when H fixes ℤ_d pointwise every block orbit lies in one fibre, so x(a) = X(fixed point) + |H| · (sum of orbit
signs) ≡ ±1 (mod 2|H|). The compressed equation then has very few solutions, and a row pattern either has none (no
solution of the full problem) or concentrates all its solutions on a few compression tuples.

**Calibration** (exact joins of `ssxx.py`, N = 93 = 3 · 31, H = ⟨64⟩ of order 5, trivial on ℤ₃):

| rows | solutions | plain Gaussian | ℤ₃ tuples | conditioned |
|---|---:|---:|---:|---:|
| 7, 9; 11, 11 | 5,616 | 10^2.54 | 9 | 10^3.94 |
| 7, 11; 9, 11 | 9,504 | 10^2.53 | 45 | 10^4.62 |
| 9, 11; 1, 13 | 0 | 10^2.08 | 0 | 0 |
| 9, 11; 7, 11 | 2,808 | 10^2.05 | 9 | 10^3.95 |
| 9, 13; 1, 11 | 21,384 | 10^2.56 | 18 | 10^4.42 |
| 11, 11; 3, 11 | 0 | 10^2.09 | 0 | 0 |
| 11, 11; 7, 9 | 1,296 | 10^2.05 | 9 | 10^3.94 |
| 11, 13; 1, 9 | 0 | 10^2.57 | 0 | 0 |

The plain estimate is 12–60 times low on the live patterns and blind to the dead ones; the conditioned estimate
(`search/comp_est.py`: every achievable compression per block with its exact count, the feasible tuples, then per
tuple a multivariate lattice Gaussian of the class sums from blocks sampled with those compressions — exact Hermite
form and Gram determinant for the lattice, pseudo-inverse covariance on its span) is 1.2–7 times high and exact on the
zeros. Bordered (N = 83, ℤ₂ × ℤ₄₁, |H| = 5, row 2 in A): 17,536 solutions, plain 10^0.83, conditioned 10^5.21.
Compressions to quotients on which H acts non-trivially (ℤ₃₁ for N = 93) keep every live pattern alive (no false
zeros) but do not see the dead ones, so the triage takes the smallest conditioned estimate over the quotients with
at most four fibre classes. The order-426 matrix has ℤ₃ compressions (15, −1, −1), (15, 1, 1), (15, −1, −1),
(1, 13, 1), one of the predicted tuples (all ≡ ±1 mod 14).

**Consequences.**
- **171 (order 342), rows 15, 15; 3, 15 under ⟨68⟩ or ⟨11⟩ is empty.** On ℤ₉ (fibre classes {0}, {3, 6}, ℤ₉*) the
  compressions are (x₀, x₃, x₁) with x₀, x₃ ≡ ±1 (mod 6) (±(5, 7, 17, 19) for A and B) and x₁ odd, and the shifts 0
  and 3 give Σ_X (x₀ − x₃)² = 4N = 684. The row sums force x₀ ≡ x₃ (mod 6) in every block, so y = (x₀ − x₃)/2 is a
  multiple of 3, Σ (y/3)² = 19, and y_A/3, y_B/3 ∈ {0, 2, 4} leave 19, 15, 11 or 3 for two squares. The two low-energy
  runs of this pattern (logE 2.0 and 2.1 by the plain model) could not succeed; they were stopped.
- **183 (order 366), rows 1, 23; 9, 11 under ⟨34⟩** lives on a single compression class (×9 by translations and
  reversals of C and D): A ≡ (−1, 1, 1), B ≡ (1, 11, 11), C ≡ (−1, −1, 11), D ≡ (−11, 11, 11); conditioned logE 6.1
  against 5.0. The unconditioned run kept the 2 · 10⁷ blocks of lowest energy out of 2.6 · 10¹⁰, centred on typical
  compressions, and found nothing in 5,900 split values. Conditioned: pools 18,480, 15,840, 6.76 · 10⁸ and
  4.97 · 10⁸ (lowe_enum leaf counts equal to the compression DP), lists of 1.9 · 10⁷ and 2 · 10⁷ (retention 0.055 and
  0.086, 10^2.8 predicted inside, 20 minutes to build), and kt4 found 8 solutions in 22 minutes. The first is the
  complex Hadamard matrix of order 366 ([C366.json](../complex-hadamard/C366.json)): compressions exactly the class
  above, blank in the 1992 table, 365 = 5 · 73 not a prime power, 183 not among the near-Williamson orders of 2026.
- **227 (order 454), bordered, ⟨49⟩ of order 7 on ℤ₂₂₆ = ℤ₂ × ℤ₁₁₃.** 227 = 2 · 113 + 1 is a safe prime, so only the
  bordered array has a multiplier group. The fibre sums over ℤ₂ are ≡ ±1 (mod 14) and the compressed equation is
  Σ_X (x₀ − x₁)² = 4N = 908; for rows 0, 0; 2, 0 it leaves 8 tuples (4 up to symmetry), all with one block at
  (∓15, ±15). Plain estimate 10^0.6, conditioned 10^6.0. The class ((1, −1), (1, −1), (1, 1), (−15, 15)) has pools
  4900, 4900, 1.66 · 10⁸ and 1.31 · 10⁸ (lowe_enum's leaf counts agree with the compression DP); lists of 2 · 10⁷
  (retention 0.22 and 0.24 by the conditioned model, 10^3.8 solutions predicted inside) were built in 5 minutes,
  and kt4 found 64 solutions in the first split values (6 minutes). Five were checked by bsz.build, and the
  first is the complex Hadamard matrix of order 454 ([C454.json](../complex-hadamard/C454.json), border
  (1, 1, 1, 1, −i, −i)): ΣPAF = −4 off 0, A and B symmetric, all four blocks ⟨49⟩-invariant, compressions as
  predicted. Blank in the 1992 table; 453 = 3 · 151 is not a prime power; no quaternary Legendre pair of length
  226 is known to us. The first matrix from the bordered array of an open order.

**The 2q-bordered family is decided by its ℤ₂ compressions.** For open N < 1000 with N − 1 = 2q, q prime, and every
odd-order H ≤ ℤ_q* with 18–50 orbits on ℤ₂ × ℤ_q, the fibre sums over ℤ₂ are ≡ ±1 (mod 2|H|) and Σ_X (x₀ − x₁)² = 4N.
Only N = 203 (|H| = 5, conditioned logE 8.0) and N = 227 (|H| = 7, 6.2) have a feasible row pattern. For N = 179
(|H| = 11), 219 (9), 263 (13), 275 (17), 303 (15), 383 (19), 467 (29), 479 (17), 483 (15), 503 (25), 563 (35),
623 (31), 699 (29), 795 (33), 803 (25), 839 (19), 863 (43), 915 (19, 57), 923 (23), 927 (21, 33) and 983 (35, 49)
no pattern is feasible, so these bordered structures are empty (exact; `python3 search/comp_scan.py bordered2q`).
Likewise N = 223 (order 446) on ℤ₂₂₂ = ℤ₆ × ℤ₃₇ with H = ⟨127⟩ of order 9 (H trivial on ℤ₆): the ℤ₂ compressions rule
out all four row patterns, although the ℤ₃ compressions alone allow 36–168 tuples per pattern.

**Propus (symmetric Hadamard): conditioning is not enough.** The same estimator with weights (A + 2B + D) matches
v = 57 = 3 · 19, ⟨7⟩, rows 3, 5; 13 with free D (864 solutions, conditioned 10^3.07, three patterns correctly
dead), but it predicts 10^1.8 for two dead patterns of v = 93, ⟨64⟩ (free D) and 10^2.25 for 65/⟨16⟩, rows 11, 3,
11, which §30 proves empty (no solution modulo 16). The doubled block enters the equations only from modulo 16 on,
and those 2-adic obstructions are invisible to compressions to ℤ₃ or ℤ₃₁ and to lattice cosets (the mod-8 bit
distributions of A and D even favour solutions at v = 93: 2.5–4 times the uniform rate). So for propus the scan
(`comp_scan.py propus`) ranks structures but a lead needs an exact or modular check before a search.

**Tools.** `lowe_enum.c` takes fibre-class constraints (per class, the weighted sum of the orbit signs over it is
fixed; exact reachability pruning per class) next to the fibre filters of §55; `search/lowe_comp.py` runs the
low-energy k-tree on one compression class at a time (classes reduced by translations by H-fixed points, reversal of
C and D, and multipliers); `search/comp_triage.py` triages every open N (plain and `--bordered`).

Also fixed: `kt4.c` leaked its hash table on every split value with fewer than 2¹⁸ pairs (the cause of the
171 crashes); allocations are now checked, and `lowe_sz.py --resume --skip k` (or `--auto`) continues a join.

## 43. Conditioned scans of every open complex order: order 430; twisted and non-cyclic conditioning

**Plain product structures.** `comp_scan.py plain` over all open N < 1000 (N = d · q, d = 3, …, 15, q prime, H ≤ ℤ_q*
of odd order fixing ℤ_d, 15–60 orbits per arbitrary block): 100 structures.
- 68 have a listable compression space. 62 of them have no row pattern with a feasible ℤ_d compression tuple, so they
  are empty. 171 (|H| = 9), 207 (11) and 921 (17) have feasible tuples but no tuple with an estimate. Three are live:
  215 (|H| = 7, conditioned 10^5.03), 453 (|H| = 15, 10^2.12) and 633 (|H| = 15, 10^2.98).
- The other 32 were skipped (the old `space()` bound ∏(orbits + 1) also undercounted: a class whose orbits have
  unequal weights takes up to 2^k values; it now counts the values). A necessary test, the reachable pairs
  (Σ_a x(a)², Σ_a x(a)) per block by a DP over the fibre classes, against Σ_X Σ_a x_X(a)² = 4N and the row pattern,
  empties 24 of them: 209, 299, 341, 403, 583, 603, 611, 671, 707, 711, 721, 737, 755, 767, 781, 869, 909, 913, 923,
  927, 955, 963, 973, 981 (for several, no row pattern exists at all). Left: 253, 319, 333, 345, 407, 553, 655, 763.
  On the coarser quotient ℤ₃, 333 (|H| = 9, 45 orbits) has conditioned estimates 10^0.4–10^1.1 per pattern, spread
  over 18–108 tuples with pools near 10⁷ each.

**Order 430.** N = 215 = 5 · 43, H = ⟨41⟩ of order 7 (trivial on ℤ₅), rows 15, 15; 11, 17: 450 feasible ℤ₅ tuples,
3 classes up to translations, reversals and multipliers, each with about 10^2.3 predicted solutions and pools of
1200, 1200, 2.4 · 10⁶ and 2.4 · 10⁶. `lowe_comp.py 215 41 --rows 15,15,11,17 --d 5 --only 0` lists the first class
completely (seconds) and kt4 found 5 solutions in its first split values (514 s in all, most of it the estimate).
All five rebuild as complex Hadamard matrices; the first is [C430.json](../complex-hadamard/C430.json). The plain
Gaussian of this pattern is 10^−0.8: the compressions put every solution into a few classes of tiny pools.

**Twisted structures** (`search/tsz_comp.py`). The full and half twists (B = A∘μ, D = C∘μ; A, B on the L′-orbits) have
compressions x_D = x_C ∘ μ̄, so C contributes P_x(s) + P_x(μ̄ s) to the compressed equation. Per block: exact
compressions (on the L- or L′-orbits for the hashed side), the feasible groups, and the lattice Gaussian of the
twisted class values. Calibration against exhaustive counts (`tsz.py test`):

| N | twist, H | plain Gaussian | conditioned | solutions |
|---:|---|---:|---:|---:|
| 39 | half, ⟨25⟩, rows 5, 9, 5 | 10^1.01 | 10^2.34 | 144 |
| 51 | half, ⟨16⟩, rows 9, 11, 1 | 10^1.41 | 10^2.44 | 288 |
| 15 | half, H = 1, rows 3, 7, 1 | 10^2.14 | 10^2.24 | 480 |
| 21 | half, H = 1, rows 3, 5, 5 | 10^2.61 | 10^2.82 | 1008 |
| 27 | bordered half, H = 1 | 10^4.42 | 10^5.80 | 406,224 |
| 31 | bordered half, H = 1 | 10^3.83 | 10^4.19 | 61,440 |
| 65 | full, ⟨16⟩, rows 11, 3 | 10^3.30 | 10^3.86 | 12,560 |

The exact search is restricted to the feasible compressions: the hashed side is the explicit list of the A (and B)
whose compressions lie in a feasible group, and `tsz.c` now takes per-class fibre constraints for the streamed C (one
C compression per translation class by H-fixed points), so the search is exhaustive for the pattern at a fraction of
the cost. The bordered half twist is often killed on ℤ₂: with r_A = 0, r_B = 2, r_C = r_D = 0 the compressed equation
reads x₀(A)² + (x₀(B) − 1)² + 2x₀(C)² = N, and when L′ acts freely on ℤ_q* with an odd number n of orbits of size k,
|x₀(A)| ≥ k − 1. At N = 203 (k = 20, n = 5) this exceeds √203, while the plain twisted Gaussian rated that structure
10^0.24. N = 219 (order 438), half twist H = ⟨82⟩ of order 6, μ = 46: conditioned 10^0.18 and 10^0.16 for rows
3, 25, 11 and 3, 23, 13 (plain 10^−2.5), one feasible group each, C streams of 6.8 · 10⁸ per translation class
against a full pool of 2.6 · 10¹⁰. Rows 3, 25, 11: no solution (exhaustive).

**Non-cyclic bordered groups** (`search/gcomp.py`): the same conditioning for bsz_group's spin structures on
ℤ₂² × ℤ_m (N = 4m + 1), with Q = ℤ₂² and fibre classes the orbit images in Q. Calibration: N = 21 (ℤ₂² × ℤ₅, trivial
symmetry), rows 2, 0, 0, 0: conditioned 10^7.84 against 4.6 · 10⁷ solutions (plain 10^7.02); N = 13: 10^4.18 per
pattern against 12,672. It does not see the obstruction that empties N = 17 over ℤ₂² × ℤ₄ (conditioned 10^6.4), which
lives on a finer quotient.

Also: `tsz.write_input` streams the (A, B) pairs in chunks (a 3 · 10⁷-pair list took 5 GB).

## 44. Order 342; finer quotients, live split values, and the estimates against the searches

**A coarse quotient can hide an obstruction.** `comp_scan.py plainall` (every H ≤ ℤ_N* with −1 ∉ H, 15–60 orbits)
rated N = 171 = 9 · 19 with H = ⟨31⟩ of order 6 (order 3 on ℤ₉, trivial on ℤ₃) at 10^3.4 for rows 9, 23; 5, 7 on
ℤ₃, with a best tuple of 10^2.0 and pools 1512, 1512, 6.0 · 10⁶, 6.0 · 10⁶. The tuple is empty: sampled from the
lists, the class sums of classes 7, 10 and 15 never vanish together (0 joint hits where independence predicts 73 in
1.5 · 10⁶ samples), although every pair of them does (about 5,300 times in 4 · 10⁶). On ℤ₉, which H does not fix
pointwise, the pattern keeps 252 feasible tuples (7 classes), conditioned 10^2.8 in all and about 10^0.5 per tuple.
The scans now use every listable quotient ℤ_d (d | N, d < 64), not only those H fixes pointwise, and every
multiplier group (even orders without −1 included; the bordered scan had skipped them).

**Order 342.** With all quotients, the strongest structure of N = 171 is H = ⟨11⟩ of order 6 (on ℤ₉ it is the whole
unit group, on ℤ₁₉ it has order 3; its element of order 2 is (−1 mod 9, 1 mod 19), so −1 ∉ H), rows 11, 11; 9, 19:
conditioned on ℤ₉ 10^4.37 (13 feasible tuples, all inequivalent), best tuple 10^3.60 with pools 1134, 378, 2.77 · 10⁷,
2.77 · 10⁷. `lowe_comp.py 171 11 --rows 11,11,9,19 --d 9 --only 0` kept 2 · 10⁷ of each large block (retention 0.93
and 0.92, 10^3.53 predicted inside), found 49,652 live split values exactly, and kt4 returned 12 solutions in 263 s;
all five checked rebuild as complex Hadamard matrices of order 342 ([C342.json](../complex-hadamard/C342.json)).

**Live split values.** kt4 enumerated split values on the product of per-class lattice steps. For complete lists of
very different sizes (1512 against 6 · 10⁶) almost none of those points is reachable from both sides: the 1,097 first
values of the run above had pairs on at most one side. `lowe_sz.sigma_file` now (a) computes the sumsets
{a_S + c_S} and {T_S − b_S − d_S} of the distinct merge-class values exactly when they are small enough and passes
only their intersection, ordered by the solution law, so that the join is exhaustive when all of them are processed,
and (b) otherwise passes the sampled live split values and the lattice points of the (A + C) coset nearest the
solution law that also lie in the (B + D) coset. On ℤ₉, tuple 0 of 171/⟨31⟩ (pools 1134,
756, 4.3 · 10⁵, 4.3 · 10⁵) has 38,086 live split values; all of them, 3.4 · 10⁸ + 2.2 · 10⁸ pairs, give no solution,
against 10^0.48 predicted.

**Calibration from the searches.** kt4 stops at its solution cap, so the N = 215 run says little about the size of the
class beyond its density: 8 solutions in the first 574 split values (1.5 · 10⁷ pairs per side), consistent with the
10^2.29 predicted. N = 453 (pools 2520, 2520, 1.6 · 10⁷, 1.6 · 10⁷, 10^1.58 predicted): every split value was live on
both sides, and all 2 · 10⁵ of them (3.8 · 10⁹ pairs per side, about 9.5% of the A × C pairs, the inner part of the
solution law) gave no solution. Sampling the class sums of that tuple shows no pair or triple of classes that never
vanishes jointly, so no low-order obstruction explains it; the finer quotient ℤ₁₅₁ (H of order 15 acting on it) is not
listable by the scan and is the next suspect.

## 45. The character-bordered Szöllősi array: order 214

**Construction** (`search/tbsz.py`, [TBORDERED.md](../complex-hadamard/TBORDERED.md)): the complex analogue of §18.
On G of order N − t with φ: G → Q and a conjugation-closed set T of t characters of Q, the core condition is
Σ_k PAF_k(g) = −4 Σ_{χ∈T} χ(φ(g)) off 0, and the border is a pair of quaternary bases (R̂ rows, Ĉ columns) of the
2t-dimensional space V_T of Q-functions (two halves) with Fourier support in T such that Xb = −R̂ Mc* Ĉ/2q is
quaternary. Q = ℤ₄ with T = {1, χ, χ̄} needs N ≡ 3 (mod 4) and makes the order-N − 3 groups available to the primes
where ℤ_N and ℤ_(N−1) have no multipliers.

**Engine.** `TStruct` (any G = ℤ_{m₁} × … with automorphism matrices, the target per ⟨H, −1⟩-class, row patterns:
Σ r² = 4t when 1 ∈ T and 4N otherwise) and `TComp` (compressions to Q, feasible tuples, conditioned lattice Gaussian
with per-class targets, as `gcomp.py`); the search reuses `lowe_enum` with fibre constraints and kt4 with the target
vector (`lowe_sz.sigma_file` and `_exact_live` now take a target vector); every hit is completed by the border search
and checked in integer arithmetic. `tbsz.py rank N --family z4cyc|v4n|v4t` ranks the structures of one family.

**N = 107.** `rank 107 --family z4cyc` (G = ℤ₈ × ℤ₁₃, Q = ℤ₄): ⟨(1, 2)⟩ (order 12) and ⟨(1, 4)⟩ (order 6) are empty
on ℤ₄; ⟨(1, 3)⟩ (order 3; 40 and 21 orbits, 20 classes) has conditioned 10^5.42 for rows 2, 0; 2, 2 and 10^4.36 for
2, 2; 2, 0 (plain Gaussian 10^4.6 for both). The search took tuple 0 (y_A = (0, 2, −2, 2), y_B = (2, −2, 2, −2),
y_C = (−4, 4, −4, 6), y_D = (0, 0, 2, 0) on ℤ₄; pools 5040, 2520, 1.97 · 10⁷, 1.92 · 10⁸; 10^3.79):
- lists: A, B, C complete; D 2.94 · 10⁷ of 1.92 · 10⁸ (lowe_enum 45 s, retention 0.50); 10^3.48 predicted inside;
- split values: 52,762, exact (sumset intersection on merge classes 17, 2, 12, 7, 18);
- kt4: 5 solutions (its cap) in the first 256 split values, 1.3 · 10⁸ + 8.9 · 10⁷ pairs, 338 s from the start.

All five complete, and the first is [C214.json](../complex-hadamard/C214.json). With 10^3.48 predicted and the
first 0.5% of the split values (by the solution law) giving 5, the estimate is consistent with the run.

**A convention pitfall.** `spence_lib.AbGroup.index` reverses all axes of a 3-dimensional coordinate array, so blocks
developed through it are X[g, h] = x(g − h). Such matrices are still complex Hadamard (the same construction
for the reflected sequences, with the border solved for the reflected compressions), but a certificate read in the
standard convention by `verify_complex.py` fails. `tbsz.devidx` therefore flattens the coordinates first, and
`certify` solves the border from the core sequences.

**Other N ≡ 3 (mod 4).** See [NOTES.md](../NOTES.md) §16: N = 167 is empty in both families on ℤ₄ × ℤ₄₁ and ℤ₂² × ℤ₄₁
(odd fibres of 41 points); N = 179 on ℤ₁₆ × ℤ₁₁ lives only for ⟨(1, 3)⟩ at 10^−1.9; N = 239 is empty in both families.
No result for N = 223 in these families is recorded.

**After order 214** (`tbsz.py rank`, families z4cyc, z3cyc, z6s/z6t, v4n/v4t, all3; `--cond` conditions a search on a
finer quotient that refines Q):
- Q = ℤ₆ with T = all characters but one (t = 5, any N ≡ 5 mod 6) has 400 quaternary vectors per half of V_T. The
  randomised border search (`solve_border2`, 200 paired column bases per tuple) found no border for the feasible
  compression tuples of N = 161, H = ⟨29⟩; the exhaustive search of §46 finds one for each of them at once. The
  family is the one that gave order 322 (§46).
- N = 249 (Q = ℤ₃, ℤ₂₄₆): H = ⟨43⟩ of order 20 (18 orbits) has 9 feasible tuples for rows 2, 0; 2, 2, and in all
  of them the merge-class sumsets of the complete lists (72, 144, 864, 864 members) do not meet: no solution.
  H = ⟨37⟩ of order 5 (54 orbits) rates 10^5.7 (best tuple 10^5.0) but with C and D pools of 8.5 · 10¹².
- N = 237 (ℤ₂₃₄): H = ⟨43⟩ of order 6 (46 orbits), best tuple 10^2.45 with pools 2.4 · 10⁵, 3.3 · 10⁵, 4.7 · 10¹⁰,
  4.7 · 10¹⁰; at L = 3 · 10⁷ the retention of C and D is 0.06 and 0.04 (10^−0.4 left in the lists).
- N = 207 (ℤ₂₀₄): H = ⟨13⟩ of order 4 rates 10^2.5 with pools 3 · 10¹⁴; N = 203 (ℤ₈ × ℤ₂₅, Q = ℤ₄): the |H| = 4
  groups rate up to 10^3.3 with pools near 10¹²; every |H| = 10 group (37 to 45 orbits) has feasible tuples but no
  estimable one (plain Gaussian near 10^−500).
- With about 40 orbits (complete lists) the estimates fall from 10^5.7 at N = 107 to about 10^1–10^2 near
  N = 230: the gain per class is 6 − log₂(√(2π · 6m)/4) bits, which shrinks with m. The remaining leads need finer
  conditioning, which splits the pools.

## 46. The ℤ₆ character-bordered array: order 322

*Status of the orders of §§41–46.* None of 322 and of 330, 342, 430 and 454 (§§41–44) is open: 4N − 1 is a prime power
for each, and Ito's Hadamard groups with Egan 2019, Theorem 4.5, give them ([NOTES.md](../NOTES.md) §18,
[PALEY.md](../complex-hadamard/PALEY.md)). The searches and their calibrations are unaffected. The 12 orders 214, 238,
302, 326, 366, 382, 426, 494, 518, 534, 754 and 1558 constructed by the searches were obtained concurrently by Zaitsev
et al. ([NOTES.md](../NOTES.md) §20).

**Equivalences.** Multiplying the sequences by λ ∘ φ (λ a real character of Q) maps the family of T to that of λT,
border included, so z4n ≡ z4cyc, v4n ≡ v4t and z6s ≡ z6t (the small-group core counts agree: 32768, 49152 and
1,419,264 for both T on ℤ₄). t is odd in every family with a quaternary border (row sums for 1 ∈ T, zero-sum border
rows for 1 ∉ T), and WLOG 1 ∈ T.

**Exhaustive border search** (`border_exh.c`, `tbsz.solve_border_exh`; now the first choice in `complete` and in the
per-tuple border check of `search`). Candidates: the normalised quaternary vectors of V_T ⊕ V_T (40,000 for ℤ₆).
The caller passes P = 2 v Mc* (Gaussian integers), and the C code builds two n × n bit matrices: B (u Mc* c ∈ 2q · μ₄,
exact integer test) and O (orthogonality from exponent differences), OpenMP over rows. A row DFS (rows ordered,
each orthogonal to the earlier ones) carries the bitset of columns compatible with all chosen rows and prunes when
fewer than 2t remain; at depth 2t it looks for a 2t-clique of O among them. For ℤ₆ the bit matrices take about a
minute (2 · 10¹⁰ integer operations) and the search 0.0–0.1 s: compatibility degree 1100–1200, orthogonality degree
3432 per candidate, solutions after 137–1898 nodes. The first row found is the all-ones vector of both halves in
each case. Neither the rows nor the columns of the borders found are of the paired form (v_i, ±c w_i) that
`solve_border2` tries, which is probably why it failed. Re-certifying order 214 with the exhaustive search gives the
same border and the same SHA-256.

**Structures beyond cyclic G** (`tbsz.py rankg N --group m_1,... --quot z6:i,j|z4:i|v4:i,j|z3:i`). G is given by
prime-power factors, every automorphism acting on each Sylow subgroup (all matrices with the divisibility pattern of
Hom(ℤ_(p^a), ℤ_(p^b))) is tried as a generator of H, kept if it preserves the target and −1 ∉ H, and deduplicated by
the orbit partition and then by an invariant (orbit sizes with their Q-images, for H and ⟨H, −1⟩). The invariant can
merge two inequivalent groups, so the ranking is not a proof of absence.

**N = 161.** Cyclic ℤ₁₅₆ (`rank 161 --family z6s`): ⟨29⟩ of order 6 is the only group alive, 10^0.34 (rows 4, 2;
0, 0), 10^−0.29 and 10^−0.76. On ℤ₂ × ℤ₂ × ℤ₃ × ℤ₁₃ with φ = (x₀, x₂ mod 3) (11 inequivalent H with 30 to 66
orbits), M = (x₀, x₁, x₂, x₃) ↦ (x₀, x₀ + x₁, 2x₂, 3x₃) of order 6 (35/20 orbits, 19 classes) rates 10^1.10 (0, 0; 4,
2), 10^1.94 (2, 0; 4, 0) and 10^1.86 (4, 0; 2, 0); plain Gaussian 10^−3.7. The search of 2, 0; 4, 0 (8 feasible
tuples, 10^2.51) took tuple 0:
- y_A = (0, 2, −2, 2, −2, 2), y_B = (−2, 2, −2, 2, −2, 2), y_C = (2, 0, 0, 2, 0, 0), y_D = (2, −4, 4, −2, 4, −4);
- complete lists 864, 432, 8,232,000, 1,317,120 (2 minutes, most of it the border check of the tuple);
- 1367 exact live split values on merge classes 5, 2, 8, 11, 14; kt4 found 1024 solutions at split value 246
  (79 s; 2.5 · 10⁸ AB and 1.7 · 10⁷ CD pairs up to there) and stopped at its cap of 5 reported;
- the 5 complete with an exhaustive-search border and verify: [C322.json](../complex-hadamard/C322.json)
  (`verify_complex.py 322`). The 1024 solutions in one split value (2¹⁰) suggest ten independent sign symmetries of
  the class sums.

**Before the fix.** The first run skipped tuple 0 with "no border found": `search` still checked borders with
`solve_border2`. With the exhaustive check the same run found the matrices.

**Two long runs stopped without a solution.**
- N = 203 (order 406), bordered array on ℤ₂₀₂ with ⟨95⟩, rows 2, 0; 0, 0 (`lowe_comp.py`): tuple 0 kept lists of
  63,504, 63,504, 3 · 10⁷ and 2.65 · 10⁷ (retention 0.008 and 0.004 of pools of 2.8 · 10¹⁰), 10^2.27 predicted inside.
  kt4 went through the first 3028 of 200,000 split values in order of the solution law, about 7% of the predicted
  mass (≈ 13 solutions), and found none.
- N = 237 (order 474), Q = ℤ₃ on ℤ₂₃₄ with ⟨43⟩, rows 2, 0; 2, 2, conditioned on ℤ₉: tuple 0 kept 1.3 · 10⁵,
  2.7 · 10⁵, 3 · 10⁷ and 3 · 10⁷ (retention 0.25 and 0.22), 10^−0.17 predicted inside, 254,175 live split values.
  Each split value pairs 10⁷ AB with 5 · 10⁷ CD candidates (5 s); the first 72 gave nothing, and the run was stopped
  as too slow for its expectation.

**167 and 179.** Odd fibres (27, 29). The bare ℤ₆ compression equation has 59 and 82 PAF-class solutions, but every
cyclic structure in range is empty on ℤ₆, as is the first on ℤ₂ × ℤ₉ × ℤ₉. For 179 with t = 3 on ℤ₄ × ℤ₄ × ℤ₁₁ (Q =
ℤ₂², H of order 10) the single feasible tuple of rows 2, 2; 2, 0 (complete lists 192, 192, 12096, 725760) has disjoint
merge-class sumsets: no solution. The rankings with t = 7 over Q = ℤ₂³ (167; fibres of 20) exceed memory in
`TComp.compressions` (8 fibre classes of up to 21 values each) and need a coarser conditioning quotient.


## 47. Refined compression tuples, and a biased lattice in the estimators

**Refinement** (`tbsz.py search --refine K [--pcond ...] --cond ...`). `search --cond` used to rank the compression
tuples on one quotient and list the blocks of the best, keeping the 3 · 10⁷ lowest-energy members of a larger pool.
Now the best coarse tuples (on Q or on `--pcond`) are refined one at a time on the finer quotient `--cond`:
`TComp.compressions` enumerates only the fine compressions whose sums over each coarse fibre equal the coarse values
(classes ordered by the coarse points they meet; a point is checked exactly once all its classes are set, and bounded
by the remaining range before that), and the fine tuples of all refined coarse tuples are ranked together. Coarse
tuples are tried until K of them have a live refinement. Borders (t > 1) are solved per coarse tuple before refining
and cached by the compressions on Q. `--by-rate` takes the fine tuples in order of estimate per join work (the product
of the C and D pools plus a fixed cost) instead of estimate. Check: for an N = 209 tuple, 3000 sampled sign vectors
per block with the coarse values all have fine compressions in the refined set except those of energy
Σ_a y(a)² ≥ 784 = tgt(0), which the exact energy bound (one block ≤ tgt(0) − 3 · 16 for odd fibres of 13) excludes.

**Lookup join.** `TComp.feasible` joined the sorted A + B keys against all C + D pairs. When one pair side is small
it now iterates over that side's pairs, vectorises over a third list and looks up the fourth (cost n_A n_B n_C instead
of n_C n_D): for an N = 167 refinement (7, 7, 29,809 and 33,459 PAF groups) 460 s become 1 s, with the same 41
feasible groups.

**Biased lattice (fixed).** `comp_est.block_stats`, `gcomp`, `tsz_comp` and `TComp.stats` took the lattice of a
block's class values from `echelon((U[1:] − U[0])[:120])`, U = `np.unique` of the samples. Those are the
lexicographically first 120 vectors, which share their leading coordinates, so their differences can span a proper
sublattice. On the N = 167 structure below (fine level) this gave index 2 and index 6 sublattices for blocks A and B
and once a lattice of rank 17 instead of 18. A sublattice makes `in_lattice` fail where the target is attainable (a
false −∞) and inflates the covolume factor where it is not. `comp_est.diff_lattice` now takes up to 400 differences of
distinct samples in random order from a random base point; all four modules use it. Estimates printed before this fix
can be off in either direction: the ℤ₂ × ℤ₃⁴ structure of N = 167 rated 10^2.41 on ℤ₆ before and 10^1.2–10^1.4 after.
Exhaustive results are unaffected, but the v = 65 ⟨21⟩ explanation of §30 ("the target is off the coset of the
lattice", from sampled differences) has to be re-derived.

**Refined ranking** (`tbsz.py rankg ... --refine K`). For each structure, `invariant_refinements` lists the
coordinate-wise quotients ∏ ℤ_(d_i) whose kernel is H-invariant and inside ker φ; the best K coarse tuples (one per
feasible group, C/D and A/B swaps identified) are refined on the smallest of them and the refined estimate is printed
with the number of tuples that stay alive. Calibration: the N = 161 structure of order 322 refines on ℤ₂² × ℤ₃ from
10^1.89 to 10^2.66 (rows 2, 0; 4, 0), every refined tuple alive.

**N = 209 (order 418), bordered array on ℤ₄ × ℤ₄ × ℤ₁₃, H = ⟨(x, y, z) ↦ (x + 2y, 2x + 3y, 3z)⟩ of order 6, rows
2, 0; 0, 0.** The four best tuples on ℤ₄ × ℤ₂ (10^1.80, 10^1.77, 10^1.36, 10^1.31, each with its C/D swap) have no
feasible refinement on ℤ₄ × ℤ₄, so no solution has those compressions; the earlier low-retention run on the first was
stopped. Each of the next coarse tuples (10^1.2–10^1.6) refines into 1536–4096 fine tuples with complete lists (864,
3456, 1.8 · 10⁷, 1.2 · 10⁷), rated about 10^−1.3 each (with the biased lattice). Five fine tuples were joined
exhaustively (kt4 over all live split values, about 14 minutes each) without a solution. Taking A and B on the coarse
quotient to cover several fine tuples per join does not pay: a fine (C, D) pair goes with only about two fine (A, B)
pairs, and the larger A and B lists made one join three times slower per fine tuple.

**N = 167 (order 334), ℤ₆ family on ℤ₂ × ℤ₃⁴, H = ⟨M⟩ of order 4** (M = (x₀, 2x₁, x₄, x₃, 2x₂), 48/28 orbits, 27
classes; refined on ℤ₂ × ℤ₃ × ℤ₃ via x₀, x₁, x₃). The coarse tuples with y_A = (3, 5, −3, −3) for rows 4, 2; 0, 0
(10^0.92–10^0.94) have 41 feasible fine groups (10,008 fine tuples), but every fine tuple tested has its target
outside the lattice of attainable class sums. This holds with 20,000 samples per block and 3,000 random
differences, so it is a congruence obstruction at the ℤ₂ × ℤ₃² level (shown on sampled lattices, not
proved). In the refined ranking all four row patterns stay alive: 10^1.06 (0, 0; 4, 2), 10^1.41 (2, 0; 4, 0),
10^1.34 (4, 0; 2, 0) and 10^1.88 (4, 2; 0, 0). For rows 2, 0; 4, 0 the first coarse tuples refine into 3000 live fine
tuples each (about 10^0.8 per coarse tuple).

**Symmetries of the fine tuples.** Joins are exact, so a tuple needs to be joined only once per orbit of the
symmetries that map core solutions to core solutions. `search` now builds these orbits (`tuple_orbit`) from four
kinds of move:
- negating a block of row sum 0;
- swapping C and D, or A and B, when their row sums agree;
- automorphisms of G given by `--sym` that normalise H and preserve the target (here x₃ ↦ −x₃);
- translating a block by an element of Fix(H). This keeps its PAF and its H-invariance; for the symmetric blocks
  A and B only translations of order ≤ 2 are allowed.

On the N = 167 structure Fix(H) = {(x₀, 0, 0, x₃, 0)} ≅ ℤ₂ × ℤ₃. One coarse tuple's 16,128 fine tuples fall into
224 orbits. An orbit counts as searched only after one of its tuples has been joined. A tuple whose border check
fails does not mark its orbit, because negations change the compressions on Q.

The first runs lacked these orbits. They joined 353 fine tuples of rows 2, 0; 4, 0 and 77 of 4, 2; 0, 0
(10^−2.0 and 10^−2.2 each, about 40 s and 3 min per join after the kt4 fix below), all without a core solution. Since
translates sit next to each other in the estimate order, these joins cover far fewer distinct orbits than their
number.

**kt4 overhead.** With complete lists of 10⁵–10⁶ members, a join has 10⁴–10⁵ live split values of a few hundred
pairs each. kt4 opened two OpenMP parallel regions per split value, and with threads spinning on a loaded machine
this cost about 16 ms per split value: 11 minutes per join. The regions now run in parallel only above 65,536 pairs,
and `search` starts kt4 with a passive wait policy: about 50 s per join.

**Calibration on the solved structure (N = 161, order 322).** The same pipeline (rows 2, 0; 4, 0, three coarse
tuples refined on ℤ₂² × ℤ₃, orbits as above, `--count-only --maxsol 100000`) joins 12 orbit representatives of
36 live fine tuples:

| Estimate per tuple | Representatives | Core solutions each |
|---|---|---|
| 10^0.69–10^0.71 | 4 | 96 |
| 10^0.33–10^0.39 | 8 | 0 |

So the estimates do not overstate the counts: 384 core solutions against about 39 predicted. The solutions are
clumped: a tuple with any solution has 96 of them. These are probably images under automorphisms that fix its
compressions, such as the multipliers of ℤ₁₃ that commute with H, together with signs. The chance that a tuple holds
a solution is therefore its estimate divided by the clump size rather than the estimate itself. For N = 167 the
refined mass per row pattern (about 10^1.5) spread over orbits of 16–72 tuples and such clumps leaves well under one
expected solution class per pattern. That is consistent with the 1,000+ empty joins, and it means a structure for 334
needs about ten times more refined mass than this one.

**Refined sweep (fixed lattice).** `rankg --refine 3` was run over 21 group and family choices for N = 167, 179, 203,
207 and 209. Two structures looked strong on their first refinement level.
- 209, bordered array on ℤ₂ × ℤ₈ × ℤ₁₃ with H = ⟨diag(1, 5, 5)⟩ of order 4 (60/38 orbits): 10^4.43 for rows
  2, 0; 0, 0.
- 203, ℤ₂² family on ℤ₂³ × ℤ₅² (37 classes): 10^3.70 for rows 2, 0; 2, 2. The ℤ₅² block of H is conjugate to
  diag(1, 3), so ℤ₂³ × ℤ₅ is a coordinate refinement.

Followed one level further down (`search --chain`, multi-level refinement), the mass falls.
- 203: from ℤ₂³ to ℤ₂³ × ℤ₅ along the axis where H acts by ×3, the two best coarse tuples carry 10^−0.12 in total,
  with C and D pools of 10^10.4–10^11.1.
- 209: the best 40 tuples on ℤ₂ × ℤ₂ refine into 58,920 live tuples on ℤ₂ × ℤ₄ (10^2.47 in total, all tuples
  counted). The best 40 of these refine on ℤ₂ × ℤ₈ (fibres of 13) into 7,360 live orbit representatives that
  carry 10^−1.21 in total. The best is rated 10^−4.75, with complete C and D pools of 7.9 · 10⁸ each, and lists
  of 3 · 10⁷ would keep 10^−6.3 of it. This structure was not searched.

The Gaussian on a coarse quotient does not see the congruences that bind on the finer ones. The level that matched
reality in the calibration above is the last one, one step from G: fibres of 13 for 161, as for 209 on ℤ₂ × ℤ₈.
Structures have to be compared there.

## 48. Paley plus three blocks over the census primes v ≡ 3 (mod 4)

Đoković and Kotsireas (arXiv:1511.08734) build skew-Hadamard matrices of order 4v, for v ≡ 3 (mod 4) prime, from
Paley's block (−1 on the nonzero squares, PAF −1 off zero) and a (v; k₁, k₂, k₃; λ) difference family X₁, X₂, X₃:
three ±1 sequences with PAF sum +1 off zero. `skew_paley3.py` searches this class with the three blocks unions of
orbits of a subgroup H of the squares, of odd order d. That leaves m = (v − 1)/d orbits besides {0}, and m ≡ 2 (mod
4).

**Method.** `exhaust v d` works as follows.

- **Encoding.** The PAF at a shift class {±gᵏH} is a quadratic form cᵀNₖc in the orbit signs. The PSD is bounded by
  3v − 1, using Gauss periods.
- **Symmetry.** The sign of {0} is fixed in each block (negation). The first block is reduced modulo the squares,
  which permute the orbits and fix Paley's block.
- **Join.** The third block is joined by hashing.
- **Validation.** It finds the factors of [skew/9424](../skew/9424) (v = 31, ⟨2⟩) and [skew/13072](../skew/13072)
  (v = 43, ⟨4⟩), and v = 71 under the subgroup of order 5. Under the subgroup of order 17 it finds nothing for v =
  239, where Đoković and Kotsireas used order 7.

**Row sums.** A block has row sum c₀ + d·e with e ≡ m (mod 2), and r₁² + r₂² + r₃² = 4v − 1. This removes most
subgroups: of the structures with 22 to 58 orbits over the census primes, eight survive.

**Exhaustive and empty.** Every census-undecided prime v ≡ 3 (mod 4) and every subgroup with 6, 10, 14 or 18 orbits
gives nothing. That is 97 structures over 65 primes:

- 49 with 6 orbits, 20 with 10 and 13 with 14 (`sweep 6 14`);
- 15 with 18 orbits, for v = 487, 523, 739, 883, 919, 991, 1063, 1423, 1567, 1747, 1783, 1999, 2143, 2251 and 2287
  (`sweep 15 20`).

None of these structures has 2, 4, 8, 12 or 16 orbits, so these four sizes cover the whole range up to 20. The other
39 undecided primes ≡ 3 (mod 4) have no subgroup in this range.

**Counts.** `estimate v/d` uses a Gaussian model:

- the PAF vector of each block is normal, with mean and covariance sampled at each row sum;
- the values lie on a lattice of spacing 4;
- the direction fixed by the row sums (Σₛ PAF(s) = r² − v) is removed.

For 31/5 and 43/7 the model predicts 10³·⁰⁹ ordered triples, and the exact count is 1008 (`p3count`, scratch).

| v/d | orbits | log₁₀ E | log₁₀ density per row-feasible triple |
|---|---:|---:|---:|
| 239/7 (Đoković–Kotsireas) | 34 | 10.0 | −18.6 |
| 1123/33 | 34 | 4.4 | −24.1 |
| 1483/39 | 38 | 3.7 | −28.4 |
| 251/5 | 50 | 14.2 | −28.7 |
| 883/21 | 42 | 6.6 | −29.4 |

All other undecided primes up to 1500 give densities below 10⁻³⁴.

**Local search does not reach them.** Simulated annealing and tabu search over the orbit signs (C, scratch) solve 71/5
within a second. They do not solve 239/7 in about ten CPU-minutes: the best energy was 48, a few classes off. The
undecided structures are at least 10⁵·⁵ sparser. The 2 + 1 join with sampled lists costs about density^(−2/3), which
is 10¹⁶ for 1123/33. So this class does not settle any census value without a better matcher.

## 49. Symmetric matrices of order 4q: inputs, census omissions, and no odd multiplier

**Inputs of the 4q family** ([../symmetric/PROPUS_MIYAMOTO.md](../symmetric/PROPUS_MIYAMOTO.md)). Of the 58
census-undecided prime powers q ≡ 1 (mod 8), 44 take a symmetric Hadamard matrix of order (q − 1)/2 from Paley I or II
and Sylvester, and 1753 takes symmetric/876. Three take published sporadic matrices, transcribed from the LaTeX
sources and rebuilt by `verify.py` in the propus array of `build_propus` (input kind `propus`): order 116 (Di Matteo,
Đoković and Kotsireas, arXiv:1503.04226v2, the first of their fifteen SDSs), 236 (Đoković, arXiv:2101.05429v2, the
first PDF for v = 59, X₁ = X₁′ ∪ −X₁′) and 604 (Balonin and Đoković, arXiv:1803.08787, Table 1, orbit representatives
under {1, 8, 19, 59, 64}). Every published solution of those papers checked (15 for v = 29, 6 for v = 59, 1 for v =
151). The last ten values need a symmetric matrix of order 4m with m = 101, 107, 119, 149, 179, 239, 249, 269, 299 or
305, all undecided.

**Seberry–Balonin with conference cores** (`sym_propus_sb.py`). Mathon's conference matrix of order 442 comes from
SageMath's `MathonPseudocyclicStronglyRegularGraph(2)` (29 s; the same matrix on a rerun), the one of order 226 from
the skew-Hadamard matrix of order 16 (Paley 8, doubled) as Q = S ⊗ S + I ⊗ J − J ⊗ I. With Yamada's circulants for
q = 881 and 449 they give 3532 and 1804. A mutated core is rejected by its SHA-256, and a changed ξ by the primitivity
test. The variant with 2n − 1 an odd square instead of a prime power would need new symmetric circulant pairs
(M, N) with M² + N² = s²I: among the census values it would reach only v = 1227 (n = 613, s = 35) and 2027 (n = 1013,
s = 45), and the lattice count for multiplier-invariant pairs, about 2 · (4c/√n)^k for k orbits, is far below one.

**No odd t multiplies the propus array** (`sym_tmult.py`). Cooper–Wallis multiplies four Williamson-type matrices by
T-matrices; the symmetric analogue would write X_k = Σ_j M_kj ⊗ Y_j (k, j ∈ {a, b, d}) with symmetric commuting (0,
±1) M_kj of disjoint supports in each row and Σ_k c_k M_kj M_ki = t c_j δ_ij for c = (1, 2, 1). At the trivial
character the row sums give an integer M with Mᵀ diag(1, 2, 1) M = t diag(1, 2, 1), so det(M)² = t³ and t must be a
square. Exhaustive searches over the symmetric elements of ℤ₃, ℤ₅, ℤ₉ and ℤ₃² find no solution (the control t = 1
finds the 16 signed permutations). Four distinct symmetric blocks do not help either: a symmetric orthogonal design of
order 4 in four variables would need four anticommuting real symmetric involutions of order 4, and their Clifford
algebra M₂(ℍ) has no module of dimension 4. Twisting by the reversal does not rescue it (`sym_arrays.py`). The blocks
XR and RX = XᵀR are symmetric for every circulant X, but they do not commute, so the Clifford argument does not apply
to arrays built from them. An exact symbolic check settles it: every symmetric Latin square of order 4 (96), every
orientation and sign, each of the four circulant variables symmetric (blocks ±X, ±XR) or arbitrary (±XR, ±RX), H²
compared with (Σ XXᵀ) ⊗ I by cancelling monomials. No array with four distinct variables exists, whichever variables
are symmetric. With X₃ = X₂ there are 768 arrays when only X₁ is symmetric (Balonin's among them), 1152 when X₁ and X₄
are, 6144 when X₁, X₂ and X₄ are, and none when X₁ is arbitrary. So every symmetric Hadamard matrix of Goethals–Seidel
type over ℤ_v needs a repeated block and a symmetric block, as in the propus family. Twisting by any involutive
multiplier does not help either (`--klein`). For v with two odd prime factors, ℤ_v* contains u ≠ ±1 with u² = 1, and X
P_u is symmetric when X is invariant under the multiplier −u. Over all 35 type assignments in {1, −1, u, −u}, none of
the arrays has four distinct variables. The twisted repeated-block arrays that do exist (448 with D of type u) impose
the same autocorrelation equation as propus with extra symmetry, so they reach nothing new.

**Corollary 27(i) with a skew core, and length 37** (`sym_skewcore.py`, `sp_mitm.c`, `l26_scheme.py`). The proof of
Seberry–Yamada's Corollary 27(i) uses of Paley's skew core only a skew core B with B R symmetric for the reversal R
(NOTES.md §29). A skew-propus quadruple of length ℓ gives, in the Balonin array, a skew-Hadamard matrix S = I + NP of
order 4ℓ with the involutive anti-automorphism P, which fixes two points. Normalised at one of them, its core,
relabelled on ℤ_n (n = 4ℓ − 1) so that P acts as x ↦ −x, is such a B. With Yamada's circulants for the prime
q = 2n − 1 = 8ℓ − 3, Lemma 26 with U = (I, 0, BR, BR) and V = (M, N, R, R) gives symmetric Williamson-type matrices of
order q + 2 with W₃ = W₄, and the propus array a symmetric Hadamard matrix of order 4(q + 2) (verify.py construction
`symmetric-propus-skewcore`; the lengths 13 and 19 give orders 412 and 604, checked). The census values in reach are
v = 8ℓ − 1 with q prime and a₃² + 2a₁² = 4ℓ − 1 solvable: ℓ = 37 (v = 295), 85 (679), 103 (823), 133 (1063), 187
(1495), 235 (1879), 277 (2215) and 295 (2359). Of these, 1879 and 2215 are also FRONTIER orders, through Kiyasu from
skew 940 and 1108.
- Length 37 has no skew-propus quadruple. `sp_mitm.c` is exhaustive. The row sums must be (7, 7, 1). A3 runs over its
  1768 classes under the multipliers and A0 over all 2¹⁸ skew-type sequences. The 2,937,404 pairs that pass the
  power-spectrum test |Â3|² + |Â0|² ≤ 148 are tabled by the autocorrelation A1 would need. All 3.65 · 10⁹ sequences
  A1 with row sum 7 and A1(0) = −1 (every A1 has such a shift) are looked up: no match (3 cores, 14 minutes). A second
  run with the work split over the first two −1 positions of A1 and a Bloom filter in front of the table also finds
  none (4 minutes).
  - Controls: the program finds quadruples at every length 3 ≤ n ≤ 35 for which Abuzin–Balonin–Đoković–Kotsireas
    (2019, Table 1) give one, and at n = 41 (row sums 1, 9; their (41; 20, 16, 16, 20; 31)).
  - It finds none for n = 25 with row sums (7, 5) or (1, 7), and none for n = 31 with (11, 1), the patterns for
    which their Table 3 has only a fourth block that is neither symmetric nor skew.
  - For 37 their Table 3 also has only such a family, under ⟨10⟩; no family with X₄ symmetric exists. So v = 295 is
    not reached this way, and neither are the amicable order 148 (16 · 851 in NOTES.md §28) nor the m(m − 1) order
    148 · 147.
  - The counts are Poisson-like. A random model (2^{2ℓ} sign choices, each of the ℓ − 1 independent autocorrelation
    classes 0 with probability 4/√(2π · 8ℓ)) predicts about 10^2.5 solutions per row pattern for 13 ≤ n ≤ 37. The
    exhaustive counts are about 10^3 when non-zero, that is one to four classes under multipliers, shifts of A1 and
    A0 ↦ A0ᵀ, and zero for four of the fourteen patterns.
- The longer lengths are out of reach. The same model gives 10^−1.2 solutions in all at ℓ = 85, 10^−8 at 133 and
  10^−28 at 235, and the only odd-order multiplier group of ℤ₈₅ is trivial, so these lengths need an algebraic
  family.
- The cyclotomic form of Lemma 26 with a repeated block (`l26_scheme.py`: U_k, V_k in the cyclotomic scheme of GF(p)
  of index e with −1 in the subgroup, v = 2p + 1) gives only p = 5 (v = 11) and p = 13 (v = 27) for p ≤ 400 and
  e ≤ 8, and nothing for 200 ≤ p ≤ 1250 and e ≤ 4.
- **Afterwards ([NOTES.md](../NOTES.md) §31).** The repeated block, and with it the skew-propus quadruple, turned out
  to be unnecessary: Yamada's data in quaternary form give a symmetric matrix for any skew core. Of the eight values
  above, 295, 679, 823, 1063 and 1495 are now built that way (skew inputs of orders 148, 340, 412, 532 and 748); 1879,
  2215 and 2359 wait for skew-Hadamard matrices of orders 940, 1108 and 1180.

## 50. Complex orders 342, 426 and 1558: the low-energy k-tree and the end of the exact sweep

Orders 342, 366, 430 and 454 were also found by the compression-conditioned route of §42–§44. This section and §51–§56
reach the same orders by the low-energy k-tree (fibre filters enter in §55 and are used only for the bordered
order-366 run of §56). Order 426 was also reported in §41 (same structure, 72 solutions). The certificates in
[complex-hadamard/](../complex-hadamard) are those of §41 (426), §42 (366 and 454), §43 (430) and §44 (342).

**Split values.** kt4 pairs its lists 0 and 1 (A + C) against its lists 2 and 3 (B + D), so a solution's split value
is a value of a + c, on the lattice of a + c (at N = 183, a fixed residue class modulo 20 on the first merge class).
`sigma_file` in `lowe_sz.py` draws the split values from the law of a + c. With this order the N = 131 control finds
416 of the 2560 solutions inside its lists in the first 3000 split values, against 400 in kt4's own order.

**Order 342.** N = 171 = 9 · 19, H = ⟨11⟩ of order 6 (33 orbits; ⟨H, −1⟩ has 17 classes), rows (1, 11; 11, 21).
- The exact join is out of reach: A and B have 7128 and 6864 members, C and D 1.5 · 10⁸ and 2.6 · 10⁸.
- The filter kept 3 · 10⁷ members of C and of D, retaining 49% and 38% of a solution's members under the model.
  Estimate inside the lists: 10^0.94 (structure 10^1.67).
- kt4 found 12 solutions at the 814th split value. [C342.json](../complex-hadamard/C342.json) has the rows (11, 11; 9,
  19) of §44; a rebuild as matrices is recorded for those rows only.

**Order 426.** N = 213 = 3 · 71, H = ⟨37⟩ of order 7 (33 orbits), rows (13, 17; 13, 15). A and B have 6435 members; C
and D share one pool of 4.36 · 10⁸ (one enumeration), of which 3 · 10⁷ were kept. kt4 found 72 solutions at the 40th
split value. Estimate inside the lists: 10^1.1. The rows (15, 17; 13, 13) of the same structure gave nothing in 11,000
split values.

**Order 1558.** N = 779 = 19 · 41, H = ⟨49⟩ of order 60 (order 3 modulo 19, order 20 modulo 41), rows (5, 7; 39, 39).
The exact sweep's lists of 120, 120, 18480 and 18480 hold 96 solutions. The per-class Gaussian estimate of this
structure is 10^−13.6. Like 742 (⟨13⟩ of order 26 on 7 · 53) and 754 (⟨9⟩ of order 42 on 13 · 29), it is a large
diagonal cyclotomic group on a product of two primes, with tiny lists and many solutions.

**The exact sweep is complete** for the open N < 1000 (pair work up to 2 · 10¹⁰ per structure and pattern); 779 was
its only hit in the range 583–995.

## 51. Accuracy of the Gaussian estimate; orders 430 and 454

The per-class estimate of `ssxx.gauss` and `bsz.gauss` multiplies one Gaussian density per class (with the class's
lattice spacing). Against exact counts (complete lists, `join4x`):

| N | H | rows | exact | per-class | full covariance | hybrid |
|---:|---|---|---:|---:|---:|---:|
| 163 | ⟨38⟩ | 17, 19; 1, 1 | 10^3.92 | 10^3.64 | 10^3.92 | 10^3.53 |
| 267 = 3 · 89 | ⟨4⟩ | 1, 1; 21, 25 | 10^2.76 | 10^−1.31 | 10^3.17 | 10^2.79 |
| 143 = 11 · 13 | ⟨3⟩ | 3, 9; 11, 19 | 10^1.20 | 10^−1.47 | 10^−0.60 | 10^0.34 |
| 247 = 13 · 19 | ⟨9, 35⟩ | 17, 19; 7, 17 | 10^1.20 | 10^−3.26 | 10^−1.99 | 10^0.72 |
| 119 = 7 · 17 | ⟨8⟩ | 1, 5; 15, 15 | 10^3.67 | 10^−4.94 | 10^−0.87 | 10^0.68 |
| 185 = 5 · 37 | ⟨9⟩ | 3, 3; 19, 19 | 10^2.11 | 10^−5.12 | 10^−2.68 | 10^0.34 |
| 377 = 13 · 29 | ⟨9⟩ | 27, 27; 5, 5 | 10^2.11 | 10^−8.33 | 10^−4.93 | 10^−0.24 |
| 371 = 7 · 53 | ⟨13⟩ | 3, 3; 25, 29 | 10^1.98 | 10^−9.25 | 10^−5.25 | 10^−3.15 |
| 779 = 19 · 41 | ⟨49⟩ | 5, 7; 39, 39 | 10^1.98 | 10^−13.61 | 10^−7.98 | 10^−4.58 |

For prime N the per-class estimate is good. For composite N it is low, by up to 15 orders of magnitude. The classes of
shifts of small order (the compressions of ℤ_N) carry sums of a few large terms, repeated |H| times, so their joint
distribution is lumpy and correlated, and a product of marginal Gaussians understates its mass at the target.

**Hybrid estimate** (`search/gauss_mc.py`). The summed class vector is sampled from random combinations of uniform
block samples (about 10⁶ combinations of 10⁴ samples per block). Classes are ranked by how much their Monte Carlo hit
rate at the target exceeds the Gaussian one and taken into an exact Monte Carlo estimate while at least 40 samples hit
all of them (two to four classes); the rest get the full-covariance Gaussian, conditioned on those classes, times the
covolume of their joint lattice (`search/gauss_lat.py`, Hermite normal form of sample differences). It is within an
order of magnitude when |H| is moderate (the first four rows), and still far low when |H| ≥ 26 — those structures are
small enough to search exactly anyway. Two more calibrations:
- **N = 93 = 3 · 31, H = ⟨4⟩** of order 5 (21 orbits; the analogue of 183 ⟨34⟩): exact counts of the eight patterns
  10^3.75, 10^3.98, 0, 10^3.45, 10^4.33, 0, 10^3.11, 0; hybrid 10^3.65–10^4.36 for all eight; per-class
  10^2.06–10^2.57.
- **Bordered, N = 31** (ℤ₃₀): ⟨19⟩ exact 10^3.36, 10^3.36, 10^4.36, 10^4.36 (2 in A, B, C, D), hybrid 10^3.39–10^3.67,
  per-class 10^2.83–10^2.90. ⟨11⟩: 0 with the 2 in A or B, 51200 = 10^4.71 with the 2 in C or D; hybrid 10^4.2–10^4.3
  for all four. ⟨15⟩ on ℤ₃₂ (N = 33): 0 for all four patterns, hybrid 10^0.9–10^2.1.

So where solutions exist the hybrid estimate is right to within half an order, but some patterns are empty for
arithmetic reasons that no Gaussian model sees. For N = 31 ⟨11⟩ with the 2 in A, each single compression (to ℤ_d,
d | 30) is satisfiable on its own (joins of the compressed complete lists); the obstruction is joint.

**Results.** The frontier scan (`search/frontier_scan.py`: every structure of every open N with 20–75 orbits and pools
up to 10¹³) re-ranked by the hybrid estimate (`search/hybrid_rank.py`) put the structures below at the top. They were
run with `lowe_sz.py` (3 · 10⁷ members kept per big block, at most 10⁴ split values):
- **Order 430.** N = 215 = 5 · 43, H = ⟨11⟩ of order 7 (trivial on ℤ₅, order 7 on ℤ₄₃), 35 orbits, 17 classes, rows
  (11, 17; 15, 15). Per-class estimate 10^−1.18 (10^−2.35 inside the lists), hybrid 10^4.09. A and B have 6435
  members; C and D share one pool of 1.45 · 10⁹ (one enumeration). kt4 found 200 solutions at the 2967th split value.
- **Order 454**, the first order from the bordered array ([BORDERED.md](../complex-hadamard/BORDERED.md)).
  N = 227 = 2 · 113 + 1 is a safe prime; over ℤ₂₂₆ = ℤ₂ × ℤ₁₁₃, H = ⟨49⟩ of order 7, 34 orbits, rows (0, 0, 2, 0).
  Per-class estimate 10^0.63, hybrid 10^4.4. C and D: pools of 6.0 · 10⁸ and 1.2 · 10⁹, 3 · 10⁷ and 2.77 · 10⁷ kept.
  kt4 found 512 solutions at the 851st split value.

**Counter-examples.**
- N = 183 = 3 · 61, H = ⟨34⟩ of order 5, rows (17, 19; 1, 9): hybrid 10^6.69, per-class 10^4.98 (10^2.04 inside lists
  of 2 · 10⁷ of 2.7 · 10¹⁰). No solution after 6833 split values, which cover 17% of the pair-product mass Σ_σ N_AC(σ)
  N_BD(t − σ). The hybrid estimate is too high here (it would put 10^4.4 solutions in the lists).
- N = 183, H = ⟨20⟩ of order 10 (it contains the element that is −1 on ℤ₃ and 1 on ℤ₆₁): exhaustively searched,
  no solution in any of its 14 patterns, against 10^0.4–10^1.0 each (hybrid) and 10^−0.3–10^−0.6 (per-class).

## 52. Retention of the low-energy filter on exact solution sets

`search/retention_test.py` lists all four blocks of a structure, finds every solution with `join4x`, and keeps the L
members of lowest energy of each big block as `lowe_sz.py` would:

| N | H | solutions | kept (C and D) | model retention (C, D) | solution members kept (C, D) | solutions kept |
|---:|---|---:|---|---|---|---:|
| 163 | ⟨38⟩ | 8352 | 5% | 0.093, 0.098 | 0.079, 0.080 | 0 (product of marginals: 53) |
| 163 | ⟨38⟩ | 8352 | 21% | 0.30, 0.29 | 0.34, 0.34 | 832 (model 728) |
| 119 | ⟨8⟩ | 4704 | 4.8% | 0.32, 0.32 | 0, 0 | 0 |
| 119 | ⟨8⟩ | 4704 | 19% | 0.74, 0.74 | 1, 1 | 4704 |
| 371 | ⟨13⟩ | 96 | 4.9% | 0.28, 0.28 | 1, 0.5 | 48 |
| 267 | ⟨4⟩ | 576 | 13% of C, 40% of D | 0.41, 0.68 | 0, 1 | 0 (model 162) |

The marginal retention follows the model on average, but a structure's solutions share few distinct members (267:
576 solutions use 8, 8, 48 and 48 members of A, B, C and D), so each structure is a lottery over those members, and
small kept fractions lose them together. Under the Gaussian model the joint retention equals the product of the
marginals (`search/joint_retention.py`: ratio 1.05–1.07 at N = 163); the measured shortfall at 5% comes from this
discreteness.

**Where the solutions sit** (`search/retention_diag.py`, N = 267). The solution members of C lie at the 31st
percentile of the energy and those of D at the 25th. Their values on the compression classes are fixed by the
arithmetic: on the class of the shift 89 (order 3) every solution member sits 1.21 standard deviations from the
filter's centre, where random members average 0.56. With the energy restricted to the unit classes (shifts prime to
N) they move to the 21st and 2.4th percentiles. The compression classes are better handled exactly than by the
Gaussian filter.

**Unit-class energy** (`lowe_sz.py --energy units`: the filter's energy uses only the classes of shifts prime to the
group order; `retention_test.py --energy units`). Solutions kept out of the exact sets, full energy against unit
classes:

| N | H | kept per big block | full energy | unit classes |
|---:|---|---|---:|---:|
| 119 | ⟨8⟩ | 2000 of 105105 (1.9%) | 0 of 4704 | 4704 of 4704 |
| 119 | ⟨8⟩ | 5000 (4.8%) | 0 | 4704 |
| 371 | ⟨13⟩ | 500 of 20592 (2.4%) | 24 of 96 | 48 of 96 |
| 93 | ⟨4⟩, rows 9, 13; 1, 11 | 15000 of 1.4 · 10⁵ (11%) | 1224 of 21384 | 1368 |
| 93 | ⟨4⟩, rows 7, 11; 9, 11 | 15000 of 1.3 · 10⁵ (11%) | 864 of 9504 | 252 |

At 3.8% kept both versions lose every solution of the two 93 patterns, and at 5% the full energy loses all of 163's
8352 solutions although each marginal keeps 8%: the joint retention of C and D falls far below the product of the
marginals at small kept fractions. Both filters ask the two blocks to absorb the same deviation of A + B from its
mean, and a solution's A + B is one of a few discrete values. The runs that found solutions kept 2–20% of their pools;
the 183 run keeps 0.07%. The two energies select different members, so each gives a separate chance.

**Split-value order.** For the 183 run above, the 6833 split values processed (kt4's order by the Gaussian law of a
solution's split) cover 17.1% of the pair-product mass; the best order for that cost, by N_AC N_BD/(N_AC + N_BD),
covers 19.0%. The order is close to optimal.

**Memory.** `lowe_sz.py` keeps only the masks and 2 · 10⁵-row subsamples of the kept lists while kt4 runs.

**Bordered sweep.** `bsz_sweep.py` was run for N = 107 and 155; no outcome is recorded for the sweep over the other
210 open N. `bsz_group.search` handles every pattern of a spin structure, so `plan()` makes one job per spin structure
instead of four.

## 53. Difference-map search (negative) and the decoded matrices of Alpöge et al.

**RRR (Elser's difference map)** for four ±1 sequences with a flat summed power spectrum (`search/rrr.c`,
`search/rrr.py`). Constraint A: ±1 vectors with fixed row sums (and symmetry), projected by sorting. Constraint B: per
frequency, Σ_k |X̂_k(ω)|² equals the target, projected by rescaling. It solves plain Goethals–Seidel quadruples at
m = 19 in about 2 · 10⁵ iterations and at m = 31 in about 3 · 10⁷ (1 of 8 runs of 4 · 10⁶; β = 0.5). The cost grows
by about 1.5 per unit of m, like local search (§26), so it is no route to unstructured lengths above about 40.

**The matrices of Alpöge et al. (August 2026)**, read from the public mirror github.com/bbeartheancient/hoa64
(commit d825201). All twelve are stored normalised (first row and column positive).
- Undoing the normalisation (in a circulant block, x(j) = c · H′[1, j] · x(j − 1)) recovers the Goethals–Seidel
  quadruples of 892, 1132, 1244, 1948 and 1964 (four sequences on ℤ_N with Σ PAF = 0). Under a border of 4 rows the
  same recovers the Spence seeds of 668, 716, 1676 and 1772 (four sequences on ℤ_{N−1}, Σ PAF = −4, row sums (2, 0, 0,
  0)); for even m each block is fixed only up to the character (−1)^t, and the twists were chosen to make
  Σ PAF constant.
- No block of any of them is symmetric or skew about any centre, so none gives a Szöllősi or bordered Szöllősi
  quadruple.
- The Goethals–Seidel blocks are invariant up to translation under multiplier groups of orders 3, 3, 5, 9 and 7, and
  the Spence seeds of 1676 and 1772 under groups of orders 15 and 12.
- The seeds of 668 and 716 (lengths 166 = 2 · 83 and 178 = 2 · 89) have no affine symmetry, no relation between
  blocks, no relation between their CRT slices, no slice correlated with a Legendre sequence beyond random size
  (largest |correlation| 30 of 83 and 89), and per-block spectra only as flat as the sum condition forces. Their ℤ₂
  splits u = (a₀ + a₁)/2, v = (a₀ − a₁)/2 satisfy Σ PAF_u = −2 and Σ PAF_v = 0 on ℤ_q with weights 2q − 1 and 2q + 1,
  as every solution must. They are unstructured solutions of length 166 and 178, far beyond local search or RRR.
- 1388, 1436 and 1916 fit neither layout.

**Safe primes.** The same split applies to the bordered Szöllősi array over ℤ₂ × ℤ_q: it is exact and the parity of
the v-stage is automatic, but with A and B symmetric the v-stage at q = 53 (order 214) has about 82 free signs for 95
bits of conditions, so a two-stage search needs some 10⁴ u-solutions per success, each a length-53 four-sequence
problem. Not feasible with these tools.

## 54. Dead split values and exact compression checks

**Dead split values.** kt4 splits a solution's values on the five merge classes as σ = a + c = t − b − d and, for
each σ, joins the (A, C) pairs with a + c = σ against the (B, D) pairs with b + d = t − σ. `sigma_file` took the σ
nearest the Gaussian law of a solution's split, on the per-class lattice (step and residue of a + c in each class).
Recounted from kt4's progress files, most of them had no pair at all on one side:

| N | H | rows | split values processed | with no pair on one side |
|---:|---|---|---:|---:|
| 341 | ⟨194⟩ | 19, 21, 11, 21 | 20740 | 20740 (100%) |
| 203 | ⟨36⟩ | 3, 25, 3, 13 | 10000 | 10000 (100%) |
| 215 | ⟨11⟩ | 11, 17, 15, 15 | 2967 | 2960 (99.8%) |
| 275 | ⟨4⟩ | 15, 15, 5, 25 | 10000 | 9683 (97%) |
| 183, bordered | ⟨3⟩ | 0, 0, 2, 0 | 642 | 393 (61%) |
| 183 | ⟨34⟩ | 17, 19, 1, 9 | 9259 | 0 |

Part of this is arithmetic. At N = 203 = 7 · 29 the values of A + C on the merge classes span a lattice of
determinant 7168 = 7 · 4⁵, of index 7 in the per-class lattice (step 4 in each class): the compression to ℤ₇ ties
the classes together, and 8574 of the 10000 split values lie off the coset. The other 1426 are on the coset but
still have no pair on one side.

`lowe_sz.py` samples 2 · 10⁷ random (A, C) pairs and 2 · 10⁷ random (B, D) pairs, keeps the merge-class values that
occur on both sides (`live_splits`) and orders only those by the Gaussian law. At N = 215 there are 20000 such values,
and the one holding the 200 solutions of order 430 ranks 6th (2967th in the old order). At N = 203 none occurs, so the
overlap Σ_σ p_AC(σ) p_BD(σ) of the two sides' distributions is below about 3/M² ≈ 10^−14 (M = 2 · 10⁷), and the join
is skipped. The 341 and 203 runs above therefore tested nothing. 341 was searched again with live split values on ⟨2⟩
and ⟨95⟩ (§57–59): no solution in 97% and 95% of the collision mass; ⟨194⟩ with rows (19, 21, 11, 21) was not searched
again.

**Exact compression checks** (`search/compress_check.py`). For every prime power d dividing the group order the
four blocks' compressions to ℤ_d are listed with the symmetry the multiplier group forces on them and joined exactly
(level 1 of `cchain4.py`). An empty join rules the structure/pattern out. Neither the Gaussian nor the hybrid
estimate sees these obstructions. Checked on exact solution sets:
- N = 93, ⟨4⟩. Rows (9, 11, 1, 13) and (11, 11, 3, 11) have no solution (complete joins of lists 126, 126,
  145860, 43758 and 126, 126, 48620, 131274) and both fail at d = 3. Rows (9, 13, 1, 11), with 21384 solutions, pass.
- Bordered 227, ⟨49⟩, rows (0, 0, 2, 0) (order 454) passes: 8 level-1 solutions at d = 2.
- Bordered 31, ⟨11⟩ with the 2 in A has no solution but passes at d = 2, 3 and 5. That obstruction is of another
  kind.
- All 291 exact counts of the calibration runs (`calib_bsz.py`; 27 structures with solutions): no structure with
  solutions fails. Of the 264 without, levels with up to 2 · 10⁵ states per block rule out 186 of 238 cyclic and 24 of
  26 bordered ones, and levels with up to 3 · 10⁶ states (at most 1.4 GB) 220 of 238 and 24 of 26.

For the bordered array (v = N − 1 even) the level d = 2 has a closed form (`bsz.z2_ok`). Σ PAF = −4 gives Σ_k
|X̂_k(ω)|² = 4v + 4 = 4N at every ω ≠ 0. At the character (−1)^t this reads Σ_k e_k² = N, where 2e_k is block k's
alternating sum. The multiplier group fixes the e_k modulo a small number.

At N = 165 (v = 4 · 41, H = ⟨37⟩: trivial modulo 4, order 5 modulo 41) the orbits force:
- e_A, e_B ≡ ±2 (mod 20) for the symmetric blocks of row sum 0;
- e_C ≡ ±1 (mod 10);
- e_D ≡ 0, ±2 (mod 10).

The sum of the squares then lies in {9, 13, 29, 33, 49, 53, 69, 73} modulo 80, and 165 ≡ 5 (mod 80) is not in that
set.

Run on the 24 queue entries of the hybrid ranking (§51), the check removed 15:
- ruled out at d = 2 (bordered): 165 ⟨37⟩ with the 2 in A or C; 269 ⟨9⟩ with the 2 in A, B or C; 223 ⟨7⟩ with the 2 in
  A or C.
- ruled out at d = 25: every 275 pattern in the queue (⟨4⟩: 15, 15, 5, 25 and 5, 15, 15, 25; ⟨14⟩: 15, 25, 5, 15 and
  5, 15, 15, 25; ⟨9⟩: 15, 25, 5, 15; ⟨104⟩: 5, 25, 15, 15). 275 ⟨4⟩ (15, 15, 5, 25) had already run: 0 solutions at
  retention 0.99.
- ruled out at d = 3: 393 ⟨52⟩, rows (27, 29, 1, 1).
- ruled out at d = 5 with 3 · 10⁶ states: 335 ⟨76⟩ in all three patterns and 395 ⟨21⟩.
- passed: 183 ⟨34⟩ in its three patterns, 341 ⟨194⟩, 171 ⟨11⟩, 213 ⟨37⟩, 203 ⟨36⟩, 453 ⟨4⟩, bordered 183 ⟨3⟩ and
  ⟨87⟩, and 403 ⟨16⟩ with rows (9, 9, 15, 35).

Over all 143 structures of the hybrid ranking the check rules out 111. Among the 32 left, apart from the solved 215
and 227, are 203 ⟨36⟩ (three patterns, hybrid 10^2.9–10^3.4), 453 ⟨4⟩, 403 ⟨16⟩ (two patterns), 219 ⟨82⟩ and ⟨10⟩,
455, 291 and 407.

`lowe_sz.py` runs the check for d ≤ 64 (up to 3 · 10⁶ states) before enumerating and stops on an empty level.

## 55. Fibre filters: the compressions inside the low-energy search

A solution's four compressions to ℤ_d form a level-1 solution (§54), so each block's compression lies in the
projection of the level-1 join. `lowe_sz.py` lists these sets for every prime power d ≤ 32 dividing the group order
whose states fit (`fibre_filters`). `lowe_enum.c` drops every member whose compression is not among them, before
computing its energy: a 64-bit key of the compression, linear in the orbit signs, is updated with each flip and looked
up at the leaves. The Gaussian statistics, the energy quantile and the retention model are taken over the allowed
sub-pool (rejection sampling).

Fractions of each block's pool that are allowed, exactly from the multiplicities of the compression states, and the
fraction of the quadruple space that the level-1 solutions cover:

| N | H | rows | d | allowed (A, B, C, D) | quadruples |
|---:|---|---|---:|---|---:|
| 453 | ⟨4⟩ | 29, 31, 1, 3 | 3 | 0.39, 0.39, 0.10, 0.10 | 1.6 · 10⁻³ |
| 371 | ⟨13⟩ | 3, 3, 25, 29 | 7 | 0.40, 0.40, 0.13, 0.13 | 9.4 · 10⁻⁵ |
| 119 | ⟨8⟩ | 1, 5, 15, 15 | 7 | 0.40, 0.40, 0.026, 0.026 | 2.3 · 10⁻⁵ |
| 267 | ⟨4⟩ | 1, 1, 21, 25 | 3 | 0.45, 0.45, 0.11, 0.33 | 7.5 · 10⁻³ |
| 267 | ⟨4⟩ | 1, 1, 21, 25 | 89 | 0.88, 0.88, 0.78, 0.79 | 2.0 · 10⁻⁶ |
| 203 | ⟨36⟩ | 3, 25, 3, 13 | 29 | 0.97, 0.995, 0.98, 0.74 | 1.2 · 10⁻³ |
| 183, bordered | ⟨3⟩ | 0, 0, 2, 0 | 2 | 0.20, 0.20, 0.28, 0.14 | 1.5 · 10⁻³ |

Solutions kept out of exact solution sets (`retention_test.py --fibre`), full energy:

| N | H | rows | kept per big block | without fibres | with fibres |
|---:|---|---|---:|---:|---:|
| 119 | ⟨8⟩ | 1, 5, 15, 15 | 2000 | 0 of 4704 | 4704 |
| 371 | ⟨13⟩ | 3, 3, 25, 29 | 500 | 24 of 96 | 56 |
| 371 | ⟨13⟩ | 3, 3, 25, 29 | 1000 | 48 of 96 | 96 |
| 93 | ⟨4⟩ | 9, 13, 1, 11 | 5000 | 0 of 21384 | 432 |
| 93 | ⟨4⟩ | 9, 13, 1, 11 | 15000 | 1224 of 21384 | 6912 |
| 93 | ⟨4⟩ | 7, 11, 9, 11 | 5000 | 0 of 9504 | 432 |
| 93 | ⟨4⟩ | 7, 11, 9, 11 | 15000 | 864 of 9504 | 1512 |
| 267 | ⟨4⟩ | 1, 1, 21, 25 | 5 · 10⁴ | 0 of 576 | 0 |

At 119 both levels together leave 1345 members of C and of D (1.3% of 105105), so 2000 keep them all. With fibres
`lowe_sz.py` finds the solutions of 371 at L = 1000 (order 742, verified).

**267, and an energy that failed.** At d = 89 the per-block projections keep 78–88% of each pool, while the join
keeps 2 · 10⁻⁶ of the quadruple space. The information is in the combination, and a per-block filter cannot use it.
The compressions also fix linear functionals of the class values: the PAF of the compression to ℤ_d at j is the sum
of PAF(g) over the shifts g ≡ j (mod d). A trial energy conditioned the Gaussian on these 5–7 directions (of 10–17
classes). It was worse everywhere: 267 kept 0 of 576 even at L = 2 · 10⁵, and 119 without fibres 0 of 4704 at
5000. Those directions carry most of the filter's discriminating power, so the option was dropped. A second trial
used the joint fibres directly. Each C member got a mixture energy: the Gaussian density of t − x, summed over the
level-1 solutions that complete its compression, with the other blocks' means taken per fibre (exact, from the
enumerated pools). At d = 89 this ranks 267's solution members at the 18.6th (C) and 31.8th (D) percentiles, no
better than the energies above. The 576 solutions use 48 members of C and 48 of D but only 8 distinct class-value
vectors in each block. A structure's solutions are a handful of distinct points, and no refinement of a Gaussian
filter makes their survival more than a lottery.

**The exact chain cannot take these structures either.** For 453 the joint state space of the levels 3 and 151
exceeds 5 · 10⁶ (and 2 · 10⁷ ran out of memory). At level 3 alone each of the 3 level-1 branches has fibres of about
1.6 · 10⁷ members, above `cchain4.py`'s listing limit of 3 · 10⁶. The fibre filter reaches the same structure from
the other side: its allowed pools are 4.8 · 10⁷ (C) and 1.6 · 10⁷ (D), so lists of 5 · 10⁷ keep every allowed member,
and the low-energy search becomes a complete search of its lists, limited only by the split values kt4 processes.

**kt4 threads.** kt4 opens two or three OpenMP parallel regions per split value over about 10⁵ pairs; on an
oversubscribed machine (4 cores, load 11) the barriers dominate: on the 453 lists one thread processes 128 split
values per second, two threads 57, and four 4. `lowe_sz.py` therefore runs kt4 with one thread (`--kt4-threads`); the
earlier runs of this section and of §50–§52 covered about 30 times fewer split values than their running time allowed.

**Where a solution's split value lies.** At N = 203 (⟨36⟩, rows 3, 25, 3, 13) the kept lists had no live split value
at all (§54). The filter's centring explains it (simulation on uniform samples, `lowe_sz` filter at 3.6% kept).
- On two classes the target lies 2.66 standard deviations (of the four-block sum) away from Σ μ.
- The filter centres each big block at z_k = t − Σ_{j≠k} μ_j, assigning the whole offset to that block, and it does so
  for C and for D at once. On those classes the kept C and D lie 5.5 of their own standard deviations from where a
  solution's C and D lie.
- The kept A + C and t − B − D sums are 2 pooled standard deviations apart, so random pairs never meet.
- Centring at the solution law instead (c_k = μ_k + Σ_k Σ⁻¹(t − Σ μ), consistent by construction) does not help: the
  separation becomes 2.6. On these classes A and B carry most of the variance, so a solution's split a + c lies in
  the tails of both sides, where sampling never looks.

The split values are therefore also taken from the lattice. The values of A + C on the merge classes lie in a coset of
Λ_A + Λ_C (a Hermite basis of the differences within each list; for 371 it has rank 4 on 5 merge classes). The points
of that coset in a box around the solution law's mean that also lie in the coset of t − B − D are added to the sampled
live values, and all are ordered by the law. kt4 discards a split value without pairs in about a millisecond.

**Other centres for the filters** (`retention_test.py --center`).
- Two alternatives to z_k were tried: splitting the offset t − Σ μ between C and D by their covariances, and the
  solution law's mean c_k. Neither helped on the exact solution sets.
- 163 (prime): 0 of 8352 at 5% kept for all three centres, about 850 at 21%.
- 93: 1224 of 21384 at 15000 with z_k, 0 with the other two.
- 371: 24 and 48 of 96 at 500 and 1000 with z_k, 0 with the others.
- 267: 0 for all three.

z_k ranks members by the likelihood ratio, which is optimal for each block alone. The joint losses come from the
few discrete values of a solution's a + b, not from the centre.

## 56. Order 366 over ℤ₁₈₂, and all levels at once

N = 183 = 3 · 61, bordered array over ℤ₁₈₂ ([complex-hadamard/BORDERED.md](../complex-hadamard/BORDERED.md)),
H = ⟨3⟩ of order 6 (trivial on ℤ₂, order 6 on ℤ₇, order 3 on ℤ₁₃), 36 orbits, rows (0, 0, 2, 0).
- **Check.** It passes the exact compression check at d = 2, 7 and 13 (16, 8 and 60216 level-1 solutions).
- **Estimates.** The hybrid estimate is 10^3.18 solutions; with fibre filters the model's retention at L = 3 · 10⁷ is
  10^−0.13 (program not included).
- **Lists.** The fibres at d = 2, 7 and 13 allow 6.9% of A and B (3880 of 55968), 2.1% of C (2.68 · 10⁷ of
  1.28 · 10⁹, all kept) and 3.6% of D (9.35 · 10⁷ of 2.57 · 10⁹, of which 3 · 10⁷ were kept).
- **Join.** kt4 (one thread) took 32423 split values in the law's order: 12599 sampled live values and the rest
  lattice points. It found 1536 solutions at the 146th, after 1.7 · 10⁹ (A, C) and 7.6 · 10⁹ (B, D) pairs. That split
  value was among the sampled live values, so the lattice points did not matter here.
- **Earlier runs.** Without fibre filters the same structure gave no solution in three runs with four threads: full
  energy (642 split values), its continuation (694) and unit-class energy (3358), in part on dead split values (§54).
  The lists kept 2.3% of C and 1.2% of D.
- **Certificate.** [C366.json](../complex-hadamard/C366.json) is the cyclic quadruple over ℤ₁₈₃ (⟨34⟩, row sums 1, 23;
  9, 11) of §42, checked by `verify/verify_complex.py`. The bordered solutions found here are not stored. Order 366 is
  blank in Seberry–Yamada's Table 11.2; Zaitsev's page has a matrix of it (Concurrent,
  [complex-hadamard/README.md](../complex-hadamard/README.md)).

The hybrid estimate of 10^3.2 for the structure and retention near 1 for the filtered lists put about 10^3
solutions in the lists; 1536 at a single split value is consistent with that. For 183 ⟨87⟩, a different subgroup of
order 6 with the same estimate, no result is recorded.

**All levels at once.** Each level's fibre filter is a necessary condition, so a solution's member of each block must
pass all of them together. When the allowed sub-pool is too small for rejection sampling, `lowe_sz.py` lists it
exactly (`lowe_enum` with the filters and no energy cut). An empty list rules the structure out, even though every
single level passed.
- Bordered 273 (ℤ₂₇₂ = ℤ₁₆ × ℤ₁₇, multiplier groups of order 8), ⟨155⟩ with rows (2, 0, 0, 0): the levels 2, 4, 8 and
  17 each allow 19%, 10%, 6.4% and 9.4% of A, but no member of A passes all four. The per-level rates sampled from
  the pool match the masses from the compression states (0.185 against 0.187 at d = 2).
- The same holds for ⟨155⟩ with rows (0, 0, 2, 0) and ⟨53⟩ with rows (2, 0, 0, 0).
- ⟨53⟩ with rows (0, 0, 2, 0) survives with 2.4% of A and B and 0.17% and 0.25% of C and D (3.2 · 10⁷ and 4.7 · 10⁷
  members): lists of 5 · 10⁷ would be complete.

## 57. Re-ranking the open complex orders with the checks and the fibres

**Screen.** The structures of the frontier map (`frontier_scan.py`, open N up to 581) with a per-class estimate of at
least 10⁻⁶ and pools of at most 10¹¹ make 1416 (structure, pattern) records. Level-1 joins with up to 2 · 10⁵ states
rule out 853 of them, among them all 40 bordered records. Levels with 3 · 10⁶ states rule out 73 of the remaining
563. Bordered records with per-class estimates between 10⁻¹² and 10⁻⁶ add 132, of which 108 fall.

**Settled by complete lists.** When the fibres leave fewer members than L in every block, the lists hold every
allowed member, so an empty join settles the structure.
- Bordered 287 = 7 · 41, ℤ₂₈₆ with H = ⟨3⟩ of order 15. With the 2 in C the lists are 424, 424, 31520 and 191712;
  with the 2 in A they are 180, 324, 577472 and 577472. An exact join (`bnw.join4`) of both gives 0 solutions. The
  patterns with the 2 in B or D are images of these.
- Bordered 273 over ℤ₂₇₂, multiplier groups of order 8. ⟨155⟩ and ⟨19⟩ in both patterns, and ⟨53⟩ and ⟨189⟩ with the 2
  in A, have no member of A that passes all levels (§56). ⟨53⟩ and ⟨189⟩ with the 2 in C survive, with complete lists
  of about 3 · 10⁷ and 5 · 10⁷ for C and D, but their hybrid estimates are 10^−1.8 and 10^−1.9.

**203 ⟨36⟩.** All three patterns ran with fibres (d = 29 only, which allows 74–100% of each pool) and lattice split
values. No random pair met on either side, and the 190512–193167 lattice points within Mahalanobis radius 3.4–4.0
of the solution law gave no solution. With 2–4% of C and D kept, the hybrid estimate of 10^3.4 is not borne out.
That is either a failure of the estimate on these classes or the joint retention of §52; the level that would
decide, d = 7, has too many states to list.

**341 = 11 · 31.** Eight multiplier groups carry per-class and hybrid estimates of 10^1.8–10^3.0, with pools of
2 · 10⁹–9 · 10⁹:
- ⟨2⟩, ⟨35⟩, ⟨95⟩ and ⟨194⟩, with fibres at d = 11;
- ⟨15⟩, ⟨27⟩, ⟨58⟩ and ⟨170⟩, with fibres at d = 11 and 31.

With fibres the expected number of solutions in lists of 3 · 10⁷ is 10^0.2–10^1.1, highest for rows (1, 21, 9, 29).
The seven best were queued; the joins of ⟨2⟩ and ⟨95⟩ found no solution (§59).

**Pools above 10¹¹.** 677 records with pools of 10¹¹–10¹³ (per-class estimate at least 10⁻⁴, or 10⁻¹⁰ for the
bordered array); 199 fall at the cheap levels. For the survivors the exact per-level fibre masses give a cheap
allowed fraction; the hybrid estimate was then computed for the best:
- **Bordered 203 = 7 · 29 over ℤ₂₀₂ = ℤ₂ × ℤ₁₀₁, H = ⟨87⟩ of order 5, rows (2, 0, 0, 0): hybrid 10^6.47.** Only d = 2
  gives a fibre (41% of C and D, 34% of A and B), so C and D keep 3 · 10⁷ of 1.1 · 10¹¹ allowed members each. With
  the 2 in C the structure fails at d = 2. The listing takes 2.8 · 10¹¹ leaves per big block; the result is in §58: no
  solution in the part of the join completed.
- 463 (prime), ⟨15⟩, rows (21, 21, 21, 23): hybrid 10^4.37, but C and D have 5.1 · 10¹¹ members each and no fibre.
- 165 = 3 · 5 · 11, ⟨16⟩: hybrid 10^−0.3–10^0.5. Fibres allow 0.8–2.6% of D but 5–36% of C.
- Bordered 235 over ℤ₂₃₄: ⟨185⟩ and ⟨29⟩ (per-class 10^−3.1 and 10^−3.4) fail at d = 3 in every pattern. ⟨43⟩ and
  ⟨121⟩ pass, with hybrid 10^0.4–10^1.6, but about 8% of 1.6 · 10¹² members stay allowed in C and D.
- Cyclic 273 (pools 10¹²) and 455 (hybrid at most 10^−1.6) are out of reach.

**Pruning the enumeration.** Each orbit fills one or two coordinates of a block's compression to ℤ_d (two when ±r
lie in different H-orbits of ℤ_d, for the symmetric blocks).
- `lowe_sz.py` orders the orbits so that whole coordinates are assigned early, most restrictive level first.
- At each depth where a coordinate becomes final, `lowe_enum` compares the partial compression with the allowed
  ones (64-bit keys of the finished coordinates) and skips the subtree on a miss.
- On eleven blocks of eight structures the kept members and their class values are identical with and without the
  checks.
- Leaves visited: 287 bordered D 2.3 · 10⁵ against 1.5 · 10⁷; 273 bordered A 2824 against 104676; 183 bordered C
  1.6 · 10⁸ against 1.28 · 10⁹ (4.3 times faster); 341 A 19305 against 51480; 219 C 1.2 · 10⁷ against 6.0 · 10⁷.
- For bordered 203 ⟨87⟩ the only level is d = 2. The even orbits come first and the check at depth 21 keeps 4
  values of the ℤ₂-compression, about 41% of the leaves.

**Runs** (fibres, live and lattice split values, one kt4 thread):
- 403 ⟨16⟩, rows (9, 9, 15, 35). The fibres at d = 13 leave 4.7% of C and 5.6% of D, so C and D are complete
  (2.04 · 10⁷ and 2.31 · 10⁷). 1.3 · 10⁶ split values (all 575692 sampled live values and lattice points) give no
  solution; the hybrid estimate was 10^0.72.
- 219 ⟨10⟩, rows (11, 25, 7, 9) and (9, 11, 7, 25). The fibres at d = 3 leave C and D complete (1.2 · 10⁷ and
  6 · 10⁶). The first 3 · 10⁵ of about 1.9 · 10⁶ live split values in the law's order give no solution; the hybrid
  estimate was 10^−1.
- Bordered 161 (ℤ₁₆₀): all four multiplier groups of order 4 (per-class estimates 10^−0.7–10^−1.5) fail at d = 2 in
  every pattern, so order 322 stays out of reach. The map has no structure for N = 179 in its orbit range.

## 58. Resumable joins; N = 453 (order 906): no solution in 91% of the collision mass

**Resumable runs.** `lowe_sz.py --budget SECONDS --auto` stops by itself before the budget ends (exit code 75), and
the same command continues:
- Big blocks (pools above 10⁸) are listed in 64 chunks of lowe_enum's parallel prefixes (`lowe_enum ... LO HI`), each
  kept under `$HADAMARD_WORK/lowe/chunks` with a fingerprint of its input; the lowest-energy members are merged at the
  end.
- The lists are saved by atomic renames, and the join records how many split values it finished, and on which merge
  classes, in `kt4_*.bin.state.json`.
- The continuation starts at the first unfinished chunk or the first unjoined split value. A run killed outright loses
  only its current chunk or split value: the join's count is reconciled with the progress file, which holds one line
  per finished split value.

**Parallel joins.** `--kt4-procs P` runs P single-thread kt4 processes on the interleaved parts rest[p::P] of the
unfinished split values and counts as done the prefix every part has finished; a killed parallel run is folded in from
the parts' progress files when the join continues. On the lists of bordered 203 (63504, 2.78 · 10⁷, 127008, 2.78 · 10⁷
members) four single-thread processes join about three times as many pairs per second as one process with four threads
(91 · 10⁶ against 30.5 · 10⁶). On bordered 183 ⟨3⟩ the parallel join again found the 1536 solutions near split value
140.

**Test.** On bordered 183 ⟨3⟩ (order 366), run from scratch in steps of 150 s:
- The chunked lists are identical to the unchunked run's: D passed 30300152 members at τ = 6.607, as before.
- The join ran over six interruptions and found the 1536 solutions at split value 146.
- With the orbit order and prefix checks of §57, C took 2.7 · 10⁷ leaves instead of 1.28 · 10⁹, and D 9.5 · 10⁷
  instead of 2.57 · 10⁹.

**453 ⟨4⟩, rows (29, 31, 1, 3).** The join on the fibre-complete lists (all 2520 allowed members of A and of B,
4.8 · 10⁷ of C and 1.6 · 10⁷ of D) finished all 4329956 sampled live split values, 6.95 · 10¹⁰ (A, C) and
2.31 · 10¹⁰ (B, D) pairs, without a solution.
- Coverage was measured with fresh independent samples (2 · 10⁷ pairs per side; program not included). The processed
  split values hold 62% of the split values met on both sides and 90.6% of the collision mass Σ_σ n_AC n_BD. The
  heaviest missed split values carry 10–16 sample pairs per side.
- So at most about one solution in ten could have escaped. The hybrid estimate of 10^1.8 (about 60) is contradicted
  by at least 1.5 orders of magnitude.
- Either the estimate fails here, as it did for 203 ⟨36⟩, or the compression to ℤ₁₅₁ (too many states to list)
  rules the structure out.

**Bordered 203 ⟨87⟩, rows (2, 0, 0, 0)** (hybrid 10^6.47, per-class 10^3.61).
- Lists. The d = 2 fibre allows 41% of C and D, and the listing visited 1.13 · 10¹¹ leaves of C over 64 chunks,
  keeping 2.78 · 10⁷ members (τ = 1.926). D has the same pool and filter, so it shares C's list. A and B are complete
  (63504 and 127008).
- Join. kt4 ordered 300000 split values (1845983 sampled live) by the solution law. The first 17505 gave no
  solution. By a fresh sample they hold 36.8% of the collision mass, close to the third of the law that the r⁵
  scaling of a 5-dimensional ellipsoid gives.
- Counting. Every solution appears about 64 times among list members (negations of B, C and D, reflections of C and
  D, C ↔ D), so the in-list estimates mean about 100 distinct solutions (hybrid) or 0.15 (per-class).
  - No solution in a third of the law contradicts the hybrid (probability about 10⁻¹⁸), unless the retention of
    every solution is near 0 (the lottery of §52).
  - It is consistent with the per-class estimate.
- 46729 split values are joined, without a solution. They hold 58.6% of the collision mass and, by the r⁵ scaling,
  about 54% of the law. The remaining 253271 (up to about 92% of the law) would take about ten more hours with four
  processes.
  - With a log-uniform prior between the two estimates, updated on this result, the chance that they hold a
    solution is about 20% (less if the lists hold none at all).
  - The search was stopped there; `runq.py` continues it from split value 46730.
- Split values with the coarse class (shift 101, step 40) at ±40 carry much collision mass but come late in the
  law's order. That class is PAF(101) = 202 − 4|supp W| for the ℤ₂ split X = (X₀, X₁), W = X₀ − X₁, and the split
  value fixes how the 203 positions of Σ_k |supp W_k| = 203 divide between (A, C) and (B, D).

## 59. Testing the estimate: value multiplicity; Zaitsev's matrices; order 682

**Joint compression check.** `compress_check.py --joint`: when every level passes, each pair of coprime levels d₁, d₂
is checked together. A block's states are then its compressions to ℤ_{d₁} and ℤ_{d₂} at once, realised by one member,
and the four blocks are joined on the PAF conditions of both levels. On the 291 structures with exact counts
(levels up to 64, as in §54):
- none of the 27 with solutions is ruled out;
- of the 264 without solutions, 244 fall at single levels and 15 more at a joint pair: twelve at N = 45 (eleven at
  d = 3, 5 and one at d = 5, 9) and three at N = 51 (d = 3, 17);
- five pass everything: bordered 31 ⟨11⟩ with the 2 in A or in B (hybrid estimate 10^4.3), 51 ⟨38⟩ (1, 11, 1, 9),
  143 ⟨3⟩ (3, 11, 9, 19) and 247 ⟨9, 35⟩ (17, 17, 7, 19). Over ℤ₃₀ the next check, all three levels at once, is the
  whole problem.

The calibration stops at N = 779, so how often the joint check bites on large structures is not measured. On the
twelve 341 records whose level d = 31 can be listed (H = ⟨15⟩, ⟨27⟩, ⟨58⟩, ⟨170⟩, three patterns each) it passes,
with 120–168 joint solutions.

**How good is the hybrid estimate?** Over the 195 calibration structures with a hybrid estimate (§51–54):
- With solutions (18): log₁₀(exact) − hybrid has median 0.09 and standard deviation 0.43, all within ±1. The
  per-class estimate is low by a median of 1.3.
- Without solutions: hybrid estimates up to 10^4.3. Of the 40 structures with hybrid ≥ 1, 22 have no solution. The
  single levels rule out 16 of these and the joint check 3 more. So among structures with hybrid ≥ 1 that pass the
  per-level checks, 6 of 24 have no solution, and 3 of 21 after the joint check.

The estimate is right in magnitude when solutions exist, but a structure can have none at any estimate. At small N the
compression checks see most of these zeros; what remains is invisible to every check used here.

**Value multiplicity** (`value_mult.py`). Members of a block with equal class values have the same PAF conditions
and the same energy, so the lists keep or drop them together and every solution comes with all such combinations.
For each block, n Σ_v p_v² (pool size n, p_v the share of class-value vector v) is measured from collisions among
uniform samples. The product over the four blocks predicts the number of solutions per distinct value quadruple:

| structure | predicted | measured |
|---|---:|---:|
| 119 ⟨8⟩ (1, 5, 15, 15) | 203 | 196 |
| 267 ⟨4⟩ (1, 1, 21, 25) | 37.1 | 36.0 |
| 93 ⟨4⟩ (7, 11, 9, 11) | 44.4 | 40.6 |
| 93 ⟨4⟩ (9, 13, 1, 11) | 45.6 | 47.5 |
| 163 ⟨38⟩ (17, 19, 1, 1) | 4.0 | 4.0 |
| 371 ⟨13⟩ (3, 3, 25, 29) | 4.1 | 4.0 |

A and B (symmetric blocks) have multiplicity 1.0 in all six. For C and D the trivial part is 2 (reversal), or 4 when
negation keeps the row sum (row sum 0). Composite N add accidental coincidences: 14 per block at N = 119.

On the big structures of §50–58 (10⁶ samples, C and D):

| structure | hybrid | C, D | distinct C, D values | outcome |
|---|---:|---|---:|---|
| 203 bordered ⟨87⟩ (2, 0, 0, 0) | 10^6.47 | 3.9, 3.9 | 10^5.3 | none in 54% of the law |
| 183 ⟨34⟩ (17, 19, 1, 9) | 10^6.69 | 5.8, 6.3 | 10^5.1 | none |
| 453 ⟨4⟩ (29, 31, 1, 3) | 10^1.79 | 6.0, 6.1 | 10^0.2 | none in 91% of the mass |
| 203 ⟨36⟩ (3, 25, 3, 13) | 10^3.38 | 14.1, 13.9 | 10^1.1 | none |
| 213 ⟨37⟩ (13, 17, 13, 15) | 10^4.82 | 6.1, 6.0 | 10^3.3 | 426 |
| 227 bordered ⟨49⟩ (0, 0, 2, 0) | 10^4.40 | 4.0, 4.0 | 10^3.2 | 454 |
| 171 ⟨11⟩ (1, 11, 11, 21) | 10^2.87 | 2.0, 2.0 | 10^2.3 | 342 |
| 215 ⟨11⟩ (11, 17, 15, 15) | 10^4.09 | 10.1, 10.1 | 10^2.1 | 430 |
| 183 bordered ⟨3⟩ (0, 0, 2, 0) | 10^3.18 | 4.2, 4.1 | 10^2.0 | 366 |
| 341 ⟨2⟩ (1, 21, 9, 29) | 10^2.94 | 1.8, 1.8 | 10^2.4 | open (below) |

The five successes lie between 10^2.0 and 10^3.3 distinct values; the four failures lie outside, on both sides.
- 453 and 203 ⟨36⟩ have one to ten distinct value quadruples in the whole pool, so whether the low-energy lists hold
  one is a lottery (§52), and the raw estimate overstated their chances by the multiplicity.
- 203 bordered and 183 ⟨34⟩ would have 10⁵ distinct quadruples. Their failures are not explained by multiplicity or
  retention; both are most likely structures without solutions that pass every check.

**Concurrent constructions.** A literature search for the complex orders found V. Zaitsev's page
(https://www.cs.uleth.ca/~zaitsevv/quaternary.html, last updated 14 October 2026, unchanged on 22 October). It has
matrices of all 12 orders of the first table of [complex-hadamard/README.md](../complex-hadamard/README.md), found
independently, and of 105 more blank orders of the 1992 table (62 of them by a Paley-type construction). Its matrices
of the 12 orders were checked (K K* = nI). It also counts 1018, 1714 and 1766 as known: Đoković's Williamson-type
matrices of orders 509, 857 and 883 (arXiv:1008.2043, from theorems of Miyamoto and Yamada). With these three and the
page's matrices, `complex_status.py` leaves exactly the page's list of unknown orders below 2000 (115 orders, 113
after orders 1286 and 1718 were constructed here; details in
[complex-hadamard](../complex-hadamard/README.md#status)): 334, 358, 418, 438, 446, 478, 498, 502, …
(`data/complex_zaitsev_2026.txt`).

**The open frontier.** Of the frontier records (§57) on the open orders, every one with a hybrid estimate of at least
1 is at N = 341 = 11 · 31 (order 682): 24 records over eight multiplier groups of order 10, with estimates of
10^2.8–10^3.0 for ⟨2⟩, ⟨35⟩, ⟨95⟩, ⟨194⟩ (fibres at d = 11) and 10^1.7–10^2.0 for ⟨15⟩, ⟨27⟩, ⟨58⟩, ⟨170⟩ (fibres at
d = 11 and 31). The next are 303 ⟨95⟩ (10^0.44, pools of 8 · 10¹⁰), 209 ⟨87⟩ (10^0.3, pools of 10¹¹), 219 (at most
10^−0.5) and 455 (at most 10^−1.6). For 341 ⟨2⟩ the per-class (10^2.84) and hybrid (10^2.94) estimates agree, and its
multiplicity is the trivial 2 per block, so it is the best-founded candidate.

**Order 682: the joins of §57–58 in collision-mass order.** Fresh samples (2 · 10⁷ pairs per side) measured how
much of the collision mass the earlier joins covered: 341 ⟨2⟩, 14196 split values, 44.6%; ⟨95⟩, 8514, 32.6%.
- The Gaussian law misorders the remaining mass: on ⟨2⟩, 40000 more split values in the law's order would reach
  78% of it, and in order of sampled mass 86%. The heaviest split values the law left out all have the lumpy first
  merge class at 2.
- `lowe_sz.py --order mass` (new joins) and `--reorder` (continued joins) order the split values by sampled
  collision mass n_AC n_BD, adding sampled split values missing from the list. On 93 ⟨4⟩ (9, 13, 1, 11) the first
  solutions came at split value 17, against 114 in the law's order.
- Both joins continued in mass order with two kt4 processes each (`runq.py`). ⟨2⟩ reached 126731 split values and ⟨95⟩
  96840, without a solution. By fresh samples (seed 2) these hold 97.3% and 94.8% of the collision mass.
- §57 expected up to about 10 solutions inside these lists (rows (1, 21, 9, 29) were the highest), that is 3–4
  distinct value quadruples after the multiplicity 3.2. If that were right, finding none in 95–97% of the mass
  would have a chance of 2–5% for each structure; an expectation ten times lower would make it likely. So the
  estimate fails here too: either these structures have no solutions, as 203 bordered and 183 ⟨34⟩ most likely,
  or the lists' retention is far below the model's.
- The other six groups of 341 are of the same kind (⟨35⟩ and ⟨194⟩ with the same pools and estimates, the four with
  fibres at d = 31 ten times lower). Nothing else on the open list has an estimate above 10^0.5, so the Szöllősi
  arrays, plain and bordered, have no strong structure left on the open orders. The joins were stopped at the split
  values stated above and can be continued from them (`--resume`).

## 60. The largest remaining symmetric classes: what fails

After the constructions of [NOTES.md](../NOTES.md) §30–§32, 499 census values v remain undecided for symmetric
Hadamard matrices of order 4v. The two largest algebraic classes are the 57 prime powers v ≡ 5 (mod 8) and the 35
values v = q + 2 with q ≡ 1 (mod 8) a prime power. No construction was found for either class or for the composite
cases. "Proved" means a complete CP-SAT run or enumeration. The programs for these searches are not part of this
repository.

- **q ≡ 5 (mod 8), order 4q.** The propus-support V-part of [NOTES.md](../NOTES.md) §30 (symmetric V of order 4n, V² =
  qI, support [[I, J, J, 0], [J, 0, I, J], [J, I, 0, J], [0, J, J, I]]) exists for n = 2, 4 (q = 5, 9) and is
  infeasible for n = 3, 6 and 10 (q = 7, 13, 21): this support needs a Hadamard matrix of order n
  ([NOTES.md](../NOTES.md) §33), and none exists for n = 3, 6 or 10; a complete CP-SAT run or enumeration agrees.
  Cyclotomic propus quadruples (Legendre-type A, D, B constant on quartic classes) exist only for q = 5 and 13 among
  the primes q ≡ 5 (mod 8) below 400, and the quaternary 2-circulant cyclotomic ansatz of order 2q only for q = 5, 13.
  A ℤ₂²-refined "skeleton" (two complementary-support symmetric weighing matrices plus a conference matrix of order
  q + 1, which builds symmetric matrices of orders 20 and 36) has no instance at q = 13 in any of its splits (proved,
  using that every 12 × 8 and 12 × 10 partial Hadamard matrix extends to H₁₂ and that Aut(H₁₂) is 5-transitive on
  columns). The known symmetric matrices of orders 4q in this class show no structure to generalise: 52 is Paley II,
  and the propus matrix of order 116 has an automorphism group of order 2.
- **v = q + 2, q ≡ 1 (mod 8).** In the quaternary template of [NOTES.md](../NOTES.md) §31 the cross-block
  orthogonality forces the core to be a unit times a real matrix, so a real symmetric conference matrix of order
  (v + 1)/2 is necessary and sufficient there; this is Seberry–Balonin's condition, and the complex Yamada–Kiyasu
  matrices with quaternary cores (q = 17, 73, 193 tested) have no anti-automorphism. Only 8 of the 35 values (235,
  523, 731, 955, 1011, 1491, 1603, 2379) have (v − 1)/2 a sum of two squares, so for the other 27 no real conference
  matrix of that order exists. A ℤ_{q+1}-circulant two-orbit family finds symmetric quaternary matrices at n = 5, 7, 9
  by annealing (at n = 9 one that comes from a non-Paley symmetric conference matrix of order 38) and none for n = 15;
  orbit-reduced CP-SAT with multiplier groups proved infeasibility in many small structures, including n = 117
  (v = 235) with four subgroups of order 12.
- **Composite q.** Theorem A's U-part, Yamada–Kiyasu's (M, N) and Seberry–Balonin's V-part are all equivalent to a
  symmetric conference matrix C of order q + 1 with a signed permutation R, R² = −I, RCRᵀ = −C (then C = [[M, N], [N,
  −M]] in a basis (e_x, Re_x)); Paley's matrices have one (x ↦ g/x, g a non-square). The three constructions go
  through for any such C (checked for random R at q = 5 to 41), so composite values need a conference matrix with R.
  Of the values undecided at this point, 39 are of Theorem-A type (v ≡ 1 mod 8 composite and a sum of two squares:
  153, 305, 369, …; 65 is settled by order 260) and 33 of Yamada–Kiyasu type (v = q + 2 with q ≡ 5 mod 8 composite and
  a sum of two squares: 119, 223, 247, …), and the Seberry–Balonin values would follow as well. No non-Paley
  conference matrix with R is known at those orders.
- **Symmetrization of quaternary arrays.** The array [[X, Y], [−Y*, X*]] of a complementary pair (x, y) becomes
  symmetric by a monomial diag(Q₁, Q₂) exactly when XQ₁, X\*Q₂ are symmetric and YQ₂ = −Q₁ᵀȲ; this gives the propus
  type (x symmetric, y real, which settled 260: [NOTES.md](../NOTES.md) §32), Hermitian or skew-Hermitian y with
  arbitrary x, and α-twisted pairs for an involutive automorphism α. All were checked exactly on small groups; no
  published pair of a remaining census order is of these types beyond Verma's.
- **Zaitsev's matrices.** Of the 27 quaternary matrices on Zaitsev's page whose order 2v has v among the remaining
  census values, nauty (side-swapping automorphisms of the phase graph, then the column phases) finds no symmetric
  form for QH130, QH214, QH238, QH326, QH426, QH470, QH490, QH494, QH518, QH534, QH538, QH566, QH610, QH686, QH754,
  QH806, QH890, QH1142, QH1462 and QH1526. QH538 and QH1018 have 16 side-swaps, none of which gives a symmetric matrix
  (QH1018 was checked without conjugation only); QH1306 and QH1466 timed out at 25 minutes, and QH1558, QH1690, QH1706
  and QH1882 were not tested ([NOTES.md](../NOTES.md) §31).

## 61. Symmetric conference matrices of orders 118 and 154: block-circulant and block-permuting searches

Brouwer's table of strongly regular graphs (read 17 October 2026) lists conference graphs on 117 and 153 vertices as
existing, credited to J. Y. Z. Ho and to C. Wang respectively as private communications of October 2026, and Verma has
a symmetric conference matrix of order 86 (arXiv:2610.08893: a 12 × 12 array of circulants of order 7 with two border
rows, automorphism group C₇ ⋊ C₃ that fixes six blocks and permutes the other six in two 3-cycles). A symmetric
conference matrix of order 118 would give, by Seberry–Balonin's Corollary 1 with q = 233, the symmetric order 940
(v = 235) and, by the Miyamoto blocks of [hadamard/MIYAMOTO_BLOCKS.md](../hadamard/MIYAMOTO_BLOCKS.md),
[FRONTIER.md](../FRONTIER.md) t = 2823 (order 11292); with an R-structure ([NOTES.md](../NOTES.md) §31) it would also
give v = 119 by the Yamada–Kiyasu route (skew input of order 60). None of these matrices is public, so searches for
such matrices were run:

- **Block-circulant search** (ℤ_m acting with k regular orbits and F fixed points, every circulant block constant on
  the orbits of a multiplier group H, row-by-row meet in the middle on the exact correlation conditions, then the
  fixed rows; program not included). Complete, with no solution: order 118 with m = 29, k = 4, F = 2 and |H| = 4, 7,
  14; m = 13, k = 9, F = 1 and |H| = 4, 6; m = 39, k = 3, F = 1 and H = ⟨4⟩, ⟨10⟩, ⟨5⟩, ⟨14, 22⟩. Order 154 with
  m = 17, k = 9, F = 1 and |H| = 8; m = 19, k = 8, F = 2 and |H| = 6, 9; m = 51, k = 3, F = 1 and H = ⟨2⟩.
  Calibration: the tool finds Paley's matrices (orders 6, 10, 14, 30) at once, but neither order 66 (m = 13, k = 5,
  F = 1, |H| = 2, 3, 4, 6) nor Mathon's order 46 (m = 15, k = 3, F = 1, |H| = 2, 4; m = 9, k = 5, |H| = 2) has a
  structure of this kind, so the 2026 matrices need either the block-permuting groups of Verma's type or no
  multipliers at all. Smaller multiplier groups (|H| ≤ 3 on ℤ₁₃ with k = 9 and on ℤ₃₉ with k = 3, |H| ≤ 4 on ℤ₁₇, ℤ₁₉
  and ℤ₅₁) were not completed (time or memory limits).
- **Block-permuting groups** (G = ⟨ℤ_m, (σ, μ)⟩ with σ permuting the blocks, brute force when the free blocks have at
  most 18 sign bits; program not included): order 118, m = 29, k = 4, μ = 2 and σ a 4-cycle (12 bits) or two
  transpositions (14 bits), and order 66, m = 13, k = 5, μ = 7 and σ a 5-cycle or a 4-cycle: none.
- **Goldberg's matrices and R.** Theorem A for v = 1225 = 35² (the only census value v ≡ 1 mod 8 with a known
  conference matrix of order v + 1) needs an R-structure on Goldberg's matrix of order 1226 built from a skew-Hadamard
  matrix of order 36. On the small members the R-structure fails: the matrix of order 10 (from skew H(4)) has one, the
  one of order 50 (from Paley's skew H(8)) has anti-automorphisms but no R with R² = −I (its whole signed
  (anti-)automorphism group, of order 14 112, was scanned), and the one of order 226 (from the doubled skew H(8)) has
  no anti-automorphism at all (nauty). For 1226 nauty did not finish within 25 minutes.

## 62. Bordered propus quadruples and quaternary Legendre pairs of even length; automorphisms of conference matrices

Searches for the smallest open symmetric values and for the class q ≡ 5 (mod 8) through quaternary pairs found
nothing. The programs for these searches are not part of this repository. A *bordered propus quadruple over ℤ_L* below
is (p, A, A, q) with p, q symmetric, Σ PAF = −4 off 0 and rows (±2, 0, 0) or (0, 0, ±2), equivalently a propus-type
quaternary Legendre pair (A real, B symmetric) of length L ([NOTES.md](../NOTES.md) §32).

- **Exhaustive census for L ≤ 32** (meet in the middle of balanced A against symmetric (p, q), every candidate
  rechecked; L = 20 and 12 recounted independently). Normalised counts (p₀ = q₀ = A₀ = 1, unordered {p, q}): L = 4: 4;
  8: 64; 12: 48; 16: 0; 20: 960; 24: 2688; 28: 0; 32: 18432. A second, independent exhaustive search over all even
  L ≤ 32 (A up to rotation; checked against a brute force at L = 16, 22, 24) finds A-necklaces with a partner for
  L = 8, 10, 12, 14, 18, 20, 24, 26, 30, 32 and none for L = 16, 22, 28. So there is no such quadruple over ℤ₁₆ or
  ℤ₂₈, and hence no uniform algebraic family of length q − 1 for q ≡ 1 (mod 4): q = 17 and 29 already fail (29 ≡ 5
  mod 8), and so does every cyclotomic or Sidelnikov decimation of index e of GF(Q) landing on L = 16 or 28 (e = 1, 3,
  4, 5, 6, 7, 10, 12, …). The zeros at L = 16 and 28 come from the antipodal compression to ℤ_{L/2}, which is
  infeasible there and at no other L ≤ 148 with small divisors; the zero at L = 22 is not explained by it (its
  compression is feasible). Solutions at L = 20, 24 are mostly unstructured, as is Verma's pair (arXiv:2609.25573) at
  L = 64.
- **Jedwab–Pender** (arXiv:2408.08472) Theorem 2(ii) (length (Q − 1)/2) and Theorem 3 (length 2p, Q = 2p − 1) are of
  propus type (checked through `propus-bordered` for Q ≤ 101 and p = 3, 7); they reach v = (Q + 1)/2 and v = 2p + 1,
  none of them a remaining census value.
- **Closed character-sum families** on GF(q)* = ℤ_{q−1} (Möbius-type sums λ(x)ψ(x + 1), class tables of order 2 and 4
  in (x, x + 1), A-first joins against all symmetric (p, q)): no quaternary Legendre pair at q = 13 (positive control
  only), 29, 37, 53, 61, 101, 109, 149. Symmetry forces B = const · ψ₄^{±1}(x)χ(x + 1), whose row sum has |Σ B|² of
  the order of q, incompatible with |Σ A|² + |Σ B|² = 2.
- **Cyclotomic orbit structures** on ℤ_{kq} (k ≤ 12, q ≤ 199) and on K × GF(Q)⁺ for v = 101, 163, 251, 325, 485, 487,
  491, 589, 973 (at most 18–22 orbits): nothing beyond v ≤ 21.
- **Hermitian, skew-Hermitian and α-twisted symmetrizations** (the types of §60). The exact monomials are worked out
  (twist (α, s): x(−αg) = x(g), y(αg) = s·ȳ(g), Q = diag(P_α, −sP_α); s = −1 is equivalent to s = +1), with a bordered
  version giving symmetric order 4(l + 1) (checked for Hermitian l = 2, 6, 10, 14 and twisted ℤ₈, ℤ₁₂, ℤ₂₀). Bordered
  Hermitian pairs over groups of order ≡ 0 (mod 4) do not exist for any abelian group of order 4, 8, 12 or 16
  (exhaustive); if that holds for all orders ≡ 0 (mod 4) it excludes every v ≡ 1 (mod 4), in particular all 57 values
  q ≡ 5 (mod 8). At l ≡ 2 (mod 4) the Hermitian counts equal the propus counts. α-twisted bordered pairs exist at
  l ≡ 0 (mod 4) (ℤ₈, α = 5: 1024; ℤ₁₂: 2304 and 3072; ℤ₂₀, α = 11: 61440). Cyclotomic Hermitian pairs on ℤ_q for the
  18 census primes q ≡ 5 (mod 8) below 800 (e = 4, 12): none (hits only at q = 5, 13, 61); Paley-product twisted
  ansätze on ℤ_m × ℤ_p: no family; multiplier-invariant sweeps at census v ∈ [101, 300] (206 structures): none.
  Lattice-corrected Gaussian estimates leave only v = 101 (bordered ℤ₁₀₀, α = 51: about 2⁴·⁵ expected solutions) and
  v = 119 (unbordered ℤ₁₁₉, α = 69: 2⁴·²) with E > 1, in spaces of 2²⁰⁰–2²⁵⁰ candidates.
- **Unstructured search** at the smallest values. A Markov chain over A inside the spectral region, with a cascade of
  compression tables and an exact lift of the partners (p, q), finds bordered propus quadruples without multiplier
  symmetry at every L = 34 … 52 tried (orders 140 … 212, all decided; 21 of the 22 found at L = 40 … 48 have trivial
  automorphism group, so solutions are plentiful and unstructured), and it recovers Verma's (p, q) from his A at
  L = 64; but its cost grows by about 1.8–1.95 per unit of L (≈ 800 CPU-s per hit at L = 52; Verma's pair was not
  re-found in 1000 CPU-s at L = 64), which puts L = 88, 100, 106 and 118 (v = 89, 101, 107, 119) at 10⁴–10¹⁴
  CPU-years. Simulated annealing at L = 100, 106, 118 and on unbordered propus over ℤ₁₀₁, ℤ₁₀₇ stays far from a
  solution.
- **Propus chain sweep** (§28, §40) over the composite census values from 445: complete with no solution (every
  structure with 1–36 orbits) for 445, 473, 485 and 493 (451 and 473 are settled elsewhere:
  [symmetric/1804](../symmetric/1804) and Đoković, arXiv:2508.17141); 455 and 475 stopped at the 50-minute limit per
  value; the sweep ended before 513.

**Block-permuting groups and an automorphism lemma.** The engine (orbit classes of (i, j, x) under transposition and
h = (σ, μ); exact correlation equations; meet in the middle on the base rows, forward-checking CSP on the σ-fixed
rows, DFS on the moving cycles; program not included) rediscovers Verma's order-86 matrix from his group alone (m = 7,
k = 12, F = 2, μ = 2, σ = (6 7 8)(9 10 11)): once as the negative of the published matrix up to signed permutation,
and in a second run as a matrix equivalent to the published one; its other border pattern is exhaustively empty. With
σ swapping blocks in pairs it finds Paley's matrices of orders 30, 42, 54, 62 and a non-Paley matrix of order 38
(which has an R-structure). It agrees with an independent CP-SAT model on 64 small configurations.

- **Lemma** (proved, and checked on about 7000 automorphisms of matrices of orders 6–86). If π is a permutation
  automorphism of order r of a symmetric conference matrix of order n, then for every d | r the number of cycles of π
  of length divisible by d is even unless √(n − 1) ∈ ℚ(ζ_d) (the trace of C on the ζ-eigenspace is (a − b)√(n − 1) ∈
  ℚ(ζ)). For n = 118 and 154 an involution has f ≡ 2 (mod 4) fixed points and an element of order 3 has f ≡ 4
  (mod 6); for n = 118 order 5 needs f ≡ 8 (mod 10) and order 7 needs f ≡ 6 (mod 14). This excludes every ℤ_m-
  translation structure with 3 | m and one fixed point ((m, k, F) = (39, 3, 1), (9, 13, 1), (3, 39, 1), (117, 1, 1)
  for 118; (51, 3, 1), (9, 17, 1), … for 154), Verma-like layouts with σ of cycle type 1⁶3 or 3³ and μ of order 3, and
  most (σ, μ) pairs (201 of 226 for (118; 13, 9, 1), 23 of 36 for (118; 29, 4, 2)).
- **Exhaustive, no solution:** for n = 118, 35 of the 37 admissible groups with at most 48 free bits in the first row,
  among them the Verma-type group (13, 9, 1), μ = 3, σ = (3 4 5)(6 7 8) (none of its 31 380 canonical first rows has a
  compatible second row; CP-SAT agrees), all groups on (29, 4, 2), (59, 2, 0) and (58, 2, 2); open: (13, 9, 1) with
  μ = 5, σ of type 1·2·2·4 and with μ = 3, σ = id. For n = 154, 144 of 322 admissible groups are settled (no
  solution), 178 not completed. So no conference matrix of order 118 or 154 with this much symmetry (about 20 or more
  bits of reduction) exists in the searched families; Ho's and Wang's matrices, if they are of this kind, use weaker
  symmetry.

## 63. Involutions with 2q − 2 fixed points; bordered propus over the norm-one torus

Two constructions for classes that Theorems A and B do not reach. The first (an involution with 2q − 2 fixed points)
gives a symmetric H(28) and a symmetric W(q + 1, q) for every q ≡ 3 (mod 4), but its Q-part is absent at q = 11 (§66);
the second gives nothing beyond q = 5 for primes q ≤ 47.

**An involution with 2q − 2 fixed points** (`sym_bigf.py`). In Theorem A the involution has 4 fixed points, the field
part is the U-part and the weighing matrix V of order 2(q − 1) is field-free. The roles can be reversed for q ≡ 3
(mod 4). Let H be symmetric of order 4q with a permutation automorphism σ of order 2 that has f fixed points and t
two-cycles. In the basis (fixed points, first points of the pairs, second points), H = [[F, B, B], [Bᵀ, P + V, P − V],
[Bᵀ, P − V, P + V]]. Here V is a symmetric weighing matrix of order t and weight q, and P is a symmetric signed matrix
on the zero set of V. With t = q + 1 the zero set of V is an involution π, P is a signed permutation matrix on it, and
H² = 4qI becomes

    V² = qI,   BᵀB = (2q − 2) I,   F B = −2 B P,   F² + 2 B Bᵀ = 4q I.

Paley's skew core K and R: x ↦ −x, with −1 at ∞, give V = KR, which is symmetric with π fixing 0 and ∞. The Q-part
(F, B, P) is field-free: it depends on q only through its size. Changing the sign of P on a pair is the same as
changing the sign of one column of B, so up to equivalence only the signs at the fixed points of π matter. A search
with B fixed must still enumerate all of them.
- q = 3 (order 12): Q-parts exist for π with 2 fixed points (24 of 48 signed patterns) and with 4 (8 of 16). Only
  the first complete: a symmetric W(4, 3) cannot have zero diagonal.
- q = 7 (order 28): B is 12 × 8 with orthogonal columns. Such a B is equivalent to 8 columns of Paley's H(12). An
  independent check (program not included) classified the k × 12 partial Hadamard matrices exhaustively, confirmed by
  a mass formula against a clique count. The class counts for k = 1, …, 12 are 1, 1, 1, 2, 1, 2, 1, 1, 1, 1, 1, 1, and
  every class with k ≥ 5 completes; this is Verheiden's completion of the last seven rows. The check also enumerated
  the 190080 automorphisms of H(12), whose 95040 column permutations act sharply 5-transitively. With B fixed, each
  row of F lies in the affine set {x ∈ {±1}¹²: Bᵀx = −2(BP)ᵢ}. Exhaustive over every involution π on the 8 columns and
  every sign pattern of P (sign changes of P are column sign changes of B, so with B fixed none can be dropped):
  - Q-parts exist only for π with 0 fixed points (24 of 1680 signed patterns) or 4 fixed points (24 of 13440). There
    are none for π with 2, 6 or 8 fixed points, so Paley's V never fits. `sym_bigf.py` enumerates every sign pattern
    and reproduces these counts.
  - A symmetric W(8, 7) with zero set π exists exactly when π has 0 or 2 fixed points (exhaustive DFS and CP-SAT).
  - Hence the completions are exactly those with π fixed-point-free: symmetric Hadamard matrices of order 28 with an
    involutory permutation automorphism with 12 fixed points exist. An independent check found them by four routes,
    including two representatives in Valtonen's classification data (arXiv:2608.06920, Zenodo
    10.5281/zenodo.22671922), and reached at least 8 of the 91 symmetric classes. It also found order 20 with 8 fixed
    points (q = 5).
  - So the reversed decoupling fails with V = KR (π fixing 0 and ∞), but not with every Paley-type V. For every prime
    power q ≡ 3 (mod 4), V_xy = −χ(1 + xy) (x, y ∈ GF(q)), V_x∞ = χ(x), V_∞∞ = −1 is a symmetric W(q + 1, q) whose
    zero set is the fixed-point-free involution x ↦ −1/x (0 ↔ ∞). This is Paley's core twisted by an element of PSL(2,
    q) with the signs χ(x), which turn it into an anti-automorphism. Checked for q = 3, 7, 11, 19, 23, 31, 43. So for
    q ≡ 3 (mod 4) only the field-free Q-part with π fixed-point-free remains: at q = 7 it exists, F = K(−2P ⊕ 2N)Kᵀ/12
    with K = H(12) and N a symmetric integer 4 × 4 matrix with N² = 7I. Later sections report on larger q.
- The other admissible shapes for q ≡ 3 (mod 4) need t ≡ 0 (mod 4), q + 1 ≤ t ≤ 2q − 2, f = 4q − 2t. A symmetric
  weighing matrix of order t ≡ 2 (mod 4) and weight q needs q to be a sum of two squares. The case f = 0 would be two
  symmetric W(2q, q) with complementary supports, which this rules out for q ≡ 3 (mod 4) prime. Paley gives V only
  for t = q + 1.
- With an order-3 automorphism in addition (the squares of GF(7)* acting on the pairs and every orbit type on the
  fixed points, with π = Paley's and free signs on the pair orbits), CP-SAT finds no solution for q = 7, consistent
  with the above (program not included).

**Bordered propus over the norm-one torus** (`border_torus.py`). §27 builds symmetric H(4v) from a bordered propus
quadruple over any abelian group of order v − 1. For v = q + 2 the norm-one torus T ⊂ GF(q²)*, cyclic of order q + 1,
is the natural group. It acts regularly on the projective line and puts Paley's matrices in two-block circulant
form. Candidate functions: χ(Tr(αt) + c) and η(t − β), with χ and η the quadratic characters of GF(q) and GF(q²),
both values taken at a zero.
- Reductions: η(t − β) = χ(1 + N(β) − Tr(β̄t)), and every χ(Tr(αt) + c) is ± a shift of one with α ∈ {1, z} (z
  primitive), because GF(q)*·T has index 2. A must be symmetric, so all its shifts are tried.
- Triples (A symmetric, B, D up to shift) with |Â|² + 2|B̂|² + |D̂|² = 4v exist only at q = 3 and 5 (80 and 32). There
  are none for q = 7, 11, 13, 17, 19, 23.
- An independent check over the full families, all α, c, β and fills, with or without symmetry of A, also finds none
  for every prime 7 ≤ q ≤ 47. It built the §27 matrices of orders 20 and 28 from all the small triples.
- §27 found the same for character blocks on GF(q)* (v = q), so single quadratic characters do not give bordered
  propus quadruples on either torus.

## 64. Skew analogues of the decoupling

Theorems A and B split a symmetric matrix along a pair-swapping involution. This section looks for skew analogues:
amicable pairs, automorphisms of order 3 and a direct skew Yamada–Kiyasu. Each claim was checked independently. No
census value is settled. The programs for these searches are not part of this repository.

**No involutions.** A skew-type Hadamard matrix H = I + S admits no signed permutation P with P² = I and PHPᵀ = H
other than ±I. A two-cycle (a, b) with sign product +1 would give H_ba = H_ab, and a diagonal P with all signs equal
is ±I. So Aut(H) contains a unique involution, −I, and its Sylow 2-subgroups are cyclic or generalised quaternion.
Examples: Q₁₆ for Paley's H(8), Q₈ for the Goethals–Seidel matrices over ℤ₇, and only {±I} for some Goethals–Seidel
skew H(36) over ℤ₉ (which refutes a claimed Q₈ in every such array). The pair decoupling of Theorems A, B and of
Miyamoto's matrices can never be skew. An automorphism of odd order d cannot decouple supports either: on every cycle
of length l > 1 a ±1 row has an odd number of entries, so its cycle sums never vanish. For an automorphism of odd
prime order p with t cycles of length p, t is even unless p ≡ 3 (mod 4) and 4v − 1 = pk² (character trace argument;
checked by CP-SAT for p = 3, 5, 7 at orders up to 36).

**Amicable pairs.** Skew H(4v) = [[I + S, A], [−A, I + S]] holds iff A is symmetric ±1, S is skew ±1 with zero
diagonal, AS = −SA and A² + SSᵀ = (4v − 1)I (Wallis–Turyn; classical).
- *Circulant normal form.* For 2 × 2 blocks of circulants over an abelian group of odd order v, every pair has the
  form A = [[a, b], [b*, −a]], S = [[s, t], [−t*, −s]] with aa* + bb* + ss* + tt* = 4v − 1 and 2as = bt* − b\*t. The
  traceless form is forced: at each character write the blocks in the Pauli basis. The scalar case would put √(4v − 1)
  or i√(4v − 1) into ℚ(ζ_e), e | v. That is impossible because the odd primes of the square-free part of 4v − 1 ramify
  there.
- *Counts over ℤ_v.* v = 3, 5, 7, 9, 11, 13: 48, 160, 112, 0, 880 and 832 solutions; none over ℤ₃². All counts are
  divisible by 4v, and every solution (2032 in all) builds a skew H(4v), checked exactly.
- *Paley plus a diagonal.* A = C + Δ with C the Paley conference matrix of order q + 1 works only for q = 5 and 9
  among q ≤ 17. The search is exhaustive over the balanced Δ for q = 13 and 17, and Δ must be balanced since tr A = 0.
  At q = 9 the admissible Δ form one orbit of 72: the index-2 cosets of the non-split tori ℤ₁₀. In torus
  coordinates the ansatz is empty for q = 13, 17, 25, 29, 37 and 41. The trivial character needs t₀² = q − 2m₀,
  which fails for q = 73, 89, 97, …
- *Twisted pairs* H = LHR, L ≠ R signed involutions. They decouple H into two eigenspace blocks with MMᵀ = VVᵀ = NI,
  with entries 0, ±1, ±√2, ±2 where fixed points occur. In detail:
  - cycles have even length, holonomy −1, and equal 2-adic valuation within one layout;
  - quartets of half-length w need 2w − 1 to be a square when w is odd;
  - a path of c vertices has a forced diagonal block of squared norm 1/sin²(π/2c), so c < (π/2)√N.
  The layouts with exactly four fixed points on each side are four odd paths of one sign: N = 12, 20, 28 and no
  larger (w = 9 by exhaustive search, w ≥ 11 by the norm bound). Other twisted layouts do exist beyond 28 (12 paths
  of 3 vertices at N = 36). None reaches the census.
- *A two-sector arrangement for v = q + 2.* Coefficient matrices over {I, J, M, N, S}, with Yamada's (M, N) and any
  skew core S, would need W = λXY with λ² = −7/16. That is impossible in ℚ(i), so this fails for every q ≥ 87.
- *Symmetric matrices are not skew in disguise.* The symmetric matrices of Theorem A (q = 9, 17, 25) and Theorem B
  (q = 5, 13, 29), Miyamoto's (q = 5, 13, 17) and Paley II (orders 28, 60) are not monomially equivalent to skew-type
  matrices. Complete searches (program not included) and CP-SAT in independent forms. The control, Paley II of order
  12, is equivalent.

**Order 3** (the Z₃-gluing theorem). Let H be skew-Hadamard of order 4v with a permutation automorphism σ of order 3,
f fixed points and t three-cycles, N = f + 3t. Write Z′_ab for the sum of the circulant block (a, b) and
V_ab = ½ Σ_c h_ab(c) ω^c. Then:
- either |Z′_ab| = 3 and V_ab = 0, or Z′_ab = ±1 and V_ab = −Z′_ab ω^{k_ab}, with k_ab the position of the minority
  entry;
- Q = [[F, √3B], [−√3Bᵀ, Z′]] has Q + Qᵀ = 2I and QQᵀ = 4vI; V has V + V* = I and VV* = vI, so V² = V − vI;
- k_ba = −k_ab for a ≠ b, and conversely F, B ±1 with these relations rebuild H.
- Each row of V has exactly v units, so f ≤ v and f ≡ v (mod 3).
- Let m₊, m₋ be the multiplicities of the eigenvalues (1 ± i√(4v − 1))/2 of V (m₊ + m₋ = t). Since h_aa(0) = 1 and
  h_aa(2) = −h_aa(1), V_aa = (1 + i√3 h_aa(1))/2, so tr V = t/2 + (i√3/2) S with S = Σ_a h_aa(1) over the t cycles;
  comparing with the eigenvalues gives (m₊ − m₋)√(4v − 1) = √3 S. If t is odd then S is odd, so m₊ ≠ m₋ and
  4v − 1 = 3s²: the centred-hexagonal gate of §10 and §14, now for every order-3 automorphism. If 4v − 1 is not of
  that form, then m₊ = m₋, t is even and the diagonal of V is balanced (S = 0).
- For f = 4, F is a skew H(4) and BBᵀ = tI.

So the invariant part fixes the signs of the ω-part. Unlike Theorem A, no input can sit freely in V. Further
results:
- Theorem A's block shape (three full blocks and a diagonal one per block row) forces the constant pairs into four
  cliques K_m.
- The borders cannot be constant on blocks, by a column-sum budget.
- Lemma D: V cannot have the form V_kl = D_kl U_l with diagonal D_kl and one free unitary U_l per block column,
  once m ≥ 3.
- If the class counts R vanish, or the three off-diagonal blocks have rank-one sign patterns, then m ≤ 4.
- Skew matrices of orders 28 and 52 with the four-clique shape exist (CP-SAT plus a SAT lift of the phases; rebuilt
  independently). The shape is specific to Z₃; it is not the support of Theorem A's V, as first claimed.
- Translation-invariant Goethals–Seidel matrices with a multiplier of order 3 (f = 4) exist for q = 7, …, 43, and
  their constant pairs never form cliques.
- The ω-part of a fixed invariant part is not free in practice. Exact counts of completions (raw / up to gauge),
  reproduced by an independent DFS and a direct CP-SAT model:
  - Paley's skew matrix of order p + 1 with its multiplier of order 3: 2/2, 24/8, 774/86 and 3834/142 for
    p = 7, 19, 31 and 43;
  - Paley(27) with the Frobenius map: 816/272, which form 15 Hadamard classes (11 if H ~ Hᵀ);
  - the order-52 clique invariant part: at least 57 up to gauge.
  None of these is parametrised by a Hadamard-type object.

**Direct skew Yamada–Kiyasu.**
- *Literature.* Yamada's J. Austral. Math. Soc. 46 (1989), Theorems 3–5 (Theorem 4 is Kiyasu's) give ordinary
  matrices only. No skew version appears in the sources read: Cati–Pasechnik's implementation of every skew
  construction to order 1208, Momihara–Xiang (arXiv:1801.08776), Karoui (arXiv:2602.16089) and Đoković
  (arXiv:2301.02751). The known cyclotomic skew families give the orders 2(q^e + 1) and 4(q^e + 1).
- *Yamada's trimming.* A skew matrix with a 4 × 4 corner and block-constant borders needs v even.
- *Z[C₃] array.* The pattern with the core in the constant part of one block has a solution at v = 7 (order 28), none
  at v = 13 (three independent proofs, among them an unsatisfiable mod-2 system), and is open at v = 19, 25, where the
  mod-2 system is satisfiable.
- *Quaternary route.* A quaternary Hadamard K of order 2v gives a skew H(4v) by Turyn's map when a monomial W makes
  M = KW twisted-Hermitian, M* = iM + (1 − i)I, and only then. M then has eigenvalues ((1 + i) ± (1 − i)√(4v − 1))/2,
  each of multiplicity v, and det M = (2iv)^v. Hence a necessary invariant: arg det K ∈ (π/2)ℤ. Zaitsev's QH322 and
  QH330 are already in this form up to column phases. Their orders, 644 and 660, were already settled; both skew
  matrices were rebuilt independently. QH214, QH238 and C214, C238 of this repository have none (complete DFS). The
  invariant excludes QH490, QH610, QH890, QH1066 and QH1690 (v = 245, 305, 445, 533, 845). QH382 and QH470 timed out,
  and 18 more at open v were not run; they are scanned in §66.

**Open case.** 106 open skew census values are v = (q + 1)/2 with q ≡ 1 (mod 8) a prime power, i.e. skew H(2(q + 1)),
the classical gap between Szekeres (q ≡ 5 mod 8) and Whiteman (q = p^t, p ≡ 5 mod 8, t ≡ 2 mod 4). §66 treats it.

## 65. Symmetric orders 404 and 596 from a periodic Golay pair (Theorem F)

The structure of §63 (H = [[F, B, B], [Bᵀ, P + V, P − V], [Bᵀ, P − V, P + V]]) settles two census values for q ≡ 5
(mod 8): [symmetric/404](../symmetric/404) and [symmetric/596](../symmetric/596), by the construction of
[symmetric/PERIODIC_GOLAY.md](../symmetric/PERIODIC_GOLAY.md) (`sym_involution.py`, verify.py construction
`symmetric-involution`). Both matrices were rebuilt with separate code, and the `verify/verify.py` builder is a third
independent implementation. The programs for the other searches of this section are not part of this repository.

**The construction (Theorem F).** Let q ≡ 1 (mod 4) be a prime power and K = [[P, Q], [Q, P]] a Hadamard matrix of
order q − 1.
- V is Paley's matrix twisted by the anti-automorphism x ↦ ν/x: V_xy = χ(xy − ν), with zero set the fixed-point-free
  involution x ↦ ν/x.
- F = A ⊗ J₂ + Z ⊗ D₂, with A[x, y] = χ(x + y) and Z the permutation x ↦ −x on GF(q)*.
- B = [𝟙 ⊗ (1, 1)ᵀ | χ ⊗ (1, 1)ᵀ | k_j ⊗ (1, −1)ᵀ], and P = [[0, 1], [1, 0]] ⊕ (−S).

The proof uses only A² = qI − 𝟙𝟙ᵀ − χχᵀ, A𝟙 = −χ, Aχ = −𝟙 (from Paley's C₀² = qI), ZKp = KpS and V² = qI. A sharper
form (Theorem F′) needs three objects:
- N, a symmetric W(q + 1, q) whose zero involution has exactly two fixed points (here N = C₀ composed with x ↦ −x,
  with A its restriction to the other points);
- V, any symmetric W(q + 1, q) whose zero set is a fixed-point-free involution;
- the input K.
No relation between N and V is needed, so prime powers are used only to obtain N and V.

**The input and its obstruction.**
- *Lemma H.* Let K be a Hadamard matrix of order 4w, w odd, with RK = KC for a fixed-point-free symmetric signed
  involution R and any symmetric signed involution C. Then w is a sum of two squares. The rational forms
  m(⟨1⟩^{g} ⊕ ⟨2⟩^{p}) and ⟨1⟩^{f} ⊕ ⟨2⟩^{p′} must be isometric, and for fixed-point-free R this reduces to
  ⟨w⟩ⁿ ≅ ⟨1⟩ⁿ, n = 2w. That fails at any prime ≡ 3 (mod 4) dividing w to an odd power, by its Hasse invariant.
- For K = [[P, Q], [Q, P]] the lemma is the classical condition on W(2w, w), with X = (P + Q)/2 and Y = (P − Q)/2
  complementary weighing matrices.
- Of the 57 census prime powers q ≡ 5 (mod 8), 44 fail it.
- Every Paley-type combination of N and V (from C₀ with k = q + 1 fixed points, V₀ with k = 0, V₂ with k = 2, and,
  for square q, the field-automorphism matrices with k = √q + 1) was listed. For non-square q they need either this
  K or a Theorem-A input (a Hadamard matrix of order (q − 1)/2), so for q ≡ 5 (mod 8) there is no way around Lemma H
  with Paley-type parts.
- Escaping it would need a symmetric W(q + 1, q) whose zero involution has 4 or 6 fixed points. None exists at q = 13
  (exact, using the uniqueness of the conference matrix of order 14), and none occurs among Paley-type matrices. Also:
  no symmetric W(q + 1, q) has q − 1 fixed points for q ≡ 1 (mod 4) (Lemma V1; exhaustive for q = 5 and 9, no proof is
  recorded here).
- The zero-involution spectra of V = C₀G (G a signed automorphism) are q = 5: {0, 2, 6}; 9: {0, 2, 4, 10};
  13: {0, 2, 14}; 17: {0, 2, 18}.

**K from periodic Golay pairs.**
- A periodic Golay pair (a, b) over an abelian group Γ of even order n gives K. In [[A, B], [−Bᵀ, Aᵀ]] the
  translation by an element of order 2 is a fixed-point-free involutory automorphism.
- Conversely, every K found whose blocks are group matrices gives a periodic Golay pair over ℤ₂ × G from its top block
  row, with a = (A + C) + (A − C)t and b = (B + D) + (B − D)t.
- Found (exact checks):
  - K for w = 5, 9, 13, 25, 37, 41, 61, 113 by class-constant searches over GF(w), with A, B, C, D ternary and
    constant on cosets of a subgroup of GF(w)* (for w = 37 with H = ⟨10⟩ of order 3, 35328 pairs);
  - K for w = 37 also from Đoković–Kotsireas's published pair of length 74.
- Census lengths: 2w = 50 and 74 settle q = 101, 149. 194 (q = 389) is undecided in Đoković–Kotsireas's list. 386,
  410, 554, 746, 810, 866, 986, 1098, 1154 and 1170 are beyond published searches, and cyclic pairs of length 810 are
  excluded by Arasu–Xiang.
- Exhaustive class-constant negatives, all with no rows:
  - w = 97, 193, 277, 373, 433 and 577 at the class indices e ∈ {2, 3, 4, 6, 12} that divide w − 1 (and e = 8 for 97);
  - multiplier-orbit structures with at most 13 orbits for n = 205, 405, 493, 549 and 585, and 14–15 orbits for
    n = 205 and 405.
- Paley's matrices give no K. Paley I and II have no fixed-point-free involutory automorphism pair except for orders 4
  and 8 (nauty, full signed groups). In particular, for the Paley conference matrices of orders 6, 10, 14, 18, 26 and
  50, no element of PΣL(2, q′) lifts to a signed automorphism R with R² = −I, which would make R ⊗ [[0, 1], [−1, 0]]
  such an involution of Paley II. The check covers the 15600 and 117600 elements for q′ = 25 and 49.
- A quick search found no cyclic pair of length 50 invariant under multiplier subgroups of order 4, 5, 10 or 20.

**Other results for q ≡ 1 (mod 4).**
- The τ-reduction (F = [[A + Z, A − Z], [A − Z, A + Z]] with A of order q − 1) built H(20), H(36) and H(68). It cannot
  reach q ≡ 5 (mod 8), q > 5:
  - Z = I is impossible for q ≡ 1 (mod 4), by a mod-4 integrality argument on (J − I − S)/2;
  - a Hadamard matrix whose columns are ±1 eigenvectors of a signed involution Z needs Z = ±I or Z fixed-point-free,
    and the latter needs a Hadamard matrix of order (q − 1)/2.
- Parity lemma: every automorphism eigenvalue ζ with √q ∉ ℚ(ζ) has even multiplicity on the Q-part. (The integer
  matrix [[F, B], [2Bᵀ, 2P]] squares to 4qI and commutes with the action, so its trace on the ζ-eigenspace, 2√q times
  the difference of the multiplicities of the eigenvalues ±2√q, lies in ℚ(ζ); hence the two multiplicities are equal.)
  In particular no abelian group of order q − 1 acts regularly on the fixed rows.
- At q = 13 there is no Q-part with a compatible V among the F invariant under free two-orbit actions of ℤ₁₂, ℤ₂ × ℤ₆,
  D₁₂ and A₄ (Dic₃ admits no F), nor under any regular action of the 15 groups of order 24.
- At q = 13 a non-commuting symmetric ternary A′ with A′² = 13I − G (G of rank 2) exists, but none of its 64 sign
  patterns gives a frame. Signature (6, 8, 0) is undecided (unknown at 2400 s). The general q = 13 problem is open.

**q ≡ 3 (mod 4).** The V-part exists for every prime power (§63: V_xy = −χ(1 + xy)). The Q-part exists at q = 7 and is
absent at q = 11 for every π, so no family is known for q ≡ 3 (mod 4) (§66).

## 66. Negative results for symmetric and skew structures, with their scope: none settles a census value

None of these searches settled a census value. Each is recorded with its exact scope. "None" means exhaustive within
the stated structure. Partial and timed-out runs are marked. The programs for these searches are not part of this
repository.

**Decoupled symmetric H(4q) with a regular group on the blocks.** This is Theorem A's shape: four fixed points and
four blocks of n = (q − 1)/2 pair points. A group G of order q − 1 acts regularly on each block, and B is constant on
blocks. Only two types occur: Type I (G = ℤ₂ × ℤ_n, V circulant) and Type II (G = ℤ_{2n}, V negacirculant). A mixed
type would give a row of Z of weight ≥ 2n.
- For non-square q and n ≥ 3:
  - bbᵀ = 4I and F² = 4I;
  - every nontrivial character block has trace 0, so up to the stabiliser of b, b is Sylvester's H(4) and F is fixed;
  - at each real character, Ĥ(χ)/2 is a symmetric integer matrix with square qI.
- Solutions that are not Paley II exist:
  - a symmetric H(44) (Type II, |Aut| = 160);
  - five inequivalent H(52) besides Paley II(25): one of Type I (|Aut| = 576) and four of Type II (384, 192, 96, 96).
- Paley II with 2q − 1 a prime power is of Type II (checked for q = 5, 7, 13, 19, 37, 61).
- q = 17: neither type (SAT; Type I with a DRAT certificate).
- q = 29, Type I: impossible by a finite parity argument, confirmed by CP-SAT and a full SAT model.
  - The three real-character matrices must be among 16 symmetric integer 4 × 4 matrices with M² = 29I.
  - These fall into 4 classes mod 2, and no three of them sum to J mod 2 as the cell parities require.
- q = 29, Type II: undecided. 96 level-1 classes survive, 69 of 74 symmetry templates are UNSAT, and the menu of
  quadratic-cyclotomic sequences over GF(29) is UNSAT.
- Not covered: signed actions (negacyclic Z), and invariance of the Q-part alone.

**The structure with 2q − 2 fixed points, q ≡ 3 (mod 4).** §63 left q = 7 as the only Q-part with π fixed-point-free.
- q = 11: there is no Q-part for any involution π, any sign pattern of P and any B.
  - The 20 × 12 matrices B with orthogonal columns were classified into 235 classes. Every level k = 1, …, 20 was
    checked against a mass formula with nauty group orders.
  - Two independent programs (Python and C++) find no F for any class (22478 leaves).
  - So no symmetric H(44) has an involutory permutation automorphism with 20 fixed points.
- Exhaustive counts of sign patterns of P that admit an F, by the number of fixed points of π (patterns / (P, F)
  pairs, over the classified B):

  | q | 0 | 2 | 4 | 6 | 8 | 10 |
  |---:|---|---|---|---|---|---|
  | 3 | none | 24/24 | 8/8 | | | |
  | 5 | 32/64 | 24/48 | 24/48 | 8/16 | | |
  | 7 | 24/96 | none | 24/96 | none | none | |
  | 9 (30 classes of B) | 108/4576 | 700/5184 | 828/3294 | 434/4534 | 120/1322 | 102/438 |
  | 11 (235 classes) | none | none | none | none | none | none |

- The q = 7 Q-part is unique up to equivalence, and its automorphism group has order 4 (with −I).
- *Parity theorem.* For q ≡ 3 (mod 4), no fixed-point-free Q-part is invariant under a group of odd order with as
  many row orbits as column orbits. This rules out the natural GF(q)-equivariant ansatz (ℤ_{(q−1)/2} with four
  regular row orbits) for every q ≡ 3 (mod 4).
- Q-parts with an involutory automorphism of the q = 7 type: none at q = 11, 15 or 19 (exact).
- So the Q-part exists at q = 7 and is absent at q = 11. It is also absent in the restricted forms above (the
  equivariant ansatz for every q ≡ 3 (mod 4); an involutory automorphism of the q = 7 type for q = 11, 15, 19). The
  general problem is undecided for q = 19, 23, …, so no family is known.

**Order-3 automorphisms.** These are symmetric H(4v) with a permutation automorphism of order 3, f fixed points and t
three-cycles. In the decoupled form the invariant part is a Q-part and the ω-part is a Hermitian U.
- Every parity-admissible (v, f) with 7 ≤ v ≤ 13 occurs, with exact matrices (all known orders): (7, 4),
  (9, 0, 3, 6 or 9), (11, 2), (11, 8), (13, 4), (13, 10).
- *Lemma R.* If U is gauge-real, then f ≤ v − √v. So t = v + 1 needs a genuinely complex U once v ≥ 10, as at
  (11, 8) and (13, 10).
- The f = 2 shape invariant under ℤ_q × ℤ₃ (a symmetric Legendre-type pair of length 3q):
  - it exists for q = 3, 5, 7 and 11 (orders 20, 32, 44, 68);
  - none at q = 15 (exhaustive, two exact lifters);
  - at q = 19 only 28 of 82 orbit representatives were completed, and none of those lifts.
- Index-6 cyclotomic Hermitian U-parts over ℤ_q: only the real (Paley) type, for all 59 primes q ≡ 1 (mod 6) up to
  691.
- A rational (Hasse) test of the Q-part form removes no parity-admissible (v, f) for any of the 499 census v undecided
  at the time.

**Skew H(2(q + 1)) for q ≡ 1 (mod 8).** 106 skew census values have v = (q + 1)/2 with q ≡ 1 (mod 8) a prime power.
The two-block families of Szekeres, Whiteman and Momihara–Xiang need q = pⁿ with p ≡ 5 (mod 8), and none of the 106
has that form.
- Szekeres-type pairs (a skew, both a and b unions of cosets of an odd-order subgroup of GF(q)* of index e, in the
  bordered two-block array), e ∈ {8, 16, 24}:
  - in the range searched they exist only for q = p^{2k} with p ≡ 5 (mod 8);
  - exhaustive for every prime power q ≤ 5000, every prime q ≤ 119993 and every composite prime power q ≤ 40000;
  - this covers the census q with v₂(q − 1) = 3 or 4 (51 and 27 values).
- e = 32: none for the 22 prime powers q ≡ 33 (mod 64) up to 5000. These include all 15 census values with
  v₂(q − 1) = 5.
- e = 64 with b = a ∘ μ (the form of every known solution): none for q = 449, 577, 961.
- Torus subgroups for the 8 open composite q (a ninth, 4913, is settled by skew/9828): none.
- For q ≡ 1 (mod 8), the Bose-quotient relative difference set behind Spence's construction has symmetric halves.
  So that route gives no skew first block.
- Not covered: the other values with v₂(q − 1) ≥ 6, larger e, and non-multiplier groups.
- Heuristic: a count makes two-block pairs over GF(q) unlikely for q above about 82 unless they are algebraic, while
  Goethals–Seidel over ℤ_v has a large expected count.

**Twisted-Hermitian quaternary matrices.** §64's quaternary route was applied to all 27 open skew census values that
have a published or repository quaternary Hadamard matrix of order 2v:
- Zaitsev's QH at v = 107, 119, 191, 235, 245, 269, 283, 295, 305, 445, 509, 571, 653, 655, 679, 731, 733, 779, 787,
  823, 853, 857, 883, 941, 943;
- the repository's C214, C238 and C382, and C1618 and C1906 (v = 809, 953).

None is monomially equivalent to a twisted-Hermitian M.
- 24 values were settled by a complete scan of the pairing equation (every c₀ and both parities e₀).
- v = 245, 305 and 445 were settled by invariants: det K must be a unit times (2v)^v, and rank_𝔭 K = v at every prime
  𝔭 above each prime dividing v.
- An independent numpy solver agrees on the eight matrices of order ≤ 590 it was run on.
- The positive controls QH322 and QH330, scrambled, are found.

**K = [[P, Q], [Q, P]] for the eleven remaining q of §65.** Theorem F needs K of order 4w, w = (q − 1)/4 ∈ {97, 193,
205, 277, 373, 405, 433, 493, 549, 577, 585}.
- *Reduction.* Let G be a group of order w. A K whose blocks are 2 × 2 arrays of G-circulants exists iff two ternary
  rows (A, B) with AA* + BB* = w and complementary supports exist. That holds iff a periodic Golay pair over ℤ₂ × G
  exists. For prime w, ℤ_w-invariance forces this 2 × 2 shape.
- *Lemma S.* For a prime p ≡ 1 (mod 4) there are no symmetric ternary A, B over ℤ_p with AA* + BB* = p.
  - Proof. At a nontrivial character, z = Â + iB̂ lies in L = ℚ(ζ_p)⁺(i) and z z̄ = p. Write p = ππ̄ in ℤ[i]. Both
    primes of ℚ(i) above p are totally ramified in L.
  - Modulo them, z reduces to s_A ± ιs_B, where s_A, s_B are the sums with s_A² + s_B² = p. Exactly one of the two
    vanishes, so (z) = πO_L or π̄O_L.
  - Hence z/π or z/π̄ is a unit of absolute value 1 in a CM field, so a root of unity, so in μ₄.
  - So Â is a rational constant c with |c| ∈ {a, b} (p = a² + b²) at every nontrivial character, and
    A = cδ₀ + κJ. That is not ternary for p > 5. p = 5 was checked directly.
  - Exhaustive searches for p ≤ 41 agree. The lemma has not been compared with the literature.
- Consequence of Lemma S: a block constant on the cosets of H ∋ −1 is symmetric. Every subgroup of even order of
  GF(p)* contains −1, so a class-constant K for prime w needs a multiplier group H of odd order. This explains all the
  negatives of §65.
- *Lemma T2.* A block constant on the cosets of H has sum ≡ A(0) ∈ {0, ±1} (mod |H|). So with w = a² + b², both a
  and b are in {0, ±1} mod |H|. For w = 97, 193, 277, 373, 433 and 577 only |H| ∈ {1, 3} survives.
- w = 97 (q = 389), the twisted shape b(x) = (−1)^x a(gx) over ℤ₂ × ℤ₉₇ (a constant on the cosets of H of order 3,
  g = r⁸):
  - exhaustive meet-in-the-middle (3.24·10⁷ table entries, 1.04·10⁹ probes): no solution;
  - the same code finds twisted pairs at w = 5, 13, 17, 29, 37, 41, 61, 101, 113;
  - the untwisted shape with |H| = 3 (2¹³² candidates) is out of reach, and a heuristic count suggests 10⁰ to 10³
    solutions.
- w = 205 (q = 821): none of the 24 multiplier subgroups of ℤ₂₀₅* with at most 25 orbits admits a ternary row, and the
  central class functions of ℤ₄₁ ⋊ ℤ₅ give none.
- Kronecker and quaternion products of smaller inputs never give ternary rows of odd composite length.
- Every known periodic Golay pair of length 2p with multiplier symmetry (Đoković–Kotsireas, lengths 74, 82, 122, 202,
  226) has |H| odd. Periodic Golay pairs of length 194 are undecided. 106 is the smallest length for which the
  existence of periodic Golay pairs is undecided (Lumsden–Kotsireas–Bright, arXiv:2408.15611).

**Symmetric multipliers of the propus array for t = 25 (§49).**
- A direct CP-SAT search over the symmetric parts of ℤ[ℤ₂₅] and ℤ[ℤ₅²] is unknown after 2400 s each.
- No factorised multiplier M = M₂M₁ through ℤ₅ exists. The trivial character forces det(M₁)² c′₁c′₂c′₃ = 250 for the
  intermediate form C′, so c′₁c′₂c′₃ ∈ {10, 250}. For none of these C′ is there an M₁ over ℤ[ℤ₅] with M₁ᵀC′M₁ = 5C
  (program not included).

**Paley's skew set as the repeated propus block, q ≡ 3 (mod 4).** In the second propus form of §17 (only A
symmetric), put B = C = I + Q, with Q Paley's skew core over ℤ_q (q ≡ 3 mod 4 prime). The known matrix with v = 43
has this B.
- BBᵀ = (q + 1)I − J, so the condition becomes PAF_A + PAF_D = 2 off 0 with A symmetric: a D-optimal pair over ℤ_q
  with one symmetric member.
- The row sums satisfy s_A² + s_D² = 4q − 2, so 2q − 1 must be a sum of two squares. That leaves 35 census prime
  powers q ≡ 3 (mod 4): 163, 167, 223, 239, 283, 343, 491, ….
- Exact counts with A on the orbits of ⟨H₀, −1⟩ and D on those of H₀ (|H₀| = m odd; program not included):
  - solutions: q = 7, m = 1: 56; q = 19, m = 1: 912; q = 19, m = 3: 48; q = 31, m = 5: 24; q = 43, m = 7: 24. The
    matrices of orders 76, 124 and 172 were built and checked.
  - none for q = 23 (m = 1, exhaustive), 31 (m = 3), 43 (m = 3), 59 (29), 79 (3, 13), 103 (17), 127 (7, 9),
    131 (13), 139 (23).
- None in the 64 structures with D on at most 19 orbits for the 35 census values.
- For q = 163 and 167, D-optimal designs of orders 326 and 334 are themselves open (§25). So there the route is at
  least as hard as those problems.

**Quadratic characters on the norm-one torus give no periodic Golay pairs.** On the torus T ⊂ GF(p²)* of order p + 1
(p = 17, 73, 89, 97), the sequences χ(Tr(αt) + c) and χ(N(t − β)), with every fill of their zeros, contain no
periodic Golay pair (exact matching of autocorrelations). §63 found the same for bordered propus. So the torus gives
no K for Theorem F through a prime q′ = 2w − 1, although q′ is prime for w = 97, 205, 405, 549 and 577.

## 67. Theorem G (four fixed points, q ≡ 3 mod 4), symmetric orthogonal designs, PGP(194), composite w for Theorem F

No census value is settled here.

**Involutions of the symmetric H(28).** Valtonen's 91 classes, expanded over their polarities, give 168 conjugacy
classes of symmetric matrices (1958 polarities, matching his file).
- Permutation involutions have f = 4 (24 classes) or f = 12 (24 classes), never both in one class.
- Elements of order 3 always act as 1⁴3⁸. Orders 7 and 13 occur only in a few classes (Paley II(13) and the symmetric
  form of Paley I(27)).
- No element of projective order 4, 5, 9 or 21 occurs, and no signed automorphism squares to −I.

**Theorem G (four fixed points).** Let n be odd, q = 2n + 1, and Q the core of a conference matrix of order n + 1 (Q
symmetric if n ≡ 1, skew if n ≡ 3 (mod 4)). Let A, B be commuting symmetric matrices of order n, A with entries ±1, B
with zero diagonal and ±1 elsewhere, and A² + B² = (q − 2)I. Then
- P = [[0, a₁Q, eI, a₂Q], [a₁Qᵀ, d₁I, a₃Q, 0], [eI, a₃Qᵀ, 0, a₄Q], [a₂Qᵀ, 0, a₄Qᵀ, d₃I]] with fixed signs has P² =
  qI − 2 I₄ ⊗ J;
- V = [[A, I, B, I], [I, −B, I, −A], [B, I, −A, −I], [I, −A, −I, B]] has V² = qI and the complementary support;
- with B₀ a Hadamard matrix of order 4 and F = −½B₀PₐB₀ᵀ, the decoupling of §63 (f = 4, t = 2q − 2) gives a symmetric
  Hadamard matrix of order 4q.

Facts about it:
- Built and checked exactly for q = 7, 11, 15, 19, 27, 31, 39, 63, 75 (the matrices were re-checked). The "sporadic"
  decoupled H(44) of §66 is the case q = 11.
- (A, B) is a symmetric conference matrix of order q − 1 with a skew anti-automorphism (Yamada's circulant pairs when
  q − 2 is a prime ≡ 1 (mod 4)). So for n ≡ 3 (mod 4) the inputs are those of the symmetric Yamada–Kiyasu construction
  ([symmetric/YAMADA_KIYASU.md](../symmetric/YAMADA_KIYASU.md)), and for n ≡ 1 (mod 4) those of Seberry–Balonin's
  conference-core corollary.
- It reaches no census value: the 40 undecided census values v ≡ 3 (mod 4) with v − 2 a prime power and no symmetric
  matrix in this repository all lack the core of order (v + 1)/2.
- Necessity within this skeleton: any V completing such a P contains a conference matrix of order q − 1, so q − 2 must
  be a sum of two squares. With Paley II (2q − 1 a sum of two squares), 39 of the 97 open prime powers q ≡ 3 (mod 4)
  fail both: 107, 179, 251, 311, 347, 359, 419, 599, ….
- Other proved obstructions:
  - A bipartite P (for example the Jacobsthal core of GF(q)* doubled) needs a Hadamard matrix of order q − 1, so the
    GF(q)-natural Q-part cannot be completed.
  - A 4 × 4 array of blocks in span{I, J, Q_GF(q)} exists only for q = 3.
  - An array of eight blocks in span{I, J, Q_GF(n)} is impossible for n ≥ 9.
- At q = 23 the ℤ₂₂-regular structure with any of 32 involutory symmetry templates is UNSAT. The generic structure
  there is unknown.
- Symmetric H(4q) with a fixed-point-free anti-automorphism R (R² = −I, RHRᵀ = −H) are exactly the Turyn images of
  complex-symmetric quaternary Hadamard matrices of order 2q. 74 of the 168 classes of order 28 have one.

**Symmetric orthogonal designs OD(4t; t, 2t, t).** These would serve as symmetric plug-ins of §32 for odd t.
- t = 3: none (28 CNFs, each with a verified DRAT proof; already the mod-4 relaxation is UNSAT).
- t = 5: none, by an exhaustive case split over the 256 conjugacy classes of symmetric W(20, 5). 248 cases are
  DRAT-certified; classes 174, 188, 189, 190, 253 and 255 rest on single kissat runs, and 252 and 254 on cube splits
  (classes numbered in the enumeration of the 256 conjugacy classes of symmetric W(20, 5); program not included).
- t = 4: an exact OD(16; 4, 8, 4) was found.
- No rational obstruction exists for (5, 10, 5), so the obstruction is integral.
- Among the other symmetric OD(12; a, b, c), only (9, 2, 1), (8, 2, 2) and (5, 5, 2) exist (two encoders agree; no
  proofs).

**PGP(194), the untwisted |H| = 3 family** (§66). Still unknown.
- The symmetries force the sums (13, 5, 13, −5), and the family has 31 independent equations in 132 bits.
- Kissat (sorting-network encoding) and CP-SAT already fail on the 52-variable analogue (37, 3), which tabu search
  solves in 0.2 s. Tabu search solves analogues up to (113, 7). At w = 97, a 3600 s kissat run, six cube runs and
  3.3·10⁸ tabu steps (best energy 224) found nothing; both methods extrapolate to 2⁶³–2⁶⁷ steps.
- A heuristic count predicts 1 to 6 equivalence classes. The decisive route is an exhaustive row enumeration with
  canonical rotations (about 6.5·10¹² nodes, 8–14 core-days, then the support join). The same scheme would settle
  (73, 3) in about a core-hour.

**Composite w for Theorem F.** No K = [[P, Q], [Q, P]] was found for q = 821, 1621, 1973, 2197 or 2341 (w = 205, 405,
493, 549, 585).
- Scope: every abelian G of order w and every conjugacy class of subgroups H ≤ Aut(G) acting with at most N orbits. In
  that scope no H-invariant ternary row (A, B) has AA* + BB* = w. So no partner rows exist, and no K of the
  group-circulant X/Y shape.
- N per group: ℤ₂₀₅ 33; ℤ₄₀₅ 33; ℤ₂₇ × ℤ₃ × ℤ₅ 26 (and most classes with 27 to 32 orbits); ℤ₉ × ℤ₉ × ℤ₅ 25;
  ℤ₉ × ℤ₃ × ℤ₃ × ℤ₅ 23 (lattice truncated at 36 orbits); ℤ₃⁴ × ℤ₅ 26 (seven overgroup families only); ℤ₄₉₃ 43; ℤ₅₄₉
  and ℤ₃ × ℤ₃ × ℤ₆₁ 36; ℤ₅₈₅ 31; ℤ₃ × ℤ₃ × ℤ₅ × ℤ₁₃ 29.
- Methods, all exact: the trivial-sum test (407 classes), compression onto coprime quotients with an exact join
  (4386 classes) and full row enumeration with exact spectral pruning (1397 classes). The calibrations at w = 9, 25
  and 45 have rows and K in many classes.
- The smallest open classes (ℤ₂₀₅ and ℤ₄₀₅ with 42 orbits, ℤ₅₄₉ with 39, ℤ₄₉₃ with 44) need more than 3·10⁷ rows in
  their first branch. The non-abelian group ℤ₄₁ ⋊ ℤ₅ and twisted shapes over non-cyclic G were not tried.

## 68. Theorem K: a symmetric conference matrix as the outer design; q = 13 without K; chain sweep; other routes

**Theorem K** ([symmetric/CONFERENCE_PLUG.md](../symmetric/CONFERENCE_PLUG.md), `search/sym_conference_plug.py`,
`verify/verify.py` construction `symmetric-conference-plug`). Let C = [[0, eᵀ], [e, Q]] be a symmetric conference
matrix of order n and m = n − 1. Put X = [[J, J], [J, −J]] and Y = [[Q + I, Q − I], [Q − I, −Q − I]]. Then H =
I ⊗ X + C ⊗ Y is a symmetric Hadamard matrix of order 2n(n − 1): X and Y anticommute, and X² + mY² = 2m(m + 1)I.
- Theorem K is Corollary 2 of E. Spence (An infinite family of Williamson matrices, J. Austral. Math. Soc. Ser. A 24
  (1977) 252–256) at r = 1, in another form. The construction here was found independently; whether its matrices are
  equivalent to Spence's was not checked.
- It is Turyn's symmetric image of Kharaghani–Seberry's complex Hadamard matrix I ⊗ J + C ⊗ (Q + iI) (Seberry–Yamada
  1992, Theorem 11.14, Corollary 11.15), which is complex-symmetric.
- It gives symmetric Hadamard matrices for the census values v = q(q + 1)/2 for q = 17, 25, 29, 49, 61 (orders 612,
  1300, 1740, 4900, 7564). These are Known: Spence's Corollary 3 covers every prime power q ≡ 1 (mod 4), and the
  census series omit it. q = 37 gives a second matrix of order 2812.
- In this scheme (commuting design, anticommuting plugs) the plug order 2k must satisfy k ≥ n − 1 unless k = 1
  (Paley II). So nothing below v = n(n − 1)/2 is reachable.
- The dual scheme (anticommuting design, commuting plugs) is the m(m − 1), m(m + 3) plug-in of §31–32.
- The builder in `verify/verify.py` is an independent implementation; for prime q its output has the SHA-256 of the
  first construction, for q = 25 and 49 it differs only in the field encoding.

**q = 13 without K.** The structure of §63 with f = 24 and t = 14 exists at q = 13, where Lemma H forbids Theorem F's
input.
- Classification of the 24 × k matrices B with orthogonal columns: complete without pruning for k ≤ 9 (1, 1, 1, 4, 5,
  31, 274, 7990, 165596 classes). From k = 9 on it uses two proved hereditary conditions on the rows of B (agreement ≤
  12, and a triple-type bound ≤ 6), and the mass formula is exact at every level.
- 2972745 classes at k = 14 pass these conditions, and 2833 pass a clique test on the possible rows of F.
- Exact Q-part search over them (calibrated on the §66 table for q = 5, 7, 9, 11; a second solver for the k = 14 case
  and a separate solver for k = 2 agree): Q-parts exist for k = 14 (V = Paley's C₀, P diagonal) and k = 2 (V = C₀R).
  There are none with k = 0, the exact shape of Theorem F.
- That gives 284 Q-part classes and symmetric H(52) with an involutory automorphism with 24 fixed points.
- They are sporadic: the automorphism groups have order at most 12 (A₄), and there is no field structure. The only
  recurring symmetry is ℤ₃, acting on PG(1, 13) as ⟨3⟩. The ℤ_m ansatz with m = (q − 1)/4 that it suggests for the 44
  census values is untested beyond m = 3.
- The matrices were rebuilt from the construction data by independent code (all 284 classes assemble, with no
  failure), the pruning conditions were re-derived and one B class was re-solved exactly. The pruned classification
  itself was not re-run; it rests on the exact mass-formula logs.
- Q-parts with k = 4, 8, 10 and 12 also exist, but no symmetric W(14, 13) has such a zero involution (exactly k ∈
  {0, 2, 14}, assuming the uniqueness of the conference matrix of order 14).

**Chain sweep over composite v from 495 to 2499.** Scope: freeD propus, −1 ∉ H, at most 40 orbits per block; 9262
structures over the 316 composite census values (6231 cyclic, 69 needing three generators, 2962 over ℤ_p × ℤ_p × ℤ_m).
- 4793 structures have an empty row-sum pool (proved), and 1588 were joined completely with no solution.
- 101 values are excluded in this scope: 97 by pools alone, plus 553, 933, 1363 and 2413.
- Every one of the 100 best-ranked structures (estimate down to 10⁻⁷·¹) was run: 89 completed, 11 are partial.
- The controls 532 and 876 were re-found.
- Deferred, capped or not run: 55 patterns above the work cap (mostly 513 and 763), and 1518 structures with estimates
  ≤ 10⁻⁸·².

**Other routes.** No census value was settled.
- *Skew to symmetric.* A skew H with an involutory signed anti-automorphism T (THᵀT = H) gives the symmetric matrix
  HT. 58 open symmetric v have a known skew matrix of the same order, but none of 29 cyclic Goethals–Seidel skew
  certificates has a T of the form E ⊗ ρ. The skew 3(q + 1) matrices at q = 7, 11 and 19 have none at all.
- *Jedwab–Pender half-length completion.* It reduces to amicable ternary pairs on ℤ_m with autocorrelation sum −2.
  These exist only for prime m ≤ 23; m = 25 is unknown.
- *Cyclotomic plug pairs* for P ⊗ X + N ⊗ Y: none for primes r ≤ 850 with e ≤ 12, or r ≤ 500 with e ≤ 16.
- *Class-constant symmetric propus* over ℤ_m × GF(q), q ≡ 3 (mod 4), is empty at (3, 7), (3, 11), (5, 11) and (5, 19).
  Two-prime generalised cyclotomy of order 2 gives symmetric propus only at (3, 5), and skew propus nowhere among 491
  prime pairs.
- *Twisted Theorem A* (f = 4) needs the same K as Theorem F.
- *Theorem A for composite v* reduces to a ℤ_n pair (u₁, c) with u₁ symmetric. Such pairs exist exactly for the prime
  powers tested (n ≤ 28) and for none of the composite v there. For v = 153 there are none with any odd-order
  multiplier.
- *A symmetric analogue of the Spence / Whiteman / Leung–Momihara H4 families* would reach 28 to 54 open v. Its
  existence is unknown: the published blocks are not of propus form, and the primary papers of these families were not
  read.

## 69. Larger plugs for Theorem K

Theorem K needs symmetric ±1 matrices X and Y of order K = 2k with XY + YX = 0 and X² + mY² = (m + 1)KI. Any such pair
with m = q (Paley) gives v = (q + 1)k/2. 266 open census values have such a representation with k > m (k odd).

**Structure for m < k < 2m (proved).** Let E = Y² − KI.
- Then X² = KI − mE, E has zero diagonal, and its entries lie in {0, ±2}.
- E_ij ≠ 0 exactly when rows i and j have the same parity of weight. So the support of E is two cliques.
- Each clique block of E/2 is a Seidel matrix with largest eigenvalue below 2, hence of the form D(I − J)D. So E =
  2((I − J_k) ⊕ (I − J_k)) after a signed-permutation conjugation.
- On the span of the two clique indicators, Y acts as ±[[1, −1], [−1, −1]] and X as [[ρ, ρ], [ρ, −ρ]] up to signs,
  with ρ² = k(m + 1) − m and ρ odd.
- Hence (k − 1)(m + 1) + 1 must be a square. k = m + 2 is impossible for every m.
- Of the 43 open (m, k) with m < k < 2m, 37 fail this test. The values 403, 611, 817, 1159, 1273, 1349, 1387, 2263 and
  2449 lose every representation in this scheme.

**k = m + 4** passes (ρ = m + 2) and is the only survivor among open values: v = 119 (13, 17), 377 (25, 29), 779
(37, 41), 1325 (49, 53) and 1539 (53, 57). 189 also passes, but it is already settled.
- Here X = [[J − 2M, J − 2I], [J − 2I, −(J − 2M)]] with M an involution on k points with c two-cycles (proved).
- The problem P(k, c) for Y does not depend on m otherwise: Y is the core of a symmetric Hadamard matrix of order
  2k + 2 that anticommutes with Z = [[M, I], [I, −M]].
- c = 0 is impossible: it forces a Hadamard core of order k + 1 ≡ 2 (mod 4).

**First plug pair with k > m** (exhaustive at k = 9; CaDiCaL and kissat agree).
- P(9, 3) has exactly 1536 labelled solutions in 6 classes, and P(9, c) has none for c ≠ 3. This gives a symmetric
  H(108) (m = 5, k = 9, not a census value), checked in plain integers.
- One class has a ℤ₃ symmetry that is translation by GF(3) when M is the Frobenius map of GF(9).

**Not settled.**
- At (9, 13): c = 1 has no solution, and there is none with a nontrivial cyclic symmetry for c = 2 to 5.
- P(17, c) (v = 119) is unknown; a kissat run on c = 7 did not finish. A calibrated heuristic predicts about 2⁻²⁶
  solutions per c at k = 13 and 2⁻⁹⁹ at k = 17.
- Paley-core Y with (semi-)affine Z gives nothing for q = 9, 13, 17.
- A GF(q²)-Frobenius family (k = q², m = q² − 4) would first reach v = 1127 (q = 7, with Mathon's conference matrix of
  order 46), which is not a census value.
- The regime k ≥ 2m, which holds the other 251 open values with a representation, is not analysed. There E may have
  entries ±4 and the clique lemma fails.

## 70. The ℤ_m ansatz beyond q = 13, and symmetric Spence analogues

No census value is settled here.

**ℤ_m-invariant Q-parts with 2q − 2 fixed points** (§68; the ansatz suggested by q = 13).
- *Lemma 1 (proved).* If ℤ_m (m an odd prime) has a regular column orbit, at most 2(q − 1)/(m + 1) rows are fixed. So
  the ansatz with 2m fixed rows needs m ≤ 3. q = 13 is the equality case.
- *Lemma 2 (proved).* If √q ∉ ℚ(ζ_m), the row and column multiplicities of every character differ by an even number.
  For m = (q − 1)/4 prime ≥ 5 this forces the semiregular shape: 8 row orbits, 4 column orbits plus 2 fixed columns.
- *The semiregular shape* (P diagonal or k = 2, exact orbit-reduced solution of the trivial-character system,
  cross-checked by CP-SAT): there are solutions for m = 1, 9 and 13 and none for m = 3, 5, 7, 11, 15, 17, 19, 21.
  Heuristic: the pattern suggests that both m and 4m + 1 must be sums of two squares, which excludes the census class.
- At q = 29, ℤ₇ (two structure cases), ℤ₅ (three cases) and ℤ₁₃ are infeasible. At q = 61, ℤ₃₁ and semiregular ℤ₁₅ are
  infeasible. The ℤ₅ case at q = 61 with 2 fixed columns and 12 column orbits is unknown.
- *The ℤ₃ family with (q − 1)/2 fixed rows* is the only live route. At q = 29 it has 480 solutions of the ω-part
  (inside a BH(7, 6) tensor ansatz). None was completed to a Q-part: about 160 trivial-level completions are
  infeasible and about 240 hit the 25 s limit.
- At q = 13 the ω-part is a completion of the Butson matrix BH(6, 6) with the profile of Haagerup's F₆. It has no
  GF(13) or quartic-character description, and the H(52) has |Aut| = 12 (it is not Paley II).

**Symmetric analogues of Spence, Whiteman and Leung–Momihara type families.**
- The blocks of these families (regenerated for 12 values of q) are never of propus shape: the four autocorrelations
  differ.
- Four distinct commuting symmetric blocks admit no symmetric 4 × 4 array (a Clifford argument, as in §49).
- *Theorem E (proved).* For a Spence q (m even) there is no symmetric analogue with the extremal pattern. One form
  needs a negacirculant Hadamard matrix of order m, and none exists for 3 ≤ m ≤ 26. The other needs m prime and 2m − 1
  a sum of two squares.
- *Theorem W (proved, built for p = 3, 5, 7, 13, 19, 31, 37, 41, 61, 79, 97).* There is a symmetric analogue of
  Whiteman's construction of order 8p + 4. It uses Paley's conference matrix of order 2p and the torus of order p in
  PSL(2, 2p − 1). It reaches no open census value for p ≤ 1249, since those v = 2p + 1 have 2p − 1 not a prime power.
- The class-constant cyclic scan over ℤ_{v−1} for the Spence-type q = 277, 317, 349, 389, 397, 461, 541, 557, 709,
  821, 1021 found nothing (2 to 20 orbits, cost cap 3·10⁹, groups up to order 1548).
- A sporadic symmetric bordered propus over ℤ₂₈ (order 116, a known order) was found.
- The primary papers on the Spence and Whiteman families were not read.

## 71. Split outer designs and plug pairs for Theorem K; the f = 0 decoupling; twisted Hermitian conference matrices

No census value is settled here. The split framework below was not found in the literature searched and works at a
small order (H(140), a known order).

**Plug pairs with a fixed-point-free ℤ_k (k prime).** In the traceless block-circulant form X = [[a, b], [bᵀ, −a]], Y
= [[e, f], [fᵀ, −e]] (a, e symmetric), X² and Y² are α ⊗ I₂ and β ⊗ I₂ for group-ring elements α, β with α + mβ =
2k(m + 1). Off zero both are sums of two periodic autocorrelations, so they are ≡ 2 (mod 4) and nonzero, and |α_t| =
m|β_t| ≥ 2m. This is the analysis of §69 in circulant form.
- At every nontrivial character the 2 × 2 blocks must be traceless. The singular alternative needs √(2k(m + 1)) ∈
  ℚ(ζ_k)⁺, which fails for prime k ≡ 1 (mod 4).
- Exhaustive meet in the middle over all circulants (program not included): k = 9, 11 and 13 with m = 5, and k = 13
  with m = 9, give no solutions. k = m = 5 gives 80, which is Theorem K.
- Cyclotomic version (blocks constant on the orbits of the index-e subgroup, e ≤ 12, e | k − 1; program not included):
  - covers primes k = 5, 13, 17, 29, 37, 41, 47, 53, 61, 71 and 73 (73 only up to e = 9), with m ∈ {5, 9, …, 61},
    m ≤ k;
  - solutions occur only at m = k;
  - a rerun after fixing a hashing bug (−0.0) gives the same counts.
- So the k ≥ 2m regime of §69 has no solution with this much symmetry for these k.

**Three commuting classes.** Take a commuting design I, C₁, C₂: symmetric weighing matrices with zero diagonal,
disjoint supports and Cₜ² = wₜ I. The traces of C₁, C₂ and C₁C₂ vanish, so the four joint eigenspaces have equal
dimension n/4. Hence n ≡ 0 (mod 4) (proved). With even plugs (n ≡ 2 mod 4, the Theorem K setting) the commuting-design
/ anticommuting-plugs scheme therefore has exactly one non-identity class. For n ≡ 4 (mod 8) with odd plugs the
argument says nothing.

**Split outer designs** (proved; built at n = 14).
- *Construction.* Let C = C₁ + C₂ be a symmetric conference matrix of order n with C₁C₂ = −C₂C₁ and Cₜ² = wₜ I
  (w₁ + w₂ = n − 1). Equivalently, x₁C₁ + x₂C₂ is a symmetric orthogonal design OD(n; w₁, w₂) with zero diagonal.
- Let X₀, X₁, X₂ be symmetric ±1 plugs of order 2k with X₀ anticommuting with X₁ and X₂, X₁X₂ = X₂X₁, and
  X₀² + w₁X₁² + w₂X₂² = 2nk I.
- Then H = I ⊗ X₀ + C₁ ⊗ X₁ + C₂ ⊗ X₂ is a symmetric Hadamard matrix of order 2nk. The cross terms are Cₜ ⊗ {X₀, Xₜ}
  and C₁C₂ ⊗ [X₁, X₂].
- X₁ = X₂ is Theorem K. The point is that the rigidity k ≥ n − 1 of Theorem K (§68) no longer applies.
- *Built.* Paley's matrix of order 14 has such a split with (w₁, w₂) = (4, 9), invariant under multiplication by the
  fourth powers of GF(13). In it C₁ = W(8, 4) ⊕ W(6, 4), a signed K₄,₄ on {∞, 0} ∪ C₀ ∪ C₁ (quartic classes) plus an
  octahedron on C₂ ∪ C₃.
- With X₀ = P ⊗ J₅ (P = [[1, 1], [1, −1]]), CP-SAT found X₁, X₂ of order 10 with X₁ ≠ ±X₂. The result is a symmetric
  H(140) with plug order 10 < 2(n − 1) = 26, checked in integers (program and data not included). v = 35 is not a
  census value; this only shows that the framework escapes the rigidity.

**Plug theory for X₀ = P ⊗ J** (proved).
- Every symmetric X anticommuting with P ⊗ J has the form I ⊗ A + P ⊗ B + M ⊗ D + iσ_y ⊗ δ, with M = [[1, −1],
  [−1, −1]], A𝟙 = B𝟙 = δ𝟙 = 0 and D𝟙 = d𝟙. On span{e ⊗ 𝟙} it acts as dM with d ∈ ℤ.
- The solution above has X₁² = I₂ ⊗ (8I + 2J) and X₂² = I₂ ⊗ (12I − 2J). Call a plug with Xₜ² =
  I₂ ⊗ ((2k − βₜ)I + βₜJ) a β-plug.
  - Then βₜ ≡ 2 (mod 4) and βₜ ≥ −2, and the split needs w₁β₁ + w₂β₂ = −2k.
  - With β₂ = −2 and β₁ = 4j + 2 this gives k = w₂ − (2j + 1)w₁.
  - The action on e ⊗ 𝟙 gives the integrality condition d² = k + (β/2)(k − 1). For β = 2 it needs 2k − 1 to be a
    square: k = 5, 13, 25, 41, …. This excludes k = 7, 9 and 11, so CP-SAT's unknown status on P(7, 2, −2) is settled
    negatively.
- β = −2 plugs: Theorem K's Y = P ⊗ Q + M ⊗ I for any symmetric conference core Q of order k. CP-SAT found one at k =
  7 as well, where no symmetric conference matrix of order 8 exists.
- β = 2k plugs: X₁ = M ⊗ J. It anticommutes with P ⊗ J and commutes with every Y = P ⊗ Q + M ⊗ I, so it is free. The
  split must then satisfy n = (k + 1)(w₁ + 1), with w₂ = k(w₁ + 1).
- Commuting pairs P(k, β₁, β₂) by CP-SAT: P(5, 2, −2) and P(5, 10, −2) exist; P(7, 2, −2) is impossible (above).

**Splits: the Pfaffian lemma** (proved; it supersedes an earlier reach list that assumed splits for every pair of sums
of two squares).
- *Lemma.* If n ≡ 2 (mod 4) and C = C₁ + C₂ is a split, then w₁w₂ is a perfect square.
- *Proof.* T = C₁C₂ is an integer matrix with Tᵀ = C₂C₁ = −T and TTᵀ = w₁w₂ I, so det T = Pf(T)² = (w₁w₂)^{n/2}. Since
  n/2 is odd, w₁w₂ is a square.
- *Consequences.*
  - Order 10 has no split at all, and order 18 has none at (4, 13) or (8, 9).
  - The routes for 299 (n = 46 with (16, 29) or (8, 37)), 765 (n = 90 with (4, 85), n = 170 with (16, 153)), 795
    (n = 318), 1677 (n = 258) and 2125 (n = 170 with (72, 97)) are dead for every conference matrix.
  - So are the routes for the generalised scheme X₀ = P ⊗ W below.
- *Survivors in the constant-β scheme* (both weights squares): 475 (n = 38, (1, 36), k = 25, β = 22), 783 (n = 54,
  (4, 49), k = 29, β = 10), 2035 (n = 110, (9, 100), k = 37, β = 14), 2057 (n = 242, (16, 225), k = 17, β = 26) and
  2197 (n = 338, (81, 256), k = 13, β = 6).
- *Relaxed list.* Only the necessary conditions dₜ ≡ k (mod 2), 0 ≤ dₜ ≤ k and w₁d₁² + w₂d₂² = k(n − k) on the
  S-sector, plus the lemma, at Paley orders ≤ 600 and n = 46. 14 open values pass: 247, 475, 575, 783, 851, 969, 1073,
  1125, 1653, 1845, 2035, 2057, 2125, 2197.
- An independent check: every one of the 72 438 admissible non-square pairs with n ≤ 5001 has an explicit local
  obstruction at a prime ≡ 3 (mod 4).

**Splits: structure and searches.**
- *Matching splits* (w₁ = 1, proved). They are exactly C₁ = [[0, I], [I, 0]] and C₂ = [[A, B], [−B, −A]] with A =
  X + Xᵀ and B = X − Xᵀ. Here X is a signed regular tournament weighing matrix: order m = n/2, zero diagonal, XXᵀ =
  ((m − 1)/2) I, and X + Xᵀ ±1 off the diagonal. So m = 2s² + 1.
  - m = 3 works, m = 9 does not (this settles order 18 completely).
  - m = 19 is order 38 (v = 475). There X has no automorphism of prime order 5, 7, 11, 13 or 17, is not circulant, and
    has no ℤ₉ or ℤ₃ with five 3-cycles. ℤ₃ with six or four 3-cycles and the trivial group are unknown.
- *Clifford parity* (proved). For any monomial automorphism group of a split, every irreducible multiplicity of the
  monomial representation is even. This prunes group actions before any search, and it forces the two-circulant form
  for abelian two-orbit groups.
- *The order-14 split is a tetrahedral object.*
  - All solutions found, from 9 saved files and 7 free CP-SAT runs, fall into one equivalence class. A₄ acts as
    vertices + faces + edges (4 + 4 + 6), with K₄,₄ joining vertices to faces and the octahedron on the edges.
  - The sign character on the edge orbit is necessary: brute force over all 4¹⁰ orbit assignments gives 32 splits, all
    equivalent.
  - The Platonic analogues do not exist in the scanned group actions: S₄ at n = 26, A₅ at n = 62, and the second A₄
    family at n = 30.
- *No anticommuting partner for the octahedron.* No 6 × 6 ±1 matrix anticommutes with W(6, 4), so a weight-4 C₁ has at
  most one octahedral component.
- *Exhaustive negatives with prescribed groups* (Paley matrix fixed, subgroup classes of PSL(2, q)):
  - q = 25 (9, 16) and (5, 20): all 20 classes with |G| ≥ 3;
  - q = 29 (4, 25): all classes with |G| ≥ 3;
  - q = 37 (1, 36): all classes with |G| ≥ 3;
  - q = 53 (4, 49): all classes with |G| ≥ 4;
  - q = 41 (16, 25): all classes with |G| ≥ 8;
  - q = 61 (25, 36): all classes with |G| ≥ 10, A₅ included;
  - q = 109 (9, 100): all classes with |G| ≥ 10.
  - Unknown: the small groups at these orders, n = 26 and n = 30 in general, and n = 110 with 11 Paley-10 blocks under
    ℤ₁₁ or ℤ₄.
- So the splits at 6 and 14 still look sporadic.

**Plug pairs.**
- *Found, k = 13.* A commuting pair (β₁, β₂) = (2, −2) with X₂ = Theorem K's Y for the Paley core of GF(13). X₁ is
  invariant under x ↦ 3x.
  - Complete enumeration in that ansatz gives exactly 8 X₁ up to sign (4 classes); all have B = 0.
  - Re-enumerated independently: 16 = 8 × 2 signs.
  - It would serve 299 at n = 46, which the Pfaffian lemma now forbids.
- *Design form of a β-plug* (proved). X = −M ⊗ J + 2Z with Z = [[A, −C], [−Cᵀ, −B]]. Here A, B, C are a-regular 0/1
  matrices, a = (k − d)/2, A² + CCᵀ = sI + (2a − s)J and AC = CB, with s = (2k − β)/4.
  - Since M ⊗ J commutes with every Theorem K plug Y, [X, Y] = 0 iff [Z, Y] = 0.
  - a = 1 gives permutation plugs with β = 2k − 8 for every k. (13, 2) is a 2-(13, 4, 2) design and (13, 6) an STS(13)
    split into two classes.
- *Regular translation groups* (proved, character argument): no pair of the census targets can be invariant under a
  regular abelian group on both halves. All useful symmetry has fixed points.
- *Rank of X₀* (proved). The only rank-2 X₀ is the balanced P ⊗ J; unbalanced [[J_a, J], [J, −J_b]] admits no plugs.
  For odd v, X₀ = R ⊗ J with R Hadamard of order 4 is impossible.
- *Negatives* (CP-SAT, exhaustive in the stated groups):
  - (17, 26, −2) with Paley Y: infeasible with no symmetry assumption.
  - (13, 6, −2): infeasible for every symmetric ansatz tried except trivial symmetry and ℤ₃ with one fixed point.
    Single β = 6 plugs are plentiful (200 classes found); the commuting partner is the obstruction.
  - k = 25, 29, 37: all groups of order ≥ 6 infeasible (Paley Y, and both plugs free); lower symmetry unknown.
  - Single β = 2 plugs at k = 25 are infeasible for ℤ₁₂ and ℤ₆.
- *Generalised X₀ = P ⊗ W* (W² = (k − b₀)I + b₀J, e.g. W = J − 2A for a polarity of PG(2, q)).
  - Plugs again act on span{e ⊗ 𝟙} as −dM with the same d-condition, and the split condition becomes w₁β₁ + w₂β₂ =
    −2b₀.
  - The candidates it adds (1677, 1785, 2353) all fail the Pfaffian lemma.

So the split framework reaches no census value yet. The live cases need both a split with square weights and a plug
pair at low symmetry. 475 needs the tournament W(19, 9) and P(25, 22, −2); 783 needs (4, 49) at n = 54 and
P(29, 10, −2).

**Quaternary decouplings.** No census value.
- *f = 0 decoupling.* [[A, B], [B, A]] is symmetric Hadamard iff (A ± B)/2 are complementary symmetric W(2q, q). This
  is classical in disjoint-weighing form (de Launey 1993).
  - The symmetric case is infeasible at q = 5, by two independent SAT encodings. The non-symmetric form exists there.
  - A symmetric H(36) of this form exists (q = 9, unstructured). q = 13 is unknown.
  - No census route.
- *Twisted Hermitian conference matrices.* L = I ⊗ J + C ⊗ (Q + iI) with C Hermitian over μ₄ becomes complex-symmetric
  after a symmetric monomial twist S with CS = S C̄ (mechanism proved, checked at w = 6).
  - It needs a complex-symmetric plug as well: a complex-symmetric quaternary Hadamard matrix of order w with an
    imaginary diagonal entry.
  - At w = 58 and 70 (v = 1653, 2415) both inputs are unknown.

## 72. Symmetric conference matrices of order 66; two regular two-graphs outside McKay's set

Symmetric conference matrices of order 66 exist: Gritsenko (arXiv:2102.05432, 2021) constructed an
srg(65, 32, 15, 16); the graph printed there was extracted from the paper and checked in exact arithmetic. Brouwer's
table lists (65, 32, 15, 16) as existing ("+") and OEIS A000952 contains 66. The smallest open order was 86
(Goryainov–Haemers–Konstantinova–Li, arXiv:2512.16653) until Verma's matrix (arXiv:2610.08893). Theorem K at order 66
gives v = 2145, which is not a census value.

**Classification.** All matrices are exact-checked, and every solution was re-verified and re-classified by
independent code.
- *Gritsenko's group.* Gritsenko's matrix has signed automorphism group Q₃₂ × ⟨−I⟩, with point orbits 1, 1, 32, 32.
  After switching it is a 4 × 4 array of circulants of order 16 with a 2-point border.
  - The block engine of §62 needed a "twisted" block permutation (h² = g⁸) to see this group.
  - With that twist, the complete enumeration under Q₃₂ gives 4096 matrices in exactly 5 classes up to
    signed-permutation equivalence and negation. Gritsenko's is one of them.
- *A second family.* A pure Verma-type group (m = 8, k = 8, F = 2, μ = 7, σ = (0 1 2 3)(4 5 6 7); the group is C₈ ⋊
  C₄, not Q₃₂) gives 4608 matrices in exactly 2 further classes.
- *Comparison with McKay's database.* The 7 classes have 48 pairwise non-isomorphic descendants srg(65, 32, 15, 16).
  - Brendan McKay's database of srg(65, 32, 15, 16) (his file sr65321516some.g6, not part of this repository; 32
    graphs, marked "maybe incomplete") is exactly the set of descendants of the 5 Q₃₂ classes.
  - The 16 descendants of the 2 C₈ ⋊ C₄ classes are not in it.
  - So two regular two-graphs on 66 points and sixteen srg(65, 32, 15, 16) are not in McKay's set of 32 graphs
    (marked incomplete); no other enumeration was located ([NOVELTY.md](../NOVELTY.md)).
  - Data: `search/data/conference66_a07eee330cbfe97f.txt` and `conference66_8e2d7e6f507ef47a.txt` (one matrix from
    each of the two C₈ ⋊ C₄ classes) and `search/data/srg65_32_15_16_new16.g6` (their 16 descendants); checker
    `search/check_conf66.py`.
- *Calibration.* The engine re-found Verma's order 86 from its group alone, and Paley matrices of orders 30–62. It
  found order-46 matrices in 2 classes with cover automorphism group of order 20; equality with Mathon's matrix was
  not tested.
- *Prime-order automorphisms of srg(65, 32, 15, 16)* (proved). √65 lies in no ℚ(ζ_p), so the number of p-cycles is
  even.
  - Orbit matrices are infeasible for p ∈ {7, 11, 13, 17, 19, 23, 29, 31}, and for p = 5 except (f, c) = (5, 12),
    which was not lifted.
  - p = 2 and p = 3 are undecided.
  - Among the semiregular groups of order 32 swept, only Q₃₂ and C₈ ⋊ C₄ admit such graphs (44 of 49 groups swept to
    completion; 2⁵ excluded by a moment argument).

## 73. Skew split plug-ins, the order-38 split, non-abelian propus, a skew Theorem K

No census value is settled.

**Order 38 (v = 247, 475).**
- *247 is dead in the split framework* (proved, for any symmetric X₀). With the matching split (1, 36) and k = 13, the
  off-diagonal entries of the three squares are even and lie in [−26, 26], with x₀ + x₁ + 36x₂ = 0. So x₂ = 0 and X₂
  would be a symmetric Hadamard matrix of order 26.
- *475 reduces to P(25, 22, −2)* (proved for X₀ = P ⊗ J₂₅), with d = (17, 1).
  - X₁ is a 2-(25, 4, 1) design with a polarity-like pairing: A² + CCᵀ = 7I + J, CᵀC + B² = 7I + J, AC = CB.
  - Its existence is unknown (CaDiCaL, 40 CPU-minutes). The same model finds the k = 13 analogue.
- *The tournament W(19, 9) has no symmetric support* (exhaustive, two independent re-enumerations).
  - The 0/1 support S satisfies SSᵀ = I (mod 2), equivalently S³ = I over GF(2). All scores of out- and
    in-neighbourhoods are then even.
  - Regular tournaments on 19 vertices with this parity and an automorphism of prime order: none for p = 5, 7, 11, 13
    and 17; only QR(19) for p = 19; 28 classes for p = 3.
  - None of the 28 classes (QR(19) and the second doubly regular tournament included) admits a signing
    (CaDiCaL and CP-SAT).
  - A mod-4 lifting test (a 171 × 171 GF(2) system) fails for 60 of 60 random even tournaments.
  - Heuristic: about 2⁻¹³⁰ expected solutions with trivial automorphism group.

**Order-38 projective-plane guess** (program not included). Both known splits (orders 6 and 14) have weights (a²,
(a + 1)²) with n = 2(a² + a + 1). The natural support guess would give conference matrices at n = 146
(a = 8, an open order): C₁ on the non-incident point–line pairs of PG(2, a) (degree a²), C₂ on everything else (degree
(a + 1)²). It is infeasible at a = 2 (n = 14). a = 3 (n = 26) is unknown (undecided after 1200 s). Not pursued.

**Skew split.**
- The skew analogue of §71 is verified. With S₁, …, S_r anticommuting skew weighing matrices partitioning J − I, H is
  skew-Hadamard of order nK when:
  - X₀ is skew-type;
  - the X_t are symmetric and pairwise commuting;
  - X₀X_t is symmetric;
  - X₀X₀ᵀ + Σ w_t X_t² = nK I.
- It is classical: with X₁ = J it is Seberry's 1978 theorem (skew/1900, 2772, …), and the general case is the amicable
  OD plug-in.
- *Built* (known orders, v ≤ 65, exact checks): H(36), H(60), H(84), H(132) and H(260). H(260) uses OD(20; 1, 2, 8, 9)
  with a plug of order 13 that is not of Seberry type.
- *n = 12 table.* Skew zero-diagonal OD(12; ·) exist for (1, 10), (2, 9), (3, 8), (5, 6) and for (1, 1, 9), (1, 2, 8),
  (1, 5, 5), (2, 3, 6). (4, 7) does not exist, by the argument below (a sketch; the condition was not found in the
  literature read): for n ≡ 4 (mod 8) the quaternion algebra (−w₁, −w₂) must be the Hamilton algebra (ramified exactly
  at 2 and infinity), and (−4, −7) is ramified at 7 and infinity.
  - *Sketch of the condition.* For n ≡ 4 (mod 8), S₁ and S₂ make ℚⁿ a free module of odd rank n/4 over the definite
    quaternion algebra D = (−w₁, −w₂), and the standard form is isometric to ⟨a₁, …, a_{n/4}⟩ ⊗ N_D up to a scalar
    factor, with N_D the norm form of D. If D is ramified at an odd prime p, N_D has anisotropic binary residue forms
    at p, and each term contributes the same nonzero Witt class of order 2 to the second residue form. An odd number
    of such classes cannot cancel, whereas the second residue form of n⟨1⟩ is 0.
- *Circulant palindromic plugs.* Over ℤ_K the solutions lie on lines in weight space: (t, 3t − 1) for K = 3, and
  (t − 1, t, 2t) for K = 5 and 7. The classification is exhaustive for n ≡ 4 (mod 8) only.
- *Reach.*
  - The (t − 1, t, 2t) line reaches the open skew values 1141, 1407, 2023 (K = 7) and 1965, 2255 (K = 5). They would
    need skew OD(4t; 1, t − 1, t, 2t) with t = 163 to 451; none is in hand.
  - The Seberry type covers 126 open skew values in principle. It needs large amicable/OD inputs; the smallest cases
    (473, 1139) die by a divisibility lemma.

**Non-abelian groups.**
- *Transplanting the arrays.* The literal transplant of the Goethals–Seidel and propus arrays (R = inversion) works
  only for class-function blocks. For ℤ_p ⋊ ℤ_r and Heis(p) it never gives a Hadamard matrix (Plancherel count).
- *The working version* replaces R by an involutory automorphism φ (blocks constant on ψ(g) = φ(g)⁻¹).
- *Skew:* impossible whenever φ fixes a prime-order quotient. This closes 87 of the 122 open skew values that admit a
  non-abelian group (201, 305 and 385 included).
- *Symmetric:* the φ-twisted propus is not vacuous. It gives symmetric H(84) over ℤ₇ ⋊ ℤ₃ (exactly 24 normalised
  solutions, exhaustive; confirmed by an independent meet in the middle) and H(108) over Heis(3).
  - Both orders are known (Paley II). The construction type was not found in the literature searched.
  - The solutions are sporadic, and far sparser than the abelian scheme: 24 against 1 068 984 at v = 21.
  - Searches with Aut(G)-subgroup symmetry at census sizes are infeasible or unknown.
- 100 of the 492 open symmetric values admit a non-abelian group.

**Skew Theorem K.**
- The real analogue H = I ⊗ X + S ⊗ Y (S skew conference matrix, X skew-type, Y symmetric, odd plug order K) is rigid
  (proved):
  - K = 1 + (n − 1)t with t even and 1 + nt a square;
  - rank(X − I) = K − 1;
  - t = 2 and t = 4 are impossible, so K ≥ 6n − 5.
- 15 open skew values pass the necessary conditions (smallest 575, with K = 115), all far beyond search.
- *Complex route.* No tournament-type complex conference matrix of even order ≥ 4 exists. The Kharaghani–Seberry /
  Theorem K complex matrix is never monomially tournament-type (proved for n ≥ 10). So Theorem K has no skew twist for
  153, 325, 435, 1225 and 1891, which are all open in the skew census.
- The n = 2 case equals "skew H with a fixed-point-free signed automorphism of square −I". Paley's skew matrices with
  q ≡ 3 (mod 8) have one (classical).

## 74. Kharaghani–Seberry plug pairs for q > p, complex split triples, routes to open complex orders

No open complex order is settled in this section (115 open orders 2N, N odd, N < 1000; Zaitsev's page lists the same
115; §78 settles 1718 and §83 settles 1286).

**The setting.** Seberry–Yamada 1992, Theorem 11.14: a skew-type complex Hadamard matrix of order p + 1 and a plug
pair (A, B) of order q with AB* = BA* and AA* + pBB* = q(p + 1)I give a complex Hadamard matrix of order q(p + 1).
- 65 open orders factor as (p + 1)q with a known skew-type order p + 1 ≡ 2 (mod 4) and odd q ≥ p.
- The published tables use only the pair A = J, B = I − iQ (q = p), and its recursion q = p^j.
- Proved: A = J forces q = p.

**Plug pairs.**
- *Rigidity, p < q < √5p* (proved; computer-assisted classification of μ₄-gain matrices with λ_max ≤ q/p, reproduced
  independently).
  - After scaling, AA* = (q − p)I + pJ and BB* = (q + 1)I − J.
  - So pq + q − p must be a square, which is the parameter equation of a symmetric (q, (q − s)/2, ·) design.
  - This kills the factorizations 42·43, 38·43, 38·51, 26·35, 30·51, 30·53, 26·47, 26·49 and 26·53. 1222, 1378 and
    1634 lose every Kharaghani–Seberry route.
- *q = p + 2 is impossible* for every p ≥ 3 (proved twice independently). 1806 = 42·43 is out.
- *Theorem 2* (proved). For q = p + 4, (A, B) = (J − 2P, −NP) with P an involution is a pair iff [[1, 1], [1, N]] is a
  Hermitian quaternary Hadamard matrix of order p + 5.
  - Found at p + 5 = 6, 10 and 14, giving pairs at (5, 9) and (9, 13): orders 54 and 130, both known.
  - (13, 17), (17, 21) and (25, 29) need Hermitian matrices of orders 18, 22 and 30: unknown.
- *Circulant pairs* exist only at q = p (orbit-constant sweeps over all 89 open factorizations).
- *The regime q ≥ √5p* (80 of 89 factorizations) is unclassified. The normal-form candidates 970, 1242, 1386, 1530 and
  1854 need design-like objects such as symmetric (97, 33, 11) or (153, 57, 21) designs with a Hermitian coupling.

**Complex split triples.** G = I ⊗ B₀ + C₁ ⊗ B₁ + C₂ ⊗ B₂ with the real splits of orders 6 and 14.
- *Reduction.* G is unitarily equivalent to I_{n/2} ⊗ M with M = [[B₀ + √w₁B₁, √w₂B₂], [√w₂B₂, B₀ − √w₁B₁]].
- *Calibration.* K = 1 and 5 work: orders 30 and 70. At n = 14, K = 5 is below the Kharaghani–Seberry bound 13. K = 3
  is infeasible for both splits.
- *Lemmas* (proved, re-derived independently):
  - group-developed (constant row sum) plugs need K ≡ 1 (mod 4);
  - the Hermitian-circulant class is impossible.
- *Searches.* Circulant searches up to K = 13 and multiplier scans (K ≤ 257 prime, several composites) found nothing.
- *The 574 candidate.* The regular-Gram candidate at (n, K) = (14, 41) was excluded by a parity argument. So no
  complex split triple is known beyond the Paley sizes.

**Literature and routes.**
- *Literature.* No accessible source uses plug pairs with q ≠ p or a 3-variable complex plug-in for orders ≡ 2
  (mod 4). Kharaghani–Seberry 1990/1993, Wallis 1973 and the 2020 book were inaccessible.
- *Quaternary conference matrices of orders 22 and 34.* Zero diagonal, units elsewhere, CC* = (n − 1)I, non-Hermitian;
  no real conference matrix exists at these orders. They are built from the zero-diagonal Szöllősi array K = [[X, Y],
  [Y*, −X*]] (A, B symmetric circulant with a₀ = b₀ = 0, C, D circulant, AAᵀ + BBᵀ + CCᵀ + DDᵀ = (4m − 2)I). Checked
  exactly by two independent codes; novelty not checked.
- *Complex Yamada–Kiyasu.* Implemented independently; it reproduces Zaitsev's 470, 566 and 1142.
- *Six open orders, one missing core each.* 1046, 1190, 1286, 1718, 1910 and 1958 need a quaternary conference matrix
  of order 262, 298, 322, 430, 478 or 490 respectively. The template side (Yamada circulants for q = 521, 593, 641,
  857, 953, 977) is verified. Cheap structured searches found nothing; this is the best-ranked constructive route. The
  core of order 430 (order 1718) is found in §78, and that of order 322 (order 1286) in §83.
- *Route map.* No classical complex route (Cor. 11.15, 11.20, Kronecker, Golay, Miyamoto, Seberry–Zhang) settles any
  open order.

## 75. Zero-diagonal cores for the complex Yamada–Kiyasu targets

No core was found, so 1046, 1190, 1286, 1718, 1910 and 1958 stay open in this section (§78 finds the core of order
430, hence 1718, and §83 that of order 322, hence 1286). Each needs a quaternary conference matrix of order c = 262,
298, 322, 430, 478 or 490 (§74). Three independent implementations were run, and a calibration of the estimators
followed.

**The array.** The zero-diagonal Szöllősi array of §74 over an abelian group of order m = c/2:
- A, B are symmetric with A(0) = B(0) = 0; C, D are ±1.
- Σ PAF = (4m − 2)δ and a² + b² + c² + d² = 4m − 2.
- Blocks are constant on the orbits of a multiplier group H with −1 ∉ H: A, B on ⟨H, −1⟩-orbits, C, D on H-orbits.

Two independent implementations rebuilt the pipeline end to end: core 142 (m = 71) gives order 566, and core 286
(m = 143, H = ⟨3⟩) gives 1142; both exact.

**Compression filters.** For a quotient ℤ_d, the fibre sums of the four blocks must solve Σ PAF = (4m − 2)δ on ℤ_d.
The feasible tuples are listed exactly, and a lattice-Gaussian estimate is conditioned on each tuple.

**Exhaustive negatives.**
- *m = 245, cyclic ℤ₂₄₅*, all 19 subgroups with −1 ∉ H, even orders included:
  - every |H| ≥ 7 is dead on ℤ₅ or ℤ₇, except ⟨18⟩ (|H| = 12);
  - ⟨18⟩: 9 of 12 row patterns die on ℤ₅, and the other 3 are empty by exact joins over all reduced ℤ₃₅ tuples;
  - for H = ⟨36⟩ (|H| = 7), an independent second implementation finds no solution of the ℤ₇-compressed equation for
    any row sums.
- *m = 245, ℤ₅ × ℤ₇²*:
  - odd H: |H| = 63, all four classes of order 21, and the unipotent |H| = 7 are dead on line, plane or ℤ₅ quotients;
  - even H: every sampled class with |H| ≥ 24 has no live pattern, and the classes nc:2054 (|H| = 8) and nc:6609
    (|H| = 12) are empty by exact joins;
  - caveat: the subgroup list comes from random closure sampling (two samples, no new class in the second); it is not
    a certified class list, and the labels nc:k index it (the list is not part of this repository).
- *m = 215*:
  - H = ⟨11⟩ (|H| = 7): the single pattern (14, 18, 13, 13) has no solution of the ℤ₅-compressed equation (two
    independent computations);
  - ⟨4⟩ and ⟨51⟩ (|H| = 14): complete lists joined, no solution;
  - |H| = 21, 28, 42 and 84: no row pattern.
- *m = 161*: H of order 11, 22, 33 and 66 has no row pattern.
- *Primes m = 131, 149, 239*: no multiplier structure. Odd-order H admits no row pattern, and every even-order
  subgroup of the cyclic ℤ_m* contains −1.

**The bordered zero-diagonal array** (proved, exact at small orders; novelty not checked).
- *Setting.* G abelian of even order n = m − 1; A, B symmetric with zero at 0; C, D ±1; Σ PAF = −4 off 0.
- *Row pattern.* Parity forces the pattern (±1, ±1, 0, 0), so a pattern always exists, unlike the unbordered array.
- *Construction.* Put α = (1 + i)/2, ρ = α r_A + ᾱ r_B, X′ = αX_A + ᾱX_B and Y′ = αX_C + ᾱX_D R. With the border x =
  0, s = ρ²t̄, y = −ρt̄v, u = ρ²t̄²v (t, v units), the matrices X = [[0, s eᵀ], [t e, X′]] and Y = [[y, u eᵀ],
  [v e, Y′]] give a quaternary conference matrix K = [[X, Y], [Y*, −X*]] of order 2m.
- *Checks.* Exact at n = 10, 12, ℤ₂ × ℤ₆, 16 and 26.
- *t-character version.* For t odd it generalizes with a conference matrix of order 2t in the border, checked at t =
  3, order 22.
- *Searches at N = 131* (core 262, order 1046):
  - H = ⟨101⟩ and ⟨9⟩ (|H| = 6) are empty, exhaustively over all ℤ₂ tuples;
  - for H = ⟨99⟩ (|H| = 4), low-energy lists of 10⁷ give 0 hits against about 0.9 expected (not conclusive);
  - H = ⟨21⟩ (|H| = 4, estimate 10^1.5, pools 4 · 10⁹) was not run.
- *N = 149 and 239*: no live, searchable bordered cyclic structure.

**Live structures** (|H| = 6; none searchable exhaustively):
- *m = 215, H = ⟨49⟩*, which acts as −1 on ℤ₅ and with order 3 on ℤ₄₃:
  - 12 row patterns survive ℤ₅;
  - pattern (14, 18, 13, 13) with ℤ₅ tuple ((−6,5,5), (−6,5,7), (1,13,−7), (1,−1,7)) has conditioned estimate 10^4.1,
    pools 4.3 · 10⁴, 4.3 · 10⁴, 2.1 · 10¹⁰ and 3.5 · 10¹⁰.
- *m = 245*: ⟨31⟩, ⟨79⟩ and the non-cyclic twins nc:36 and nc:842 over ℤ₅ × ℤ₇² (pair work about 2 · 10¹² per tuple).
- *m = 161*: H of order 3 (log₁₀ E ≈ 10, no listable quotient).

**Which estimate to trust.**
- *The discrepancy.* Two Gaussian estimators disagree by 10^3.6 on m = 215, ⟨49⟩, (14, 18, 13, 13):
  - per-class independent densities give 10^0.2;
  - the full covariance with the joint lattice gives 10^3.8.
- *Test family.* m = 5p with H of order 6 acting as −1 on ℤ₅ and with order 3 on ℤ_p, the shape of the live structure.
  Exact counts:
  - m = 35: 0 solutions in all 5 patterns;
  - m = 65: 352 solutions in (0, 4, 11, 11), 0 in the other 9;
  - m = 95: 0 solutions in all 14 patterns, although the full estimator predicts about 10³ in total.
- *Result.* The ℤ₅ compression has no feasible tuple for exactly the empty patterns at m = 65 and 95. On the live
  pattern, the ℤ₅-conditioned estimate is 10^2.86 against 10^2.55 exact.
- *Reading.* The unconditioned estimators are meaningless here. The ℤ₅-conditioned one is calibrated to within a
  factor of 2 on this family, which supports the 10^4.1 at m = 215.
- *Counterexample.* The ℤ₅ × ℤ₇² class nc:2054 (m = 245, |H| = 8) had conditioned estimate 36 and an exact count of 0.

**The m = 215, ⟨49⟩ join** (program not included).
- *Setup.* Lists A, B complete (42,875 each), C 2 · 10⁷ (kept from 60% of prefixes), D 4.03 · 10⁶ (from 20%). kt4 with
  merge classes 22, 18, 13, 20, 16; split values in order of the Gaussian law of a solution's split.
- *Coverage.* The 200,000 most likely split values (Mahalanobis radius ≤ 2.97, P(χ²₅ ≤ 8.82) = 0.88): 47% of the A × C
  pairs and 49% of the B × D pairs.
- *Result.* 0 hits, against about 2.5 · 0.88 ≈ 2.2 expected (Poisson P(0) ≈ 0.11). Unknown, but the conditioned
  estimate for these partial lists now looks too high by a factor of several. The first 20,000 values gave 0 against
  0.95.

## 76. A quaternary conference matrix of order 850, and H(6796)

**Result.** The unbordered zero-diagonal Szöllősi array (§74–75) has solutions over the non-cyclic group G = ℤ₅² × ℤ₁₇
with the symmetry group H = ⟨(M, 8)⟩ of order 24. Here M = [[0, 1], [1, 2]] acts on ℤ₅², M³ = 2I, so M is
multiplication by an element of order 12 of GF(25)*, and 8 has order 8 modulo 17. −1 ∉ H; there are 21 orbits.
- *Search.* The exhaustive compression-conditioned join (one row pattern, (0, 40, 7, 7), two feasible ℤ₁₇ tuples)
  finds 1280 solutions in 77 s; with z ↦ 2z instead, 512.
- *Chain.* A solution gives a quaternary conference matrix of order 850, the complex Yamada–Kiyasu template (q = 1697)
  a quaternary Hadamard matrix of order 3398, and Turyn's doubling H(6796).
- *Settles.* n = 1699 of Table 4 ([hadamard/6796](../hadamard/6796)).
- *The estimate.* The Gaussian estimate gave 10^0.65 and 10^0.52, against 10^3.1 and 10^2.7 exact: an underestimate by
  a factor of 300.
- *Checks.* `verify/verify.py` (construction `complex-yamada-kiyasu`) rebuilds the core from the four sequences,
  Yamada's circulants and K. The certification pipeline (an independent implementation; program not included), a
  bit-packed integer checker and a third rebuild from the sequences agree.
- *The template.* It needs only q a prime power ≡ 1 (mod 4) and any quaternary conference matrix of order (q + 3)/2:
  no symmetry, no normal form (proof in the builder's docstring). So it reaches the twelve
  [FRONTIER.md](../FRONTIER.md) entries n = 643, 859, 955, 1051, 1291, 1699, 2083, 2659, 2731, 2803, 2899 and 2971 (q
  = n − 2 ≡ 1 mod 8, core of order (n + 1)/2 ≡ 2 mod 4). With q ≡ 5 (mod 8) it also reaches 1879, 2215 and 2383, from
  cores of orders 940, 1108 and 1192.

**Structures tried.** Only the generic group search found a core.
- *Core 322 (n = 643).*
  - *The cpf filter* (proved). For a subgroup K ≤ G and a character ψ ≠ 1 of K, Parseval on G/K gives Σ_k Σ_cosets
    |Σ_{x∈K} X_k(g + x) ψ(x)|² = 4m − 2 (unbordered) or 4n + 2 (bordered). Bounding each coset term by its exact
    zonotope maximum over the orbit-constant signs gives a necessary condition.
    - Calibration: 96 small structures with exact counts, 0 false kills, 72 of 85 empty ones killed.
    - Over ℤ₁₆₁, all 8 subgroups with |H| ≥ 6 and −1 ∉ H are dead by it, among them ⟨24⟩ and ⟨114⟩, whose estimates in
      §75 were unconditioned. The survivors are ⟨93⟩, ⟨22⟩ and ⟨139⟩.
  - *−1 ∈ H is allowed.* All-symmetric solutions exist (624, 896 and 6624 at m = 9, 15, 21; cores 18, 30, 42 built).
  - *⟨93⟩* reduces exactly to two stages over ℤ₂₃. Its conditioned estimate is 10^10, but on the m = 21 analogue the
    same estimator predicts 85 and 17 for two empty patterns. Not searched.
  - *Bordered, |G| = 160.* Every H with 5 | |H| on ℤ₂⁵ × ℤ₅ is dead (exhaustive compression to G/ℤ₅). The live
    structures (ℤ₁₆ × ℤ₂ × ℤ₅ with |H| = 8; ℤ₂⁵ × ℤ₅ and ℤ₄² × ℤ₂ × ℤ₅ with |H| = 6) are uncalibrated and unsearched.
  - *Side results.* m = 215 ⟨4, 39⟩ and m = 323 ⟨8, 18⟩ (all-symmetric) are empty.
- *Character-bordered arrays.* The character-bordered zero-diagonal construction is proved: Xb Xb* = (2t − 1)I holds
  automatically, and the border exists iff the entries of Xb are 0 or units.
  - Every target has N ≡ 5 (mod 6), so the ℤ₆ family (t = 5, n = N − 5) applies to all of them.
  - Exhaustively empty at N = 161: ℤ₂² × ℤ₃ × ℤ₁₃ with |H| = 6 (the structure that gave the complex order 322) and
    ℤ₁₅₆ with |H| = 6.
  - ℤ₂ × ℤ₇ × ℤ₁₁ (Q = ℤ₇, t = 7, |H| = 5): 0 solutions in 15 tuple orbits, 6 of them with complete lists.
  - Over structured H the Q-conditioned estimate overshot every exact count (about 60 predicted in all, 0 found).
- *Direct constructions.*
  - *No group-developed object.* No group-developed quaternary Hadamard matrix of order 2n exists, since 2n is not a
    sum of two squares. A quaternary periodic complementary pair is a Goethals–Seidel quadruple with extra conditions.
    Bordered group-developed cores of order c are impossible for 13 of the 15 core orders, by self-conjugacy of an
    unramified prime.
  - *No reduction to smaller objects.* Kronecker, plug-in and doubling constructions cannot give c ≡ 2 (mod 4). Paley
    and Mathon parameters do not match.
  - *Relaxed cores* (proved, exact at orders 38 and 102). A relaxed core, with D = αI + βR in the template, also
    works. It comes from a bordered array with the zero moved from B to C, and is rarer than the standard one.
- *Generic group search.* A generic engine over any abelian group, any automorphism group and any invariant quotient
  (exact hash-and-verify joins, validated on brute-force counts; program not included) gave the result above.
  - *Exhaustively empty* (complete joins): bordered ℤ₂₃₈ ⟨81⟩; unbordered ℤ₆₆₅ ⟨4, 9⟩, ⟨62⟩, ⟨8, 11⟩ and ⟨39, 46⟩;
    ℤ₇₂₅ ⟨54⟩; four structures over ℤ₅² × ℤ₂₉ (|H| = 56 to 140).
  - *Dead by exact compressions:* most of the remaining structures with at most 72 orbits for m = 215, 323, 425, 665,
    725 and the bordered N = 215 to 743.
  - *Live and unsearchable:* bordered N = 263 over ℤ₂₆₂ with H = ⟨53⟩ (estimate 10^8.9, C and D pools about 10¹⁴);
    ℤ₂₁₅ ⟨136⟩ and ℤ₃₂₃ ⟨69⟩, whose estimates are unreliable (−1 on a factor).
  - *Z₂ lemma* (proved). In the bordered array with n = 2p and odd |H| = h, the Z₂-compression gives
    d_A² + d_B² + d_C² + d_D² = 8p + 2. This kills N = 215 bordered, and N = 263 except h = 5.
  - *Caveat on the tool of §75.* The compression-tuple generator of §75 packs PAF vectors into a 64-bit integer that
    wraps for d ≥ 8. Dead verdicts stay sound; live ones may rest on spurious tuples.
- *Certification pipeline* (program not included). It turns a core into a certified directory in under two minutes,
  even at c = 1486. Validated on 83 cases, including cores 142 → 1132 and 286 → 2284. A sweep confirmed Yamada's lemma
  for all 232 prime powers q ≡ 1 (mod 4) below 3000.

## 77. Non-cyclic structures for the other Yamada–Kiyasu targets

The order-850 effect (§76) was tried at the other targets over non-cyclic groups, unbordered and bordered. No core was
found. Every verdict below is an exact compression infeasibility, an exact border test or a complete join; none rests
on the Gaussian estimate.

**Which targets admit a non-cyclic group.** An abelian group of square-free order is cyclic. Among the unbordered
targets only m = 596 (core 1192), 725 (1450) and 245 (490, complex only) have a non-cyclic group, besides 425.

**Unbordered.** GAP listed every conjugacy class of subgroups of Aut(G). Class labels (c138, ...) and structure
numbers ([61], ...) below index lists that are not part of this repository.
- *m = 596, ℤ₂² × ℤ₁₄₉* (proved).
  - Parity for even m: a and b are odd, c and d even, and exactly one of c, d is 2 (mod 4).
  - The classes with 37 | |H| have no row pattern: the signed orbit sums are 0 (mod 74).
  - Every other class has |H| ≤ 6 and at least 150 orbits; out of reach. The cyclic ℤ₅₉₆ fares the same way.
- *m = 725, ℤ₅² × ℤ₂₉* (complete list: 400 classes without −1).
  - 248 classes have no solution: no row pattern, cpf-dead, or dead on a quotient.
  - Classes c138 and c143 (|H| = 28) were searched completely: 0 solutions.
  - The exact analogue of the 850 group (a Singer element of order 24 times an element of order 4) is dead on a
    quotient.
  - 21 live classes remain, the cheapest at 5 · 10¹¹ pair work (classes c139 and c142, |H| = 28 and 35 orbits,
    searched in §80).
- *m = 245, ℤ₇² × ℤ₅* (272 classes).
  - 220 classes have no solution; 9 more were searched completely, with 0 solutions.
  - 25 live classes remain, the cheapest at 3.8 · 10¹².
- *Defect.* The row-pattern generator of the earlier engine keeps only c > 0, so at even m it drops the patterns with
  c = 0. Earlier exhaustive verdicts at even m that used it are incomplete; a corrected version exists
  (program not included).

**Bordered and character-bordered.** The engine was extended to any target function, with the border tested exactly on
every tuple before any list is built.
- *Two character-bordered families* (not found in the literature searched) admit Singer-type symmetry:
  - Q = ℤ₃² with T = {1} plus two lines (t = 5), for N ≡ 5 (mod 9): 239, 743, 725, 554, 131, 149; order 28 checked
    exactly;
  - Q = ℤ₄² with t = 7, for N ≡ 7 (mod 16): 215, 263, 743.
- *Twisted family.* ℤ₂ × ℤ₃² can never carry a border (only 4 vectors in V_T).
- *Exhaustively empty.*
  - Bordered: N = 665 over ℤ₂³ × ℤ₈₃ and ℤ₂ × ℤ₄ × ℤ₈₃; N = 161 over ℤ₁₆ × ℤ₂ × ℤ₅ (the three |H| = 8 structures left
    live in §76).
  - The ℤ₃² family: one structure at N = 239, two at N = 131, and 25 at N = 725 with |H| = 48 to 192.
  - Most structures of the ℤ₄² family at N = 215, 263, 743.
- *Live, at 4.6 · 10¹⁰ to 10¹² pair work.* N = 725 in the ℤ₃² family (five structures, [61], [37], [74], [76] and [38]
  of the sweep's list; those numbered 96–286 were not triaged); N = 215, ℤ₄² family: [42] (a tuple exceeded the memory
  cap) and [36]; N = 239: [16] and [17] (§80, §82).
- *Observation* (heuristic). In the ℤ₃² family no tuple of the patterns (1, 1, 0, 4) or (1, 3, 2, 2) has had a border.
- *Out of the exact border solver's reach:* Q = ℤ₂⁴ (t = 7), ℤ₂ × ℤ₃² (t = 17; it admits the full GF(9)* Singer cycle
  and applies to N ≡ 17 mod 18, among them 161, 215, 323, 521, 665, 683, 701), ℤ₅² (t = 13), ℤ₁₁² (t = 21), ℤ₁₅
  (t = 11).

**Reading.** The 850 effect did not repeat. The Singer structures of the other targets die on row sums or on
quotients, and what survives has small or scalar H. §78 turns instead to a two-level sampler (the FLC method) on the
strongest cyclic structures.

## 78. Two-level sampling settles n = 859 and 1051

**The method.** A two-level sampler: compress, sample compressed instances, complete each exactly through its fibres,
as in the FLC classification method ([README.md](../README.md) reference [17]). FLC classifies circulant Williamson,
good and best matrices by:
- compressing each block by a small factor;
- listing the compressions that pass the PSD bound;
- joining the four lists on a hash linear in the compressed PAF (class partitions, a one-word Bloom filter, SIMD pair
  filters against the partner lists' minima);
- completing every instance (compressed quadruple) through the fibres of its blocks, with an (A, D) × (B, C) hashed
  join on the full PAF.

The two-level split is no shortcut for exhaustive search, because level-1 instances outnumber solutions (§75). When
solutions are plentiful, random instances completed exactly can cost far less per solution than any join of
full-sequence lists, which is what §75 to §77 used.

**Core 526, n = 1051** ([hadamard/4204](../hadamard/4204); program not included).
- *Structure.* The bordered array over ℤ₂₆₂ = ℤ₂ × ℤ₁₃₁ with H = ⟨53⟩ (order 5).
- *Split* (proved). With c = x(0, ·) + x(1, ·) and d = x(0, ·) − x(1, ·) on ℤ₁₃₁, PAF_x(0, s) = (PAF_c + PAF_d)/2 and
  PAF_x(1, s) = (PAF_c − PAF_d)/2. So the solutions are exactly the pairs (instance, completion):
  - level 1: Σ PAF_c = (8p + 2) δ − 8;
  - level 2: Σ PAF_d = (8p + 2) δ, with, away from 0, supp d the complement of supp c.
- *Necessary condition* (proved; not found in the literature searched). For p ≡ 3 (mod 4), 2(4h − 1)² ≤ 8p + 2. It
  explains every empty small case.
- *Calibration.*
  - Exact at p = 13, 41, 31; the last two reproduce earlier exact counts, 132096 at p = 41 and 764880 at p = 31.
  - p = 43, h = 3, a residue no earlier calibration covered: E = 10^7.91 against the Gaussian 10^7.87.
- *Population at p = 131.* N₁ = 1.3 · 10¹³ instances (Monte Carlo); about 6 · 10⁻⁵ completions per instance at the
  Gaussian E.
- *Zero-count stratification* (the gain).
  - The zero-orbit count of a solution's C or D is close to Binomial(26, 1/2). Exactly, 2z_A + 2z_B + z_C + z_D = 52,
    which pins z_A + z_B to (52 − z_C − z_D)/2.
  - The uniform C/D population peaks at z = 8 to 10.
  - Sampling only the stratum z_C = z_D = 18 raised the yield per pair test by about 10^3.6.
  - First solution after 103 s; 19 verified (17 distinct) in about 8 CPU-minutes.
  - Measured E ≈ 10^8.0–8.3, a factor 4 to 8 below the Gaussian 10^8.9.

**Core 430, n = 859** ([hadamard/3436](../hadamard/3436); program not included).
- *Structure.* ℤ₂₁₅ with ⟨49⟩, the model of §75.
  - Level 1: the compression to ℤ₄₃, g₀ = f₀ + 2f₁ + 2f₂, conditioned on the ℤ₅ tuple.
  - Level 2: the signs of v on the columns with |g₀| = 1, with Σ PAF_u + 5 PAF_v = 858 δ and Σ CCF(u, v) + CCF(v, u)
    = 0.
- *Key-group completion* (not found in the literature searched). All level-1 pairs with the same PAF sum are completed
  at once: the A × C fibre pairs are hashed and the B × D fibre pairs probe them. Per instance this would cost n_AC ·
  n_BD completions per group.
- *Calibration* (exact).
  - m = 65: 1408 solutions (§75 counts 352 in another normalisation; an independent brute force agrees).
  - m = 155 class 0: 1320 per ordered tuple from 7.7 · 10⁷ instances.
  - m = 35, 95: empty.
- *Trend.* Completions per instance: 0.92, 1.7 · 10⁻⁵, 4.2 · 10⁻⁷ and 3 · 10⁻⁸ at p = 13, 31, 37, 43.
- *Hits.* At p = 43, class 6 (rows (14, 18, 13, 13)): 3 hits with 1/64 of C, 1/64 of D and a quarter of the hash
  classes, the first after about 2 CPU-minutes. Class 15 (rows (6, 10, 19, 19)): 1 hit.
- *Estimate.* About 5 · 10⁴ solutions per ordered class-6 tuple (90% interval [1.4 · 10⁴, 1.3 · 10⁵]).
- *Why the join of §75 failed.* Its kt4 join multiplied Gaussian retention guesses for mass-selected full-sequence
  lists; a uniform prefix sample with exactly known coverage calibrates reliably.

**Core 322, n = 643** (program not included): no core; a sampled run, so unknown (§83 finds one).
- *Structure.* ℤ₁₆₁ = ℤ₇ × ℤ₂₃ with ⟨93⟩.
  - Level 1: the ℤ₂₃ quadruple of X₀ = c₀ + 6u.
  - Level 2: the signs of the V columns, with Σ PAF_Z = 642 δ for Z = R + i√7 E.
  - Exact at m = 21 and 35; cores of orders 154 and 182 built as controls.
- *Estimates* (heuristic).
  - The level-1 Gaussian is accurate to about 2×. The level-2 Gaussian overestimates completions 5 to 10×, more as p
    grows.
  - At m = 161: N₁ ≈ 10^18.4 and E ≈ 10^11, so cores should be plentiful. Completions per instance ≈ 10⁻⁷: about 56
    free signs against 22 conditions.
  - 900 instances completed (1 to 10 s each) gave 0.
- *Lesson.* The split pays when the fibres stay small (FLC: at most about 2²⁰ per block), as at 1051 and 859. With
  fibres of 2⁴⁸ to 2⁶⁰ it does not.

**Checks.** For both orders:
- the cores were rebuilt with independent code and K K* = nI was checked exactly;
- the certification pipeline (program not included) and a bit-packed checker confirmed H Hᵀ = nI;
- `verify/verify.py` rebuilds both from the four sequences. Construction `complex-yamada-kiyasu` supports the bordered
  form and the reflection R[g, c − g]; any c works, since R X R = Xᵀ.

The quaternary matrices of orders 1718, 2102 and 3398 have certificates in `complex-hadamard/`;
`verify/verify_complex.py` calls the same builder.

**Remaining structures of the ℤ₃² family at N = 725** (§77). Five structures were completed with 0 cores; a sixth has
no tuple with a border.

## 79. The two-level sampler over ℤ₂ × ℤ_p for composite p

No core was found, so 955, 1291, 2731 and 2971 stay open. The sampler of §78 was generalised to the bordered array
over ℤ₂ × ℤ_p for any odd p and any H ≤ ℤ_p*, and every subgroup with −1 ∉ H was triaged for N = 239, 323, 683 and 743
(p = 119, 161, 341, 371).

**The sampler** (program not included). Orbits may have unequal sizes (fixed points included), so row sums and zero
counts are weighted. Each weighted zero-count stratum is sampled exactly uniformly by a dynamic program.
- Validated against the sampler for prime p (p = 13: 1792; p = 31: 764880 for h = 3, 0 for h = 5), against an
  independent brute force without the split (p = 15 ⟨4⟩: 34560, ⟨11⟩: 19200; p = 13 with −1 ∈ H: 576) and against the
  generic group engine of §76 on 15 empty structures at p = 21 to 39.
- Positive control: at p = 21, H = ⟨8⟩, it found 104, 38417 and 16570 verified solutions; one was taken through the
  certification pipeline to H(684).

**Proved.**
- *Split for any odd p.* The identities of §78 hold with fixed points; c_A(0) = c_B(0) = 1 always, while c(0) = 0 is
  forced for C and D only when the nonzero H-orbit sizes have a common factor.
- *Weighted zeros.* Z_A + Z_B + Z_C + Z_D = 2p, with Z_k the number of zeros of c_k on ℤ_p.
- *Z_C ≡ Z_D (mod 4)* in every solution (not found in the literature searched). From the ℤ₂ lemma Σ D_k(1)² = 8p + 2
  with D_A ≡ −1 + 2Z_A (mod 8) and D_C = 2(Z_C mod 2), reduced mod 16 and combined with the zero count. Checked on
  121643 solutions at p = 15, 21, 31, 37 and 131: no violation. At p = 21 the joint stratum matters: (13, 13) gave 0
  completions in 16858 instances, (13, 9) gave 7169.
- *The bound 2(4h − 1)² ≤ 8p + 2 of §78 needs equal orbit sizes.* It fails for composite p: N = 743, H = ⟨15⟩
  (h = 13 > 9) has 16 feasible ℤ₂ tuples.

**Triage** (exact ℤ₂ and ℤ_d compressions, then the ℤ₂-conditioned Gaussian, log₁₀ E).
- *N = 239.* ⟨4⟩ (⟨81⟩ in ℤ₂₃₈, §76; exhaustively empty) −16.4; ⟨13⟩ +0.85; ⟨64⟩ (⟨183⟩ in ℤ₂₃₈) +3.22; ⟨18⟩ +8.3;
  ⟨50⟩ +14.3. Most larger H die on ℤ₂, ℤ₇ or ℤ₁₇.
- *N = 323.* Everything with |H| ≥ 6 is dead on ℤ₂ or has E < 10⁻⁴⁰, except ⟨93⟩ (+7.9, 69 orbits).
- *N = 683, 743.* Every live structure has E ≤ 10^−4.8.

**Calibration** (heuristic). The Gaussian is accurate to 0.9 in log₁₀ on composite p where H acts as −1 on a factor
(p = 15, 21). On the shape of N = 239 ⟨64⟩ (trivial on one prime factor, cyclic of order 4 containing −1 on the other)
it overestimates by at least 10²: p = 39 and 51 are exactly empty against 10^0.7 and 10^1.4 predicted, and p = 65 gave
0 completions over 88.6% of the hash classes against 10^2.4. A cost model checked on the sampler of §78
(optimistic by about 10×) puts N = 239 ⟨64⟩ at about 4 · 10⁵ CPU-years to a first core, ⟨13⟩ at 3 · 10⁶.

**Runs** (6 to 12 minutes each, 0 level-1 instances, as the model predicts): N = 239 ⟨64⟩ in the stratum Z_C = Z_D =
79 and uniformly; N = 743 ⟨15⟩; N = 683 ⟨5⟩. Exhaustive quotient joins cost 10²³ to 10⁴³ pair tests.

**Reading.** The bordered cyclic ℤ₂ × ℤ_p route is closed for these four targets: each structure is below one expected
solution or about 10⁵ CPU-years away. The two N = 239 structures left live in §77 could not be run here: the sweep
finds no usable quotient for them (both are closed in §82). Two structures of one N = 725 family have no bordered
tuple in the patterns (1, 1, 0, 4) and (1, 3, 2, 2); the pattern (3, 3, 0, 0) (about 3 · 10¹¹ pair work) was not
completed within an hour.

## 80. Why the order-850 solutions were plentiful; the bordered split for even or non-cyclic G′

No core was found.

**Anatomy of the order-850 solutions.**
- *Tied solutions.* All 1280 (H = ⟨(M, 8)⟩) and 512 (⟨(M, 2)⟩) solutions satisfy D = C ∘ n⁻¹ for an element n of the
  normaliser of H of odd–odd type: a non-square on GF(25) and a non-residue on ℤ₁₇. Also A ∘ n⁻¹ = −A and B ∘ n⁻¹ = B.
  The untied count is 0. Under N(H)/H (order 32), the reflections and swap of C and D, and A ↦ −A, they form 10 and 4
  free orbits of size 128; one core per orbit (14 in all) was built and checked.
- *A and B are explicit.* They are index functions: A(x, 0) = ±η₂₅(x), A(0, z) = ±η₁₇(z), and inside A(x, z) =
  f(ind₂₅ x ± ind₁₇ z mod 4) with f of type (+, +, −, −). B is 1 on both axes. The values |χ(X)|² on the axis
  character orbits are the same integers in every solution.
- *Tied join* (proved). With D = C ∘ n⁻¹, V_D = P V_C for a class permutation P, and the equation reads
  V_A + V_B + (I + P)V_C = T. Stage 1 joins A and B on the annihilator of Im(I + P), stage 2 scans C. It recovers all
  1280 solutions in 0.8 s (generic join 71 s). An independent brute force agrees on 850, on ℤ₈₅ (CD ties 4032, AB
  ties 9152) and on 510 small structures.
- *Why the estimate failed.* For the generic stratum the old estimate (4.6 solutions) is right; the whole count lives
  on the tied stratum, where (I + P)V_C has rank 4 instead of 9. A rank-reduced lattice Gaussian gives 187, and 456
  with the measured mod-8 residue bias, against 1280 (heuristic; 2.8× low instead of 285×).
- *Census.* Tied and symmetric-sum searches over every class with at most 60 orbits of the open groups
  (unbordered and bordered) found nothing; the plug-in shapes were the classes the engine could not handle.
- *Caveat on earlier negative verdicts.* The tied-join search (program not included) silently drops every n that does
  not preserve the kernel of the chosen compression (83% of the (n, tuple) pairs for class c139 of m = 725), and
  caches fibre permutations by object identity, which goes stale in long sweeps. So its negative verdicts are not
  exhaustive for classes with reducible H, and the conjecture that the quartic-twist family is absent from the open
  groups is not established.

**The bordered split for even or non-cyclic G′.**
- *Estimator.* The estimate table of §79 returned "null" on float64 underflow, not for a missing estimate. With
  log-sum-exp, N = 701 has E = 10^−1908 (|H| = 60) and 10^−943 (|H| = 30). Its Singer-type |H| = 24 structures have no
  feasible ℤ₂ tuple (exact).
- *Closed in this setting.* N = 701, 521, 665 and 725: E ≤ 10^−40 or exactly dead, except N = 725 with |H| = 9 (E =
  10^2.3 but about 10¹² CPU-years). N = 665: all seven structures with at most 11 A-orbits are exhaustively empty.
- *N = 161.* The two structures with E above 10 have |H| = 3 acting trivially on ℤ₅. That shape is unreliable for the
  estimator: its small analogues are exactly empty against E = 10^2.3 to 10^6.3. The next two candidates (G′ =
  ℤ₁₆ × ℤ₅ with H = ⟨(5, 2)⟩ and ⟨(3, 2)⟩) are dead on the ℤ₄ compression, which 1.7 hours of exhaustive stratum runs
  (0 level-1 instances) are consistent with.
- *Calibration.* On 181 small structures with exact counts, the estimate was never null or negative where solutions
  exist, but it gave false alarms up to 10^4.5 on empty structures.

**m = 725 exact joins** (n = 2899). Classes c139 and c142 of §77 (|H| = 28, 35 orbits), the two cheapest live ones,
are exhaustively empty: all six live row patterns of each, about 7.4 · 10¹¹ pair tests per class, 0 solutions. Their
tied strata are empty too, so class c141 (untied stratum only, 1.3 · 10¹² pair tests) was not run.

**n = 955, character-bordered ℤ₃² family** (§77). The two structures left live there are not unusable: the sweep
admits only quotients that refine Q unless that restriction is lifted. Without it both run on the quotient ℤ₃₉, but
`lowe_enum` then dies (SIGSEGV at the 3.5 GB cap) while listing one block (§82). Unknown at this point; §82 decides
both.

## 81. Paley plug-in cores

No core was found.

**The plug-in class.** G = K × F_q with H = {1} × S, S the squares of F_q* (H acts trivially on K). Each block is
three functions on K: x₀ on K × {0}, x₁ on K × S, x₂ on K × N. The class came from a calibration census run alongside
§80 (not tabulated there), where the order-153 class ℤ₃² × ℤ₁₇ has tens of thousands of solutions against 5 predicted.
In all of them x₂ = −x₁ at all but two of the 36 block positions: the blocks are x₀ ⊗ δ + x₁ ⊗ (Legendre symbol) plus
a few *defects* (x₁ = x₂).

**The K-level system** (proved; re-derived independently and checked against a separate brute force).
- With s = x₀ − e, t = (x₁ − x₂)/2, e = (x₁ + x₂)/2, each block is X = sI + tQ + eJ, Q the Paley matrix of F_q. The
  bordered array has the same form with its own right-hand sides. The unbordered condition Σ PAF = (4m − 2)δ_G is
  equivalent to
  - (I) Σ(s s̃ + q t t̃) = (4m − 2)δ_K,
  - (J) Σ(q e ẽ + s ẽ + e s̃ − t t̃) = 0,
  - (Q) Σ(t s̃ + ε s t̃) = 0, with ε = 1 for q ≡ 1 and −1 for q ≡ 3 (mod 4).
- *Defect identity.* The identity coefficients give e₊(q + 1) + e₋(q − 3) + e₀(q − 1) = 4|K| (unbordered), or
  4(|K| − 1) (bordered). Here e₊ and e₋ count defects with x₀x₁ = ±1, and e₀ ≤ 2 counts those at A(0), B(0).
- *q ≡ 3 (mod 4).* Symmetry of A and B forces x₁(0) = x₂(0), so e₀ = 2, plus a defect of A and of B at every 2-torsion
  point of K.
- *Bordered row lemma.* q ≥ |K| + 3 forces an even number of defects in every block.
- *Further exact filters.* The PE bound (defect differences must cancel between blocks when (q − 1)² is large), a
  row-sum lemma, and a quotient test that includes the identity coefficients.

**Results.**
- *n = 1291, both plug-in classes exhaustively empty.* ℤ₁₉ × F₁₇ (H of order 8, 57 orbits): only (e₊, e₋, e₀) =
  (1, 3, 1) is possible, giving three canonical defect configurations. Two die on a mod-8 row-sum lemma; the third has
  no solution of (I) + q(J) (a separate search over all 2⁵³ sign choices agrees). ℤ₁₇ × F₁₉ (order 9): only (0, 2, 2);
  one configuration, 128 × 21632 pairs, 0 hits (reproduced independently, 346112 t-quadruples).
- *Bordered, proved empty by the identity and the lemmas, and for 955, 2971 and 1046 also by exhaustive runs:* n = 955
  (ℤ₁₄ × F₁₇), n = 2971 (ℤ₁₄ × F₅₃), n = 2731 (ℤ₂₂ × F₃₁, the 66-orbit class no earlier sweep reached; e₀ = 2 leaves
  24 = 32a + 28b), and the complex-only order 1046 (ℤ₁₀ × F₁₃). The engine kills the same classes for 955, 2971 and
  2731 on exact compressions.
- *Impossible by the identity:* unbordered ℤ₅² × F₂₉ (2899), ℤ₇ × F₂₃ (643) and every even-m shape (1879, 2215, 2383);
  bordered N = 323, 701, 725, 665, 149, 215, 245.
- *n = 2899, ℤ₂₉ × GF(25).* Two of three configurations die on row sums. The third has no solution of (I) + q(J) with
  at most one exceptional coordinate (4 · 10⁷ candidates); the Gaussian expectation of a full solution is 10⁻¹¹
  (heuristic).
- *n = 2083, |K| = 40, q = 13:* 12 to 14 defects; the estimate is 10^−3.3 (ℤ₄₀) or lower; unknown.
- *n = 2659, ℤ₃₅ × F₁₉:* only (2, 4, 2) is possible; estimates 10⁻³⁷ and below; unknown.

**Existence trend.**
- *The plain plug-in runs out of freedom.* At q = 17 the plug-in has solutions for |K| = 9, 11, 13 (43392, 80 and 1344
  canonical) and none for |K| = 15, 17, 19.
- *The one plentiful target shape.* n = 643 (ℤ₂₃ × F₇, H = ⟨93⟩) has K-level estimate 10^11.9, but needs 12 to 22
  defects.
- *Ties.* Only 2.65% of the ℤ₃² × F₁₇ class is tied, so the plug-in class contains, but is not, the tied family of
  §80. A tie D = C ∘ n⁻¹ needs e₊ and e₋ even, so the (1, 3, 1) shapes admit none.

**Positive controls.** A core of order 442 over ℤ₁₃ × F₁₇ (plug-in, defects (2, 0, 1)) was taken through the
certification pipeline (program not included) to H(3532). A bordered core of order 42 over ℤ₄ × F₅ was taken to
H(332). Neither order is new.

**A parity lemma that fails.** A lemma stating that one configuration type has no defect-free solution of (I) + q(J).
It is false when (|K| − 1)/2 is odd (an explicit clean assignment exists at |K| = 19, and clean solutions exist at |K|
= 11 and 15). The ℤ₁₉ verdict stands on the exhaustive searches alone.

## 82. A sampler for n = 643, index-r plug-ins, the remaining n = 955 structures

No core was found. The only corrections are the three sub-statements noted below.

**n = 643: the plug-in class ℤ₂₃ × F₇** (H = ⟨93⟩ on ℤ₁₆₁, the only target plug-in class with a large expected count).
- *Defect structure* (proved; re-derived independently). For q = 7 the defect types are p (s = 0), m and z (the two
  forced zero positions). They satisfy 2n_p + n_m = p − 3 = 20 with n_z = 2, so E = n_p + n_m + 2 lies in [12, 22].
  There are 909 feasible allocations of defects to the four blocks. The K-level estimate is 10^12.2 ordered solutions
  (10^11.95 with the low-E correction).
- *Exact filters* on level-1 instances (not found in the literature searched): Σ PAF_{x₀} ≡ 0 (mod 3), and two row-sum
  identities (R0, R1). R1 forces Q = τ_C² + τ_D² ≡ 2 (mod 8) for the chosen allocation.
- *Tools* (levels 1 and 2 are `plug7_gen.c` and `plug7_l2.c`, see §83; the checks and the certification pipeline are
  not included).
  - Level 1 (x₀ quadruples sampled from mod-3 buckets, defect-pattern lists, class-restricted 4-list hash join) costs
    0.3 s per instance.
  - Level 2 (Gray-code fibres, AVX2 PSD pair filter, bucketed (A, D) × (B, C) hash join) costs 0.45 s per instance.
  - These are 5–10× and 8–14× faster than §78. Both equal brute force on complete sets at p = 3 and 5 and per instance
    at p = 11 to 23 (checked independently).
  - Certification was exercised end to end at p = 13 (core 182 → H(1452)) and on 24 p = 17 solutions
    (cores of order 238).
- *Calibration* (exact). p = 13: 1.2 · 10⁶ instances, 21914 completions over 87 allocations. p = 17, allocation
  (3, 3, 4, 4): 275421 instances, 48 completions (24 cores). By row-sum class Q: Q = 10 6/3713, Q = 18 0/11028, Q = 26
  10/23503, Q = 50 22/104337.
- *Rate at p = 23* (heuristic, three independent estimators). A partial-key collision profile of the level-2 join
  (count of 4-tuples matching the first j key coordinates, linear in j on a log scale, calibrated at p = 13 and 17)
  gives about 10⁻⁵ cores per instance for allocation (3, 3, 4, 4), i.e. about 23 CPU-hours per core. The class Q = 10
  is about twice as efficient: about 12 CPU-hours per core (range 4 to 30).
- *Run.* 2.0 CPU-hours over a mix of classes: 8790 instances, 0 completions (0.09 expected).

**Index-r plug-ins.** H = {1} × S_r with S_r of index r in F_q*, |H| = L = (q − 1)/r.
- *The system* (proved, validated on 27 shapes). Each block is r + 1 functions on K (one per cyclotomic class), with
  cyclotomic structure constants. The balance identity is Σ σ(Lσ + 2x₀) = 4r(|K| − ε) (ε = 0 unbordered, 1 bordered),
  where σ is the sum of the r class values at a position.
- *Corollaries.* Odd r needs L ≤ r + 1 (bordered also |K| ≥ r + 1). An exact profile/row-sum DP filter (two
  independent re-implementations, 0 false kills) proves empty, among others, r = 3 at n = 1291 (ℤ₁₇ × F₁₉), 2659 and
  2899, every r ≤ 12 of the bordered n = 2899 shapes over ℤ₄ and ℤ₂², and r = 2, 3 of bordered n = 2803 and 2731.
- *Three sub-statements were corrected.* The claim that the r + 1 conditions are Galois conjugates fails for prime
  powers q = p^e with gcd(r, (q − 1)/(p − 1)) > 1 (counterexample q = 9, r = 4). The bordered y = 0 equation has
  right-hand side −4, not 0. And "F real iff L even" fails for q = 2^e.
- *Existence trend.* r = 2 dies at target sizes, but r = 4 at q = 17 has solutions where r = 2 has none
  (exact 742400 at K = ℤ₅). The tilted-Gaussian estimator (within about 1 decade on non-rigid exact counts; unreliable
  in rigid shapes, in both directions) predicts 10^8.8 for n = 1291 (K = ℤ₁₉, r = 4) and 10^14.8 for bordered n = 2083
  (K = ℤ₄₀, q = 13, r = 4). It also predicts 10^16.7 for bordered n = 2803 (ℤ₂₈ × GF(25), r = 8) and 10^6.6 for
  complex 1958 (ℤ₅ × GF(49), r = 8).
- *The obstruction is search, not existence.* Blocks have (r + 1)|K| bits (95 for r = 4, K = ℤ₁₉). The s/t decoupling
  of r = 2 has no analogue. Simulated annealing stalls at density 10⁻¹⁵, and exact joins reach only |K| ≤ 5 (r ≥ 3).

**The remaining n = 955 structures.**
- *The `lowe_enum` crash.* The block list on ℤ₃₉ has 20⁶ = 6.4 · 10⁷ members. `lowe_enum` (`search/lowe_enum.c`,
  `leaf()`) does not check `realloc`: an allocation of 3.1 GB returned NULL and the program died with SIGSEGV. A fixed
  copy (checked realloc, streamed output, chunks) was validated; it is not included.
- *What decided the case* is the quotient Q = ℤ₃² itself, which every earlier sweep had dropped for its compression
  space (1.1 · 10⁷ > 3 · 10⁶). Of the two structures left live in §77, [16] has no Q-tuple in any row pattern and [17]
  is exhaustively empty: about 1.1 · 10¹¹ pair tests over 45 orbit representatives; an independent pipeline agrees.
- *The whole ℤ₃² family for N = 239* (34 classes, census confirmed independently): 19 classes die at the Q level;
  three more are closed (one in §77, and the two above); 12 remain, the cheapest at 2.6–5.8 · 10¹¹ pair tests per
  orbit representative, with lists above the memory cap.
- *Observation* (heuristic, 327 Q-tuples in 8 structures): rows (1, 1, 0, 4) and (1, 3, 2, 2) never carry a border,
  (3, 3, 0, 0) always does.
- *Other families at the Q level.* N = 743: 24 of 27 structures die. N = 725: 58 of the 287 structures of one family
  die, and none of the 46 of the other.
- *The bordered index-4 sub-plug-in over ℤ₁₄ × F₁₇* (|H| = 4, 70 orbits). The K-level system and the identity 832 =
  4(q − 1)(|K| − 1) were derived and checked independently, as was the lemma |Σ ε| ≤ 1 per block (so some position
  has |e| = 2; not found in the literature searched). Existence is unknown; the small bordered analogues K = ℤ₄, ℤ₆
  are exactly empty.

## 83. A core of order 322: n = 643 settled

The sampler of §82 (`plug7_gen.c`, `plug7_l2.c` and `plug7_sol.py`, in this directory) ran on four cores in the
row-sum class τ_C² + τ_D² = 10 of the allocation (3, 3, 4, 4).
- *Hit.* After 42451 instances of that class (about 13.5 CPU-hours, plus the 2 CPU-hours of §82), one instance gave 2
  completions, i.e. one core up to the sign of t. The calibrated estimate was about 12 CPU-hours per core.
- *The core.* ℤ₁₆₁ = F₇ × ℤ₂₃, H = ⟨93⟩. Row sums (−20, −4, 15, −1). 14 defects: A at k = 0, 10, 13; B at 0, 4, 19; C
  at 9, 10, 15, 19; D at 0, 12, 19, 20 (k the ℤ₂₃ coordinate). The core is neither symmetric nor Hermitian.
- *Checks.* `verify.py` rebuilds the core from the four sequences (reflection 0), then the template and the doubling
  (q = 641, r = 3, ξ = 3 + w), to H(2572), SHA-256 ebf232c3…; `verify_complex.py` gives the quaternary matrix of order
  1286, SHA-256 5e2c874c…. A second implementation (not included) produces the same hash. A bit-packed check of all
  3306306 row pairs of the CSV finds every pair orthogonal.
- *Novelty* ([NOVELTY.md](../NOVELTY.md)). No construction of H(2572), of a quaternary H(1286) or of a quaternary
  conference matrix of order 322 was found.

**Reading.** The Paley plug-in turned a class that round 23 had sampled without success (completions per instance
about 10⁻⁷ at 1–10 s per instance) into one that yields a core in half a CPU-day. The gain came from four sources:
the defect identity fixes the defect counts; R0, R1 and the mod-3 filter cut the level-1 population exactly; the
row-sum class is chosen by calibration at p = 13 and 17; and both levels are 5–14× faster.

## 84. Rounds 28–29: bordered plug-ins, the limits of Yamada–Kiyasu-type templates, the transfer theorem (October 21–22)

No core was found. A container restart cut off two lenses (closed forms, sampler speed-up) and the judges of the
template and transfer lenses, so the claims below marked PROVED are the lenses' own, not yet independently checked.
Scratchpad `r28/`, `r29/`.

**Bordered plug-ins over ℤ₂ × ℤ_p × F₇** (n = 1291, p = 23; n = 955, p = 17; `r28/bplug/`).
- *Checked parts.* The exact bordered K-level system, its row-sum identities, the defect identity, the forced defects
  and the ℤ₂ (c/d) split were checked against G-level PAF. A new exact sampler `bp9.c` was validated against brute
  force on 32 planted cases.
- *Positive controls.* The first bordered q = 7 plug-in cores, of orders 198 (K = ℤ₁₄) and 254 (K = ℤ₁₈), were
  certified; the judge verified both. p = 3 and 5 are exhaustively empty.
- *Infeasible.* The classes have Gaussian expectations of 10^11.7 (p = 17) and 10^12.3 (p = 23), but only one
  level-1 instance in 10^11 to 10^12.6 (p = 17), or 10^18.5 to 10^21 (p = 23), completes. The ℤ₂ split, ties,
  strata and extra symmetry give no gain of the 10⁶ to 10⁸ needed. Every solution found at p = 7, 9 lies in a
  degenerate class (one block nearly constant) that cannot persist at p ≥ 17. HEURISTIC; existence UNKNOWN.

**Yamada–Kiyasu-type templates** (`r29/I2/`; not yet judged).
- *The skeleton.* Complex Yamada–Kiyasu is a ℤ₂-split Goethals–Seidel/Szöllősi array. Its W-channel is a "G-object"
  of order m = (n − 1)/2, i.e. a weighing matrix W(2m, 2m − 1) of the form [[M, N], [−N, M]]; the S-channel is the
  core. The template works for any G-object, circulant or not (PROVED, verified numerically).
- *The obstruction.* By Bruck–Ryser/Raghavarao, n − 2 must be a sum of two squares when n ≡ 3 (mod 4). When
  n ≡ 1 (mod 4) the core would have odd order, which is impossible. So 79 of the 124 open n have no template of this
  type, among them 599, 719, 739, 751, 839 and 863. 34 more would need a non-Paley symmetric conference matrix of
  order n − 1 with this structure.
- *Complex Miyamoto* (Lemma 8.6) is the same skeleton with the two channels swapped. For q ≡ 3 (mod 4) it fails
  only because its V-part would need a quaternary Hadamard matrix of odd order (PROVED).
- *Further results.* A generalisation n = q + 2w was derived but not tested beyond w = 1. T-matrix multiplication of
  the complex-YK matrices fails, because the core blocks do not commute (REFUTED on the core of order 322).

**Transfer theorem and K-symmetric plug-ins** (`r29/I1/`; the judge did not finish).
- *The theorem* (PROVED, verified numerically). The plug-in system and its defect identity depend only on
  (K, q, sign). Any symmetric conference core of order q can replace Paley's. A skew core works iff it has an
  involution P with PQP = Qᵀ, which has exactly one fixed point; none exists for q = 15.
- *Verified with a non-Paley core.* A quaternary conference matrix of order 950 = 2 · 19 · 25 was built from a
  Latin-square conference graph not isomorphic to Paley(25).
- *No new target.* The identity forces q ≤ 4|K| + 3, and every non-prime-power core order up to 55 (45; skew 15, 35,
  39, 51, 55) fails it at every target.
- *K-symmetry does not rescue the large-K targets* (2659, 2899, 2803, 2083, 2971): E < 1 below 35 orbits, and
  ρ ≤ 10⁻¹⁷ where E ≥ 1 (HEURISTIC).
- *Unification.* The core of n = 859 (§78) is itself a K-symmetric plug-in: ℤ₂₁₅ = ℤ₄₃ × F₅ with ⟨49⟩ = (order 3
  on ℤ₄₃) × (the squares of F₅). Both 859 and 643 sit where the level-2 sign count n_t is 32 to 56 and the completion
  probability ρ ≥ 10⁻⁷.
- *Scaling law* (HEURISTIC). At fixed |G| a larger q makes the sampler much cheaper. At q = 7, E peaks near
  |K| = 41 while ρ falls by about 10^−0.6 per unit of |K|.
- *One lead.* The complex order 1958 (core 490) via K = ℤ₃₅, q = 7, H_K = ⟨11⟩: E ≈ 10^3.7, n_t = 33, ρ ≈ 10⁻⁷,
  about 30 to 100 CPU-hours with a naive sampler.

## 85. Round 30: Paley plug-ins in the plain Szöllősi array; the complex order 418 (October 22)

The quaternary Hadamard matrix of order 418 was found. Scratchpad `r30/` (census and order-1958 lenses).

**The plain-array plug-in** (`search/plugh`). As §81, but in the plain array (all entries ±1, A and B symmetric, Σ PAF
= 4|G|δ), which gives quaternary Hadamard matrices of order 2|G|.
- *System.* The K-level equations (I), (J), (Q) of §81 hold with 4|G|δ_K on the right of (I). The defect identity is
  e₊(q + 1) + e₋(q − 3) = 4|K|; for q ≡ 3 (mod 4) A and B have a defect at 0.
- *Checked.* Against G-level brute force on 289 small shapes (26 nonzero, e.g. ℤ₁₃ × F₃ with H = ⟨5, −1⟩: 14592 both
  ways). The C solvers `l1.c` (level 1: exact A, B lists against random C, D configurations, class-sorted hash join)
  and `pl2.c` (level 2: row-sum filter, PSD channels, meet-in-the-middle on the signs t) reproduce the exact counts of
  five shapes (ℤ₁₉ × F₁₉, H = ⟨7⟩: 3072; ℤ₁₉ × F₇, ⟨7⟩: 159744; ℤ₁₃ × F₇, ⟨3⟩: 27648; ℤ₁₃ × GF(9), ⟨3⟩: 0).

**Census of the open orders** (Gaussian/lattice estimator, then exact or sampled runs; 38 orders ranked by cost).
- *Estimator.* Normal regime (q < |K|, H = 1): estimate/measured completion rate per level-1 instance 0.6 (ℤ₁₃, q =
  7), 0.15 and 0.23 (ℤ₁₃, q = 9), 0.2 (ℤ₁₅, q = 7), > 4.4 (ℤ₁₇, q = 11: 0 of 21586) and 0.025 (ℤ₁₉, q = 11). So it is
  good to about 1.5 decades either way, with no systematic correction. Rigid regime (q > |K|, or few orbits): six of
  six shapes with estimates 10^2.3 to 10^4 are EXHAUSTIVELY empty: ℤ₁₁ × F₂₃ (506), ℤ₁₁ × GF(25) (550), ℤ₁₅ × F₂₃
  (690), ℤ₄₁ × F₇ with ⟨10⟩ (574), ℤ₁₉ × F₁₁ with ⟨7⟩ (418), ℤ₃₁ × F₁₁ with ⟨2⟩ (682). There the level-1 count is
  underestimated 20 to 30 times and full solutions do not exist.
- *Small K-symmetric shapes.* 470 more (4 to 13 orbits per block, |K| ≤ 120, over the 85 orders that pass the
  identity) are EXHAUSTIVELY empty; almost all die at level 1.
- *Order 418* (ℤ₁₉ × F₁₁, H = 1, n_t = 51, estimate 10^5.7 solutions): the only order within reach. Two workers found
  solutions in their first rounds (0.6 CPU-hours, 7641 level-1 instances, 40 completions, five certified distinct).
  The certificate has defects (e₊, e₋) = (5, 2) in the allocation (1, 1, 2, 3), row sums 1, 21, −15, 13 (1 + 441 + 225
  + 169 = 836). `verify_complex.py`, the sampler's own check and a separate numpy check (with the Turyn double H(836))
  give K K* = 418 I exactly.
- *Out of reach of this level 1:* 506 (ℤ₂₃ × F₁₁), 598 (ℤ₂₃ × F₁₃), 550 (ℤ₂₅ × F₁₁), about 10² to 10⁶ CPU-hours by the
  estimates. For |K| ≥ 23 the exact A, B lists and the random C, D lists are both too large; a structured level 1
  (random background, exact join on the defects, as `plug7_gen.c` for n = 643) is needed. 1026 (ℤ₁₉ × GF(27), ⟨7⟩) and
  1258 (ℤ₁₇ × F₃₇) were stopped without a result (UNKNOWN).

**Complex order 1958** (conference core of order 490 over ℤ₃₅ × F₇, H_K = ⟨11⟩, §84). Exact engine validated
exhaustively at K = ℤ_m, m ≤ 3, and on m = 4 (counts 40, 64, 198, ≥ 128 canonical). A bug in its canonical test (the
block-sign gauge) made the first production run incomplete; it was fixed and the run repeated. Production: 7.3 · 10⁶
canonical level-1 instances by 07:08 UTC, 0 solutions; the probe-based expectation of the unprocessed part is about
0.4 (HEURISTIC). The small-m solutions all lie in one or two structured allocations that the generic estimator does
not see.

## 86. Rounds 31–32: real Goethals–Seidel plug-ins, a structured level 1, the order-506 campaign (October 22)

No new order. Scratchpad `r31/`, `r32/`.

**Real plug-ins** (four free blocks: unbordered GS over K × F_r for composite n, Spence over a group of order n − 1,
bordered GS over n − t). The K-level system is that of §81 with no symmetric block; the defect identity is
(r + 1)e₊ + (r − 3)e₋ = 4|K| + φ(0) with φ(0) = 0 (GS), −4 (Spence), −4t (bordered).
- *Checked.* Against independent G-level brute force on about 550 small shapes (r = 3 to 13, cyclic and non-cyclic K,
  H_K trivial, nontrivial and twisted, Q = ℤ₂ to ℤ₅ and ℤ₂²), 0 mismatches; certificates of known orders (36 to 836)
  pass `verify.py`.
- *Census of the 124 open n.* 6695 shapes (4 to 28 orbits per block); the exact filters leave 2729 and the estimator
  gives at least one solution for 14 shapes on 12 n. All of them are now decided and empty: Spence 599, 947 (two
  shapes), 1427 and GS 1981 by exhaustive search; Spence 1031 and 2711 by a character-sum test; GS 955 (ℤ₁₉₁ × F₅,
  ⟨5⟩) because no allocation satisfies the trivial-character identities; Spence 947 (ℤ₂ × ℤ₄₃ × F₁₁, ⟨(1, 4)⟩) by an
  exhaustive search of all 466 allocations; 739, 2315 and ten further survivor shapes likewise. An independent
  judge re-derived the identities and re-ran 947 and 955 with its own code. The estimator gave up to 10^6.9 solutions
  on exactly empty small shapes and over-predicted the 947 level-1 count by 3 to 5 decades.

**A structured level 1** (`r31/sl1`): the background is fixed by an exact congruence (Σ PAF of the backgrounds ≡ 0
mod q − 1 at every nonzero shift), then every defect pattern of the four blocks is joined exactly. It is 44
times faster than round 30's level 1 at ℤ₁₉ and about 2400 times at ℤ₂₃. A symmetry-reduced exhaustive mode
(translations and reflections of C and D, units, A ↔ B) settles whole allocations.
- *Checked.* Equal to the exact level-1 counts of six shapes; sampling rate × recall reproduces exact solution counts
  (27136 ± 735 against 27648, 159252 ± 1128 against 159744). An adversarial audit proved each symmetry and compared
  the reduced search with its own complete enumeration on 188 (N, q, allocation) cases: 0 missing solutions.
- *Results.* At ℤ₁₉ × F₁₁ the allocation (1, 1, 2, 3) holds exactly one solution orbit (the order-418 core of §85) and
  (1, 1, 2, 4) a second, different one; at ℤ₁₇ × F₁₁ one orbit gives an order-374 matrix (not open). Round 30's empty
  ℤ₁₇ calibration row was a sampler failure; the exact completion rate there is 2 · 10⁻³.
- *Calibration.* In orbit units the round-30 estimator over-predicts by a factor κ ≈ 6 (2 to 16) for prime |K|, and
  under-predicts for composite |K| (ℤ₁₅).
- *A faster level 2.* `pl2f` is 4 to 5 times faster than `pl2` with identical output on 181,000 instances (two
  independent comparisons, including every known core and brute-forced small shapes); the production build avoids
  AVX-512.

**Complex order 506** (ℤ₂₃ × F₁₁, H = 1, 70 feasible allocations). Exhaustively empty (symmetry-reduced):
(1, 1, 3, 4), (1, 3, 2, 3), (1, 1, 2, 5), (1, 1, 2, 4), which carried about half of the estimator's mass. Running:
(1, 3, 3, 3), (1, 3, 2, 4), (1, 1, 3, 5). The remaining expectation is about 0.4 solution orbits (HEURISTIC).

**The other open complex orders** (re-ranked with the new costs). For eight of them (334, 358, 694, 718, 958, 1438,
1678, 1726) the defect identity rules out every plug-in shape. Newly EXHAUSTIVELY empty: 1258 (ℤ₁₇ × F₃₇), 1026
(ℤ₁₉ × GF(27)), 1806 (ℤ₂₁ × F₄₃), 690 (ℤ₁₅ × F₂₃), 682 (ℤ₃₁ × F₁₁ with H of order 2 or 5) and 574 in the bordered
array. Only 550 (ℤ₂₅ × F₁₁, about 40 to 60 CPU-hours per allocation, a few per cent each) is within reach after 506;
it is running.

## 87. Round 33: skew Paley products with many w-positions; skew 476 and symmetric 7612 (October 22–23)

Scratchpad `r33/`. The skew-Hadamard matrix of order 476 (v = 119, census-undecided) was found, and from it the
symmetric Hadamard matrix of order 7612 (v = 1903, census-undecided).

**The skew plug-in** (`search/skewplug`). §15's Paley products in the K-level form of §81: G = K × F_q with q ≡ 3 (mod
4), H = {1} × S, blocks sI + tQ + eJ. B, C and D are free; the skew block A has s and e odd under k ↦ −k off 0, s(0) =
1, e(0) = 0 (no w-position at 0, the others in pairs ±k of the same type) and t even. The defect identity is (q + 1)e₊
+ (q − 3)e₋ = 4|K|, §15's (q + 1)N_w − 4N₋ = 4m. The solution set is invariant under translations of B, C and D, their
negations, the reflection of each block (A ↦ 2δ − A), and the units of K: a group of order 2⁷|K|³φ(|K|) (1.0 · 10⁷ for
K = ℤ₁₇), times 6 for the permutations of B, C, D.
- *Checked.* G-level brute force equals K-level brute force on 31 shapes (12 nonzero, e.g. ℤ₅ × F₇: 74,880,000
  quadruples), and the two-level pipeline equals both on 28 (10 nonzero). A judge's own brute force reproduced all 31
  counts. Each symmetry was tested on 2000 random solutions per map from the complete solution sets of ℤ₃ × F₇
  (176,256) and ℤ₃ × F₁₁ (110,592); negating A or swapping it with another block fails, and translating A mostly
  fails.

**v = 119 = 17 · 7.** 2e₊ + e₋ = 17, so 9 to 17 w-positions: the regime §15 could not search. The first round of the
sampler (allocation (2, 2, 2, 3), the exact list of skew blocks against sampled lists of the others) gave the matrix
([skew/476](../skew/476)). Three more solutions came from the allocations (2, 2, 3, 4), (0, 2, 3, 5) and (4, 1, 3, 4),
and the judge found a fifth with its own sampler and level 2. One allocation slot alone has about 5 · 10¹³ level-1
instances, and 15 of 20 sampled ones complete at level 2; the class holds an estimated 10¹⁹ to 10²⁰ ordered quadruples
(HEURISTIC). For q = 7 the compression entries are ±1, ±5, ±7, and (q − 2)² = 25 is within the noise of the
autocorrelation sums, so no cancellation between big entries is forced and level 1 is huge; for large q the regime is
rigid (below).

**The cascade.** Doubling gives skew 952 = 8 · 119, below Cati–Pasechnik's 119(4), and with q = 1901 ≡ 5 (mod 8) the
symmetric Yamada–Kiyasu theorem gives the symmetric Hadamard matrix of order 7612
([symmetric/7612](../symmetric/7612)). The other six Yamada–Kiyasu targets (v = 1711, 1879, 2215, 2359, 2383, 2391)
need skew orders 856, 940, 1108, 1180, 1192, 1196, none of which has a usable Paley-product shape.

**The other skew census values with a many-w-position shape** (|K| > q). 783 = 29 · 27 and 1271 = 41 · 31 fail the
identity for every q ≥ 7. Every K-symmetric shape that could be run to the end is EXHAUSTIVELY empty: 589, 779, 575
(also over ℤ₅ × ℤ₅), 627, 759, 943, 473, 817 and 989 (for the last three only H_K of order 7 and 21; order 3 is
undecided), and about thirty further undecided values with at most 13 orbits. A bound on the products of big
compression entries proves 361 = 19 · 19 and 529 = 23 · 23 empty for H_K = 1 (the estimator had ranked an empty slot
of 361 first, at 10^8.6); 361 with H_K = ⟨7⟩ is empty by an exhaustive run. Shapes with H_K = 1 and |K| ≥ 25 have more
than 21 free t-orbits per block, beyond the present level 2; 575 = 25 · 23 is the best of them (HEURISTIC). Sampling
245, 385, 473, 623 and 805 found nothing (HEURISTIC).

**Symmetric matrices.** For group-developed blocks the propus array (A, B, B, D) is the only symmetric array
(`search/sym_arrays.py`). Its Paley plug-in was checked against G-level brute force on 230 shapes (27 nonzero); the
positive controls ℤ₁₁ × F₇ and ℤ₁₃ × F₇ have one and three solution orbits (orders 308 and 364, known). At v = 119 all
245 allocations are EXHAUSTIVELY empty (symmetry-reduced, 0 solutions), as the prior estimate of 10⁻² orbits suggested
(HEURISTIC). The small-N solutions use special blocks (a Legendre first block, at N = 13 the Singer (13, 4, 1)
difference set); their analogues at N = 17 gave nothing.

## 88. Round 34: a multiplier on K, negacyclic arrays; skew 1892, 1528 and 2216 (October 23)

**A multiplier on K.** §87's skew Paley product with H = H_K × S, H_K ≤ K* of odd order acting on every block.
*Twisted generators are impossible* (EXACT): if H contains (h_K, ν) with ν a non-square of F_q, then by Euler's criterion
it contains (h_K^{(q−1)/2}, −1), which maps (0, f) to (0, −f), so an H-invariant block has A(−g) = A(g) there,
against the skew condition. The same holds with field automorphisms for q = pᵉ, e odd (the norm of a non-square is a
non-square). A brute-force check agrees: none of 2262 generators with a non-square ν (nine small K) admits an
invariant skew block, and over 1512 generators of seven small K those that do are exactly those of odd order. `search/skewplug` gained the production loop `prod_ks.py` (all level-1 instances
of a round go to level 2; certification on a hit) and the level 2 `pl2f2.c` (`skf2`), equal to `gl2.c` in completions
and printed solutions on 7500 level-1 instances of seven shapes, 1.8 to 6.7 times faster.

**v = 473 = 43 · 11, H_K = ⟨6⟩ of order 3** (§87: undecided). 15 orbits of ℤ₄₃, 45 per block; (q + 1)N_w − 4N₋ = 172.
The tenth round of the loop gave [skew/1892](../skew/1892), allocation (0, 4, 7, 10), N_w = 21, N₋ = 20. Four workers
found five solutions with different row sums or allocations (inequivalence not tested; row sums (1, −9, 17, −39), (1, 21, −19, −33), (1, −19, −3, −39), (1, 17, −9, −39),
(1, 35, 15, 21)) in 133,196 level-2 runs. The completing level-1 instances of a round form one cluster (2 · 2 · 2
combinations of two A, two B and two C with one D, each with one t-solution up to t ↦ −t). The Gaussian estimate of
N1 at the hit slot is within a factor 1.5 of the measured one, but every level-2 estimate (Gaussian, and an
instance-level local limit over the exact fibres) is 3 to 4 decades too low, here and at 119: no estimator available
can discard a K-symmetric shape (HEURISTIC). For H_K = 1 and |K| ≥ 25 the exact level 2 needs 2⁴⁴ (ℤ₂₅) to 2⁷⁷ (ℤ₄₃)
pair evaluations per level-1 instance after the single-block filters (measured fibres); nothing was run there.

**Other K-symmetric and rigid shapes** (`search/skewplug` plus a pattern pipeline: a big–big filter, an exact
necessary condition on the positions of the w-positions from the products of large compression entries, sound on
the known solutions of seven shapes; symmetry reduction; an exhaustive per-pattern level 1 and level 2, equal to the allocation
pipeline on three nonzero shapes and recovering [skew/1100](../skew/1100) from its pattern).
- *EXHAUSTIVELY empty:* 589 (ℤ₃₁ × F₁₉, ⟨5⟩), 627 (ℤ₅₇ × F₁₁, ⟨4⟩), 759 (ℤ₆₉ × F₁₁, ⟨4⟩), 979 (ℤ₈₉ × F₁₁, ⟨2⟩),
  1539 (ℤ₈₁ × F₁₉, ⟨4⟩), 623 (ℤ₈₉ × F₇, ⟨2⟩: 374 allocations, 1.6 · 10⁸ level-1 instances, none completes), 989
  (ℤ₄₃ × F₂₃, ⟨6⟩), 1209 (ℤ₃₉ × F₃₁, ⟨16⟩), 817 (ℤ₄₃ × F₁₉, ⟨6⟩: 457 pattern classes after the big–big filter, 317 after the
  row-sum filter), 805 (ℤ₃₅ × F₂₃, ⟨11⟩); with
  non-cyclic K or matrix H_K: 539 (ℤ₇² × F₁₁, order 21), 1323 (ℤ₇² × GF(27), order 21), 1029 and 2401 (ℤ₇³ × F₃,
  F₇, Singer subgroup of order 171), 1953 (ℤ₃² × ℤ₇ × F₃₁, order 9); and 575 (ℤ₂₅ × F₂₃, ⟨21⟩).
- *EXACTLY empty:* 851 (ℤ₃₇ × F₂₃, ⟨10⟩: the orbit weights cannot realise the identity), 897 (ℤ₃₉ × F₂₃, ⟨16⟩: all
  7332 patterns fail the big–big bound) and 1333 (ℤ₄₃ × F₃₁, ⟨6⟩); 1323 (ℤ₇² × GF(27), H_K of order 7 or 3: the
  row-sum filter kills every pattern that passes the big–big bound); 1859 and 2375 with their Singer groups (weights);
  for H_K = 1 every pattern dies at 361, 529, 899, 961, 1763 and 1269 (ℤ₄₇ × GF(27)).
- *The N_w = 2 lemma* (EXACT). For q = 2m + 1 the identity allows two w-positions only as allocation (0, 0, 1, 1), one
  of each type. Translating both to 0, the level-1 equation at a lag k ≠ 0 reads q ε_C a_k + (q − 2) ε_D b_k + SS(k) = 0
  with a_k = s_C(k) + s_C(−k), b_k likewise and |SS(k)| ≤ 4m − 4 < 2(q − 2); so a_k = 0 ⟺ b_k = 0 and ε_C a_k =
  −ε_D b_k, whence |σ_C| = |σ_D + 2ε_D| for the row sums σ of s. Of the exact row-sum patterns 2 of 7 survive at m = 11
  (exactly those with level-1 instances), 0 of 4 at m = 23 and 0 of 6 at m = 29: **1081 = 23 · 47 and 1711 = 29 · 59
  have no solution in this product for any H_K**. 1081 was re-checked by exhaustive class lists and an independent join.
- *Sampled without a hit* (UNKNOWN): 245 (ℤ₃₅ × F₇, ⟨11⟩: 5.3 · 10⁵ level-1 instances), 385 (ℤ₃₅ × F₁₁, ⟨11⟩ and
  ℤ₅₅ × F₇, ⟨16⟩), 777, 781, 957, 1749; 627 (ℤ₅₇ × F₁₁, ⟨7⟩), 1045, 1141, 1243, 1267, 1393, 1397, 1617 and 2189 gave
  no or almost no level-1 instances.

**Negacyclic Goethals–Seidel arrays** ([skew/NEGACYCLIC.md](../skew/NEGACYCLIC.md); `search/negaskew`). For v odd,
ℤ[x]/(x^{2v} + 1) ≅ ℤ[i][ℤᵥ] (x ↦ i · g), a ±1 negacyclic block is (1 + i)U with U quaternary on ℤᵥ, and the array with
negacyclic blocks is Hadamard of order 8v exactly when Σ U_X U_X* = 4v δ. A palindromic first block (skew type) is
u(−g) = i · conj(u(g)) for g ≠ 0, u(0) ∈ {1, −i}. Under a multiplier group H of odd order a skew block has 2^{n_H} choices and a
free block 4^{n_H}; the real skew class with the same H has 2^{(n_H − 1)/2} and 2^{n_H}, with half the equations.
- *Validation.* In real mode `qskew` reproduces the exact counts 176,256 (ℤ₃ × F₇), 110,592 (ℤ₃ × F₁₁) and 74,880,000
  (ℤ₅ × F₇) of §87; quaternary solutions over ℤ₃ and ℤ₅ give skew-Hadamard matrices of orders 24 and 40.
- *Sweep.* All multiplier groups of odd order with n_H ≤ 13 orbits at the v ≤ 999 with Table 5 entries v(m), m ≥ 4,
  and the v ≤ 2499 undecided in the census: 373 shapes (371 on the present repository, which settles 473 and 933),
  2 hits, 370 EXHAUSTIVELY empty, ℤ₂₁₉₇ unfinished; GF(pᵏ)
  shapes at 361, 841, 1369, 1849 and 2197 empty. The real class is empty on all 510 shapes with n_H ≤ 17 (508 now).
- *Hits.* v = 191, H = ⟨180⟩ of order 19, 11 orbits, 4 s: [skew/1528](../skew/1528); v = 277, H = ⟨273⟩ of order 23,
  13 orbits, 3 minutes: [skew/2216](../skew/2216). Table 5: 191(226) and 277(5) become 191(3) and 277(3).
- *Row sums.* With k = |H|, the sum of a block is u(0) + kS with S a Gaussian integer, and for the skew block
  S = m(1 + i) (each pair of orbits ±O contributes w + i · conj w = ±(1 + i) per element). In both hits three sums are
  units and the fourth is u(0) + k(±1 ± i): R_A = −19 − 20i at 191, R_D = 24 − 23i at 277. So 4v − 3 = k² + (k + 1)²,
  that is v = 1 + k(k + 1)/2. The prime v of this form with k ≡ 3 (mod 4) are 7, 29, 67, 191, 277, 379, 631, 947,
  1129, 1327, … (H the subgroup of order k, (k + 3)/2 orbits); whether all of them carry solutions is open.
- *Paley products in this ring.* The skew block always has a half-defect at (0, f) (x₂(0) = i · conj x₁(0)), of
  weight (q + 1)/2 or (q − 3)/2; for q ≳ 3|K| its products with the other big entries cannot cancel unless another
  block locks onto it, so single-defect shapes are empty (413) and the rigid candidates are 473, 517, 1139 and 1207
  (UNKNOWN). Spin shapes (C = B ∘ μ, D = B ∘ μ², μ³ ∈ H) with at most 19 orbits fail the row-sum condition or are
  EXHAUSTIVELY empty (817), as is 589 with ⟨574⟩ and μ = 5 at 21 orbits.
- *Coverage of Table 5.* Of its 172 entries with m ≥ 4, 463 is not census-undecided and 29 had a repository matrix of
  order 4v. Of the other 142, 191 and 277 are the hits; of the remaining 140 (473 and 933 are settled by now), 77 have
  a multiplier shape with at most 13 orbits (all empty), 18 only shapes with 15 to 21 orbits, beyond the four-list join (4^{n_H − 1} per list; simulated annealing and
  CP-SAT do not even recover the 191 solution), and 45 none with at most 21.

## 89. Rounds 35–36: Paley plug-ins in the negacyclic ring; skew 2392 and 4712; the octic ring (October 23)

The commands that reproduce the results of §88 and §89, with their expected output, are in the sections
"Reproduction of search/NOTES.md §88–89" of [skewplug/README.md](skewplug/README.md) and
[negaskew/README.md](negaskew/README.md) (`patterns/`, `estimates/`, `spin/` and `octic/` hold the programs).

**The plug-in.** §88's quaternary family with H = H_K × S over G = K × F_q (q ≡ 3 mod 4, S the squares, H_K ≤ ℤ_|K|*
of odd order). The algebra of §15 holds with Gaussian-integer coefficients: a block is s ⊗ δ + t ⊗ Q + e ⊗ J, with
x₀ at (k, 0), x₁ = e + t on the squares and x₂ = e − t on the non-squares, s = x₀ − e; Q* = −Q and Q² = J − qδ give
U U* = (ss* + q tt*) δ + (ts* − st*) Q + (q ee* + se* + es* − tt*) J. The compressions w = x₀ + (q − 1)e satisfy
Σ w w* = 4|K| q δ (level 1), and t solves Σ(ss* + q tt*) = 4|K| q δ and Σ(ts* − st*) = 0 (level 2). Positions are of
three kinds: x₂ = −x₁ (e = 0, t = x₁ a free unit), full defects x₂ = x₁ (t = 0, cost (q − 1) + 2 Re(x₀ ē) ∈ {q + 1,
q − 1, q − 3}) and half defects x₂ = ±i x₁ (|e|² = 1/2, t = ±i e, cost (q − 1)/2 ± 1); the costs add up to 4|K|, and
the skew relation x₂(0) = i · conj x₁(0) forces a half defect in A at k = 0.
- *Engine* (`search/negaskew/qp.c`, `qshape.py`, `plugcert.py`). Level 1 joins the complete list of skew blocks and
  the list of all other blocks with total cost ≤ 4|K| and |ŵ(χ)|² ≤ 4|K|q at every character (a class-bucketed hash
  join on the exact compressed autocorrelations); level 2 is an exact meet in the middle over the t parts. All filters
  are necessary conditions. `real=1` restricts the alphabet to ±1: the real skew class of §87.
- *Validation* [EXACT]. `qp` equals `qskew` (§88, the plain orbit join over ℤ_v) on ℤ₃ × F₇ (79,164,864 normalised
  quadruples), ℤ₃ × F₁₁ (3,753,216), ℤ₇ with ⟨2⟩ × F₁₁ (3840) and on three empty shapes; in real mode it gives the
  exact counts of §87 (22,032, 13,824 and 9,360,000 normalised). An independent reimplementation of both levels
  returns at ℤ₁₃ with ⟨3⟩ × F₂₃ exactly the same 13,056 solutions as a set, and of level 1 the same 261,984 instances
  at ℤ₁₉ with ⟨7⟩ × F₃₁. A separate pattern-class pipeline found another solution of the 299 class.

**Hits.** v = 299 = 13 · 23, H_K = ⟨3⟩ (5 orbits on ℤ₁₃, 15 per block): 64 skew blocks, 2303 others, 36,000 level-1
instances, 13,056 solutions (exhaustive, 11 s) — [skew/2392](../skew/2392). v = 589 = 19 · 31, H_K = ⟨7⟩ (7 orbits on
ℤ₁₉, 21 per block): 261,984 level-1 instances, and on a 1% sample 84 completed (168 solutions) —
[skew/4712](../skew/4712). Table 5: 299(4) → 299(3), 589(5) → 589(3) (589(4) through [skew/9424](../skew/9424)). In both
certified solutions every defect is at k = 0: 299 has A's half defect (cost 10) and full defects of cost 20 and 22 in C
and D; 589 has 16 in A and 32 and 28 in B and D. The Gaussian estimate put 299 at 2^−1.3 ≈ 0.4 expected solutions:
estimates must not be used to discard plug-in shapes (HEURISTIC, as in §88).

**Other plug-in shapes** (quaternary, Table 5 targets with m ≥ 4).
- *EXHAUSTIVELY empty:* 539 (ℤ₄₉, q = 11, ⟨2⟩ of order 21: 10,856,160 level-1 instances, none completes), 667
  (ℤ₂₉, 23, ⟨7⟩), 707 (ℤ₁₀₁, 7, ⟨5⟩: 290,880 instances), 851 (ℤ₃₇, 23, ⟨7⟩), 899 (ℤ₂₉, 31, ⟨7⟩), 589 (ℤ₃₁, 19, ⟨2⟩),
  817 (ℤ₄₃, 19, ⟨4⟩ and ℤ₁₉, 43, ⟨7⟩), 989 (ℤ₄₃, 23, ⟨4⟩) and 413 (ℤ₇, 59, H_K = 1 and ⟨2⟩): all but 539 and 707 have
  no level-1 instance. 667 and 851 also die by the row sums in every defect pattern.
- *EXACTLY empty:* no defect pattern has cost 4|K| at 295, 335, 515, 655, 695, 995 (K = ℤ₅), 927 (ℤ₉), 611 (ℤ₁₃,
  q = 47); in the census beyond Table 5, 1363 and 2291 are EXHAUSTIVELY empty and 1403 dies by the
  row sums.
- *Sampled without a hit* (UNKNOWN): 737 and 869 (q = 11: every sampled level-1 instance fails the row-sum test of level
  2), 721 and 973 (q = 7), and with 9 orbits on K 943, 779, 759 (ℤ₃₃, 23) and 627 (ℤ₃₃, 19). Not run (lists of 10⁶ to
  10⁹): 575, 597, 623, 959, 627 (ℤ₅₇ and ℤ₂₀₉), 759 (ℤ₆₉ and ℤ₂₅₃), 803 and 979.
- *Origin-defect shapes.* The Table 5 targets where 4|K| is a sum of A's half-defect cost and at most three other
  origin costs are 299, 413, 517 and 589 (the others, 209, 253, 341, 345, 391, 473, 903, are settled in the
  repository).

**The real class** (`qp real=1`, skew order 4v). Exhaustive within the time limits: Table 5 targets with 3 to 10 orbits
on K, 85 shapes: 81 empty, 4 timeouts (623, 627, 759, 959 with 9 orbits); with 11 to 15 orbits: 517, 649, 737, 869
and 913 (K = ℤ₁₁, H_K = 1), 299 and 611 (ℤ₁₃, H_K = 1), 589 (ℤ₃₁, ⟨5⟩) and 851 (ℤ₃₇, ⟨10⟩) empty, the rest timeouts; census-undecided v
with 7 or 9 orbits on K, 139 shapes: 128 empty, 11 timeouts; with 11 to 15 orbits, 151 shapes: 65 empty, 86 timeouts.
No hit. (The plain real shapes with at most 5 orbits were swept in §88.)

**Census shapes of the 473 type** (real, a multiplier on K; `search/skewplug`). EXHAUSTIVELY empty: 627 (ℤ₃₃ × F₁₉,
⟨4⟩), 759 (ℤ₃₃ × F₂₃, ⟨4⟩), 1539 (ℤ₅₇ × GF(27), ⟨4⟩), 1953 (ℤ₆₃ × F₃₁, ⟨4, 22⟩), 2001 (ℤ₈₇ × F₂₃, ⟨7⟩), 2349 (ℤ₈₇ × GF(27),
⟨7⟩), 781 (ℤ₇₁ × F₁₁, ⟨20⟩) and 803 (ℤ₇₃ × F₁₁, ⟨2⟩); every pattern fails the big–big bound at 1083 (ℤ₅₇ × F₁₉, ⟨4⟩)
and 1705 (ℤ₅₅ × F₃₁, ⟨16⟩), and no slot has a row-sum pattern at 1199, 1397 and 2071. The raw level-1 counts of §88's sampled
shapes are mostly symmetric copies: each block may be reflected, and B, C, D negated and translated by the points fixed
by H_K, a group of order 2⁷ |Fix|³ that overstates distinct instances up to a thousandfold when |Fix| > 1, so those
zero results are weaker than they look. The best remaining real shapes are the three of 539 with Fix = {0}
(ℤ₄₉ × F₁₁ with ⟨18⟩, ℤ₇² × F₁₁ with ⟨(2, 2)⟩ and ⟨(2, 4)⟩), each run for under an hour; the second gave skew 2156 in
round 37 (§90). A rerun of the 473 shape with
symmetry-reduced instances gave a sixth solution of order 1892 (row sums 1, −21, −9, 37).

**Spin shapes** (C = B ∘ μ, D = B ∘ μ², 21 orbits). EXHAUSTIVELY empty: 737, 817 (two groups), 589 (two), 869, and
the census values 1273 (two), 1397, 1661, although the spin estimate gave each of them several bits of surplus
(`search/negaskew/spin/spin_est.py`).

**The octic ring** (order 16v). ℤ[x]/(x^{4v} + 1) ≅ ℤ[ζ₈][ℤ_v] for v odd; a ±1 sequence of length 4v becomes a
sequence on ℤ_v over 16 letters c₀ + c₁ζ + c₂ζ² + c₃ζ³, of squared moduli 4 ± 2√2, and Σ U U* = 16v δ has a peak
equation (the coefficient of √2 at shift 0) with no analogue in the order-8v family. Exhaustive counts of ordered skew
quadruples: 6144 (v = 1), 2,057,674,752 (v = 3), 75,652,141,056,000 (v = 5), with real matrices of orders 16, 48 and 80
verified. Quadruples made of two pairs (X, m · θX), with θ: x ↦ −x and m a monomial, are exactly the images of
order-8v quadruples (a, b, c, d) with b centre-palindromic (b_k = b_{2v−1−k}); they are 7.0% of the solutions at v = 3
and 0.044% at v = 5, and at ℤ₆₁ with |H| = 15 the
octic class has 15,360 normalised solutions where the order-8v family has none. For the Table 5 entries with m ≥ 5
the class is out of reach: a block constant on the orbits of H with |H| = h has row sum U(0) + hS with S of even
coordinates, so Σ|R_X|² = 16v needs (2h − 1)² + 15 ≤ 16v. Every multiplier shape with at most 14 orbits of the 61
targets is therefore infeasible, and the feasible ones (17 orbits at 337, 21 at 101 and 821, 23 at 947, …) need lists
of at least 2^57 blocks. No order 16v was found.

## 90. Round 37: skew 2156, the census value 539 (October 23)

**The hit.** §88's real skew Paley product with a multiplier on K, over K = ℤ₇ × ℤ₇ (non-cyclic) and q = 11, with
H_K = ⟨(2, 2)⟩ the scalars of order 3: 17 orbits on K, 51 per block, H of order 15. The identity of §15 reads
12 N_w − 4 N₋ = 196. The solution has no w-position in A and 6, 7, 7 in B, C, D (N_w = 20, N₋ = 11), row sums 1, −29,
15, 33 — [skew/2156](../skew/2156). 539 = 7² · 11 is undecided in the Đoković–Seberry census, so this settles a census
value; Table 5 had 539(4).

**How it was found.** Round 35 had ranked three K-symmetric shapes of 539 close to the 473 shape of §88 (ℤ₄₉ × F₁₁ with
⟨18⟩, ℤ₇² × F₁₁ with ⟨(2, 2)⟩ and ⟨(2, 4)⟩), but its slot weights came from Gaussian estimates of 16 sampled slots per
shape, and the ℤ₇² shapes ran for 25 minutes. `est_extend.py` estimates every one of the 2265 slots of the top energy
window in minutes; the best slots then scored logE = +0.94 (ℤ₇² with ⟨(2, 2)⟩), +0.81 (ℤ₇², ⟨(2, 4)⟩) and +0.66 (ℤ₄₉),
against +0.67 for the 473 shape with the same tool (whose two hit slots rank 6th and 40th of 161). With these weights
five seeded jobs of `prod_dedup.py` (level-1 instances reduced modulo the exact symmetries of the K-level system) found
the matrix after 35,188 distinct level-2 classes of the three shapes (13,334 from round 35), in round 7 of the job with
seed base 5, slot (0, 6, 7, 7) split (0, 3, 4, 4) (logE +0.79, third of 266). Per class the completion rate is about
five times lower than at 473, whose level-1 spaces are five times smaller (HEURISTIC). The logged state before the hit
replays the certificate exactly (`search/skewplug/data/repro539`, 11 s).

**Which shapes the sampler can reach** (HEURISTIC, calibrated on the hits at 473 and 539). A slot needs both
logE ≳ 0 (expected solution classes) and logρ = logE − logN1 ≳ −8.5 (expected completions per level-1 class); the
loop runs 10⁴ to 10⁵ classes per core-hour. Best slots (logE / logρ): 473: +0.67 / −7.0; 539: +0.94 / −7.9; 781
(ℤ₇₁ × F₁₁, ⟨5⟩): −1.9 / −7.6; 707 (ℤ₁₀₁ × F₇, ⟨36⟩): +0.14 / −16.5; 1243, 1267, 1397 and 637: logE ≤ −4.6. No other
census shape of round 35's list is in the reachable regime.

**An exhaustive real engine with exact symmetry reductions** (`search/negaskew/real37`, `qp37r.c`). `qp.c` with `real=1`
plus existence reductions of the K-level system: the units of K act on all four blocks at once (A canonical, `uA`); each
of B, C, D may be negated, reflected (k ↦ −k) and translated by a point of K fixed by H_K (`tF`: canonical = the
lexicographically least compression vector over this group, which replaces the normalisation x₀(0) = 1); C ↔ D (`sym`)
and class(B) ≤ class(C) ≤ class(D) (`s3`); also `dry=1` (exact level-1 probe totals) and `part=i nparts=n`. VALIDATED by
orbit coverage, the group orbits of the reduced solutions being exactly the full set (`symtest.py`, `symtest2.py`): ℤ₃ ×
F₇ 22,032 ← 436, ℤ₃ × F₁₁ 13,824 ← 384, ℤ₅ × F₂₃ 240,000 ← 2800 (80 with `tF`), ℤ₅ × F₁₉ 384,000 ← 4480 (128), ℤ₇⟨2⟩ ×
F₃₁ 576 ← 12, ℤ₁₃⟨3⟩ × F₇ 116,352 ← 706; the judge rebuilt all six, and `part` splitting sums to the unsplit count.
Effect: ℤ₁₃₇ × F₇ with ⟨16⟩ (9 orbits, estimated at over 30 minutes) takes 1.6 s, and 539 = ℤ₄₉ × F₁₁ with ⟨8⟩ (13
orbits, Fix = 7) needs 3.37 · 10¹⁰ probes (exact dry run, 276 skew blocks, 507,636 others) instead of 10¹³. Positive
controls (found in seconds, certified with `verify.py`): ℤ₁₁ × F₂₃ (skew 1012), ℤ₁₅ × F₁₉ (1140), ℤ₁₅ × F₂₃ (1380), ℤ₁₅
× F₃₁ (1860), ℤ₁₀₃ × F₃ with |H_K| = 17 (1236: q = 3, 526,344 level-1 instances), ℤ₅² × F₁₁ with a Singer subgroup of
order 3 (the class of [skew/1100](../skew/1100)).

**Verdicts of the real class** (EXHAUSTIVE for K, q, H_K as written, H = H_K × squares, blocks ±1, solutions up to the
symmetries above; `real37/data/RESULTS.txt`: 66 jobs, 65 EMPTY and 1 TIMEOUT, 201 = ℤ₆₇ × F₃, settled in §91; no
solution; the positive controls below are separate runs).
- The timeouts of §89 that the reductions make finite: 623 (ℤ₈₉ × F₇), 627 (ℤ₅₇ × F₁₁ and ℤ₂₀₉ × F₃), 759 (ℤ₆₉ × F₁₁ and
  ℤ₂₅₃ × F₃: 2.3 · 10⁸ level-1 instances), 781 (ℤ₇₁ × F₁₁, ⟨20⟩), 817, 897, 959, 1045, 1407, 1631, 1659 (ℤ₅₅₃ × F₃: 1.9
  · 10⁸), 1705, 1947, 1953, 2001, 2185 (ℤ₁₁₅ × F₁₉), 2277, 2365 and 2413; with 11 orbits and q = 3: 573, 1473, 1713,
  1893, 2433 (no compatible class pair, and independently by a row-sum argument: every H-orbit pairs with its negative,
  so A has row sum 1 and a free block has the signed sum of its orbit sizes, and 1 + R_B² + R_C² + R_D² = 4v has no
  solution; the same argument kills 927, 1719, 2151, 361, 529, 961 with ⟨7⟩, 1849 with ⟨9⟩ and 2209, EXACT).
- 539 = ℤ₄₉ × F₁₁ with ⟨8⟩: **EXHAUSTIVELY empty** (two halves of 25 CPU-minutes: 28,086 level-1 instances, every one
  through level 2, probes 1.687 · 10¹⁰ + 1.684 · 10¹⁰, no solution; the judge replicated both halves, 4.4 minutes each
  on a free core). So 2156 has no real skew matrix with this symmetry; the remaining 539 shapes (17 to 21 orbits) are
  sampled only. ℤ₇² × F₁₁ with the diagonal group of order 9 (9 orbits, agrees with the §88 sweep) and ℤ₇₇ × F₇ with ⟨4⟩
  are empty.
- Families never swept before (K cyclic with q | m, q = 27 or 343, non-cyclic K with diagonal H_K): 29 shapes with at
  most 11 orbits for 22 census values (361, 529, 783, 927, 961, 1029, 1083, 1269, 1323, 1539, 1587, 1593, 1647, 1719,
  1849, 1917, 2107, 2151, 2209, 2303, 2401, 539) and 8 shapes with 13 to 15 orbits and q ≥ 27: all EMPTY, seconds each.
- K = ℤ_p² with a matrix H_K (a Singer subgroup): 575 = ℤ₅² × F₂₃ (order 3: no compatible class pair), 1183 = ℤ₁₃² × F₇
  (order 21: 5806 level-1 instances, every one fails the DC test of level 2), 1083 = ℤ₁₉² × F₃ (order 45: no compatible
  class pair).
- UNKNOWN (time, not memory): 917, 1057, 1477, 2317 (q = 7, 11 orbits: 6 · 10⁹ to 2 · 10¹⁰ probes), 2185 = ℤ₉₅ × F₂₃,
  the rigid shapes with 15 orbits whose level 2 needs 2¹⁴ to 2¹⁵ sign fibres per block (2485, 1645, 1209, 2109, 2093,
  805), and every shape with 13 or more orbits and q = 3, 7, 11 except 539 with ⟨8⟩ (lists of 10⁹ to 10¹¹). 201 = ℤ₆₇ ×
  F₃ (7 orbits) was UNKNOWN here (78% of its probes, no solution) and is settled in §91.
- *A defect of `qshape_mat.py`* (found by the judge, fixed in `qshape37m.py`): it wrote the unit action of `uA=1` for
  all diagonal pairs (u₁, u₂) of ℤ_p², but diag(u₁, u₂) permutes the H_K-orbits only if it normalises H_K, which fails
  for a Singer group when u₁ ≠ u₂. The three matrix verdicts above were re-run unreduced (575: |A| = 16, no instance;
  1183: 1,999,872 level-1 instances, all fail; 1083: |A| = 1116, |L| = 686,754, no instance) and stand; at 1083 the
  reason is 'no compatible class pair', not 'empty skew list'. No other shape is affected (`unitcheck.py` on 169 shape
  files).

**The quaternary plug-in beyond Table 5** (`search/negaskew/quat37`, `qp37q.c`: `qp.c` with `usym` (A canonical under
the units of K and the conjugation), `tsym` (B, C, D under translation by the fixed points of H_K), `sym`, exact
weights). Shapes: K = ℤ_m, q ≡ 3 (mod 4) prime with q ∤ m and q ≤ 8m + 3, H_K every subgroup of odd order, 3 to 9
orbits, for the census-undecided v in (999, 2499] with 4v not in the repository: 452 values, 204 shapes, 130 realisable
by the defect identity (46, 25, 36, 23 with 3, 5, 7, 9 orbits). VALIDATED: weighted counts equal the unreduced ones (ℤ₃
× F₇ 79,164,864 with and without `tsym`, ℤ₇⟨2⟩ × F₁₁ 3840, 299: 13,056 with 36,000 level-1 instances; the judge also).
EXHAUSTIVELY empty: all 46 with 3 orbits (42 by an empty skew list); all 25 with 5 orbits (including the origin-defect
shapes 1363, 1711, 2291 = ℤ₂₉ × F₄₇, F₅₉, F₇₉, which have 696, 540, 576 level-1 instances and no level-2 solution, and
1903, 1603); 11 of the 36 with 7 orbits (1273, 1333 over ℤ₄₃ × F₃₁, 1397, 1529, 1661, 1793, 1817, 1917, 2077, 2369,
2449); with 9 orbits 1271 and 1763; the prime-power q = 27 shapes 1269, 1323, 1593, 1647, 1917. **943 = 41 · 23, H_K =
⟨10⟩ of order 5 (Table 5: 943(4)): EXHAUSTIVELY empty** — 6.0047 · 10¹⁰ probes (class-level estimate and measurement
agree to 10⁻⁵), no level-1 instance, one hour of CPU. UNKNOWN: 1333 over ℤ₃₁ × F₄₃, 1133 (stopped at class 768 of 1911),
869, 2189, the 22 shapes with q = 7, 3 and 7 orbits and the 21 census shapes with 9 orbits (and 779, 759, 627 among the
targets): lists of 10⁷ and more entries. Of the empty shapes' v, 17 have 8v − 1 a prime power (Paley I), so their
verdicts are not about new orders. The judge re-ran 35 of 36 cheap shapes unreduced with the engine of §89 (all agree;
1661 timed out). A hit would have been a skew order 8v with the symmetry; none was found.

**The origin-defect shapes** (`search/negaskew/origin`; quaternary class with every non-origin position of type E).
Writing a block as s ⊗ δ + t ⊗ Q + e ⊗ J, e_X = ε_X δ₀ (ε = 0, a unit, or (1 ± i)/2 times a unit). EXACT (derived from
Q* = −Q, Q² = J − qδ, QJ = 0; checked to 10⁻¹² on 1005 + 220 + 68 solutions and re-derived by the judge): Σ(ss* + q tt*)
= 4mq δ (d), Σ(ts* − st*) = 0 (Q), Σ(q ee* + se* + es* − tt*) = 0 (J); with n/2 = Σ|ε_X|², c_s = Σ|s|² = qn/2, c_t = 4m
− n/2, so q = 2c_s/n is fixed by the defect pattern. PROVED: at a lag k ≠ 0 Σ_X P_{s_X} = −q Σ_X Φ_X with Φ_X = ε_X conj
u_X(k) + conj ε_X u_X(−k) ∈ ℤ[i], |Σ P_s| ≤ 4(m − 2) + 2Σ|S_X|, so q above this bound forces ΣΦ = 0 (the 'Lock': for the
rigid pair (A, M) u_B = λ u_A off the origin with λ = −2 conj(ε_A) ε_B; the judge found λ ∈ {1, −i} in the formal-q
sets); the row sums satisfy Σ(|R^s|² − |S|²) = 4(m − 1) + qZ, Σ(|R^t|² − |τ|²) = 4(m − 1) − Z with Z ∈ ℤ, which for Z =
0 gives h ≤ 4p + 5 (h = |H_K|, p = number of pairs ±O of orbits); every K-level solution with q = 2c_s/n a prime power ≡
3 (mod 4) gives a skew matrix of order 8mq. Shapes with a single big defect (n = 1: q = 8m − 1 or 8m + 3) are EXACTLY
empty for every K and H_K (413, 2363).
- Members, EXHAUSTIVE per shape in the restricted class (`d1=1`, or `d1=2` with H_K trivial): (m, q) = (3, 7), (3, 11)
  (3,753,216 solutions), (5, 7), (5, 11), (7, 11), (7, 19), (7, 27), (7, 31) with H_K = 1; (13, 19), (13, 23), (13, 27)
  with ⟨3⟩ (1536, 13,056 and 2304 solutions; the (13, 27) set checked at G level over ℤ₁₃ × F₂₇, every eighth solution);
  (19, 31) with ⟨7⟩ (589, about 17,000 solutions, sampled). The certified matrices of orders 280, 616, 1736 and 1976 (8
  · 35, 77, 217, 247) are not new (those v are census-decided; the judge rebuilt them). 299 has 136 classes of 96 under
  the symmetry group (multiplier ℤ₁₃*/H_K, x₁ ↔ x₂, conjugation, reflection, S₃ on B, C, D).
- Empty (EXHAUSTIVE per shape in the class): every shape with 3 orbits (109 shapes, m up to 383) except (7, 11), whose
  3840 solutions have defects off the origin; the shapes with 5 orbits at m = 29, 37, 49, 53, 61, 101, 109, 121, 149,
  157 and 13 × 103, 13 × 107; with 7 orbits 19 × 79, 27 × 71, 31 × 43, 43 × 59, 43 × 71; the prime powers 11, 19, 23, 29
  × 27, 61 to 211 × 243 and 43 × 343. The judge's own enumeration of the nu = 3 identity-feasible pairs found 225 (the
  109 are a subset: q = 3, 7, 11 with many defects were not tested). The judge noted that 'EXHAUSTIVE for nu = 3' is
  therefore for the shapes tested only.
- **517 = 11 · 47** (Table 5: 517(4)): the quaternary Paley plug-in class over ℤ₁₁ × F₄₇ is EXHAUSTIVELY empty **for
  every multiplier group H_K** (H_K = 1 included). The defect budget 4m = 44 leaves A's origin half defect of cost 22
  and one more half defect of cost 22 (hand proof, re-derived by the judge), which can be translated to the origin:
  pattern (A, M | E, E); the decoupling bound 46.32 < q = 47 gives the Lock and a K-level equation 2P_C + 2P_D =
  −(P_{s_A} + P_{s_B}) that `rigid_s.py` solves exactly (all A, C over all 4¹⁰ sequences, D modulo translation and unit:
  95,326 classes): no solution in any of the four origin variants. VALIDATED: three independent implementations
  (`rigid_s.py`; `rigid_s_check.py`: FFT and full D, no linear keys, all 1024 A of one variant in 1831 s; the judge's
  `rigid_j.py`), a planted-solution test (5/5), and the Lock counts m = 3: 144, m = 7: 30,576 per origin variant (m = 5,
  9, 11: 0; the judge's `rigid_j.py` re-run for the port gives 144, 0, 30,576 and, at (11, 47), 0 for all 1024 A of each
  of the four variants, 7 min CPU; at m = 3 it also counts two variants with a full defect on B, outside the pattern (A,
  M | E, E)). The class is empty at ℤ₁₁ × F₄₇ whatever H_K; other multiplier groups on ℤ₅₁₇ are not covered.
- HEURISTIC: a Gaussian model E(m, p) ≈ (3.4 · 10⁴/m³)ᵖ per origin variant, calibrated on 299, 589, 247 (and a held-out
  check at (13, 27): 24 classes against the predicted 80), predicts no members for m ≳ 30: the family is expected to be
  finite. The model fails for rigid shapes (at m = 11 it predicts 3 · 10³ to 4 · 10⁴ (Lock + S) solutions, the exact
  count is 0). The uniform cyclotomic conjecture (class values independent of m) is FALSIFIED (EXACT numerics: the
  13,056 solutions at (13, 23) extended to m' = 29, 37, 53, 61, 101, 109 and the 589 sample to 31, 43, 67, 79, 103 leave
  a residual that grows like √m). A CP-SAT model of the K-level system does not find the 299 solution.
- The real hit 2156 is not an origin-defect solution: its defect budget 4|K| = 196 is spread over 20 positions of 6
  orbits (cost 12 on 9 and 8 on 11 positions), the other regime of the same K-level system (many small defects, small q,
  large |K|).

**What the judge checked** (the rebuild of 2156 and the corrections of the production statistics and of the hit slot's
rank are already in §90 and skew/2156/README.md). Not re-run in full by the judge: the verdicts with 13 to 15 orbits
(its own builds later reproduced most: 627/ℤ₂₀₉, 759/ℤ₂₅₃, 1045, 1407, 1659, 1947, 2185/ℤ₁₁₅, 2277, then 817, 897, 1323,
1705, 1953, 2001, 2365, 2413 and the q = 3 ones). The judge's `unitcheck.py` (169 shape files) and `dc_check.py` are in
`real37`; his `rigid_j.py` is in `origin`.

## 91. Round 38: the full-slot re-scan of the census shapes; skew 1540, the census value 385 (October 23)

Scratchpad `r38/`.

**The re-scan** (HEURISTIC, the exact-dead claims pending a judge). 641 undecided census values give 339 pairs (v, q)
with q ≡ 3 (mod 4) a prime power and a solution of the coarse identity, and 554 K-symmetric real skew shapes: K every
abelian group of order |K|, H_K the odd-order subgroups of Aut(K) up to conjugacy for cyclic K and ℤ_p², torus and
Singer subgroups for ℤ_pᵏ (k ≥ 3), diagonal units otherwise, at most 25 orbits on K. 249 are decided by the exhaustive
verdicts of rounds 34–37; 9 die by exact arguments (the lemmas for N_w = 1 and 2 at 1139, 1207 and 1081, the big-big
filter of §87 at 361 and 529 with H_K = 1, the row sums over all energies at 1539, 1603, 1919 and 2401); 241 were
scanned with `est_extend.py` (44 over every slot of the window, 197 sampled). None meets the criterion of §90 with a
payable level 2. The best margins are those of 575 = 25 · 23 with H_K = 1 (over ℤ₂₅: logE +0.37, logρ −7.31; over ℤ₅²:
+0.17, −7.46; both with 2²⁴ sign fibres per block at level 2), then 245 = ℤ₃₅ × F₇ and 385 = ℤ₃₅ × F₁₁ with H_K = ⟨11⟩,
and 781. The judge found that the scan's row-sum cache was keyed by object identity, so a worker that ran several
shapes marked some feasible slots infeasible (575, 627, 957, 385 over ℤ₅₅ × F₇); refilling them moves no margin above
these. Its estimate files cover only the top energy window E ≥ E_max − 5 for most shapes (for 245, E = 30 to 35 of
18 to 35); only 385, 575, 1083, 1159, 1645 and 2485 were estimated over the whole range. 1793 and 2189 have no slot at
any energy (the orbit weights cannot meet the identity); 917 has 161 feasible slots outside the window. The q = 3 shape
with the strongest signal, 201 = ℤ₆₇ × F₃ with H_K of order 11 (logE +2.01), has no solution in the round-37 run
resumed in round 38 (2.3 · 10⁹ level-1 and 2.5 · 10⁸ level-2 instances, the same engine; not rerun independently), so
the estimate overstates q = 3 shapes by at least two decades. Not enumerated: subgroups of Aut(ℤ_pᵏ) with unipotent
elements and non-diagonal automorphisms of non-elementary p-groups (33 shapes flagged).

**An enumerating level 1.** On 385 and 245, `prod_dedup.py` ran about 1650 classes in 8 minutes without a completion:
each class of level-1 instances recurred about 2300 times, and 46% of the rounds drew a pattern with no instance.
`exh38.py` generates, per block, defect weight, m-weight and row sum, the block list until it covers that space, keeps
one block per orbit of the exact block symmetries (those of `prod_dedup.py`'s keys), reduces the skew block by the
units of K, and joins every row-sum pattern of a slot with `gl1.c`; the new classes go to the exact level 2. It runs
about 100 classes per second per core. Checked: on the slot (2, 2, 4, 7) / (2, 1, 3, 4) of 385 it finds 6184 classes,
among them every class of the two sampler runs, and the unit-reduced run gives 1546 classes with the same unit orbits.

**The hit.** Swept in the order of the estimated completions per class, the seventh slot of 385, allocation (6, 1, 1, 4)
with m-split (0, 0, 0, 1) (N_w = 12, N₋ = 1, so 12 · 12 − 4 = 140 = 4 · 35), holds 100 classes, and 2 complete, 45 s
after the launch: two skew-Hadamard matrices of order 1540, row sums 1, −7, 31, 23 ([skew/1540](../skew/1540)). 385 = 5
· 7 · 11 is undecided in the Đoković–Seberry census; Table 5 had 385(3). The slot is the lowest energy of the shape (E =
12 of 12 to 17), and its completion rate, 2 in 100 classes, is five decades above its estimate (logρ −7.49; logE −1.54,
1.6 decades below the 539 hit and outside the criterion of §90). For an enumerating level 1 the relevant quantity is the
number of distinct classes of a slot, 10² to 10⁴ for the leading slots of these shapes, and level 2 (3 ms per class) is
not the bottleneck. The second matrix has the same row sums but a different level-1 class (C's defect at 7 instead of
21, D's at 1, 11, 16, 21); its equivalence to the first was not tested. An independent judge rebuilt both from their
certificates and replayed the slot. The sweep of the hit slot replays in 14 s (`search/skewplug/data/repro385`).

**The other shapes** (HEURISTIC; a finished slot covers 60 to 90% of each raw block list, so no verdict is exhaustive).
About 22,700 classes reached level 2 in total: 245 (11,447 classes, 154 of 2247 slots), 1477 = ℤ₂₁₁ × F₇ with ⟨101⟩ and
⟨100⟩ (4036 and 1251), 1267 = ℤ₁₈₁ × F₇ (2277), 1057 = ℤ₁₅₁ × F₇ (885), 781 = ℤ₇₁ × F₁₁ (762), 1645, 2485, 605 and 385 =
ℤ₅₅ × F₇ (a few each), and no completion outside the hit slot; the 29 lowest-energy slots of 245 gave 5331 classes and
none. By the summed Gaussian estimates 245 holds the most solutions (about 10⁴, HEURISTIC), and the lowest-energy slots
of every shape, which the estimates did not cover, are the natural continuation.

**Details of the re-scan and the production loop.**

- *The enumeration* (`scan38/enum38.py`, 106 s, 552 shapes when re-run now because 385 has since been settled; 554 at
  the scan, the annotated list is kept in `scan38/data/shapes38.json`): of the 641 undecided values, 235 have at least
  one K-symmetric shape, 129 of them still an undecided one, 106 only decided or EXACT-dead shapes, and 406 values have
  no shape in this family. `decided38.py` collects 642 shape records of earlier exhaustive verdicts (r36 sweep 452, r34
  80, r37 64, r35 46); `prep38.py` annotates: 296 decided (47 of them also dead by the orbit weights), 258 undecided.
- *Calibration of the reachability criterion* (HEURISTIC). The solved shapes with H_K = 1 span a wide range of the best
  slot (logE / logρ): 253: +5.78 / −0.57, 285: +6.81 / −2.11, 345: +6.30 / −2.33, 465: +5.08 / −1.95, 391 = ℤ₁₇ × F₂₃:
  +6.12 / −3.01, 275 = ℤ₂₅ × F₁₁: +8.74 / −8.68 (the class of skew/1100), 119 = ℤ₁₇ × F₇ (the 476 hit): +12.35 / −5.45
  (logN1 17.8); the weakest hits, 473 and 539, were found by a sampler after 5 · 10³ to 5 · 10⁴ classes with the
  estimator 3 to 4 decades low there. The criterion logE ≳ 0, logρ ≳ −8.5 is the lower edge of this range. A third
  quantity is the level-2 cost, tbits = max over B, C, D of ν(1 − nₓ/|K|), the logarithm of the number of t-sign fibres
  per block (hits: 473 12.9, 539 15.6; the exact level 2 handles about 2¹⁵ in seconds, and 2⁴⁴ pair evaluations per
  instance at H_K = 1 and |K| = 25): 575 = 23 · 25 with H_K = 1 meets the first two (+0.37 / −7.31 over ℤ₂₅; ℤ₅²: +0.17
  / −7.46 after the correction of the judge) but has tbits = 24, so level 2 is out of reach (r34: UNKNOWN; the big–big
  filter does not kill it).
- *Coverage bias of a sample* (measured by bootstrap on the complete scans 539, 781, 1243): a uniform sample of a
  fraction c of the slots understates the best logE by (median / 90th percentile) 0.12 / 0.3 decades at c = 0.3, 0.3 /
  0.95 at 0.1, 0.8 / 2.0 at 0.03, 1.15 / 2.4 at 0.01 (reproduced by the judge: 0.08 / 0.31, 0.31 / 0.89, 0.78 / 1.82,
  1.14 / 2.49); `table38.py` reports the margin min(logE, logρ + 8.5) of the best slot, corrected by these shortfalls
  (adjMed, adjP90). For |K| > 150 the windows are sampled lazily and the estimate is degenerate (logE ≈ −5 · 10⁴): those
  shapes are heuristically dead but not enumerated.
- *Estimator biases.* The Gaussian estimate knows nothing of the exact lemmas: 1207 = ℤ₁₇ × F₇₁ (one w-position, logE
  +2.06), 1139 (N_w = 1 lemma of §15), 1081 (N_w = 2 lemma of §88) and 361 = ℤ₁₉ × F₁₉ (+4.15 / −4.0, killed by the
  big–big filter) score high and are dead; shapes with N_w ≤ 8 and large q need the filter, not the estimate. For q = 3
  the identity gives no upper bound on E, so the top energy window is the all-defect corner; the window is E in [N, 1.6
  N] (calibrated on the solvable 309 = ℤ₁₀₃ × F₃: 68,866 slots, 426 row-sum feasible, best +1.34 / −5.67, top-40 sum
  219), and the 13 decided-empty q = 3 shapes of §90 show no feasible slot, but the estimator is no discriminator there:
  201 scored +2.01 / −4.98, better than the control, and is EXHAUSTIVELY empty (two runs with identical lists, rows 0 to
  19 and 19 to 53 of 54, 2.28 · 10⁹ level-1 and 2.5 · 10⁸ level-2 instances).
- *The exact all-energy row-sum scan* (`rowscan38.py`: every allocation of every energy, every m-split; ALL_INFEASIBLE =
  EXACT-dead): 1183 (79,895 slots), 1539 = ℤ₃⁴ × F₁₉ with ⟨5⟩ (92), 1603 = ℤ₂₂₉ × F₇ with ⟨19⟩ (1,114,897), 1859 and
  2375 (Singer subgroups), 1919 = ℤ₁₀₁ × F₁₉ with ⟨5⟩ (388), 2401 = ℤ₇³ × F₇ with order 57 (75,747) and 171. Alive at
  other energies: 917 = ℤ₁₃₁ × F₇ with ⟨13⟩ (161 feasible slots, the first at E = 118). 1793 and 2189 have no slot at
  any energy (the orbit weights cannot realise the identity; the first scan listed them as alive by admitting m-splits
  that are not sums of orbit weights); 2443 = ℤ₃₄₉ × F₇: 2.46 · 10⁶ slots, none feasible (PARTIAL).
- *The production loop* (`exh38.py`; `scan38/rotate.py` rotates it over shapes with a time slice proportional to the
  measured yield of new classes per second, `supervise.sh`, `mklow.py` writes the low-energy slot files). Per shape,
  distinct classes sent to level 2 / completions: 385: 2015 (35 of 2747 slots, in the order of estimated completions per
  class) plus about 4000 of the sampler, **2** (the hit); 245 = ℤ₃₅ × F₇: 13,814 in 168 slots of the top window and 5331
  in the 29 lowest-energy slots (E = 18), 0; 1477 ⟨101⟩ 4036 and ⟨100⟩ 1251, 1267 2277, 1057 885, 781 762 + 902 low
  energy, 1645 39, 2485 4, 0 each. 605 and 385 over ℤ₅₅ × F₇ sent 2 and 4 classes; 627, 737, 803, 1141, 1243 and 1393
  sent none (the skew-block lists exceed `PROD_MAXA`, or every join exceeds `EX_JOINCAP`): UNKNOWN. The slots marked
  finished cover 60 to 90% of each raw block list (`EX_COV`), so no verdict is EXHAUSTIVE. The measured facts that
  motivate the next sweep (HEURISTIC): the level-1 spaces are tiny after symmetry (a class of level-1 instances appeared
  about 2300 times in the sampler; the number of classes of a slot is N1/(2(4|Fix|)³ times the unit orbits), 10² to 10⁴
  for the leading slots), so enumeration beats sampling 50- to 400-fold and level 2 (3 ms per class) is never the
  bottleneck; the Gaussian estimate is no guide to completions (the hit slot has logE −1.54, logρ −7.49 and 2
  completions in 91 classes, none of 1900 other classes of the first 34 slots of 385); the scan's top-window files leave
  the low energies of most shapes unestimated, and the hit sits at the lowest energy.
- *Conventions of the port.* The shape indices of the tables (`scan38/data/TABLE38.txt`, `SUMMARY38.txt`,
  `rowscan38.jsonl`, the `S<idx>_…` files of `est/`) refer to the annotated list `scan38/data/shapes38.json` (554
  records); `enum38.py` re-run now gives 552 (385 is settled), with other indices, so every command that takes an index
  takes that file. Two defects of the code of round 38 were found and fixed while porting: `rowfast2.py` cached block
  classes under `id(sh)` and reused them across shapes of one worker process (51 slots of a control: 1 mismatch against
  the exact row patterns before, 0 after; the estimates of 575 over ℤ₅² moved from −0.47 / −7.38 to +0.17 / −7.46, no
  margin moved above the others), and the control of 539 now runs from `scan38/ctl539.py` (22 s, 2265 slots, 266
  feasible, best +0.94 / −7.87 at alloc (6, 4, 6, 7), the hit slot third at +0.79 / −7.75, top-40 sum 107.2; HEURISTIC)
  because `est_extend.py` of the repository runs 29 slots of this shape in 240 s.
- *Checks by the judge* (1540; see §91): independent rebuild (H Hᵀ = 1540 I, H + Hᵀ = 2I, row sums 1, −7, 31, 23), two
  bit-exact replays (12 s), the estimates of 781 (all 1735 slots), 1645, 1057, 2485, 1267, 385, 245, 575 and the
  controls reproduced with `est_extend.py`, the nine exact-dead shapes rederived (1139, 1207, 1081 by the lemmas; 361,
  529 by `bbhk1.py`; 1539, 1603, 1919, 2401 by row sums), 1207 also by the pattern pipeline (EXHAUSTIVELY empty, 102 s).
  Not rerun: 201 (100 CPU-minutes), 2443.

## 92. Round 39: the 2026 literature rechecked, the triangular family, the orders 8p (October 23–24)

Scratchpad `r39/`.

**The literature.** The survey of Egan, Kotsireas, London and Ó Catháin (arXiv:2610.27175, 19 October 2026) credits
this repository (commit 2d45dbb of 18 October) with 2148, 2284, 2636, 3964, the thirteen Paley-block orders and 6108, and
states that 3732, 4492, 5020, 5788 and 5932 "were found by the present authors", by Goethals–Seidel arrays with
multiplier groups of orders 31, 33, 50, 39 (5788 = 4(4 · 19² + 3) by its Theorem 7 over ℤ₂² × F₁₉²). Cooley's
Middletownbooks/Hadamard committed 2092 and 3732 on 14 October. Ours of 2092, 4492 and 5932 were committed on 11 October
and 3732 and 5020 on 13 October; the five are relabelled Concurrent (our 3732 is skew, neither other one is stated to
be), 1447 leaves [FRONTIER.md](../FRONTIER.md), and the survey's six open orders 8p (p = 599, 719, 739, 839, 1031, 1087)
are recorded there ([NOVELTY.md](../NOVELTY.md)). Its 95 open orders up to 10 000 agree with FRONTIER.md except for our
2572, 3436, 4204 and 6796, which it predates.

**Spence's construction for 2396 with H_K = 1.** §86's census estimated the Spence plug-in over ℤ₂ × ℤ₁₃ × F₂₃ (the
group of order 598 = n − 1 for n = 599) with trivial H_K at logE 16.3 and left it untouched because its level 2 has 99
free signs. It is EXACTLY empty. Identity (I) at 0, Σ(s² + 23 t²) = 4 · 598 + 4 = 2396 with 24 per z-position and s² ∈
{0, 4} at a w-position, gives 24 N_w − 4 N₂ = 100 with N₂ ≤ N_w, so five w-positions, each with s = −2e and compression
u = 21e = ±21. Two of them lie in one of the four blocks, at a distance d; at the shift d their product contributes ±441
to Σ PAF_u, which must be −92, while all other terms there are bounded by 8 · 21 (big against small) plus 95 (small
against small): |−92 ∓ 441| ≥ 349 > 263.

**The triangular family** (§88: v = 1 + k(k + 1)/2 prime, k ≡ 3 (mod 4), H ≤ ℤ_v* of order k, e = (k + 1)/2 cosets;
tools and commands in [negaskew/tri/README.md](negaskew/tri/README.md), proofs P1 to P7 in
[negaskew/tri/PROOFS.md](negaskew/tri/PROOFS.md)).
- *Row sums* (PROVED). Each pair of cosets ±C_j contributes ±(1 + i) to the skew block's coset sum, so S_A = m(1 + i)
  with m ≡ e/2 (mod 2); the identity 4v − 3 = k² + (k + 1)² is met by one block with S = ±1 ± i and three with S = 0,
  and for k ≥ 7 this is the only pattern. So A is the big block iff k ≡ 3 (mod 8) (67, 191, 379, 631, 947, 1327), a free
  block iff k ≡ 7 (mod 8) (29, 277, 1129, 2017).
- *The skew spectrum* (PROVED, [negaskew/tri/PROOFS.md](negaskew/tri/PROOFS.md) P4). At every t, t = 0 included, the
  transform U of the twisted-skew block satisfies conj U = (1 + i) − iU, so U = 1 + (1 + i)r with r real (for u(0) = 1
  and for u(0) = −i) and |U|² = 2r² + 2r + 1: 2|U|² − 1 = (2r + 1)² is the square of a real element of ℚ(i, ζ_v) =
  ℚ(ζ_{4v}), whose only real quadratic subfield is ℚ(√v). Three flat free blocks (|U|² = v + 1 at every t ≠ 0, PAF −1
  off the origin, e.g. three Legendre sequences when v ≡ 3 mod 4) would give, for k ≡ 3 (mod 8), |U_A(t)|² = 4v − 3(v +
  1) = v − 3, that is (2r + 1)² = 2v − 7 = k² + k − 5, which is neither a square ((2k + 1)² − 4n² = 21 only for k = 2,
  5) nor v times a square (v | 7) for k ≥ 7; for k ≡ 7 (mod 8) the big free block is flat at v − 3, |U_A|² = v + 1 and
  (2r + 1)² = 2v + 1, excluded below; so no solution has three flat free blocks. No skew block of a solution is flat
  either: its flat value would be v − 3 (|R_A|² = 4v − 3) or v + 1 (|R_A|² = 1), and 2v + 1 is neither a square nor v
  times a square. The exception is v = 7 (2v − 7 = v). A single flat free block does occur, at v = 29 and 67 (960 of
  5952 and 576 of 6192 solutions), and not among the solutions found at 191 (all 99,360), 277 (all 646 reduced) and 379
  (the sample of 448).
- *Exact counts* (normalised: B, C, D with u(0) = 1). v = 29: 5952; 67: 6192; 191: 99,360 (105 classes under the
  symmetry group ℤ_e × tilde³ × S₃ × (A ↦ −A) of order 96e; `qskew`, `tri_join` and `tri2` agree; `tri2` needs 7
  CPU-seconds at 191, 1.8 ns per probe against about 70 ns for `tri_join`); 277: 359,424 (EXHAUSTIVE in three complete
  runs of `tri2` ([negaskew/tri/tri2.c](negaskew/tri/tri2.c)): two shift-reduced with different moduli and seeds, one
  unreduced (7.29 · 10¹¹ probes), logs in `negaskew/tri/data`; the reduced run was repeated (4.6 CPU-minutes) and is
  byte-identical to its log; and in the judge's own unreduced join; 313 classes: 311 of 1152 solutions and 2 of 576).
- *v = 379 has solutions* (EXACT): `tri_join 379 27 2 262144 170258 170259 1000` finds 32 in 2 minutes, the images of
  one solution under the tilde of B, C, D, the swap of C and D and A ↦ −A; `tri_cert.py` gives a skew-Hadamard matrix of
  order 3032 = 8 · 379 (verify.py and an independent numpy check; not a new order, since 1516 = 4 · 379 is known). My
  partial run of 14 even residues of M = 2¹⁸ found nothing because a solution without a stabiliser brings 31 others into
  its residue (the tilde of B, C, D, the swap of C and D and A ↦ −A keep both pair sums of the join: a group of order
  32, PROOFS.md P5; each of the 14 hit partitions of a random sample of 1484 holds exactly 32): E/(32 · 2¹⁷) · 14 = 0.3
  to 0.5 independent clusters were expected there (E = 9.5 · 10⁴ to 1.6 · 10⁵), 1.1 with the Gaussian estimate. At v =
  277 (B the big free block, S = 1 + i) the group has order 16 and the hits per partition are multiples of 8, 16 for
  2322 of the 2996 hit partitions. HEURISTIC: about 10⁵ normalised solutions at 379 (9.5 · 10⁴ ± 35% from the sample of
  partitions; about 70 inequivalent solutions) against the Gaussian 3.3 · 10⁵ (ratio 0.28; 0.77 at 29, 0.97 at 67, 0.49
  at 191, 0.41 at 277). An earlier pooled value, 1.25 · 10⁵ (ratio 0.38), is inside the error of the logs but none of
  the three ways of pooling them (plain, inverse-variance, unreduced plus reduced) reproduces it.
- *No explicit structure* (EXACT on the solution sets of 29, 67, 191, and 277, 379 where possible): no two small blocks
  equal up to unit and tilde beyond v = 29, no spin or shift relation beyond e = 6, no sparse coset DFT, no flat or
  Legendre-type free block beyond 67, no character-type block (c_j = i^(aj+b), `charblk.py`) beyond 67 (at 29 a = 1 and
  a = 2 occur, at 67 a = 2; none at 191, 277, 379).
- *Open members* (UNKNOWN). With the measured ratio 0.28 the expected numbers of inequivalent solutions are 17 (631),
  6.6 (947, Table 5's 947(6)), 14 (1129), 0.9 (1327), 0.12 (2017) (HEURISTIC; `est_fast.py`, E ρ/(96e)): the family is
  probably finite. With random sub-lists of each free list (fraction f, E ρ f³ = 1), the four-list join needs about 6 ·
  10¹⁵ probes at 631 and between 4 · 10²⁰ and 10²² at 947, depending on the memory for the sub-lists (the plain scan
  costs |L|² = 2.5 · 10²³ there). No construction for the larger members is known.

**The orders 8p** (p = 599, 719, 739, 839, 1031, 1087: the survey's six open orders 4792, 5752, 5912, 6712, 8248, 8696,
n = 2p, so Goethals–Seidel arrays over groups of order 2p and Spence's bordered array over groups of order 2p − 1 =
1197, 1437, 1477, 1677, 2061, 2173; programs, data and every command in [eightp/README.md](eightp/README.md)). No matrix
was found and none is claimed.
- *Plain structures* (EXACT kills, HEURISTIC unreachability). Real Goethals–Seidel over ℤ₂ₚ and the negacyclic array
  over ℤ_p (quaternary blocks) with a multiplier group H ≤ ℤ_p* must meet the ℤ₂-compression identities Σ(s₀² + s₁²) =
  8p, Σ s₀s₁ = 0 with s_a ≡ ±1 (mod 2|H|), and Σ|u₀ + |H|S|² = 4p. Every group is dead except |H| = 2 (every p), 13 at p
  = 599 (quaternary only: 47 orbits), 3, 6 and 9 at 739, 5 and 10 at 1031 (10: quaternary only), 3 and 6 at 1087: every
  live pool has at least 2⁹⁴ sign vectors per block (2¹⁶⁶ in the real array). (The first census listed the odd orders
  only; the judge's census of all orders agrees on the rest.)
- *The theorems* (EXACT from the survey's hypotheses). Paley I and II (8p − 1 is never a prime power), Paley blocks (4p
  − 1 = 5 · 479, 5³ · 23, 3 · 5 · 197, 5 · 11 · 61, 7 · 19 · 31, 3³ · 7 · 23; (s, t) = (2, p) needs T-matrices of order
  ≥ 599), Miyamoto III (2p − 3 = 1195, 1435, 1475, 1675, 2059, 2171 are composite), Seberry–Yamada 8q (WT(299), WT(359),
  WT(369), WT(419) and WT(515) are in no table; WT(543) is, with key wi, and the survey excludes 543 = 3 · 181
  explicitly, but that route is valid and gives 8696: §95), Cooper–Wallis, T-sequences, Spence, Mathon, Yamada and
  Agaian: none applies. The Craigen–Holzmann–Kharaghani decomposition 8(x + y), x + y = p, would exist at five of the
  six p if the complex Golay lengths were closed under products; they are not: there is no complex Golay pair of length
  9 = 3 · 3 (EXHAUSTIVE over 4⁹ sequences, by the judge).
- *Spence's construction over an odd group* (VALIDATED). It does not need an even group: four ±1 sequences on an abelian
  G of order v with Σ PAF = −4 off 0 have row sums (±1, ±1, ±1, ±1) (Σ r² = 4, r ≡ v mod 2); with S the 4 × 4 Hadamard
  matrix of the row sums, P = Q = S, K = −S and M the group-developed Goethals–Seidel array, [[K, Q ⊗ 1ᵀ], [P ⊗ 1, M]]
  is a Hadamard matrix of order 4(v + 1) (H Hᵀ = nI checked for ℤ_v, v = 3 to 13, orders 16 to 56, and by the judge for
  v ≤ 11). `verify.py` accepts only the row sums (2, 0, 0, 0), so a hit needs the patch in `eightp/oddspence/`; sampler
  hits at |G| = 21, 27, 35 (orders 88, 112, 144) pass the patched builder. Consequence (EXACT): the row sums allow
  balanced blocks, so for odd |G| the multiplier group is not bounded by a row-sum condition as in the plain structures.
- *The plug-in census* (EXACT filters). G = K × F_r of order 2p − 1, blocks sI + tQ + eJ constant on the orbits of H_K ×
  squares(F_r*) (cyclic and product H_K, twists, 2 to 40 orbits), the defect identity (r + 1)e₊ + (r − 3)e₋ = 4(|K| − 1)
  and the row-sum identity Σρ² = 4 in orbit units: of 464, 3, 25, 169, 69, 22 shapes at p = 599, 719, 739, 839, 1031,
  1087 the live ones are 337, 1, 7, 49, 26, 0; at p = 1087 (2p − 1 = 2173 = 41 · 53) no allocation of orbit-unit defects
  exists. The exact τ filter (the same identities with the exact t-sums of one block, joined over four blocks; 0 false
  kills on the 35 odd-|G| shapes with a nonzero count, 91 of the 124 empty ones killed) was applied to the 213 live
  shapes with at most 20 orbits (212) or 21 orbits (one): p = 599: 87 alive (8 to 21 orbits), 75 dead, 3 undecided (DP
  memory); 719: dead; 739: one alive (ℤ₂₁₁ × F₇, |H_K| = 15, 15 orbits); 839: 11 alive (9 to 20 orbits); 1031: all 20
  dead. So p = 719, 1031 and 1087 have no plug-in shape with at most 20 orbits. 207 of the 208 live shapes with 21 or
  more orbits are not swept (UNKNOWN).
- *The exhaustive engine* (VALIDATED, EXHAUSTIVE on 14 shapes). The K-level solver of §86 against an independent G-level
  brute force on 166 Spence shapes of odd |G| (K up to ℤ₁₅, ℤ₃², ℤ₅ × ℤ₃; r = 3 to 9; twisted and untwisted H_K): 159
  MATCH (35 nonzero, the largest 741,520,000 ordered quadruples), 0 MISMATCH, 7 time-limited; the judge re-counted 75 at
  G level with its own code and re-ran eight with the engine. All allocations, exact lists and exact level 2 give no
  solution on 14 distinct shapes: p = 599: ℤ₇ × ℤ₁₉ with r = 9 (six shapes, 6 or 8 orbits), ℤ₉ × ℤ₁₉ with r = 7 (6
  orbits), ℤ₃² × ℤ₇ with r = 19 and ℤ₃² × ℤ₁₉ with r = 7 and a sign twist (8 orbits); p = 719: ℤ₄₇₉ × F₃, H_K = ⟨2⟩ (3
  orbits; the judge's G-level count agrees); p = 739: ℤ₂₁₁ × F₇ with ⟨4⟩ and ⟨5⟩ (3 and 7 orbits); p = 1031: ℤ₂₂₉ × F₉
  with ⟨3⟩ and ⟨4⟩ (5 and 7). Thirteen of the fourteen are also τ-dead. Not finished (UNKNOWN): three shapes at p = 599
  (8 and 9 orbits, 300 s limit, no solution so far), three at p = 839 (ℤ₁₃ × ℤ₄₃, r = 3: memory of the list generator),
  six at p = 1031 (ℤ₃ × ℤ₂₂₉, 687 points, over the 640 of the engine).
- *Plain cyclic Spence* (EXHAUSTIVE per shape, VALIDATED). Four ±1 sequences on ℤ_v, v = 1197, 1437, 1477, 1677, 2061,
  2173, constant on the orbits of every subgroup H ≤ ℤ_v*, all row sums +1, Σ PAF = −4 off 0, the count by a bucketed
  hash join: with at most 22 orbits there are 55, 5, 33, 45, 26, 22 shapes, of which 12, 2, 6, 2, 2, 4 have a sign
  vector with row sum +1 (the others are empty for that reason), and all have count 0; the same for the shapes with 23
  and 24 orbits at v = 1197 (20 shapes), 1677 (7) and 1477 (3). The judge re-counted every shape with a +1 vector at 22
  orbits or fewer, and those of 1197 and 1677 at 23 to 24, with its own joins and subgroup census. The first solver lost
  matches when it split a bucket of more than 4 · 10⁷ pairs (a hash residue modulo a number that is not a power of two;
  found by a test with a tiny bucket limit); the shapes it could have affected (1197 and 1677 at 24 orbits, 1477 at 21,
  the non-cyclic sample below) were re-solved with the corrected solver, EMPTY again. Non-cyclic ℤ₃² × ℤ₁₃₃: 40,000
  random 2-generator subgroups of GL(2, 3) × ℤ₁₃₃*, 583 distinct partitions, 283 with a +1 vector, all empty (a sample:
  HEURISTIC coverage). v = 1437 is complete (no subgroup of ℤ₁₄₃₇* has 10 to 479 orbits). UNKNOWN: v = 2061 at 23 orbits
  (three shapes with a +1 list of 291,720 vectors; the exact join was stopped after 3 h 13 min of CPU with no shape
  finished; days per shape by extrapolation, HEURISTIC) and v = 2173 at 23 orbits (154,000 vectors), and every shape
  with 25 or more orbits (1477: 28 orbits, 4.5 · 10⁶ vectors; 1677: 27 to 30 orbits).
- *Block-symmetric structures* (EXHAUSTIVE, VALIDATED). Blocks (A, B, μB, μ²B), (X, μX, μ²X, μ³X) and (X, μX, Y, μY)
  over ℤ_v with A (or X) invariant under K = ⟨H, μ⟩ and B under H of index 3, 4 or 2, all (H, μ) with at most 30 orbits
  of H: 655 pairs (238, 61, 356 for d = 3, 4, 2; 208 with a +1 list), count 0 for all six v (d = 3 exists only for v =
  1197, 1477, 1677, 2061: 3 does not divide |ℤ_v*| for 1437 and 2173). Validated against a brute force on 126 pairs (d =
  2, 4: 46 nonzero) and, by the judge, on 50 cases of d = 3 (7 nonzero; the engine counts pairs (B, distinct A vector)).
- *Bordered arrays* over |G| = n − t (t = 2, 4, 5, 6, 7): the exact filters leave 54 shapes, all with t = 2 (ℤ₄ or ℤ₂² ×
  ℤ_m × F_r; 21, 21, 12 at p = 599, 739, 1031). 52 have Gaussian estimates, all below 10⁻²⁵ except one (10⁷·⁵, p = 1031,
  ℤ₄ × ℤ₁₀₃ × F₅, 28 orbits, 76 t-sign bits); two 8-orbit shapes at p = 1031 (ℤ₄ × ℤ₁₀₃ and ℤ₂² × ℤ₁₀₃, r = 5) have none
  (singular covariance). None was searched (UNKNOWN).
- *Sampling and cost* (HEURISTIC). The two-level sampler with exact row patterns ran 360 s on each of the eight smallest
  alive shapes (seven completed; 372 of the 374 allocation tokens of ℤ₃² × ℤ₇ × F₁₉ exhausted, no instance); no round of
  any sampler gave a complete level-1 instance at p = 599, 739, 839, so the sampled coverage of the alive space is about
  0 (not evidence of absence). The modelled cost of the cheapest alive shapes is 10⁵·⁸ s (p = 599, ℤ₉ × ℤ₁₉ × F₇, |H_K|
  = 18, 13 orbits) and 10⁶·¹ s (p = 739, ℤ₂₁₁ × F₇, |H_K| = 15, 15 orbits) on one core; the Gaussian estimator is off by
  up to +2.8 decades (over) and, on six of 35 nonzero shapes, by 3.7 to 315 decades (under); it gives up to 10⁷ on
  exactly empty small shapes. Open (UNKNOWN): every alive shape (8 to 21 orbits at p = 599, 9 to 20 at 839, 15 at 739)
  and the unswept ones, plain cyclic Spence beyond the orbit counts above, non-cyclic groups beyond the sample, the
  bordered shapes.

## 93. Round 40: the lowest energies of the census shapes; skew 980 and 3124, the census values 245 and 781 (October 24)

Scratchpad `r40/`. The lesson of §91, applied to every open real skew census shape: sweep the slots of the lowest
energies with `exh38.py`, in increasing order of their number of level-1 classes.

**The sweep** (HEURISTIC coverage: lists of 60 to 99% of each raw block space, EX_COV). Of the 241 shapes scanned in
round 38, 35 were inventoried through E_min + 3, 23 through level 0 or 1 or partly at level 0, 76 had no feasible slot
within a 15 to 40 s cap, and 107 (63 of them q = 3 shapes with 17 or more orbits) were not inventoried. "Level" is the
energy above the lowest row-sum-feasible one, which for 781, 1057, 1477 and 2317 lies above the lowest energy allowed by
the identity. 2833 slots were swept (961 enumerable and 395 heavy ones at level 0, 1299 at level 1, 178 at level 2),
317,775 classes of level-1 instances went to the exact level 2 (133 of them, at 627 and 1539, only in a re-run: §94),
and 11 completed: the 2 of the control 385, 8 at 245 and 1 at 781. The control re-found 1540 in the 19th of the 34
enumerable lowest-energy slots of 385, after 474 classes. A positive control for non-cyclic K with a matrix H_K (the
slot of [skew/1100](../skew/1100)) gave the known matrices.

**The hits.**
- 245 = ℤ₃₅ × F₇ with ⟨11⟩ (group ℤ₃₅ × ℤ₇, not cyclic): the slot (6, 3, 4, 6) / (0, 3, 0, 0) at E = 19, one above
  the lowest feasible energy (N_w = 19, N₋ = 3, 8 · 19 − 4 · 3 = 140), holds 11,812 classes, 3 of which complete; other
  slots of the shape gave six more, and the judge's re-run a seventh (slot (6, 4, 4, 5) / (0, 0, 1, 2), row sums
  1, 3, 21, 23). Ten skew-Hadamard matrices of order 980 = 4 · 245 (a lower bound), pairwise inequivalent as
  Goethals–Seidel quadruples under Aut(G) and the block symmetries (row sums such as 1, −9, −27, 13 and 1, 3, 31, 3)
  ([skew/980](../skew/980)). 245 is undecided in the census; Table 5 had 245(4).
- 781 = ℤ₇₁ × F₁₁ with ⟨25⟩ (ℤ₇₈₁ cyclic): the slot (0, 6, 10, 11) / (0, 0, 5, 5) at the lowest feasible energy E = 27
  (N_w = 27, N₋ = 10, 12 · 27 − 40 = 284), a heavy slot (up to 6 · 10⁶ raw configurations per block specification),
  1841 classes, 1 completion: a skew-Hadamard matrix of order 3124 = 4 · 781, row sums 1, −51, 9, 21
  ([skew/3124](../skew/3124)). 781 is undecided in the census; Table 5 had 781(3).
Both slots replay alone in 1 to 3 minutes (`search/skewplug/data/repro245`, `repro781`). An independent judge rebuilt
the nine matrices of order 980 then known and both forms of 3124 with its own code (H Hᵀ = nI, H + Hᵀ = 2I, exact), ran
`verify.py` on a fresh copy, replayed both slots and redid the novelty checks.

**What the hits say** (HEURISTIC). The three census values settled by the enumerating engine (385, 245, 781) all came
from the lowest one or two energies of shapes that the Gaussian criterion of §90 rated unreachable (logE −1.5 to −2 at
their best slots), and at 781 from a slot that the round-38 production had skipped as too heavy. At the lowest energies
the block pools are small and a class completes with probability 10⁻⁴ to 10⁻², far above the estimates.

**The F1 family at low defect** (second lens of round 40, scratchpad `r40/f1low/`). The quaternary negacyclic plug-in of
§88 (skew 8v) has no defect pattern in `qp.c`; it was sliced by E, the number of positions of K carrying a defect of
positive cost (costs F₊ q + 1, Fᵢ q − 1, F₋ q − 3, M₊ (q + 1)/2, M₋ (q − 3)/2, Σ ℓ c = 4|K|). The sliced engine `qpe10`
reproduces the unsliced counts exactly (299: 13,056; ℤ₃ × F₇: 79,164,864; ℤ₇⟨2⟩ × F₁₁: 3840; ℤ₁₀₁⟨5⟩ × F₇: the level-1
weight 145,440), and both known hits, 299 and 589, lie at the lowest energy E = 3 (589 is found there in 0.5 s). Over
1000 (shape, E) slots, lowest first, gave no solution:

- EXACT within the engine class (q prime ≥ 7, K cyclic or ℤ_p², H = H_K × squares, ν ≤ 9): of the 128 census values with
  shapes, 94 have every shape empty at every feasible energy (among them 943, 1133, 1679, 2189 and 2303), 22 are Paley-I
  orders, 34 are partial.
- Partial: 779 = ℤ₄₁⟨10⟩ × F₁₉ (Table 5's 779(10)) is empty at E = 9, 12, 13, 14 (lists up to 4 · 10⁶, 6.5 · 10¹⁰
  probes) and E = 17 to 19 were not run; 759, 627, 623, 575, 2107, 1387 and 1539 are empty at their lowest slots, while
  their heavier slots (lists of 10⁶ to 10⁷, 10¹¹ to 10¹² probes) are UNKNOWN. 627 = ℤ₅₇⟨4⟩ × F₁₁ at E = 21 and 1075 at E
  = 3 were sampled at level 2 only (HEURISTIC).
- q = 3: the skew list is empty for all 24 shapes with ν = 3, for 1083 and 1587 (ν = 5) and for thirteen shapes with ν =
  7. For an odd number of orbit pairs there is a parity reason: R_A = w₀ + (1 + i)S with S odd, and |R_A|² ≤ 4|K|q fails
  for large ℓ; ν = 5 rests on the computation. The other q = 3 shapes (1047, 1527, 1671, 1983, 2271, 597 and nine with ν
  = 9) need a list builder that treats the zero-cost defects F₋, M₋ as free.

The judge rebuilt `qpe10` against the repository's `qp`, reproduced the validation counts, more than 25 exhaustive-empty
slots with identical list sizes and probe counts (all slots of 943, those of 779 up to E = 14, 1679, 1133, 2189, 1539 at
E = 17, the low slots of 759, 627, 623, 575, 1387) and the ν = 3 verdicts with the independent `qskew`. Not rederived:
the per-value tally, the ν = 5 verdicts beyond the earlier runs (no independent engine) and the level-2 samples.

**Corrections of the round-40 judge.** A re-run of the other slots of 245 gave a tenth matrix of order 980 (above). The
slot count of the sweep is 2833 (a draft said 2755). The Gaussian estimate against the measurement, over all 1535 swept
slots with an estimate: 1228 slots had no class; measured/estimated 0.92 for estimates of 10 to 100 classes, 0.36 for
100 to 1000, 0.34 for 1000 to 10⁴; 39 of the 950 slots estimated below 0.1 classes had a class (HEURISTIC).

## 94. Rounds 41–43: the remaining skew slots, symmetric propus and complex plug-ins at the lowest energies (October 24)

Scratchpad `r41/`, `r42/`, `r43/`. No new matrix. The judges of rounds 41 and 42 did not run (a usage limit stopped
them), so their verdicts below are those of the searches, with the validations each search ran on its own engine; the
code is not yet in the repository.

**Round 41: the rest of the skew census shapes.**

- *An exhaustive list engine* (`exh42.py`, `gl1e.c`; VALIDATED). It enumerates every configuration of a block
  specification once instead of sampling: the enumerated sets equal round 40's covering sets byte for byte on 32
  specifications of 8 shapes, the class counts equal round 40's on 16 slots of 245 and 385, and the hit slots replay
  (385 and 781: the same matrices; 245: the same 11,812 classes and 3 completions). Lists are built 3 to 50 times
  faster.
- *Heavy slots* (one share of the shapes). 113 slots swept, 96 EXHAUSTIVE for the engine's reductions and 17 with
  sampled lists; 2 level-1 classes (957 = ℤ₈₇ × F₁₁), 0 completions. Out of reach of the engine: 1235, 1679 and most of
  2323 (block specifications of 10⁸ to 10¹⁰ configurations), and the four best slots of 707 = ℤ₁₀₁ × F₇ (estimates up to
  6.8 · 10⁴ classes, blocks of 3.5 · 10⁹ to 1.2 · 10¹⁰ raw configurations: they need a streamed join).
- *The whole energy range* (the other share). An inventory over the full identity range takes seconds per shape. Of the
  76 shapes without a feasible slot in round 40, 57 have one, 13 have none (EXACT) and 6 were not reached. Feasible
  energies are sparse sums of orbit weights (917 = ℤ₁₃₁ × F₇ has only E = 66, 92 and 118; round 40 stopped at E_min +
  3). Swept: 1424 enumerable lowest-energy slots (17,834 classes) and 329 heavy or higher-energy slots (574,359 classes:
  917 at E = 92 and 118, 1057 at E = 121 and 151, 1267, 1477); 0 completions, and the level-2 row filter rejected every
  class of the higher energies. A heavy slot of 627, 737 or 803 costs 2 to 3 CPU-hours.
- *The estimator* (HEURISTIC, 2200 slots of round 40). The Gaussian slot estimate is right to a factor of about 3
  between 1 and 100 classes, too low by 5 to 60 times below 0.1 class and too high by 2 to 3 times above 100; all ten
  completions of round 40 lie in slots with logE ≥ −2, and the slot of [skew/3124](../skew/3124) ranks first of the 51
  lowest-energy slots of 781. Almost all slots left unswept have logE ≤ −6.

**Round 42: symmetric propus plug-ins** (the propus array of §17 over K × F_r with a multiplier group on K: A symmetric,
B twice, D free, weights (1, 2, 1) in the defect identity; §87 had run H_K = 1 only).

- *Engine* (VALIDATED). Raw counts equal round 33's K-level and G-level counts on 176 shapes; at ℤ₁₃ × F₇ with H = 1,
  slot E = 13, allocation (1, 3, 6), round 33's independent engine and the new one both give 843,648 level-1 instances
  and 86,528 solutions; the controls 308 (ℤ₁₁ × F₇) and 364 (ℤ₁₃ × F₇, ⟨3⟩) are found; the certificate path passes
  verify.py.
- *Results* (EXHAUSTIVE within the engine class for every finished shape; 0 completions everywhere). (a) The 654
  K-symmetric shapes of round 33's list (ν ≤ 22), every energy: 120 have row-sum-feasible slots, 14,888 slots, 4151
  level-1 classes; this includes 209 (ℤ₁₉ × F₁₁, ⟨7⟩), 259 (ℤ₃₇ × F₇, ⟨10⟩), 287, 341 and 247, the best-ranked shapes
  that round 33 never counted. (b) 630 more shapes with ν ≤ 13 and r ≥ 5, every energy: 368,404 slots, 1373 classes. (c)
  2143 shapes with ν = 7 to 22 and |K| ≤ 900, the lowest four energies: 148,234 slots, 163,416 classes. (d) 24 shapes
  with r = 3, every energy (for example 267 = ℤ₈₉ × F₃ with ⟨2⟩: 741,080 classes): 95,203 slots, 4,842,692 classes.
- *Why the shapes die* (EXACT, the 534 shapes of (a) without a feasible slot). 234 force the symmetric block to be a
  defect at every position, 225 have no allocation of orbit weights that meets the defect identity, 73 have no integral
  energy, 2 have allocations but no row-sum pattern.

**Round 42: complex plug-ins with 14 to 24 orbits per block** (the plain Szöllősi array of §85 for the 112 open orders
2N; round 30 had exhausted 4 to 13 orbits for cyclic K with |K| ≤ 120).

- *A spectral filter* (PROVED, VALIDATED). The level-1 equation Σ_X PAF_{u_X} = 4N δ gives Σ_X |û_X(ψ)|² = 4N for every
  character ψ of K, so every block of a solution has |û(ψ)|² ≤ 4N; applied to each list before the join it removes 28 to
  99% of the large lists. With orderly generation of the defect sets (one per orbit of the block symmetries) the cost
  moved from the join to the enumeration. Validated against round 30's counts (3072, 159,744, 27,648, 0), ten of its
  empty verdicts, the Γ-class sets of `exhaust_c` on nine shapes and the order-418 control.
- *Results.* 892 shapes inventoried: 564 have no row-sum-feasible slot in the examined energy window (EXACT), 191 have
  slots, 137 are UNKNOWN (identity energy beyond the orbit DP, or the inventory timed out). 32,380 slots of 160 shapes
  swept (32,371 EXHAUSTIVE), 7903 level-1 classes, 0 completions. Not run: about 900 slots of 10 to 100 s each and the
  heavier ones.

**Round 43: a streamed join for the heavy slots.** The best unswept skew slots of rounds 40–41 failed only because one
block's list (up to 1.2 · 10¹⁰ raw configurations) could not be stored. Round 43 (scratchpad `r43/`) keeps the three
small blocks in a hash table and streams the heavy block's configurations through an exact spectral filter (Σ_X
|û_X(ψ)|² = 4|G| at every character ψ of K, so each block has |û(ψ)|² ≤ 4|G| minus the least value of the other blocks),
at 2.3 · 10⁸ configurations per second, then joins and deduplicates the survivors as before. It needed a new level-2
solver, `skf4`, that fits under the memory limit.

- *Validation* (VALIDATED; the round's judge re-ran the parts marked). The three hit slots replay class for class with
  the same matrices (245: 11,812 classes, 3 completions; 385: 100 and 2; 781: 1841 and 1, sha256 of
  [skew/3124](../skew/3124)) (judge); the class sets equal the stored join's on 15 slots of 10 shapes; the judge's own
  enumerator gives the same filtered lists for 15 block specifications of up to 1.2 · 10¹⁰ raw configurations and the
  same 14 and 68 classes on two slots of 707; `skf4` agrees with `skf2` on 3009 instances and with `sl2` on the Z₅₇ and
  707 shapes.
- *Results* (EXHAUSTIVE within the engine class, 0 completions). 707 = ℤ₁₀₁ × F₇: the four slots of the lowest energy E
  = 51 (14, 68, 5731 and 6577 classes; the Gaussian estimate was 10 to 99 times too high) and two slots of E = 53 (1805
  and 615 classes); 1083 = ℤ₅₇ × F₁₉ and 1159 = ℤ₆₁ × F₁₉: 33 slots, 6 of them with 1 to 21 classes, the others with
  none. Unfinished: one slot of 707 at E = 56 (2 of 6 patterns, 2851 classes, 0 completions), 1235 [4, 3, 3, 3] (1 of
  133 patterns) and about 720 queued slots; one slot of 707 at E = 58 exceeds the memory of the cell table (UNKNOWN).
- *A silent failure, and its correction.* The level-2 solver `skf2` needs about 2.9 GB for one class of the Z₅₇ shapes;
  in a batch it ran out of memory under `ulimit -v` (exit 3, no output), and `exh38.py` counted the batch as run without
  completion. In round 40 this hit 627 = ℤ₅₇ × F₁₁ (112 of the 115 classes of the slot [18, 0, 0, 1] / [0, 0, 0, 0]) and
  1539 = ℤ₅₇ × GF(27) (18 classes in two slots); every other level-2 verdict of rounds 38 and 40 has a result for each
  class, and round 41 had level-1 classes only at 957. The 133 classes of the three slots were re-run with `skf4` (14 of
  them also with `sl2`): all decided, 0 completions, so the statements of §93 stand. `skf_io.py` and `sl2io.py` now mark
  every instance without a result line as L2FAIL.

**Reading** (HEURISTIC). In the skew array the lowest energies gave completions at rates far above the estimates (§93);
in the propus and plain complex arrays five million classes at the lowest energies gave none. The excess of the skew
family does not carry over to these arrays at these sizes; what is left there is either heavy (streamed joins) or has
estimates far below one.

## 95. Round 44: a conference matrix of order 118 (symmetric 940, n = 2823); the orders 8 · 1087 and 8 · 1619 (October 24–25)

Scratchpad `r44/`. Five lenses (algebra, literature, routes, closure, algorithm), each proposal with adversarial judges.
The disk filled during the round: the workflow journal was truncated, and the synthesis and some judge verdicts were
lost. What follows was rebuilt from the surviving journal lines and the judges' files, and every positive result was
re-run.

**Closure.** A catalogue of the objects known in three tiers (published by 2025, published in 2026 by others, this
repository) and of the recursive rules of the tables (Paley, Miyamoto, Yamada–Kiyasu, Cooper–Wallis, Seberry–Yamada's
8q, the Corollary 25/27/29/31 family for Williamson-type matrices, doubling), closed under the rules, against the
FRONTIER list, the six open orders 8p and the skew and symmetric census values.

- *Reached:* one target, 8696 = 8 · 1087. Seberry–Yamada's Lemma 8.12 with Corollaries 8.14 and 8.18 gives a Hadamard
  matrix of order 8q (q ≡ 3 mod 4 prime) from symmetric Williamson-type matrices of order (q − 1)/2. Those of order
  543 come from Corollary 27(i) (p = 541, (p + 1)/2 = 271 prime), which `symwt.py` built since §26. §92 had set the
  theorem aside because the 2026 survey excludes order 543; Table A.1 records it with key wi, a valid route. Built
  ([hadamard/8696](../hadamard/8696), `sy8q.py`, verify.py construction `seberry-yamada-8q`); two judges rebuilt it,
  one byte for byte, one with another primitive root, and checked H Hᵀ = 8696I exactly. **Known** (Table A.2 and
  Cati–Pasechnik's 1087(3)); the survey's list of open orders 8p should drop it. Five remain: 4792, 5752, 5912, 6712
  and 8248.
- *Beyond the 8p list:* 12952 = 8 · 1619 from Williamson-type matrices of order 809, which Corollary 8.8(ii) gives from
  [symmetric/404](../symmetric/404) (`sy8q.py c29`). Table 4's 1619(3) rests on Table A.1's key wl (Corollary 8.8(i)),
  so **First explicit construction** ([hadamard/12952](../hadamard/12952)). In the same way
  [symmetric/596](../symmetric/596) gives Williamson-type matrices of order 1193, which repairs the routes of 3579 and
  5965 (Corollary 8.8(i) list of NOTES.md §7: 373, 613, 757, 1109, 1301 and 1493 keep no route).
- *Near misses* (a single missing object unlocks the most targets): Williamson-type matrices of orders 599, 223, 419
  and 359 and a symmetric conference matrix of order 478 (four targets each), then T-matrices of order 223 and
  Williamson-type matrices of orders 719, 739, 839, 863, 947, 955, 191, 971, 167, 573 and 179 (three each). WT(223),
  for example, gives n = 669, 1115 and 2899 by Cooper–Wallis with T-matrices of orders 3, 5 and 13. The skew census
  part found no direct family hit and no improvement of Wallis's products.

**A conference matrix of order 118.** A symmetric conference matrix of order 118 (a conference graph SRG(117, 58, 28,
29)) gives, by Corollary 27, symmetric Williamson-type matrices of order 235 with a repeated block, hence the census
value v = 235 (Seberry–Balonin) and FRONTIER n = 2823 (Miyamoto blocks, §61). The literature lens proposed Mathon's
block structure (Canad. J. Math. 30 (1978): order pq² + 1 with blocks indexed by a field of order p and hyperplanes of
AG(2, q), for p = q + 2). Here 117 = 13 · 3², p = 13 ≠ q + 2, so Mathon's theorem does not apply; the judge made the
structure non-circulant and solved it as an exact cover ([search/conference118](conference118/README.md)).

- *Ansatz.* Points ℤ₁₃ × AG(2, 3) and ∞; diagonal blocks J − I − 2R (R the graph of two of the four directions), the
  block (i, j) equal to −χ(i − j)(J − 2M) with M a hyperplane of AG(2, 3)². On the characters of ℤ₃², C² = 117I splits
  into Paley's core of order 13 (trivial character) and, for the nontrivial ones, a cubic "flag graph" on the 52 pairs
  (point, direction) whose components must satisfy W² = 13I for W = D + 2Φ over ℤ[ω] (D = ±1 by direction class, Φ
  the signed cube-root phases of the edges).
- *Search.* Solvable template classes (`tmpl.c`): K₄ 48, prism 64, K₃,₃ 64, cube 256. Their placements on ℤ₁₃ with
  Paley signs (`inst.cpp`): 585, 17160, 624 and 196560. An exact cover of the 78 pairs by K₄ components alone
  (`xc.cpp`) succeeds after 13,589 nodes; 13 K₄ covering the pairs of ℤ₁₃ are the lines of PG(2, 3). `build.py`
  assembles C and checks C = Cᵀ, zero diagonal and C² = 117I exactly. `run.sh` reproduces the stored matrix
  (`data/conference118_5fe566972f21790b.txt`) byte for byte in about 10 s; no randomness.
- *The matrix.* Its signed automorphism group (pynauty, on the signed graph) is {±I}: the matrix has no symmetry
  beyond negation. The same ansatz over ℤ₁₇ gives conference matrices of order 154 (two
  solutions checked, C² = 153I); they give no new order in range (154 = (q + 3)/2 for q = 305 and (q − 1)/2 for
  q = 309, neither a prime power; 2 · 154 · 153 = 4 · 11781 lies beyond the census), and Ho and Wang reported graphs
  on 153 vertices as well (§61).
- *Consequences.* (a) `symwt.py 'cor27c(233, …)'`: symmetric Williamson-type matrices of order 235 with W₃ = W₄,
  checked (one judge's concern that they might not commute was settled by the exact check). (b) The Miyamoto blocks
  with q = 941 = 4 · 235 + 1 and T-matrices of order 3: [hadamard/11292](../hadamard/11292), n = 2823, **New**; one
  judge rebuilt it independently. (c) Seberry–Balonin's Corollary 1: [symmetric/940](../symmetric/940), v = 235,
  **Concurrent**: Brouwer's table cites J. Y. Z. Ho (email, October 2026) for a conference graph on 117 vertices, not
  public, and that graph gives 940 through the published corollary (as for 2284 and Zaitsev's matrix). The conference
  matrix itself is Concurrent with Ho's. (d) `cor29r(941, …)`: symmetric Williamson-type matrices of order 941, hence
  Cooper–Wallis orders 941 · t (not yet harvested).
- *Labels.* 2823 is New and not Concurrent: Ho's graph alone does not give it through published theorems (the
  Miyamoto blocks and the repeated-block arrangement are ours, §7, §26).

**Refuted or not pursued** (judge verdicts that survived).

- Singer–Paley difference families at q = 25: none of 1.1 · 10¹² triples passes at index 4 (EXACT within the family).
- A quaternary Paley-II doubling: trivial, no new order.
- Egan et al.'s Theorem 7 over q = 25 (n = 2503): estimated about 2% chance per complete run; not run.
- Spence's construction at 839 and 1187: dead by the slice equation (§3.1); the proposed flat filter is not new.
- Decoding Zaitsev's quaternary matrices into real routes: decodable, but no cheap sweep.
- Tied-stratum cores: about 1%, a repeat of §80. Decoupling of PGP(194): a repeat of §66.
- Non-circulant Williamson-type matrices of order 369 for 5912: expected count far below 1.
- A flat-pool k-tree for n = 1103: refuted by the one judge whose verdict survived.

**Next.** Harvest 941 · t and the orders beyond the tables that Williamson-type matrices of orders 235 and 941 give;
look for an R-structure in the solutions of the exact cover (it would give v = 119 by Yamada–Kiyasu, §61); run the same
ansatz for p = 29, 53 and 89 (orders 262, 478 and 802, the symmetric conference matrices that would give the census
values 523, 955 and 1603, and with 478 also FRONTIER 955).

## 96. Round 45: the flag equation for conference matrices of order 9p + 1; FLC's engines at order 65 (October 25)

Scratchpad `r45/`. Two lines: the ansatz of §95 for every prime p ≡ 1 (mod 4), aimed at the conference graphs on 9p
vertices that Brouwer's table (read 25 October 2026) lists as open: 261, 333, 369, 477, 549, 657 and 801 (p = 29, 37,
41, 53, 61, 73, 89; 117 and 153 exist, by J. Y. Z. Ho and Chen Wang). Their consequences: C(262) gives the census value
523 (Seberry–Balonin, q = 521) and, through Yamada–Kiyasu with q = 521, the complex order 1046; C(334) gives the
complex order 334 (C + iI), the smallest open one; C(478) gives FRONTIER n = 955 and 2865 (Corollary 27 with q = 953,
then Cooper–Wallis with T-matrices of order 3), the census value 955 and the complex orders 478 and 1910; C(550)
gives the complex order 550 and the census value 1099; C(658) the complex order 658; C(802) the census value 1603.
And FLC (github.com/vatskishan03/FLC), whose branches `w65-engine` and `good-engine` add classifiers with a five-fold
compression.

**The flag equation** ([search/flag9p](flag9p/README.md)). On the nontrivial characters of ℤ₃², C² = 9p I becomes
Ψ² = p I for Ψ = D + 2Φ on the 4p flags (point, direction), D = +1 on two directions and −1 on the other two, Φ with one
entry χ(j − i)ω^c in each 4 × 4 block (i, j), i ≠ j (EXACT; §95 is the case p = 13). Each component of the flag graph
satisfies it separately; its trace forces an even number of flags, half with D = +1. Bipartite components are X X* =
d I (d = (p − 1)/4); a complete bipartite K_{a,a} needs D constant on each side and a Butson matrix BH(a, 6) with the
Paley sign pattern. Quaternary non-Hermitian variants (for the complex Yamada–Kiyasu core of order 478) add nothing: the
D-terms of ΨΨ* force every flag edge to have its reverse with Φ_vu = ±conj Φ_uv, and with single-class edges the
supports must be symmetric, so d is even (p ≡ 1 mod 8) (EXACT).

- *Cyclic, bipartite* (`bipcyc.c`; EXHAUSTIVE for the ansatz, both sign conventions): p = 13 none (as the general cyclic
  search of round 44); p = 17: 55,296 structures (27,648 with the first row canonical); every one built and checked
  is a conference matrix of order 154 (like the cyclic ones of §95; concurrent with Chen Wang's graph on 153 vertices).
  The decoupled ones (one orientation per row, 1728 with the first row canonical) place a complementary pair of length
  2 on translates of {1, 2} by 2, 8 and 10. p = 29: none. p = 37 and p = 41 (decoupled rows):
  running (`r45/alam`).
- *Alamouti sub-family* [[P₁, P₂], [−λP₂*tᵏ, λP₁*tᵏ]] (`alam.py`; supports forced into 2k·{1, …, (p − 1)/2}): p = 13
  none, p = 17 twelve, p = 29, 37, 41 none.
- *Odd multiplier groups* (`mulcyc.c`; all edge types, bipartite or not, phase twist when 3 | k; EXHAUSTIVE): p = 5,
  k = 1 gives 576 (Mathon's order 46, the control); none for (p, k) = (13, 3), (29, 7), (37, 9), (53, 13), (61, 15),
  (109, 27) and (73, 9) with the first sign convention (the second is running). Without a twist, X(1) ∈ kℤ[ω] needs
  k² | d.
- *Non-cyclic, clique and biclique components.* Four K_{d+1} and six K_{d,d} account for exactly the p(p − 1)/2 pairs and
  4p flags for every p, and D-balance needs d odd (p ≡ 5 mod 8, the class where the cyclic ansatz fails). Templates
  (`kn.c`, `kg.c`): K₄ 2 normalised solutions, K₆ 12, K₈ (p = 29) 1056, all K₈ in one switching class with D-labels
  (Seidel spectrum ±√13, (±√5)³); on ℤ₂₉ there are 4466 placements (12 affine orbits, one D-labelling each). BH(7, 6):
  in the 804,033 normalised matrices enumerated before the time cap (not complete), the 5040 distinct sign patterns form
  one class, J − 2P; no pair of 7-sets of ℤ₂₉ has a Paley sign matrix in it, so the clique and biclique covers of ℤ₂₉
  need other component types (UNKNOWN). With K₈'s through one point and bicliques between their other points, D-balance
  fails for every p (EXACT).

**FLC's engines** (digest of the branches, a Sonnet agent; the trees in `flcwt/`). w5 classifies circulant symmetric
Williamson matrices and w5g good matrices of odd order n = q·m (q = 3, 5, 7, m ≤ 33): compression mod m, lists by the
power spectrum, a 2 × 2 hash join of the compressed quadruples (step 4), then the fibres of each instance joined on the
full autocorrelation (step 5). Both test suites pass on x86 (AVX-512), and the published instance fingerprints are
reproduced. They reach no open Hadamard or skew order (4n ≤ 372), but they settle open cells of the classifications:

- *Good matrices of order 55:* exactly one inequivalent class (w5g, q = 5, 10,773 instances, 297 CPU-s; the same class
  without grouping). Row sums 1, 7, −13, −1; checked independently (own decoding and Goethals–Seidel array: A + Aᵀ =
  2I, B, C, D symmetric, Σ X Xᵀ = 220 I, a skew-Hadamard matrix of order 220). Good matrices of order 55 were not in
  the literature we read (Đoković–Kotsireas, arXiv:1802.00556: orders ≤ 39, 43, 45, 127; Bright et al.,
  arXiv:1811.05094: multiples of 3 to 69), but the vxvware Goethals–Seidel challenge board (GitHub, 2026) lists one
  class at 53, 55 and 61 from participants; the census (exactly one class) is new.
- *Williamson matrices of order 65* (open: Kharaghani et al., arXiv:2605.08661, and the board) and *good matrices of
  order 65* (the board has an unverified emptiness claim): running (`r45/flc`, PRODUCTION.md; 70,699 and 67,437
  instances, fingerprints c7000a01e7b84192 and d2617b1ba87e7899 as published; chunks of 1500 instances, 0 hits so far).

## 97. Round 46: the AG(2, 3) lift for Hadamard matrices, which reduces to a known plug-in; the order 4 · 1719 (October 25)

Scratchpad `r46/`. The idea was to turn the ansatz of §95 and §96 from conference matrices to Hadamard matrices: for an
Hadamard matrix K of order N, put H_ij = K_ij(J − 2M_ij) with M_ij a hyperplane {⟨a, x⟩ + ⟨b, y⟩ = c} of AG(2, 3)². If
the construction worked in general, it would give H(36w) from H(4w). That would settle FRONTIER n = 1719 = 9 · 191 from
Đoković's H(764), and the 29 census values v = 9w whose 4w is a known symmetric order: 369, 513, 621, 657, 693, 765,
783, 891, 909, 927, 1017, 1089, 1125, 1143, 1269, 1341, 1359, 1539, 1557, 1593, 1647, 1719, 1773, 1845, 1953, 1989,
2133, 2313 and 2349. It does not work in general. Its usable case is a known plug-in, and that plug-in needs an
input we do not have ([search/x9split](x9split/README.md)).

**The lift** (EXACT).
- *Trivial character.* The block J − 2M has row sum 3, so the trivial character gives 3K, and (3K)(3K)ᵀ = 9N I holds
  for every K.
- *Nontrivial characters.* Each block maps one character class at j to one at i with modulus 6, as in §96. Every point
  has N blocks over 4 directions.
- *The condition.* With real phases and one orientation, H Hᵀ = 9N I is equivalent to the following. The matrix W of
  order 4N, with W[(i, f_ij), (j, g_ij)] = K_ij and one entry in each 4 × 4 block, satisfies W Wᵀ = (N/4) I. With
  cube-root phases, W becomes a unit weighing matrix over the sixth roots.
- *The ℤ₃²-invariant case.* Take M_ij = {⟨a, x − y⟩ = 0}, so f = g. Then W is block diagonal and the condition reads:
  K = K₀ + K₁ + K₂ + K₃ with disjoint weighing matrices W(N, N/4), and H = Σ K_a ⊗ Y_a.
- *The blocks Y_a = J − 2L_a.* They satisfy Y_a Y_bᵀ = J for a ≠ b and Σ Y_a Y_aᵀ = 36I (checked). So the cross
  terms collapse to (K Kᵀ − Σ K_a K_aᵀ) ⊗ J = 0.
- *Status: Known.* The Y_a are a regular 4-set of order 9 (Seberry–Whiteman, Graphs Combin. 4 (1988), from Mathon's
  construction; Seberry–Zhang, Graphs Combin. 9 (1993)). The plug-in of a regular s-set into an Hadamard matrix that is
  a sum of s disjoint weighing matrices is theirs. A Baumert–Hall array, for example from T-matrices of order w,
  supplies such a split. In that case the order 36w is also Cooper–Wallis with Williamson matrices of order 9.
- *Symmetry.* H is symmetric if every K_a is (c = 0). A skew H is impossible: every diagonal block is ±Y_a, whose
  diagonal is −1.

**Other orders for the plug-in** (PROVED). Use a regular s-set of odd order m on a split of an Hadamard matrix of
order 4w (w odd) into s disjoint W(4w, 4w/s). Then s | 4w, and W = Σ Y_a satisfies W Wᵀ = smI + s(s − 1)cJ with c odd.
- *mod 4:* if s is even, W ≡ 0 (mod 2), so s(s − 1)c ≡ 0 (mod 4) and 4 | s. In particular there is no 2-split plug-in
  of odd order (no order 12w from a 2-split).
- *Determinant:* det(W)² = (sm)^m (1 + (s − 1)c), so sm(1 + (s − 1)c) must be a square.
- *s = 4:* this leaves only the odd squares among m ≤ 51 (9, 25, 49). So there is no ×3, ×5 or ×7 analogue; for
  example 1713 = 3 · 571, 2855 = 5 · 571 and 1981 = 7 · 283 are not reachable.
- *Search over s ≤ 15, m ≤ 51:* checked against the 122 FRONTIER orders, every (s, m) that passes both tests reaches
  only n = 1719 (s = 4, m = 9).
- Brute force finds no regular s-set of order 3 for s = 2, 3, 4 and 6.

**Splits found and not found.**
- *Goethals–Seidel quadruples, circulant splits* (`csplit.py`, CP-SAT, one label per position of each sequence; each
  layer needs Σ PAF = wδ).
  - Found: Williamson order 3, and three random quadruples of order 11, in under 21 s each. The resulting H(108) and
    H(396) were checked exactly.
  - With SA (`sa_split.c`: 4 restarts of 3 · 10⁷ swaps, 120 s) on quadruples from `gsgen.c`: all three of order 13
    split, and none of the three at each of orders 17, 19, 23 and 29.
  - Hall's Williamson quadruples of orders 19 and 23: CP-SAT UNKNOWN after 600 s, and SA stops at energies 5 and 9.
  - So neither method finds splits beyond order about 15, although the count below predicts many at order 19.
- *Propus quadruples, symmetric circulant splits:* none for w = 5, 7, 9, 11 (three quadruples each; CP-SAT proves
  infeasibility in at most 13 s).
- *Explicit symmetric matrices of order 12:* no split of Paley II (q = 5), and no symmetric 4-split of any symmetric
  H(12) (labels and signs free). Both UNKNOWN after 1200 s and 900 s. Symmetric W(4, 3), W(8, 3), W(12, 3), W(8, 5)
  and W(20, 5) exist (CP-SAT), so there is no obstruction at the level of the parts.

**n = 1719.** It needs a split of an Hadamard matrix of order 764. The only ones known are Đoković's two
Goethals–Seidel quadruples (arXiv:math/0703312): ℤ₁₉₁, H = ⟨39⟩ of order 5, supports of sizes 85, 90, 90 and 100.
Both were reconstructed (`dk764.py`) and checked: Σ PAF = 764δ.
- *H-invariant splits.* Each layer holds exactly one of the four zero positions (weight 191 ≡ 1 mod 5) and 38 cosets.
  CP-SAT ran 600 s with no answer. The best SA energy was 153 over the 76 shift classes.
- *Free circulant splits* (764 labels, 380 equations): the best SA energy was 985.
- *Count* (HEURISTIC). A layer's PAF at a shift has about 48 terms. Parity costs one bit per equation (the support
  autocorrelation over GF(2)), and the size of a nonzero even sum costs about 3.1 bits.
  - H-invariant: 2^292 balanced labellings against 76 · 4.1 bits gives about 2^−20 splits per quadruple.
  - Free: about 2^−40 per quadruple.
  - The same count gives about 2^60 splits at w = 19, where none were found.
  - A joint search (the quadruple and its split, H-invariant) has about 2^49 solutions at a density of 2^−419.
- *Not pursued.*

**Other results since §96.**
- The exact maximum packing of K₈ placements on ℤ₂₉ (`k8clique.py`, max clique on 4466 vertices and 2,492,028
  compatibility edges) is 4. That is what the four K₈ of a clique-and-biclique cover need, but the bicliques have no
  placement (§96).
- The paused searches were not completed:
  - `bipcyc3 37` (full bipartite, p = 37) and `mulcyc 73 9 1` are stopped; `bipcyc3d 41` (decoupled, p = 41) is
    running again.
  - The FLC Williamson-81 probe hit its time cap while stopped and was not resumed.
- FLC at order 65: 16 of 47 Williamson chunks and 12 of 45 good chunks done, 0 hits (running).

## 98. Round 46b: the sign classes of the bicliques and Lindsey's bound (October 25)

Scratchpad `r46/cover`, programs in [search/flag9p](flag9p/README.md). §96 left the clique-and-biclique covers of
ℤ₂₉ with one unchecked assumption, that the BH(7, 6) enumeration (stopped at a time cap) had found every sign class.
This round removes it, does the same for p = 37 (C(334), whose C + iI is the smallest open complex order), and adds a
spectral filter.

**Sign classes** (EXACT).
- *Why the inner products are ±3.* Seven sixth roots that sum to 0 form one triangle and two antipodal pairs. So every
  row and column inner product of a BH(7, 6) sign pattern is ±3.
- *BH(7, 6)* (`bh7signs.py`). The patterns with all inner products ±3 are 7 normalised matrices, all with Gram
  4I + 3J. They form one class, J − 2P, and it admits cube-root phases.
- *BH(9, 6)* (`bh9signs.py`). Nine sixth roots summing to 0 are three triangles, or one triangle and three pairs, so the
  inner products are ±3 or ±9. The 1009 normalised patterns fall into exactly three classes, each admitting phases:
  - rank one (cut norm 81);
  - one class with row inner products 3 and 9 (cut norm 63);
  - H′ ⊗ J₃ with H′ = [[1, 1, 1], [1, 1, −1], [1, −1, 1]] (cut norm 45).

**Lindsey's bound** (PROVED).
- *The bound.* For a conference core S of order p (eigenvalues ±√p and 0) and ±1 vectors x, y on disjoint d-sets,
  |xᵀSy| ≤ d√p. So a K_{d,d} component whose sign class is M needs ‖M‖_{∞→1} ≤ d√p.
- *J − 2P* (cut norm d(d − 2)) survives at p = 29 (35 ≤ 37.7) and fails for every d ≥ 9.
- *p = 37:* only H′ ⊗ J₃ survives (45 ≤ 54.7).
- *Rank one* needs d² ≤ p. Direct search (`rank1bic.py`) finds 144 such bicliques in the Paley two-graph of ℤ₁₃ (0 on
  one side) and none at p = 29 or 37.

**Paley cores** (EXACT within the clique-and-biclique family).
- *Possible counts.* With components K_{d+1} and K_{d,d} only, the flag count (d + 1)a + 2db = 4p and the pair count
  reduce to the same equation:
  - p = 29: (a, b) = (4, 6) or (11, 2);
  - p = 37: (4, 6) or (13, 1).
  Every option contains a biclique.
- *Placements* (`k99place.py`). Paley(37) has no K₉,₉ in the class H′ ⊗ J₃, and Paley(29) has none in J − 2P (§96).
- *Conclusion.* Neither Paley core has a cover by cliques and bicliques.
- *The combinatorics is not the obstruction* (`kkcover.py`, CP-SAT). D-balanced decompositions of K_p into four
  K_{d+1} and six K_{d,d} exist for p = 13, 29 and 37 (1, 27 and 58 s); none is of the one-point type of §96. The
  obstruction is in the signs.

**Other cores.**
- `srg_k77.py` asks whether any SRG(29, 14, 6, 7) contains a K₇,₇ in the class J − 2P. CP-SAT gave UNKNOWN after
  1800 s; the p = 13 control is found at once.
- The 41 conference graphs on 29 vertices and the 6760 on 37 could not be fetched (Spence's and Brouwer's pages
  returned 503 on 25 October), so the non-Paley cores are untested.

**Multipliers.** (73, 9) with the second sign convention has no solution either (0 of 254,803,968), so the odd
multiplier groups of §96 are empty for every (p, k) listed there. `bipcyc3d 41` (decoupled bipartite, p = 41) is still
running.

## 99. Round 46c: the interval family behind p = 17, and symmetric quaternary matrices for n = 2571 (October 25)

**The p = 17 solutions are one family.** The 1728 decoupled bipartite solutions of §96 (first row canonical) have
exactly four support partitions:

{1, 2, 3, 4} | {5, 6, 7, 8},  {1, 2, 4, 7} | {3, 5, 6, 8},  {2, 3, 5, 7} | {1, 4, 6, 8},  {3, 4, 5, 6} | {1, 2, 7, 8}.

They are the images under the multipliers 1, 3, 5 and 9 of one interval structure, with 432 solutions each. Row A1 has
atoms at the offsets 1, …, d and row A2 at −(d + 1), …, −2d. Each row is then a sixth-root complementary pair whose two
supports split an interval, with the Legendre signs prescribed.

- *Reduction* (EXACT). Every condition depends only on index differences (or sums) and on products χ(i)χ(j). So the
  multiplier and the sign convention drop out, and one run per prime decides the family.
- *Search* (`flag9p/intcyc.c`). The atoms are assigned alternately from both ends, so the long shifts close first.
  - p = 17: 27 gauge-fixed solutions, rebuilt by `build_bip.py` into a conference matrix of order 154.
  - p = 41 and 73: no row satisfies the row condition (0 s and 179 s).
  - p = 89 and 97: not decided within the 900 s cap.

So p = 17 looks like a small case: the Legendre pattern on [1, d] has no sixth-root complementary split at p = 41 or
73.

**Symmetric quaternary Hadamard matrices and n = 2571.** The literature agent's leading near miss, checked here
against [hadamard/MIYAMOTO_BLOCKS.md](../hadamard/MIYAMOTO_BLOCKS.md#symmetric-williamson-type-matrices):
- *The route.* Lemma 26 (Corollary 29) turns a symmetric quaternary Hadamard matrix of order (q − 1)/2 into symmetric
  Williamson-type matrices of order q. For q = 857 that order is 428, and Cooper–Wallis with T-matrices of order 3
  then gives n = 2571.
- *Known sources do not apply.*
  - Turyn's C ± I needs a symmetric conference matrix of order 428, impossible since 428 ≡ 0 (mod 4).
  - A real symmetric H(428) is the open census value v = 107.
  - A Kronecker product would need order 214, and Turyn's C + iI does not exist there: 213 is not a sum of two
    squares.
  - The repeated-block source would need symmetric Williamson-type matrices of order 107. A quadruple (Y, Z, X, X) of
    odd order m gives a symmetric quaternary matrix of order 2m through M₃ = [[Y, X], [X, −Z]] and
    M₄ = [[Z, X], [X, −Y]], but that order is again 2 (mod 4) (see below).
- *A searched structure* (`flag9p/symq_mitm.py`, exhaustive): L = [[X, WR], [WR, −X̄]], with X a symmetric quaternary
  circulant of order m and W real. L is symmetric, and L L̄ = 2mI reduces to PAF_X + PAF_W = 2mδ.
  - Solutions exist for m = 3, 5, 7, 9, 13, 15 and 19 (24 to 5016 of them, checked exactly). There are none for m = 11
    and 17, exactly the m where 2m − 1 is not a sum of two squares.
  - Every solution has an imaginary part of weight 1 or m − 1. So the structure is Turyn's C + iE, and Belevitch's
    condition explains the pattern.
  - Circulant repeated-block quadruples are a special case. If (Y, Z, X, X) is such a quadruple, then
    A = (Y + Z)/2, B = (Y − Z)/2 and W = X satisfy PAF_A + PAF_B + PAF_W = 2mδ. A direct search finds none for
    m = 11 and 17.
  - So every structured source of symmetric quaternary Hadamard matrices of order n ≡ 2 (mod 4) that we know needs
    n − 1 to be a sum of two squares. Whether this holds in general is open here (CONJECTURE; an unrestricted CP-SAT
    run at order 22 was stopped undecided).
  - Calibrated on m = 19 (predicted 2^12 solutions, found 5016), the same count gives about 2^−74 for the order-428
    version of this structure, so it was not run.
- *What remains.* The route to 2571 needs a symmetric quaternary Hadamard matrix of order 428 directly; a real one is
  the census value v = 107.

Agent reports of this round (literature, combinations): no new verified route. The near misses match §95. The
leading one is the order-107 object above, which settles n = 2571 and the census value 107. A further reference to
check is Seberry–Finlayson (JCMCC 43, 2002, skew orders 4q < 4000), which was not reachable.
