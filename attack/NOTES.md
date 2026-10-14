# Attack notes: structure, obstructions and new searches

Working notes for the push on the unresolved orders in [FRONTIER.md](../FRONTIER.md).
Everything here is either proved below, or measured with the code in this directory.

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

### 3.1 Spence's construction and the slice character (FRONTIER correction)

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
  E ≈ 10^−1 to 10^−5. These are lottery tickets that can be searched exhaustively
  (`find4.c`, `gs_exhaustive.py`, validated by constructing H(148) and H(172)).
  The best is n = 933 (order 3732), |H| = 62, 22 orbits, E ≈ 10^−0.9.
- Legendre pairs (ℓ = 2n − 1) have no coarse character except the trivial one,
  and the border absorbs it. They are 2-list problems, linear in the pool. The
  repository never searched them. `lp_search.py` recovers the classical pairs and
  many others for prime ℓ. For composite ℓ, the only case open orders produce,
  small-structure pairs are rare.

## 5. Results of the lottery searches

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
neither small Goethals–Seidel nor small Spence structures leave a useful lottery
ticket. The only exception is n = 933 with |H| = 31, which has about two expected
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
an interval. Validation: for each of the three ways to pair up main's T(113) blocks, the exact join on
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
above. This matches where main's T-matrices stop (t ≤ 197, at most 36 orbits).

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

## 8. Second pass: twisted Legendre pairs, local search, base sequences

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
  magnitude beyond this machine.
- **The d = 3 structure for 7572, run here** (`d3/prep1893.py`, `d3/walk.c`). This
  is the METHOD.md structure over ℤ₁₈₉₃: H = ⟨8⟩ of order 30, K = ⟨8, 2⟩ of
  order 90, 30 and 86 orbits, 15 classes, row sums (87; 1, 1, 1).
  - The complete first-block pool has 24,608,948 sequences, all stored in an exact
    hash table.
  - The model gives log₁₀ E = 2.15, and log₁₀(success per exactly centred lookup)
    = −14.7, matching METHOD.md.
  - The tabu walk over the second block checks every row-sum-preserving swap
    neighbour: 6.5·10⁶ lookups/s on 4 threads. That is at least 900 days per
    expected solution on this machine.

## 9. Structural ideas checked in the third pass

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

- **Structure.** Goethals–Seidel over ℤᵥ, all four blocks invariant under H ≤ ℤᵥ* with −1 ∉ H, and the
  first block of skew type (A(0) = 1, A(−x) = −A(x)). Negation pairs the H-orbits, so there are
  2^{(O−1)/2} skew first blocks, and the row sums satisfy 1 + b² + c² + d² = 4v.
- **Triage** (`skew_triage.py`, `skew_rank.py`) uses the Gaussian model of §2 with the skew A pool.
  **For composite v it is not to be trusted at all**: it is off by 10²–10⁵ in both directions, which is
  the same lesson as main's order 3964. Composite v are searched exactly (`chain.py`, below); the model only
  orders work for prime v.
  - Under: 209 with ⟨26⟩ was predicted at 3 solutions and has 288.
  - Over: 201 with ⟨22⟩ (|H| = 11, rows 1, 1, 19, 21) was predicted at 10³ and has none
    (exhaustive `find4b`).
  - Under: 377 with ⟨16⟩ was predicted at 0.1 and gave 64 in the first τ; 469 with ⟨9⟩ was predicted
    at 0.2 and has 192; 889 with ⟨2, 25⟩ came from the exact sweep, not from the model.
- **Matching.** `match4.c` is `search/match.c` with four distinct lists: left pairs (A, largest pool),
  right pairs from the other two. With lists of very different sizes, the yield per τ grows with the
  smaller pair list, so the clamp S is set as loose as memory allows (pair lists ≤ 2.5·10⁷).
  Pools beyond 4·10⁶ are replaced by uniform samples (`--maxpool`), at a cost of L/N per list.
- **New orders**: 836, 1036, 1100, 1148, 1236, 1364, 1412, 1508, 1612, 1852, 1876, 2404, 3556,
  3676 and 1372, i.e. v = 209, 259, 275, 287, 309, 341, 343, 353, 377, 403, 463, 469, 601, 889, 919. 889 came from the exact sweep, on
  H = ⟨2, 25⟩ of order 63 with 21 orbits. Each took seconds to minutes on 4 cores.
  - 259 used sampled lists: 3.9·10⁶ of 3–4·10⁷ per pool.
  - 353 (⟨22⟩) also used sampled lists: 1.2·10⁷ of 5.7·10⁸ per pool, 2 solutions in 33 values of τ.
  - 403 (⟨100⟩) used sampled lists of 1.05·10⁷ of 4–6·10⁷: 4 solutions in 154 values of τ on
    rows (1, 5, 25, 31). A 119 run with low-energy lists (8·10⁶ per role) found nothing.
  - 469 (⟨9⟩, predicted 0.2) gave 192 solutions in the exact matcher. The k-tree had missed them
    because 36 values of τ with the top-6 values per class do not cover the solutions. For pools
    of ≤ 10⁵, run the exact matcher.
- **Spin structure** (`spin_skew.py`, `spin_ls.c`), the type of Đoković's skew 4·397.
  - A is invariant under K = ⟨H, μ⟩ with μ of order 3, and C = B∘μ, D = B∘μ².
  - The condition P_A + P_B + P_B∘μ + P_B∘μ² = 0 lives on the ⟨K, −1⟩-classes, a third as many as
    the ⟨H, −1⟩-classes. It is a 2-list problem: B ranges over 10¹³–10¹⁵ choices, and a tabu
    search on (A, B) in C does about 5·10⁴ moves/s.
  - The row sums are forced: 1 + 3b² = 4v. For open v ≤ 999 that leaves v = 469, 721, 817 and 919.
  - 919 (|H| = 17, μ = 52, 9 classes, predicted 2·10⁵ solutions) was solved in 37 s.
  - 721 (|H| = 17, μ = 46, predicted 200) is still open after 6·10⁶ moves; the best residual is
    two classes off by 4.
  - 817: an earlier version required μ³ = 1 and found no structure with E > 10⁻³. That was too
    narrow, because μ³ ∈ H suffices: the blocks B, B∘μ, B∘μ² close up since B is H-invariant. With
    μ³ ∈ H, H = ⟨11⟩ (order 21) and μ ∈ {6, 9, 17, 24} give spin structures with 8 classes and model
    E ≈ 10³·². Tabu searches on them are running.
  - Since a skew block has row sum 1, spin needs 4v = 1 + 3b², i.e. v = 3m² + 3m + 1, a centred
    hexagonal number. Đoković's 397 and the DGK order 631 are of this form, as are the four values
    above.
  - `spin_skew.py exhaust` joins all A against all B exactly when the pools are small. For 817
    (|H| = 63, μ = 49 or 552) and 721 (|H| = 51) it finds no solutions; most of these
    structures have no B with row sum b at all.
- **General local search.** The same C engine handles plain skew GS: B, C, D are stacked into one
  vector with a block-diagonal form, and swaps stay inside a block (`spin_skew.py gsearch`). It
  solves 209/⟨26⟩ in under a second.
- **Compression obstruction** (`skew_compress.py`). This explains the Gaussian model's misses.
  - For d | v, the compressions to ℤ_d are complementary: Σ_X PAF_x(k) = 4v·[k = 0] for k ∈ ℤ_d. This
    is an identity of integers, equivalent to Σ_X |X̂(χ)|² = 4v at every character. The compressions
    are also constant on the H-orbits of ℤ_d, and the skew block has x₀ = 1 and x₋ⱼ = −xⱼ.
  - Enumerating each role's achievable compressed vectors and meeting the integer equations is exact
    and cheap for small d. These are the coarse characters of §1 again.
  - **Correction.** The first version compared the norms |X̂(χ)|² after rounding them to integers.
    For d ≥ 5 they are not rational integers (they lie in ℚ(ζ_d)⁺), so that test was not exact; only its
    d = 3 verdicts were sound. Every verdict below was re-checked with the integer test.
  - It proves there are no solutions for:
    - 201/⟨22⟩, rows (1, 1, 19, 21): d = 3, predicted 10³;
    - 301/⟨4⟩, rows (1, 1, 19, 29): d = 7, predicted 4;
    - 265/⟨16⟩, rows (1, 1, 23, 23): d = 5.
  - The solved structures pass. 209/⟨26⟩ has 11 and 21 feasible combinations at d = 11 and 19;
    469/⟨9⟩ has 1 at d = 7.
  - Any triage of composite v should apply this test before the Gaussian estimate.
- **Which tool when.** The k-tree enumerates every left pair (A, D) and every right pair (B, C) once
  over all values of τ. Its work per solution is therefore (N_A·N_D + N_B·N_C)/E′, with E′ the yield
  of the (sampled) lists. At about 10⁷ pairs/s on this machine:
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
- **Second-tier sweep** (pair work 1.5·10⁹–5·10¹⁰, ≤ 37 orbits). The first list of verdicts came
  from the rounded-norm test and is **withdrawn**. The lists are being recomputed with the integer test,
  and the structures are covered exactly by the compression-chain sweep below.
- **Compression-chain engine** (`chain.py`, `chain_sweep.py`), after main's 3964 search. For v = p·q
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
skew-Hadamard matrix of order (q + 3)/2 exists, then a Hadamard matrix of order 4(q + 2) exists. Main's
`verify.py` has the `yamada` builder, and `attack/make_yamada.py` nests one of our skew certificates in it,
after k doublings. For a skew order s (≡ 0 mod 4), q = 2s − 3 is automatically 5 (mod 8), so only its
primality matters.

| skew order s | q = 2s − 3 | Hadamard order 4(q + 2) | status of that order (agent audit, Oct 2026) |
|---|---|---|---|
| 836 | 1669 | 6684 (n = 1671) | known; n ≤ 2999 not in the open list; built as a pipeline test |
| 1036, 1100, 1148, 1372 | 2069, 2197 = 13³, 2293, 2741 | 8284, 8796, 9180, 10972 | known (n ≤ 2999, not open) |
| 1564 | 3125 = 5⁵ | 12508 (n = 3127 = 53 · 59) | known: T-sequences of length 59 × Williamson-type of order 53 (Handbook "."); built as a test of the prime-power case |
| 1612 | 3221 | 12892 | known: T-matrices of order 11 × Williamson-type of order 293 |
| 1852 (Đoković 1994) | 3701 | 14812 | known: T(7) × Williamson-type 529; also Kiyasu from Đoković's 1852 |
| 2200 = 2·1100 | 4397 | 17596 | known: T-sequences of length 53 × Williamson-type 83 |
| 3016 = 2·1508 | 6029 | 24124 | implied since 2026-08-30 by T-matrices of order 163 (JD-Jones-ASES/Hadamard-T) × Williamson 37, but stated nowhere we found; main's lead |
| 3556 | 7109 | 28444 | known: Paley II (q = 14221) |
| 3676 (Schumer 2026) | 7349 | **29404** (n = 7351) | implied since 2026-09-17 by Schumer's skew 3676, stated nowhere we found; Handbook Table 1.53 marks 7351 ".", apparently in error (every twin-prime n ≡ 3 mod 4 without a Paley route between 4000 and 9999 is marked "."); built from our skew 3676: [29404](../29404) |
| 4808 = 2·2404 | 9613 | 38460 | known: Paley I (q = 38459) |
| **6448 = 4·1612** | **12893** | **51580 (n = 12895 = 5·2579)** | **probably open**: Paley I/II fail (51579, 25789 composite); no Williamson-type family reaches 2579 or 12895; Seberry–Yamada 1992 Table A.2 gives only 8·2579; n > 10000 is beyond Handbook Table 1.53; no 2026 deposit mentions it |
| 7112, 7408, 14816 | 14221, 14813, 29629 | 56892, 59260, 118524 | known (Paley I, Paley II, Williamson-type routes) |

So the new skew order 1612 (4·403), doubled twice, gives the Hadamard order 51580, which is probably new, and
Schumer's (and our) skew order 3676 gives 29404, whose Handbook entry we could not trace.

**Prime powers.** `build_yamada` and `make_yamada.py` also accept q = p^k: GF(q) is given by a primitive polynomial
over GF(p) in the certificate (fields `p`, `poly`), and r and ξ are encoded by base-p digits. Checked on q = 13, 29
and 5³ (orders 60, 124, 508) and on q = 5⁵ with skew 1564 (order 12508, H Hᵀ = 12508 I, 42 s).

`verify/verify.py 51580` was run separately (2026-10-12, 40 minutes, about 3.5 GB): `51580 ok`, with the recorded digest.

## 14. The census of undecided skew orders: leads and calibration

The Đoković–Seberry census lists 737 odd v < 2500 with no skew-Hadamard matrix of order 4v (May 2025;
`attack/data/skew_census_undecided.txt`). Leaving out the 13 orders of this repository and Schumer's 919,
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
- a CP-SAT model with the products of signs as XNOR literals (OR-Tools 9.12): UNKNOWN after 600 s on 2
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
  matrix of order 4 · 6103 = 24412 exists. It is **probably new**, see [24412](../24412):
  - Handbook Table 1.53 lists 6103 with t = 3, read from an unofficial copy and checked against known
    entries. Paley I gives 8 · 6103.
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
  - `split_join` turned out not to be needed: the exhaustive joiner `attack/join4x.c` makes the single-level
    final joins cheap (next item), and the extra level is now off by default (`--split`).
- **The exhaustive joiner `join4x.c`.** It partitions both pair sides exactly on the three classes with the
  largest spread and probes small open-addressing tables keyed by a 64-bit hash of all classes. Partition keys
  are processed in lexicographic chunks that bound the build side; since lexicographic order is translation
  invariant, the probe cells that complete a chunk form an interval, found by binary search. The build phase
  inserts in parallel with compare-and-swap. It agrees with `find4b` on random instances and reproduces the
  288 solutions of 209/⟨26⟩. On the heaviest 623 branch (fibres 1,120, 3.4·10⁵, 2.7·10⁵, 2.2·10⁵) it
  finishes in 103 s; `find4b` had timed out at 900 s.
- **Exhaustive results with the new joiner** (`python3 attack/chain.py v gens rows --budget 30000`):
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

- **N_w = 3.** Two w-positions in one block (or a pair ±a in A) put a term of size about q² into the
  autocorrelation of the compressions y_k = u_k + q w_k, which only another such pair could cancel. So for N_w ≤ 3
  each w-position sits alone in its own non-A block. Without any decoupling assumption the system is then:
  stage 1, the compressions y_k (±1 with one entry q or q − 2 at 0 in each w-block) have Σ PAF_y = 0; stage 2, the
  z_k have Σ PAF_z(s) = R(s) = Σ_{w-blocks} (u_k(s) + u_k(−s)) and zero cross terms. This general search reproduces
  the u-quadruple counts of the decoupled N_w = 2 search exactly (9,520 and 209,536).
- **Precedent.** The counting identity is classical in another setting. Seberry Wallis (J. Combin. Theory A, 1976),
  Lemma 8, fills an orthogonal design of order 2ᵗ with J, J − 2I and a Paley-core block under
  a(v + 1) + b(v − 3) = 2ᵗ, which is our identity with a = N_w − N₋ and b = N₋. Seberry (Ars Combin., 1978),
  Theorems 6 and 7, give skew orders 2q(q − 1) and mn(n − 1) from orthogonal designs, and Wallis (1969, 1970) uses
  the same J, J − 2I, Paley-core blocks. In every version we found, the outer structure is an orthogonal design or
  a Kronecker product. Here it is a Goethals–Seidel quadruple over ℤ_m, and no located theorem produces the orders
  below (research-agent audit, 2026-10-12).
- **Prime powers.** Nothing above uses that q is prime: with χ the quadratic character of GF(q), q = pᵏ ≡ 3 (mod 4),
  the group becomes ℤ_m × GF(q)⁺ (non-cyclic for k > 1). The census has one such case, 513 = 19 · 27 (N_w = 3, N₋ = 2);
  `--field 3,3` lifts to ℤ₁₉ × ℤ₃³ and verifies with the `goethals-seidel-group` builder.

**`attack/paley_config.py`** solves the reduced problem exactly. Stage 1 lists the u-quadruples (or, in the general
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

- **1140 was within classical reach.** The audit of 1140 found that Seberry's theorem (Ars Combin. 6 (1978),
  Thm 7: an OD(mn; 1, m, mn − m − 1) and amicable designs of order n give a skew order mn(n − 1)) produces it with
  m = 3, n = 20 from an OD(60; 1, 3, 56), which a Goethals–Seidel search over ℤ₁₅ finds in 0.1 s; we rebuilt that
  matrix and it verifies. The design was not recorded anywhere we looked, and every list had 1140 as unknown. In our
  terms this is the case where (q + 1)/4 divides m (285 = 15 · 19 with (19 + 1)/4 = 5). For 1012, 1380, 1564, 1860
  and 6004 it does not divide m, and none of Seberry's theorems applies.
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
    κ: Σ (Da* − aD*) + Re J′(Db* − bD*) + Im J′(Eb* − bE*) = 0; λ: Σ (Ea* − aE*) + Im J′(Db* − bD*) − Re J′(Eb* − bE*) = 0.
    δ + q·G is again the compression condition. A position with s = 0 is either χ-type (b = ±1, D = E = 0) or
    ψ-type (b = 0, D, E = ±1); both add 2 + 2q to the δ identity coefficient. Skewness forces A to be ψ-type at
    a = 0. If that were the only ψ-type position in all four blocks, the κ-equation would need
    |±Re J′ ± Im J′| = 1, which is impossible, so ψ-type positions occur elsewhere too.
  - Stage 2 has six channel patterns per position against two in the Paley case. A direct low-energy test on 1341
    (O = 45) needed 8 GB in the exact join and was stopped; a dedicated search is still to be written.


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
the search uses the two compressions as exact intermediate levels instead (`attack/cyclo_product.py`, `join3.c`):

1. stage 1: the compressions to ℤ_m, up to negation, translation and reversal of the blocks;
2. level B: the compression to GF(q) (summing over ℤ_m), whose autocorrelations on the classes C_t, t < e/2, are
   the group key;
3. level C: the autocorrelations at the shifts (b, C_c), 1 ≤ b ≤ (m − 1)/2; stage 1 + B + C is the full condition.

**Single-spike family.** One spike with s = ±2 and f = −sign(s) costs 4F − 4, which equals 4me iff F = m e + 1, i.e.
q = m e² + e + 1. For m = 3: q = 53 (e = 4, v = 159), 311 (e = 10, v = 933), 991 (e = 18, v = 2973), 1753 (e = 24,
v = 5259). For m = 9: 149 (e = 4, v = 1341), 331 (e = 6, v = 2979).

**933 = 3 · 311, e = 10 (order 3732).** The single spike leaves one compression quadruple up to symmetry, rows
(1, 61, 3, 1). The lists have 8,064 (skew block), 13,335,840 (spike block) and 16,003,008 members; level B groups
them into 512, 55,915, 57,892 and 57,917 groups. join3 found 24 solutions among the first 6,978 of about 6 · 10⁷
matching group quadruples (5 minutes on 3 threads): of the order of 10⁵ skew solutions in all. A Gaussian model over
the 14 independent equations conditioned on the compressions predicts 10^4.2 to 10^4.6 (skew) and 10^8.1 (plain).
The lottery estimate in §5 for the same subgroup was 10^0.3 with 10¹⁷ pair operations. This settles the order 3732,
which FRONTIER.md listed as unresolved; see [3732](../3732).
  - Since 3733 is prime and ≡ 1 (mod 4), Miyamoto's theorem gives H(4 · 3733) = H(14932) from H(3732). The repository's
    Miyamoto builder rebuilds and checks it (H Hᵀ = 14932 I) in about 4 minutes.

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
the ℤ_m coordinate) to H. For plain structures K_m may contain −1. For 1255, H = {h ≡ ±1 (mod 5), h a 10th power
mod 251} has order 50 and 33 orbits per block; this is the |H| = 50 structure that the lottery ranking of §5 rated
at log₁₀ E = −5.0. With F = 25 a spike costs 96 or 104 at position 0 and twice that on the orbits {±1}, {±2}, so the
identity (Σ = 200) forces one spike of each kind at position 0 of two blocks, and stage 1 has a single class (rows
1, −1, −53, −47). The lists have 1.3–1.6 · 10⁷ members, but level B splits them into about 10⁶ groups, so a complete
join is out of reach. Keeping the 8,000 largest groups of every list (about 8 % of the members) gave 32 solutions in
17,273 group quadruples (6 minutes); see [5020](../5020). An FFT-based Gaussian estimate with the K_m symmetry gives
log₁₀ E ≈ 1.5 for the whole structure, which the run shows to be far too low.

**Sweep of the plain structures for the open orders.** For every composite n of FRONTIER.md (47 orders) and every
product Z_m × GF(q), even index e with (q − 1)/e ≥ 3 and at most about eight spikes, and every cyclic multiplier
group K_m ≤ ℤ_m* (including those containing −1) whose lists have at most 2 · 10⁷ members per block, stage 1 was
run exactly: 239 structures. Only 1255 (e = 10, K₅ = {±1}) has a compression quadruple; 933 had been solved
already. Every other structure has no solution. Structures with larger lists (for instance 2973 = 3 · 991 with
e = 18, a single spike and 4.9 · 10⁴ patterns per position) are not covered.

**1081 = 23 · 47 (Paley, N_w = 2).** `paley_config.py 1081 --maxu 3` went through all 6¹¹ joint (u₁, u₂) options in
2,759 s: stage 1 is empty, so 1081 has no solution in the Paley product with N_w = 2 (the only admissible N_w).

**Skew sweep with a multiplier on ℤ_m (`attack/skew_kmult_sweep.py`, 2026-10-13).** Every census value v = m q with q prime,
e | q − 1 even and (q − 1)/e odd, and every cyclic K_m ≤ ℤ_m* without −1 such that the structure has at most 22
position orbits and at most 2·10⁷ members per block: stage 1 run for 210 structures, 8 with stage-1 classes, all 8
joined (`cyclo_product.py`, `--sub 8000`, K = 1 per class). One hit: **v = 1417 = 13 · 109**, e = 4, K_m = ⟨3⟩,
H = ⟨1405, 328⟩ of order 81: skew-Hadamard matrix of order 5668 ([skew/5668](../skew/5668)). The structure is tiny
(lists of 144, 7,776, 5,184 and 5,184 members, 12 stage-1 classes); an exhaustive join of all 12 classes gives
8 solutions in each of the classes 2, 3, 6 and 11 and none in the others, 1 s per class. No hit for 1313 = 13 · 101
and 481 = 13 · 37 (e = 4, K_m = ⟨3⟩), 759, 261 and 1333.

**1341 = 9 · 149 (skew, e = 4), partial.** `cyclo_product.py 1341 --q 149 --e 4 --config 0 --K 3 --threads 3 --bkey 3
--sub 8000` (level-B key extended by the shifts (3, C_c)): lists of 5,184, 10,077,696, 10,077,696 and 6,718,464 members,
of which the 8,000 largest groups kept 100 %, 13 %, 22 % and 3.3 %; 2,604,080 group quadruples, 2.0 · 10¹⁰ member
operations, 1,961 s, no solution. This covers about 10⁻³ of the configuration's product space and rules nothing out.

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
(and D symmetric for the first form). `attack/propus.py` searches multiplier structures exhaustively: A on the orbits
of ⟨H, −1⟩, B and D on the orbits of H (or of ⟨H, −1⟩), every row-sum pattern r_A² + 2r_B² + r_D² = 4v, exact
three-list joins with join4x, and every hit is rebuilt as a matrix and checked (H = Hᵀ, H Hᵀ = 4v I).

**Checks.** v = 13: thousands of solutions in both forms (orders 52, all verified). v = 29: with A and D symmetric
and B arbitrary there is no cyclic solution at all (both row-sum patterns (3, 7, 3) and (7, 3, 7) exhausted, pools of
2–3 · 10³ symmetric and 3.5–6.8 · 10⁷ arbitrary sequences), so the known symmetric matrix of order 116 is not a
cyclic instance of the first form. It is an instance of the second: the family of Di Matteo, Đoković and Kotsireas
(Spec. Matrices 3 (2015) 227–234, arXiv:1503.04226) with −1 entries A = {0, 4, 5, 6, 7, 9, 13, 16, 20, 22, 23, 24, 25},
B = C = {0, 1, 2, 5, 6, 8, 11, 12, 14, 16, 22}, D = {0, 1, 2, 3, 5, 8, 10, 13, 14, 15, 18, 22, 25} gives a symmetric
Hadamard matrix of order 116 in `symmetric_matrix_v3` (checked).

**Literature (survey of October 2026).** Every computer-found symmetric matrix of the census is of the second kind:
one non-repeated block symmetric, the repeated block B = C neither symmetric nor skew (except v = 43, where it is the
skew Paley set), the fourth block arbitrary (Seberry–Balonin's propus construction, Australas. J. Combin. 69 (2017),
and Balonin–Đoković's "propus difference families", arXiv:1803.08787). Both blocks A and D symmetric is known only for
v = 13 and Turyn's series (arXiv:1708.05098). Cyclic propus searches are exhaustive for odd v ≤ 63 only. Since the
census, Đoković (arXiv:2508.17141, August 2025) proved that a symmetric Hadamard matrix of order q(q + 1) exists for every
prime power q ≡ 3 (mod 8), which settles 4 · 189 (q = 27) and 4 · 473 (q = 43); the paper also builds q = 67, which
gives 4 · 1139, an order on the census list. Nothing settles 65, 89, 101, 107, 119, 133, 149, 153, 163 or 167; the
symmetric orders 668, 716, 892, … remain open (the matrices of Alpöge et al. are not symmetric).

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
main's order 3964 use Q = ℤ₂ × ℤ₂ with T the three nontrivial characters (n ≡ 3 mod 4).

**The family Q = ℤ_q, T = all characters.** Here Γ = 4q δ_K: the targets are −4q on K \ 0 and 0 off K, the border has
size 4q and X, R̂, Ĉ are Hadamard of order 4q. For q = 3 it applies to every n divisible by 3 (|G| = n − 3), in
particular to n ≡ 1 (mod 4), where Alpöge's family does not exist. The compressions to ℤ₃ satisfy
Σᵢ PAF_{yᵢ} = 12 δ; for even |K| they are three single ±2 entries in three different blocks (4 · 27 · 8 = 864
patterns). `attack/bordered.py patterns` (C solver `attack/border3.c`) finds a border for all 864 patterns (4 min).
`attack/bordered.py validate` assembles complete examples: orders 36 (G = ℤ₃ × ℤ₂) and 60 (G = ℤ₃ × ℤ₄ and
ℤ₃ × ℤ₂ × ℤ₂), each with H Hᵀ = nI. Negation on the ℤ₃ factor seems unusable as a symmetry: the structures tried have no solution (n = 15, 21, 27,
exhaustive) or no row pattern whose lattice classes contain the targets (n = 573, all single generators).

**2148 = 4 · 537 in this family.** G = ℤ₅₃₄ = ℤ₆ × GF(89) is forced (t = 3). As a cyclotomic product (multiplier
classes of index 8, F = 11) the compression identity becomes Σ over the 24 positions of (11 s² + 2 f s) = 96: exactly
two spikes, ±23 at (0, 0) and ∓23 at (0, 1) of the block whose ℤ₃-compression is zero (`attack/bgs.py`, one stage-1
class). No further multiplier is compatible: −1 on ℤ₃ forces equal compression values at ±1 and breaks the delta
condition, and twisted multipliers force s ∈ 4ℤ at the spikes. The lists then have 1.2 · 10¹¹ members per block, too
many for the joins of §16. (The diagonal Gaussian estimate of 10⁷ given here earlier ignores the strong correlations
between the classes; the full-covariance estimate is 10^13.2.) §19 solves this structure: order 2148. The other orders 537 − t with a usable quotient (t = 9: |G| = 528, t = 15: 522)
give much smaller estimates. `bgs.py scan` (cyclic G = ℤ_m × GF(r), m ≤ 30, t ≤ 30, at most four spikes) found no
other order with a feasible stage 1 among the FRONTIER orders up to n = 1423.

**Q = ℤ₂² × ℤ₃, T = (nontrivial on ℤ₂²) × (all on ℤ₃), t = 9.** This product of Alpöge's set and the ℤ₃ set applies to
every n ≡ 9 (mod 12), i.e. to all open n ≡ 1 (mod 4) divisible by 3, and allows Alpöge's order-3 automorphism α of
ℤ₂² coupled with multipliers of K. Gaussian estimates (`attack/bgs_gauss.py`, summed over all row patterns): 789 (K = ℤ₅ × ℤ₁₃,
σ = (α, ×2, ×4), 84 orbits, 47 classes) 10^8.1; 789 with (α, ×4, ×2) 10^3.4; 537 at best 10^−3.1 (G = ℤ₂⁴ × ℤ₃₃);
669, 1437, 2661 below 10^−10. The border problem becomes: two orthogonal bases of the 36-dimensional space V_T made of
"fibre-character" vectors (one signed character of ℤ₂² on each of the 12 fibres (block, a)), with all mutual inner
products ±2 (in units where the norm is 12) for the single-spike compressions, i.e. a pair of unbiased weighing
matrices W(36, 12) with this fibre structure. The twisted-product ansatz (rows h_j(f) e_{π_f(u)} from a Hadamard
matrix of order 12) cannot work: it would need four columns of a Hadamard matrix of order 12 with constant product.
A meet-in-the-middle search over all 6¹² candidate columns found none for random bases of this kind. Open.

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
The other U-classes are (0, C_t) for t < 4, (1, 0) and (1, C_t) for t < 8, with two linear relations from stage 1
(12 independent equations); the W-condition has 11 independent classes, and its values lie on a parity lattice of index
2¹¹ once the supports are fixed. Given the U's, the W of each slice is a sign vector on the zero cells of U with a sum
fixed by stage 1 (C(4, 2) = 6 choices for four zero cells), so each U-solution leaves a small exact four-list problem
of about 200 vectors per block. This is the four-list split of T-matrices (§6) applied inside one block.

**Calibration on n = 177** (G = ℤ₆ × GF(29), e = 4, F = 7: the same stage-1 pattern with spikes ±13):
- the structure holds exactly 26,095,424 solutions in its stage-1 class (`join3.c` on the full lists of 46,656 and
  20,736 members, 46 minutes); the full-covariance Gaussian estimate is 10^6.87, the diagonal one 10^4.17;
- 157,482,016 U-solutions (`join4x.c`), with 0.166 W-completions each on average;
- reconstructed solutions give Hadamard matrices of order 708 (H Hᵀ = 708 I, border from `border3.c`).

**The search for 537** (`attack/bgs_uw.py`, `attack/kt4.c`):
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
5. `python3 attack/bgs_uw.py search 89 8 --L 10000000 --sigmas 2 --maxu 300000 --threads 4 --seed 2`: 284,292
   U-solutions, **one completion**, verified as a Hadamard matrix of order 2148 (197 s on 4 shared cores); see
   [2148](../2148). The observed rate 3.5 · 10⁻⁶ matches the lattice-corrected estimate.

So the structure that §18 judged out of reach takes minutes once the ℤ₂-split separates the conditions: the U-problem
has 12 equations and very many solutions, and each U-solution leaves a W-problem small enough to solve exactly.

**Generalization** (`UWG` in `attack/bgs_uw.py`). For any even m = 2m′ the pairs (α, α + m′) give
PAF_x(β + j m′, d) = 2 (PAF_U(β, d) + (−1)^j PAF_W(β, d)), with PAF_U cyclic and PAF_W negacyclic over ℤ_m′, so the
targets of U and W are (t(β) ± t(β + m′))/4. A multiplier group K_m = ⟨k⟩ on ℤ_m maps pairs to pairs (k is odd), with
W changing sign when a slice lands on the second position of its image pair; an orbit that maps a slice to itself with
a flip forces W = 0 there. Validation: on small structures with K_m (n = 24 with ℤ₄ × GF(5), n = 30 with ℤ₄ × GF(7))
the U/W counts equal exhaustive counts over all x-quadruples (36 and 128), and for n = 177 the class `UWG` reproduces
the completion rates of `UW`.

**Other open orders in this shape.** `python3 attack/bgs_uw.py scan <the 148 n of FRONTIER.md>` looks at cyclic
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

**The tool** (`attack/cp_lowe.py`, `attack/cp_enum.c`).

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

**Sweeps.** `attack/cp_hunt.py` runs this estimate-then-search loop over every census composite and, in plain mode,
every composite n of FRONTIER.md (`attack/cp_sweep.py` records the estimates only).

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
- 1123 (e = 34, F = 33, rows 1, 1, 1, 67, one class): 10^6.63 for the class and 10^2.24 inside lists of
  6 · 10⁷, 10⁷ and 10⁷ members (all 131,072 members of the skew block). Predicted 10^11.8 pairs per solution. Running.

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
`attack/lp_lowe.py` implements it (stage 1 on ℤ_m with spike costs summing to 2(m − 1)e, PAF sums −2q; level B target
−2m, level C −2; exact hash join) and `attack/lp_hunt.py` runs it over the FRONTIER orders. It finds Legendre pairs
of length 143 (ℤ₁₁ × GF(13), e = 4, K₁₁ = ⟨3⟩: the twin-prime class). But the count is against it:
- An H-invariant pair has two blocks of O orbits each, so about 2 · 0.3 · O decimal digits of freedom, against about
  O/2 conditions that each cost log₁₀(√(4πl)/4). That gives log₁₀ E ≈ O (0.6 − 0.5 log₁₀(√(4πl)/4)), which is
  negative for l ≳ 320. Goethals–Seidel structures have twice the freedom per condition (four blocks of length n
  against two of length 2n − 1), which is why their FRONTIER estimates are 10²⁰ and more.
- Measured: l = 1501 (n = 751), ℤ₇₉ × GF(19), e = 6, K₇₉ of order 13, the only stage-1 class: blocks of
  7.2 · 10⁸ and 2.2 · 10⁸ members, Gaussian estimate 10^−24.8 (rank 20).
So the route is closed for these lengths except through algebraic families.

**Prime-power fibres** (`attack/cp_field.py`, `cp_lowe.py --field p,k`). The cyclotomic product ℤ_m × GF(p^k) lives
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

**The 1123 run.** With lists of 5.8 · 10⁷, 10⁷ and 10⁷ (10^2.23 solutions expected inside), the old `kt4` processed
72 of 6000 split values in 830 s without a solution before it was killed for memory (another job was running). It
is being rerun with lists of 3 · 10⁷, 10⁷, 10⁷ (10^1.98) and the new `kt4`.

**Bordered structures through the same machinery** (`attack/bcp.py`: `bgs.Structure` seen as a cyclotomic product,
level-B/level-C targets −(m/q) Σ_α Γ(α) and −Γ(b mod q)). On the 2148 structure (ℤ₆ × GF(89), e = 8, Q = ℤ₃) the class
estimate is 10^12.6 and the estimate inside low-energy lists of 10⁷ members per block is 10^3.98, so the direct
four-list search would find it without the U/W split. (In-list estimates need about 10⁶ samples per block when the
lists keep less than 10⁻⁴ of a block; with 6 · 10⁴ samples the 2148 figure came out as 10^−56.)
- Total Gaussian estimates (all row patterns and stage-1 classes; random H-invariant sequences, PAF at the orbits of
  ⟨H, −1⟩ on nonzero shifts) for 573: with Q = ℤ₃ (border 12) every structure with at most 64 orbits per block is below
  10^−0.6; Spence-type (t = 1) ℤ₅₂ × GF(11), e = 2, K₅₂ of order 3: 10^10.9 over 60 orbits, spread over very many
  stage-1 classes.
- Larger totals need character sets whose border problem is not known to be solvable: 573 with ℤ₈ × GF(71), e = 10,
  Q = ℤ₈, T = characters of order 1 and 8 (t = 5): 10^16.6 (55 orbits); 643 with ℤ₁₂ × GF(53), e = 4, Q = ℤ₁₂, T of
  orders 1, 3, 12 (t = 7): 10^16.3 (60 orbits); 599 with ℤ₁₆ × GF(37), Q = ℤ₈, orders 1, 4, 8: 10^7.6. For Q = ℤ₈ and
  T = {1} ∪ {order 8}, a ±1 row of R̂ is on each block either constant or antiperiodic (f(α + 4) = −f(α)), since its
  compression to ℤ₄ must be constant. Whether such borders exist is open here.

**Reach of the low-energy four-list search.** In the tail regime a list of L members drawn from a block of N members
keeps a share of about g·L/N of the solutions (g ≈ (4/3)^{r/2} for rank r and equal covariances; the factor is
(1/2)^{−r/2} for the two-list case). So the expected number inside the lists is E_class · Π_k g L_k/N_k, while the
join costs about |A| L_B + L_C L_D pairs. Enlarging a structure (bigger N, larger r) at a fixed list size therefore
lowers the yield; the useful structures are those whose blocks can be kept nearly whole, N ≲ 10¹⁰, with
E_class ≳ 10³. The census map (`census_reach.py`, necessary norm and row-sum condition only, blocks of 10⁵–10¹¹ per
class) finds such structures for 44 of the undecided composite v up to 1089 — 119, 225, 319, 325, 369, 429, 473, 475,
481, 539, 553, 575, 589, 609, 623, 627, 637, 657, 689, 693, 707, 725, 737, 777, 781, 783, 793, 817, 845, 851, 891, 897,
899, 925, 927, 933, 943, 949, 957, 989, 1073, 1075, 1081, 1089 — and for none of the others. Many of these fail the
autocorrelation part of stage 1 or have class estimates below 1 (e.g. 481: 10^1.5; 589: at most 10^−2.6), so the
method has a natural frontier in the census; most undecided v admit no multiplier structure of searchable size at all.
- Non-cyclic K_m (two generators, no twist, m ≤ 120) add a single census structure that passes the norm and
  row-sum condition with blocks of 10⁵–10^11.5 per class: 2223 = 117 · 19, e = 2, |K| = 9, 25 position orbits.


## 22. Bordered arrays over ℤ₄, pair-spin structures, and calibrations (October 13, evening)

**The skew 4492 run is negative.** `cp_lowe.py 1123 --q 1123 --e 34 --config 0 --L 3e7,1e7,1e7 --nsig 6000 --threads 3
--K 5` (seed 1): lists of 131,072 (the whole skew block), 2.9 · 10⁷, 1.0 · 10⁷ and 1.0 · 10⁷ members, 6000 split values
holding a predicted share of 0.0746 of the 10^1.98 solutions inside the lists, 7.4 · 10⁸ AB and 2.3 · 10¹¹ CD pairs in
4932 s: no solution. With about 7 expected, the prime-field estimate inside low-energy lists is optimistic by a
factor of at least 2.4 (95 %). For comparison, the four blocks of the known plain solution of 4492 (same H and rows
1, 1, 1, 67) lie at energy quantiles 2 · 10⁻⁴ to 1 · 10⁻³ of their plain blocks (ranks 5 · 10⁵ to 2 · 10⁶ of 2.3 · 10⁹),
well inside lists of 10⁷. That solution was found by a search that favours low energy, so this does not test the
estimate.

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

**Pair-spin structures** (`attack/pairspin.py`). With μ ∉ H, μ² ∈ H and K = ⟨H, μ⟩, a pair of blocks (X, X ∘ μ)
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

**Border problems over ℤ_q** (`attack/border_probe.py`). The rows of R̂ and the columns of Ĉ are ±1 vectors of
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
  ℤ_{n−3} (and over ℤ₄ × K), where the ℤ₂ × ℤ₂ family needs a non-cyclic group. A census of these structures for the
  FRONTIER orders n ≡ 3 (mod 4) is the next step.
- Q = ℤ₈, T = {1} ∪ {characters of order 8} (t = 5, the 573 lead of §21): the compressions have row sums (±4, ±2, 0,
  0) and an antiperiodic part with five ±1 cells and Σ NPAF = 0 off 0 (48,640 patterns). Random bases give S_R = ∅
  (each R̂ is compatible with about 12.5 % of all ±1 vectors of V_T per row). Unresolved.

**The large bordered totals of §21 were artifacts.** The totals of `btotal_scan.py` average over random sequences and
ignore the exact compression budget Σ (F s² + 2 f s) = 4(m − t)e of `bgs.stage1`:
- 573, ℤ₈ × GF(71), e = 10, T = {1} ∪ {order 8}: stage 1 is empty for K_m = ⟨3⟩ and ⟨7⟩ (the 10^16.6 structures);
  with K_m trivial (88 orbits per block) it is too large for the exact join.
- 643, ℤ₁₂ × GF(53), e = 4, T of orders {1, 3, 12} or {1, 4, 12}: stage 1 empty (10^16.3).
- 599, ℤ₄ × GF(149), e = 4 (T of orders {1, 4}): the budget is 16 while a spike costs at least 4F − 4 = 144, so stage 1
  is empty (the total was 10^3.8).

**Prime-power fibres.** `cp_hunt.py --fibre-only --estimate-only` over the 76 census values with a fibre GF(p^k):
GF(25) and GF(49) structures (F = 3) come out between 10^−67 and 10^−0.04 per class; the best lead is 845 = 5 · 13²,
ℤ₅ × GF(169), e = 8, F = 21, K₅ trivial, class 1 of the sampled stage 1, rows (1, 43, −3, −39): 10^4.54 for the
class and 10^1.48 inside lists of 5 · 10⁷, 2.5 · 10⁷ and 2.5 · 10⁷ members (blocks of 1.3–1.7 · 10⁹). Queued after 319.

## 23. Skew orders 1276 and 3380, a symmetric order 532, and calibrations (October 13, night)

**Skew-Hadamard matrix of order 1276** ([skew/1276](../skew/1276)). The 319 run of §21–22 (ℤ₁₁ × GF(29), e = 4,
F = 7, H = ⟨45⟩ of order 7, stage-1 class 0 of the sampled three, rows 1, 25, 19, 17; lists of 31,104 (skew block),
60,466,176 (the whole second block), 9,966,434 and 10,012,823 members; estimate 10^1.64 inside the lists, share 0.09
in 6000 split values) found a solution at the 321st split value of the resumed run (the first run had searched 267
values; 131 recurred), about 450 distinct split values in all. With about 4 solutions expected over 6000 values and the
densest values first, this is in line with the estimate for composite v (unlike the prime-field 1123 run of §22).
Census: undecided; Cati–Pasechnik Table 5: 319(3); the audit found no earlier skew 1276 (its README).

**Skew-Hadamard matrix of order 3380** ([skew/3380](../skew/3380)), the lead of §22: ℤ₅ × GF(169), e = 8, F = 21,
H = ⟨w⁸⟩ of order 21, stage-1 class 1, rows 1, 43, −3, −39. Lists of 78,400 (skew block), 50,002,409, 25,105,989 and
25,135,175 members (estimate 10^1.48 inside them); four solutions at the 26th split value, 732 s in all on three
shared cores. The group is non-cyclic; the certificate uses the `goethals-seidel-group` builder with the automorphism
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
- **The Gaussian estimate misses by a factor 10⁹ here.** `propus_rank.py` gives log₁₀ E = −7.33 for exactly this
  structure and row pattern (pools 10^3.0, 10^6.5, 10^6.4; 12 independent classes), against 72 solutions (about two
  up to the multipliers and reversals that preserve the structure). The other row patterns of the structure (also
  about 10^−7.3) have none. The calibration below already showed factors of 10³–10⁷; propus structures must be
  searched exhaustively, not ranked.
- Sweeps (`--maxo 30`) are now running over the symmetric census values 89 < v ≤ 400 and 400 < v ≤ 1000
  (`$HADAMARD_WORK/logs/sym_queueA.log`, `sym_queueB.log`). v = 89 has no structure with a non-empty pool pattern
  within the cap.

**Symmetric Hadamard matrices: the propus estimate is far too pessimistic.** `attack/propus_rank.py` (Gaussian estimate
of the propus structures of `propus.py`, full covariance and lattice of sampled class values) gives, for the structures
in which Balonin and Đoković found propus difference families: v = 67, H of order 3 (A symmetric, D free) 10^−2.9;
v = 109, H of order 3, 10^−6.6; v = 151, H of order 5, 10^−6.5. Solutions exist in all three (six, three and one
families in their Table 1).
- v = 65 (the smallest open symmetric order, 260): every cyclic propus structure scores at most 10^−4.7 (H of order 3,
  25 orbits); the exhaustive sweep `propus.py 65 --sweep --maxo 30 --maxpool 1e8` (all structures up to 30 orbits,
  207 row patterns, including the 25- and 26-orbit structures beyond the earlier cap of 24) finds no solution.
- Larger structures (v = 163 with H of order 3: 28/55/55 orbits, 10^−12.1; v = 193: 10^−15.5) need a low-energy
  three-list search, not written yet.

**Parallel tempering for the spin structures** (`attack/spin_pt.c`, same input as `spin_ls.c`; replicas at geometric
temperatures, swap or pair-flip moves in O(C), replica exchange after every sweep). On 631 it needs 5–23 s per
solution (tabu: about 1 s); on 919 none in 10 CPU-minutes with 24 replicas between T = 20 and 20,000. It does not
beat the tabu search with restarts.

**ℤ₄ bordered families (§22) for 599, 643 and 647** (`attack/z4triage.py`, both character sets T, every multiplier
group with 14–44 orbits): 30, 437 and 192 structures, stage 1 empty in every one.

**Census hunt** (`cp_hunt.py skew`, composite census values up to 627 in the order of the list, `--mine 0.5`): beyond
301, 319 and 553 (already run) no structure reaches the threshold; 481 (ℤ₁₃ × GF(37), e = 4) has 10^1.46 for its
best class but only 10^−3.95 inside lists of 10⁷; most other values have empty stage-1 samples or blocks above the
enumeration cap (627: 1.15 · 10¹¹ members). Stopped at 627 to free the cores for the symmetric sweeps.

## 24. Symmetric and non-cyclic sweeps (October 13, late night)

**Tools.**
- `attack/propus_group.py`: the propus search over non-cyclic groups G = ℤ_p² × … × ℤ_m (group-developed blocks,
  R[g, h] = [g + h = 0]; the arrays of §17 hold over any abelian group). Structures: subgroups of GL(2, p) of any order
  (Singer, diagonal, torus, unipotent, Borel) on each split prime times any subgroup of ℤ_m*. Check: symmetric
  matrices of order 100 over ℤ₅² (Singer subgroup of order 3).
- `attack/propus.py --plain`: Goethals–Seidel quadruples (A, B, B, D) on H-orbits for ordinary Hadamard matrices
  (hits in `gs_propus_hits.jsonl`); `--maxwork` skips row patterns whose pair work exceeds a cap, the two smallest
  lists now form the pairs, and `--only-freeD` sweeps only the variant with A symmetric and B, D free (the form of
  every computer-found symmetric matrix of the census, §17).
- `attack/skew_group_sweep.py`: every skew structure over a non-cyclic group (`skew_group.py`, odd-order subgroups)
  with 8–30 orbits and every row pattern, exactly (find4b) when the pair work is at most 3 · 10¹⁰ and otherwise by the
  k-tree of `skew_ktree.run` with 2000 split values (it finds the 1100 solutions at the 810th).
- Speed: `Structure.class_values` is one matrix product per block of rows (7 times faster, identical values) and
  propus pools are enumerated by meet in the middle (`pool_exact`; 400 times faster at 26 orbits). The 133 pattern
  that holds the 532 solutions takes 26 s instead of 175 s.

**Symmetric Paley products** (`attack/propus_paley.py`). For v = m q with q ≡ 1 (mod 4) a prime power and
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
four-block lottery of §5.

## 25. Compressed propus searches, D-optimal designs and Legendre pairs (October 13, night)

**`attack/propus_chain.py`: compressions make large symmetric structures cheap.** For d | v the compressions of the
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
- Running: all composite symmetric census values v ≤ 1000, every structure with at most 44 orbits (freeD form);
  the composite open orders of FRONTIER.md with `--plain` and 20–44 orbits.

**D-optimal designs and Legendre pairs** (`attack/two_circ.py`, `attack/two_chain.py`). Two ±1 sequences with
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

**Williamson matrices of orders 65 and 67** (`attack/williamson.py`; the only odd orders below 71 for which
Kharaghani, Mohammadian, Tayfeh-Rezaie and Zaitsev, arXiv:2605.08661, list Williamson matrices as not yet found):
four symmetric ±1 sequences on the orbits of ⟨H, −1⟩, all row patterns, exact four-list joins. Checked on n = 13, 25
and 37. Running.

## 26. Local search without symmetry, and the skew hunt beyond 627 (October 13, late night)

**Why unstructured searches are tempting.** A crude Gaussian count for cyclic families without multipliers (each
equation P(s) = 0 met with probability about 4/(σ√(2π)), σ² the variance of P(s); about 10 bits for the row-sum
pattern) gives log₂ E ≈ v (2.5 − ½ log₂(0.63 √(7v))) − 10 for propus triples (A symmetric, B and D free), that is
about 2³³ triples for v = 65 and 2³⁵ for v = 89, positive up to v ≈ 370; for skew Goethals–Seidel quadruples
log₂ E ≈ v (3.5 − ½ log₂(0.63 √(5v))) − 13, about 2¹³⁴ for v = 89. Under the same count a multiplier group of order 2
leaves about 2⁻⁴ propus triples for v = 65: symmetry discards almost all of the expected solutions, and for the
smallest open symmetric order (260) and skew order (356) only the unstructured space is large. The question is
whether a solution can be *found* there.

**It cannot, by local search** (tabu search, every move evaluated exactly, restart after 20,000 steps without a new
best; `attack/propus_ls.c`, the general `attack/gs_ls.c` for the types skew / symmetric / free with weights, and
`attack/orbit_ls.c` with `orbit_ls.py` for multiplier-orbit structures). Steps per solution:

| v | 25 | 29 | 33 | 35 | 37 | 41 | 45 |
|---|---:|---:|---:|---:|---:|---:|---:|
| propus (A symmetric, B, D free) | 3.0 · 10⁴ | 3.8 · 10⁵ | 9.4 · 10⁵ | 3 · 10⁶ | 5.1 · 10⁶ | 1.0 · 10⁸ | none in 9 CPU-min |
| skew Goethals–Seidel | | 5 · 10³ | | | 1.0 · 10⁶ | | none in 7 CPU-min |

Every solution is checked (the propus ones as matrices: v = 29 gives a symmetric Hadamard matrix of order 116 with
H = Hᵀ and H Hᵀ = 116 I). The cost grows by about a factor 1.7 per unit of v for propus and 2 for skew quadruples,
which puts v = 65 near 10¹³ steps (about three CPU-years at the rate reached for v = 41) and v = 89 far beyond.
Two refinements do not help:
- **Fixed compressions** (`attack/comp_ls.c`, `comp_margins.py`): for v = p q every block is kept on the arrays with
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
solution. The v = 99 and 123 D-optimal runs that died in the DP are being repeated.

**Skew census hunt beyond 627** (`cp_hunt.py --trivial-k`, then all K_m). Of the 410 remaining composite census values
in 629..2499, 213 have cyclotomic-product structures within the member caps (477 structures, 89 with K_m trivial). The
trivial-K structures have large class sizes (F = 11 … 45) and their sampled stage 1 is empty except for 1565
(ℤ₅ × GF(313), e = 8, F = 39: 10^1.93 for the class but 10^−2.45 inside lists of 10⁷, blocks of 10^9.2), 1313
(10^−1.86) and 927 (a block of 3.8 · 10¹¹ members, above the enumeration cap). With F > me a spike alone exceeds
the budget Σ y² = 4v, so large-F structures are infeasible as soon as F ≳ me. The structures with K_m nontrivial
are running.

**Operational note.** The container schedules by autogroup: `nice` acts only inside a session. Priorities between
the background searches are set with `echo N > /proc/PID/autogroup`.

## 27. Bordered propus: symmetric Hadamard matrices of order 4(m + 1) from groups of order m (October 14)

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

**Order 260 (v = 65, the smallest open symmetric order).**
- ℤ₆₄, direct joins: every multiplier structure with at most 25 orbits is empty (both row patterns).
- ℤ₆₄, compressions: ⟨31, 33⟩ (25 orbits) has no compatible ℤ₁₆-compression triple for either pattern (ruled out);
  ⟨17⟩ (17/32 orbits): pattern (0, 0, 2) has no compatible ℤ₈ triple, pattern (2, 0, 0) has 640 (pair work
  5.2 · 10¹⁰, being completed); ⟨31⟩ (25/33) running. ⟨33⟩ (25/48 orbits) is out of reach: 55,616 compatible ℤ₈
  triples with pair work 3 · 10¹⁸, and H acts trivially on ℤ₁₆ and ℤ₃₂, so no finer compression stays small.
- Non-cyclic groups of order 64, cyclic automorphism subgroups, at most 30 orbits: ℤ₂ × ℤ₃₂ (6 structures),
  ℤ₄ × ℤ₁₆ (12), ℤ₂² × ℤ₁₆ (14), ℤ₂ × ℤ₄ × ℤ₈ (21), ℤ₈², ℤ₄³ (17), ℤ₂³ × ℤ₈, ℤ₂² × ℤ₄² — no solution. The
  2-generated subgroups (≤ 30 orbits) and the quotient chain (`propus_gchain.py`, 26–46 orbits) are running.

**All census values v ≤ 1000** (197 values): ℤ_{v−1} with at most 30 orbits (direct) and 31–44 orbits
(compressions), and the 319 non-cyclic groups of order v − 1 whose endomorphisms number at most 2 · 10⁶, are
queued (`border_queue.sh`, `bgroup_census.sh`).

## 28. A symmetric Hadamard matrix of order 876, and the lottery of small product structures (October 14)

**Order 876 = 4 · 219** ([symmetric/876](../symmetric/876)). The compressed propus sweep of the composite census
values (`propus_chain.py v --sweep --mino 1 --maxo 36`, smallest structures first) found it in ℤ₂₁₉ with H = ⟨4⟩ of
order 9: H fixes the ℤ₃ coordinate of ℤ₃ × ℤ₇₃ and acts on ℤ₇₃* as the eighth powers, so the structure is the
cyclotomic product ℤ₃ × GF(73) of index 8 (27 orbits for B and D, 14 for A, 13 classes). Row pattern (3, 17, 17):
compressions to ℤ₇₃ give 768 compatible triples (pair work 3.4 · 10⁷ against 6.9 · 10⁹ for the direct join); the
573rd holds 18 solutions, and all 768 together hold 288, in 16 triples (the 16 cosets of H in ℤ₂₁₉* permute them).
Verified by `verify.py`, by `check_matrix` and by an independent rebuild; reproduced from scratch with the same
sequences. Novelty audit (its README): 219 is undecided in the Đoković–Seberry census of 2025-05-28 and in no series
of it; no paper, database or 2026 deposit we could read has a symmetric matrix of order 876, and every classical
route fails by arithmetic; the 2018 Balonin–Seberry tables (n < 400) could not be read.

**The Gaussian model misses by 10^10.6.** `propus_rank.py` gives 10^−8.1 solutions for this structure and pattern
(pools 10^3.0, 10^6.9, 10^6.9; 13 classes, full covariance), against 288. For 532 the miss was 10^9.6 (§23). Both
structures are products ℤ_m × GF(q) in which H meets GF(q)* in a cyclotomic subgroup: for 876, H is the eighth powers
and fixes ℤ₃; for 532 (ℤ₇ × GF(19)), H = {(h, φ(h))} is the graph of a homomorphism φ: ℤ₇* → ⟨7⟩ ≤ GF(19)* of order 3,
so its orbits on ℤ₇* × GF(19)* are indexed by GF(19)* (an index-18 cyclotomic product with a twist). There the class values
are linear in a few quadratic invariants of short ℤ_m-sequences (the f-layer autocorrelations and the
cross-correlations of the class sequences, weighted by cyclotomic numbers), so they sit on a sparse lattice that a
Gaussian with the right covariance badly underweights near 0. Practical consequence: rank nothing by the Gaussian in
these families; exhaust them, smallest first. Running on the same principle:
- the compressed propus sweep of the composite census values ≤ 1000 up to 36 orbits (now with priority), then
  37–44 orbits, then the composite values 1000 < v < 2500 up to 44 orbits;
- the bordered sweeps of §27 over ℤ_{v−1} and the non-cyclic groups of order v − 1;
- for skew orders, the exact chain (`screen_chain.py 3600 36 0 --order orbits`) over the 102 census-screen survivors
  with at most 36 orbits that no chain run had touched (§17–§22 ranked them by the Gaussian, which is unreliable
  here).

## 29. ℤ_q-bordered propus: symmetric Hadamard matrices of order 4(m + q) (October 14)

**Construction.** Let G be abelian of order m, φ: G → ℤ_q a surjection with kernel K, and A (symmetric), B, D ±1
functions on G with

    PAF_A + 2 PAF_B + PAF_D = 4m at 0,  −4q on K \ 0,  0 off K.

Let P be the propus array of §17 over G, y_X(c) = Σ_{φ(g)=c} X(g) the compressions to ℤ_q, Mc the propus array over
ℤ_q of (y_A, y_B, y_D) (Mc² = 4q I), and R̂ a Hadamard matrix of order 4q with R̂ Mc R̂ᵀ = −4q X for a ±1 matrix X.
With R[:, (i, g)] = R̂[:, (i, φ(g))],

    H = [[X, R], [Rᵀ, P]]

is a symmetric Hadamard matrix of order 4(m + q): P² = I₄ ⊗ (4(m + q) I − 4q [g − h ∈ K]) and RᵀR = I₄ ⊗ 4q [φ(g) = φ(h)]
cancel; R Rᵀ = (m/q) R̂ R̂ᵀ = 4m I; compressing a group-developed block (also its transpose and its product with the
reflection g ↦ −g) over the fibres of φ gives the ℤ_q-developed block of the compression, so R P = (R̂ Mc)∘φ and
X R + R P = 0 from X R̂ = −R̂ Mc; X² = 4q I. q = 1 is §27. The row sums satisfy a² + 2b² + d² = 4q, all even when
v = m + q is odd. A census value v is thus reached through every abelian group of order v − q with a ℤ_q quotient,
for each odd q dividing v.
- **Compressions.** z = y/2 is integral with PAF_zA + 2 PAF_zB + PAF_zD = q δ on ℤ_q (z_A symmetric). q = 3: 48
  patterns, single spikes (z_A = ±δ₀ with one B spike, or one B spike and one D spike). q = 5: 120 patterns, and in all
  of them z_B = 0, so B̂ vanishes on the characters of G that are trivial on K.
- **The border problem** is a clique problem: 4q mutually orthogonal ±1 rows r_i (first entry +1) with
  r_i Mc r_jᵀ = ±4q for all i, j. q = 3: 2048 candidate rows, all 48 patterns solved at once (`propus.qborder_solve`;
  Goethals–Seidel arrays of ℤ₃-circulants also work). q = 5: 91,648 candidate rows, mean degree ≈ 3,200; a plain
  depth-first search stalls, a colouring-bounded clique search on each vertex neighbourhood
  (`attack/qborder_clique.c`, Tomita's MCQ) solves the patterns met so far in about 150 s (the first after 34 start
  vertices). Goethals–Seidel borders do not exist for the single-spike ℤ₅ pattern (z_A = δ₀, z_D = 2δ₀): none among
  1536 block conjugates × 120,000 GS quadruples. A second version avoids the global adjacency (random adjacent pairs
  (s, t), MCQ on N(s) ∩ N(t), |N(s) ∩ N(t)| ≈ 250): every ℤ₅ pattern met so far is solved in about a second.
  Simulated annealing over the signed column permutations of a fixed Paley matrix of order 20 does not work (cost
  stuck far from 0). q = 7: 17.4 million candidate rows, |N(s)| ≈ 370,000, |N(s) ∩ N(t)| ≈ 15,000; 329 start pairs on
  the pattern (z_A, z_B, z_D) = (δ₀, δ₀, (1, 1, 0, 1, 0, 0, −1)) found no clique (open; a longer run is going). Hits
  without a border are kept in `qborder_pending.jsonl`.
- **Checks.** ℤ₃: every one of the 576 triples on ℤ₆ gives a symmetric Hadamard matrix of order 36; the chain finds and
  verifies orders 108 (ℤ₂₄), 132 (ℤ₃₀) and 36 on ℤ₂ × ℤ₃; ℤ₁₂ and ℤ₁₈ have no triple (exhaustive, ℤ₁₂ also by brute
  force). ℤ₅: ℤ₁₀ and ℤ₂₀ have none (exhaustive), ℤ₃₀ has several (six compression patterns), and each pattern gives a
  verified symmetric Hadamard matrix of order 140. ℤ₇: sequences exist on ℤ₁₄ and ℤ₂₈ (55 compression patterns), but
  no border is known. Existence depends on m, as for §27.
- **Tools.** `propus.py --qborder q`, `propus_chain.py --qborder q` (per-class targets in stage 1 and in join4x),
  `propus_bgroup.py --qborder q [--phi c₁,…]` (automorphisms that preserve ker φ), `verify.py` construction
  `propus-qbordered` (fields q, quotient = coefficients of φ, R), `border_certify.py` (re-solves the border in the
  verifier's convention: the reflection g ↦ −g, where `propus.py` uses the back diagonal).
- **Literature.** We have not seen this bordering in the propus papers we read (Seberry–Balonin 2017,
  Balonin–Đoković arXiv:1803.08787, Đoković arXiv:2101.05429); it is the symmetric analogue of the bordered
  Goethals–Seidel arrays of §18 with T = all characters of ℤ_q, with the extra requirement that the border be R̂ and R̂ᵀ.
- **Runs.** The 117 census values divisible by 3 (`qborder_pool.sh symQ3.txt … 3`: direct joins up to 30 orbits, the
  chain for 31–44) and the 98 divisible by 5; v = 65 (order 260, the smallest open symmetric order) on ℤ₆₀ (chain,
  1–44 orbits) and on ℤ₂ × ℤ₃₀ (`propus_bgroup.py 2,30 --qborder 5 --phi 0,1 --pairs`). On ℤ₆₀ every structure up to
  24 orbits is empty (direct joins), and so far the chain has found no compatible compression triple in any structure
  except one of 30 orbits (60 triples, no solution).

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
