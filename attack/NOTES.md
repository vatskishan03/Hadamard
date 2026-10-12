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
- **Triage** (`skew_triage.py`, `skew_rank.py`) uses the Gaussian model of §2 with the skew A pool. It
  is only a rough guide here.
  - Under: 209 with ⟨26⟩ was predicted at 3 solutions and has 288.
  - Over: 201 with ⟨22⟩ (|H| = 11, rows 1, 1, 19, 21) was predicted at 10³ and has none
    (exhaustive `find4b`).
  - Under: 377 with ⟨16⟩ was predicted at 0.1 and gave 64 in the first τ.
- **Matching.** `match4.c` is `search/match.c` with four distinct lists: left pairs (A, largest pool),
  right pairs from the other two. With lists of very different sizes, the yield per τ grows with the
  smaller pair list, so the clamp S is set as loose as memory allows (pair lists ≤ 2.5·10⁷).
  Pools beyond 4·10⁶ are replaced by uniform samples (`--maxpool`), at a cost of L/N per list.
- **New orders**: 836, 1036, 1148, 1236, 1364, 1508, 1852, 1876, 2404, 3676, i.e.
  v = 209, 259, 287, 309, 341, 377, 463, 469, 601, 919. Each took seconds to minutes on 4 cores.
  - 259 used sampled lists: 3.9·10⁶ of 3–4·10⁷ per pool.
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
  - 817 has no structure with predicted E > 10⁻³.
- **General local search.** The same C engine handles plain skew GS: B, C, D are stacked into one
  vector with a block-diagonal form, and swaps stay inside a block (`spin_skew.py gsearch`). It
  solves 209/⟨26⟩ in under a second.
- **No hits:**
  - 201 with ⟨22⟩ and 301 with ⟨4⟩ (both exhaustive; 301 was predicted at 4);
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
