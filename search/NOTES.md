# Search notes: structures, obstructions and searches

Working notes on the searches for the unresolved orders of [FRONTIER.md](../FRONTIER.md) and for the skew and
symmetric orders, in the order the work was done. Everything here is either proved below or measured with the code in
this directory. [NOTES.md](../NOTES.md) collects the structural results.

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
  - 721 (|H| = 17, μ = 46, predicted 200) is still open after 6·10⁶ moves; the best residual is two classes off by 4.
  - 817: an earlier version required μ³ = 1 and found no structure with E > 10⁻³. That was too narrow, because μ³ ∈ H
    suffices: the blocks B, B∘μ, B∘μ² close up since B is H-invariant. With μ³ ∈ H, H = ⟨11⟩ (order 21) and μ ∈ {6, 9,
    17, 24} give spin structures with 8 classes and model E ≈ 10³·². Tabu searches on them are running.
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
  - **Correction.** The first version compared the norms |X̂(χ)|² after rounding them to integers. For d ≥ 5 they are
    not rational integers (they lie in ℚ(ζ_d)⁺), so that test was not exact; only its d = 3 verdicts were sound. Every
    verdict below was re-checked with the integer test.
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
- **Second-tier sweep** (pair work 1.5·10⁹–5·10¹⁰, ≤ 37 orbits). The first list of verdicts came
  from the rounded-norm test and is **withdrawn**. The lists are being recomputed with the integer test,
  and the structures are covered exactly by the compression-chain sweep below.
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
estimate in §5 for the same subgroup was 10^0.3 with 10¹⁷ pair operations. This settles the order 3732, which
FRONTIER.md listed as unresolved; see [3732](../hadamard/3732).
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
   [2148](../hadamard/2148). The observed rate 3.5 · 10⁻⁶ matches the lattice-corrected estimate.

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

**The 1123 run.** With lists of 5.8 · 10⁷, 10⁷ and 10⁷ (10^2.23 solutions expected inside), the old `kt4` processed
72 of 6000 split values in 830 s without a solution before it was killed for memory (another job was running). It
is being rerun with lists of 3 · 10⁷, 10⁷, 10⁷ (10^1.98) and the new `kt4`.

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
The census map (`census_reach.py`, necessary norm and row-sum condition only, blocks of 10⁵–10¹¹ per class) finds such
structures for 44 of the undecided composite v up to 1089 — 119, 225, 319, 325, 369, 429, 473, 475, 481, 539, 553,
575, 589, 609, 623, 627, 637, 657, 689, 693, 707, 725, 737, 777, 781, 783, 793, 817, 845, 851, 891, 897, 899, 925,
927, 933, 943, 949, 957, 989, 1073, 1075, 1081, 1089 — and for none of the others. Many of these fail the
autocorrelation part of stage 1 or have class estimates below 1 (e.g. 481: 10^1.5; 589: at most 10^−2.6), so the
method has a natural frontier in the census; most undecided v admit no multiplier structure of searchable size at all.
- Non-cyclic K_m (two generators, no twist, m ≤ 120) add a single census structure that passes the norm and row-sum
  condition with blocks of 10⁵–10^11.5 per class: 2223 = 117 · 19, e = 2, |K| = 9, 25 position orbits.


## 22. Bordered arrays over ℤ₄, pair-spin structures, and calibrations (October 13, evening)

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
  ℤ_{n−3} (and over ℤ₄ × K), where the ℤ₂ × ℤ₂ family needs a non-cyclic group. A census of these structures for the
  FRONTIER orders n ≡ 3 (mod 4) is the next step.
- Q = ℤ₈, T = {1} ∪ {characters of order 8} (t = 5, the 573 lead of §21): the compressions have row sums (±4, ±2, 0,
  0) and an antiperiodic part with five ±1 cells and Σ NPAF = 0 off 0 (48,640 patterns). Random bases give S_R = ∅
  (each R̂ is compatible with about 12.5 % of all ±1 vectors of V_T per row). Unresolved.

**The large bordered totals of §21 were artifacts.** The totals of `btotal_scan.py` average over random sequences and
ignore the exact compression budget Σ (F s² + 2 f s) = 4(m − t)e of `bgs.stage1`:
- 573, ℤ₈ × GF(71), e = 10, T = {1} ∪ {order 8}: stage 1 is empty for K_m = ⟨3⟩ and ⟨7⟩ (the 10^16.6 structures); with
  K_m trivial (88 orbits per block) it is too large for the exact join.
- 643, ℤ₁₂ × GF(53), e = 4, T of orders {1, 3, 12} or {1, 4, 12}: stage 1 empty (10^16.3).
- 599, ℤ₄ × GF(149), e = 4 (T of orders {1, 4}): the budget is 16 while a spike costs at least 4F − 4 = 144, so stage
  1 is empty (the total was 10^3.8).

**Prime-power fibres.** `cp_hunt.py --fibre-only --estimate-only` over the 76 census values with a fibre GF(p^k):
GF(25) and GF(49) structures (F = 3) come out between 10^−67 and 10^−0.04 per class; the best lead is 845 = 5 · 13²,
ℤ₅ × GF(169), e = 8, F = 21, K₅ trivial, class 1 of the sampled stage 1, rows (1, 43, −3, −39): 10^4.54 for the class
and 10^1.48 inside lists of 5 · 10⁷, 2.5 · 10⁷ and 2.5 · 10⁷ members (blocks of 1.3–1.7 · 10⁹). Queued after 319.

## 23. Skew orders 1276 and 3380, a symmetric order 532, and calibrations (October 13, night)

**Skew-Hadamard matrix of order 1276** ([skew/1276](../skew/1276)). The 319 run of §21–22 (ℤ₁₁ × GF(29), e = 4, F = 7,
H = ⟨45⟩ of order 7, stage-1 class 0 of the sampled three, rows 1, 25, 19, 17; lists of 31,104 (skew block),
60,466,176 (the whole second block), 9,966,434 and 10,012,823 members; estimate 10^1.64 inside the lists, share 0.09
in 6000 split values) found a solution at the 321st split value of the resumed run (the first run had searched 267
values; 131 recurred), about 450 distinct split values in all. With about 4 solutions expected over 6000 values and
the densest values first, this is in line with the estimate for composite v (unlike the prime-field 1123 run of §22).
Census: undecided; Cati–Pasechnik Table 5: 319(3); the audit found no earlier skew 1276 (its README).

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
- **Gaussian estimate (corrected in §29).** `propus_rank.py` gave log₁₀ E = −7.33 for exactly this structure and row
  pattern (pools 10^3.0, 10^6.5, 10^6.4; 12 independent classes), against 72 solutions. That version dropped the
  lattice factor of the class values (about C log₁₀ 4); with it the estimate is 10^0.5, and 72 is a factor of about
  10^1.4 above it. The other row patterns of the structure have none. [The conclusion first drawn here, that propus
  structures cannot be ranked, rested on the faulty estimator.]
- Sweeps (`--maxo 30`) are now running over the symmetric census values 89 < v ≤ 400 and 400 < v ≤ 1000
  (`$HADAMARD_WORK/logs/sym_queueA.log`, `sym_queueB.log`). v = 89 has no structure with a non-empty pool pattern
  within the cap.

**Symmetric Hadamard matrices: the propus estimate seemed far too pessimistic** [it was the estimator; §29].
`search/propus_rank.py` (Gaussian estimate of the propus structures of `propus.py`) gave, for the structures in which
Balonin and Đoković found propus difference families: v = 67, H of order 3 (A symmetric, D free) 10^−2.9; v = 109, H
of order 3, 10^−6.6; v = 151, H of order 5, 10^−6.5. Solutions exist in all three (six, three and one families in
their Table 1). These numbers lack the lattice factor (§29) and are too low by roughly C log₁₀ 4 for C classes.
- v = 65 (the smallest open symmetric order, 260): every cyclic propus structure scores at most 10^−4.7 (H of order 3,
  25 orbits); the exhaustive sweep `propus.py 65 --sweep --maxo 30 --maxpool 1e8` (all structures up to 30 orbits, 207
  row patterns, including the 25- and 26-orbit structures beyond the earlier cap of 24) finds no solution.
- Larger structures (v = 163 with H of order 3: 28/55/55 orbits, 10^−12.1; v = 193: 10^−15.5 with the faulty
  estimator, both far above 1 after the correction of §29) need a low-energy three-list search, not written yet.

**Parallel tempering for the spin structures** (`search/spin_pt.c`, same input as `spin_ls.c`; replicas at geometric
temperatures, swap or pair-flip moves in O(C), replica exchange after every sweep). On 631 it needs 5–23 s per
solution (tabu: about 1 s); on 919 none in 10 CPU-minutes with 24 replicas between T = 20 and 20,000. It does not
beat the tabu search with restarts.

**ℤ₄ bordered families (§22) for 599, 643 and 647** (`search/z4triage.py`, both character sets T, every multiplier
group with 14–44 orbits): 30, 437 and 192 structures, stage 1 empty in every one.

**Census hunt** (`cp_hunt.py skew`, composite census values up to 627 in the order of the list, `--mine 0.5`): beyond
301, 319 and 553 (already run) no structure reaches the threshold; 481 (ℤ₁₃ × GF(37), e = 4) has 10^1.46 for its
best class but only 10^−3.95 inside lists of 10⁷; most other values have empty stage-1 samples or blocks above the
enumeration cap (627: 1.15 · 10¹¹ members). Stopped at 627 to free the cores for the symmetric sweeps.

## 24. Symmetric and non-cyclic sweeps (October 13, late night)

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

## 25. Compressed propus searches, D-optimal designs and Legendre pairs (October 13, night)

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
- Running: all composite symmetric census values v ≤ 1000, every structure with at most 44 orbits (freeD form);
  the composite open orders of FRONTIER.md with `--plain` and 20–44 orbits.

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
solution. The v = 99 and 123 D-optimal runs that died in the DP are being repeated.

**Skew census hunt beyond 627** (`cp_hunt.py --trivial-k`, then all K_m). Of the 410 remaining composite census values
in 629..2499, 213 have cyclotomic-product structures within the member caps (477 structures, 89 with K_m trivial). The
trivial-K structures have large class sizes (F = 11 … 45) and their sampled stage 1 is empty except for 1565
(ℤ₅ × GF(313), e = 8, F = 39: 10^1.93 for the class but 10^−2.45 inside lists of 10⁷, blocks of 10^9.2), 1313
(10^−1.86) and 927 (a block of 3.8 · 10¹¹ members, above the enumeration cap). With F > me a spike alone exceeds
the budget Σ y² = 4v, so large-F structures are infeasible as soon as F ≳ me. The structures with K_m nontrivial
are running.

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

**Order 260 (v = 65, the smallest open symmetric order).** Settled on October 17 by a quadruple with no multiplier
symmetry, decoded from Verma's quaternary Legendre pair of length 64 (NOTES.md §32); the structured searches below
had no chance to find it.
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

## 28. A symmetric Hadamard matrix of order 876, and small product structures (October 14)

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

**The Gaussian model seemed to miss by 10^10.6** [corrected in §29: the estimator dropped the lattice factor; the
corrected estimate is 10^0.6, a factor 10^1.9 below the 288 found, and 10^0.5 against 72 for 532]. `propus_rank.py`
gave 10^−8.1 solutions for this structure and pattern (pools 10^3.0, 10^6.9, 10^6.9; 13 classes, full covariance),
against 288. For 532 the apparent miss was 10^9.6 (§23). Both structures are products ℤ_m × GF(q) in which H meets
GF(q)* in a cyclotomic subgroup: for 876, H is the eighth powers and fixes ℤ₃; for 532 (ℤ₇ × GF(19)), H = {(h, φ(h))}
is the graph of a homomorphism φ: ℤ₇* → ⟨7⟩ ≤ GF(19)* of order 3, so its orbits on ℤ₇* × GF(19)* are indexed by
GF(19)* (an index-18 cyclotomic product with a twist). There the class values are linear in a few quadratic invariants
of short ℤ_m-sequences (the f-layer autocorrelations and the cross-correlations of the class sequences, weighted by
cyclotomic numbers), so they sit on a sparse lattice that a Gaussian with the right covariance badly underweights near
0. Practical consequence: rank nothing by the Gaussian in these families; exhaust them, smallest first. [§29 withdraws
most of this paragraph: with the lattice factor restored the excess is 10^1.4–10^1.9, not 10^10, and the estimator
matches exact counts within ±0.3 on small v.] Running on the same principle:
- the compressed propus sweep of the composite census values ≤ 1000 up to 36 orbits (now with priority), then 37–44
  orbits, then the composite values 1000 < v < 2500 up to 44 orbits;
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
  (z_A, z_B, z_D) = (δ₀, δ₀, (1, 1, 0, 1, 0, 0, −1)) found no clique (open; a longer run is going). Hits without a
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

**The Gaussian estimator was wrong, and most of the apparent yield was an artifact.** `propus_rank.estimate` dropped
the lattice factor of the class values whenever the sampled covariance was singular, and for propus it always is: the
row sums fix Σ_c |c| · val_c. The estimates were therefore too low by about Σ_c log₁₀ g_c ≈ C log₁₀ 4. Corrected (the
density on the hyperplane w^⊥, w = (|c|), times the covolume ∏ g_c · |w| / gcd(g_c w_c) of the lattice points there),
the estimator matches exhaustive counts on ℤ_v with trivial H within ±0.3 in log₁₀ (v = 7, 9, 11, 13, every row
pattern; one arithmetic zero at v = 7, rows (1, 1, 5), predicted 10^1.8). Revised figures:

| structure | old estimate | corrected | found |
|---|---|---|---|
| 876: ℤ₂₁₉, ⟨4⟩, rows (3, 17, 17) | 10^−8.1 | 10^0.6 | 288 |
| 532: ℤ₁₃₃, ⟨26⟩, rows (1, 1, 23) / (1, 11, 17) | 10^−7.4 | 10^0.5 | 72 (rows (11, 5, 19)) |
| v = 65, H trivial, rows (11, 3, 11) | 10^−8.3 | 10^10.9 | (not searchable) |

So the two finds sit in structures where about one solution was expected, with a structural excess of 10^1.4–10^1.9,
and the small structures (≤ 30 orbits) we swept mostly have expectations below one. The large ones are where solutions
are plentiful: the propus ranking can now be trusted to choose them, and the problem becomes search technology for
structures of 45–80 orbits (multi-level compression for composite v, low-energy lists for prime v), as it was for the
skew orders. The skew estimators (`cp_lowe.estimate_in_lists`, census screen) are separate code and are not affected
by this fix; they are checked separately.

**Calibration of the corrected estimator on the swept structures.** For v = 133 and 219 (all structures up to 30
orbits searched exhaustively) the corrected estimates stay below 10^0.9, and the two structures that hold solutions
are among the highest (133/⟨26⟩: 10^0.77, 72 found; 219/⟨4⟩: 10^0.57, 288 found; the best empty one, 133/⟨4⟩ with H
of order 9, 10^0.81). So in the swept range each structure holds at most about one expected solution, and a sweep
succeeds or fails by chance. v = 65 is the exception:

| v = 65 structure | orbits A, B/D | corrected estimate (best rows) | searched | found |
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

## 30. Solution orbits, and dedicated row searches: order 260 over ℤ₆₅ is out (October 14)

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
- 2-adic check (scratch `local_density.py`): for 65/⟨16⟩, rows (11, 3, 11), the fraction of pool triples that solve
  the class equations modulo 2^k, against the uniform 2^−(C−1) per level, is 0.53, 0.35, 0.90 for k = 1, 2, 3 (a
  factor of about 1/6) and 0 modulo 16. But modulo 16 the class-value maps are nearly injective (170,948 distinct
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
  - 183/⟨34⟩ (e = 12): 8 triples, B ≈ 4·10⁶ and D ≈ 4–6·10⁸ per triple, about 2·10⁹ probes each; running.
  - 213/⟨25⟩ (e = 14): about 3,000 rows per row sum, B ≈ 3·10¹⁰ before the filters per triple: beyond this join.

**Order 260 by the bordered families**, in expected orbits: ℤ₅-bordered on ℤ₆₀, ⟨49⟩ 0.37 and ⟨11⟩ 0.19 (both empty,
above), the rest ≤ 0.02; ℤ₁₃-bordered on ℤ₅₂, ⟨9⟩ 0.22 (20 orbits, but the border is a clique problem with rows of
length 52 that the solver of §29 cannot handle yet); bordered (q = 1) on ℤ₆₄, at most 10⁻⁵ (the character obstruction
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
expected orbits each) never appeared in any sweep log; they are being run with the current chain (no solution so far).

## 31. Symmetric Hadamard matrices of order m(m − 1), and orders 5700 and 732 (October 14)

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
  none for n = 9 (35 is not a² + 2b²). These existence results were published before (found by the audit of 17
  October 2026): Abuzin, Balonin, Đoković and Kotsireas (Inform. Control Syst. 2019, no. 5, 2–9, Table 1, "X1 skew,
  X2 = X3, and X4 symmetric") list such families for every odd v < 50 with 4v − 1 = a² + 2b² except 37, 43 and 47,
  for which their Table 3 has only families with X4 neither symmetric nor skew; for 37 there is none with X4
  symmetric (`sp_mitm.c`, exhaustive, §49). Their families of lengths 19, 25
  and 31 are equivalent to ours (a multiplier, a shift of A1, and for 31 a reversal and sign change of A1).
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
verify.py rebuilds it). The new amicable orders 52, 76, 100 (and 124 from n = 31) are ≡ 4 (mod 8) but not of the
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
triples. The corrected estimate of `propus_rank.py` for the structure was 0.69 expected solution orbits (|Γ| = 432),
0.06 for this row pattern. The remaining triples were not finished; the rerun of the fifth triple alone (commit
`913f7b5`) finishes it with exactly the two solutions.

**prow.c part keys (fix).** The D parts of `join3`/`join4` were chosen by the additive key
(F₀ − 3)/4 + 7 (F₁ − 3)/4 mod np of the first two x = 0 classes. Those class values span only a few hundred residues,
so with np in the thousands almost all D landed in a few parts: the 213/⟨25⟩ run (np = 1811, table sized for the
average part) stopped with "part 0 overflows its table" after its B passes. The key is now Σ_k μ_k (F_k − 3)/4 over
all e/2 classes with fixed pseudo-random μ_k (μ₀ = 1), which is still additive over the rows of D and splits into an
A part and a B part for a needed vector. The planted tests (join3 and join4) recover their solutions; 213/⟨25⟩ was
restarted with it.

**Skew audit (5700, 2652, 9900).** The skew audit (statuses in the READMEs) found no earlier construction for 5700 and
9900, and none for amicable orders 52, 76, 100: Seberry's tables of amicable orders 2ᵗq (1978, 1991, 1992) give
t = 3 for q = 13, 19, 25. The 1992 summary line "AV 2(q + 1), q ≡ 5 (mod 8); J. Wallis" is the |t| = 1 case
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
| 52 (new) | **1989 (3): done, [skew/7956](../skew/7956)** |
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
| 1989 | 39 | 3 | 52 (new) | {1, 16, 22} | three singletons | **2688 solutions → skew 7956** |
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
| 891 | 81 | 27 | 12 | ⟨4⟩, order 27 | no template | (order 9: running) |
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

## 32. The symmetric weighted plug-in in Balonin form; generic group engines (October 14, afternoon)

**The plug-in.** A symmetric OD(4t; μ, 4t − μ) xP + yN (P, N symmetric with complementary supports, P² = μI,
N² = (4t − μ)I, PN + NP = 0) and a symmetric Hadamard matrix of order h with 4t = μh give the symmetric Hadamard matrix
P ⊗ J_{h−1} + N ⊗ D of order μh(h − 1), D the core of the symmetric matrix (§31: the J terms cancel exactly when the
weight of P is μ = 4t/h). No amicable pair is needed, only a symmetric matrix of order h = 4k, so on paper the census
values v = μk(4k − 1) are in reach: 125 symmetric-undecided values in 139 splits, 122 in 134 without the values
already settled here (k = 1, h = 4 alone covers every multiple of 3), e.g. 475 (t = 25, μ = 5, h = 20), 891 (t = 81, μ = 27, h = 12), 1989 (t = 39, μ = 3, h = 52).

**Balonin form.** With A symmetric and B, D any circulants (or group-developed matrices),

    [ −A   BR   BR   DR ]
    [ BR   RD   A   −RB ]
    [ BR   A   −RD   RB ]
    [ DR  −RB   RB   A  ]

is symmetric and orthogonal whenever AAᵀ + 2BBᵀ + DDᵀ is scalar (the off-diagonal blocks cancel identically, also
with variables, provided both parts of A are symmetric). Splitting each block into an x part and a y part, the design
conditions are PAF(A_x) + 2PAF(B_x) + PAF(D_x) = μδ, the same for y with 4t − μ, and X_A + 2X_B + X_D = 0 for the
cross terms X(s) = Σ_g x(g)[y(g + s) + y(g − s)]. Summing the cross condition over s gives the row-sum screen
r_Aρ_A + 2r_Bρ_B + r_Dρ_D = 0 with r_A² + 2r_B² + r_D² = μ and ρ_A² + 2ρ_B² + ρ_D² = 4t − μ
(`search/symod_reach.py`). It removes 1989 (μ = 3), 693 (μ = 21), 1235, 1353, 2013, 2093 (μ = 1, as in §31) and most of the
n = 4 splits; it keeps 475, 891, 1089, 2185, 2375, 1617, 1749, 1815, 1947, 2277, 513, 657, 669 and 723.

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
too), μ = 1 skew-propus designs, and the ℤ₂₅ negative. Fixed on the way: the signed row-sum size of a skew pair is 0,
not its weight. verify.py's Theorem-7 "od" form takes an optional "group" (blocks over ∏ ℤ_n^k, R the negation).

**Theorem 7 with n = 4.** A skew OD(4t; 1, t, 3t − 1) gives a skew-Hadamard matrix of order 12t, so t = v/3 reaches
every undecided skew v divisible by 3 (172 of them). Exhaustively over ℤ_t (H trivial, Goethals–Seidel form):
- t = 3, 5, 9: designs exist;
- t = 7: impossible, since 7 is not a sum of three squares (Eades), as for every t ≡ 7 (mod 8).
- t = 11: a design exists (template from a capped sample of 2 · 10⁶ b quadruples).
- t = 13: running.
In every design found the b parts sit in a short window around 0, where periodic and aperiodic autocorrelations agree.
They are four short ternary complementary sequences with u₁ antisymmetric. For t = 11:
- u₁ = (1, 0, −1) on positions (−1, 0, 1);
- u₂ = −(δ₀ + δ₁ + δ₂);
- u₃ = (1, −1, −1) on (−1, 0, 1);
- u₄ = (−1, 1, −1) on (−2, −1, 0).

Their aperiodic autocorrelations at shifts 1 and 2 sum to 0, and the weights sum to 2 + 3 + 3 + 3 = 11. So b parts
of any weight are plentiful (from complementary ternary sequences); the c part is the hard half.

The census sweep (`n4sweep.py` in the work directory) covers every cyclic H of odd order with at most 15 orbits for
the 120 targets with t a sum of three squares: 3-orbit groups (squares) and most 5–9-orbit groups have no b
template at all, and no c part was found in any run (results in the table below when complete).

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
- Orbit annealing (`osa.c`) finds the 1716 design. For 1089 with H of order 3 (55 H-orbits, about 52 variables per block) the energy
  never drops below about 10⁵ (random level), at T₀ = 40 or 3000. This agrees with §26: local search is not
  competitive above about 40 free positions per block.

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
with row sum −1 at all, nA = 0). The Kiyasu routes are empty in this form for the largest groups: 235 (ℤ₄₆₉, |H| = 33), 277 (ℤ₅₅₃,
|H| = 39) and 149 (ℤ₂₉₇, |H| = 45). verify.py has a construction "two-block" (kinds skew and symmetric, any abelian
group) ready for certificates.

## 33. Plug-in balances, two weight-1 symmetric variables, and cyclotomic plugs (October 14, evening)

**Theorem 7 outside Goethals–Seidel form (negative).** Eades' condition (§31) only rules out 4-type designs. A
skew-type OD(n; 1, m, c) of order n ≡ 4 (mod 8) in any form gives a rational family, which exists only if
⟨1, m, c, mc⟩ ≅ 4⟨1⟩, i.e. the Hilbert symbols (−m, −c)_p are 1 at every odd prime p (the transpose is the canonical
involution of the quaternion algebra (−m, −c) generated by the two skew variable matrices, so the trace form is a
multiple of its norm form). For the twelve Theorem-7 census targets where Eades fails this test fails as well:
627 (p = 13), 759 (7), 1045 (13), 1221 (7, 29), 1235 (3, 41), 1323 (47), 1683 (7), 1815 (151), 2145 (7, 17), 2185 (109),
2375 (3, 79). So no array of any kind (two-block, non-abelian groups, bordered) can give these designs. The test passes
for 429, 475, 693, 891, 957, 1089, 1353, 1617, 1749, 1947, 1989, 2013, 2277, 2365, 2409 and 2457, exactly Eades'
list. (`/tmp` script `hilbert.py`; the same test for symmetric designs is not valid: their transpose is an
orthogonal involution, and the symmetric OD(76; 1, 75) of §31 would fail it.)

**No symmetric OD(4n; 1, 1, 4n − 2) for n > 1.** In the symmetric plug-in (§31) the balance s_J r + s_P (r − 4) = s_core
with weight 1 on J, weight 1 on J − 2I and the rest on cores gives r = 2n + 1, i.e. symmetric Hadamard matrices of
order 4n(2n + 1) (census: 903, 1953) from a symmetric OD(4n; 1, 1, 4n − 2) and any symmetric Hadamard matrix of order
2n + 2. No such design exists. Let P₁, P₂ be its weight-1 matrices (symmetric signed permutations with underlying
involutions τ₁, τ₂, which commute because P₁P₂ = −P₂P₁) and N the rest, N anticommuting with both. Evaluate
NP₁ + P₁N at (i, τ₂(i)) and NP₂ + P₂N at (i, τ₁(i)): with x = N[i, τ₁τ₂(i)] and y = N[τ₁(i), τ₂(i)] they give
x s₁(τ₂i) = −s₁(i) y and x s₂(τ₁i) = −s₂(i) y, while P₁P₂ = −P₂P₁ gives s₁(i)s₁(τ₂i) = −s₂(i)s₂(τ₁i). So x = 0. Since N
vanishes only on the supports of P₁ and P₂, every i is fixed by exactly one of τ₁, τ₂. On X₁ = Fix(τ₂), P₂ is a
diagonal ±1 matrix D₂ with D₂(τ₁i) = −D₂(i), and NP₂ = −P₂N forces N[i, j] = 0 whenever D₂(i) = D₂(j); but N[i, j] ≠ 0
for all j ∉ {i, τ₁(i)} in X₁, so |X₁| ≤ 2, likewise |X₂| ≤ 2, and the order is at most 4 (where the design exists).
This also explains the empty Balonin-form searches for this type (§31).

**Weighted cyclotomic plugs (negative).** Generalising Theorem 7's plug, a skew-type OD(4t; 1, s₂, …, s_k) with
I ⊗ A + Σ X_i ⊗ Y_iR, A group-developed of skew type a and Y_i of symmetric y_i over GF(p), is skew-Hadamard of order
4tp iff PAF_a + Σ s_i PAF_{y_i} = 0 off 0 (R-type symmetric blocks commute with each other and are amicable with any
group-developed A, so nothing else is needed). With a, y_i constant on cyclotomic classes of order e (−1 ∈ C_{e/2})
the condition is e/2 linear equations in the class PAF vectors; a meet in the middle over all sign patterns
(`/tmp` script `wplug.py`, primes p ≤ 700, e = 4, 6, 3- and 4-variable types of orders 4–28 with square products)
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

**Faster spin tabu (`spin_ls2.c`, `SPIN_LS=spin_ls2 spin_skew.py csearch ...`).** The same moves and acceptance
rule as `spin_ls.c` (identical trajectories for a given seed: the best residual after 60 restarts on 1657 is 26 496
in both), but the swap scan evaluates, for each + orbit i, all partner orbits j at once in single precision with
the j-dependent terms (2g_cj + 2M_cjj and the row i of M_c) laid out contiguously, so the compiler vectorises it; the
exact integer objective is recomputed after every move. 2.75 times faster: 60 restarts of 3000 steps on 1657 take 4 s
instead of 11 s. Measured restarts per solution on 919 (|H| = 17, μ = 8, 9 classes): 1361, 3753, 8848 (steps 3000),
328 and 1695 (steps 8000), 14 543 (steps 1000), i.e. one to five minutes per solution on a loaded machine, against
twenty minutes before; on 631 (|H| = 7, μ = 2, 7 classes) 7–50 restarts (under a second). Two variants did not
help on 631: jumping A to the best of the 2^NP pool every 10 or 50 steps, and restarting half the time from the best
state so far perturbed by six swaps. The census queue (817 and 721 on their 8-class structures, then 1657, 1801 and
1951 in two-hour slices) runs with it.

## 34. Spence's skew 4(1 + q + q²) construction, made explicit; the missing ingredient is a cyclic CDS (October 14, evening)

Spence (Canad. J. Math. 27 (1975) 555–560) builds a skew-Hadamard matrix of Goethals–Seidel type of order
4(1 + q + q²) from a cyclic projective plane of order q and complementary difference sets (CDS) in ℤ_v,
v = 1 + q + q². The plane part, as implemented in `search/spence_plane.py`:
- Q = q², N = Q² + Q + 1 = v w with w = q² − q + 1 (coprime to v). D = {i mod N : Tr(αⁱ) = 0} is the Singer
  (N, Q + 1, 1) difference set of PG(2, Q); it is fixed by the multiplier q.
- f(r) = #{d ∈ D : d ≡ r mod v} takes values 0, 1, 2 (the fibres are orbits of the subgroup of order w, arcs).
  Every difference of D occurs once, so Σ_r f(r) f(r + s) = w for s ≢ 0, and u = f − 1 has PAF_u = w − 2(Q + 1) + v = 0
  off 0: a perfect ternary sequence.
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

## 35. Three skew blocks in one multiplier orbit (negative for the census) (October 14, night)

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

## 36. Spin structures: an exact fibre test, and the landscape of the tabu search (October 14, night)

**Fibre test (`spin_fibre.py`).** Over a composite ℤ_v, take a shift s in a proper subgroup dℤ_v. Its
autocorrelation splits over the d cosets a + dℤ_v (the fibres): P_X(s) = Σ_a P_{X_a}(s), where
X_a(y) = X(a + dy). The fibre orbits under a group G of multipliers are independent, and the fibre g·a contributes
the term of fibre a at the pulled-back shift g⁻¹s. So a DP over fibre orbits gives exactly:
- the attainable values of P_A(s) for skew K-invariant A (fibres a and −a tied by the skew condition);
- the attainable values of T_B(s) = Σ_j P_B(μʲs) for H-invariant B, together with the row sum of B.

The spin condition P_A(s) + T_B(s) = 0 needs the two sets to meet on every subgroup class. The DP agrees with brute
force (all 256 A of 721) and contains every sampled T_B.

The two composite hexagonal structures that had been in the tabu queue are **infeasible**:
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

## 37. Near-Williamson quadruples by exhaustive search: order 191 and the complex orders (October 14, night)

**The route.** Cooper–Wallis needs Williamson-type matrices of order w and T-matrices of order t to give H(4tw).
Near-Williamson quadruples (A, B, C, D) — circulant, B, C, D symmetric, ΣPAF = 4wδ — give the Williamson-type
matrices (AR, B, C, D) ([near-williamson](../near-williamson)), and since three blocks are symmetric they are also
Szöllősi quadruples, i.e. complex Hadamard matrices of order 2w. For w = 191, T-matrices of orders 5, 9 and 15
(from base sequences BS(3, 2), BS(5, 4), BS(8, 7)) would settle three entries of FRONTIER.md at once: 3820
(n = 955), 6876 (n = 1719) and 11460 (n = 2865). §6 and williamson-type/README.md had ruled out T-matrices of order
191 and the bordered form on ℤ₁₉₀; the plain form on ℤ₁₉₁ had not been ranked.

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
6876 and 11460 now needs a quadruple without multiplier symmetry. 183 (complex order 366) and 213 (426) are open complex orders; their best
structures are empty too. Single-precision sums of spectra are compared with 4w + 1, so rounding cannot reject a
quadruple; the 191 run was repeated with this margin (the first run used 4w + 10⁻⁶, also empty).

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

## 38. Anti-automorphisms of the Theorem 7 designs: no symmetric 7956 this way (October 14, night)

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
  equal, and P₁² = I needs k_i² = 1. When Stab(A) has odd order (here {1, 16, 22}) the only such k is −1: reflections, as
  in `gs_anti.py`. The b part (one entry at 0 in each of
  blocks 2–4, so ±R in every off-diagonal block of B) forces c_i + c_j to be the same for all i ≠ j, hence one common
  centre, and P₁ is block-monomial in ±R: an involution π and signs ε, 80 candidates. For every sign pattern of the b
  entries, the signed union-find of P₁GP₁ = Gᵀ (G the design in the normal form of `od13.c`) leaves 48 consistent
  configurations (computed for t = 5, 7 and 39; the form does not depend on t). In every one of them one of W₂, W₃,
  W₄ is symmetric and the other two are tied, W_l(g) = ±W_k(±g). The bc condition of §31 asks that the symmetric
  parts of W₂, W₃, W₄ add to 0 on every pair {s, −s}; the symmetric block contributes ±1 and the tied pair S_k + S_l
  with S_l = ±S_k, which is 0 or ±2. So **no Goethals–Seidel OD(4t; 1, 3, 4t − 4) has a block-monomial involutive
  anti-automorphism** when the stabiliser of the first block has odd order, for any t. For the design of the certificate, the direct test over all reflection
  configurations with arbitrary centres also finds none (`od_anti.py cert`).

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

## 39. Spin structures over ℤ₇² × ℤ₄₃: the order 8428 = 4 · 2107 is out of reach this way (October 15)

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

## 40. What the rankings say about the remaining symmetric and FRONTIER leads (October 15)

**The symmetric chain sweep beyond v ≈ 440.** The compressed propus sweep of §28 had covered the composite census
values up to 435 (a coverage map rebuilt from the chain logs) before it was stopped; it is running again on the 85
composite values from 445 to 995 (orbits 1–36). The corrected Gaussian (§29) puts every ranked structure with at most
45 orbits for v ≥ 440 below 10^−1.5 (best: 589/⟨28⟩ 10^−1.5, 451/⟨5⟩ 10^−2.3, 445/⟨16⟩ 10^−3.1), so the sweep mostly
produces exact negatives. The estimates were low by 10^1.4–10^1.9 for the two hits of this family (532, 876), so a
hit is not excluded.

**163 (order 652).** The estimates disagree with the exact count where the count is available. With H = ⟨38⟩ (order
9; A on 10 orbits, B and D on 19, 9 classes) the corrected estimate is 10^1.47 for the row patterns (17, 1, 19) and
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

## 41. Paley products for the complex orders: order 330 (October 15)

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
(7, 6). N_w = 5 (34 skeletons) has no solution. Among the 2455 skeletons of N_w = 6 the run had found three verified
solutions when it ran out of memory (it kept every stage-1 orbit as a tuple; now they are packed bytes, and the
exhaustive count is being rerun). This gives the complex Hadamard matrix of order 330 ([C330.json](../complex-hadamard/C330.json), rows
7, 9; 1, 23). The order is blank in Seberry–Yamada's Table 11.2, 329 = 7 · 47 is not a prime power, and no
quaternary Legendre pair of length 164 is known (neither family of Jedwab–Pender covers it).

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

**A bug in `lowe_sz.py`, and order 426.** Its ordered split values were computed for the pair (A, B), while kt4
joins (A, C) and (B, D): at N = 183 every split value missed the lattice of A + C (residue 18 against 2 modulo 20 on
the first merge class), and the first join examined no pair at all. With the fix, the first structure ranked after it
by the mid-size triage, N = 213 = 3 · 71 with H = ⟨37⟩ of order 7 (33 orbits for C and D, 16 classes; rows 13, 17;
13, 15; estimate 10^2.2 in the structure and 10^0.95 inside lists of 6435, 6435, 3 · 10⁷ and 3 · 10⁷ kept from pools
of 4.4 · 10⁸), gave 72 solutions within the first split values (six minutes in all), all checked. This is the complex
Hadamard matrix of order 426 ([C426.json](../complex-hadamard/C426.json)); 425 = 5² · 17 is not a prime power and no
quaternary Legendre pair of length 212 is known. The 183 join is being repeated from its saved lists.

## 42. Compressions decide which row patterns live: orders 366 and 454 (October 15)

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

**Tools.** `lowe_enum.c` takes optional fibre-class constraints (per class, the weighted sum of the orbit signs over
it is fixed; exact reachability pruning per class); `search/lowe_comp.py` runs the low-energy k-tree on one
compression class at a time (classes reduced by translations by H-fixed points, reversal of C and D, and multipliers);
`search/comp_triage.py` triages every open N (plain and `--bordered`).

Also fixed: `kt4.c` leaked its hash table on every split value with fewer than 2¹⁸ pairs (the cause of the
171 crashes); allocations are now checked, and `lowe_sz.py --resume --sig-from k` continues a join.

## 43. Conditioned scans of every open complex order: order 430; twisted and non-cyclic conditioning (October 15)

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
  On the coarser quotient ℤ₃, 333 (|H| = 9, 45 orbits) has conditioned estimates 10^0.4–10^1.1 per pattern, spread over
  18–108 tuples with pools near 10⁷ each.

**Order 430.** N = 215 = 5 · 43, H = ⟨41⟩ of order 7 (trivial on ℤ₅), rows 15, 15; 11, 17: 450 feasible ℤ₅ tuples,
3 classes up to translations, reversals and multipliers, each with about 10^2.3 predicted solutions and pools of
1200, 1200, 2.4 · 10⁶ and 2.4 · 10⁶. `lowe_comp.py 215 41 --rows 15,15,11,17 --d 5 --only 0` lists the first class
completely (seconds) and kt4 found 5 solutions in its first split values (514 s in all, most of it the estimate).
All five rebuild as complex Hadamard matrices; the first is [C430.json](../complex-hadamard/C430.json). The plain
Gaussian of this pattern is 10^−0.8: the compressions put every solution into a few classes of tiny pools.

**Twisted structures** (`search/tsz_comp.py`). The full and half twists (B = A∘μ, D = C∘μ; A, B on the L′-orbits) have
compressions x_D = x_C ∘ μ̄, so C contributes P_x(s) + P_x(μ̄ s) to the compressed equation. Per block: exact
compressions (on the L- or L′-orbits for the hashed side), the feasible groups, and the lattice Gaussian of the twisted
class values. Calibration against exhaustive counts (`tsz.py test`):

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
symmetry), rows 2, 0, 0, 0: conditioned 10^7.84 against 4.6 · 10⁷ solutions (plain 10^7.02); N = 13: 10^4.18 per pattern
against 12,672. It does not see the obstruction that empties N = 17 over ℤ₂² × ℤ₄ (conditioned 10^6.4), which lives
on a finer quotient.

Also: `tsz.write_input` streams the (A, B) pairs in chunks (a 3 · 10⁷-pair list took 5 GB).

## 44. Order 342; finer quotients, live split values, and the estimates against the searches (October 15)

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
and (b) otherwise filters candidates by membership in both lattice cosets. On ℤ₉, tuple 0 of 171/⟨31⟩ (pools 1134,
756, 4.3 · 10⁵, 4.3 · 10⁵) has 38,086 live split values; all of them, 3.4 · 10⁸ + 2.2 · 10⁸ pairs, give no solution,
against 10^0.48 predicted.

**Calibration from the searches.** kt4 stops at its solution cap, so the N = 215 run says little about the size of
the class beyond its density: 8 solutions in the first 574 split values (1.5 · 10⁷ pairs per side), consistent with
the 10^2.29 predicted. N = 453 (pools 2520, 2520, 1.6 · 10⁷, 1.6 · 10⁷, 10^1.58 predicted): every split value was live
on both sides, and all 2 · 10⁵ of them (3.8 · 10⁹ pairs per side, about 9.5% of the A × C pairs, the inner part of
the solution law) gave no solution. Sampling the class sums of that tuple shows no pair or triple of classes that never vanishes jointly,
so no low-order obstruction explains it; the finer quotient ℤ₁₅₁ (H of order 15 acting on it) is not listable by the
scan and is the next suspect.

## 45. The character-bordered Szöllősi array: order 214 (October 15)

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

**A convention bug, caught by the verifier.** `spence_lib.AbGroup.index` reverses all axes of a 3-dimensional
coordinate array, so the first version of `tbsz.szarray` and `tbsz.build` developed every block as X[g, h] = x(g − h).
Its matrices were still complex Hadamard (the same construction for the reflected sequences, with the border solved
for the reflected compressions), but the certificate, read in the standard convention by `verify_complex.py`, failed.
`tbsz.devidx` now flattens the coordinates first, `certify` re-solves the border from the core sequences, and the
certificate verifies.

**Other N ≡ 3 (mod 4).** See [NOTES.md](../NOTES.md) §16: N = 167 is empty in both families on ℤ₄ × ℤ₄₁ and ℤ₂² × ℤ₄₁
(odd fibres of 41 points); N = 179 on ℤ₁₆ × ℤ₁₁ lives only for ⟨(1, 3)⟩ at 10^−1.9; N = 223 and 239 are being ranked.

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

## 46. The ℤ₆ character-bordered array: order 322 (October 15)

*Correction (October 16).* 322 was not open, nor were 330, 342, 430 and 454 of §§41–44: 4N − 1 is a prime power
for each, and Ito's Hadamard groups with Egan 2019, Theorem 4.5, give them ([NOTES.md](../NOTES.md) §18,
[PALEY.md](../complex-hadamard/PALEY.md)). The searches and their calibrations stand. The orders 214, 366 and 426
of §§41–45, and the other eight found by search, were constructed concurrently by Zaitsev et al.
([NOTES.md](../NOTES.md) §20).

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
cyclic structure in range is empty on ℤ₆, as is the first on ℤ₂ × ℤ₉ × ℤ₉. For 179 with t = 3 on ℤ₄ × ℤ₄ × ℤ₁₁ (Q = ℤ₂²,
H of order 10) the single feasible tuple of rows 2, 2; 2, 0 (complete lists 192, 192, 12096, 725760) has disjoint
merge-class sumsets: no solution. The rankings with t = 7 over Q = ℤ₂³ (167; fibres of 20) exceed memory in
`TComp.compressions` (8 fibre classes of up to 21 values each) and need a coarser conditioning quotient.


## 47. Refined compression tuples, and a biased lattice in the estimators (October 15)

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

**Biased lattice (fixed).** `comp_est.block_stats`, `gcomp`, `tsz_comp` and `TComp.stats` took the lattice of a block's
class values from `echelon((U[1:] − U[0])[:120])`, U = `np.unique` of the samples. Those are the lexicographically
first 120 vectors, which share their leading coordinates, so their differences can span a proper sublattice. On the
N = 167 structure below (fine level) this gave index 2 and index 6 sublattices for blocks A and B and once a lattice of
rank 17 instead of 18. A sublattice makes `in_lattice` fail where the target is attainable (a false −∞) and inflates
the covolume factor where it is not. `comp_est.diff_lattice` now takes up to 400 differences of distinct samples in
random order from a random base point; all four modules use it. Estimates printed before this fix can be off in either
direction: the ℤ₂ × ℤ₃⁴ structure of N = 167 rated 10^2.41 on ℤ₆ before and 10^1.2–10^1.4 after. Exhaustive results
are unaffected, but the v = 65 ⟨21⟩ explanation of §30 ("the target is off the coset of the lattice", from sampled
differences) has to be re-derived.

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

## 48. Paley plus three blocks over the census primes v ≡ 3 (mod 4) (October 17)

Đoković and Kotsireas (arXiv:1511.08734) build skew-Hadamard matrices of order 4v, for v ≡ 3 (mod 4) prime, from
Paley's block (−1 on the nonzero squares, PAF −1 off zero) and a (v; k₁, k₂, k₃; λ) difference family X₁, X₂, X₃:
three ±1 sequences with PAF sum +1 off zero. `skew_paley3.py` searches this class with the three blocks unions of
orbits of a subgroup H of the squares, of odd order d. That leaves m = (v − 1)/d orbits besides {0}, and m ≡ 2 (mod
4).

**Method.** `exhaust v d` works as follows.

- **Encoding.** The PAF at a shift class {±gᵏH} is a quadratic form cᵀNₖc in the orbit signs. The PSD is bounded by
  3v − 1, using Gauss periods.
- **Symmetry.** The sign of {0} is fixed in each block (negation). The first block is reduced modulo the squares, which
  permute the orbits and fix Paley's block.
- **Join.** The third block is joined by hashing.
- **Validation.** It finds the factors of [skew/9424](../skew/9424) (v = 31, ⟨2⟩) and [skew/13072](../skew/13072)
  (v = 43, ⟨4⟩), and v = 71 under the subgroup of order 5. Under the subgroup of order 17 it finds nothing for v = 239,
  where Đoković and Kotsireas used order 7.

**Row sums.** A block has row sum c₀ + d·e with e ≡ m (mod 2), and r₁² + r₂² + r₃² = 4v − 1. This removes most
subgroups: of the structures with 22 to 58 orbits over the census primes, eight survive.

**Exhaustive and empty.** Every census-undecided prime v ≡ 3 (mod 4) and every subgroup with 6, 10, 14 or 18 orbits
gives nothing. That is 97 structures over 65 primes:

- 49 with 6 orbits, 20 with 10 and 13 with 14 (`sweep 6 14`);
- 15 with 18 orbits, for v = 487, 523, 739, 883, 919, 991, 1063, 1423, 1567, 1747, 1783, 1999, 2143, 2251 and 2287
  (`sweep 15 20`).

None of these structures has 2, 4, 8, 12 or 16 orbits, so these four sizes cover the whole range up to 20. The other 39
undecided primes ≡ 3 (mod 4) have no subgroup in this range.

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

## 49. Symmetric matrices of order 4q: inputs, census omissions, and no odd multiplier (October 17)

**Inputs of the 4q family** ([../symmetric/PROPUS_MIYAMOTO.md](../symmetric/PROPUS_MIYAMOTO.md)). Of the 58
census-undecided prime powers q ≡ 1 (mod 8), 44 take a symmetric Hadamard matrix of order (q − 1)/2 from Paley I or II
and Sylvester, and 1753 takes symmetric/876. Three take published sporadic matrices, transcribed from the LaTeX sources
and rebuilt by `verify.py` in the propus array of `build_propus` (input kind `propus`): order 116 (Di Matteo, Đoković
and Kotsireas, arXiv:1503.04226v2, the first of their fifteen SDSs), 236 (Đoković, arXiv:2101.05429v2, the first PDF
for v = 59, X₁ = X₁′ ∪ −X₁′) and 604 (Balonin and Đoković, arXiv:1803.08787, Table 1, orbit representatives under
{1, 8, 19, 59, 64}). Every published solution of those papers checked (15 for v = 29, 6 for v = 59, 1 for v = 151).
The last ten values need a symmetric matrix of order 4m with m = 101, 107, 119, 149, 179, 239, 249, 269, 299 or 305,
all undecided.

**Seberry–Balonin with conference cores** (`sym_propus_sb.py`). Mathon's conference matrix of order 442 comes from
SageMath's `MathonPseudocyclicStronglyRegularGraph(2)` (29 s; the same matrix on a rerun), the one of order 226 from
the skew-Hadamard matrix of order 16 (Paley 8, doubled) as Q = S ⊗ S + I ⊗ J − J ⊗ I. With Yamada's circulants for
q = 881 and 449 they give 3532 and 1804. A mutated core is rejected by its SHA-256, and a changed ξ by the primitivity
test. The variant with 2n − 1 an odd square instead of a prime power would need new symmetric circulant pairs
(M, N) with M² + N² = s²I: among the census values it would reach only v = 1227 (n = 613, s = 35) and 2027 (n = 1013,
s = 45), and the lattice count for multiplier-invariant pairs, about 2 · (4c/√n)^k for k orbits, is far below one.

**No odd t multiplies the propus array** (`sym_tmult.py`). Cooper–Wallis multiplies four Williamson-type matrices by
T-matrices; the symmetric analogue would write X_k = Σ_j M_kj ⊗ Y_j (k, j ∈ {a, b, d}) with symmetric commuting (0, ±1)
M_kj of disjoint supports in each row and Σ_k c_k M_kj M_ki = t c_j δ_ij for c = (1, 2, 1). At the trivial character the
row sums give an integer M with Mᵀ diag(1, 2, 1) M = t diag(1, 2, 1), so det(M)² = t³ and t must be a square. Exhaustive
searches over the symmetric elements of ℤ₃, ℤ₅, ℤ₉ and ℤ₃² find no solution (the control t = 1 finds the 16 signed
permutations). Four distinct symmetric blocks do not help either: a symmetric orthogonal design of order 4 in four
variables would need four anticommuting real symmetric involutions of order 4, and their Clifford algebra M₂(ℍ) has no
module of dimension 4. Twisting by the reversal does not rescue it (`sym_arrays.py`). The blocks XR and RX = XᵀR are
symmetric for every circulant X, but they do not commute, so the Clifford argument does not apply to arrays built from
them. An exact symbolic check settles it: every symmetric Latin square of order 4 (96), every orientation and sign, each
of the four circulant variables symmetric (blocks ±X, ±XR) or arbitrary (±XR, ±RX), H² compared with (Σ XXᵀ) ⊗ I by
cancelling monomials. No array with four distinct variables exists, whichever variables are symmetric. With X₃ = X₂
there are 768 arrays when only X₁ is symmetric (Balonin's among them), 1152 when X₁ and X₄ are, 6144 when X₁, X₂ and X₄
are, and none when X₁ is arbitrary. So every symmetric Hadamard matrix of Goethals–Seidel type over ℤ_v needs a repeated
block and a symmetric block, as in the propus family. Twisting by any involutive multiplier does not help either
(`--klein`). For v with two odd prime factors, ℤ_v* contains u ≠ ±1 with u² = 1, and X P_u is symmetric when X is
invariant under the multiplier −u. Over all 35 type assignments in {1, −1, u, −u}, none of the arrays has four distinct
variables. The twisted repeated-block arrays that do exist (448 with D of type u) impose the same autocorrelation
equation as propus with extra symmetry, so they reach nothing new.

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
  A1 with row sum 7 and A1(0) = −1 (every A1 has such a shift) are looked up: no match (3 cores, 14 minutes). A
  rerun with the work split over the first two −1 positions of A1 and a Bloom filter in front of the table also finds
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
- **Afterwards (NOTES.md §31).** The repeated block, and with it the skew-propus quadruple, turned out to be
  unnecessary: Yamada's data in quaternary form give a symmetric matrix for any skew core. Of the eight values above,
  295, 679, 823, 1063 and 1495 are now built that way (skew inputs of orders 148, 340, 412, 532 and 748); 1879, 2215
  and 2359 wait for skew-Hadamard matrices of orders 940, 1108 and 1180.


## 50. The largest remaining symmetric classes: what fails (October 17, exploration round 2)

After NOTES.md §30–32 the census has 499 undecided values; the two largest algebraic classes are v ≡ 5 (mod 8) prime
powers (57 values) and v = q + 2 with q ≡ 1 (mod 8) a prime power (35 values). Four exploration agents attacked them
and the composite cases; none found a construction. Their exact findings (scripts kept only in the session scratchpad;
"proved" means a complete CP-SAT run or enumeration):

- **q ≡ 5 (mod 8), order 4q.** The propus-support V-part of NOTES.md §30 (symmetric V of order 4n, V² = qI, support
  [[I, J, J, 0], [J, 0, I, J], [J, I, 0, J], [0, J, J, I]]) exists for n = 2, 4 (q = 5, 9) and is infeasible for
  n = 3, 6 and 10 (q = 7, 13, 21; proved). Cyclotomic propus quadruples (Legendre-type A, D, B constant on quartic
  classes) exist only for q = 5 and 13 among the primes q ≡ 5 (mod 8) below 400, and the quaternary 2-circulant
  cyclotomic ansatz of order 2q only for q = 5, 13. A ℤ₂²-refined "skeleton" (two complementary-support symmetric
  weighing matrices plus a conference matrix of order q + 1, which builds symmetric matrices of orders 20 and 36) has
  no instance at q = 13 in any of its splits (proved, using that every 12 × 8 and 12 × 10 partial Hadamard matrix
  extends to H₁₂ and that Aut(H₁₂) is 5-transitive on columns). The known symmetric matrices of orders 4q in this
  class show no structure to generalise: 52 is Paley II, and the propus matrix of order 116 has an automorphism group
  of order 2.
- **v = q + 2, q ≡ 1 (mod 8).** In the quaternary template of NOTES.md §31 the cross-block orthogonality forces the
  core to be a unit times a real matrix, so a real symmetric conference matrix of order (v + 1)/2 is necessary and
  sufficient there; this is Seberry–Balonin's condition, and the complex Yamada–Kiyasu matrices with quaternary cores
  (q = 17, 73, 193 tested) have no anti-automorphism. Only 8 of the 35 values (235, 523, 731, 955, 1011, 1491, 1603,
  2379) have (v − 1)/2 a sum of two squares, so for the other 27 no real conference matrix of that order exists. A
  ℤ_{q+1}-circulant two-orbit family finds symmetric quaternary matrices at n = 5, 7, 9 by annealing (at n = 9 one
  that comes from a non-Paley symmetric conference matrix of order 38) and none for n = 15; orbit-reduced CP-SAT with
  multiplier groups proved infeasibility in many small structures, including n = 117 (v = 235) with four subgroups of
  order 12.
- **Composite q.** Theorem A's U-part, Yamada–Kiyasu's (M, N) and Seberry–Balonin's V-part are all equivalent to a
  symmetric conference matrix C of order q + 1 with a signed permutation R, R² = −I, RCRᵀ = −C (then C = [[M, N],
  [N, −M]] in a basis (e_x, Re_x)); Paley's matrices have one (x ↦ g/x, g a non-square). The three constructions go
  through for any such C (checked for random R at q = 5 to 41), so composite values need a conference matrix with R:
  40 values of Theorem-A type (65, 153, 305, …), 18 of Yamada–Kiyasu type (119, 223, 247, …) and the Seberry–Balonin
  values would follow; no non-Paley conference matrix with R is known at those orders.
- **Symmetrization of quaternary arrays.** The array [[X, Y], [−Y*, X*]] of a complementary pair (x, y) becomes
  symmetric by a monomial diag(Q₁, Q₂) exactly when XQ₁, X*Q₂ are symmetric and YQ₂ = −Q₁ᵀȲ; this gives the propus
  type (x symmetric, y real, which settled 260: NOTES.md §32), Hermitian or skew-Hermitian y with arbitrary x, and
  α-twisted pairs for an involutive automorphism α. All were checked exactly on small groups; no published pair of a
  remaining census order is of these types beyond Verma's.
- **Zaitsev's matrices.** None of QH130 … QH1526 of the remaining census orders is equivalent to a symmetric matrix
  (NOTES.md §31).

## 51. Conference matrices for the Seberry–Balonin and Yamada–Kiyasu routes; symmetric structures (October 18)

Brouwer's table (read 17 October 2026) records conference graphs on 117 and 153 vertices (J. Y. Z. Ho, and C. Wang
for 153, private communications of October 2026), and Verma has a symmetric conference matrix of order 86
(arXiv:2610.08893: a 12 × 12 array of circulants of order 7 with two border rows, automorphism group C₇ ⋊ C₃ that fixes
six blocks and permutes the other six in two 3-cycles). A symmetric conference matrix of order 118 would give, by
Seberry–Balonin's Corollary 1 with q = 233, the symmetric order 940 (v = 235) and, by the Miyamoto blocks of
hadamard/MIYAMOTO_BLOCKS.md, FRONTIER 2823; with an R-structure (NOTES.md §31) it would also give v = 119 by the
Yamada–Kiyasu route (skew input of order 60). None of these matrices is public, so we searched for our own (scripts in
the session scratchpad `conf/`).

- **Block-circulant search** (`bconf2.py`: ℤ_m acting with k regular orbits and F fixed points, every circulant block
  constant on the orbits of a multiplier group H, row-by-row meet in the middle on the exact correlation conditions,
  then the fixed rows). Complete, with no solution: order 118 with m = 29, k = 4, F = 2 and |H| = 4, 7, 14; m = 13,
  k = 9, F = 1 and |H| = 4, 6; m = 39, k = 3, F = 1 and H = ⟨4⟩, ⟨10⟩, ⟨5⟩, ⟨14, 22⟩. Order 154 with m = 17, k = 9,
  F = 1 and |H| = 8; m = 19, k = 8, F = 2 and |H| = 6, 9; m = 51, k = 3, F = 1 and H = ⟨2⟩. Calibration: the tool finds
  Paley's matrices (orders 6, 10, 14, 30) at once, but neither order 66 (m = 13, k = 5, F = 1, |H| = 2, 3, 4, 6) nor
  Mathon's order 46 (m = 15, k = 3, F = 1, |H| = 2, 4; m = 9, k = 5, |H| = 2) has a structure of this kind, so the
  2026 matrices need either the block-permuting groups of Verma's type or no multipliers at all. Smaller multiplier
  groups (|H| ≤ 3 on ℤ₁₃ with k = 9 and on ℤ₃₉ with k = 3, |H| ≤ 4 on ℤ₁₇, ℤ₁₉ and ℤ₅₁) exceeded the tool's time or
  memory.
- **Block-permuting groups** (`gconf.py`: G = ⟨ℤ_m, (σ, μ)⟩ with σ permuting the blocks, brute force when the free
  blocks have at most 18 sign bits): order 118, m = 29, k = 4, μ = 2 and σ a 4-cycle (12 bits) or two transpositions
  (14 bits), and order 66, m = 13, k = 5, μ = 7 and σ a 5-cycle or a 4-cycle: none.
- **Goldberg's matrices and R.** Theorem A for v = 1225 = 35² (the only census value v ≡ 1 mod 8 with a known
  conference matrix of order v + 1) needs an R-structure on Goldberg's matrix of order 1226 built from a skew-Hadamard
  matrix of order 36. On the small members the R-structure fails: the matrix of order 10 (from skew H(4)) has one,
  the one of order 50 (from Paley's skew H(8)) has anti-automorphisms but no R with R² = −I (its whole signed
  (anti-)automorphism group, of order 14 112, was scanned), and the one of order 226 (from the doubled skew H(8)) has
  no anti-automorphism at all (nauty). For 1226 nauty did not finish in 25 minutes.

## 52. Propus-type quaternary pairs, Hermitian types and the chain sweep: round 3 (October 18)

Three exploration agents (search, character sums, Hermitian types) attacked the smallest open symmetric values and the
q ≡ 5 (mod 8) class through quaternary pairs; nothing was found. Scripts in the session scratchpad `round3/`. A
*bordered propus quadruple over ℤ_L* below is (p, A, A, q) with p, q symmetric, Σ PAF = −4 off 0 and rows (±2, 0, 0) or
(0, 0, ±2), equivalently a propus-type quaternary Legendre pair (A real, B symmetric) of length L (NOTES.md §32).

- **Exhaustive census for L ≤ 32** (meet in the middle of balanced A against symmetric (p, q), every candidate
  rechecked; L = 20 and 12 recounted independently). Normalised counts (p₀ = q₀ = A₀ = 1, unordered {p, q}): L = 4: 4;
  8: 64; 12: 48; 16: 0; 20: 960; 24: 2688; 28: 0; 32: 18432. A second, independent exhaustive search over all even
  L ≤ 32 (A up to rotation; checked against a brute force at L = 16, 22, 24) finds A-necklaces with a partner for
  L = 8, 10, 12, 14, 18, 20, 24, 26, 30, 32 and none for L = 16, 22, 28. So there is no such quadruple over ℤ₁₆ or ℤ₂₈, and
  hence no uniform algebraic family of length q − 1 for q ≡ 5 (mod 8): q = 17 and 29 already fail, and so does every
  cyclotomic or Sidelnikov decimation of index e of GF(Q) landing on L = 16 or 28 (e = 1, 3, 4, 5, 6, 7, 10, 12, …).
  The zeros come from the antipodal compression to ℤ_{L/2}, which is infeasible at L = 16 and 28 only (among L ≤ 148
  with small divisors). Solutions at L = 20, 24 are mostly unstructured, as is Verma's pair at L = 64.
- **Jedwab–Pender** (arXiv:2408.08472) Theorem 2(ii) (length (Q − 1)/2) and Theorem 3 (length 2p, Q = 2p − 1) are of
  propus type (checked through `propus-bordered` for Q ≤ 101 and p = 3, 7); they reach v = (Q + 1)/2 and v = 2p + 1,
  none of them a remaining census value.
- **Closed character-sum families** on GF(q)* = ℤ_{q−1} (Möbius-type sums λ(x)ψ(x + 1), class tables of order 2 and 4
  in (x, x + 1), A-first joins against all symmetric (p, q)): no quaternary Legendre pair at q = 13 (positive control
  only), 29, 37, 53, 61, 101, 109, 149. Symmetry forces B = const · ψ₄^{±1}(x)χ(x + 1), whose row sum has |Σ B|² of
  the order of q, incompatible with |Σ A|² + |Σ B|² = 2.
- **Cyclotomic orbit structures** on ℤ_{kq} (k ≤ 12, q ≤ 199) and on K × GF(Q)⁺ for v = 101, 163, 251, 325, 485, 487,
  491, 589, 973 (at most 18–22 orbits): nothing beyond v ≤ 21.
- **Hermitian, skew-Hermitian and α-twisted symmetrizations** (the types of §50). The exact monomials are worked out
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
  automorphism group, so solutions are plentiful and unstructured), and it recovers Verma's (p, q) from his A at L = 64;
  but its cost grows by about 1.8–1.95 per unit of L (≈ 800 CPU-s per hit at L = 52; Verma's pair was not re-found in
  1000 CPU-s at L = 64), which puts L = 88, 100, 106 and 118 (v = 89, 101, 107, 119) at 10⁴–10¹⁴ CPU-years.
  Simulated annealing at L = 100, 106, 118 and on unbordered propus over ℤ₁₀₁, ℤ₁₀₇ stays far from a solution.
- **Propus chain sweep** (§28, §40) resumed on the composite census values from 445: complete with no solution (every
  structure with 1–36 orbits) for 445, 473, 485 and 493; 455 and 475 stopped at the 50-minute limit per value; the
  job was then stopped by the session's time limit before 509.

**Round 4 (October 18): block-permuting groups and an automorphism lemma.** A dedicated engine (orbit classes of
(i, j, x) under transposition and h = (σ, μ); exact correlation equations; meet in the middle on the base rows,
forward-checking CSP on the σ-fixed rows, DFS on the moving cycles; scratchpad `round4/blk/`) was calibrated first:
from Verma's group alone (m = 7, k = 12, F = 2, μ = 2, σ = (6 7 8)(9 10 11)) it rediscovers his order-86 matrix (equal
to the negative of the published one up to signed permutation, and a rerun found one equivalent to the published
one itself; its other border pattern is exhaustively empty), and with σ swapping
blocks in pairs it finds Paley's matrices of orders 30, 42, 54, 62 and a non-Paley matrix of order 38 (which has an
R-structure). It agrees with an independent CP-SAT model on 64 small configurations.

- **Lemma** (proved, and checked on about 7000 automorphisms of matrices of orders 6–86). If π is a permutation
  automorphism of order r of a symmetric conference matrix of order n, then for every d | r the number of cycles of π
  of length divisible by d is even unless √(n − 1) ∈ ℚ(ζ_d) (the trace of C on the ζ-eigenspace is (a − b)√(n − 1) ∈
  ℚ(ζ)). For n = 118 and 154 an involution has f ≡ 2 (mod 4) fixed points and an element of order 3 has f ≡ 4
  (mod 6); for n = 118 order 5 needs f ≡ 8 (mod 10) and order 7 needs f ≡ 6 (mod 14). This excludes every ℤ_m-
  translation structure with 3 | m and one fixed point ((m, k, F) = (39, 3, 1), (9, 13, 1), (3, 39, 1), (117, 1, 1)
  for 118; (51, 3, 1), (9, 17, 1), … for 154), Verma-like layouts with σ of cycle type 1⁶3 or 3³ and μ of order 3, and
  most (σ, μ) pairs (201 of 226 for (118; 13, 9, 1), 23 of 36 for (118; 29, 4, 2)).
- **Exhaustive, no solution:** for n = 118, 35 of the 37 admissible groups with at most 48 free bits in the first
  row, among them the Verma-type group (13, 9, 1), μ = 3, σ = (3 4 5)(6 7 8) (none of its 31 380 canonical first rows
  has a compatible second row; CP-SAT agrees), all groups on (29, 4, 2), (59, 2, 0) and (58, 2, 2); open: (13, 9, 1)
  with μ = 5, σ of type 1·2·2·4 and with μ = 3, σ = id. For n = 154, 144 of 322 admissible groups are settled (no
  solution), 178 unresolved (the sweep was cut by machine load, not by a finding). So no conference matrix of order
  118 or 154 with this much symmetry (about 20 or more bits of reduction) exists in the searched families; Ho's and
  Wang's matrices, if they are of this kind, use weaker symmetry.

## 53. Large fixed sets for q ≡ 3 (mod 4), and bordered propus over the torus (October 18)

Two constructions aimed at classes that Theorems A and B do not reach. Neither works as first posed; the first one
survives with a non-Paley V (corrected the same day).

**An involution with 2q − 2 fixed points** (`sym_bigf.py`). In Theorem A the involution has 4 fixed points, the
field part is the U-part and the weighing matrix V of order 2(q − 1) is field-free. The roles can be reversed for
q ≡ 3 (mod 4). Let H be symmetric of order 4q with a permutation automorphism σ of order 2 that has f fixed points
and t two-cycles. In the basis (fixed points, first points of the pairs, second points),
H = [[F, B, B], [Bᵀ, P + V, P − V], [Bᵀ, P − V, P + V]]. Here V is a symmetric weighing matrix of order t and weight q,
and P is a symmetric signed matrix on the zero set of V. With t = q + 1 the zero set of V is an involution π, P is a
signed permutation matrix on it, and H² = 4qI becomes

    V² = qI,   BᵀB = (2q − 2) I,   F B = −2 B P,   F² + 2 B Bᵀ = 4q I.

Paley's skew core K and R: x ↦ −x, with −1 at ∞, give V = KR, which is symmetric with π fixing 0 and ∞. The Q-part
(F, B, P) is field-free: it depends on q only through its size. Changing the sign of P on a pair is the same as
changing the sign of one column of B, so up to equivalence only the signs at the fixed points of π matter. A search
with B fixed must still enumerate all of them.
- q = 3 (order 12): Q-parts exist for π with 2 fixed points (24 of 48 signed patterns) and with 4 (8 of 16). Only
  the first complete: a symmetric W(4, 3) cannot have zero diagonal.
- q = 7 (order 28): B is 12 × 8 with orthogonal columns. Such a B is equivalent to 8 columns of Paley's H(12). An
  independent check (scratchpad `verify53/extension/`) classified the k × 12 partial Hadamard matrices exhaustively,
  confirmed by a mass formula against a clique count. The class counts for k = 1, …, 12 are 1, 1, 1, 2, 1, 2, 1, 1,
  1, 1, 1, 1, and every class with k ≥ 5 completes; this is Verheiden's completion of the last seven rows. The check
  also enumerated the 190080 automorphisms of H(12), whose 95040 column permutations act sharply 5-transitively.
  With B fixed, each row of F lies in the affine set {x ∈ {±1}¹²: Bᵀx = −2(BP)ᵢ}. Exhaustive over every involution π
  on the 8 columns and every sign pattern of P (sign changes of P are column sign changes of B, so with B fixed none
  can be dropped):
  - Q-parts exist only for π with 0 fixed points (24 of 1680 signed patterns) or 4 fixed points (24 of 13440). There
    are none for π with 2, 6 or 8 fixed points, so Paley's V never fits.
  - A symmetric W(8, 7) with zero set π exists exactly when π has 0 or 2 fixed points (exhaustive DFS and CP-SAT).
  - Hence the completions are exactly those with π fixed-point-free: symmetric Hadamard matrices of order 28 with an
    involutory permutation automorphism with 12 fixed points exist. An independent check found them by four routes,
    including two representatives in Valtonen's classification data (arXiv:2608.06920, Zenodo
    10.5281/zenodo.22671922), and reached at least 8 of the 91 symmetric classes. It also found order 20 with 8 fixed
    points (q = 5).
  - The first version of this section said no such matrix of order 28 exists. That version fixed the pair signs of P
    to +1, which is invalid when B is fixed, and so searched only 1/16 of the patterns. `sym_bigf.py` now enumerates
    every sign pattern and reproduces the counts above.
  - So the reversed decoupling fails with V = KR (π fixing 0 and ∞), but not with every Paley-type V. For every
    prime power q ≡ 3 (mod 4), V_xy = −χ(1 + xy) (x, y ∈ GF(q)), V_x∞ = χ(x), V_∞∞ = −1 is a symmetric W(q + 1, q)
    whose zero set is the fixed-point-free involution x ↦ −1/x (0 ↔ ∞). This is Paley's core twisted by an element of
    PSL(2, q) with the signs χ(x), which turn it into an anti-automorphism. Found by a round-7 agent; checked here for
    q = 3, 7, 11, 19, 23, 31, 43 (an earlier version of this paragraph wrongly argued from unsigned permutations that
    Paley's core has no such V). So for q ≡ 3 (mod 4) only the field-free Q-part with π fixed-point-free remains: at
    q = 7 it exists, F = K(−2P ⊕ 2N)Kᵀ/12 with K = H(12) and N a symmetric integer 4 × 4 matrix with N² = 7I. Later
    sections report on larger q.
- The other admissible shapes for q ≡ 3 (mod 4) need t ≡ 0 (mod 4), q + 1 ≤ t ≤ 2q − 2, f = 4q − 2t. A symmetric
  weighing matrix of order t ≡ 2 (mod 4) and weight q needs q to be a sum of two squares. The case f = 0 would be two
  symmetric W(2q, q) with complementary supports, which this rules out for q ≡ 3 (mod 4) prime. Paley gives V only
  for t = q + 1.
- With an order-3 automorphism in addition (the squares of GF(7)* acting on the pairs and every orbit type on the
  fixed points, with π = Paley's and free signs on the pair orbits), CP-SAT finds no solution for q = 7, consistent
  with the above (scratchpad `bigf/qinv2.py`).

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

## 54. Skew analogues of the decoupling (October 18, exploration round 4)

Theorems A and B split a symmetric matrix along a pair-swapping involution. Three agents looked for skew analogues:
amicable pairs, automorphisms of order 3, and a direct skew Yamada–Kiyasu. Each claim was checked by two or three
independent refuters. No census value is settled. Scripts are in the session scratchpad (`skewwf/`).

**No involutions.** A skew-type Hadamard matrix H = I + S admits no signed permutation P with P² = I and
PHPᵀ = H other than ±I. A two-cycle (a, b) with sign product +1 would give H_ba = H_ab, and a diagonal P with all
signs equal is ±I. So Aut(H) contains a unique involution, −I, and its Sylow 2-subgroups are cyclic or generalised
quaternion. Examples: Q₁₆ for Paley's H(8), Q₈ for the Goethals–Seidel matrices over ℤ₇, and only {±I} for some
Goethals–Seidel skew H(36) over ℤ₉ (a refuter's counterexample to a claimed Q₈ in every such array). The pair
decoupling of Theorems A, B and of Miyamoto's matrices can never be skew. An automorphism of odd order d cannot
decouple supports either: on every cycle of length l > 1 a ±1 row has an odd number of entries, so its cycle sums
never vanish. For an automorphism of odd prime order p with t cycles of length p, t is even unless p ≡ 3 (mod 4) and
4v − 1 = pk² (character trace argument; checked by CP-SAT for p = 3, 5, 7 at orders up to 36).

**Amicable pairs.** Skew H(4v) = [[I + S, A], [−A, I + S]] holds iff A is symmetric ±1, S is skew ±1 with zero
diagonal, AS = −SA and A² + SSᵀ = (4v − 1)I (Wallis–Turyn; classical).
- *Circulant normal form.* For 2 × 2 blocks of circulants over an abelian group of odd order v, every pair has the
  form A = [[a, b], [b*, −a]], S = [[s, t], [−t*, −s]] with aa* + bb* + ss* + tt* = 4v − 1 and 2as = bt* − b*t. The
  traceless form is forced: at each character write the blocks in the Pauli basis. The scalar case would put
  √(4v − 1) or i√(4v − 1) into ℚ(ζ_e), e | v. That is impossible because the odd primes of the square-free part of
  4v − 1 ramify there; a refuter supplied this step.
- *Counts over ℤ_v.* v = 3, 5, 7, 9, 11, 13: 48, 160, 112, 0, 880 and 832 solutions; none over ℤ₃². All counts are
  divisible by 4v, and every solution (2032 in all) builds a skew H(4v), checked exactly.
- *Paley plus a diagonal.* A = C + Δ with C the Paley conference matrix of order q + 1 works only for q = 5 and 9
  among q ≤ 17. The search is exhaustive over the balanced Δ for q = 13 and 17, and Δ must be balanced since tr A = 0.
  At q = 9 the admissible Δ form one orbit of 72: the index-2 cosets of the non-split tori ℤ₁₀. In torus
  coordinates the ansatz is empty for q = 13, 17, 25, 29, 37 and 41. The trivial character needs t₀² = q − 2m₀,
  which fails for q = 73, 89, 97, …
- *Twisted pairs* H = LHR, L ≠ R signed involutions. They decouple H into two eigenspace blocks with MMᵀ = VVᵀ = NI,
  with entries 0, ±1, ±√2, ±2 where fixed points occur. Corrected by a refuter:
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
  matrices. Complete searches: `skeweq.c`, and CP-SAT in independent forms. The control, Paley II of order 12, is
  equivalent.

**Order 3** (the Z₃-gluing theorem). Let H be skew-Hadamard of order 4v with a permutation automorphism σ of order 3,
f fixed points and t three-cycles, N = f + 3t. Write Z′_ab for the sum of the circulant block (a, b) and
V_ab = ½ Σ_c h_ab(c) ω^c. Then:
- either |Z′_ab| = 3 and V_ab = 0, or Z′_ab = ±1 and V_ab = −Z′_ab ω^{k_ab}, with k_ab the position of the minority
  entry;
- Q = [[F, √3B], [−√3Bᵀ, Z′]] has Q + Qᵀ = 2I and QQᵀ = 4vI; V has V + V* = I and VV* = vI, so V² = V − vI;
- k_ba = −k_ab for a ≠ b, and conversely F, B ±1 with these relations rebuild H.
- Each row of V has exactly v units, so f ≤ v and f ≡ v (mod 3).
- If t is odd then 4v − 1 = 3s². This is the centred-hexagonal gate of §35–36, now for every order-3 automorphism.
  Otherwise t is even and the diagonal of V is balanced.
- For f = 4, F is a skew H(4) and BBᵀ = tI.

So the invariant part fixes the signs of the ω-part. Unlike Theorem A, no input can sit freely in V. Further
results:
- Theorem A's block shape (three full blocks and a diagonal one per block row) forces the constant pairs into four
  cliques K_m.
- The borders cannot be constant on blocks, by a column-sum budget.
- Lemma D: V cannot have the form V_kl = D_kl U_l with diagonal D_kl and one free unitary U_l per block column,
  once m ≥ 3.
- If the class counts R vanish, or the three off-diagonal blocks have rank-one sign patterns, then m ≤ 4.
- Skew matrices of orders 28 and 52 with the four-clique shape exist (CP-SAT plus a SAT lift of the phases; a refuter
  built its own). The shape is specific to Z₃; it is not the support of Theorem A's V, as first claimed.
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
- *Z[C₃] array.* The pattern with the core in the constant part of one block has a solution at v = 7 (order 28),
  none at v = 13 (three independent proofs, among them a mod-2 system UNSAT in 0.1 s), and is open at v = 19, 25,
  where the mod-2 system is satisfiable.
- *Quaternary route.* A quaternary Hadamard K of order 2v gives a skew H(4v) by Turyn's map when a monomial W makes
  M = KW twisted-Hermitian, M* = iM + (1 − i)I, and only then. M then has eigenvalues
  ((1 + i) ± (1 − i)√(4v − 1))/2, each of multiplicity v, and det M = (2iv)^v. Hence a necessary invariant:
  arg det K ∈ (π/2)ℤ. Zaitsev's QH322 and QH330 are already in this form up to column phases. Their orders, 644 and
  660, were already settled; a refuter rebuilt both skew matrices. QH214, QH238 and our C214, C238 have none
  (complete DFS). The invariant excludes QH490, QH610, QH890, QH1066 and QH1690 (v = 245, 305, 445, 533, 845).
  QH382 and QH470 timed out, and 18 more at open v were not run: round 6 (below) scans them.

**Lead.** 106 open skew census values are v = (q + 1)/2 with q ≡ 1 (mod 8) a prime power, i.e. skew H(2(q + 1)),
the classical gap between Szekeres (q ≡ 5 mod 8) and Whiteman (q = p^t, p ≡ 5 mod 8, t ≡ 2 mod 4). Round 5 attacks
it.
