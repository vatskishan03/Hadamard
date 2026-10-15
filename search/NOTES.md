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
  none for n = 9 (35 is not a² + 2b²).
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
Hadamard matrix of order 124 · 123 = 15252 = 4 · 3813 (built and checked by `sym_mm1.py build`). Paley II (q = 7625),
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

## 32. Complex orders 342, 426 and 1558: the low-energy k-tree and the end of the exact sweep (October 14)

**Two fixes.** `sigma_file` in `lowe_sz.py` built the law of a solution's split value from kt4's lists 0 and 2
(A + B), but kt4 pairs its lists 0 and 1 (A + C). Its split values were off the lattice of a + c, so kt4 formed no
pairs at all. With the fix, the N = 131 control finds 416 of the 2560 solutions inside its lists in the first 3000
split values, against 400 in kt4's own order. In `kt4.c` the small-join branch (fewer than 2¹⁸ pairs per split value)
never freed its hash arrays: about 9 MB/s at N = 213, past 10 GB in one run.

**Order 342.** N = 171 = 9 · 19, H = ⟨11⟩ of order 6 (33 orbits; ⟨H, −1⟩ has 17 classes), rows (1, 11; 11, 21).
- The exact join is out of reach: A and B have 7128 and 6864 members, C and D 1.5 · 10⁸ and 2.6 · 10⁸.
- The filter kept 3 · 10⁷ members of C and of D, retaining 49% and 38% of a solution's members under the model.
  Estimate inside the lists: 10^0.94 (structure 10^1.67).
- kt4 found 12 solutions at the 814th split value (62 s).

**Order 426.** N = 213 = 3 · 71, H = ⟨37⟩ of order 7 (33 orbits), rows (13, 17; 13, 15). A and B have 6435 members;
C and D share one pool of 4.36 · 10⁸ (one enumeration), of which 3 · 10⁷ were kept. kt4 found 72 solutions at the
40th split value; the whole run took 274 s. Estimate inside the lists: 10^1.1. The rows (15, 17; 13, 13) of the same
structure gave nothing in 11,000 split values.

**Order 1558.** N = 779 = 19 · 41, H = ⟨49⟩ of order 60 (order 3 modulo 19, order 20 modulo 41), rows (5, 7; 39, 39).
The exact sweep's lists of 120, 120, 18480 and 18480 hold 96 solutions (1.8 s). The per-class Gaussian estimate of
this structure is 10^−13.6. Like 742 (⟨13⟩ of order 26 on 7 · 53) and 754 (⟨9⟩ of order 42 on 13 · 29), it is a large
diagonal cyclotomic group on a product of two primes, with tiny lists and many solutions.

**The exact sweep is complete** for the open N < 1000 (pair work up to 2 · 10¹⁰ per structure and pattern): the two
halves of the range 583–995 met at 833. 779 was its only hit in that range.

## 33. How good is the Gaussian estimate? (October 14)

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
  10^3.75, 10^3.98, 0, 10^3.45, 10^4.33, 0, 10^3.11, 0; hybrid 10^3.65–10^4.36 for all eight; per-class 10^2.06–10^2.57.
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
  members; C and D share one pool of 1.45 · 10⁹ (one enumeration). kt4 found 200 solutions at the 2967th split value
  (1107 s in all).
- **Order 454**, the first order from the bordered array ([BORDERED.md](../complex-hadamard/BORDERED.md)). N = 227 =
  2 · 113 + 1 is a safe prime; over ℤ₂₂₆ = ℤ₂ × ℤ₁₁₃, H = ⟨49⟩ of order 7, 34 orbits, rows (0, 0, 2, 0). Per-class
  estimate 10^0.63, hybrid 10^4.4. C and D: pools of 6.0 · 10⁸ and 1.2 · 10⁹, 3 · 10⁷ and 2.77 · 10⁷ kept. The first
  run's parent process was killed for memory at the 850th split value; the resumed join found 512 solutions at the
  851st.

**Counter-examples.**
- N = 183 = 3 · 61, H = ⟨34⟩ of order 5, rows (17, 19; 1, 9): hybrid 10^6.69, per-class 10^4.98 (10^2.04 inside lists of
  2 · 10⁷ of 2.7 · 10¹⁰). No solution after 6833 split values, which cover 17% of the pair-product mass
  Σ_σ N_AC(σ) N_BD(t − σ). The hybrid estimate is too high here (it would put 10^4.4 solutions in the lists).
- N = 183, H = ⟨20⟩ of order 10 (it contains the element that is −1 on ℤ₃ and 1 on ℤ₆₁): exhaustively searched,
  no solution in any of its 14 patterns, against 10^0.4–10^1.0 each (hybrid) and 10^−0.3–10^−0.6 (per-class).

## 34. Retention of the low-energy filter, measured on exact solution sets (October 14)

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
mean, and a solution's A + B is one of a few discrete values. The runs that found solutions kept 2–20% of their
pools; the 183 run keeps 0.07%. The two energies select different members, so they are separate lottery tickets.

**Split-value order.** For the 183 run above, the 6833 split values processed (kt4's order by the Gaussian law of a
solution's split) cover 17.1% of the pair-product mass; the best order for that cost, by N_AC N_BD/(N_AC + N_BD),
covers 19.0%. The order is close to optimal.

**Memory.** The `lowe_sz.py` parent kept the full kept value arrays (about 1 GB per list) while kt4 ran, and the
out-of-memory killer took two parents during this work (their kt4 children ran on). The parent now keeps only the
masks and 2 · 10⁵-row subsamples.

**Bordered sweep.** `bsz_sweep.py` had been run only for N = 107 and 155. It is now running over the other 210 open N.
`bsz_group.search` handles every pattern of a spin structure, so `plan()` now makes one job per spin structure
instead of four.

## 35. Iterated maps, and the matrices of August 2026 (negative, October 14)

**RRR (Elser's difference map)** for four ±1 sequences with a flat summed power spectrum (`search/rrr.c`,
`search/rrr.py`). Constraint A: ±1 vectors with fixed row sums (and symmetry), projected by sorting. Constraint B: per
frequency, Σ_k |X̂_k(ω)|² equals the target, projected by rescaling. It solves plain Goethals–Seidel quadruples at
m = 19 in about 2 · 10⁵ iterations and at m = 31 in about 3 · 10⁷ (1 of 8 runs of 4 · 10⁶; β = 0.5). The cost grows
by about 1.5 per unit of m, like local search (§26), so it is no route to unstructured lengths above about 40.

**The matrices of Alpöge et al. (August 2026)**, read from the public mirror github.com/bbeartheancient/hoa64
(commit d825201). All twelve are stored normalised (first row and column positive).
- Undoing the normalisation (in a circulant block, x(j) = c · H′[1, j] · x(j − 1)) recovers the Goethals–Seidel
  quadruples of 892, 1132, 1244, 1948 and 1964 (four sequences on ℤ_N with Σ PAF = 0). Under a border of 4 rows the same
  recovers the Spence seeds of 668, 716, 1676 and 1772 (four sequences on ℤ_{N−1}, Σ PAF = −4, row sums (2, 0, 0, 0));
  for even m each block is fixed only up to the character (−1)^t, and the twists were chosen to make Σ PAF constant.
- No block of any of them is symmetric or skew about any centre, so none gives a Szöllősi or bordered Szöllősi
  quadruple.
- The Goethals–Seidel blocks are invariant up to translation under multiplier groups of orders 3, 3, 5, 9 and 7, and
  the Spence seeds of 1676 and 1772 under groups of orders 15 and 12.
- The seeds of 668 and 716 (lengths 166 = 2 · 83 and 178 = 2 · 89) have no affine symmetry, no relation between blocks,
  no relation between their CRT slices, no slice correlated with a Legendre sequence beyond random size (largest
  |correlation| 30 of 83 and 89), and per-block spectra only as flat as the sum condition forces. Their ℤ₂ splits
  u = (a₀ + a₁)/2, v = (a₀ − a₁)/2 satisfy Σ PAF_u = −2 and Σ PAF_v = 0 on ℤ_q with weights 2q − 1 and 2q + 1, as every
  solution must. They are unstructured solutions of length 166 and 178, far beyond local search or RRR.
- 1388, 1436 and 1916 fit neither layout.

**Safe primes.** The same split applies to the bordered Szöllősi array over ℤ₂ × ℤ_q: it is exact and the parity of
the v-stage is automatic, but with A and B symmetric the v-stage at q = 53 (order 214) has about 82 free signs for 95
bits of conditions, so a two-stage search needs some 10⁴ u-solutions per success, each a length-53 four-sequence
problem. Not feasible with these tools.

## 36. Dead split values, and exact compression checks before a search (October 14)

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

`lowe_sz.py` now samples 2 · 10⁷ random (A, C) pairs and 2 · 10⁷ random (B, D) pairs, keeps the merge-class values
that occur on both sides (`live_splits`) and orders only those by the Gaussian law. At N = 215 there are 20000 such
values, and the one holding the 200 solutions of order 430 ranks 6th (2967th in the old order). At N = 203 none
occurs, so the overlap Σ_σ p_AC(σ) p_BD(σ) of the two sides' distributions is below about 3/M² ≈ 10^−14 (M = 2 · 10⁷),
and the join is skipped. The 341 and 203 runs above therefore tested nothing; 341 is queued again on its saved
lists (`--resume`).

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
  solutions fails. Of the 264 without, levels with up to 2 · 10⁵ states per block rule out 186 of 238 cyclic and 24
  of 26 bordered ones, and levels with up to 3 · 10⁶ states (seconds, at most 1.4 GB) 220 of 238 and 24 of 26.

For the bordered array (v = N − 1 even) the level d = 2 has a closed form; `bsz.z2_ok` already applied it in the exhaustive sweep (`bsz_sweep.py`), but `lowe_sz.py`, `frontier_scan.py` and `hybrid_rank.py` did not, so impossible bordered structures were ranked and queued. Σ PAF = −4 gives
Σ_k |X̂_k(ω)|² = 4v + 4 = 4N at every ω ≠ 0. At the character (−1)^t this reads Σ_k e_k² = N, where 2e_k is block
k's alternating sum. The multiplier group fixes the e_k modulo a small number.

At N = 165 (v = 4 · 41, H = ⟨37⟩: trivial modulo 4, order 5 modulo 41) the orbits force:
- e_A, e_B ≡ ±2 (mod 20) for the symmetric blocks of row sum 0;
- e_C ≡ ±1 (mod 10);
- e_D ≡ 0, ±2 (mod 10).

The sum of the squares then lies in {9, 13, 29, 33, 49, 53, 69, 73} modulo 80, and 165 ≡ 5 (mod 80) is not in that set.

Run on the queue of hybrid-ranked candidates (§33), the check removed 15 of 24 tickets:
- ruled out at d = 2 (bordered): 165 ⟨37⟩ with the 2 in A or C; 269 ⟨9⟩ with the 2 in A, B or C; 223 ⟨7⟩ with the 2
  in A or C. Three of them were already running and were stopped.
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

`lowe_sz.py` now runs the check for d ≤ 64 (up to 3 · 10⁶ states) before enumerating and stops on an empty level.

## 37. Fibre filters: the compressions inside the low-energy search (October 14)

A solution's four compressions to ℤ_d form a level-1 solution (§36), so each block's compression lies in the
projection of the level-1 join. `lowe_sz.py` now lists these sets for every prime power d ≤ 32 dividing the group
order whose states fit (`fibre_filters`). `lowe_enum.c` drops every member whose compression is not among them,
before computing its energy: a 64-bit key of the compression, linear in the orbit signs, is updated with each flip
and looked up at the leaves. The Gaussian statistics, the energy quantile and the retention model are taken over the
allowed sub-pool (rejection sampling).

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

**kt4 on a shared machine.** For each split value kt4 opens two or three OpenMP parallel regions over about 10⁵
pairs. With several joins and sweeps on 4 cores (load 11) their barriers wait for descheduled threads. On the 453
lists a single thread processes 128 split values per second, two threads 57, and the default four 4. The earlier
runs of this section and of §32–§34 thus covered about 30 times fewer split values than their running time allowed.
`lowe_sz.py` now runs kt4 with one thread (`--kt4-threads`).

**Where a solution's split value lies.** At N = 203 (⟨36⟩, rows 3, 25, 3, 13) the kept lists had no live split value
at all (§36). The filter's centring explains it (simulation on uniform samples, `lowe_sz` filter at 3.6% kept).
- On two classes the target lies 2.66 standard deviations (of the four-block sum) away from Σ μ.
- The filter centres each big block at z_k = t − Σ_{j≠k} μ_j, assigning the whole offset to that block, and it does so
  for C and for D at once. On those classes the kept C and D lie 5.5 of their own standard deviations from where a
  solution's C and D lie.
- The kept A + C and t − B − D sums are 2 pooled standard deviations apart, so random pairs never meet.
- Centring at the solution law instead (c_k = μ_k + Σ_k Σ⁻¹(t − Σ μ), consistent by construction) does not help: the
  separation becomes 2.6. On these classes A and B carry most of the variance, so a solution's split a + c lies in
  the tails of both sides, where sampling never looks.

The split values are therefore now also taken from the lattice. The values of A + C on the merge classes lie in a
coset of Λ_A + Λ_C (a Hermite basis of the differences within each list; for 371 it has rank 4 on 5 merge classes).
The points of that coset in a box around the solution law's mean that also lie in the coset of t − B − D are added
to the sampled live values, and all are ordered by the law. kt4 discards a split value without pairs in about a
millisecond.

**Other centres for the filters** (`retention_test.py --center`).
- Two alternatives to z_k were tried: splitting the offset t − Σ μ between C and D by their covariances, and the
  solution law's mean c_k. Neither helped on the exact solution sets.
- 163 (prime): 0 of 8352 at 5% kept for all three centres, about 850 at 21%.
- 93: 1224 of 21384 at 15000 with z_k, 0 with the other two.
- 371: 24 and 48 of 96 at 500 and 1000 with z_k, 0 with the others.
- 267: 0 for all three.

z_k ranks members by the likelihood ratio, which is optimal for each block alone. The joint losses come from the
few discrete values of a solution's a + b, not from the centre.

## 38. Order 366 (October 14)

N = 183 = 3 · 61, bordered array over ℤ₁₈₂ ([complex-hadamard/BORDERED.md](../complex-hadamard/BORDERED.md)),
H = ⟨3⟩ of order 6 (trivial on ℤ₂, order 6 on ℤ₇, order 3 on ℤ₁₃), 36 orbits, rows (0, 0, 2, 0).
- **Check.** It passes the exact compression check at d = 2, 7 and 13 (16, 8 and 60216 level-1 solutions).
- **Estimates.** The hybrid estimate is 10^3.18 solutions; with fibre filters the model's retention at L = 3 · 10⁷ is
  10^−0.13 (`fibre_rank`).
- **Lists.** The fibres at d = 2, 7 and 13 allow 6.9% of A and B (3880 of 55968), 2.1% of C (2.68 · 10⁷ of
  1.28 · 10⁹, all kept) and 3.6% of D (9.35 · 10⁷ of 2.57 · 10⁹, of which 3 · 10⁷ were kept).
- **Join.** kt4 (one thread) took 32423 split values in the law's order: 12599 sampled live values and the rest
  lattice points. It found 1536 solutions at the 146th, about 440 s into the join, after 1.7 · 10⁹ (A, C) and
  7.6 · 10⁹ (B, D) pairs. That split value was among the sampled live values, so the lattice points did not matter
  here.
- **Earlier runs.** Without fibre filters the same structure was run three times: full energy, a resume of it, and
  unit-class energy. The lists kept 2.3% of C and 1.2% of D, and the runs processed 642, 694 and 3358 split values
  (four threads, in part dead split values) without a solution.
- **Certificate.** C366.json (construction `szollosi-bordered`, border x = s = t = y = 1, u = v = −i), checked by
  `verify/verify_complex.py`. Order 366 is blank in Seberry–Yamada's Table 11.2 (label **New**, NOVELTY.md).

The hybrid estimate of 10^3.2 for the structure and retention near 1 for the filtered lists put about 10^3
solutions in the lists; 1536 at a single split value is consistent with that. For 183 ⟨87⟩, a different subgroup of
order 6 with the same estimate, the run is in progress.

**All levels at once.** Each level's fibre filter is a necessary condition, so a solution's member of each block
must pass all of them together. When the allowed sub-pool is too small for rejection sampling, `lowe_sz.py` now lists
it exactly (`lowe_enum` with the filters and no energy cut). An empty list rules the structure out, even though every
single level passed.
- Bordered 273 (ℤ₂₇₂ = ℤ₁₆ × ℤ₁₇, multiplier groups of order 8), ⟨155⟩ with rows (2, 0, 0, 0): the levels 2, 4, 8 and
  17 each allow 19%, 10%, 6.4% and 9.4% of A, but no member of A passes all four. The per-level rates sampled from
  the pool match the masses from the compression states (0.185 against 0.187 at d = 2).
- The same holds for ⟨155⟩ with rows (0, 0, 2, 0) and ⟨53⟩ with rows (2, 0, 0, 0).
- ⟨53⟩ with rows (0, 0, 2, 0) survives with 2.4% of A and B and 0.17% and 0.25% of C and D (3.2 · 10⁷ and 4.7 · 10⁷
  members): lists of 5 · 10⁷ would be complete.

## 39. Re-ranking the open orders with the checks and the fibres (October 14)

**Screen.** The structures of the frontier map (`frontier_scan.py`, open N up to 581) with a per-class estimate of at
least 10⁻⁶ and pools of at most 10¹¹ make 1416 (structure, pattern) records. Level-1 joins with up to 2 · 10⁵ states
rule out 853 of them, among them all 40 bordered records. Levels with 3 · 10⁶ states rule out 73 of the remaining
563. Bordered records with per-class estimates between 10⁻¹² and 10⁻⁶ add 132, of which 108 fall.

**Settled by complete lists.** When the fibres leave fewer members than L in every block, the lists hold every
allowed member, so an empty join settles the structure.
- Bordered 287 = 7 · 41, ℤ₂₈₆ with H = ⟨3⟩ of order 15. With the 2 in C the lists are 424, 424, 31520 and 191712;
  with the 2 in A they are 180, 324, 577472 and 577472. An exact join (`bnw.join4`) of both gives 0 solutions. The
  patterns with the 2 in B or D are images of these.
- Bordered 273 over ℤ₂₇₂, multiplier groups of order 8. ⟨155⟩ and ⟨19⟩ in both patterns, and ⟨53⟩ and ⟨189⟩ with
  the 2 in A, have no member of A that passes all levels (§37). ⟨53⟩ and ⟨189⟩ with the 2 in C survive, with complete
  lists of about 3 · 10⁷ and 5 · 10⁷ for C and D, but their hybrid estimates are 10^−1.8 and 10^−1.9.

**203 ⟨36⟩.** All three patterns ran with fibres (d = 29 only, which allows 74–100% of each pool) and lattice split
values. No random pair met on either side, and the 190512–193167 lattice points within Mahalanobis radius 3.4–4.0
of the solution law gave no solution. With 2–4% of C and D kept, the hybrid estimate of 10^3.4 is not borne out.
That is either a failure of the estimate on these classes or the joint retention of §34; the level that would
decide, d = 7, has too many states to list.

**341 = 11 · 31.** Eight multiplier groups carry per-class and hybrid estimates of 10^1.8–10^3.0, with pools of
2 · 10⁹–9 · 10⁹:
- ⟨2⟩, ⟨35⟩, ⟨95⟩ and ⟨194⟩, with fibres at d = 11;
- ⟨15⟩, ⟨27⟩, ⟨58⟩ and ⟨170⟩, with fibres at d = 11 and 31.

With fibres the expected number of solutions in lists of 3 · 10⁷ is 10^0.2–10^1.1, highest for rows (1, 21, 9, 29).
The seven best are queued.

**Pools above 10¹¹.** 677 records with pools of 10¹¹–10¹³ (per-class estimate at least 10⁻⁴, or 10⁻¹⁰ for the
bordered array); 199 fall at the cheap levels. For the survivors the exact per-level fibre masses give a cheap
allowed fraction; the hybrid estimate was then computed for the best:
- **Bordered 203 = 7 · 29 over ℤ₂₀₂ = ℤ₂ × ℤ₁₀₁, H = ⟨87⟩ of order 5, rows (2, 0, 0, 0): hybrid 10^6.47.** Only d = 2
  gives a fibre (41% of C and D, 34% of A and B), so C and D keep 3 · 10⁷ of 1.1 · 10¹¹ allowed members each. With
  the 2 in C the structure fails at d = 2. Running (2.8 · 10¹¹ leaves per big block).
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

**Runs of this round** (fibres, live and lattice split values, one kt4 thread):
- 403 ⟨16⟩, rows (9, 9, 15, 35). The fibres at d = 13 leave 4.7% of C and 5.6% of D, so C and D are complete
  (2.04 · 10⁷ and 2.31 · 10⁷). 1.3 · 10⁶ split values (all 575692 sampled live values and lattice points) give no
  solution; the hybrid estimate was 10^0.72.
- 219 ⟨10⟩, rows (11, 25, 7, 9) and (9, 11, 7, 25). The fibres at d = 3 leave C and D complete (1.2 · 10⁷ and
  6 · 10⁶). The first 3 · 10⁵ of about 1.9 · 10⁶ live split values in the law's order give no solution; the hybrid
  estimate was 10^−1.
- Bordered 161 (ℤ₁₆₀): all four multiplier groups of order 4 (per-class estimates 10^−0.7–10^−1.5) fail at d = 2 in
  every pattern, so order 322 stays out of reach. The map has no structure for N = 179 in its orbit range.

## 40. Resumable runs; 453 settled to 91% (October 15)

**Restart.** The container was reclaimed about five minutes after the session went idle, and every background job
died: the 341 joins at their 8164th and 8514th split values, and bordered 203 ⟨87⟩ in the middle of listing block C.
`lowe_sz.py --budget SECONDS --auto` now makes a run resumable:
- Big blocks (pools above 10⁸) are listed in 64 chunks of lowe_enum's parallel prefixes (`lowe_enum ... LO HI`),
  each kept under `WORK/lowe/chunks` with a fingerprint of its input; the lowest-energy members are merged at the end.
- The lists are saved by atomic renames, and the join records how many split values it finished, and on which merge
  classes, in `kt4_*.bin.state.json`.
- A call stops by itself before its budget ends (exit code 75), and the same command continues: from the first
  unfinished chunk, or from the first unjoined split value.

The test was bordered 183 ⟨3⟩ (order 366), from nothing, in steps of 150 s.
- The chunked lists are identical to the unchunked run's: D passed 30300152 members at τ = 6.607, as before.
- The join ran over six interruptions and found the 1536 solutions at split value 146.
- With the orbit order and prefix checks of §39, C took 2.7 · 10⁷ leaves (12 s) instead of 1.28 · 10⁹, and D
  9.5 · 10⁷ instead of 2.57 · 10⁹.

**453 ⟨4⟩, rows (29, 31, 1, 3).** The join on the fibre-complete lists (all 2520 allowed members of A and of B,
4.8 · 10⁷ of C and 1.6 · 10⁷ of D) finished all 4329956 sampled live split values, 6.95 · 10¹⁰ (A, C) and
2.31 · 10¹⁰ (B, D) pairs, without a solution.
- Coverage was measured with fresh independent samples (2 · 10⁷ pairs per side, `coverage.py` in the scratchpad). The
  processed split values hold 62% of the split values met on both sides and 90.6% of the collision mass Σ_σ n_AC n_BD.
  The heaviest missed split values carry 10–16 sample pairs per side.
- So at most about one solution in ten could have escaped. The hybrid estimate of 10^1.8 (about 60) is contradicted
  by at least 1.5 orders of magnitude.
- Either the estimate fails here, as it did for 203 ⟨36⟩, or the compression to ℤ₁₅₁ (too many states to list)
  rules the structure out.
