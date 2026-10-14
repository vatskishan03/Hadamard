# Notes on structured searches (October 2026)

What we learned while searching for the unresolved orders in [FRONTIER.md](FRONTIER.md), including the
searches that found nothing. Section 1 is a mechanism, sections 2–4 and 6 are negative results, section 5 is
the literature as of October 2026.

## 1. Why some structures have many solutions

Take a group G = M × ℤ_p, p prime, and a symmetry group S of the four sequences. Let C = S ∩ ℤ_p* act on
the ℤ_p coordinate, f = (p − 1)/|C|, and t = |π_p(S)|/|C|, where π_p is the projection to ℤ_p*. On every
fibre {a} × ℤ_p a sequence is constant on the cosets of C, so it is x₀(a) on 0 and v₀(a), …, v_{f−1}(a) on
the f cyclotomic classes. Expanding the class vector in the characters of ℤ_p*/C ≅ ℤ_f, every nonzero
mode pairs with a Gauss sum of absolute value √p, and the Goethals–Seidel condition splits by the
difference of the modes.

- **f = 2 (a Paley product).** Every fibre is either constant ("w-type") or ±(Legendre sequence)
  ("z-type"). The identity coefficient forces (p + 1)N_w − 4N₋ = 4|M| for the numbers of w-type fibres,
  and the equations split by powers of p. This is what makes the skew Paley-product orders on the skew
  branch easy.
- **f = 6 (order 3964).** For our matrix of order 3964, M = ℤ₂² × ℤ₁₃, p = 19, C = ⟨11⟩ of order 3,
  so f = 6 and t = 3. Its group of order 27 also contains (α², 3, 1), which pairs the order-3
  automorphism α of ℤ₂² with an order-3 multiplier of ℤ₁₃. Of its 208 fibres, 75 are ±Legendre (36%,
  against 2 of 64 patterns, or 3%, for a random fibre), 7 are constant and 126 are other cubic-class
  patterns.
- **f = 8 (order 1916).** Alpöge et al.'s matrix of order 1916 has p = 17, C = {±1}, f = 8 and t = 2.

With f = 2 the counting identity can only be met by a few fibres of energy ≈ p, which the border's flat
coset sums forbid (section 2). With f = 6 or 8 the fibre sums move in steps of 2, so many small terms
can balance it.

As a rough model, the Gaussian estimate of the solution count misses a factor of about √p for each
class whose character is nontrivial on ℤ_p. On the skew branch's exact count for v = 253 this brings
the prediction to within 10^0.8–10^1.8 of the truth; the uncorrected model is off by 10^6–10^9.

The shape is necessary but not sufficient. 9,677 of the structures in section 3 have it (p ≥ 11,
4 ≤ f ≤ 12, t ≥ 2), and 108 of those are exactly empty.

## 2. Paley products do not reach the unresolved orders

- **Goethals–Seidel over ℤ_n.** A Paley product needs n = m·q with q ≲ 4m + 3. Every composite
  unresolved n ≤ 2999 has its prime factors split as q ≤ 13 and m ≥ 167, or the reverse. That gives
  N_w ≥ 48, and the equations no longer separate.
- **Bordered, with K ⊇ ℤ_p.** The admissible coset sums are c + 2e on K and c elsewhere, with
  |c| ≲ 2√v/q. A w-type fibre contributes about ±p to its coset, so N_w ≥ 1 is impossible.
- **Bordered, with K = M × {0} or K = K_M × {0}.** This needs p | n (family F) or p | n + 1 (family U)
  with p ≈ 2√(n/N_w). No unresolved n has such a p in family F. In family U, flatness then forces
  several w-type fibres into one block, and their p²-terms cannot be cancelled.
- **No w-type fibres.** Taking N_w = 0 forces n = m(p + 1), where Kronecker-type constructions already
  give the order.

## 3. Twisted structures: what was searched

**Census.** The family of 1916 and 3964 is G = ℤ₂² × ℤ_κ (family U, q = 4) with
S = ⟨(α, u₁), (1, u₂)⟩, where u₁ runs over all units. We also used ℤ₂³ × ℤ_κ (q = 8), with
automorphisms of order 7 (Singer) and 3. Built by `twinbuild.py` with 30–64 orbits, this gives 36,444
structures for 174 unresolved orders: 114 with n ≡ 3 (mod 4) and 60 with n ≡ 7 (mod 8).

**Results.** About 5,300 were run:

| Outcome | Count |
|---|---:|
| exactly empty | 857 |
| skipped: 10^8–10^12 compressed vectors at the first level | 4,274 |
| timed out | about 240 |
| exceeded memory caps (the level-2 lifts at n = 647 need 50–80 GB) | about 50 |

**Exact analogues of 3964.** These have a cubic-class fibre on ℤ_p and α paired with an order-3
multiplier of ℤ_m. All are empty at the fibre-side quotient:

| n | κ | structures |
|---|---|---:|
| 1039 | 7 · 37 | 10 |
| 1979 | 26 · 19 | 150 |

Other hosts of this type are larger and untested: 1823 (35 · 13), 1879 (7 · 67), 2215 (7 · 79),
2815 (19 · 37) and 2887 (7 · 103).

## 4. Search engine notes

- **Order of the levels.** For 3964 the fibre-side quotient G/ℤ₁₃ has 405 compressed vectors per
  role, against 3.2–4.6·10⁵ for G/ℤ₁₉. Joining the small side first makes level 1 instant but moves
  the cost into the lifts. The original search took ℤ₁₉ first (304 branches).
- **join4c.** It materialises every pair sum: about 1.5·10¹¹ per side for 3964's level 1, hours on
  two cores. A cap on the number of level-1 solutions gives an early exit, but solutions are rare per
  partition, so it helps little here.
- **Memory.** Container memory includes /dev/shm. Per-run address-space caps are needed when lifts
  can take tens of GB.
- **Threads.** Pin OPENBLAS_NUM_THREADS per process, or every run starts one BLAS thread per visible
  core.

## 5. Literature (search of October 2026)

- **Alpöge et al., 12 orders.** They were posted with a decoder and no paper. Speyer's analysis of 668
  on MathOverflow shows a bordered Goethals–Seidel array with a border of size 4 over ℤ₁₆₆.
- **Legendre pairs** have been computed up to length 147 (Kotsireas–Koutschan 2021). The largest
  Goethals–Seidel / SDS searches reached n = 631 (Đoković–Golubitsky–Kotsireas 2013, about 10¹⁴ pair
  sums).
- **Cyclotomic SDS families** (Lehmer 1976, Yamada 1992, Xia–Xia–Seberry–Wu 2004, Leung–Ma–Schmidt
  2006, Leung–Momihara–Xiang 2019) need conditions such as 4n − 1 a prime power or n = 2q² + 1, which
  the unresolved orders checked do not meet.
- **Solution counts.** We found no published model of the number of solutions of multiplier-invariant
  difference families.
- **Cocyclic matrices.** For order 4p with p ≡ 3 (mod 4), every cocyclic Hadamard matrix is of
  Williamson type or a transposed Ito matrix (Barrera Acevedo, Ó Catháin, Dietrich 2019).

## 6. The two-spike shape and T-matrix obstructions (October 13)

**Two-spike bordered shape** (the structure of order 2148 on the progress branch). It is a bordered
Goethals–Seidel array over ℤ₂q × GF(r) with Q = ℤ_q, T = all characters, index e and F = (r − 1)/e odd.

- The budget identity Σ(F s² + 2 f s) = 4qe and the target −4qr at the spike distance q force
  F ∈ {qe/2 − 1, qe/2, qe/2 + 1}, so n = q(2r + 1) with r = eF + 1.
- Among the unresolved n ≤ 2999 only 537 (q = 3, e = 8) and 2515 (q = 5, e = 10, r = 251, F = 25) have this form.
- For 2515, an exact meet-in-the-middle over the compressions to ℤ₁₀ finds no stage-1 class. It covers all eight
  sign choices of the spikes ±49, ±51 at distance 5 and all ±1 completions. The same code returns the single class
  of the branch for 537 and for its analogue 177.
- So the shape does not reach any other unresolved order up to 2999.

**Spence over ℤ₂ × GF(419) (n = 839).**
- Through the ℤ₂-split, the w-sums satisfy Σ s_w² = n, with s_w = ε + F a in the three blocks whose zero cell is a
  w-cell and s_w = F d in the fourth.
- This gives F A + 2B = 2e with A = Σ a² and |B| ≤ A. It has no solution for F ∈ {11, 19, 22, 38}.
- So no Spence structure over ℤ₈₃₈ with a multiplier group of order at least 11 exists.

**T-matrices.** By Cooper–Wallis, T-matrices of order t give H(4tw) for every Williamson order w.

- T(191) would settle 573, 955, 1719 and 2865. With Miyamoto's theorem, 573 also gives 2293.
- T(223) would settle 669, 1115 and 2899, and 669 gives 2677.
- T(263) would settle 789, 1315, 1841 and 2893.

For prime t and a multiplier group H ≤ ℤ_t^*, one T-matrix holds the point 0 and every other orbit has |H|
points, so the row sums satisfy (±1 + |H|a)² + |H|²(b² + c² + d²) = t. For t ≡ 3 (mod 4) the left side is
≡ 1 (mod 4) whenever |H| is even, so H has odd order. Hence:

- 191 admits no nontrivial multiplier group (|H| = 5, 19 and 95 also fail the equation);
- 263 admits none (|H| = 131 fails);
- 223 admits only |H| = 3 (rows −5, 9, 9, 6; 74 cosets, 37 classes).

## 7. Audits and Williamson-type routes (October 14)

**1823 is not open.** Seberry and Whiteman (Graphs Combin. 4 (1988) 355–377) construct symmetric conference matrices of
order 5·9^(2t+1) + 1. For t = 1 the order is 3646 = 2·1823, and [[C+I, C−I], [C−I, −C−I]] is a Hadamard matrix of
order 7292 = 4·1823; their abstract states this matrix.
- Seberry–Yamada's 1992 Table A.2 lists 1823 with t = 3, and Cati–Pasechnik and FRONTIER.md inherited the entry.
- We rebuilt C from the 3×3-cell iteration in Balonin's catalogue (one sign of the printed formula has to be flipped).
- Exact checks: C = Cᵀ, zero diagonal, C² = 3645 I and H Hᵀ = 7292 I.
- SHA-256 of the canonical CSV of H: `4f2289bc0e2ce8f8212457c3e1044b63396eca8179b0b0703debf0b0b365b9ef`.

**1747.** Handbook Table 1.53 (2007) prints t = 2. Seberry–Yamada (1992) give t = 3, from Miyamoto's 8q class, and no
construction for t = 2 could be traced. It is probably a table error, so 1747 stays on the list.

**Corollary 29(i) of Seberry–Yamada (1990).** It states: if q ≡ 1 (mod 4) is a prime power and WT((q−1)/4) or
H((q−1)/2) exists, then WT(q) exists.
- **What the proof forces.** The proof uses their Lemma 26 with Yamada's matrices as U₁ and U₂. This forces U₃ = U₄ = 0,
  V₁ = I and V₂ = 0. So V₃ and V₄ must be symmetric, commuting ±1 matrices with V₃² + V₄² = (q−1) I: a symmetric
  quaternary Hadamard matrix of order (q−1)/2.
- **Repeated block.** A WT((q−1)/4) with a repeated block (Y, X, X, Z) supplies V₃ = [[Y, X], [Xᵀ, −Z]] and
  V₄ = [[Z, X], [Xᵀ, −Y]]; X need not be symmetric. Checked for q = 29, 37, 53 and 61: WT(q) and the Williamson array
  H(4q) verify exactly.
- **Four distinct blocks.** For every class of Williamson matrices of order m ≤ 19, a signed 2×2 arrangement of the
  blocks exists exactly when a block repeats.
- **Order 67.** The three known WT(67) (near-Williamson, see below) give no arrangement. The repeated-block structure at
  67 with H = {1, 29, 37} is empty, by an exhaustive search. Other multiplier orders fail on row sums. So
  1883 = 7·269 is not reached this way.
- **Consequences.** Đoković (Int. Math. Forum 6 (2011)) settled 1527 = 3·509, 2545 = 5·509 and 2571 = 3·857 through
  Cor. 29(i), using WT(509) from good matrices of order 127, and H(428). Under the Lemma 26 proof these need a symmetric
  BH(254, 4) and a symmetric BH(428, 4); we found neither.
- **Open point.** We could not read Miyamoto's 1991 paper (J. Combin. Theory A 57), which is credited with the
  Williamson version. Unless it uses a different construction, these three orders are only conditionally settled.

**Spence-W: Williamson-type matrices from Spence sequences.** Take four ±1 sequences on ℤ₂ × ℤ_n with ΣPAF = −4 off zero
(Spence's condition for order 4(2n+1)), one invariant under a multiplier group H and three under ⟨H, −1⟩.
- **What they give.** Through the ℤ₂-split they give the U and V of Lemma 26, hence WT(2n+1). With any
  OD(4w; w, w, w, w) this gives a Hadamard matrix of order 4(2n+1)w.
- **Why it helps.** The symmetry comes from ℤ_n, so the route applies to primes 2n+1 whose own multiplier groups are
  useless. For example, T-matrices of order 191 or 263 admit no multiplier group (section 6).
- **Checks.** The Williamson array H(4(2n+1)) verifies exactly for each of:
  - WT(21) and WT(25), with n even and no multiplier;
  - WT(83), with n = 41 and |H| = 5.
- **269.** The structure with n = 134 and |H| = 11 (15 classes) is exhaustively empty, so it does not give 1883.

**Orders above 4·2999 from near-Williamson matrices.** Kharaghani, Mohammadian, Tayfeh-Rezaie and Zaitsev
(arXiv:2605.08661, 2026) give near-Williamson matrices: circulant ±1 matrices A, B, C, D with B, C, D symmetric and
ΣXXᵀ = 4w I. Then (AR, B, C, D) are Williamson-type matrices. They exist for w = 47, 59, 65, 67, 103 and 133, among
others. Williamson matrices do not exist for w = 35, 47, 53 and 59 (Đoković 1993; Holzmann, Kharaghani and
Tayfeh-Rezaie 2008) and have not been found for w = 65 and 67 (the paper's §2.3: "either known not to exist or have
not yet been found, namely n ∈ {53, 59, 65, 67}").
- **Five new orders.** Combined with the T-matrices of orders 103, 151 and 163 by Cooper–Wallis, they give 24308,
  27604, 35636, 38468 and 39260. See their READMEs.
- **Further candidates** (by the same audit):
  - 4·3149 = 4·47·67, 4·4757 = 4·67·71, 4·4841 = 4·47·103, using T-sequences of lengths 47 and 71;
  - 4·8549 = 4·83·103, which needs only T103 and the classical WT(83).
- None of these products lands on an unresolved n ≤ 2999.
- All four candidates are now built: [12596](12596) and [19364](19364) with the bordered WT(47) of
  [williamson-type](williamson-type), [34196](34196) with the bordered WT(83), and [19028](19028) with T-sequences of
  length 71 (from Turyn-type TT(24)) and the near-Williamson WT(67).
