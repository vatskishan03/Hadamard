# Notes on structured searches (October 2026)

What we learned while searching for the unresolved orders in [FRONTIER.md](FRONTIER.md), including the searches
that found nothing. Section 1 is a mechanism, sections 2–4 and 6 are negative results, and section 5 is the literature
as of October 2026. Section 7 collects corrections to the tables and the Williamson-type routes, sections 8–12 the
complex Hadamard matrices, section 13 the Spence spin structures behind order 2292, and section 14 the complex orders 330 and 426. The notes on the individual
search programs are in [search/NOTES.md](search/NOTES.md).

## 1. Why some structures have many solutions

Take a group G = M × ℤ_p, p prime, and a symmetry group S of the four sequences. Let C = S ∩ ℤ_p* act on
the ℤ_p coordinate, f = (p − 1)/|C|, and t = |π_p(S)|/|C|, where π_p is the projection to ℤ_p*. On every
fibre {a} × ℤ_p a sequence is constant on the cosets of C, so it is x₀(a) on 0 and v₀(a), …, v_{f−1}(a) on
the f cyclotomic classes. Expanding the class vector in the characters of ℤ_p*/C ≅ ℤ_f, every nonzero
mode pairs with a Gauss sum of absolute value √p, and the Goethals–Seidel condition splits by the
difference of the modes.

- **f = 2 (a Paley product).** Every fibre is either constant ("w-type") or ±(Legendre sequence)
  ("z-type"). The identity coefficient forces (p + 1)N_w − 4N₋ = 4|M| for the numbers of w-type fibres,
  and the equations split by powers of p. This is what makes the skew Paley-product orders of [skew](skew) easy.
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
class whose character is nontrivial on ℤ_p. On the exact count for v = 253 ([search/NOTES.md](search/NOTES.md) §15)
this brings the prediction to within 10^0.8–10^1.8 of the truth; the uncorrected model is off by 10^6–10^9.

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
- **Memory.** Files in /dev/shm count against the memory limit. Per-run address-space caps are needed when lifts
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

**Two-spike bordered shape** (the structure of order [2148](hadamard/2148)). It is a bordered
Goethals–Seidel array over ℤ₂q × GF(r) with Q = ℤ_q, T = all characters, index e and F = (r − 1)/e odd.

- The budget identity Σ(F s² + 2 f s) = 4qe and the target −4qr at the spike distance q force
  F ∈ {qe/2 − 1, qe/2, qe/2 + 1}, so n = q(2r + 1) with r = eF + 1.
- Among the unresolved n ≤ 2999 only 537 (q = 3, e = 8) and 2515 (q = 5, e = 10, r = 251, F = 25) have this form.
- For 2515, an exact meet-in-the-middle over the compressions to ℤ₁₀ finds no stage-1 class. It covers all eight
  sign choices of the spikes ±49, ±51 at distance 5 and all ±1 completions. The same code returns the single
  stage-1 class of the 2148 search for 537, and the class of its analogue 177.
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

**1823 is not open.** Seberry and Whiteman (Graphs Combin. 4 (1988) 355–377) construct symmetric conference matrices
of order 5·9^(2t+1) + 1. For t = 1 the order is 3646 = 2·1823, and [[C+I, C−I], [C−I, −C−I]] is a Hadamard matrix of
order 7292 = 4·1823; their abstract states this matrix.
- Seberry–Yamada's 1992 Table A.2 lists 1823 with t = 3, and Cati–Pasechnik and FRONTIER.md inherited the entry.
- We rebuilt C from the 3×3-cell iteration in Balonin's catalogue (one sign of the printed formula has to be flipped).
- Exact checks: C = Cᵀ, zero diagonal, C² = 3645 I and H Hᵀ = 7292 I.
- SHA-256 of the canonical CSV of H: `4f2289bc0e2ce8f8212457c3e1044b63396eca8179b0b0703debf0b0b365b9ef`.

**1747.** Handbook Table 1.53 (2007) prints t = 2. Seberry–Yamada (1992) give t = 3, from Miyamoto's 8q class, and no
construction for t = 2 could be traced. It is probably a table error, so 1747 stays on the list.

**Corollary 29(i) of Seberry–Yamada (1990).** It states: if q ≡ 1 (mod 4) is a prime power and WT((q−1)/4) or
H((q−1)/2) exists, then WT(q) exists.
- **What the proof forces.** The proof uses their Lemma 26 with Yamada's matrices as U₁ and U₂. This forces U₃ = U₄ =
  0, V₁ = I and V₂ = 0. So V₃ and V₄ must be symmetric, commuting ±1 matrices with V₃² + V₄² = (q−1) I: a symmetric
  quaternary Hadamard matrix of order (q−1)/2.
- **Repeated block.** A WT((q−1)/4) with a repeated block (Y, X, X, Z) supplies V₃ = [[Y, X], [Xᵀ, −Z]] and V₄ = [[Z,
  X], [Xᵀ, −Y]]; X need not be symmetric. Checked for q = 29, 37, 53 and 61: WT(q) and the Williamson array H(4q)
  verify exactly.
- **Four distinct blocks.** For every class of Williamson matrices of order m ≤ 19, a signed 2×2 arrangement of the
  blocks exists exactly when a block repeats.
- **Order 67.** The three known WT(67) (near-Williamson, see below) give no arrangement. The repeated-block structure
  at 67 with H = {1, 29, 37} is empty, by an exhaustive search. Other multiplier orders fail on row sums. So 1883 =
  7·269 is not reached this way.
- **Consequences.** Đoković (Int. Math. Forum 6 (2011)) settled 1527 = 3·509, 2545 = 5·509 and 2571 = 3·857 through
  Cor. 29(i), using WT(509) from good matrices of order 127, and H(428). Under the Lemma 26 proof these need a
  symmetric BH(254, 4) and a symmetric BH(428, 4); we found neither.
- **Open point.** We could not read Miyamoto's 1991 paper (J. Combin. Theory A 57), which is credited with the
  Williamson version. Unless it uses a different construction, these three orders are only conditionally settled.

**Spence-W: Williamson-type matrices from Spence sequences.** Take four ±1 sequences on ℤ₂ × ℤ_n with ΣPAF = −4 off
zero (Spence's condition for order 4(2n+1)), one invariant under a multiplier group H and three under ⟨H, −1⟩.
- **What they give.** Through the ℤ₂-split they give the U and V of Lemma 26, hence WT(2n+1). With any OD(4w; w, w, w,
  w) this gives a Hadamard matrix of order 4(2n+1)w.
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
- All four candidates are now built: [12596](hadamard/12596) and [19364](hadamard/19364) with the bordered WT(47) of
  [williamson-type](williamson-type), [34196](hadamard/34196) with the bordered WT(83), and [19028](hadamard/19028)
  with T-sequences of length 71 (from Turyn-type TT(24)) and the near-Williamson WT(67).

## 8. Complex Hadamard matrices, the order-3 twist and quaternary Legendre pairs (October 14)

**Szöllősi's two-symmetric array.** If A, B, C, D are circulant ±1 matrices of odd order N with ΣXXᵀ = 4N I and only A
and B symmetric, Theorem 4 of arXiv:2603.09572 gives a quaternary Hadamard matrix of order 2N. Williamson-type
matrices (all four symmetric up to R) are not needed. With a multiplier group H ∌ −1, C and D are H-invariant, A and B
are ⟨H, −1⟩-invariant, and the condition lives on the ⟨H, −1⟩-classes (`search/ssxx.py`, exact four-list joins).
- **N = 71**, H = ⟨5⟩ of order 5. The Gaussian model gives 10^3.8 per row pattern; the exact joins find 5992 and 6076
  solutions. This gives order 142, which no Williamson-type matrices of order 71 provide: Table A.1 of Seberry–Yamada
  1992 is blank at 71, and the 2026 near-Williamson results stop at 69.
- **N = 77**, H = ⟨15⟩ of order 5, with 15 ≡ 1 (mod 7) and the squares mod 11. Here the model is wrong by four
  orders: it predicts 10^−1.7 to 10^−2.1, and the exact joins find 2352, 2352, 2352 and 1176 solutions in four
  patterns. These are Paley-product structures, X = u ⊗ 1 + w ⊗ J + z ⊗ χ over ℤ₇ × ℤ₁₁.
  [search/NOTES.md](search/NOTES.md) §15 shows that the per-class densities are strongly correlated for them.
- **The sweep** (`search/ssxx_sweep.py`, exact joins in order of cost; `search/ssxx_sample.py`, uniform samples of the
  large lists when they cannot be listed). For prime N the Gaussian estimate of the hits inside the samples is
  accurate: at 71 it predicts 1445 and the join finds 1430, and at 131 it predicts 10^4.0 and the join finds 12618.
  Results so far, with the exact count in the structure that gave each order:
  - 238 (N = 119, H = ⟨8⟩, 4704);
  - 262 (N = 131, H = ⟨53⟩, sampled);
  - 286 (N = 143, H = ⟨3⟩, 16);
  - 326 (N = 163, H = ⟨38⟩, 8352).
- **Turyn's theorem** gives a real Hadamard matrix of order 4N from a complex one of order 2N. So the sweep cannot
  succeed at an unresolved real order (573, 599, …) without settling that order. It succeeds where real H(4N) is
  known by constructions that do not produce complex matrices (Paley I, Miyamoto, products).
- **Safe primes are out of reach.** For N with (N − 1)/2 prime (107, 167, 179, 227, 263, …), the only subgroup of
  odd order is the squares. Then A and B are constant on ℤ_N*, C and D are of Legendre type, and their PAFs cannot
  sum to 0.
- **Known orders.** The order 2N counts as known if it comes from one of:
  - Williamson-type matrices of order N;
  - a Paley conference matrix (2N − 1 a prime power);
  - complex Golay sequences;
  - a quaternary Legendre pair of length N − 1;
  - Szöllősi's 94.

  The sweep over the remaining odd N < 1000 is `search/ssxx_sweep.py`.

**The order-3 twist on the bordered form** (`search/wt_twist.py`, `search/twist_join.c`). The twist of
arXiv:2605.08661 §4.4, type (ii), is (A; B, B∘u, B∘u²) with u of order 3. It can be placed on ℤ_(N−1), where the
bordered construction of [williamson-type](williamson-type) applies; for prime N with 3 ∤ N − 1 there is no order-3
multiplier on ℤ_N itself. The B-list is streamed in Gray order against a hash table of the A-list, with an O(1) hash
update per step, so the search of a structure is exact. Exact counts agree with brute force (N = 31 plain: 180 and
510; N = 19 bordered: 72), and the hits rebuild verified Williamson arrays. The twist cannot reach N ≡ 2 (mod 3):
- **Bordered.** The rows are (±2, 0, 0, 0), and the reduction ℤ_m → ℤ₂ gives d_A² + 3d_B² = 4N with d = (even sum) −
  (odd sum). A symmetric B with row sum 0 has d_B = 2e_B, so N = (d_A/2)² + 3e_B². This holds for m ≡ 2 (mod 4), where
  e_B is odd, and for m ≡ 0 (mod 4), where e_B is even.
- **Plain.** Here r_A² + 3r_B² = 4N.
- **Consequence.** In both forms N must be a norm from ℤ[ω], so N ≢ 2 (mod 3).

This excludes 71, 77, 107, 119, 131, 143, 155, 161, 167, 179, 185 and 191, among them 191 (573) and 263 (789).
**Order 111.** The plain twist at N = 111 (u = 10, K′ = ⟨38⟩, 13 classes) predicted 10^−0.16 and holds 36 quadruples
with rows (21, 1, 1, 1): near-Williamson matrices of order 111
([near-williamson/NW111.txt](near-williamson/NW111.txt), 23 minutes for the 2^38-step B-stream). They give
Williamson-type matrices of order 111, which are not among the 2026 near-Williamson orders and appear in
Table A.1 only with t = 3, and through them the complex Hadamard order 222. The
other two orders those FRONTIER rows need, 223 and 163, fail on cost:
- 223, which would settle 669, 1115 and 2899, has one admissible twisted structure, K′ = 1, with 10^0.2 expected
  solutions but A- and B-lists of 10^22 and 10^32.
- 163 is at most 10^−17 in every twisted structure.

**Quaternary Legendre pairs** (qLP; open lengths up to 100 after Lebedev and Verma 2026: 70, 76, 88, 92, 94 and 100).
In Gray rows A = ((p + q) + i(p − q))/2 and B = ((r + s) + i(r − s))/2, a qLP is ΣPAF = −4 over the four rows plus an
amicability condition. Amicability is automatic when each pair (p, q), (r, s) is symmetric about a common centre.
Lebedev exhausted the family with all centres at 0, at lengths 64, 70 and 76.

*Half-integer centres.* A row symmetric about −1/2 pairs each even position with an odd one, so its alternating sum
vanishes. A row symmetric about 0 has alternating sum d ≡ 2 (mod 4) or ≡ 0 (mod 4) by its row sum. The condition at
frequency ℓ/2, Σ_rows d² = 4ℓ + 4, then forces:
- All four rows centred at −1/2: impossible.
- One pair at each kind of centre: (d_p/2)² + (d_q/2)² = ℓ + 1, so ℓ + 1 must be a sum of two squares.
  - Impossible for ℓ = 70 (71), 76 (77), 92 (93) and 94 (95).
  - Possible for ℓ = 88 (89) and ℓ = 100 (101). This mixed family lies outside Lebedev's queues.

So at lengths 70 and 76 a qLP must have a pair of Gray rows that is not symmetric about a common centre.

*A binary.* Verma's length-64 pair has A binary. For binary A and symmetric (r, s) the condition is
2PAF_A + PAF_r + PAF_s = −4, with r_A = 0 and (r_r, r_s) = (2, 0). `search/qlp_bin.py` searches it with A on the
H-orbits and r, s on the ⟨H, −1⟩-orbits; it finds pairs at lengths 14 and 26 and verifies them exactly. At the open
lengths every multiplier structure has a Gaussian estimate of 10^−2.2 or less: 70: 10^−4.4 (⟨11⟩); 76: 10^−2.3 (⟨45⟩);
94: 10^−2.2 (⟨3⟩).

## 9. Bordered Szöllősi arrays and routes to the open orders (October 14)

**Which new objects would settle FRONTIER orders.** Factoring the 145 open n ≤ 2999 of FRONTIER.md as n = t·w
against the known T-matrix orders and Williamson-type orders (Table A.1 of 1992, the 2026 near-Williamson orders,
47, 83, 111) shows that no product is overlooked. A single new order would settle several entries. Every such w is a
prime ≡ 3 (mod 4): for q ≡ 1 (mod 4) Williamson-type matrices of order (q + 1)/2 are classical, and none of these
2w − 1 is a prime power.

| new object | settles (Cooper–Wallis with the other factor) |
|---|---|
| WT(263) or T(263) | 789, 1315, 1841, 2893 |
| WT(191) or T(191) | 955, 1719, 2865 |
| WT(223) or T(223) | 669, 1115, 2899 |
| WT(167), WT(179), WT(227), WT(359) | 1169 and 2171; 1969 and 2327; 1589 and 2951; 1795 and 2513 |

The bordered Williamson-type structures of [williamson-type](williamson-type) rate every multiplier group at about
10^−4 for these w. The FRONTIER rows 7516, 8860 and 9532 need skew matrices of orders 940, 1108 and 1192 (Kiyasu's
theorem as proved by Yamada, J. Austral. Math. Soc. 1989, Theorem 4: q ≡ 5 (mod 8) a prime power and a skew matrix
of order (q + 3)/2 give H(4(q + 2))). These skew orders (v = 235, 277, 149) are out of reach of the skew searches
([search/NOTES.md](search/NOTES.md) §10).

**A 3-bordered Williamson-type structure does not exist.** Try G = ℤ₂² × ℤ_m with p = 4m + 3. Take three border
points whose border rows and columns are constant on the cosets of K = 0 × ℤ_m, and group-developed blocks with
three of them symmetric. The ΣPAF condition is then the one of the 3964 construction: −12 on K \ 0 and +4 off K.
But amicability fails:
- For each nontrivial character χ of ℤ₂², Σ_k x̂_k(χ)² = 12 forces three of the four character sums to be ±2.
- Every assignment of characters to border rows then leads to a contradiction in the χ-coefficients of the
  cross-term equations. The case analysis covers one, two, three or four blocks per border coordinate.

A SAT model of the general coset-constant border (CaDiCaL, `pysat`) confirms this: it is unsatisfiable at p = 11.
The GS version escapes because the Goethals–Seidel array needs no amicability.

**Bordered Szöllősi array (`search/bsz.py`, `search/bsz_group.py`).** Spence's bordered condition carries over to
Szöllősi's complex array with only two symmetric blocks. Let A, B, C, D be ±1 sequences on an abelian group G of
order N − 1 with:
- A and B symmetric;
- ΣPAF = −4 off 0;
- hence row sums (±2, 0, 0, 0) in some order.

Put X′ = αA + ᾱB and Y′ = αC + ᾱDR (α = (1 + i)/2, group-developed matrices, R the matrix of g ↦ −g) and

    X = [[x, s eᵀ], [t e, X′]],   Y = [[y, u eᵀ], [v e, Y′]],   K = [[X, Y], [Y*, −X*]].

Then K is a complex Hadamard matrix of order 2N when the border entries x, s, t, y, u, v ∈ {±1, ±i} satisfy five
scalar equations in the eigenvalues ρ_X = αr_A + ᾱr_B and ρ_Y = αr_C + ᾱr_D of e. The derivation:
- X′ is symmetric, so X′Y′ = Y′X′.
- X′X′* + Y′Y′* = (AAᵀ + BBᵀ + CCᵀ + DDᵀ)/2 = 2NI − 2J.

For r_A = 2 one solution is x = 1, s = −t̄, u = −iyt̄, v = ity, with t and y free. Every row-sum pattern has 32
border solutions.

Checks:
- Exhaustively over ℤ₆: all 2304 quadruples give complex Hadamard matrices of order 14.
- Over ℤ₂² × ℤ_m: orders 26, 42 and 50, with 12672 and 4.6 · 10⁷ solutions at N = 13 and 21 against Gaussian
  estimates of 10^4.25 and 10^7.02.
- The model misses arithmetic obstructions: N = 17 over ℤ₂² × ℤ₄ is empty against 10^5.6.

`verify/verify_complex.py` accepts the constructions `szollosi-bordered` and `szollosi-bordered-group`.

This is the complex analogue of the hierarchy GS (no symmetric block) → Szöllősi (two) → near-Williamson (three)
on ℤ_N, one level down: Spence (none) → bordered Szöllősi (two) → bordered Williamson-type (three) on a group of
order N − 1. The point is the same as for [williamson-type](williamson-type): it moves the search to a group of
order N − 1, which helps exactly where ℤ_N has no usable multipliers. `search/bsz_sweep.py` runs every cyclic
multiplier structure and every spin structure on ℤ₂² × ℤ_((N−1)/4) whose complete lists join within 2 · 10¹⁰ pair
operations, for the open N.

**Order 214 (N = 107) is now the smallest open complex Hadamard order** in our accounting. It is blank in Table
11.2. No near-Williamson or quaternary Legendre pair of the matching size (ℓ = 106) is known to us, and 107 is a safe
prime:
- ℤ₁₀₇* and ℤ₁₀₆* have no useful multiplier groups.
- The only bordered structure the row sums admit, H = {±1}, rates 10^1.0 but has pools of 5 · 10¹⁴.

## 10. Exact compression chains and streaming joins for complex Hadamard orders (October 14)

**Compression chain (`search/cchain4.py`, `search/cchain_sweep.py`).** The method of order 3964 (and of
`search/chain.py` for skew families), generalised to four blocks with per-block symmetry. Blocks are of kind L
(⟨H, −1⟩-orbits) or A (H-orbits). Four problems are covered:

| mode | group | blocks | target ΣPAF |
|---|---|---|---|
| Szöllősi | ℤ_N | L L A A | 0 |
| bordered Szöllősi | ℤ_(N−1) | L L A A | −4 |
| near-Williamson | ℤ_N | A L L L | 0 |
| bordered Williamson-type | ℤ_(N−1) | A L L L | −4 |

For d | v the compressions to ℤ_d satisfy the same identity with target T·v/d. The steps:
- A dynamic program lists the joint compressions to ℤ_d₁ × ℤ_d₂ of each block, with multiplicities.
- The level-1 join is exact.
- For each level-1 solution, the level-2 join inside the fibres is exact.
- The fibres are expanded by meet in the middle over the orbits.
- The last join groups each list by its compressions to the prime-power levels the joint program left out,
  joins them exactly, and only then runs join4x on all classes.

Nothing is sampled. A pattern with an empty level-1 join is ruled out at once, and most are. Positive controls
(N = 77, 143, 185, 247 and bordered N = 21) recover known solutions in seconds.

*Order 310.* N = 155 = 5 · 31 with H = ⟨16⟩ of order 5 has 35 orbits and lists of 1.3 · 10⁴, 6.4 · 10³, 7 · 10⁸ and 6
· 10⁸, far beyond the exact join. The Gaussian estimate is 10^1.4 for several row-sum patterns, but the best of these,
(9, 13; 3, 19) and (7, 11; 3, 21), are ruled out by their ℤ₅ compressions. Patterns (9, 9; 13, 17) and (9, 13; 9, 17)
pass level 1 with 600 branches each, and both gave verified solutions within 5 minutes (`complex-hadamard/C310.json`).

**Streaming k-tree (`search/kstream.c`, `search/kstream.py`).** For prime N there is no compression. When A and B
(symmetric) are small but C and D have 10¹⁰ or more members, C and D are enumerated without being stored:
- a depth-first search over their orbits, with class values updated incrementally (vectorised int32 tables);
- pass 1 joins each C to the A list on k classes, for a set of targets near the centre, and stores
  (fingerprint of the other classes, A, C);
- pass 2 streams D against B on the same classes and probes the stored fingerprints;
- every match is verified on all classes, then as a complex Hadamard matrix.

A pseudo-random fraction of the depth-first prefixes may be enumerated in either pass. With the pass-1 table held
at a fixed size, the expected yield is E · (fraction of D) · P(a + c ∈ targets) and does not depend on the fraction
of C. Positive control N = 131 (⟨53⟩): 4372 verified solutions from complete passes in 4 minutes, 144 from 10% / 20%
of the two blocks in 34 s.

The target is N = 191, the only open prime with a large estimate. H = ⟨39⟩ of order 5 gives 39 orbits and
19 classes. Lists are 9.2 · 10⁴ for A and B and 3.4 · 10¹⁰, 2.9 · 10¹⁰ for C and D. The Gaussian estimate is 10^7.55
solutions, and for prime N the estimate has been accurate (71: 1445 predicted, 1430 found; 131: 10^4.0, 12618).

*Order 382.* Run: `python3 search/kstream.py 191 39 11,11,9,21 --k 7 --wide 5 --vals 2 --fracC 0.18 --fracD 0.3
--maxm1 2.2e8 --threads 4`. Pass 1 joined 18% of the C pool with A on the first 7 classes. It kept 1.73 · 10⁸
fingerprints in 36 minutes on 4 cores. Pass 2 found a verified solution 6.5 minutes into the D stream (30% of the D
pool was scheduled). This gives the complex Hadamard matrix of order 382 (`complex-hadamard/C382.json`). The order is
blank in Table 11.2, and 191 is prime, so no compression or product reaches it.

## 11. Twisted Szöllősi quadruples (October 14)

**Construction** (`search/tsz.py`, `search/tsz.c`). Let H ≤ ℤ_N^* with −1 ∉ H, L = ⟨H, −1⟩, and μ a unit with μ ∉ L
and μ² ∈ L. Every quadruple below satisfies Szöllősi's hypotheses (A, B symmetric, ΣPAF = 0 off 0), so it gives a
complex Hadamard matrix of order 2N and a Goethals–Seidel matrix of order 4N.

- *Full twist.* A is L-invariant, C is H-invariant, B = A∘μ, D = C∘μ.
  - B is symmetric because A is, and PAF_{X∘μ}(s) = PAF_X(μs).
  - So ΣPAF(s) = F(s) + F(μs) with F = PAF_A + PAF_C, and the condition holds on the orbits of ⟨L, μ⟩.
  - The row sums satisfy r_A² + r_C² = 2N, so N is a sum of two squares.
  - Parity: F ≡ 2N ≡ 2 (mod 4). So no class may be mapped to itself by μ (on such a class the condition reads F(s) =
    0).
- *Half twist.* A and B are invariant under L′ = ⟨L, μ⟩, C is H-invariant, D = C∘μ.
  - The condition PAF_A + PAF_B + PAF_C(s) + PAF_C(μs) = 0 holds on the L′-classes.
  - The row sums satisfy r_A² + r_B² + 2r_C² = 4N.
  - It needs [L′ : L] = 2, which is impossible for a prime N ≡ 3 (mod 4): there ℤ_N^* is cyclic of order 2 · odd.
    Every composite N with two distinct prime factors has such twists.
- *Bordered half twist.* The same on ℤ_{N−1} with ΣPAF = −4 and row sums (0, 2, 0)
  ([BORDERED.md](complex-hadamard/BORDERED.md)). The full twist cannot be bordered: it forces r_A = r_B.

**Cost against yield.** The unknowns are about 1.5 · O bits, where O is the number of H-orbits. The constraints are
about O/4 classes. The search is a two-list join: the A side is hashed (or, for the half twist, the pairs (A, B)),
and the C pool is streamed by the depth-first search of `kstream.c`. Its cost is about the size of the C pool, so
structures with O ≤ 40 are exhaustive in minutes. The untwisted search has twice the exponent in its estimate
(about 2^0.7O against 2^0.35O) but needs a four-list join over pools of 2^O. At equal cost the untwisted structure
wins; the twisted structures add other pairs (H, μ), reachable where the untwisted join is not.

**Calibration** (`python3 search/tsz.py test`; every hit rebuilt and checked exactly):

| N | kind | Gaussian estimate | solutions found |
|---:|---|---:|---:|
| 13 | full, H = 1 | 10^3.15 | 780 |
| 17 | full, H = 1 | 10^3.87 | 6256 |
| 29 | full, H = 1 | 10^5.91 | 989016 |
| 41 | full, H = ⟨10⟩ | 10^1.76 | 48 |
| 61 | full, H = ⟨13⟩ | 10^3.57 | 4680 |
| 21 | half, H = 1 | 10^2.61 | 1008 |
| 39 | half, H = ⟨25⟩ | 10^1.01 | 144 |
| 51 | half, H = ⟨16⟩ | 10^1.41 | 288 |
| 27 | bordered half, H = 1 | 10^4.42 | 406224 |
| 31 | bordered half, H = 1 | 10^3.83 | 61440 |

As for the untwisted array, the estimate is accurate for primes and low for composites.

**A two-block array that does not help.** With X = αA + ᾱDR and Y = αB + ᾱCR (α = (1 + i)/2, R the reversal), XX* +
YY* equals (ΣXXᵀ)/2 for *any* four circulants. Then K = [[X, Y], [−Ȳ, X̄]] is complex Hadamard iff ΣPAF = 0 and ABᵀ +
CDᵀ is symmetric. This needs neither A nor B symmetric, but it replaces 2 · (N − 1)/2 symmetry bits by (N − 1)/2
integer constraints of about 4 bits each. Its estimate is lower than Szöllősi's by about 2^−N, so it was not searched.

## 12. Low-energy k-tree for Szöllősi structures (October 14)

`search/lowe_sz.py` (with `search/lowe_enum.c` and the four-list join `search/kt4.c` of the low-energy skew searches,
[search/NOTES.md](search/NOTES.md) §20) for structures whose blocks C, D have 10⁹–10¹¹ members.
- Every member of a big block is visited by the depth-first search of `kstream.c`.
- A member is kept when its energy E_k(x) = (x − z_k)ᵀ Σ_{−k}⁺ (x − z_k) is small, where z_k = t − Σ_{j≠k} μ_j and
  Σ_{−k} is the covariance of the other three blocks' sum. Under the Gaussian model, exp(−E_k/2) is proportional to
  the probability that the other three blocks complete x to a solution, so the L members of lowest energy carry far
  more than their share L/|pool| of the solutions.
- The kept lists are joined by `kt4` with the pairing (A, C), (B, D).
- The estimate inside the lists is E times the retention of each filtered block: the probability, under the
  conditional law N(c_k, Σ_k − Σ_k Σ⁺ Σ_k) of a solution's member, that it passes the filter.

Calibration, N = 131, H = ⟨53⟩, rows (9, 9; 1, 19):
- Keeping 2 · 10⁵ of the 1.0 · 10⁷ C members and 2 · 10⁵ of the 7.7 · 10⁶ D members retains 4.9% and 5.7% of
  solutions' members, about 2.4× their share.
- The estimate inside the lists is 10^3.37.
- kt4 found 8 verified solutions within its first split values (whole run: 30 s).
- An exhaustive join of the same lists (all 3.4 · 10⁸ pairs per side, kt4 on two merge classes) counts 2560
  solutions inside them, against the 10^3.37 = 2344 predicted. The retention model is accurate.

At N = 183, H = ⟨34⟩, keeping 2 · 10⁷ of 2.6 · 10¹⁰ members retains 3–5% of solutions' members: a boost of about 50
per block. That turns the structure's 10^5.0 into 10^2.2 solutions inside the lists, against about 0.1 expected from
the streaming k-tree of §10.

## 13. Spence spin structures: orders 2292 and 9172 (October 14)

**Structure.** Spence's bordered construction over a group G of order n − 1 (ΣPAF = −4 off zero, row sums 2, 0, 0, 0),
with blocks (A, B, B∘μ, B∘μ²). Here H is a group of automorphisms, μ normalises H and μ³ ∈ H. A is invariant under K =
⟨H, μ⟩ and B under H. The condition holds on the orbits of ⟨K, −1⟩.
- **Small cases.** Exhaustive solution sets exist at n = 13 (G = ℤ₂² × ℤ₃, 720 solutions) and n = 19 (ℤ₁₈, 23328).
  Both give Hadamard matrices through `build_spence`.
- **Arithmetic gate.** A real character χ with χ∘μ in the H-orbit of χ forces n = a² + 3b². For the other n, μ must
  act without fixed points on Hom(G, ℤ₂), which needs every part size of the 2-part of G to occur an even number of
  times. That leaves 17 of the non-Eisenstein open n, among them 573, 789, 1169, 1437, 1473 and 1509.

**573.** G = ℤ₂² × ℤ₁₄₃, H = ⟨27⟩ (quadratic residues on ℤ₁₁), μ = (α, ×3) with α of order 3 on ℤ₂².
- **Size.** 27 classes, 54 K-orbits and 156 H-orbits. The Gaussian estimate is 10^7.5 solutions (10^10.9 with the full
  covariance). A plain join would cost about 10^34.
- **Exact gates.** For μ-invariant N, the compressed vectors to G/N are products of per-orbit value sets: each orbit
  lies over one orbit of G/N and meets each fibre in the same number of points. The structure passes ℤ₂², ℤ₁₁, ℤ₁₃,
  ℤ₂² × ℤ₁₁, ℤ₁₄₃ and ℤ₂² × ℤ₁₃, and a lattice test on the class values.
- **Level ℤ₂² × ℤ₁₃ (52 points).** A hash join at ℤ₁₃ is followed by an energy-matched search: A drawn at a fixed
  energy, B columns restricted to the complementary energies. This gives a verified solution every 20–60 CPU seconds;
  without the energy matching, 3·10¹⁵ pairs gave none.
- **Lift to G.** The fibre ℤ₁₁ is rigid. A compressed value ±1 forces (b₀, s·χ(ρ)) with one free sign s, and ±9, ±11
  force everything. This leaves 61–63 signs, 8 linear equations and the rest quadratic.
- **Lift solver.** An exhaustive four-list meet in the middle on the linear equations (about 10⁷ solutions per lift,
  15–40 s), followed by the full 27-class check. On planted targets it returned the planted signs every time.
- **Results.** Two of the first five level solutions lifted, against the Gaussian estimate of about 0.2 each. Both
  give Hadamard matrices of order 2292; the first is [2292](hadamard/2292). Miyamoto with q = 2293 gives
  [9172](hadamard/9172).

**Other survivors** of the same gates, not yet searched:
- 1169 over ℤ₂⁴ × ℤ₇₃ (Singer cycle of order 15 coupled with an order-3 multiplier): 17 classes, estimate 10^4–10^6,
  ℤ₇₃ gate undecided.
- 1169 over ℤ₄² × ℤ₇₃: 27 classes, estimate 10^2–10^9.
- 789 over ℤ₂² × ℤ₁₉₇: 29 classes, estimate 10^0–10^2.

Every other open n tested is excluded by the arithmetic gate, by row sums, or at a coarse quotient. Coverage is
incomplete for some non-cyclic groups (ℤ₂⁴ beyond ΓL(2,4), ℤ₃³, ℤ₅³).

## 14. The complex orders 330 and 426 (October 15)

For N = m q and H the squares of GF(q) acting trivially on ℤ_m, every block of Szöllősi's array is u ⊗ δ + w ⊗ J +
z ⊗ χ, with at each a ∈ ℤ_m either a z-position (w = 0, u, z = ±1) or a w-position (z = 0, w = ±1, u ∈ {0, −2w}).
The condition splits into three equations over ℤ_m; the compressions y = u + q w carry the first two, and their value
at 0 is the counting identity (q + 1) N_w − 4 N₋ = 4m of §1. `search/ssxx_paley.py` lists the w-skeletons up to
symmetry and joins the remaining ±1 values exactly in two stages; it reproduces the exact solution counts of
`search/ssxx.py` at N = 65 and 77.
- **N = 165 = 15 · 11** (q = 11, N_w = 6): verified solutions, hence the complex Hadamard matrix of order 330
  ([C330.json](complex-hadamard/C330.json)), blank in the 1992 table and in every family we checked.
- **Empty**: N = 171, 253, 275, 323, 345, 391, 513 and 629 (every admissible N_w, exhaustive), and 207 by the counting
  identity alone. As for the skew Paley products, small m is rich and the cases q ≈ 2m are thin
  ([search/NOTES.md](search/NOTES.md) §41).
- **N = 213 = 3 · 71**, H = ⟨37⟩ of order 7, by the low-energy k-tree of §12 after a fix: its split values had been
  computed for the wrong pair of lists, so the first 183 run examined no pair at all. Lists of 6435, 6435, 3 · 10⁷ and
  3 · 10⁷ (kept from pools of 4.4 · 10⁸), 10^0.95 solutions predicted inside them, 72 found in six minutes: the complex
  Hadamard matrix of order 426 ([C426.json](complex-hadamard/C426.json)).
