# Notes on structured searches (October 2026)

What we learned while searching for the unresolved orders in [FRONTIER.md](FRONTIER.md), including the searches that
found nothing. Section 1 is a mechanism, sections 2–4 and 6 are negative results, and section 5 is the literature as of
October 2026. Section 7 collects corrections to the tables and the Williamson-type routes, sections 8–12 the complex
Hadamard matrices, section 13 the Spence spin structures behind order 2292, section 14 the complex orders 330 and 426,
section 15 the compressions that decide which row patterns live, with the complex orders 342, 366, 430 and 454, section
16 the character-bordered array and the complex order 214, section 17 its ℤ₆ form and the complex order 322, section 18
Paley's matrix halves, section 19 Miyamoto's lemma with symmetric inputs, section 20 the page of Zaitsev et al. and
section 21 the quaternion type of Paley's matrix. The notes on the individual search programs are in
[search/NOTES.md](search/NOTES.md).

**Correction (October 16).** Sections 8–17 call the complex orders 142, 154, 222, 262, 286, 310, 322, 330, 342, 370,
430, 454 and 742 open before our searches. They are not: 4N − 1 is a prime power for each, and Ito's Hadamard groups
with Egan 2019, Theorem 4.5, give them (§18). The other 11, 214, 238, 302, 326, 366, 382, 426, 494, 518, 534 and
754, were constructed concurrently by Zaitsev et al. (§20).

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
- **Consequences (corrected October 17).** Đoković (arXiv:1008.2043; Int. Math. Forum 6 (2011)) claimed 1527 = 3·509,
  2545 = 5·509 and 2571 = 3·857 through Cor. 29(i), using WT(509) from good matrices of order 127, and H(428). The
  printed proof fails for every Williamson-type input and for every non-symmetric H (with a symmetric H it works,
  but that case is part (ii)). With U₄ = 0 and V₁ = I, X₁X₄ᵀ − X₄X₁ᵀ = diag(0, 2(V₄ᵀ − V₄) ⊗ [[1, −1], [−1,
  1]]), and the printed V₄ = [[W₃, W₄], [−W₄, W₃]] has V₄ᵀ − V₄ with off-diagonal block −(W₄ + W₄ᵀ), whose diagonal
  is ±2 for every ±1 matrix W₄. With V₃ = V₄ = H, X₁ and X₃ are amicable only if H = Hᵀ, which is part (ii). Both
  papers define Williamson-type matrices as pairwise amicable (JCMCC 7 p. 98, (4); 1992 Definition 3.3), and Lemma 26
  needs it (its condition (ii)). The H variant is Seberry–Yamada's own ("new", their Remark 30); the Williamson-type
  variant is credited to Miyamoto (JCTA 57, 1991), which we could not read (ScienceDirect, the Elsevier API, CORE and
  the Internet Archive all failed). We found no valid proof, erratum or other construction, in Đoković's later
  papers, Cati–Pasechnik, SageMath, Zaitsev's page or KMTZ. FRONTIER.md therefore lists 1527, 2545 and 2571 as
  unresolved, together with 2427 = 3·809 and 2611 = 7·373, which Seberry–Yamada's Table A.2 (code o2) builds from
  Williamson-type matrices of orders 809 and 373 that their Table A.1 obtains only from this corollary. We checked
  this with code written from the printed pages, against positive controls (Corollary 29(ii) with symmetric
  conference matrices, part (i) with a symmetric H) that do give Williamson-type matrices.
- **Beyond Table 4.** Đoković's Propositions 3.2–3.4 rest on the corollary alone for 13 entries between 3001 and 9999:
  4285, 5999, 6617, 7179, 7323, 8079, 8187, 8653, 8751, 8859, 9427, 9627 (printed as 9087) and 9671. 8079 is now reached
  by the Paley-block route ([hadamard/32316](hadamard/32316)), 7179 and 7323 by Miyamoto blocks (§23), 6617, 8653 and
  9671 with near-Williamson matrices (§25), and 8187, 8751, 8859 and 9627 with symmetric Williamson-type matrices (§26);
  4285, 5999 and 9427 (q = 857) are unresolved as far as we know. His other uses of the corollary have valid
  alternatives, checked for every q: Corollary 29(ii) with a symmetric H((q − 1)/2) or a symmetric conference matrix, or
  Turyn (e.g. 7167 through Paley's conference matrix of order 1194, 9111 through Turyn with 6073 prime).
- **Through Table A.1 (October 17).** The failure also reaches orders built by Cooper–Wallis (or Yang's products) from
  the Williamson-type matrices of Table A.1 that rest on the corollary alone: the 19 primes q = 173, 293, 373, 521, 613,
  757, 773, 809, 953, 1109, 1193, 1301, 1453, 1493, 1637, 1693, 1733, 1753 and 1993 (all of key wl, credited to
  Miyamoto's paper, which we could not read), for which Turyn, Corollary 8.8(ii) with a symmetric input known before
  2026, and Corollaries 8.7 and 8.9 all fail. (An earlier version of this list left out 173, whose wl entry comes from
  Williamson matrices of order 43; it would need a symmetric conference matrix of order 86 for part (ii), and Verma's of
  arXiv:2610.08893 (6 October 2026) now provides one; see the last bullet.) The orders n ≤ 9999 that a table records as
  known and that rest only on these primes, after checking Paley I and II, Miyamoto, Yamada–Kiyasu (with the
  skew-Hadamard matrices known to the Đoković–Seberry census or built here), Spence, the other T × WT splits and
  SageMath's development version (Unknown for each), are:
  - n ≤ 2999, omitted by Cati–Pasechnik's Table 4 (m = 2 from Seberry–Yamada's Table A.17; code o2 in the 1992 Table
    A.2): 519, 865, 1211, 1557, 1903 and 2249 (from WT(173)), 1465 and 2051 (from WT(293)), which the Paley-block route
    reaches (q = 691 and 1171) and which are built here, 2427 (above), built here by Miyamoto blocks with a periodic
    Golay pair (§23), and 2611 (above), which no route we know reaches;
  - Handbook Table 1.53 with t = 2: 3223, 3633, 3809 and 3979 (Paley blocks; 3223 also by Yamada–Kiyasu from our skew
    matrix of order 1612), 3065 (Miyamoto blocks with Williamson matrices of order 153, §23), 3903 (the same with order
    325, not built) and 3579 and 3865 (Miyamoto blocks with symmetric Williamson-type matrices, §26);
  - Đoković's set Δ (arXiv:0912.5091), whose proof takes Williamson-type matrices from Table A.1: 4981, 5017, 5363,
    5567, 5709, 7093, 7325, 7439, 8131, 8477, 9083, 9169, 9325 and 9515 (Paley blocks), 8899 and 6743 (Miyamoto blocks
    with a periodic Golay pair or Williamson matrices of order 153, §23), 9107 (the same with order 325, not built),
    5259 (below), 6341 and 8579 (near-Williamson matrices, §25), 4765, 5199, 5411, 5545, 5965, 5979, 8465, 8503, 8665
    and 9965 (symmetric Williamson-type matrices, §26), and 4479, 6813 and 8327 (no route).
  - 5259 = 3 · 1753 is now settled: 1753 ≡ 1 (mod 8), and Corollary 8.8(ii) with the symmetric Hadamard matrix of
    order 876 of this repository gives symmetric Williamson-type matrices of order 1753 (checked exactly), hence
    H(21036) by Cooper–Wallis ([hadamard/21036](hadamard/21036)). 8765 = 5 · 1753 has Paley I (35059 is prime).
  - The orders built here are labelled First explicit construction: a table records them as known, and the only
    construction we could trace is the corollary's. If Miyamoto's 1991 paper proves the Williamson-type case, they
    were known.

- **Repaired primes (October 17, later).** The list of 19 was too long. The corollary works when the Williamson-type
  matrices of order (q − 1)/4 have a repeated block (above), and Seberry–Yamada's own symmetric constructions give
  repeated blocks: Corollary 8.8(ii) with M₃ = M₄ = H, Corollary 27 (U₃ = U₄, V₃ = V₄) and Corollary 31 all have X₃ =
  X₄. Their Remark 30 derives Williamson-type matrices of orders 293, 773, 1453, 1637, 1693 and 1733 from orders 73,
  193, 363, 409, 423 and 433, and those are exactly of these kinds (73, 193, 409 and 433 by Corollary 8.8(ii) with
  symmetric Hadamard matrices of orders 36, 96, 204 and 216; 363 by Corollary 27 (their Remark 28) or Corollary 31 with
  q = 181 over GF(361), which we used; 423 by Corollary 27(i)). The repeated-block form of Lemma 26 then gives symmetric
  Williamson-type matrices of all six orders, checked exactly (`python3 search/symwt.py "cor29r(773, cor29h(193,
  kron(paley2(5), h2(), h2(), h2())))"` and so on), and also of order 2917 from 729 (Đoković's Proposition 3.3). So
  Table A.1's entries for these primes are correct, with a proof we found in none of the sources we could read
  (Seberry–Balonin, arXiv:1512.01732, Lemmas 1, 3 and 4, observe that these corollaries give propus-type matrices;
  Miyamoto's paper may contain the arrangement through its symmetric block form), and the orders resting on them were
  obtainable classically: 1465, 2051, 3223, 3809, 4981, 5567, 7325 and 9083 (q = 293, built by Paley blocks) and 3865,
  5411, 8503, 5199, 8665, 8465 and 8751 (§26). The complex order 1546 follows at once: the array [[X, Y], [−Ȳ, X̄]] of
  the symmetric Williamson-type matrices of order 773 is a complex Hadamard matrix (checked). For 173, Verma's symmetric
  conference matrix of order 86 (arXiv:2610.08893, 6 October 2026; we checked C = Cᵀ, zero diagonal and C² = 85I from
  the printed 12 × 12 array of circulants) gives symmetric Williamson-type matrices of order 173 by Corollary 8.8(ii),
  so the orders 519, 865, 1211, 1557, 1903, 2249, 3633 and 3979 and the other 173 · t entries built by Paley blocks have
  a published route since then, which no one has stated.
- **Kronecker products (audit of October 17).** Lemma 26 asks only for a symmetric quaternary Hadamard matrix
  (M₃ + M₄)/2 + i(M₄ − M₃)/2 of order (q − 1)/2 (§19). Paley's symmetric conference matrices give one (C ± I),
  symmetric Hadamard matrices give (H, H), the repeated-block arrangement gives one from (Y, Z, X, X), and Kronecker
  products of these are again of this kind (`symwt.py "cor29q(q, …)"`). This repairs 521 (260 = 10 · 26, conference
  matrices of orders 10 and 26), 953 (476 = 14 · 34, with Williamson-type matrices of order 17 from Corollary 8.8(ii)
  with a symmetric H(8)) and 1993 (996 = 6 · 166, Corollary 31 with p = 41), and Đoković's 2729, 2953 and 3209 (§26),
  all checked exactly. The primes of the list still without a route are 373, 613, 757, 809, 1109, 1193, 1301 and 1493
  (1753 was settled by H(876) above). For six of them (q − 1)/2 = 186, 306, 378, 554, 650 and 746 is twice an odd
  number, so it has no factorisation into even orders (a quaternary Hadamard matrix of order > 1 has even order), and
  no conference matrix (185, 305, 377, 553, 649 and 745 are not prime powers) or repeated-block Williamson-type matrix
  of half that order is known. 404 = 4 · 101 and 596 = 4 · 149 would need symmetric Hadamard matrices of those orders,
  which the census lists as undecided, or repeated-block Williamson-type matrices of orders 101 or 149.
  A search for circulant Williamson matrices (A, B, C, C) of orders 43, 93, 153, 189, 277, 325 and 373 over the
  multiplier groups containing −1 with at most 18 orbits (some row-sum classes skipped as too large) found none.

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

**Order 214 (N = 107) is now the smallest open complex Hadamard order** in our accounting (constructed in §16). It is
blank in Table 11.2. No near-Williamson or quaternary Legendre pair of the matching size (ℓ = 106) is known to us, and
107 is a safe prime:
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

## 15. Compressions decide which row patterns live: the complex orders 342, 366, 430 and 454 (October 15)

The Gaussian estimate of §12 treats the class sums of a structure as independent continuous variables. For d | N
the fibre sums x(a) = Σ_{g ≡ a (mod d)} X(g) obey Σ_X PAF_{x_X} = 4N δ on ℤ_d exactly, and when H fixes ℤ_d
pointwise each fibre sum is ≡ ±1 (mod 2|H|). The compressed equation then has very few solutions: a row pattern
either has none, or all its solutions sit on a handful of compression tuples. Exact counts at N = 93 (ℤ₃ × ℤ₃₁,
|H| = 5): the three patterns the ℤ₃ compressions rule out have 0 solutions where the Gaussian predicted 10^2.1–10^2.6,
and the five live ones have 1,296–21,384, 12–60 times the prediction. Conditioning on the compressions (per tuple,
a lattice Gaussian of the class sums from blocks sampled with those compressions) gets the zeros exactly and the
live counts within a factor 1.2–7; bordered at N = 83 (ℤ₂ × ℤ₄₁) it predicts 10^5.2 for 17,536, where the plain
estimate was 10^0.8 ([search/NOTES.md](search/NOTES.md) §42).
- **The pattern 15, 15; 3, 15 of N = 171 is empty** under ⟨68⟩ and ⟨11⟩ (the ℤ₉ compression, a two-square argument);
  the low-energy runs on it, rated 10^2 by the plain estimate, could not succeed.
- **N = 183 (order 366), ⟨34⟩ of order 5, rows 1, 23; 9, 11.** The ℤ₃ compressions leave one class (A ≡ (−1, 1, 1),
  B ≡ (1, 11, 11), C ≡ (−1, −1, 11), D ≡ (−11, 11, 11), up to translations and reversals of C and D). The earlier
  low-energy lists, chosen without it, held almost no blocks of this class and gave nothing in 5,900 split values.
  With lists drawn from the class (pools 18,480, 15,840, 6.8 · 10⁸, 5.0 · 10⁸; lists of 2 · 10⁷, 10^2.8 predicted
  inside) kt4 found 8 solutions in 22 minutes: the complex Hadamard matrix of order 366
  ([C366.json](complex-hadamard/C366.json)), blank in the 1992 table, 365 = 5 · 73 not a prime power.
- **N = 227 (order 454), bordered.** 227 is a safe prime: ℤ₂₂₇ has no usable multiplier group, and the bordered
  array lives on ℤ₂₂₆ = ℤ₂ × ℤ₁₁₃ with ⟨49⟩ of order 7. The plain estimate was 10^0.6; conditioned on the ℤ₂
  compressions it is 10^6.0. One compression class (A ≡ (1, −1), B ≡ (1, −1), C ≡ (1, 1), D ≡ (−15, 15)) has pools
  4900, 4900, 1.7 · 10⁸ and 1.3 · 10⁸; low-energy lists of 2 · 10⁷ (10^3.8 solutions predicted inside) gave 64
  solutions in six minutes: the complex Hadamard matrix of order 454 ([C454.json](complex-hadamard/C454.json)), blank
  in the 1992 table, 453 = 3 · 151 not a prime power, and not in any family we checked.
- **Every product structure of the open orders, conditioned** (`search/comp_scan.py plain`: N = d · q < 1000 with
  d = 3, …, 15, H ≤ ℤ_q* of odd order fixing ℤ_d, 15–60 orbits): 100 structures. Of the 68 whose compression space
  can be listed, 62 have no row pattern with a feasible compression tuple, so they are empty; 171, 207 and 921 have
  feasible tuples but none with an estimate; three are live: N = 215 (|H| = 7, conditioned 10^5.0), 453 (|H| = 15,
  10^2.1) and 633 (|H| = 15, 10^3.0). Of the 32 with too many compressions to list (ℤ₁₁ with six values per fibre,
  …), 24 are empty already by the value at 0 of the compressed equation together with the row sums, Σ x² = 4N with
  every fibre sum x ≡ ±1 (mod 2|H|): for N = 209 = 11 · 19, |H| = 9, x = ±1 + 18k with k ∈ {−1, 0, 1}, so
  44 + 288a + 360b = 836, i.e. 4a + 5b = 11, which has no solution in non-negative integers.
- **N = 215 (order 430), ⟨41⟩ of order 7, rows 15, 15; 11, 17.** Over ℤ₅ the pattern has 450 feasible compression
  tuples in 3 classes up to symmetry, each predicting about 10^2.3 solutions with pools of only 1200, 1200, 2.4 · 10⁶
  and 2.4 · 10⁶. The exact four-list join of the complete lists of the first class (A ≡ (15, 1, −1), B ≡ (15, −1, 1),
  C ≡ (15, −1, −1, −1, −1), D ≡ (1, 13, 1, 1, 1)) found 5 solutions in nine minutes: the complex Hadamard matrix of
  order 430 ([C430.json](complex-hadamard/C430.json)). Blank in the 1992 table; 429 = 3 · 11 · 13 is neither a prime
  power nor a sum of two squares (so no conference matrix of order 430 exists); ℓ = 214 lies in neither infinite
  quaternary Legendre pair family.
- **Every multiplier group, every quotient** (`comp_scan.py plainall`): H ranges over all subgroups of ℤ_N* without −1,
  and the estimate of a pattern is its smallest over all listable quotients ℤ_d, including those on which H acts. A
  coarse quotient alone can mislead: for N = 171 and H = ⟨31⟩ (order 3 on ℤ₉) the best ℤ₃ tuple, rated 10^2.0, is dead,
  since the class sums of three classes never vanish together (0 joint hits where independence predicts 73); on ℤ₉ the
  pattern keeps 10^0.5 per tuple, and two tuples searched exhaustively gave nothing.
- **N = 171 (order 342), H = ⟨11⟩ of order 6, rows 11, 11; 9, 19.** H permutes the units of ℤ₉ transitively and has
  order 3 on ℤ₁₉, and −1 ∉ H. On ℤ₉ (fibre classes {0}, {3, 6}, ℤ₉*) the pattern has 13 feasible compression tuples,
  conditioned 10^4.5 in all. The first, A ≡ (−5, 1, 5), B ≡ (−17, 3, 5), C ≡ (1, −1, 7), D ≡ (−1, 1, 7), has pools 1134,
  378, 2.8 · 10⁷ and 2.8 · 10⁷ and 10^3.5 predicted inside lists of 2 · 10⁷; the exact sumset of the merge classes gave
  49,652 live split values, and kt4 found 12 solutions in 263 s: the complex Hadamard matrix of order 342
  ([C342.json](complex-hadamard/C342.json)), the fourth smallest open order. Blank in the 1992 table; 341 = 11 · 31 is
  neither a prime power nor a sum of two squares; ℓ = 170 lies in neither quaternary Legendre pair family. (The pattern
  15, 15; 3, 15 of the same N is empty under ⟨11⟩ by the two-square argument above; the row pattern matters.)

## 16. The character-bordered Szöllősi array: the complex order 214 (October 15)

**The construction** ([complex-hadamard/TBORDERED.md](complex-hadamard/TBORDERED.md), `search/tbsz.py`). The bordered
Goethals–Seidel arrays with a character set T ([search/NOTES.md](search/NOTES.md) §18) carry over to Szöllősi's array.
G abelian of order m = N − t, φ: G → Q onto a group of order q, T a conjugation-closed set of t characters of Q, and
four ±1 sequences on G (A, B symmetric) with Σ_k PAF_k(g) = −4 Σ_{χ∈T} χ(φ(g)) for g ≠ 0. The border is three
quaternary matrices: R̂ (2t × 2q) and Ĉ (2q × 2t), whose rows and columns are orthogonal bases of the space V_T of
pairs of functions on Q with Fourier support in T, and Xb = −R̂ Mc* Ĉ / 2q, which must be quaternary (Mc is Szöllősi's
array over Q of the compressions). Then [[Xb, R], [R′, Kc]] is a complex Hadamard matrix of order 2N; t = 1 is the
bordered array of §9.
- Instances: Q = ℤ_q with T all characters (q | N); **Q = ℤ₄ with T = {1, χ, χ̄}** and Q = ℤ₂² with T the three
  nontrivial characters or {1, χ₁, χ₂}, all with t = 3, which need N ≡ 3 (mod 4). (The two choices on ℤ₂² are one
  family, see §17.)
- The negative results of §9 do not apply: they concern a 3-point border for the Williamson-type array, and for the
  complex array a border made of single characters. Here the border vectors are arbitrary quaternary vectors in
  V_T. With Q = ℤ₄, V_T has 36 quaternary vectors per half, so the border problem is a small depth-first search
  (1296 candidates, 324 after normalisation), solved for every core we tried.
- Checks: every core on ℤ₈, ℤ₂ × ℤ₄ (N = 11) and ℤ₁₂ (N = 15), samples completed and verified (orders 22 and 30);
  `python3 search/tbsz.py test`. For Q = ℤ₃ with T all characters, all 192 compression quadruples of the ℤ₃ problem
  have a border, and complete examples of orders 18 and 30 were verified.

**N = 107 (order 214), the smallest open order.** On G = ℤ₁₀₄ = ℤ₈ × ℤ₁₃ with φ = x mod 4 and H = ⟨(1, 3)⟩ of
order 3 (40 orbits, 20 classes), the pattern 2, 0; 2, 2 has 128 feasible ℤ₄ compression tuples with a positive
conditioned estimate, 10^5.7 in all; its plain Gaussian estimate is 10^4.6. Under H = ⟨(1, 4)⟩ of order 6 and
⟨(1, 2)⟩ of order 12 both patterns are dead on ℤ₄. The first tuple (pools 5040, 2520, 1.97 · 10⁷, 1.92 · 10⁸; 10^3.8
predicted) was searched with A, B, C complete and D cut to 2.9 · 10⁷ (retention 0.5): kt4 found 5 solutions in the
first 256 of 52,762 exact split values, six minutes in all. Each completes to a complex Hadamard matrix of order 214
([C214.json](complex-hadamard/C214.json)). Blank in the 1992 table; 213 = 3 · 71 is not a prime power, and ℓ = 106
lies in neither quaternary Legendre pair family (213 and 105 are not prime powers); the quaternary Legendre pairs of
2026 have lengths at most 100.

**The same family for the other primes N ≡ 3 (mod 4)** (`tbsz.py rank`, H generated by one automorphism: diag(a, b)
on ℤ_(2^k) × ℤ_o, or an element of GL(2, 2) times a unit on ℤ₂² × ℤ_o; 12 to 70 orbits):
- N = 167: every such structure on ℤ₁₆₄ = ℤ₄ × ℤ₄₁ (Q = ℤ₄, 8 groups) and on ℤ₂² × ℤ₄₁ (Q = ℤ₂², both choices of T)
  is empty on Q. The fibres have 41 points, and the orbits of H pin their sums to coarse progressions.
- N = 179 on ℤ₁₆ × ℤ₁₁ (Q = ℤ₄): only ⟨(1, 3)⟩ of order 5 lives, at 10^−1.9 and 10^−2.7. For the non-cyclic groups
  of order 176 see §17.
- N = 239: 59 is a safe prime, so ℤ₄ × ℤ₅₉ has only H of order 29 in range (12 orbits, empty on ℤ₄), and ℤ₂² × ℤ₅₉
  none.

## 17. The ℤ₆ character-bordered array: the complex order 322 (October 15)

**Two reductions** ([complex-hadamard/TBORDERED.md](complex-hadamard/TBORDERED.md)). Multiplying the four sequences by
λ ∘ φ, λ a real character of Q, turns a core for T into one for λT, and the border carries over (K becomes D K D with
D diagonal ±1). And t must be odd: for 1 ∈ T and t even the row sums are odd and Σ r² ≡ 4 (mod 8) ≠ 4t; for 1 ∉ T the
border rows have zero sums on each half, which forces q, hence m, even. A conjugation-closed T of odd size contains a
real character, so 1 ∈ T may be assumed. So Q = ℤ₄ with T = {χ, χ², χ³}, the two choices on ℤ₂², and the two on ℤ₆
are not new families; for t = 3 the families are {1, χ₁, χ₂} on ℤ₂² and {1, χ, χ̄} with χ of order 3, 4 or 6.

**The border for Q = ℤ₆.** With T = all characters but the sign (t = 5), the family reaches every N ≡ 5 (mod 6),
among them 161, 167 and 179. §45 of [search/NOTES.md](search/NOTES.md) had parked it after a randomised border search
(orthogonal bases of the form (v_i, ±c w_i)) found no border. An exhaustive search (`search/border_exh.c`) decides
the question. It takes the 40,000 normalised quaternary vectors of V_T ⊕ V_T, builds bit matrices for
compatibility (u Mc* c / 2q quaternary) and orthogonality, runs rows depth-first while keeping the common
compatible columns, and ends with a 10-clique of columns. It found a border in under a second for each of the ℤ₆
compression tuples of N = 161 we gave it, and the solutions are not of the paired form. The earlier negative was the
solver's.

**N = 161 (order 322), then taken for the smallest open order.** (It is not open: 643 = 4 · 161 − 1 is prime, so
Ito and Egan 2019 give it; see §18.) The plain array over ℤ₁₆₁ (multiplier groups with at most 64
orbits) and the bordered array over ℤ₁₆₀ (33 to 50 orbits) had nothing alive. With Q = ℤ₆:
- cyclic G = ℤ₁₅₆: the best group, ⟨29⟩ of order 6 (40 orbits), rates 10^0.34 conditioned on ℤ₆;
- G = ℤ₂ × ℤ₂ × ℤ₃ × ℤ₁₃, φ = (x₀, x₂ mod 3), with H ranging over the cyclic subgroups generated by K-preserving
  automorphisms (shears on ℤ₂², units on ℤ₃ and ℤ₁₃): H = ⟨M⟩, M = (x₀, x₁, x₂, x₃) ↦ (x₀, x₀ + x₁, 2x₂, 3x₃), of
  order 6 (35 and 20 orbits, 19 classes), rates 10^1.1, 10^1.9 and 10^1.9 for the patterns 0, 0; 4, 2 and 2, 0;
  4, 0 and 4, 0; 2, 0 (plain Gaussian 10^−3.7). Its pools are small: for 2, 0; 4, 0 the best tuple has complete
  lists of 864, 432, 8.2 · 10⁶ and 1.3 · 10⁶ sequences (10^1.8 predicted).
The exact join of that tuple (1367 live split values on five merge classes) found 1024 solutions at the 246th split
value, 80 seconds into the join and six minutes after the start. Every solution checked completes to a complex
Hadamard matrix of order 322 ([C322.json](complex-hadamard/C322.json)). The 1992 table is blank at 161. ℓ = 160 lies
in neither quaternary Legendre pair family (321 = 3 · 107 and 159 = 3 · 53 are not prime powers). There is no
conference matrix of order 322 (321 is not a sum of two squares).

**167 and 179 (orders 334 and 358).** Their ℤ₆ fibres are odd (27 and 29 points). The ℤ₆ compression equation itself
has solutions (59 and 82 PAF-class combinations, A and B symmetric, no H), but the multiplier groups tried are empty
on ℤ₆: on ℤ₁₆₂, whose unit group is cyclic, so that |H| is odd, both groups with 12 to 46 orbits (|H| = 27 and 9);
on ℤ₁₇₄ all five with 12 to 46 orbits (|H| = 4, 7, 14, 14 and 28); on ℤ₂ × ℤ₉ × ℤ₉ the first one tried, of order 9
(34 orbits). For 179, t = 3 on ℤ₄ × ℤ₄ × ℤ₁₁ and ℤ₂ × ℤ₂ × ℤ₄ × ℤ₁₁ (Q = ℤ₄ or ℤ₂², H generated by one automorphism that
preserves K, 30 to 70 orbits) rates at most 10^−0.45: ℤ₄ × ℤ₄ × ℤ₁₁ over ℤ₂² with H of order 10, rows 2, 2; 2, 0.
Its single feasible tuple has complete lists of 192, 192, 1.2 · 10⁴ and 7.3 · 10⁵ sequences, and the merge-class
sumsets of the two sides of the join do not meet: no solution.


## 18. Paley's matrix halves: 76 known complex orders (October 16)

**How it came up.** The searches for 406, 414 and 418 had stalled: their best structures lose almost all their
estimated mass at finer conditioning levels. So we looked for product and plug-in constructions instead, and started
from the equivalence behind Turyn's theorem. A complex Hadamard matrix of order 2N is a real Hadamard matrix of order
4N with a *complex structure*: a monomial automorphism (P, Q) with P² = Q² = −I. In a basis adapted to P such a
matrix is [[A, B], [−B, A]], and (A + iB)/(1 + i) is quaternary (Seberry–Yamada 1992, Lemma 11.16).

Paley's skew Hadamard matrix of order q + 1 (q ≡ 3 mod 4) admits SL(2, q) acting monomially, P_g e_u = χ(λ_u) e_(g u)
on the projective line. The centre acts as χ(−1) I = −I. So any g of order 4, for example u ↦ −1/u, gives
P_g² = −I, and P_g commutes with H. That is a complex Hadamard matrix of order (q + 1)/2 for every prime power
q ≡ 3 (mod 4) ([complex-hadamard/PALEY.md](complex-hadamard/PALEY.md)).

**Checks.** `search/paley_complex.py` builds K directly from H, with GF(p^k) by log tables, and checks K K* = m I
exactly. `verify/verify_complex.py` rebuilds K with its own polynomial arithmetic (squares and inverses found
directly). Both give the same SHA-256 digests. Tested on q = 7 … 643, including 27 = 3³ and 243 = 3⁵, and certified
for the 63 open orders.

**Counts below N = 1000.** There are 145 odd N with 4N − 1 a prime power. The 1992 table knows 66 of them and leaves 79
blank. Of the blanks:
- 70 came from Đoković's good matrices (1993), and 442 from a Mathon conference matrix;
- 13 had been constructed here by searches: 142, 154, 222, 262, 286, 310, 322, 330, 342, 370, 430, 454, 742;
- the other 63 (406, 414, 474, 526, 546, 582, 646, 666 (q = 11³), …, 1966) now have certificates from the theorem.

All 76 are known, from Ito and Egan 2019 (below). That leaves 143 open N < 1000 (115 after §19 and §20). None of the
smallest open orders (334, 358, 418, 438, 446, 478, 498, …) has 4N − 1 a prime power. The 11 orders found by search
that do not follow (214, 238, 302, 326, 366, 382, 426, 494, 518, 534, 754) were found concurrently by Zaitsev et al.
(§20).

**Novelty: known.** The first version of this section said that no source combines the automorphism group of
Paley's matrix (Kantor 1969) with the half-order lemma, and labelled the 63 orders New. Searching the literature on
Ito's Hadamard groups, which we had flagged as unread, turned up the source. Ito showed that Paley's matrix makes the
dicyclic group Q₈ₜ, 4t = q + 1, a Hadamard group: it has a central relative (4t, 2, 4t, 2t) difference set. Schmidt
(Des. Codes Cryptogr. 17 (1999), Corollary 3.6) extends this to t = 2^a 10^b 26^c m. And Egan (Des. Codes Cryptogr.
87 (2019) 67–74, Theorem 4.5) turns every such difference set into a BH(2t, 4); his Remark 4.6 builds BH(70, 4) this
way from q = 139. Our element u ↦ −1/u is the element of order 4 of the dicyclic group. So the theorem is Egan's, and
it covers 13 orders this repository had called New (142, 154 and 322 among them) as well as the 63. Neither 2026
paper cites Egan, and both take their smallest open order from the 1992 table, as we did. Ito's conjecture (Q₈ₜ is
a Hadamard group for every t) is verified for t ≤ 46 only (Schmidt), so it adds nothing in the open range.

**Quaternion type, and why it does not multiply.** For q ≡ 3 (mod 8) the Sylow 2-subgroup Q₈ of SL(2, q), generated
by u ↦ −1/u and h = [[a, b], [b, −a]] with a² + b² = −1, acts on the q + 1 points modulo its centre with all orbits of
size 4. This was checked for q = 11 … 83. So Paley's matrix is of quaternion type: a Williamson array of four blocks
of odd order (q + 1)/4 with Σ X Xᵀ = (q + 1) I. In the basis given by orbit representatives no two of the blocks are
amicable (checked for q ≤ 139). So Cooper–Wallis does not apply directly, and these blocks give no new real orders
as they stand.

**Two products that cannot work.**
1. *Complex Agaian halving.* An Agaian-type product M₁ ⊗ E + M₂ ⊗ F needs a complex Hadamard matrix [N₁ N₂] of
   order 2n whose column halves satisfy N₁ ∘ N̄₂ ∈ {±1}. Write N₁ = N₂ ∘ S with S = (s_j) a ±1 matrix of columns.
   Orthogonality of the columns j, j′, k, k′ forces each of the four sets {r : (s_j(r), s_k(r)) = (±1, ±1)} to have
   even size, because a sum of an odd number of fourth roots of unity is not 0. But s_j is balanced with n entries
   +1, and n is odd. So no such matrix exists for odd n > 1, and complex Hadamard matrices of orders 2m and 2n
   give nothing of order 2mn this way.
2. *Cooper–Wallis for Szöllősi quadruples.* Plug Szöllősi blocks (A₁, A₂ symmetric, A₃, A₄ arbitrary) into a 4 × 4
   array of (0, ±1) circulants of odd order t, with the two symmetric output blocks built from A₁ and A₂ only. The
   column weights then force the array to split into two 2 × 2 arrays. In the arbitrary half, X = T₃₃ + T₃₄ and
   Z = T₄₃ + T₄₄ would be a periodic complementary pair of odd length t, which does not exist for t > 1.

**Three side checks.**
- *The 2026 matrices.* E. Piza's record and the Hadamard-B data give the seeds of all twelve Alpöge–Voinov–
  Reynolds-Haertle matrices as bordered Goethals–Seidel arrays over abelian groups:
  - 668, 716, 1676 and 1772 have a border of one point (PAF sum −4);
  - 892, 1132, 1244, 1948 and 1964 have no border;
  - 1388, 1436 and 1916 are character-bordered, over ℤ₂ × ℤ₃ × ℤ₅₇, ℤ₂ × ℤ₂ × ℤ₈ × ℤ₁₁ and ℤ₂ × ℤ₂ × ℤ₁₁₉, with
    t = 5, 7 and 3: PAF sums {−20, 4}, {−28, 4} and {−12, 4}, the family of §18 of search/NOTES.md.

  If two seeds of one order were symmetric (after a translation), Szöllősi's array (plain, bordered or
  character-bordered) would give the complex order 2N: 334, 358, 446, …. None of the 48 seeds is symmetric about any
  centre.
- *The symmetric analogue is classical.* With R: x ↦ −x one has R S Rᵀ = Sᵀ. So D R H, with D = −1 at ∞, is a
  symmetric Hadamard matrix equivalent to Paley's (checked for q = 7 … 43). Consistently, none of the 595 undecided
  v of the symmetric census has 4v − 1 a prime power. The complex table, by contrast, leaves every half order
  blank.
- *A lead not yet built.* Seberry and Yamada's Corollary 8.8 (after Miyamoto) turns Williamson-type matrices of order
  (q − 1)/4, or a Hadamard matrix of order (q − 1)/2, into Williamson-type matrices of order q ≡ 1 (mod 4). The 2026
  inputs make three orders reachable that we then took for open. They are known from Zaitsev's explicit matrices,
  whose real blocks are the printed matrices X₁, …, X₄ of this corollary (X₁ is not amicable with X₃ or X₄), so they
  are instances of the completion of §19, not Williamson-type matrices; Zaitsev attributes 1018 and 1714 to
  Đoković's Williamson-type matrices and 538 to unpublished work (§7, §19):
  538 (q = 269, from the near-Williamson matrices of order 67 of arXiv:2605.08661, made pairwise amicable as
  (AR, B, C, D)), 1018 (q = 509, from Đoković's good matrices of order 127) and 1714 (q = 857, from H(428)). The
  survey does not give the matrices U_i, V_i of its Lemma 8.6. The natural choice (U = (C₁, C₂, 0, 0) from the
  quadratic-residue blocks of the conference matrix, V = (I, 0, K₁, K₂)) needs a symmetric complex Hadamard matrix of
  order (q − 1)/2, and 2 × 2 block patterns built from four symmetric Williamson-type blocks cannot give one. So these
  three orders are not claimed.

  Two follow-ups.
  - Seberry and Yamada's 1989 Kobe conference version of the M-structure paper (Wollongong repository,
    infopapers/1044) states the lemma (its Lemma 13) and the corollary (its Corollary 16(i), credited to Miyamoto)
    without proof. The JCMCC 7 (1990) 97–137 version numbers them Lemma 26 and Corollary 29 and prints the proof
    examined in §7. Its Corollary 36 (Corollary 23 of the Kobe version) turns a complex Hadamard matrix of order
    (q − 1)/2 into a *real* matrix of order 4q only. With the conference blocks split along the quartic residues, the
    amicability of the V's forces their 2 × 2 patterns into {Z, X} with Z and X not amicable, so that ansatz fails
    too.
  - Miyamoto's real array (verify.py) does not carry a complex structure over from K. For q = 13 and K of order 12 of
    the form [[A, B], [−B, A]], the matrix of order 52 has an automorphism group of order 8 (nauty), and no element
    squares to −I.

  Also the nauty test of the 2026 matrix of order 668 and of our symmetric 876 for a complex structure (which would
  give 334 and 438) did not finish: the Hadamard graphs of these orders ran for two hours without a group.


## 19. Miyamoto's lemma with a symmetric complex input (October 16)

**The reduction.** Let q ≡ 1 (mod 4) be a prime, n = (q − 1)/2, g a primitive root, χ the quadratic character. Take,
in Seberry and Yamada's Lemma 8.6 (after Miyamoto),

    U_1 = −[χ(g^(2i) − g^(2j))],   U_2 = [χ(g^(2i) − g^(1 − 2j))],   U_3 = U_4 = 0,
    V_1 = I,   V_2 = 0,   V_3 = M_3,   V_4 = M_4,

with M_3, M_4 symmetric ±1 matrices of order n, M_3 M_4 = M_4 M_3 and M_3² + M_4² = 2n I. In the exponents U_1 is a
symmetric circulant and U_2 a back-circulant, so they are amicable; U_1² + U_2 U_2ᵀ = q I − 2J because Paley's core
satisfies Q² = q I − J; the row sums are 1 and 0. The lemma then gives Williamson-type matrices of order q, and these
give a complex Hadamard matrix of order 2q. A pair (M_3, M_4) is the same thing as a symmetric quaternary complex
Hadamard matrix (M_3 + i M_4)/(1 + i) of order n. The inputs of Seberry and Yamada's Corollary 8.8(ii) are special
cases: a symmetric conference matrix W (M_3, M_4 = W ± I) or a symmetric Hadamard matrix (M_3 = M_4).

A second step supplies such pairs: if A, B are symmetric ±1 circulants of order m and W is any ±1 circulant with
A Aᵀ + B Bᵀ + 2 W Wᵀ = 4m I (a propus quadruple (A, W, W, B) with both outer blocks symmetric), then
M_3 = [[B, W], [Wᵀ, −A]] and M_4 = [[A, W], [Wᵀ, −B]] form a pair of order 2m. With q = 4m + 1 prime this gives a
complex Hadamard matrix of order 2q = 8m + 2.

`search/miyamoto_wt.py selftest` builds all of this and checks it exactly for q = 13, 17, 29, 37 and 61, with inputs
from Paley's symmetric Legendre pairs, from Sylvester's matrix of order 8, and from doubly symmetric propus quadruples
of orders 3, 7 and 9 found by brute force.

**What it would give.** The complex orders 2q with q ≡ 1 (mod 4) prime that we had listed as open are 538, 1018, 1306,
1466, 1706, 1714 and 1882 (q = 269, 509, 653, 733, 853, 857, 941). Williamson-type matrices of order q would need a
symmetric quaternary Hadamard matrix of order 134, 254, 326, 366, 426, 428 and 470. The complex order 2q needs only some
quaternary Hadamard matrix of that order (the completion below), and all seven exist: 134 (Kharaghani, Mohammadian,
Tayfeh-Rezaie and Zaitsev, 2026), 254 (Đoković's good matrices of order 127), 326, 366 and 426 (this repository), H(428)
and 470 (Zaitsev's page). So the completion gives these seven complex orders on its own, which were already known.
Đoković (arXiv:1008.2043, Propositions 3.2–3.4) applies Corollary 8.8(i) to get Williamson-type matrices of orders 509
(from good matrices of order 127) and 857 (from H(428)), and Corollary 8.9 (JCMCC 7 Corollary 31) gives 883 (q = 441);
the near-Williamson matrices of order 67 (2026) would give 269 the same way. Corollary 8.8(i) fails as printed (§7), but
538, 1018 and 1714 are known from explicit matrices on Zaitsev's page (checked exactly). The page attributes 1018 and
1714 to Đoković's Williamson-type matrices and gives no further construction. Their real blocks are the printed X₁, …, X₄ of
Lemma 26: X₁ and X₂ are symmetric, and X₁ is not amicable with X₃ or X₄. So they are not Williamson-type matrices but
instances of the completion below, and Williamson-type matrices of orders 269, 509 and 857 remain unestablished. 1766
is known from QH1766, which splits into explicit Williamson-type matrices of order 883. The other four, 1306, 1466, 1706
and 1882, are on Zaitsev's page (§20), so none of this form is open.

**Why we found none.**
- *Doubly symmetric propus over ℤ_m* (m = 67, 127, 163, 183, 213, 235; `search/propus_rank.py`, default structure):
  for m = 67 the multiplier groups of orders 3, 11 and 33 are excluded by the row sums. Under |H| = 3, for example,
  the symmetric blocks have row sums ±1 + 6k with k odd, W has ±1 + 6k for any k, and no such values satisfy
  a² + b² + 2w² = 268. H = {±1} rates 10^−12, and H = 1 rates 10^0.6 in a space of 2^135. The best structures for
  the other m rate 10^−1.3 (163, H = ⟨38⟩ of order 9, which we exhausted: no solution), 10^−2.4 (127), 10^−3.8 (213)
  and 10^−3.9 (183). The Szöllősi solutions we certified for 326, 366 and 426 do not have C = D, not even up to
  shift, sign and reversal.
- *m = 93 (October 17).* With Yamada's U₁, U₂, Lemma 26 needs exactly a symmetric commuting pair (M₃, M₄) of order
  (q − 1)/2 with M₃² + M₄² = (q − 1)I, i.e. a symmetric quaternary Hadamard matrix; a doubly symmetric propus
  quadruple of order m gives one of order 2m. For real orders this is the open input of four FRONTIER entries:
  2611 = 7 · 373 needs one of order 186 (m = 93), 1527 and 2545 one of order 254 (m = 127), 2427 and 2571 one of
  order 404 or 428 (or a symmetric Hadamard matrix of that order, v = 101 or 107 of the symmetric census), and 1883
  = 7 · 269 one of order 134. For m = 93 every multiplier subgroup of ℤ₉₃* of order 3 to 30 was searched
  exhaustively, A and D symmetric and B free or symmetric (`search/propus.py 93 g`, and `search/propus_chain.py 93
  25 --symD` for the 33-orbit structure): no solution. The best of them rated 10^−0.8 (H = ⟨4⟩); H = ⟨61⟩ and
  ⟨32⟩ rate 10^−8.8 and 10^−6.1 and H = 1 rates 10^−2.2 in a space of 2^140, beyond exhaustive search.
- *Bordered symmetric pairs.* A border of one point turns the problem into a symmetric Legendre pair of length
  n − 1 (PAF_a + PAF_b = −2). An exhaustive search of cyclic lengths up to 33 finds them only for 3 and the primes
  ≡ 1 (mod 4), where they come from Paley's matrix. A count of orbits explains this. Symmetric blocks have half the
  free signs of arbitrary ones but face the same number of PAF classes, so the expected number falls exponentially.
- *Symmetric Hadamard of order 428* (v = 107) is itself undecided in the Đoković–Seberry census.

**Corollary 8.8(i) (Williamson-type matrices of order (q − 1)/4 suffice) could not be reconstructed.** Miyamoto's
paper (JCTA 57, 1991) is in Elsevier's open archive, but the full text could not be retrieved from here. The 1992
Table A.1 does use the corollary: for q ≡ 5 (mod 8) it lists Williamson-type matrices of order q exactly when
(q − 1)/4 has them or 2q − 1 is a prime power (checked for every prime q < 1000). We tried the natural forms:
- U-blocks along the quartic ratio classes of Paley's core, and V-blocks taken from four commuting symmetric
  Williamson-type blocks in 2 × 2 patterns. The V-system has 80 solutions (`search/miyamoto_wt.py vpatterns`), and
  each needs a U-block made of a single quartic class in some row. Its row sum is a cyclotomic number that varies with q
  (for q = 269 the four are 6, 5, −7, −5), so the row-sum condition fails.
- Without the 2 × 2 patterns: once U_1 = −Q_SS is forced by the row sums, a pattern that puts the Q_SN part in two
  different U's fails at V_1 V_2ᵀ = V_2 V_1ᵀ. The remaining pattern needs exactly the symmetric pair above.

For q ≡ 1 (mod 8) the octic ratio classes allow more patterns. From σ_(−k) = χ(g^k) σ_k, the combinations with a
row sum that does not depend on q are {1, 7} and {3, 5} with equal signs, {2, 6} with opposite signs, and their
unions. Q_SN can then be split into the halves {1, 7} and {3, 5}, which are 4 × 4 block patterns of order (q − 1)/8.
But U_1 = −Q_SS is still forced, so V_1 is a signed permutation P and the V_j P^T must again be symmetric and
commuting. An arbitrary Hadamard matrix of order (q − 1)/2, which Corollary 8.8(i) allows, does not fit this frame
either.

So the corollary must use U-matrices that are not invariant under the quartic (or octic) subgroup, or a different
lemma. The October 17 audit (§7) settles the question for the printed proof: it fails for every Williamson-type input
and every non-symmetric H. With Yamada's U₁ and U₂, Lemma 26 needs a symmetric quaternary Hadamard matrix of order
(q − 1)/2, and no 2 × 2 block pattern in the Williamson-type blocks (signs, R, transposes) gives one unless two blocks
have equal Gram matrices (exhaustive for orders 3, 7 and 9). If W₁ = W₂ for symmetric commuting Williamson matrices W₁,
…, W₄, then V₁ = I, V₃ = [[W₃, W₁], [W₁, W₄]] and V₄ = [[W₄, −W₁], [−W₁, W₃]] satisfy every hypothesis (checked for q
= 13, 29, 37, 53, 61); this is the repeated-block case of §7.

**A completion that the printed matrices do give (October 17).** Let q ≡ 1 (mod 4) be a prime power, n = (q − 1)/2, and
C any quaternary Hadamard matrix of order n, symmetric or not (a real Hadamard matrix of order n will do). Write
(1 + i)C = M₃ + iM₄ and take U = (A, BR, 0, 0) and V = (I, 0, M₃, M₄) in Lemma 8.6, with A, BR Yamada's matrices as in the
printed proof. The X_i are then ±1 matrices with ΣX_iX_iᵀ = 4qI. They are pairwise amicable only when C is symmetric
(Corollary 8.8(ii)), but in every case X = ½((X₁ + X₂) + i(X₁ − X₂)) and Y = ½((X₃ + X₄) + i(X₃ − X₄)) give
K = [[X, Y], [−Y*, X*]], with entries ±1, ±i and KK* = 2qI.
- *Proof.* X₁ and X₂ are symmetric, and each commutes with X₂, X₃ and X₄: this uses the row and column sums (1, 0, 0,
  0) of the U's, AC = CA, U₃ = U₄ = 0, V₁ = I and V₂ = 0. Hence XY = YX, so the off-diagonal blocks of KK* vanish.
  XX* + YY* = ½ΣX_iX_iᵀ + (i/2)(N₁₂ + N₃₄) with N_ab = X_aX_bᵀ − X_bX_aᵀ. N₁₂ = 0, and N₃₄ = diag(0, 2(M₃M₄ᵀ − M₄M₃ᵀ)
  ⊗ [[1, −1], [−1, 1]]), which is 0 because CC* = nI gives M₃M₄ᵀ = M₄M₃ᵀ. The other diagonal block is
  Y*Y + X*X = ½ΣX_iᵀX_i − (i/2)(N'₁₂ + N'₃₄) with N'_ab = X_aᵀX_b − X_bᵀX_a. ΣX_iᵀX_i = 4qI because A and BR are
  symmetric with column sums 1 and 0 and M₃ᵀM₃ + M₄ᵀM₄ = 2nI, and N'₃₄ = 0 because M₃ᵀM₄ = M₄ᵀM₃; both identities
  come from C*C = nI.
- Checked exactly for q = 13 (a scrambled non-real BH(6, 4)), 17, 41 and 89 (column-scrambled real Hadamard matrices),
  and in the audit for q = 13, 29, 37, 61 and 269. We found it stated in none of the papers we read: JCMCC 7 Corollary
  36 gives only the real order 4q from a complex Hadamard matrix of order (q − 1)/2. Zaitsev's QH538, QH1018 and QH1714
  have exactly its form, with the Williamson-type inputs of the printed proof, but the page does not say so; the theorem
  explains why they are complex Hadamard matrices although their X_i are not Williamson-type. `search/miyamoto_wt.py
  complete` builds it, and `verify/verify_complex.py` (construction `miyamoto-completion`) rebuilds it independently.
- Consequences. The complex orders 746, 1514, 1546, 1618 and 1906 (q = 373, 757, 773, 809, 953), which the 1992 Table
  11.2 records only through Corollary 8.8(i), are known through the completion; without it the open count of N < 1000
  would be 120 instead of 115. 1618 and 1906 are built here from H(404) and H(476)
  ([complex-hadamard](complex-hadamard/README.md)); the other three need complex matrices of orders 186, 378 and 386,
  from Table A.1's Williamson-type matrices of orders 93, 189 and 193, which we did not build. The completion settles
  no open N < 1000 by itself.
- Real orders are unaffected: the completion gives complex 2q, hence real 4q, which Miyamoto's first theorem already
  gives (C doubles to a real H(q − 1)), and it gives no Williamson-type matrices, so the Cooper–Wallis products that
  Đoković built on WT(509) and WT(857) stay unproved.

**The 2026 matrices once more.** A complex structure on a Goethals–Seidel matrix can also come from seeds related by a
multiplier (Y = ±X∘μ up to translation, as in the twisted arrays), not only from symmetric seeds. For the nine
cyclic matrices of Alpöge, Voinov and Reynolds-Haertle (668, 716, 892, 1132, 1244, 1676, 1772, 1948, 1964) we tested
every ordered pair of seeds, every unit μ, translation and sign. Each seed is invariant under the multiplier subgroup
of its search: orders 3, 3, 5, 15, 12, 9 and 7 for the last seven, none for 668 and 716. No seed is ±μ-related to a
different seed, and no seed goes to its negative. So these matrices carry no complex structure of this kind.

## 20. The page of Zaitsev et al.: the 11 complex orders are concurrent (October 16)

V. Zaitsev's page [Quaternary Hadamard matrices](https://www.cs.uleth.ca/~zaitsevv/quaternary.html), with Kharaghani,
Mohammadian and Tayfeh-Rezaie, says "Last updated: 14 October 2026" (HTTP Last-Modified 15 October 2026, 02:50 UTC).
We found it on 16 October. It has matrices of 122 orders 2N; we downloaded every one and checked it exactly
(entries in {±1, ±i}, K K* = 2N I; [search/data/complex_zaitsev_2026.txt](search/data/complex_zaitsev_2026.txt)). It
attributes 23 orders to Szöllősi's array with two symmetric blocks, the construction of §8–14 (142, 154, 222, 238,
262, 286, 302, 310, 326, 342, 366, 370, 382, 426, 430, 494, 518, 534, 686, 742, 754, 814 and 1558), 27 to unpublished
approaches (among them 214, 406, 454, 470, 490, 1066, 1306, 1466, 1706 and 1882), 62 to "a construction utilizing
Paley-type matrices" (our 63 halves of §18 without 406, 814 and 1066, plus 322 and 330), and five, for which it
claims nothing, to Mathon's conference matrices and to Đoković's Williamson-type matrices of orders 509, 857 and 883.

So all 11 orders we had labelled New are on the page, and the page was last modified before we committed 214 and 366.
They are now labelled Concurrent. At ten of them the matrices differ as well: for each pair (a, b) of rows let n_k be
the number of columns c with K[a,c] K̄[b,c] = iᵏ; the multiset of the tuples (n_0, …, n_3), each taken up to rotation
and reversal, is unchanged by row and column permutations, unit scalings and conjugation, and we computed it for K and
Kᵀ. It differs between our matrix and theirs at 214, 302, 326, 366, 382, 426, 494, 518, 534 and 754. At 238 it
agrees, and so does the finer invariant that keeps, for each row, the multiset of its tuples; the two matrices of
order 238 may be equivalent.

**A real order as well.** Turyn's doubling turns a quaternary Hadamard matrix of order 2N into a real one of order 4N,
and we first checked the page only for complex orders. Its matrix of order 1142 = 2 · 571 doubles to H(2284), which we
had committed on 11 October as New ([hadamard/2284](hadamard/2284/README.md#status); `search/quaternary_check.py
--double` checks it). 2284 is now Concurrent too. No other order of the page doubles to an order of FRONTIER.md or of
this repository: for 2092, 2148, 2292, 3732 and 3964 the page lists 1046, 1074, 1146, 1866 and 1982 as unknown, and the
other orders are above its range.

The page also covers 24 orders that we had listed as open, 470, 490, 686 and 1306 among them. With §19 the count of open
N < 1000 falls from 143 to 115 (`search/complex_status.py`); the smallest open orders are 334, 358, 418, 438, 446, 478,
498, 502, 506 and 550, none with 4N − 1 a prime power. The page ends with its own list of the even orders below 2000
with no quaternary Hadamard matrix known; it is the same set of 115. The searches of §8–17 that found nothing at 334,
358, 418, 438 and 446 stand as negative results for those structures.

## 21. Paley's matrix is of quaternion type, not of Williamson type (October 16)

Williamson-type matrices of order 263 would settle four FRONTIER orders (789, 1315, 1841, 2893; §9), and 4 · 263 − 1 =
1051 is a prime. Ito's dicyclic Hadamard group of order 2(q + 1) contains Q₈ = ⟨x ↦ ix, x ↦ c x̄⟩ (N(c) = −1) when
q ≡ 3 (mod 4), and Q₈ acts semiregularly on the rows of Paley's matrix with its centre acting as −I. So, in a basis
adapted to the N = (q + 1)/4 orbits, Paley's matrix is a Williamson array H = Σ R_u ⊗ X_u of four blocks of order N
(`search/quaternion_type.py pattern q`; checked for q = 11, 19, 23, 43 and 1051). That is not enough:
- **The blocks are not amicable.** H Hᵀ = 4N I only forces [X_i, X_1] = −[X_k, X_j], [X_j, X_1] = −[X_i, X_k] and
  [X_k, X_1] = −[X_j, X_i], with [X, Y] = X Yᵀ − Y Xᵀ, and the commutators are non-zero. These are Hadamard matrices
  of quaternion type (with a semiregular Q₈ action), but not Williamson-type matrices.
- **No base point helps.** Changing the base point of each orbit, independently for rows and columns, mixes the blocks
  entrywise. For q = 19 (N = 5) none of the 4¹⁰ choices gives pairwise amicable blocks (`rebase 19`); for q = 11 160 of
  4096 do, but Williamson-type matrices of order 3 are classical anyway.
- **The weaker condition does not survive plugging.** In an orthogonal design Σ x_u M_u the cross terms are
  Σ M_u M_vᵀ ⊗ [X_u, X_v], so quaternion-type blocks need M_1 M_iᵀ = ±M_j M_kᵀ and its cyclic versions. The
  Cooper–Wallis Baumert–Hall arrays of T-matrices of orders 3, 5 and 7 do not satisfy this: for q = 19 none of the 384
  relabellings and signs of the variables gives a Hadamard matrix (`plug 19 t`).

So the dicyclic structure gives quaternion-type matrices of every order 4N with 4N − 1 a prime power, but not the
Williamson-type matrices of order N that Cooper–Wallis needs, and 263 stays out of reach this way. A product that needs
no amicability exists: since left and right multiplications commute, Σ_{u,v} L_u R_v ⊗ T_v ⊗ X_u is a Hadamard matrix
of order 4tN whenever Σ L_u ⊗ X_u is one and Σ R_v ⊗ T_v is a t × t orthogonal matrix over {±1, ±i, ±j, ±k} (for
circulants, a perfect sequence over Q₈). But such matrices exist only for even t: an inner product of two rows is a sum
of t units of Q₈, each contributing ±1 to a single coordinate, so it is non-zero when t is odd (brute force agrees: no
perfect Q₈ sequence of length 3, 5, 7 or 9, many of lengths 2, 4 and 6). The product therefore only reaches multiples
of 8.

**Plug-ins from association schemes (negative).** The symmetric plug-in of [search/NOTES.md](search/NOTES.md) §31 needs
commuting symmetric ±1 matrices X, D of order r with X² + (m − 1) D² = m r I. With X and D in the Bose–Mesner algebra of
a strongly regular graph (every feasible parameter set with r odd, r ≤ 800, and every X, D = J − 2(sum of classes)), the
only solutions are the known families: D the core of a regular symmetric Hadamard matrix of order r + 1 = 4u², and X = J
or J − 2I. In the cyclotomic schemes of GF(q) (q prime < 3000, e = 3, 4, 5, 6, 8 classes, −1 ∈ C₀) no pair satisfies the
equation with m ≡ 0 (mod 4).

**The structure of the matrices on Zaitsev's page.** Their monomial automorphism groups (computed with nauty on the
usual 8n-vertex graph) have orders 32 at 238, as for ours (the multiplier group ⟨8⟩ times the four scalars), 8 at 470
and 24 at 490, with two fixed rows in both. So the orders the page attributes to unpublished approaches come with
very little symmetry: an involution at 470, a group of order 6 at 490.

## 22. Paley's matrix in sixteen circulant blocks: 14 FRONTIER orders (October 16)

§21 left 263 out of reach because the quaternion-type blocks of Paley's matrix are not amicable. A different block
structure and a product that never asks for amicability settle the four orders of 263 and ten more
([hadamard/PALEY_BLOCKS.md](hadamard/PALEY_BLOCKS.md)).

- **Circulant blocks.** For q ≡ 3 (mod 8), s = (q + 1)/4 is odd. An element g of order s in a non-split torus of
  PSL(2, q) has four orbits of length s on PG(1, q), and Paley's matrix satisfies H(gx, gy) = e(x) e(y) H(x, y). Since
  g^s = 1, the products of e over the four orbits are equal; with c this common product (c^s = c, s odd), signs with
  E(gx) = c e(x) E(x) close up along every orbit, and g is an automorphism of E H E, whose blocks along the orbits are
  circulant (`search/paley_blocks.py check q`). Matrices of the same kind are in the literature: Ito's relative
  difference sets in the dicyclic group of order 8s, for every prime power 4s − 1, give Hadamard matrices in Schmidt's
  form (9), whose sixteen blocks are ±A, ±B, ±C, ±D or their transposes for circulants with A Aᵀ + B Bᵀ + C Cᵀ + D Dᵀ =
  4s I and A Bᵀ − B Aᵀ + C Dᵀ − D Cᵀ = 0 (Schmidt 1999, Theorems 2.2 and 3.3), and Yamada's generalized quaternion
  type gives them for Paley's matrix (Seberry–Yamada 1993, Theorem 1).
- **The product.** With T-matrices T₁, …, T₄ of order t, A_j = Σ_k T_k ⊗ N_kj has entries ±1 (disjoint supports)
  and Σ_j A_j A_jᵀ = Σ_k T_k T_kᵀ ⊗ 4s I = 4st I, because the block rows of E H E are orthogonal. Only that
  orthogonality is used, no commutation of the blocks (amicability). This is Seberry–Yamada (1990),
  Theorem 4, with an OD(4s; 4s) of sixteen circulant blocks. For gcd(s, t) = 1 the A_j are circulants of length st.
- **Reach.** Against FRONTIER.md, with every prime power q ≡ 3 (mod 8), every T-matrix order we know (odd t ≤ 210
  outside Seberry–Yamada's 26 exceptions, plus 97, 103, 109, 113, 131, 151, 163, 181 and 197, the Golay numbers plus
  one, 3m − 1 for TT(m) with m ≤ 44 even, 2m + 1 for base sequences with m ≤ 40) and gcd(s, t) > 1 allowed
  (group-developed Goethals–Seidel), the route reaches exactly 14 entries: s = 263 (t = 3, 5, 7, 11), 503 (3, 5),
  227 (7, 13), 563 (3, 5), 383, 431 (5) and 887, 911 (3). It also reaches six entries of Cati–Pasechnik Table 4 that
  other 2026 constructions have settled (1359, 1441, 1703, 2227, 2335, 2489). All 14 are certified and pass
  `verify/verify.py`, which rebuilds Paley's matrix from the symplectic form on PG(1, q).
- **Why it was missed.** Each step is published; we found no source that combines them. Seberry–Yamada (1993) use
  the generalized quaternion structure only in products that keep it (their Theorem 6 asks T₀T₁ᵀ − T₁T₀ᵀ + … = 0,
  which ordinary T-matrices do not satisfy) and reach 2ʳ m₁⋯mᵣ (q + 1) with Williamson orders mᵢ. Schmidt uses the
  four-way split only for a three-squares identity. Balonin and Đoković (arXiv:1508.00640, §9) compute the blocks
  for q = 139 and observe that they are not of Williamson type, so Cooper–Wallis does not apply; Theorem 4 does
  not need it. The novelty audit of the route (three literature sweeps and a critic) is summarised in
  PALEY_BLOCKS.md.
- **What does not extend.** The other FRONTIER entries n = s·t with t a T-matrix order would need a Hadamard matrix of
  order 4s with sixteen circulant blocks (for instance from a Hadamard group Q₈ₛ, as Ito conjectured for every s), for s
  = 167, 179, 191, 223, 239, 269, 283, 347, 359, 419, 443, 463, 479, 491, 571, 573, 599, 631, 643, 863, 941, 947, 955,
  971 or 991. We know none: 4s − 1 and 2s − 1 are not prime powers, and Williamson-type matrices of these orders are not
  known to us. Fewer rows would do with fewer T-matrices: T-matrices with T₄ = 0 (t = 3, 5 and every Golay number plus
  one) need only three block rows with circulant blocks, orthogonal to each other. But a Goethals–Seidel matrix has at
  most two such rows (its blocks are circulant on the diagonal and back-circulant off it, and no affine relabelling x ↦
  ax + b of each block row and column makes three rows circulant, since that would need a = −a), and two rows are
  useless because T-matrices of odd order t > 1 with two zero matrices do not exist (T₁ ± T₂ would be a periodic
  complementary ±1 pair of odd length, but PAF ≡ t (mod 4) at every shift).
- **The same object as a quaternary pair.** A matrix in Schmidt's form (9) is a quaternary periodic complementary pair
  of length s: P = (1 − i)(A + iB)/2 and Q = (1 − i)(C + iD)/2 have entries ±1, ±i, and P P* + Q Q* = 2s I is exactly
  the two conditions above. So each s of the list would also give the complex order 2s (334 for s = 167, the smallest
  open complex order). With multiplier symmetry the counting is poor: the row sums leave only |H| = 5 at s = 191 and
  only the trivial group at s = 167 and 179, and the crude count (4s bits against s − 1 real conditions of about four
  bits each) gives about one solution in a space of 2¹⁵⁶ at s = 191, |H| = 5. We did not search.
- **No skew version.** A₁ = Σ T_k ⊗ N_k1 is skew-type only if T₁ = δ₀ and every other T_k with T_k ≠ 0 is odd
  (T_k(−x) = −T_k(x)) with N_k1 a symmetric circulant. In Paley's array each block column has exactly one symmetric
  off-diagonal block (checked for q = 43, …, 227), so at most T₁ and one odd T₂ survive, and a ternary odd T₂ with
  T₂T₂* = (t − 1)δ₀ cannot exist for odd t > 1 (PAF at a nonzero shift is a sum of t − 2 terms ±1).
- **No product with arbitrary Goethals–Seidel quadruples.** X_j = Σ_k A_k ⊗ M′_kj would turn every Goethals–Seidel
  quadruple of length s into one of length st if M′ were a 4 × 4 array of (0, ±1) circulants of order t whose block
  columns are T-matrix quadruples, pairwise orthogonal. Exhaustively there is none for t = 3 (192 admissible columns,
  no 4-clique of orthogonal ones) or t = 5 (3,840 columns). In the arrangement of form (9) the condition is that
  U = T₁ + iT₂ and V = T₃ + iT₄ be complementary, and that is impossible for odd t > 1: W = U + V has unit entries,
  U − V = εW with ε = ±1, and U U* + V V* = t I would need Σ W(x) W̄(x + d) = 0 over the x with ε(x) = ε(x + d), a
  sum of an odd number of units of ℤ[i] (the number of agreements of a ±1 sequence of odd length is odd).

The tables of §9 ("no product is overlooked") missed this product; their entries WT(263) and WT(227) are no longer
needed. The novelty re-audit of October 16 found that the same route also reaches three of our orders beyond the
tables, so they were obtainable classically: 19028 = 4 · 71 · 67 (q = 283, Sawade's T-matrices of order 67), 39260 =
4 · 755 · 13 (q = 3019) and, through Miyamoto's theorem, 21748 = 4 · 5437 from 5436 = 4 · 453 · 3 (q = 1811). It also
found that SageMath's development version builds 31316 since October 3, 2026 (Miyamoto from a Cooper–Wallis matrix
of order 7828 with the T-matrices of order 103 decoded from Schneider's matrix of order 2060), nine days before our
commit, so 31316 is Known. The order READMEs say so.

Three more orders follow by Miyamoto's theorem. q = 4n + 1 is a prime for four of the fourteen; for 6037 the order
4q = 24148 is Paley II already (2q − 1 = 12073 is a prime), and the other three give 21044 = 4 · 5261
(from 5260; Handbook Table 1.53 gives 5261 with t = 6), 40244 = 4 · 10061 (from 10060) and 45044 = 4 · 11261 (from
11260), the last two beyond every table. Paley I and II, Yamada–Kiyasu and the product theorems do not reach them,
and they are certified and verified like the others. For the other ten, 4n + 1 is composite.

**Beyond Table 4 (October 17).** The route also reaches 138 of the 691 entries 3001 ≤ n ≤ 9999 of Handbook Table 1.53
that had no valid construction: 637 that it still lists with t ≥ 3 after Đoković's updates, and 54 recorded as known
only through Corollary 8.8(i) (§7: 13 settled by Đoković in 2010, 33 of his set Δ and 8 with t = 2). Of the 120 it
reaches with t ≥ 3 or from 2010 we built 117, of orders 12052 to 39916
([hadamard/PALEY_BLOCKS.md](hadamard/PALEY_BLOCKS.md#beyond-table-4)): 109 with T-matrices known before 2026 (one of
them 32316, from the 13) (orders 3 to 73, from base sequences, Golay pairs, SageMath's stored T-sequences, and
multiplier-subgroup searches for 19, 31, 41, 49 and 73), six with the 2026 T-matrices of orders 103, 151 and 163, and
16892 and 25132, which SageMath's development version builds. Two extensions were needed: the group-developed array over
ℤ_t × ℤ_s when gcd(s, t) > 1 (seven orders), and Paley's matrix over GF(3⁵) (s = 61, three orders), with the field
recorded in the certificate. Of the other four, three were already in the repository. For 28820 = 4 · 131 · 55 we built
T-matrices of order 55 from base sequences BS(6, 5) by Yang's multiplication by 5 (Seberry–Yamada 1992, Theorem 5.23,
with the Golay pair of length 2; `t-matrices/T55.txt`). The other 18 reached are among the orders recorded only through
Corollary 8.8(i), all built (§7). Nothing of this reaches FRONTIER.md. A Paley II analogue (s = (q + 1)/2) would add
nothing, because those s are Turyn's Williamson orders and Cooper–Wallis already covers them. No ordering, sign or
reversal of a Paley-block quadruple satisfies the form-(9) condition ((q, t) = (11, 5), (19, 3), (43, 3), (43, 5) and
(67, 3)), as the parity argument above predicts, so the route gives no complex orders.

## 23. Miyamoto's matrix in bordered blocks, times T-matrices: 15 orders beyond Table 4 (October 17)

The Paley-block route of §22 needs a Hadamard matrix of order 4s made of sixteen blocks in a commutative algebra with a
reversal. Miyamoto's matrix of order 4q (q ≡ 1 mod 4 prime) is one, provided its Hadamard matrix K of order q − 1 is of
2N-type, K = [[A, B], [−Bᵀ, Aᵀ]] with A, B negacyclic, i.e. a negaperiodic Golay pair of length (q − 1)/2
(Balonin–Đoković, arXiv:1508.00640). With Paley's matrix indexed by ∞, 0, g^(2j), g^(2j+1), each block is σ[[1, ±e],
[±eᵀ, J₂ ⊗ U + D₂ ⊗ V]] with U circulant and V negacyclic over ℤ_(q−1)/2. Bordered matrices [[x, u e], [u eᵀ, M]] with M
in the algebra spanned by J₂ ⊗ circulants and D₂ ⊗ negacyclics and x − r(M) = u form a commutative algebra (J₂D₂ = 0;
the borders of XY and YX agree because w(x − μ) = u(y − ν)), and the permutation that fixes the border and sends (h, j)
to (h, −j) for j = 0 and (1 − h, −j) otherwise transposes every block (on the D₂ part it is the reversal times diag(1,
−1, …, −1), which transposes negacyclic matrices). Seberry–Yamada's multiplication by T-matrices (Theorem 4 of 1990)
then goes through over this algebra: H(4qt). The observation came from the frontier-routes workflow; it was checked
there for q = 13, 29, 37 with Williamson-based K, and here independently (`search/miyamoto_blocks.py selftest`: q = 13,
29, 37, 41, 53, 61, with NG pairs from both Paley series, Williamson matrices of orders 9 and 13 and Turyn's
multiplication by the Golay pair of length 2, times T-matrices of orders 3, 5 and 7). The bordered form itself is
Seberry and Yamada's: they write Miyamoto's matrix as an M-structure (JCMCC 7 (1990), Theorem 34 and Corollary 35;
survey 1992, Theorem 8.10 and Corollary 8.11), but their M-structure multiplications need circulant or type-one blocks
(Theorem 4) or sixteen matrices T_ij with row and column orthogonality (Theorem 21), which exist for no odd t > 1: for
each entry the 4 × 4 slice (T_ik(r, c)) is a signed permutation matrix, and a sum of t of them over c cannot vanish mod
2. So they never multiply Miyamoto's matrix by T-matrices; what is new here is the commuting choice of K and of Paley's
indexing. The novelty audit of October 17 found this attribution; chapter 6 of Seberry–Yamada's 2020 book (§6.4, §6.5)
was not accessible.

- *NG pairs.* Proposition 6 of Balonin–Đoković (length p + 1, p an odd prime power; for p ≡ 1 (mod 4) their first Paley
  series) and the second Paley series (length (p + 1)/2, p ≡ 3 mod 4) come from a negacyclic conference matrix on PG(1,
  p): c_k = χ(Δ(1, w^k)/Δ(1, w)) with Δ(x, y) = x y^p − x^p y and w primitive in GF(p²); Δ(w^i, w^j) = N(w)^i Δ(1,
  w^(j−i)) with χ(N(w)) = −1 makes it negacyclic after the signs (−1)^i. Williamson (or quasi-Williamson) matrices of
  odd order w give length 2w through the array (29) of Balonin–Đoković, read in ℤ[y]/(y^(2w) + 1) ≅ ℤ[x]/(x² + 1) ⊗
  ℤ[ℤ_w] (y ↦ xz, y^j = (−1)^⌊j/2⌋ x^(j mod 2) z^(j mod w)); the CRT reindexing without these signs does not give a
  negacyclic pair. Turyn's multiplication by a Golay pair multiplies the length (their (15)–(16)).
- *Reach.* Among the Handbook entries without a valid construction (§22), with q prime and the NG pairs above, the route
  reaches 27: 11 also reached by Paley blocks (already here), 6031 (here by Yamada–Kiyasu, order 24124) and 15 new ones,
  3279 to 9991 ([hadamard/MIYAMOTO_BLOCKS.md](hadamard/MIYAMOTO_BLOCKS.md)). Eleven use T-matrices of order 3, 5, 7 or
  9: nine are obtainable classically and two are First explicit constructions (below); four use the 2026 T-matrices of
  orders 103, 151 and 163 and also follow by Cooper–Wallis. Two, 7179 = 3 · 2393 (Golay 26 × Williamson matrices of
  order 23) and 7323 = 3 · 2441 (Golay 10 × the first Paley series over GF(121)), are among the 13 that Đoković claimed
  only through Corollary 8.8(i) (§7). Prime powers q (25, 49) add 4075 and 6419, Paley-block orders already here;
  Williamson matrices of order 153 (found with a multiplier group of order 16) give 3065 and 6743 (orders 12260, 26972),
  entries recorded only through Corollary 8.8(i); those of order 325 would give 3903 and 9107 (not found by the
  program).
- *FRONTIER, negaperiodic pairs.* None: every FRONTIER entry n = qt needs an NG pair that is open (2427 is reached with
  a periodic pair, below). 1883 = 7 · 269 needs quasi-Williamson matrices of order 67, 2823 = 3 · 941 of order 235, and
  the orders reopened in §7 need 93 (2611), 127 (1527, 2545), or NG pairs of length 404 (2427) and 428 (2571). Balonin
  and Đoković list quasi-Williamson matrices of orders 47, 59, 65, 67, 73, 81, 89, 93, 101, 103, 107, 109, 113 and 119
  as open since 1999.
- *Periodic Golay pairs: FRONTIER 2427.* The algebra needs only that the sixteen V_rc commute and are transposed by the
  reversal that transposes the U_rc. For a periodic Golay pair (a, b) of length m, K = [[A, B], [−Bᵀ, Aᵀ]] with A, B
  circulant is a Hadamard matrix, every V_rc is circulant, and the plain reversal (h, j) ↦ (h, −j) transposes every
  block (J₂ and D₂ are symmetric); commutation is as before. Formulas (15)–(16) also multiply a Golay pair by a periodic
  Golay pair (Georgiou–Stylianou–Drosou–Koukouvinos 2014). For q = 809, m = 404 = 2 · 202 and Đoković–Kotsireas (Numer.
  Algorithms 69 (2015), arXiv:1310.5773) give a periodic Golay pair of length 202 (orbits of {1, 87, 95, 137, 185} on
  ℤ₂₀₂; transcribed and checked against the arXiv source). With T-matrices of order 3 this gives H(9708), n = 2427 = 3 ·
  809, verified densely: the first valid construction of an order that Table 4 omits only through Williamson-type
  matrices of order 809 from Corollary 8.8(i), and that FRONTIER.md listed as unresolved. With T-matrices of order 11 it
  gives H(35596), n = 8899 (Đoković's Δ, same dependence). The variant came out of the novelty audit of this route (an
  audit agent tested it on q = 809 and built H(9708) independently, with its own implementation of Miyamoto's array).
  Over all FRONTIER entries and the 691 Handbook entries, periodic pairs of lengths g · v (g a Golay number, v ∈ {1, 34,
  50, 58, 72, 74, 82, 122, 202, 226}) add only these two. 2611, 1527, 2545, 2571 and 1883 would need pairs of length
  186, 254, 428 or 134; periodic Golay pairs of these lengths cannot exist (2m = 4 · 93, 4 · 127, 8 · 107 and 4 · 67 are
  not sums of two squares), and negaperiodic ones are open.
- *T-matrices of order 55.* 9515 and 28820 needed them and SageMath stores none. Yang's multiplication by 5
  (Seberry–Yamada 1992, Theorem 5.23, Golay pair of length 2) applied to base sequences BS(6, 5) gives T-sequences of
  length 55: in the template X = (A, C; 0; ±A, ±C; 0; ±B*, 0), Y likewise with B, D, Z, W on the even blocks and the
  last block (0, ±D*), (0, ±C*), 64 of the 4096 sign choices work (t-matrices/T55.txt). Two attempts to read the
  multiplication off Yang's Lagrange identity with block-concatenated or interleaved supports gave no layout at all:
  with T-sequence outer factors, three of the four outer supports must be symmetric about a common centre, which no
  TS(5) satisfies.
- *Searches for the missing pairs (October 17).* (i) Near-Williamson quadruples (A, B, B, D) with A, B symmetric and D
  free (`propus.py --nwrep`, the repeated-block input that repairs Corollary 8.8) over ℤ_v for v = 93, 67, 127, 43 and
  73, every multiplier structure whose pools fit (52 exact joins): none. (ii) Quasi-Williamson matrices of orders 67,
  93, 107, 127, 235 and the other orders of Balonin–Đoković's list (31), constant on the orbits of a multiplier group of
  odd order (so not symmetric), joined on PAF sums and the antisymmetric part of the cross-correlations
  (`search/qw_mult.py`), for every such group with at most 15 orbits: none. (iii) NG pairs directly, as two sequences in
  ℤ[y]/(y^v + 1) constant up to sign on the orbits of a multiplier group y ↦ y^k, filtered by |p(ζ)|² ≤ 2v at the odd
  2v-th roots of unity and joined exactly on the NAF (`search/ng_mult.py`, which finds the Paley and other pairs of
  lengths 14 to 38 at once): for v = 134 (QW(67), FRONTIER 1883) the structures with up to 24 orbits give none (24
  orbits: 94578 of 8388608 sequences pass the spectral test, no pair), and the next structure has 46 orbits; v = 94
  (QW(47)) has no structure between 6 and 48 orbits. For v = 186 (2611), 214 (2571, through the doubling to 428) and 254
  (1527, 2545) every structure with at most 24, 6 and 20 orbits gives none; the next ones have 28, 108 and 30 orbits. A
  Gaussian count predicts about 70 pairs for the 24-orbit structure of v = 134, so multiplier-invariant NG pairs are
  much rarer than the count suggests. (iv) Williamson matrices of order 153 constant on the orbits of a multiplier group
  of order 16 exist and were found at once (12260, 26972); `williamson.py` with up to 24 orbits finds none of order 325
  (for 3903, 9107) so far.

## 24. Which Hadamard matrices T-matrices multiply; NG pairs and Yamada's bordered blocks (October 17)

Theorem 4 of Seberry–Yamada (1990) uses only the orthogonality of the block rows of H = (N_kj), plus the fact that A_j =
Σ_k T_k ⊗ N_kj can go into the Goethals–Seidel array. So it applies to any Hadamard matrix of order 4s on whose rows and
columns an abelian group of order s acts semiregularly, with four orbits each, by signed permutations. In an orbit basis
the sixteen blocks are group-developed, and the inversion of the group transposes them all. Every cocyclic Hadamard
matrix over a group with an abelian subgroup of odd order s qualifies: the expanded group splits over that subgroup
because s is odd. §22 is the case of a cyclic group from PSL(2, q); §23 is a bordered analogue, and so is Yamada's
matrix below. Three more cases were examined; none reaches FRONTIER.md.

- *Negaperiodic and periodic Golay pairs of length 2s, s odd.* With y ↦ −ix, ℤ[y]/(y^(2s) + 1) ≅ ℤ[i] ⊗ ℤ[ℤ_s]. So an NG
  pair of length 2s, equivalently quasi-Williamson matrices of order s (Balonin–Đoković, arXiv:1508.00640, §9, after Ito
  and Schmidt), is a Hadamard matrix of order 4s in sixteen circulant blocks of order s; so is a periodic Golay pair of
  length 2s. Times T-matrices of order t this gives H(4st). The second Paley series (p = 4s − 1) is the Paley-block
  route of §22. The first series (2s − 1 a prime power) is Paley II, whose s are Turyn's Williamson orders. Williamson
  matrices are the classical case. Every FRONTIER entry n = s · t with T-matrices of order t known needs s ∈ {167, 179,
  191, 223, 239, 269, 283, 347, 359, 419, 443, 463, 479, 491, 509, 571, 573, 599, 631, 643, 857, 863, 941, 947, 955,
  971, 991}. No quasi-Williamson matrices of any of these orders are known, and an NG pair of length 334 would also
  settle the open complex order 334.
  Balonin and Đoković (arXiv:1508.00640) show that NG pairs of every even length are equivalent to Ito's conjecture,
  so for s = 167, 179, 191 and 223 the route is Ito's conjecture at t = s. An exact meet-in-the-middle over NG pairs of
  lengths 334, 358, 382 and 446 invariant under every admissible multiplier group with at most 24 sign pairs (|H| =
  166, 83; 178, 89; 190, 95, 38, 19; 222, 111, 74, 37) found none; the same program finds NG pairs of lengths 26–42
  under groups of order 2 to 6 (`search/ng_mitm.py`).
- *Yamada's matrix of order 4(q + 2) (Kiyasu; search/NOTES.md §13) in bordered blocks.* Its inner blocks lie in P ⊗
  circ(n₀) ⊕ J₂ ⊗ alg(S, J), with P = [[1, −1], [−1, 1]], PJ₂ = 0, S the skew core of order n₀ = (q + 1)/2 and SJ = 0,
  so they commute. Reassign the four border columns to the block columns in the order (0, 3, 1, 2) and sign them (1, −1,
  −1, −1). Then every block restricted to span(border, e) is symmetric and lies in one two-dimensional algebra; 16 of
  the 2 · 4!² · 2⁸ choices of border order, signs and reversal parity work. So the blocks commute, and diag(1, I₂ ⊗ R)
  transposes all of them whenever R is an involutive permutation with R S R = Sᵀ, as x ↦ −x is for Paley cores. Times
  T-matrices this gives H(4(q + 2)t), checked exactly for (q, t) = (5, 3), (5, 5), (13, 3), (13, 7) and (37, 3) (Paley
  cores). Nothing is new: the same inputs, a skew core with such an R and symmetric circulants M, N of order n₀ with MMᵀ
  + NNᵀ = qI, are what Seberry–Yamada's Lemma 26 (JCMCC 7 (1990), Corollary 27(i)) turns into symmetric Williamson-type
  matrices of order q + 2, and Cooper–Wallis does the rest.
- *New cores for Lemma 26.* FRONTIER 2315 = 5 · 463 would follow from a skew-Hadamard matrix of order 232 with an
  involutive signed anti-automorphism fixing exactly two points (necessarily with signs + and −). The structured
  versions fail. (a) The back-circulant input U₃ = CR, U₄ = DR, with C, D circulants of order 231 that are amicable
  (Ĉ·conj(D̂) real) and satisfy PAF_C + PAF_D = −2: at a character of order 3, Ĉ = s v and D̂ = r v with v ∈ ℤ[ω], so
  (s² + r²) N(v) = 462. Here N(v) ∈ {1, 3, 7, 21}, and none of 462, 154, 66, 22 is a sum of two squares. (b) Two borders
  and two ℤ₁₁₅-orbits that the anti-automorphism swaps: the orthogonality of the two block rows forces KM = 0, so |M̂|²
  or |K̂|² is 231 at a character of order 5, impossible since 3 is inert in ℚ(ζ₅). (c) One border and three ℤ₇₇-orbits,
  one fixed and two swapped: |Â + 2D̂|² = 231 (or |Â|² = 231) at every non-trivial character, impossible at order 7
  since 3 is inert in ℚ(ζ₇). Bordered two-circulant skew matrices from a skew-type U and a symmetric V with PAF_U +
  PAF_V = −2 (J. Wallis, Yamada: w = s² + 4 a prime power) do have such an anti-automorphism: x ↦ −x on both orbits with
  signs + and −, borders swapped. Lemma 26 turns their cores into symmetric Williamson-type matrices of order 4w + 3
  when 4w + 1 is a prime power (checked for w = 13: order 55); the next case is w = 1093 (order 4375). Szekeres's pairs
  (Remark 18) give cores of prime-power order, which Corollary 27 already covers.
- *The other M-structures of JCMCC 7* (Corollaries 36–45, eight and sixty-four blocks) give orders 8m only, outside
  FRONTIER.md.

## 25. Miyamoto blocks with near-Williamson matrices: 1527, 1883, 2545 and 2611 (October 17)

§23 needed the four quadrants of the Hadamard matrix K of order q − 1 in Miyamoto's array to lie in a commutative
algebra transposed by one signed reversal, and looked for K among negaperiodic or periodic Golay pairs (two-block
matrices of 2N or 2C type). A third kind works. If X₁, …, X₄ are pairwise commuting symmetric ±1 matrices of order w
with Σ Xᵢ² = 4wI, their Williamson array has quadrants [[X, Y], [−Y, X]], which form a commutative algebra closed under
transposition, and the half swap transposes them. Near-Williamson matrices give such Xᵢ: A, B, C, D circulant with B, C,
D symmetric (arXiv:2605.08661v2, September 2026; good matrices and Williamson matrices are special cases), and Xᵢ = (AR,
B, C, D) with R the back-identity, since AR is symmetric and commutes with every symmetric circulant. Relabelled so that
the half swap is j ↦ 1 − j on ℤ_m, which also reverses Paley's circulant blocks, the sixteen bordered blocks commute and
the permutation fixing the border and sending (b, j) to (b, 1 − j) transposes them
([hadamard/MIYAMOTO_BLOCKS.md](hadamard/MIYAMOTO_BLOCKS.md#near-williamson-matrices)).

The judge stage of the frontier workflow of October 17 (scratch work, not in the repository) built H(7532) this way with
near-Williamson matrices of order 67 and its own code; the run ended before its result was read. The saved matrix checks
exactly (±1 entries, H Hᵀ = 7532 I). `search/miyamoto_blocks.py` now takes the near-Williamson matrices as SPEC
`nw:FILE` and stores them as "nw" in the certificate, and `verify/verify.py` rebuilds K from them (`miyamoto_nw_k`, with
the back-identity, so the matrix differs from the judge's). All four orders are verified densely:

- 7532 = 4 · 1883 (q = 269, near-Williamson matrices of order 67, T-matrices of order 7): New. Table 4 lists 1883 as
  open, and no one claimed it through Corollary 8.8(i).
- 6108 = 4 · 1527 and 10180 = 4 · 2545 (q = 509, order 127 from the database of the 2026 paper, T-matrices of order 3
  and 5): First explicit constructions. Table 4 lists both as open; Đoković claimed them in 2010 through Corollary
  8.8(i) with his good matrices of order 127, which satisfy the same hypotheses.
- 10444 = 4 · 2611 (q = 373, order 93 from the database, T-matrices of order 7): First explicit construction. Table 4
  omits it through Seberry–Yamada's Table A.2 (Williamson-type matrices of order 373 from Corollary 8.8(i)).

This does not give Williamson-type matrices of order q; Corollary 8.8(i) stays unproven in general (§7; for some q a
repeated block repairs it, §26). FRONTIER.md now has 130 entries. Of its entries of the form n = qt with q ≡ 1 (mod 4) a
prime and T-matrices of order t known, the route still misses 2571 (q = 857, w = 214), 2823 (q = 941, w = 235) and 2677
(t = 1, which needs H(2676)). Of the database's existence orders (53 to 181), those with q = 4w + 1 prime are 67, 73,
79, 93, 97 and 127. Among the Handbook entries with t ≥ 3 they reach n = 3497, 4573, 8339, 9953 (q = 269), 9827 (317),
6341, 8579 (373), 6617, 8653 and 9671 (509), now all built (hadamard/MIYAMOTO_BLOCKS.md; T-matrices of order 37 from
`search/tmat_mult.py 37 10`); 4981, 5567, 7325, 9083 (q = 293) and 9325 (373) are reached too but were built by Paley
blocks. 6341 and 8579 are in Đoković's set Δ and 6617, 8653 and 9671 among his 2010 claims, all through Corollary
8.8(i). 9827 is in Δ as well, validly: Cooper–Wallis with Williamson-type matrices of order 317 from Corollary 8.8(ii)
and Paley's conference matrix of order 158. So 39308 is a new matrix of a known order; an earlier version of this
section called it new beyond the tables, and the audit of October 17 corrected it. The other four are new beyond the
tables. For q = 509 Đoković's good matrices of order 127 (1993) satisfy the hypotheses too, so 1527, 2545 and the 509 ·
t entries were obtainable classically, and for q = 317 Turyn's Williamson matrices of order 79 would do; for q = 269 and
373 the inputs are the 2026 near-Williamson matrices (J. Wallis's Williamson-type matrices of order 93 are not
symmetric, and the propus symmetric Hadamard matrices of orders 268 and 372 do not have commuting quadrants, checked by
the audit).

## 26. Miyamoto blocks with symmetric Williamson-type matrices, and the corollary with a repeated block (October 17)

§25 used only that X₁, …, X₄ are symmetric, commute pairwise and have Σ Xᵢ² = 4wI. Seberry and Yamada's corollaries of
1990 make many such matrices: Corollary 25 (order 2uv from orders u and v), 27 (order p + 2), 29(ii) (order q from a
symmetric conference or Hadamard matrix of order (q − 1)/2) and 31 (order 2p + 1), seeded by circulant Williamson
matrices and Turyn's family of order (p + 1)/2. `search/symwt.py` builds them all, over GF(p^k) where needed (Yamada's
circulants α, β with α² + β² = pI from GF(p²), Miyamoto's lemma over GF(729), Paley II over GF(49)), and checks them.
`miyamoto_blocks.py certify q t "wt:RECIPE"` stores the four matrices (upper triangles, base64), and `verify.py` checks
them before building K (`miyamoto_wt_k`; a flipped bit is rejected).

Over the Handbook entries 3001 ≤ n ≤ 9999 with t ≥ 3, or recorded only through Corollary 8.8(i), the route reaches 148.
124 already had valid constructions (Đoković's set Δ and his 2010 propositions, or orders built here), and 8
corollary-only ones were built here before by other routes. The other sixteen are built (dense checks of verify.py:
hadamard/MIYAMOTO_BLOCKS.md), all First explicit constructions: 3579 and 5965 (q = 1193, w = 298 = 2 · 1 · 149), 3865,
5411 and 8503 (773, w = 193), 4765 (953, w = 238 = 2 · 7 · 17), 5199 and 8665 (1733, w = 433), 5545 (1109, w = 277),
5979 and 9965 (1993, w = 498 = 2 · 3 · 83), 8187 (2729, w = 682 = 2 · 11 · 31), 8465 (1693, w = 423), 8751 (2917, w =
729), 8859 (2953, w = 738 = 2 · 9 · 41) and 9627 (3209, w = 802 = 2 · 1 · 401), with T-matrices of order 3, 5, 7 or 11
(hadamard/MIYAMOTO_BLOCKS.md). Among the 85 other reach orders not built here none lies outside Δ and Đoković's 2010
propositions. FRONTIER.md gains nothing: 2571 and 2823 would need order 214 = 2 · 107 or 235, which none of the
corollaries gives.

The search for inputs showed that the corollary itself can be repaired for seven primes: Seberry–Yamada's symmetric
matrices from Corollaries 8.8(ii) with M₃ = M₄, 27 and 31 have X₃ = X₄, and the repeated-block form of Lemma 26 (§7)
gives symmetric Williamson-type matrices of orders 293, 773, 1453, 1637, 1693, 1733 and 2917, checked exactly. §7 now
says which entries this repairs (also the complex order 1546). The audit added Kronecker products of symmetric
quaternary Hadamard matrices (§7), which repair 521, 953, 1993, 2729, 2953 and 3209 as well; of the sixteen orders only
3579, 5545 and 5965 (q = 1193 and 1109) keep no Williamson-type route that we know. All sixteen keep the label First
explicit construction (NOVELTY.md): the tables' printed route fails, and the repairs were found here.

The audit of §25 (workflow of October 17: an independent rebuild of the four near-Williamson orders, which reproduced
their SHA-256 values, two novelty lenses and a claims check) found that 9827 is in Đoković's set Δ through valid
Williamson-type matrices of order 317, so 39308 is relabelled Known, and that Verma's symmetric conference matrix of
order 86 (arXiv:2610.08893, 6 October 2026) settles Williamson-type matrices of order 173 through Corollary 8.8(ii); we
checked the matrix exactly. Smaller corrections (citations of Đoković's propositions, the range of Balonin–Đoković's
list (31), the provenance of the near-Williamson files, the Table A.17 inference) are applied in the READMEs.

Negative: no circulant Williamson matrices (A, B, C, C) of orders 43, 93, 153, 189, 277, 325 or 373 in the multiplier
structures with at most 18 orbits (some row-sum classes skipped), so the repair does not reach 173, 373, 613, 757, 1109,
1301 or 1493 this way.

## 27. Skew-Hadamard matrices of order 3(q + 1) from Paley's blocks (October 17)

The T-matrix multiplication of §22 gives Hadamard matrices; we asked when its first block A_c = Σ_k T_k ⊗ N_kc can be
skew-type. Scanning the sixteen blocks of Paley's matrix for q = 11 to 331 showed the same pattern in every case: in
some block column c the diagonal block is skew-type (as always), one block is symmetric up to a shift, and the other
two are paired, N_dc(y) = ε N_bc(δ − y). After the base point of orbit d is moved by −δ the pair is N_dc = ε N_bcᵀ, and
the T-matrices of order 3, T_c = I, T_b = Z, T_d = −εZ², 0 on the fourth orbit, make A_c skew-type. The Goethals–Seidel
array with A_c first is then a skew-Hadamard matrix of order 12s = 3(q + 1).

The pattern has a short explanation, which turns the observation into a theorem for every prime power q ≡ 3 (mod 4)
([skew/PALEY_SKEW.md](skew/PALEY_SKEW.md)). On the points z^m of PG(1, q) (z primitive in GF(q²) = GF(q)(i)), Paley's
matrix is [m = n] + (−1)^m χ(Im z^(n−m)); its sixteen blocks along the orbits of z⁴ are negacyclic; and Frobenius, an
anti-automorphism (det(ū, v̄) = −det(u, v)), fixes two orbits and swaps the other two, which gives the pairing in
column 0 identically. For q ≡ 7 (mod 8), s is even, the blocks cannot be made circulant (every lift of the torus
element to SL(2, q) has order 2s), and the construction runs over ℤ₃ × N_s with negacyclic blocks.

- **Reach.** 59 census values v < 2500 have v = 3(q + 1)/4 with q ≡ 3 (mod 8) a prime power. Seven were already in this
  repository (261, 285, 429, 465, 513, 663, 2475); the other 52 are certified, from 996 = 4 · 249 to 9996 = 4 · 2499.
  879 is one of the census values that Seberry–Yamada's Table 7.1 and Cati–Pasechnik's Table 5 record as known (like
  663, [skew/2652](skew/2652)), so 3516 is a first explicit construction. For q ≡ 7 (mod 8) the orders are 8v or 16v,
  and six of them lie below the smallest power of two in Table 5: 1800 (225(4)), 4344 (543(4)), 4632 (579(5)),
  7272 (909(4)), 7656 (957(4)) and 14352 = 16 · 897 (897(5)). Five more such values (345, 369, 513, 693, 765) are
  already settled at 4v. The audit of October 17 found that Table 5 misses Wallis–Whiteman's 4(q + 1) for q = 1913 ≡ 9
  (mod 16), so 7656 was known (relabelled Known), and the family (|t| + 1)(q + 1) for 765 (8 · 765 with q = 509); and
  that 14352 equals Đoković's skew 276 times our amicable pair of order 52. It also found the pairing in print: Yamada
  (1991), restated as Theorem 1 of Seberry–Yamada (JCMCC 13, 1993), puts Paley's matrix (q ≡ 3 mod 8) in a
  sixteen-block form with A skew and B₁ = −B₀ᵀ. The new orders are therefore labelled New, but obtainable classically:
  every ingredient was published, and we found the combination nowhere.
- **Checks.** `verify.py` rebuilds all 58 matrices from Paley's matrix and checks H Hᵀ = nI and H + Hᵀ = 2I densely
  (the paley-blocks builder gained optional shifts, a column order and negacyclic blocks; the 45 earlier paley-blocks
  certificates up to order 17000 and 25132 rebuild with unchanged SHA-256). `search/skew_3q.py` builds the family from
  the closed formula, independently of `paley_blocks.py`, and passes for all 86 primes q ≡ 3 (mod 4) from 7 to 991.
- **No other multiplier.** For A_c to be skew-type in general, T_c = I and the other three T-matrices must be X, −εXᵀ
  and a skew Y, with 2 PAF(X) + PAF(Y) = (t − 1)δ; the trivial character forces t = 2a² + 1. Exhaustive searches over
  ℤ₉, ℤ₃ × ℤ₃ and ℤ₁₉ found none, so t = 3 is the only multiplier we have.
- **An orthogonal design.** Splitting the product into the identity, the rest of row c of Paley's matrix and the paired
  rows b, d gives a skew-type OD(3(q + 1); 1, q, 2(q + 1)) (the cross terms vanish because the paired part is skew;
  checked for q = 11, 19 and 43). It does not feed Seberry's Theorem 7, whose weight m must divide the order.
- **Twisted T-matrices.** The proof of t = 2a² + 1 uses an integer character of the block ring Z[x]/(x^s + 1), which
  exists only for s odd. For q ≡ 7 (mod 8) (s even) we let each point g of ℤ_t carry a unit ±x^a of the block ring
  as well as an orbit. The skew conditions then force the self-paired orbit to stay empty (a mirror pair at difference
  d contributes twice, the self-mirror pair at ±d/2 once), and everything reduces to one object: T_b supported on a half
  G of ℤ_t (G ⊔ −G = ℤ_t ∖ 0) with unit entries and T_b T_b* = (t − 1)/2 in ℤ[ℤ_t] ⊗ ℤ[x]/(x^s + 1). Units are ≡ 1
  modulo (2, x − 1), so every difference must occur an even number of times in G: t = 3, or t ≡ 3 (mod 8) with G the
  quadratic residues. t = 5 and 7 have no solution (exhaustive, q = 23 and 31); for t = 11 the congruences on the
  exponents are unsolvable for every even s (Smith form 1, 1, 1, 11, 0), and there is not even a complex solution with
  4th, 8th, 12th or 16th roots of unity; t = 19 has none with 8th roots; no half of ℤ₉ or ℤ₃² passes the parity test.
- **No symmetric analogue this way.** T_c = 0 and I on the other fixed orbit make A_c symmetric, but the propus array
  needs two equal blocks among the other three; for q = 19, 43 and 67 no two of them are related by any affine map of ℤ₃
  × ℤ_s. Equal blocks are more than the propus array needs: its rows are orthogonal as soon as A is symmetric and B Bᵀ =
  C Cᵀ, so two blocks with equal Gram matrices would do. Over every placement of δ₀, δ₁, δ₂ (T-matrices of order 3) on
  three of the four orbits, all signs and all relative orbit shifts, q = 19 and 43 give no symmetric block together with
  two equal Grams that passes the dense check (q = 11 gives a symmetric H(36), since s = 3 = t). The automorphisms of
  Paley's matrix that normalise the torus move every orbit, which is why no pair of columns shares a Gram. The same test
  on Paley II (Turyn's Williamson array (X₁, X₂, X₃, X₃), which would give symmetric 6(q + 1) for 39 census values)
  gives only q = 5 (H(36)) among q = 5, 9, 13, 17, 25 and 29. Nor do the new matrices have an involutive
  anti-automorphism P (which would make HP symmetric, as for [symmetric/7956](symmetric/7956)): by nauty (pynauty 2.8.8,
  on the signed digraph of H − I), for q = 11, 19, 43 and 59 they have no anti-automorphism at all, and their
  automorphism groups have order 4, 2, 2 and 2. Nor does any of the 384 Goethals–Seidel variants (order of the last
  three blocks, transposes, signs) for q = 11 and 19; they fall in two classes.
- **No complex structure.** For q ≥ 11 no pairing (A_a ± iA_b, A_c ± iA_d) of the four blocks has X X* + Y Y* real, and
  the element of order 4 behind Paley's complex Hadamard matrices
  ([complex-hadamard/PALEY.md](complex-hadamard/PALEY.md)) cannot pass through T-matrices of order 3: it swaps the
  orbits in two pairs, so T_{k+2} would be the image of T_k under an involution of ℤ₃, whose fixed point would lie in
  two disjoint supports. The same holds for every odd t, since an involution of a set of odd size has a fixed point.
  (An earlier version of this note said that composite t allows fixed-point-free affine involutions; that was wrong.)
  Nor can Paley's complex Hadamard matrix of order (q + 1)/2, a 2 × 2 array of circulants, be multiplied by a pair
  U₀, U₁ of quaternary sequences of odd length t with disjoint supports and U₀U₀* + U₁U₁* = tδ: units are ≡ 1 mod
  (1 + i), so for every shift the number of pairs inside one support would be even, but it is t minus an even number.
- **No iteration.** The new matrices are Goethals–Seidel arrays whose off-diagonal blocks A_jR are back-circulant. A
  sixteen-block structure of the kind the construction needs would be a signed automorphism of order 3s acting with
  four orbits; the automorphism groups above (order 2 for q = 19, 43 and 59) contain only ±I, so the construction does
  not iterate to 9(q + 1).
- **Szekeres's matrices.** For q ≡ 5 (mod 8) a sixteen-block form of a skew matrix of order 2(q + 1) with the pairing
  would give 6(q + 1), 21 census values. The Szekeres matrices from C₀ ∪ C₁ and C₀ ∪ C₃ (cyclotomic classes of order 4)
  have automorphism group of order 2|PSL(2, q)| (nauty, q = 29, 37 and 53), and its elements of order (q + 1)/2 act with
  four orbits, so they are sixteen negacyclic blocks. On the points of GF(q)² \ 0 modulo the squares of GF(q)* they are
  S(u, v) = ±1 as det(u, v) lies in C₀ ∪ C₁ or C₂ ∪ C₃, completed by I ⊗ [[1, 1], [−1, 1]] on the pairs {u, νu}, ν a
  non-square (we checked that this matrix is Hadamard and equivalent to the two-circulant one for q = 13, 29, 37 and
  53). Frobenius (determinant −1 ∈ C₂) reverses S but not the completion, and they have no anti-automorphism; in the
  T(3) test no two blocks of a column pair up (all columns, base points and signs, q = 37), and the off-diagonal blocks
  that are symmetric up to a shift admit only skew T_k, whose row sum 0 forces t = 1. The other Szekeres-type matrices
  (from C₀ ∪ C₂, when they exist) have anti-automorphisms but automorphism groups of order q(q − 1) only, with no
  element of order (q + 1)/2.
