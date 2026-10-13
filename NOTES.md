# Notes on structured searches (October 2026)

What we learned while searching for the unresolved orders in [FRONTIER.md](FRONTIER.md), including the
searches that found nothing. Section 1 is a mechanism, sections 2–4 are negative results, section 5 is the
literature as of October 2026.

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
