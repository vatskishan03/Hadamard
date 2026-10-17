# Miyamoto's matrix in bordered blocks, times T-matrices

Seberry and Yamada write Miyamoto's Hadamard matrix of order 4q as an M-structure, a 4 × 4 array of bordered blocks
(J. Combin. Math. Combin. Comput. 7 (1990), Theorem 34 and Corollary 35; 1992 survey, Theorem 8.10 and Corollary 8.11).
Their multiplication theorems for M-structures do not cover it: Theorem 4 needs circulant or type-one blocks, and
Theorem 21 needs sixteen matrices T_ij that exist for no odd order t > 1 (each 4 × 4 slice is a signed permutation
matrix, and t of them cannot be orthogonal mod 2). We use that when the Hadamard matrix of order q − 1 in Miyamoto's
construction is of 2N-type (a negaperiodic Golay pair) or made of two circulants (a periodic Golay pair), and Paley's
points are listed along a primitive root, the sixteen blocks commute and are transposed by one permutation. Then
T-matrices of every order t multiply the array, as for Paley's matrix in circulant blocks
([PALEY_BLOCKS.md](PALEY_BLOCKS.md)), and give Hadamard matrices of order 4qt. We found this nowhere.

The route gives 15 orders that Table 1.53 of the Handbook lists as unknown (n = 3279 to 9991), two of them orders that
Đoković claimed only through Seberry–Yamada's Corollary 8.8(i), and four orders that the tables record only through that
corollary ([below](#orders-recorded-only-through-corollary-88i)), among them 9708 = 4 · 2427, which
[FRONTIER.md](../FRONTIER.md) listed as unresolved. With near-Williamson matrices of 2026 in place of the Golay pair
([below](#near-williamson-matrices)) it gives three orders of Cati–Pasechnik's Table 4, n = 1527, 1883 and 2545, the
reopened order n = 2611 and ten Handbook entries.

| Order | n | q | t | negaperiodic Golay pair of length (q − 1)/2 |
|------:|---:|---:|---:|---|
| [13116](13116) | 3279 = 3 · 1093 | 1093 | 3 | Paley, second series (p = 1091) |
| [17516](17516) | 4379 = 29 · 151 | 29 | 151 | Paley, Proposition 6 (p = 13) |
| [18908](18908) | 4727 = 29 · 163 | 29 | 163 | Paley, Proposition 6 (p = 13) |
| [20028](20028) | 5007 = 3 · 1669 | 1669 | 3 | Paley, second series (p = 1667) |
| [28092](28092) | 7023 = 3 · 2341 | 2341 | 3 | Paley, second series (p = 2339) |
| [28716](28716) | 7179 = 3 · 2393 | 2393 | 3 | Golay 26 × Williamson matrices of order 23 |
| [29292](29292) | 7323 = 3 · 2441 | 2441 | 3 | Golay 10 × Paley, Proposition 6 (p = 121 = 11²) |
| [30604](30604) | 7651 = 7 · 1093 | 1093 | 7 | Paley, second series (p = 1091) |
| [33380](33380) | 8345 = 5 · 1669 | 1669 | 5 | Paley, second series (p = 1667) |
| [35780](35780) | 8945 = 5 · 1789 | 1789 | 5 | Paley, second series (p = 1787) |
| [36668](36668) | 9167 = 89 · 103 | 89 | 103 | Paley, Proposition 6 (p = 43) |
| [38660](38660) | 9665 = 5 · 1933 | 1933 | 5 | Paley, second series (p = 1931) |
| [39036](39036) | 9759 = 3 · 3253 | 3253 | 3 | Paley, second series (p = 3251) |
| [39348](39348) | 9837 = 9 · 1093 | 1093 | 9 | Paley, second series (p = 1091) |
| [39964](39964) | 9991 = 97 · 103 | 97 | 103 | Paley, Proposition 6 (p = 47) |

## Construction

**1. Miyamoto's matrix.** Let q ≡ 1 (mod 4) be a prime, m = (q − 1)/2 and g a primitive root mod q. List the points of
Paley's symmetric conference matrix C of order q + 1 as ∞, 0, g^(2j), g^(2j+1) (j ∈ ℤ_m). Then
C(g^(2i), g^(2j)) = χ(g^(2(j−i)) − 1) and the other quadratic-residue blocks are circulant in j as well. Miyamoto's
theorem (J. Combin. Theory Ser. A 57 (1991)) turns C and a Hadamard matrix K of order q − 1 into a Hadamard matrix of
order 4q; `verify/verify.py` implements it (`miyamoto`, in SageMath's block layout). Its sixteen blocks of order q are
(Seberry–Yamada 1990, proof of Theorem 34, with U and V as in their Corollary 35)

    N_rc = σ_rc [[1, ±e], [±eᵀ, J₂ ⊗ U_rc + D₂ ⊗ V_rc]],     D₂ = [[1, −1], [−1, 1]],

with σ_rc = ±1, U_rc ∈ {0, C₁, C₂, C₂ᵀ, C₄} the circulant residue blocks of C and V_rc ∈ {0, I, K_i, K_iᵀ} the
quadrants of K.

**2. A Hadamard matrix of 2N-type.** Take K = [[A, B], [−Bᵀ, Aᵀ]] with A and B negacyclic of order m. It is a
Hadamard matrix exactly when the first rows (a, b) form a negaperiodic Golay pair (NG pair) of length m: NAF_a + NAF_b
= 2m δ (Balonin and Đoković, arXiv:1508.00640, Proposition 2; Ito's "associated pairs"). Then every V_rc is
negacyclic.

**3. The algebra.** Negacyclic matrices commute and are closed under transposition, and so are circulants; J₂ and D₂
are orthogonal (J₂D₂ = 0), so the matrices J₂ ⊗ U + D₂ ⊗ V form a commutative, transpose-closed algebra in which the
all-ones vector e is an eigenvector (row sum r(M) = 2 r(U)). Bordered matrices [[x, u e], [u eᵀ, M]] with x − r(M) =
λu are closed under products (the border of XY is (xw + uν, uy + wμ), equal because w(x − μ) = u(y − ν) = λuw) and
commute; every N_rc has λ = 1. The permutation R_B that fixes the border and sends (h, j) (half h of the core, index j)
to (h, −j) for j = 0 and to (1 − h, −j) for j ≠ 0 satisfies R_B N R_B = Nᵀ: on J₂ ⊗ U it is the reversal j ↦ −j,
which transposes circulants, and on D₂ ⊗ V it is the reversal times diag(1, −1, …, −1), which transposes negacyclic
matrices. The U and V parts never interact (J₂D₂ = 0), so all that is needed of K is that its quadrants lie in a
commutative algebra transposed by a signed reversal j ↦ ±(−j): NG pairs, periodic Golay pairs
([below](#periodic-golay-pairs)), or quasi-Williamson matrices over any abelian group of order (q − 1)/4 through a
bijection that carries inversion to j ↦ −j.

**4. Multiplication by T-matrices.** With T-matrices T₁, …, T₄ of order t put A_j = Σ_k T_k ⊗ N_kj. These are ±1
matrices of order tq in a commutative, transpose-closed algebra, R = R_t ⊗ R_B (R_t: x ↦ −x on ℤ_t) satisfies
R A_j R = A_jᵀ, and Σ_j A_j A_jᵀ = Σ_k T_k T_kᵀ ⊗ 4q I = 4qt I by the orthogonality of the block rows of Miyamoto's
matrix. The Goethals–Seidel array

    [[A₁, A₂R, A₃R, A₄R], [−A₂R, A₁, A₄ᵀR, −A₃ᵀR], [−A₃R, −A₄ᵀR, A₁, A₂ᵀR], [−A₄R, A₃ᵀR, −A₂ᵀR, A₁]]

is a Hadamard matrix of order 4qt. This is Seberry and Yamada's Theorem 4 (1990) over the algebra of bordered blocks
instead of a group ring; the Paley-block route is the case of circulant blocks. Prime powers q work the same way
(checked for q = 9, 25, 49, 81 during the audit); the program implements prime q.

**NG pairs used.** All from Balonin–Đoković (arXiv:1508.00640) and checked by the program and by `verify/verify.py`:
- *Paley pairs of length p + 1* (their Proposition 6; for p ≡ 1 (mod 4) their first Paley series): for an odd prime
  power p, a negacyclic conference matrix of order p + 1 with first row (0, c₁, …, c_p) gives the pair
  (1, c₁, …, c_p), (−1, c₁, …, c_p) of length p + 1. We build the matrix on PG(1, p) along a primitive element w of
  GF(p²): c_k = χ(Δ(1, w^k)/Δ(1, w)) with Δ(x, y) = x y^p − x^p y.
- *Second Paley series* (Proposition 8; it contains Ito's series): for p ≡ 3 (mod 4), (1, c₂, c₄, …, c_(p−1)) and
  (c₁, c₃, …, c_p) of length (p + 1)/2.
- *Williamson matrices* of odd order w through their array (29), read in ℤ[y]/(y^(2w) + 1) ≅ ℤ[x]/(x² + 1) ⊗ ℤ[ℤ_w]
  (y ↦ xz): length 2w. For 7179 the symmetric circulant Williamson matrices of order 23 are found by the program, and
  for 12260 and 26972 those of order 153, constant on the orbits of a multiplier group of order 16 (`wmult:153:145,35`;
  Williamson matrices of order 153 = 17 · 18/2 are classical, key w4 of Seberry–Yamada's Table A.1).
- *Turyn's multiplication* by Golay pairs of length 10 and 26 (their Proposition 3, formulas (15) and (16)).

## Periodic Golay pairs

If (a, b) is a periodic Golay pair of length m (PAF_a + PAF_b = 2m δ), K = [[A, B], [−Bᵀ, Aᵀ]] with A, B circulant is a
Hadamard matrix of order q − 1, every V_rc is circulant, and the plain reversal (h, j) ↦ (h, −j) (border fixed)
transposes every block, since J₂ and D₂ are symmetric. Steps 3 and 4 then hold verbatim with R_B replaced by this
reversal. Formulas (15) and (16) of Balonin and Đoković multiply a Golay pair of length g and a periodic Golay pair of
length v into a periodic Golay pair of length gv (Georgiou, Stylianou, Drosou and Koukouvinos 2014, Theorem 1).

For q = 809, m = 404 = 2 · 202, and Đoković and Kotsireas (Numer. Algorithms 69 (2015), arXiv:1310.5773) give a periodic
Golay pair of length 202 as unions of orbits of H = {1, 87, 95, 137, 185} on ℤ₂₀₂. No negaperiodic Golay pair of length
404 is known (Balonin–Đoković, list (31)). The variant was found during the novelty audit of this route.

## Near-Williamson matrices

Let A, B, C, D be circulant ±1 matrices of order w with B, C, D symmetric and AAᵀ + BBᵀ + CCᵀ + DDᵀ = 4wI:
near-Williamson matrices (Kharaghani, Mohammadian, Tayfeh-Rezaie and Zaitsev, arXiv:2605.08661v2, September 2026). Good
matrices (A of skew type) and Williamson matrices (A symmetric) are special cases. With R the back-identity, X₁ = AR is
symmetric, and X₁, X₂ = B, X₃ = C, X₄ = D commute pairwise: AR · B = ABR = B · AR because R B R = B. So their Williamson
array K = [[X₁, X₂, X₃, X₄], [−X₂, X₁, −X₄, X₃], [−X₃, X₄, X₁, −X₂], [−X₄, −X₃, X₂, X₁]] is a Hadamard matrix of order
4w, and each of its four quadrants of order 2w has the form [[X, Y], [−Y, X]] with X, Y in the commutative algebra
generated by X₁, …, X₄. These quadrants form a commutative algebra, closed under transposition ([[X, Y], [−Y, X]]ᵀ =
[[X, −Y], [Y, X]]), and the half swap (h, y) ↦ (1 − h, y) transposes every one of them.

For q = 4w + 1 a prime, take this K in Miyamoto's array (m = (q − 1)/2 = 2w). Relabel the indices of each quadrant by
the pairs {j, 1 − j} of ℤ_m, (0, y) ↦ j and (1, y) ↦ 1 − j, so that the half swap becomes the reflection j ↦ 1 − j. That
reflection also transposes the circulant blocks from Paley's matrix, and the U and V parts never interact, so the
permutation R_B that fixes the border and sends (b, j) to (b, 1 − j) transposes all sixteen blocks, which commute. Steps
3 and 4 then hold as before. `search/miyamoto_blocks.py` checks the commutation and the transposition for q = 29, 37 and
53 with near-Williamson matrices whose A is not symmetric (`selftest`).

The idea came from a judge agent of this repository's frontier workflow, which built H(7532) with near-Williamson
matrices of order 67 and its own implementation; the certificates here were built again by `miyamoto_blocks.py`.

| Order | n | q | t | near-Williamson matrices of order (q − 1)/4 | Status |
|------:|---:|---:|---:|---|---|
| [6108](6108) | 1527 = 3 · 509 | 509 | 3 | 127, database of arXiv:2605.08661v2 | First explicit construction: Table 4 lists it as open; Đoković claimed it through Corollary 8.8(i) |
| [7532](7532) | 1883 = 7 · 269 | 269 | 7 | 67, Table 7 of arXiv:2605.08661v2 | **New**: Table 4 lists it as open |
| [10180](10180) | 2545 = 5 · 509 | 509 | 5 | 127, database of arXiv:2605.08661v2 | First explicit construction: Table 4 lists it as open; Đoković claimed it through Corollary 8.8(i) |
| [10444](10444) | 2611 = 7 · 373 | 373 | 7 | 93, database of arXiv:2605.08661v2 | First explicit construction: omitted by Table 4 through Seberry–Yamada Table A.2 (Corollary 8.8(i)) |

The database also holds orders 73, 79 and 97 with q = 4w + 1 a prime (293, 317, 389). Over the Handbook entries 3001 ≤ n
≤ 9999 with no valid construction, the route reaches ten more, all built (T-matrices of order 37 from
`search/tmat_mult.py 37 10`, [t-matrices/T37.txt](../t-matrices/T37.txt)):

| Order | n | q | t | Status |
|------:|---:|---:|---:|---|
| [13988](13988) | 3497 = 13 · 269 | 269 | 13 | New, beyond the tables (Handbook t = 3) |
| [18292](18292) | 4573 = 17 · 269 | 269 | 17 | New, beyond the tables (t = 3) |
| [25364](25364) | 6341 = 17 · 373 | 373 | 17 | First explicit construction: Đoković's set Δ, through Corollary 8.8(i) (Handbook t = 8) |
| [26468](26468) | 6617 = 13 · 509 | 509 | 13 | First explicit construction: Đoković 2010, through Corollary 8.8(i) (t = 4) |
| [33356](33356) | 8339 = 31 · 269 | 269 | 31 | New, beyond the tables (t = 3) |
| [34316](34316) | 8579 = 23 · 373 | 373 | 23 | First explicit construction: Đoković's set Δ (t = 3) |
| [34612](34612) | 8653 = 17 · 509 | 509 | 17 | First explicit construction: Đoković 2010 (t = 3) |
| [38684](38684) | 9671 = 19 · 509 | 509 | 19 | First explicit construction: Đoković 2010 (t = 3) |
| [39308](39308) | 9827 = 31 · 317 | 317 | 31 | New, beyond the tables (t = 3) |
| [39812](39812) | 9953 = 37 · 269 | 269 | 37 | New, beyond the tables (t = 3) |

The orders with q = 269 and 317 need the near-Williamson matrices of 2026. For q = 509, Đoković's good matrices of order
127 (J. Combin. Math. Combin. Comput. 14 (1993)) satisfy the same hypotheses, so those orders were obtainable
classically, by a combination we found nowhere. For q = 373 we know no classical input of this kind (Table A.1 records
Williamson-type matrices of order 93, but not symmetric ones).

## Orders recorded only through Corollary 8.8(i)

| Order | n | q | t | pair of length (q − 1)/2 | Recorded by |
|------:|---:|---:|---:|---|---|
| [9708](9708) | 2427 = 3 · 809 | 809 | 3 | periodic: Golay 2 × Đoković–Kotsireas PG(202) | Seberry–Yamada Table A.2; omitted by Cati–Pasechnik Table 4 |
| [12260](12260) | 3065 = 5 · 613 | 613 | 5 | negaperiodic: Williamson matrices of order 153 | Handbook Table 1.53, t = 2 |
| [26972](26972) | 6743 = 11 · 613 | 613 | 11 | negaperiodic: Williamson matrices of order 153 | Đoković's set Δ (arXiv:0912.5091) |
| [35596](35596) | 8899 = 11 · 809 | 809 | 11 | periodic: Golay 2 × Đoković–Kotsireas PG(202) | Đoković's set Δ (arXiv:0912.5091) |

The tables obtain these orders from Williamson-type matrices of order 809 or 613 that Seberry–Yamada's Table A.1
obtains only from Corollary 8.8(i), whose printed proof fails ([NOTES.md](../NOTES.md) §7). They are labelled First
explicit construction ([PALEY_BLOCKS.md](PALEY_BLOCKS.md#orders-recorded-only-through-corollary-88i) has the other
orders of this kind).

## Reach

Over the 691 entries 3001 ≤ n ≤ 9999 of Handbook Table 1.53 that had no valid construction
([PALEY_BLOCKS.md](PALEY_BLOCKS.md#beyond-table-4)), the route needs n = qt with q ≡ 1 (mod 4) a prime power, T-matrices
of order t, and a pair of length (q − 1)/2: an NG pair from the families above (lengths g · v with g a Golay number and
v − 1 or 2v − 1 a prime power, or v = 2w with Williamson or quasi-Williamson matrices of order w known), or a periodic
Golay pair (lengths g · v, v ∈ {1, 34, 50, 58, 72, 74, 82, 122, 202, 226}). With q prime and NG pairs from Paley's
series, Golay products and the Williamson orders of Balonin–Đoković it reaches 27 entries: 11 are also reached by the
Paley-block route and are in this repository already, 6031 is in it by Yamada–Kiyasu ([24124](24124)), and the other 15
are the orders of the first table. Prime powers q = 25 and 49 add 4075 and 6419, Paley-block orders already here.
Whiteman's Williamson orders p(p + 1)/2 add 3065 and 6743 (w = 153) and 3903 and 9107 (w = 325), and periodic Golay
pairs add 8899; 3065, 6743 and 8899 are built (table above), 3903 and 9107 are not yet (no Williamson matrices of order
325 found by the program so far; those of Table A.1, key w4, from J. Wallis and Whiteman, are not implemented here).

Among the FRONTIER entries the Golay-pair version reaches one, 2427, with the periodic Golay pair of length 404, and the
near-Williamson version reaches 1527, 1883 and 2545 and the reopened 2611 ([above](#near-williamson-matrices)). Of the
reopened orders of Corollary 8.8(i), 2571 = 3 · 857 remains: it would need near-Williamson matrices of order 214 or a
negaperiodic Golay pair of length 428. 2823 = 3 · 941 would need either of order 235 or length 470. Periodic Golay pairs
of length 134, 186, 254, 428 or 470 cannot exist, since 2m is not a sum of two squares.

## Status

- Of the 11 orders of the first table with T-matrices of order 3, 5, 7 or 9, nine are **New, beyond the tables, but
  obtainable classically**: every ingredient was published before 2026 (Paley, Miyamoto 1991, the Paley NG pairs,
  Williamson matrices of order 23, Golay pairs, Turyn's T-matrices), and only their combination is new. The other two,
  7179 and 7323, were claimed by Đoković (arXiv:1008.2043, Proposition 3.2) through Corollary 8.8(i) only
  ([NOTES.md](../NOTES.md) §7) and are labelled **First explicit construction**.
- The 4 orders with T-matrices of order 103, 151 or 163 (17516, 18908, 36668, 39964) are **New, beyond the tables**.
  Cooper–Wallis with the same T-matrices and Williamson-type matrices of order 29, 89 or 97 also gives them (Table A.1
  of 1992; it obtains those of order 89 through Corollary 8.8(i), but Corollary 8.8(ii) with a symmetric Hadamard
  matrix of order 44 gives them validly).
- 9708, 12260, 26972 and 35596 are **First explicit constructions**: tables record them as known, through Corollary
  8.8(i) only.
- 6108, 10180 and 10444 are **First explicit constructions** and 7532 is **New**. 7532 and 10444 need the
  near-Williamson matrices of 2026; 6108 and 10180 were obtainable classically with Đoković's good matrices of order
  127. The ten Handbook entries of the second near-Williamson table are labelled as there.
- SageMath's development version (October 17, 2026) constructs none of the 19. It returns Unknown for 14 and raises
  UnknownError in its Miyamoto test for 12260, 33380, 35780, 38660 and 39348 (no symmetric conference matrix of order
  n + 1 is known to it; with that Unknown treated as False it returns Unknown).

The sources searched are listed in [NOVELTY.md](../NOVELTY.md).

## Reproduction and verification

```
python3 search/miyamoto_blocks.py selftest                    # q = 13 ... 61 with T-matrices of order 3, 5, 7
python3 search/miyamoto_blocks.py certify 1093 3 paley2:1091   # writes hadamard/13116/certificate.json
python3 search/miyamoto_blocks.py certify 809 3 golay2:dk202   # periodic Golay pair: hadamard/9708
python3 search/miyamoto_blocks.py certify 613 5 wmult:153:145,35   # Williamson matrices of order 153: hadamard/12260
python3 search/miyamoto_blocks.py certify 269 7 nw:near-williamson/NW67.txt   # near-Williamson: hadamard/7532
python3 verify/verify.py 13116                                # rebuilds the matrix and checks H Hᵀ = 4n I
```

The certificates store q, g, the NG pair ("ng"), the periodic Golay pair ("pg") or the near-Williamson matrices ("nw")
and the T-matrices. `verify/verify.py` (construction `miyamoto-blocks`) checks the pair and the T-matrices, rebuilds
Paley's matrix, the blocks, the reversal and the Goethals–Seidel array independently of the search program (Miyamoto's
array with the same `miyamoto` function, which the search program imports from `verify/verify.py`), and checks H Hᵀ = 4n
I densely. An independent re-implementation written for the audit reproduces the SHA-256 of 13116, 17516 and 29292 and
the dense check of 9708.
