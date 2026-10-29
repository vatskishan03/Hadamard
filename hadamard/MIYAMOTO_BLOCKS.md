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
corollary ([below](#orders-recorded-only-through-corollary-88i)), among them 9708 = 4 · 2427. With near-Williamson
matrices of 2026 in place of the Golay pair ([below](#near-williamson-matrices)) it gives three orders of
Cati–Pasechnik's Table 4, n = 1527, 1883 and 2545, the reopened order n = 2611, nine Handbook entries and a new matrix
of the known order 9827. With Seberry and Yamada's symmetric Williamson-type matrices
([below](#symmetric-williamson-type-matrices)) it gives sixteen Handbook entries that the tables record only through
Corollary 8.8(i), and, with a symmetric conference matrix of order 118 found here, n = 2823 of Table 4
([11292](11292), [below](#reach)).

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
(checked for q = 9, 25, 49, 81); the program implements prime q.

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
404 is known (Balonin–Đoković, list (31)).

## Near-Williamson matrices

Let A, B, C, D be circulant ±1 matrices of order w with B, C, D symmetric and AAᵀ + BBᵀ + CCᵀ + DDᵀ = 4wI:
near-Williamson matrices (Kharaghani, Mohammadian, Tayfeh-Rezaie and Zaitsev, arXiv:2605.08661v2, September 2026). Good
matrices (A of skew type) and Williamson matrices (A symmetric) are special cases. With R the back-identity, X₁ = AR is
symmetric, and X₁, X₂ = B, X₃ = C, X₄ = D commute pairwise: AR · B = ABR = B · AR because R B R = B. So their Williamson
array K = [[X₁, X₂, X₃, X₄], [−X₂, X₁, −X₄, X₃], [−X₃, X₄, X₁, −X₂], [−X₄, −X₃, X₂, X₁]] is a Hadamard matrix of order
4w, and each of its four quadrants of order 2w has the form [[X, Y], [−Y, X]] with X, Y in the commutative algebra
generated by X₁, …, X₄. These quadrants lie in the commutative algebra {[[X, Y], [−Y, X]] : X, Y ∈ alg(X₁, …, X₄)},
which is closed under transposition ([[X, Y], [−Y, X]]ᵀ = [[X, −Y], [Y, X]]), and the half swap (h, y) ↦ (1 − h, y)
transposes every one of them.

For q = 4w + 1 a prime, take this K in Miyamoto's array (m = (q − 1)/2 = 2w, even, so that j ↦ 1 − j has no fixed point
on ℤ_m). Relabel the indices of each quadrant by the pairs {j, 1 − j} of ℤ_m, (0, y) ↦ j and (1, y) ↦ 1 − j, so that the
half swap becomes the reflection j ↦ 1 − j. That reflection also transposes the circulant blocks from Paley's matrix,
and the U and V parts never interact, so the permutation R_B that fixes the border and sends (b, j) to (b, 1 − j)
transposes all sixteen blocks, which commute. Steps 3 and 4 then hold as before (the border condition uses only the row
sums r(C₁) = r(C₄) = 1 and r(C₂) = 0 of Paley's residue blocks, from the cyclotomic numbers of order 2 for q ≡ 1 mod 4).
`search/miyamoto_blocks.py` checks the commutation and the transposition for q = 29, 37 and 53 with near-Williamson
matrices whose A is not symmetric (`selftest`).

The certificates here are built by `miyamoto_blocks.py` and checked by `verify/verify.py`.

| Order | n | q | t | near-Williamson matrices of order (q − 1)/4 | Status |
|------:|---:|---:|---:|---|---|
| [6108](6108) | 1527 = 3 · 509 | 509 | 3 | 127, database of arXiv:2605.08661v2 | First explicit construction: Table 4 lists it as open; Đoković claimed it through Corollary 8.8(i) |
| [7532](7532) | 1883 = 7 · 269 | 269 | 7 | 67, Table 7 of arXiv:2605.08661v2 | **New**: Table 4 lists it as open |
| [10180](10180) | 2545 = 5 · 509 | 509 | 5 | 127, database of arXiv:2605.08661v2 | First explicit construction: Table 4 lists it as open; Đoković claimed it through Corollary 8.8(i) |
| [10444](10444) | 2611 = 7 · 373 | 373 | 7 | 93, database of arXiv:2605.08661v2 | First explicit construction: omitted by Table 4 through Seberry–Yamada Table A.2 (Corollary 8.8(i)) |

The database also holds orders 73, 79 and 97 with q = 4w + 1 a prime (293, 317, 389). Over the Handbook entries 3001 ≤ n
≤ 9999, the route reaches ten more, all built (T-matrices of order 37 from
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
| [39308](39308) | 9827 = 31 · 317 | 317 | 31 | Known: Đoković's set Δ (Cooper–Wallis with Williamson-type matrices of order 317, Corollary 8.8(ii)); independent matrix |
| [39812](39812) | 9953 = 37 · 269 | 269 | 37 | New, beyond the tables (t = 3) |

The orders with q = 269 need the near-Williamson matrices of 2026. For q = 317, Turyn's Williamson matrices of order 79
= (157 + 1)/2 (1972) satisfy the same hypotheses, and 9827 was already known (Đoković's set Δ, through Corollary 8.8(ii)
with Paley's conference matrix of order 158). For q = 509, Đoković's good matrices of order 127 (J. Combin. Math.
Combin. Comput. 14 (1993)) satisfy the same hypotheses, so those orders were obtainable classically, by a combination we
found nowhere. For q = 373 we know no classical input of this kind: J. Wallis's Williamson-type matrices of order 93
(Linear and Multilinear Algebra 3 (1975)) are amicable, but A₁ is not symmetric and does not commute with B₁ or D₁, and
the symmetric propus Hadamard matrices of orders 372 (Đoković, arXiv:2301.02751) and 268 (Balonin–Đoković,
arXiv:1803.08787) have quadrants that do not commute (checked).

## Symmetric Williamson-type matrices

The argument for near-Williamson matrices uses only that X₁, …, X₄ are symmetric, commute pairwise and satisfy X₁² +
X₂² + X₃² + X₄² = 4wI. For symmetric matrices, amicable means commuting, so any symmetric Williamson-type matrices of
order w = (q − 1)/4 give K, with the same relabelling and the same R_B. Seberry and Yamada (J. Combin. Math. Combin.
Comput. 7 (1990)) construct many:

- Corollary 25: order 2uv from orders u and v, symmetric if the inputs are;
- Corollary 27: order p + 2 (p ≡ 1 mod 4 a prime power, with a skew or symmetric core of order (p + 1)/2 or (p + 3)/2);
- Corollary 29(ii): order q from a symmetric conference matrix or a symmetric Hadamard matrix of order (q − 1)/2;
- Corollary 31: order 2p + 1 when 2p − 1 is a prime power.

Turyn's Williamson matrices of order (p + 1)/2 (p ≡ 1 mod 4 a prime power) and circulant Williamson matrices serve as
seeds. `search/symwt.py` builds and checks all of these (Yamada's circulants α, β with α² + β² = pI over GF(p²) for
prime powers p, Paley's conference matrices over GF(p^k)). `miyamoto_blocks.py certify q t "wt:RECIPE"` stores the four
matrices in the certificate, and `verify.py` checks them before it builds K (`miyamoto_wt_k`).

Over the Handbook entries 3001 ≤ n ≤ 9999 that Table 1.53 lists with t ≥ 3 or that the tables record only through
Corollary 8.8(i), the route reaches 148 (with prime powers in Corollary 29(ii) and in Turyn's family). 124 already have
valid constructions (Đoković's set Δ of arXiv:0912.5091 and his Propositions 3.1–3.4 of arXiv:1008.2043, by routes
checked in [NOTES.md](../NOTES.md) §7, or orders built here before), and 8 that rest on the corollary were built here
by other routes (6341, 6617, 7179, 7323, 8079, 8579, 8653 and 9671). The other sixteen are built here:

| Order | n | q | t | symmetric Williamson-type matrices of order (q − 1)/4 | Recorded by |
|------:|---:|---:|---:|---|---|
| [14316](14316) | 3579 = 3 · 1193 | 1193 | 3 | 298: Corollary 25 from orders 1 and 149 (Corollary 29(ii), conference matrix of order 74) | Handbook Table 1.53, t = 2 |
| [15460](15460) | 3865 = 5 · 773 | 773 | 5 | 193: Corollary 29(ii), symmetric Hadamard matrix of order 96 | Handbook Table 1.53, t = 2 |
| [19060](19060) | 4765 = 5 · 953 | 953 | 5 | 238: Corollary 25 from Williamson matrices of orders 7 (Turyn) and 17 | Đoković's set Δ |
| [20796](20796) | 5199 = 3 · 1733 | 1733 | 3 | 433: Corollary 29(ii), symmetric Hadamard matrix of order 216 | Đoković's set Δ |
| [21644](21644) | 5411 = 7 · 773 | 773 | 7 | 193, as for 15460 | Đoković's set Δ |
| [22180](22180) | 5545 = 5 · 1109 | 1109 | 5 | 277: Corollary 29(ii), conference matrix of order 138 | Đoković's set Δ |
| [23860](23860) | 5965 = 5 · 1193 | 1193 | 5 | 298, as for 14316 | Đoković's set Δ |
| [23916](23916) | 5979 = 3 · 1993 | 1993 | 3 | 498: Corollary 25 from orders 3 (Turyn) and 83 (Corollary 31, p = 41) | Đoković's set Δ |
| [32748](32748) | 8187 = 3 · 2729 | 2729 | 3 | 682: Corollary 25 from Williamson matrices of orders 11 and 31 (Turyn) | Đoković 2010, Proposition 3.2 |
| [33860](33860) | 8465 = 5 · 1693 | 1693 | 5 | 423: Corollary 27(i), p = 421 | Đoković's set Δ |
| [34012](34012) | 8503 = 11 · 773 | 773 | 11 | 193, as for 15460 | Đoković's set Δ |
| [34660](34660) | 8665 = 5 · 1733 | 1733 | 5 | 433, as for 20796 | Đoković's set Δ |
| [35004](35004) | 8751 = 3 · 2917 | 2917 | 3 | 729: Corollary 29(ii) over GF(729), symmetric Hadamard matrix of order 364 | Đoković 2010, Proposition 3.3 |
| [35436](35436) | 8859 = 3 · 2953 | 2953 | 3 | 738: Corollary 25 from Williamson matrices of orders 9 (Turyn) and 41 | Đoković 2010, Proposition 3.2 |
| [38508](38508) | 9627 = 3 · 3209 | 3209 | 3 | 802: Corollary 25 from orders 1 and 401 (Corollary 29(ii), symmetric Hadamard matrix of order 200 over GF(49)) | Đoković 2010, Proposition 3.2 (printed as 9087) |
| [39860](39860) | 9965 = 5 · 1993 | 1993 | 5 | 498, as for 23916 | Đoković's set Δ |

These tables obtain each order from Williamson-type matrices of order q through Corollary 8.8(i): for q < 2000 by
Table A.1's key wl, which Seberry and Yamada credit to Miyamoto (1991), with Williamson-type matrices of order (q − 1)/4
(for 953, 1193 and 1993 those of orders 238, 298 and 498 from key wo), for 2917 by Đoković's Proposition 3.3 with
Williamson-type matrices of order 729, and for 2729, 2953 and 3209 by his Proposition 3.2 with an Hadamard matrix of
order (q − 1)/2. The printed proof of that corollary fails
([NOTES.md](../NOTES.md) §7). Miyamoto's paper, which we could not read, defines matrices of symmetric block form,
[[A, B], [−B, A]] with A and B symmetric; his corollary probably assumes an Hadamard matrix of order q − 1 of that form,
and the failing step may be Seberry and Yamada's restatement. All sixteen are labelled First explicit construction.
Every ingredient here was published before 2026; the combination is ours. All sixteen pass the dense check of
`verify.py` (H Hᵀ = nI, SHA-256 as in the READMEs). An independent rebuild (own Miyamoto array, relabelling and
T-matrices) checks H Hᵀ = nI densely for 14316, 15460 and 35004 and reproduces their SHA-256 in the repository's
convention.

**The corollary with valid inputs.** Lemma 26 (Seberry–Yamada's Corollary 29) needs symmetric, commuting ±1 matrices
M₃, M₄ of order (q − 1)/2 with M₃² + M₄² = (q − 1)I, that is a symmetric quaternary Hadamard matrix (M₃ + M₄)/2 +
i(M₄ − M₃)/2 ([NOTES.md](../NOTES.md) §19). Paley's symmetric conference matrices give one (C ± I, Turyn), symmetric
Hadamard matrices give one (H, H), and so do symmetric Williamson-type matrices (Y, Z, X, X) with a repeated block,
arranged as M₃ = [[Y, X], [X, −Z]], M₄ = [[Z, X], [X, −Y]]. Such propus-type inputs come out of Turyn's family and of
Corollaries 27, 29(ii) with M₃ = M₄ and 31 (Seberry–Balonin, arXiv:1512.01732, Lemmas 1, 3 and 4). The symmetric
quaternary Hadamard matrices are closed under Kronecker products. So symmetric Williamson-type matrices of order q
exist, and Cooper–Wallis gives the order as well, for

- q = 773, 1693, 1733 and 2917, from repeated-block inputs of orders 193, 423, 433 and 729 (Corollaries 29(ii) with
  M₃ = M₄ = H and 27(i)): 3865, 5411, 8503, 8465, 5199, 8665 and 8751 (also q = 293, 1453 and 1637, from orders 73, 363
  and 409);
- q = 953, 1993, 2729, 2953 and 3209, from Kronecker products: (q − 1)/2 = 476 = 14 · 34, 996 = 6 · 166, 1364 = 22 · 62,
  1476 = 18 · 82 and 1604 = 2 · 802, with Paley conference matrices of orders 14, 6, 62, 18 and 82, Sylvester's matrix
  of order 2 and the repeated-block arrangements of orders 34, 166, 22 and 802 (Williamson-type matrices of orders 17,
  83, 11 and 401): 4765, 5979, 9965, 8187, 8859 and 9627. The same gives q = 521 (260 = 10 · 26).

`symwt.py "cor29r(q, …)"` and `"cor29q(q, …)"` check all of them. We found the repeated-block arrangement in none of the
sources we could read; Miyamoto's paper may contain the equivalent statement through its symmetric block form. Only q =
1109 (order 5545) has no Williamson-type route that we know: Lemma 26 would need a symmetric quaternary Hadamard matrix
of order 554 = 2 · 277. For q = 1193 (orders 3579 and 5965) the symmetric Hadamard matrix of order 596
([symmetric/596](../symmetric/596)) gives one through Corollary 8.8(ii) (`symwt.py "cor29h(1193, …)"`, checked exactly).

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
pairs add 8899; 3065, 6743 and 8899 are built (table above), 3903 and 9107 are not built (the program has found no
Williamson matrices of order 325, and those of Table A.1, key w4, from J. Wallis and Whiteman, are not implemented
here).

Among the FRONTIER entries the Golay-pair version reaches one, 2427, with the periodic Golay pair of length 404, and the
near-Williamson version reaches 1527, 1883 and 2545 and the reopened 2611 ([above](#near-williamson-matrices)). Of the
reopened orders of Corollary 8.8(i), 2571 = 3 · 857 remains: it would need near-Williamson matrices of order 214 (in the
sense used here, B, C, D symmetric; Kharaghani et al. define them for odd orders only) or a negaperiodic Golay pair of
length 428. Periodic Golay pairs of length 134, 186, 254, 428 or 470 cannot exist, since 2m is not a sum of two
squares. Symmetric Williamson-type matrices of order 214 would do as well; the constructions of
[the section above](#symmetric-williamson-type-matrices) do not give them: 214 = 2 · 107 needs order 107, a prime ≡ 3
(mod 4) that is not p + 2 or (p + 1)/2 for a prime power p ≡ 1 (mod 4) (105 and 213 are not), nor 2p + 1 with 2p − 1 a
prime power (2 · 53 − 1 = 105).

**Order 2823.** The order 2823 = 3 · 941 needs symmetric Williamson-type matrices of order 235 or a negaperiodic Golay
pair of length 470. By Corollary 27, 235 = 233 + 2 needs (233 + 1)/2 = 117 to be a prime power or a symmetric
conference matrix of order 118. Brouwer's table of strongly regular graphs (read on 22 October 2026) lists the
conference graph SRG(117, 58, 28, 29), and so such a matrix, as existing (credited to Jose Yip Zheng Ho). Egan,
Kotsireas, London and Ó Catháin (arXiv:2610.27175) report that it was communicated to Brouwer in October 2026 and is
not publicly available. We found one independently by an exact cover
([search/conference118](../search/conference118/README.md), [search/NOTES.md](../search/NOTES.md) §95).

With it, Corollary 27 (U = (I, 0, Q, Q), V = (M, N, I, I), Yamada's circulants for q = 233) gives symmetric
Williamson-type matrices of order 235 with a repeated block (`symwt.py 'cor27c(233, …)'`). The Miyamoto blocks with
T-matrices of order 3 then give [11292](11292) (n = 2823), and Seberry–Balonin's Corollary 1 gives the symmetric
Hadamard matrix of order 940 ([symmetric/940](../symmetric/940)). The repeated block also gives symmetric
Williamson-type matrices of order 941 = 4 · 235 + 1 (`cor29r`), hence Cooper–Wallis orders 941 · t.

**Orders 4 · 941 · t.** The blocks of q = 941 with the symmetric Williamson-type matrices of order 235 take T-matrices
of every order t and give Hadamard matrices of order 4 · 941 · t (n = 941t); t = 3 is [11292](11292). Of the odd t from
5 to 53, t = 5, 7, 9, 15, 21, 25, 27, 29, 35, 37 and 49 give orders that Paley I or II or the Paley-block route reach
already (checked), and t = 39, 45 and 51 have no T-matrices here. The other eleven give orders above Table 1.53 of the
Handbook (n ≥ 10351 > 9999). Paley I and II, Miyamoto and Yamada–Kiyasu do not reach them. The Paley-block route fails
or needs T-matrices of order 941 or a multiple of it, which are not known to us. Cooper–Wallis reaches them only
through the repeated-block matrices of order 941 above, that is, through the same conference matrix of order 118. All
eleven are certified by

```
python3 search/miyamoto_blocks.py certify 941 t 'wt:cor27c(233, "conference118_5fe566972f21790b.txt")' t-matrices/Tt.txt
```

(t = 13: built-in T-matrices, no file). The range stops at order 200000 (t = 53). T-matrices of order 55 and above give
larger orders in the same way, which we did not examine.

| Order | n | t | T-matrices | SHA-256 of the matrix | Check |
|------:|---|---:|---|---|---|
| 41404 | 10351 = 11 · 941 | 11 | [T11](../t-matrices/T11.txt): Golay pair of length 10 | `68a4d200…` | H Hᵀ = 4nI checked |
| 48932 | 12233 = 13 · 941 | 13 | built in (base sequences BS(7, 6)) | `096d5a5c…` | H Hᵀ = 4nI checked |
| 63988 | 15997 = 17 · 941 | 17 | [T17](../t-matrices/T17.txt): Golay pair of length 16 | `e15b6cce…` | H Hᵀ = 4nI checked |
| 71516 | 17879 = 19 · 941 | 19 | [T19](../t-matrices/T19.txt): `tmat_mult.py 19 7` | `fc3f3d9f…` | H Hᵀ = 4nI checked |
| 86572 | 21643 = 23 · 941 | 23 | [T23](../t-matrices/T23.txt): SageMath | `d088897f…` | not densely checked |
| 116684 | 29171 = 31 · 941 | 31 | [T31](../t-matrices/T31.txt): `tmat_mult.py 31 5` | `57700063…` | not densely checked |
| 124212 | 31053 = 33 · 941 | 33 | [T33](../t-matrices/T33.txt): Golay pair of length 32 | `b82c4aac…` | not densely checked |
| 154324 | 38581 = 41 · 941 | 41 | [T41](../t-matrices/T41.txt): `tmat_mult.py 41 10` | `877c7d91…` | not densely checked |
| 161852 | 40463 = 43 · 941 | 43 | [T43](../t-matrices/T43.txt): `tmat_mult.py 43 6` | `403b581d…` | not densely checked |
| 176908 | 44227 = 47 · 941 | 47 | [T47](../t-matrices/T47.txt): SageMath | `28735ed8…` | not densely checked |
| 199492 | 49873 = 53 · 941 | 53 | [T53](../t-matrices/T53.txt): SageMath | `d9742c3a…` | not densely checked |

The standard `certify` and `verify.py` assemble the matrix in int8 and need about 28n² bytes: 3.0 GB for 41404 and 70 GB
for 199492. The first four were checked by running the standard checks (the four matrices of order 235, the T-matrices,
`gram_is_scalar`, the SHA-256 of the canonical CSV) on a row-on-demand version of the Goethals–Seidel array, within
3.5 GB (581 s to 5282 s on one core). Its rows equal the standard builder's for t = 1, 3, 5 and 7, and it reproduces the
SHA-256 of 11292.

The other seven were certified by the same command, with the SHA-256 computed by streaming, but not checked densely: the
Gram check costs about N³ operations (estimated 2.5 hours for 86572, more than 18 hours for 199492 on one core under
load). For 86572 a second row generator reproduces the SHA-256 and finds all inner products among 1017 random and
boundary rows zero. None of the eleven is stored in this repository (the standard check of the smallest needs 3 GB); the
command above rebuilds each certificate.

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
  127. The ten Handbook entries of the second near-Williamson table are labelled as there (39308 is Known: Đoković's set
  Δ).
- The sixteen orders of the symmetric Williamson-type table are **First explicit constructions**. For thirteen of them
  (all but 3579, 5545 and 5965, that is q = 1109 and 1193) Cooper–Wallis also works once the corollary is used with a
  repeated block (above). For q = 1193 (3579 and 5965), Corollary 8.8(ii) with the symmetric Hadamard matrix of order
  596 gives Williamson-type matrices of order 1193 (above). SageMath 10.8 (passagemath 10.8.13) returns Unknown for
  nine of them and raises UnknownError in its Miyamoto test for the other seven.
- 41404, 48932, 63988 and 71516 (not stored; q = 941, T-matrices of order 11, 13, 17 and 19) are **New, beyond the
  tables**: n = 941t lies above Table 1.53, and the symmetric Williamson-type matrices of order 235 need the conference
  matrix of order 118, Concurrent with J. Y. Z. Ho's ([11292](11292)), so the orders are not obtainable classically. The
  seven further orders of the table in [Reach](#reach) are certified, not densely checked.
- SageMath's development version (17 October 2026) constructs none of the 19 orders of the first table and of the table
  of orders recorded only through Corollary 8.8(i). It returns Unknown for 14 and raises UnknownError in its Miyamoto
  test for 12260, 33380, 35780, 38660 and 39348 (no symmetric conference matrix of order n + 1 is known to it; with that
  Unknown treated as False it returns Unknown).

The sources searched are listed in [NOVELTY.md](../NOVELTY.md).

## Reproduction and verification

```
python3 search/miyamoto_blocks.py selftest                    # q = 13 ... 61 with T-matrices of order 3, 5, 7
python3 search/miyamoto_blocks.py certify 1093 3 paley2:1091   # writes hadamard/13116/certificate.json
python3 search/miyamoto_blocks.py certify 809 3 golay2:dk202   # periodic Golay pair: hadamard/9708
python3 search/miyamoto_blocks.py certify 613 5 wmult:153:145,35   # Williamson matrices of order 153: hadamard/12260
python3 search/miyamoto_blocks.py certify 269 7 nw:near-williamson/NW67.txt   # near-Williamson: hadamard/7532
python3 search/miyamoto_blocks.py certify 773 5 "wt:cor29h(193, kron(paley2(5), h2(), h2(), h2()))"   # hadamard/15460
python3 search/symwt.py "cor29r(773, cor29h(193, kron(paley2(5), h2(), h2(), h2())))"   # WT(773), repeated block
python3 verify/verify.py 13116                                # rebuilds the matrix and checks H Hᵀ = 4n I
```

The certificates store q, g, the NG pair ("ng"), the periodic Golay pair ("pg"), the near-Williamson matrices ("nw") or
the symmetric Williamson-type matrices ("wt", their upper triangles in base64) and the T-matrices. `verify/verify.py`
(construction `miyamoto-blocks`) checks the pair and the T-matrices, rebuilds Paley's matrix, the blocks, the reversal
and the Goethals–Seidel array independently of the search program (Miyamoto's array with the same `miyamoto` function,
which the search program imports from `verify/verify.py`), and checks H Hᵀ = 4n I densely. An independent
re-implementation reproduces the SHA-256 of 13116, 17516 and 29292 and the dense check of 9708.
