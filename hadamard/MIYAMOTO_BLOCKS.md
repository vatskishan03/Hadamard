# Miyamoto's matrix in bordered blocks, times T-matrices

Miyamoto's Hadamard matrix of order 4q, built with a Hadamard matrix of order q − 1 of 2N-type, is a 4 × 4 array of
bordered blocks that commute and are transposed by one permutation. Like Paley's matrix in circulant blocks
([PALEY_BLOCKS.md](PALEY_BLOCKS.md)), it can therefore be multiplied by T-matrices in the Goethals–Seidel array, and
gives Hadamard matrices of order 4qt. We found the route stated nowhere. It gives 15 orders that Table 1.53 of the
Handbook lists as unknown (n = 3279 to 9991), two of them orders that Đoković claimed only through Seberry–Yamada's
Corollary 8.8(i). With a periodic Golay pair in place of the negaperiodic one ([below](#periodic-golay-pairs)) it also
gives 9708 = 4 · 2427, an order that FRONTIER.md listed as unresolved, and 35596 = 4 · 8899.

| Order | n | q | t | negaperiodic Golay pair of length (q − 1)/2 |
|------:|---:|---:|---:|---|
| [13116](13116) | 3279 = 3 · 1093 | 1093 | 3 | Paley, second series (p = 1091) |
| [17516](17516) | 4379 = 29 · 151 | 29 | 151 | Paley, first series (p = 13) |
| [18908](18908) | 4727 = 29 · 163 | 29 | 163 | Paley, first series (p = 13) |
| [20028](20028) | 5007 = 3 · 1669 | 1669 | 3 | Paley, second series (p = 1667) |
| [28092](28092) | 7023 = 3 · 2341 | 2341 | 3 | Paley, second series (p = 2339) |
| [28716](28716) | 7179 = 3 · 2393 | 2393 | 3 | Golay 26 × Williamson matrices of order 23 |
| [29292](29292) | 7323 = 3 · 2441 | 2441 | 3 | Golay 10 × Paley, first series (p = 121 = 11²) |
| [30604](30604) | 7651 = 7 · 1093 | 1093 | 7 | Paley, second series (p = 1091) |
| [33380](33380) | 8345 = 5 · 1669 | 1669 | 5 | Paley, second series (p = 1667) |
| [35780](35780) | 8945 = 5 · 1789 | 1789 | 5 | Paley, second series (p = 1787) |
| [36668](36668) | 9167 = 89 · 103 | 89 | 103 | Paley, first series (p = 43) |
| [38660](38660) | 9665 = 5 · 1933 | 1933 | 5 | Paley, second series (p = 1931) |
| [39036](39036) | 9759 = 3 · 3253 | 3253 | 3 | Paley, second series (p = 3251) |
| [39348](39348) | 9837 = 9 · 1093 | 1093 | 9 | Paley, second series (p = 1091) |
| [39964](39964) | 9991 = 97 · 103 | 97 | 103 | Paley, first series (p = 47) |

## Construction

**1. Miyamoto's matrix.** Let q ≡ 1 (mod 4) be a prime, m = (q − 1)/2 and g a primitive root mod q. List the points of
Paley's symmetric conference matrix C of order q + 1 as ∞, 0, g^(2j), g^(2j+1) (j ∈ ℤ_m). Then
C(g^(2i), g^(2j)) = χ(g^(2(j−i)) − 1) and the other quadratic-residue blocks are circulant in j as well. Miyamoto's
theorem (J. Combin. Theory Ser. A 57 (1991)) turns C and a Hadamard matrix K of order q − 1 into a Hadamard matrix of
order 4q; `verify/verify.py` implements it (`miyamoto`, in SageMath's block layout). Its sixteen blocks of order q are

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
matrices.

**4. Multiplication by T-matrices.** With T-matrices T₁, …, T₄ of order t put A_j = Σ_k T_k ⊗ N_kj. These are ±1
matrices of order tq in a commutative, transpose-closed algebra, R = R_t ⊗ R_B (R_t: x ↦ −x on ℤ_t) satisfies
R A_j R = A_jᵀ, and Σ_j A_j A_jᵀ = Σ_k T_k T_kᵀ ⊗ 4q I = 4qt I by the orthogonality of the block rows of Miyamoto's
matrix. The Goethals–Seidel array

    [[A₁, A₂R, A₃R, A₄R], [−A₂R, A₁, A₄ᵀR, −A₃ᵀR], [−A₃R, −A₄ᵀR, A₁, A₂ᵀR], [−A₄R, A₃ᵀR, −A₂ᵀR, A₁]]

is a Hadamard matrix of order 4qt. This is Seberry and Yamada's Theorem 4 (1990) over the algebra of bordered blocks
instead of a group ring; the Paley-block route is the case of circulant blocks.

**NG pairs used.** All from Balonin–Đoković (arXiv:1508.00640) and checked by the program and by `verify/verify.py`:
- *First Paley series* (their Proposition 6): for an odd prime power p, a negacyclic conference matrix of order p + 1
  with first row (0, c₁, …, c_p) gives the pair (1, c₁, …, c_p), (−1, c₁, …, c_p) of length p + 1. We build the matrix
  on PG(1, p) along a primitive element w of GF(p²): c_k = χ(Δ(1, w^k)/Δ(1, w)) with Δ(x, y) = x y^p − x^p y.
- *Second Paley series* (Proposition 8; Ito's series): for p ≡ 3 (mod 4), (1, c₂, c₄, …, c_(p−1)) and (c₁, c₃, …, c_p)
  of length (p + 1)/2.
- *Williamson matrices* of odd order w through their array (29), read in ℤ[y]/(y^(2w) + 1) ≅ ℤ[x]/(x² + 1) ⊗ ℤ[ℤ_w]
  (y ↦ xz): length 2w. For 7179 the symmetric circulant Williamson matrices of order 23 are found by the program.
- *Turyn's multiplication* by Golay pairs of length 10 and 26 (their Proposition 3, formulas (15) and (16)).

## Periodic Golay pairs

The algebra of step 3 only needs the sixteen V_rc to commute and to be transposed by the same reversal as the U_rc. If
(a, b) is a periodic Golay pair of length m (PAF_a + PAF_b = 2m δ), K = [[A, B], [−Bᵀ, Aᵀ]] with A, B circulant is a
Hadamard matrix of order q − 1, every V_rc is circulant, and the plain reversal (h, j) ↦ (h, −j) (border fixed)
transposes every block, since J₂ and D₂ are symmetric. Steps 3 and 4 then hold verbatim with R_B replaced by this
reversal. Formulas (15) and (16) of Balonin and Đoković multiply a Golay pair of length g and a periodic Golay pair of
length v into a periodic Golay pair of length gv (Georgiou, Stylianou, Drosou and Koukouvinos 2014, Theorem 1).

For q = 809, m = 404 = 2 · 202, and Đoković and Kotsireas (Numer. Algorithms 69 (2015), arXiv:1310.5773) give a periodic
Golay pair of length 202 as unions of orbits of H = {1, 87, 95, 137, 185} on ℤ₂₀₂. No negaperiodic Golay pair of length
404 is known (Balonin–Đoković, list (31)). With T-matrices of order 3 and 11 this gives:

| Order | n | q | t | pair of length (q − 1)/2 | Status |
|------:|---:|---:|---:|---|---|
| [9708](9708) | 2427 = 3 · 809 | 809 | 3 | periodic: Golay 2 × Đoković–Kotsireas PG(202) | First explicit construction |
| [35596](35596) | 8899 = 11 · 809 | 809 | 11 | periodic: Golay 2 × Đoković–Kotsireas PG(202) | First explicit construction |

Cati–Pasechnik's Table 4 omits n = 2427, and Đoković's set Δ (arXiv:0912.5091) contains 8899, both through
Williamson-type matrices of order 809 that Seberry–Yamada's Table A.1 obtains only from Corollary 8.8(i)
([NOTES.md](../NOTES.md) §7); FRONTIER.md listed 2427 as unresolved. The variant was found during the novelty audit of
this route. Over the FRONTIER entries and the 691 Handbook entries, periodic pairs (lengths g · v with v ∈ {34, 50, 58,
72, 74, 82, 122, 202, 226}) add only these two.

## Reach

Over the 671 entries 3001 ≤ n ≤ 9999 of Handbook Table 1.53 that had no valid construction
([PALEY_BLOCKS.md](PALEY_BLOCKS.md#beyond-table-4)), the route needs n = qt with q ≡ 1 (mod 4) a prime, T-matrices of
order t, and an NG pair of length (q − 1)/2 from the families above (lengths g · v with g a Golay number and v − 1 or
2v − 1 a prime power, or v = 2w with Williamson or quasi-Williamson matrices of order w known). It reaches 29 entries.
Thirteen are also reached by the Paley-block route and are in this repository already; 6031 is in it by
Yamada–Kiyasu ([24124](24124)). The other 15 are the orders above.

With negaperiodic pairs it reaches no FRONTIER order. The reopened orders of Seberry–Yamada's Corollary 8.8(i) would
need NG pairs that are open: 2611 = 7 · 373 a pair of length 186 (quasi-Williamson matrices of order 93), 1527 and 2545
one of length 254 (order 127), 2427 one of length 404, 2571 one of length 428; 1883 = 7 · 269 needs quasi-Williamson
matrices of order 67, which Balonin and Đoković list as open (their list (31)). A periodic Golay pair of length 404
exists, and gives 2427 ([above](#periodic-golay-pairs)); periodic Golay pairs of length 186, 254, 428 or 134 cannot
exist, since 2m is not a sum of two squares.

## Status

- Of the 11 orders with T-matrices of order 3, 5, 7 or 9, nine are **New, beyond the tables, but obtainable
  classically**: every ingredient was published before 2026 (Paley, Miyamoto 1991, Ito's and the Paley NG pairs,
  Williamson matrices of order 23, Golay pairs, Turyn's T-matrices), and only their combination is new. The other two,
  7179 and 7323, were claimed by Đoković (arXiv:1008.2043, Proposition 3.2) through Corollary 8.8(i) only
  ([NOTES.md](../NOTES.md) §7) and are labelled **First explicit construction**.
- 9708 and 35596 (periodic Golay pairs) are **First explicit constructions**: tables record them as known, through
  Corollary 8.8(i) only.
- The 4 orders with T-matrices of order 103, 151 or 163 (17516, 18908, 36668, 39964) are **New, beyond the tables**.
  Cooper–Wallis with the same T-matrices and Williamson-type matrices of order 29, 89 or 97 (Table A.1 of 1992) also
  gives them.
- SageMath's development version (October 17, 2026) returns Unknown for all 17.

The sources searched are listed in [NOVELTY.md](../NOVELTY.md).

## Reproduction and verification

```
python3 search/miyamoto_blocks.py selftest                  # q = 13 ... 61 with T-matrices of order 3, 5, 7
python3 search/miyamoto_blocks.py certify 1093 3 paley2:1091  # writes hadamard/13116/certificate.json
python3 search/miyamoto_blocks.py certify 809 3 golay2:dk202   # periodic Golay pair: hadamard/9708
python3 verify/verify.py 13116                              # rebuilds the matrix and checks H Hᵀ = 4n I
```

The certificates store q, g, the NG pair ("ng") or the periodic Golay pair ("pg") and the T-matrices. `verify/verify.py`
(construction `miyamoto-blocks`) checks the pair and the T-matrices, rebuilds Paley's matrix, Miyamoto's matrix, the
blocks, the reversal and the Goethals–Seidel array independently of the search program, and checks H Hᵀ = 4n I densely.
