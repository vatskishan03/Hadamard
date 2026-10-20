# Symmetric Hadamard matrices of order 4q from a periodic Golay pair of length (q − 1)/2

**Theorem.** Let q ≡ 1 (mod 4) be a prime power. If a Hadamard matrix K = [[P, Q], [Q, P]] of order q − 1 exists,
then a symmetric Hadamard matrix of order 4q exists. Such a K exists whenever there is a periodic Golay pair of length
(q − 1)/2 over an abelian group.

For q ≡ 1 (mod 8) the theorem adds nothing to [PROPUS_WEIGHING.md](PROPUS_WEIGHING.md), which needs only some Hadamard
matrix of order (q − 1)/2. For q ≡ 5 (mod 8) it is the first construction we know of symmetric Hadamard matrices of
order 4q. Of the 57 prime powers q ≡ 5 (mod 8) that the Đoković–Seberry census (2025-05-28) leaves undecided, it
settles 101 and 149 ([404](404), [596](596)). For eleven more it would need an input of order 388, 772, …, 2340 (see
[Reach](#reach)), and for the other 44 no input exists.

The matrices have an involutory automorphism with 2q − 2 fixed points
([search/NOTES.md](../search/NOTES.md) §63, §65). The construction was found in round 8 of the exploration of 18 October
2026 (its own proof and builder), re-derived independently, and verified by four adversarial checks. Its builder in
[verify/verify.py](../verify/verify.py) (construction `symmetric-involution`) is a third implementation.

## Proof

Let χ be the quadratic character of GF(q), m = q − 1 and n = m/2. Index the projective line by GF(q) ∪ {∞}.

**1. The frame.** For a symmetric (0, ±1) matrix V of order q + 1 with V² = qI, a symmetric signed permutation matrix
P with P² = I whose support is the complement of the support of V, a symmetric ±1 matrix F of order 2m and a ±1
matrix B of size 2m × (q + 1) with

    Bᵀ B = 2m I,   F B = −2 B P,   F² + 2 B Bᵀ = 4q I,

the matrix

```
    [ F    B      B     ]
H = [ Bᵀ   P + V  P − V ]
    [ Bᵀ   P − V  P + V ]
```

is symmetric with entries ±1, and H² = 4qI. The diagonal block of the last two block rows is
BᵀB + (P ± V)² + (P ∓ V)² = 2mI + 2P² + 2V² = (2m + 2 + 2q) I = 4qI. The off-diagonal block is
BᵀB + 2P² − 2V² = (2m + 2 − 2q) I = 0, since the cross terms PV and VP cancel. The first block row gives
FB + 2BP = 0 and F² + 2BBᵀ = 4qI. In the coordinates (fixed points, first points of the pairs, second points) this is
the shape of every symmetric H(4q) with an involution of 2m fixed points and q + 1 two-cycles. V carries all the field
structure. F, B and P form the Q-part.

**2. The Q-part from K.** On the m points x ≠ 0 let A[x, y] = χ(x + y) and Z[x, y] = [x + y = 0]. Then A = XZ, where X is
Paley's core restricted to GF(q)* and Z is the permutation x ↦ −x, which commutes with X. Paley's conference matrix
C₀ is [[0, 1, 1ᵀ], [1, 0, χᵀ], [1, χ, X]] in the order (∞, 0, GF(q)*), and C₀² = qI gives

    A² = X² = qI − 𝟙𝟙ᵀ − χχᵀ,   A𝟙 = X𝟙 = −χ,   Aχ = Xχ = −𝟙,

the last by the Jacobi sum Σ_y χ(y) χ(y − x) = −1. Relabel the rows of K by the points of GF(q)*: row r < n is the
point r + 1, and row n + r is −(r + 1). Then x ↦ −x acts on the rows as the block swap S, and Z Kp = Kp S, because
S K = K S. Write J₂ = [[1, 1], [1, 1]] and D₂ = [[1, −1], [−1, 1]]. Put

    F = A ⊗ J₂ + Z ⊗ D₂ = [[A + Z, A − Z], [A − Z, A + Z]]   (in the order (first copy, second copy)),
    B = [ 𝟙 ⊗ (1, 1)ᵀ | χ ⊗ (1, 1)ᵀ | k_j ⊗ (1, −1)ᵀ  (j < m) ],

with k_j the columns of Kp. Since J₂ D₂ = 0, J₂² = 2J₂ and D₂² = 2D₂:

- F is symmetric with entries ±1 (A is zero exactly where Z is one), and F² = 2A² ⊗ J₂ + 2I ⊗ D₂;
- BBᵀ = (𝟙𝟙ᵀ + χχᵀ) ⊗ J₂ + Kp Kpᵀ ⊗ D₂ = (𝟙𝟙ᵀ + χχᵀ) ⊗ J₂ + mI ⊗ D₂, so
  F² + 2BBᵀ = 2qI ⊗ (J₂ + D₂) = 4qI;
- BᵀB = 2mI, because 𝟙 ⊥ χ, the columns of K are orthogonal, and (1, 1) ⊥ (1, −1);
- F(𝟙 ⊗ (1, 1)ᵀ) = 2A𝟙 ⊗ (1, 1)ᵀ = −2 χ ⊗ (1, 1)ᵀ, F(χ ⊗ (1, 1)ᵀ) = −2 𝟙 ⊗ (1, 1)ᵀ, and
  F(k_j ⊗ (1, −1)ᵀ) = 2 Z k_j ⊗ (1, −1)ᵀ = 2 k_{S(j)} ⊗ (1, −1)ᵀ.

So FB = −2BP with P = [[0, 1], [1, 0]] ⊕ (−S). This P is a symmetric signed permutation matrix with P² = I, and it has
no fixed points.

**3. The V-part.** Let ν be a non-square. Then

    V₀[x, y] = χ(xy − ν),   V₀[x, ∞] = V₀[∞, x] = χ(x),   V₀[∞, ∞] = 1

is a symmetric weighing matrix of order q + 1 and weight q. Its zeros are y = ν/x and (0, ∞), a fixed-point-free
involution, since x² = ν has no solution. The entries of V₀² are character sums.
- A diagonal entry is the number of nonzero entries, q.
- For x ≠ x′, both nonzero, the sum is Σ_y χ((xy − ν)(x′y − ν)) + χ(xx′) = −χ(xx′) + χ(xx′) = 0. A quadratic with
  leading coefficient xx′ and nonzero discriminant ν²(x − x′)² has character sum −χ(xx′).
- The sums with row 0 or row ∞ are linear sums Σ_y χ(x′y − ν) = 0, or the quadratic Σ_y χ(xy² − νy) = −χ(x),
  cancelled by V₀[x, ∞]V₀[∞, ∞] = χ(x).

V₀ is Paley's matrix twisted by the anti-automorphism x ↦ ν/x: V₀ = C₀W, where W is the signed lift of that map, with
W C₀ Wᵀ = −C₀ and W² = −I. Relabel V₀ so that its zeros are exactly the support of P. Any bijection of 2-cycles works,
since only V² = qI and the support enter step 1. Step 1 then gives the theorem.

**4. K from a periodic Golay pair.** If (a, b) is a periodic Golay pair over an abelian group Γ of order n, i.e.
PAF_a + PAF_b = 2n δ, then K′ = [[A, B], [−Bᵀ, Aᵀ]] (A, B group-developed) is a Hadamard matrix of order 2n. Every
translation of Γ, acting on both halves, is an automorphism. If n is even, Γ has an element t of order 2, and the
translation by t is a fixed-point-free involution. Ordering the indices by (coset of ⟨t⟩, half, transversal) puts K′ in
the form [[P, Q], [Q, P]].

## Which q

**A necessary condition.** K = [[P, Q], [Q, P]] is Hadamard exactly when X = (P + Q)/2 and Y = (P − Q)/2 are weighing
matrices W(n, n/2) with complementary supports. Indeed KKᵀ = [[2(XXᵀ + YYᵀ), 2(XXᵀ − YYᵀ)], …]. For q ≡ 5 (mod 8),
n = (q − 1)/2 ≡ 2 (mod 4), and a weighing matrix of order n ≡ 2 (mod 4) has weight a sum of two squares. So
(q − 1)/4 must be a sum of two squares. The same holds for every Hadamard matrix of order q − 1 with a fixed-point-free
involutory automorphism, by a Hasse–Minkowski argument (NOTES.md §65). Periodic Golay pairs of length ℓ need ℓ to be a
sum of two squares as well.

**Reach.** Of the 57 census values q ≡ 5 (mod 8) that are prime powers (`python3 search/sym_involution.py targets`),
44 have (q − 1)/4 not a sum of two squares. The other 13 need K of order q − 1:

| q | 4q | K of order | status |
|---:|---:|---:|---|
| 101 | [404](404) | 100 | periodic Golay pair of length 50 over ℤ₂ × ℤ₅² (ours, from classes of index 8 in GF(25)) |
| 149 | [596](596) | 148 | periodic Golay pair of length 74 (Đoković–Kotsireas, arXiv:1310.5773) |
| 389 | 1556 | 388 | needs a K of order 388; periodic Golay pairs of length 194 are undecided (Đoković–Kotsireas' list) |
| 773, 821, 1109, 1493, 1621, 1733, 1973, 2197, 2309, 2341 | | 772 … 2340 | no K known (periodic Golay pairs of lengths 386, 410, 554, 746, 810, 866, 986, 1098, 1154, 1170 are not known; cyclic ones of length 810 do not exist, by Arasu–Xiang) |

Round 9 also built K of orders 36, 52, 164, 244 and 452, and symmetric H(148), H(212) of known orders (q = 37, 53). An
exhaustive multiplier search found no K for the other census orders in the structures it covered (NOTES.md §65).

## Status

Both orders are **New**: v = 101 and v = 149 are undecided in the census, and we found no symmetric Hadamard matrix of
order 404 or 596 and no statement of this construction in the sources we read (audit of 18 October 2026, NOVELTY.md).

- Miyamoto (1991, Theorem 5) claimed a Hadamard matrix of order 4q "of symmetric block form" from a conference matrix
  of order q + 1 and a symmetric Hadamard matrix of order q − 1, which would cover both orders. The census says his
  proof does not give it ([PROPUS_MIYAMOTO.md](PROPUS_MIYAMOTO.md#status)). The input here is different: K need not
  be symmetric, but must have the form [[P, Q], [Q, P]].
- Ordinary Hadamard matrices of orders 404 and 596 are classical (Miyamoto's theorem with any H(q − 1)). The symmetry
  is what is new.
- The periodic Golay pair of length 74 is Đoković and Kotsireas'. For length 50 the round-9 agent cites earlier pairs
  (Kotsireas–Koukouvinos 2008; Đoković 2011), which we did not check; the one used here is from our own search. None
  of ours is cyclic: our search found no cyclic pair of length 50 invariant under a multiplier subgroup of order 4,
  5, 10 or 20 (exhaustive within those subgroups).
- Symmetric matrices of orders 148 (q = 37) and 212 (q = 53) were also built this way; those orders were known.

## Reproduction and verification

```
python3 search/sym_involution.py targets          # the census prime powers q = 5 (mod 8) and their inputs
python3 search/sym_involution.py certify 101 149  # writes symmetric/404 and symmetric/596
python3 verify/verify.py symmetric/404            # rebuilds H and checks it
```

`verify/verify.py` (construction `symmetric-involution`) does the following.
- From the pair: it checks that (a, b) is a periodic Golay pair over the stated group and that t has order 2, and
  builds K and checks that it is of the form [[P, Q], [Q, P]] and Hadamard.
- From the Q-part: it builds F, B and P as in step 2 and checks BᵀB = 2mI, FB = −2BP, F² + 2BBᵀ = 4qI, and that P is
  a fixed-point-free symmetric signed permutation.
- From the V-part: it builds V₀, checks V₀² = qI and its zero pattern, and relabels it.
- On the result: it checks H = Hᵀ, HHᵀ = 4qI and the SHA-256 of the canonical CSV.
