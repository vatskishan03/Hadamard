# Skew-Hadamard matrix of order 7956

7956 = 4 · 1989 = 3 · 52 · 51, with 1989 = 3² · 13 · 17.

## Construction

A skew-Hadamard matrix is a Hadamard matrix H with H + Hᵀ = 2I. This one comes from Seberry's theorem (Ars Combin. 6
(1978), Theorem 7): an orthogonal design of type (1, m, mn − m − 1) in order mn and amicable Hadamard matrices of order
n give a skew-Hadamard matrix of order mn(n − 1). Here n = 52 and m = 3, and both ingredients are new: the amicable
pair of order 52 of [skew/2652](../2652), and an orthogonal design of type (1, 3, 152) in order 156 found for this
matrix.

**1. Amicable cores of order 51** (as in [skew/2652](../2652)). The skew-propus quadruple of length 13 (the three lines
of `sequences.txt`) gives a symmetric orthogonal design xP + yN of type (1, 51), hence the amicable pair S = I + PN
(skew) and M = P + N (symmetric) of order 52. Normalised at the first index r with P[r, r] = +1, they give W (order
51) skew with W1 = 0 and (I + W)(I + W)ᵀ = 52 I − J, and D symmetric with D1 = −1 and D² = 52 I − J, with
WD + DW = 0.

**2. A skew-type OD(156; 1, 3, 152).** Four circulant blocks of length t = 39 in the Goethals–Seidel array:
A₁ = aI + cW₀ and Aᵢ = bI + cWᵢ (i = 2, 3, 4), with W₀ skew (W₀[0] = 0, W₀[−g] = −W₀[g]) and Wᵢ[0] = 0, ±1 elsewhere
(`certificate.json`, "od": codes 1 = a, ±2 = ±b, ±3 = ±c). The array is an orthogonal design aI + bB + cC of type
(1, 3, 152) exactly when

- Σᵢ PAF(Wᵢ)(s) = 0 for s ≠ 0 (the c² terms),
- W₂[s] + W₂[−s] + W₃[s] + W₃[−s] + W₄[s] + W₄[−s] = 0 for every s (the bc terms; the ac terms vanish because W₀ is
  skew, the b² terms because each b part is a single entry).

Then A = I, B and C are skew, BBᵀ = 3 I, CCᵀ = 152 I and BCᵀ + CBᵀ = 0.

**3. The plug-in.**

    H = I₁₅₆ ⊗ (I + W) + B ⊗ J + C ⊗ D        (order 156 · 51).

H + Hᵀ = 2I because W, B and C are skew and J, D symmetric. The cross terms vanish (W J = J W = 0, WD + DW = 0,
DJ = JD = −J, B and C skew with BCᵀ + CBᵀ = 0), so

    H Hᵀ = I ⊗ [(I + W)(I + W)ᵀ + 3 J² + 152 D²] = I ⊗ [52 I − J + 3 · 51 J + 152 (52 I − J)] = 7956 I.

**Structure of the design (why t = 39 works).** In the Goethals–Seidel form the b part must be three single entries in
A₂, A₃, A₄: an entry pair in one block has off-peak autocorrelation that, with total weight 3, nothing can cancel, and
A₁ can only carry skew (paired) b entries. Shifting each block puts the entries at 0. The bc condition then says that
on each pair {s, −s} the symmetric parts Sᵢ(s) = (Wᵢ[s] + Wᵢ[−s])/2 of W₂, W₃, W₄ are (0, 0, 0) or a permutation of
(1, −1, 0). The row sums satisfy r₂ + r₃ + r₄ = 0 and r₂² + r₃² + r₄² = 4t − 4, so with rᵢ = 2ρᵢ,
ρ₂² + ρ₂ρ₃ + ρ₃² = (t − 1)/2: a value of the norm form of ℤ[ω], i.e. every prime ≡ 2 (mod 3) divides (t − 1)/2 to an
even power. For t = 39, 19 = 2² + 2 · 3 + 3² and the row sums are 4, 6, −10. (t = 5, 11, 13, 21, 31, 35 and 37 fail
the condition.)

## Search

`attack/od13x.c` searches these designs exhaustively among those invariant under a multiplier group H of ℤₜ. Here
H = {1, 16, 22}, the largest admissible group: since W₀ is skew, no element of H may be ≡ −1 modulo 3 or modulo 13.
The 19 pairs {s, −s} fall into 7 H-orbits. The 2⁷ possible W₀ go into a hash table keyed on their autocorrelations
at the orbit representatives, and every assignment of one of the 20 admissible patterns per orbit to (W₂, W₃, W₄)
with symmetric-part sums (−5, 2, 3) is looked up. The PAFs are updated incrementally from tables of orbit-pair
contributions. Up to permuting W₂, W₃, W₄ and negating all three, every admissible row-sum target is equivalent, so
one target suffices. The run took 6 seconds: 4 285 440 assignments, 2688 solutions (1344 triples (W₂, W₃, W₄), each
with W₀ and −W₀). The certificate uses the first one printed.

Simulated annealing over all designs (`attack/od13.c`, and `attack/od13b.c` with the row sums held fixed) solves
t ≤ 27 at once but did not reach t = 39 in about 20 minutes. The symmetric search did.

## Reproduction

| | |
|---|---|
| Code | `attack/od13x.c` (search); `attack/sym_mm1.py` (`amicable_cores`, `skew_od_from_codes`, `check_skew_od`, `skew_seberry`); `verify/verify.py` (`build_skew_amicable_plug`, "od" form) |
| Command | `gcc -O2 -o od13x attack/od13x.c && ./od13x 39 16` (the first `SOL` block is the "od" of the certificate) |
| Randomness and caps | none; the search is exhaustive over its structure |

`verify/verify.py` (construction `skew-amicable-plug`) rebuilds the amicable pair of order 52 from the skew-propus
quadruple and checks every identity of its cores. It then rebuilds the Goethals–Seidel design from the codes and
checks that it is a skew-type OD(156; 1, 3, 152) with 156 = 3 · 52. Finally it assembles H and checks H Hᵀ = 7956 I,
H + Hᵀ = 2I and the SHA-256 of the canonical CSV.

## Status

To our knowledge this is the first skew-Hadamard matrix of order 7956. The literature audit is in progress. Until it
is recorded here, read this as "found by us", not as "new". What has been checked so far:

1. **The Đoković–Seberry census** (PDF of 2025-05-28): 1989 is in its list of 737 undecided skew values
   ("… 1985, 1987, 1989, 1993, …") and in its list of undecided symmetric values.
2. **Tables of skew orders 2ᵗq**: Cati–Pasechnik (arXiv:2411.18897), Table 5, stops at n ≤ 999; Seberry–Yamada
   (1992), Table 7.1, stops at q = 999.
3. **Classical skew families** (arithmetic):
   - Paley needs 7955 = 5 · 37 · 43 to be a prime power.
   - 2(q + 1) needs q = 3977 = 41 · 97.
   - 4(q + 1) needs q = 1988.
   - Spence's q² + q + 1 skips 1989 (1981, 2071).
   - Seberry's q(q + 1) and 2ʳn(2ʳn − 1) have no integer solution.
   - (p − 1)^u + 1 needs 7955 to be a perfect power.
   - (|t| + 1)(q + 1) has no factorisation with q ≡ 5 (mod 8) a prime power.
   - Kronecker products and doubling are impossible, since 7956/2 is not a Hadamard order.
   - Wallis's 1971 theorem would need amicable orders 52 and 4 with an auxiliary skew order 39.
   - The only other splits of Seberry's theorem are n = 4 (OD(2652; 1, 663, 1988)) and n = 2 (degenerate).
4. **SageMath** (develop branch of 2026-10-11, its existence logic traced, not executed): `Unknown`. Its
   orthogonal-design branch needs OD(156; 1, 3, 152) and amicable matrices of order 52, and has neither.
5. **Texts already collected** for the audits of [skew/5700](../5700) (the census, Seberry 1978, Seberry–Yamada 1991
   and 1992, Đoković's skew papers, Cati–Pasechnik, the 2026 deposits): no mention of 7956, and 1989 only as a year or
   in the census lists.

## Files

| File | Contents |
|---|---|
| `certificate.json` | the skew-propus quadruple ("blocks", positions of −1) and its row sums; the OD(156; 1, 3, 152) ("od": type and the four blocks as codes) |
| `sequences.txt` | A3, A1, A0 of the skew-propus quadruple, one per line |
| `H7956.csv.gz` | the matrix |

## Verification

```
python3 verify/verify.py skew/7956
gzip -dc skew/7956/H7956.csv.gz | ./check_matrix
```

SHA-256 of `H7956.csv`: `85337d21f28bfaa5066fd85e480e55de293a3f195af0d88b6a2e085dc501eda6`
