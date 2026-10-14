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
even power. For t = 39, 19 = 2² + 2 · 3 + 3² and the row sums are 4, 6, −10. (Among the odd t < 39, the condition fails for
t = 5, 11, 13, 17, 21, 23, 29, 31, 35 and 37.)

## Search

`search/od13x.c` searches these designs exhaustively among those invariant under a multiplier group H of ℤₜ. Here
H = {1, 16, 22}, the largest admissible group: since W₀ is skew, no element of H may be ≡ −1 modulo 3 or modulo 13.
The 19 pairs {s, −s} fall into 7 H-orbits. The 2⁷ possible W₀ go into a hash table keyed on their autocorrelations
at the orbit representatives, and every assignment of one of the 20 admissible patterns per orbit to (W₂, W₃, W₄)
with symmetric-part sums (−5, 2, 3) is looked up. The PAFs are updated incrementally from tables of orbit-pair
contributions. Up to permuting W₂, W₃, W₄ and negating all three, every admissible row-sum target is equivalent, so
one target suffices. The run took 6 seconds: 4 285 440 assignments, 2688 solutions (1344 triples (W₂, W₃, W₄), each
with W₀ and −W₀). The certificate uses the first one printed.

Simulated annealing over all designs (`search/od13.c`, and `search/od13b.c` with the row sums held fixed) solves
t ≤ 27 at once but did not reach t = 39 in about 20 minutes. The symmetric search did.

## Reproduction

| | |
|---|---|
| Code | `search/od13x.c` (search); `search/sym_mm1.py` (`amicable_cores`, `skew_od_from_codes`, `check_skew_od`, `skew_seberry`); `verify/verify.py` (`build_skew_amicable_plug`, "od" form) |
| Command | `gcc -O2 -o od13x search/od13x.c && ./od13x 39 16` (the first `SOL` block is the "od" of the certificate) |
| Randomness and caps | none; the search is exhaustive over its structure |

`verify/verify.py` (construction `skew-amicable-plug`) rebuilds the amicable pair of order 52 from the skew-propus
quadruple and checks every identity of its cores. It then rebuilds the Goethals–Seidel design from the codes and
checks that it is a skew-type OD(156; 1, 3, 152) with 156 = 3 · 52. Finally it assembles H and checks H Hᵀ = 7956 I,
H + Hᵀ = 2I and the SHA-256 of the canonical CSV.

## Status

To our knowledge this is the first skew-Hadamard matrix of order 7956, and the OD(156; 1, 3, 152) the first
orthogonal design of that type.

**Novelty checklist** (audits of 2026-10-14 by research agents; downloaded sources were read as text, never
executed):

1. **The Đoković–Seberry census** (PDF of 2025-05-28): 1989 is in its list of 737 undecided skew values
   ("… 1985, 1987, 1989, 1993, …") and in its list of undecided symmetric values.
2. **Tables of skew orders 2ᵗq**:
   - Cati–Pasechnik (arXiv:2411.18897), Table 5, stops at n ≤ 999.
   - Seberry–Yamada (1992), Table 7.1, stops at q = 999.
   - Seberry–Yamada (Utilitas Math. 40 (1991)) promises tables for q ≤ 2000, but its printed skew table stops at
     q = 499.
3. **The orthogonal design.** By Geramita–Verner (quoted in Seberry 1978; Handbook of Combinatorial Designs, 2nd
   ed., Thm 2.74), OD(156; 1, 3, 152) exists exactly when OD(156; 3, 152) does. Neither was found.
   - The published tables of 2- and 3-variable designs in orders 4 · odd stop far below 156. They list (1, 3, 8)
     (order 12), (1, 3, 24) (order 28) and (1, 3, 32) (order 36). They list (3, 16) in order 20 as nonexistent, in
     line with the norm condition above. Sources: Georgiou–Koukouvinos–Mitrouli–Seberry, J. Statist. Plann.
     Inference 106 (2002); Koukouvinos–Mitrouli–Seberry–Karabelas, Australas. J. Combin. 13 (1996); Handbook
     Tables 2.98–2.99.
   - The appendices of Seberry's Orthogonal Designs (Springer 2017) stop at order 44 for 4 · odd (table of
     contents only).
   - The general families (good, G-, best and T-matrices; Seberry's (1, 1, 4p − 2)) give other types.
4. **Seberry's Theorem 7 with m ≥ 3**: every use we found has m = 1 (Seberry 1978, Corollaries 8 and 9; the census,
   §4; SageMath, whose OD table holds only (1, 1, 26); Cati–Pasechnik, Theorem 12). Seberry 1978 says of these
   designs that "in general little is known".
5. **Classical skew families** (arithmetic):
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
6. **SageMath** (develop branch of 2026-10-11, its existence logic traced, not executed): `Unknown`.
7. **Searches**: web search, the arXiv API (including all Hadamard and orthogonal-design postings of June–October
   2026), the Zenodo API, GitHub code search, OEIS, and the 2026 deposits already collected (Schumer, Karoui, Piza,
   Kline, schneiderlo/hadamard-2060, JD-Jones-ASES, bbeartheancient/hoa64, munimtechnologies/hadamard-orders). None
   mentions 7956, 4 · 1989 or OD(156; 1, 3, 152).

Not checked:
- the body of Seberry–Yamada's book (Wiley 2020; its "Summary of Skew Hadamard Orders");
- the bodies of Seberry's Orthogonal Designs (2017) and of Geramita–Seberry (1979);
- Đoković, J. Combin. Inform. System Sci. 20 (1995);
- MathSciNet and zbMATH.

The census of 2025, co-written by Seberry, still lists 1989 as undecided. This is a finding of nothing, not a proof
that no earlier construction exists.

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
