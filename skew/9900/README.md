# Skew-Hadamard matrix of order 9900

9900 = 4 · 2475 = 100 · 99, with 2475 = 3² · 5² · 11.

## Construction

A skew-Hadamard matrix is a Hadamard matrix H with H + Hᵀ = 2I. This one comes from Seberry's theorem (Ars Combin.
1978): amicable Hadamard matrices of order n and a skew-type orthogonal design of type (1, m, mn − m − 1) in order mn
give a skew-Hadamard matrix of order mn(n − 1). Here n = 100 and m = 1. Amicable Hadamard matrices of order 100 were
not known (the classical amicable orders ≡ 4 (mod 8) are q + 1 with q ≡ 3 (mod 8) a prime power, and 99 = 3² · 11);
they come from a skew-propus quadruple of length 25, as for [symmetric/5700](../../symmetric/5700) (length 19) and
[skew/2652](../2652) (length 13).

**1. Amicable Hadamard matrices of order 100.** A skew-propus quadruple of length 25 (the first three lines of
`sequences.txt`: A3 symmetric, A1, A0 of skew type, row sums 9, 3, 1, PAF_A3 + 2 PAF_A1 + PAF_A0 = 100 δ) gives, in
the Balonin array, a symmetric orthogonal design xP + yN of type (1, 99): P a symmetric signed permutation, N
symmetric with the complementary support, P² = I, N² = 99 I,
PN + NP = 0. Then S = I + PN is a skew-Hadamard matrix (PN is skew because P and N anticommute) and M = P + N a
symmetric Hadamard matrix, and since PNP = −N,

    S Mᵀ = (I + PN)(P + N) = P + N + PNP + 99 P = 100 P = (P + N)(I + NP) = M Sᵀ,

so (S, M) is an amicable pair.

**2. Amicable cores.** Let r be the first index with P[r, r] = +1 and e = S[r]. Conjugating S and M by diag(e) makes
row r of both all ones (because S[r, j] = N[r, j] = M[r, j] off the diagonal). With r moved first,
S = [[1, 1ᵀ], [−1, I + W]] and M = [[1, 1ᵀ], [1, D]], where W (order 99) is skew with W1 = 0 and
(I + W)(I + W)ᵀ = 100 I − J, D is symmetric with D1 = −1 and D² = 100 I − J, and WD + DW = 0 (from S Mᵀ = M Sᵀ).

**3. A skew-type OD(100; 1, 1, 98).** Four ±1 sequences of length 25 (the last four lines of `sequences.txt`;
"G-matrices"), k₀ and k₁ of skew type and g₂, g₃ symmetric, with row sums 1, 1, 7, 7 and
PAF_k₀ + PAF_k₁ + PAF_g₂ + PAF_g₃ = 100 δ. In the Goethals–Seidel array the circulant blocks aI + c(K₀ − I),
bI + c(K₁ − I), cG₂, cG₃ give aI + bB + cC with B and C skew, BBᵀ = I, CCᵀ = 98 I and BCᵀ + CBᵀ = 0.

**4. The plug-in.**

    H = I ⊗ (I + W) + B ⊗ J + C ⊗ D        (order 100 · 99).

H + Hᵀ = 2I because W, B and C are skew and J, D symmetric, and

    H Hᵀ = I ⊗ [(I + W)(I + W)ᵀ + J² + 98 D²] + (B + Bᵀ) ⊗ J − C ⊗ (WD + DW) + (BCᵀ + CBᵀ) ⊗ JD
         = I ⊗ [100 I − J + 99 J + 98 (100 I − J)] = 9900 I.

## Search

`search/sym_mm1.py propus 25` searches skew-propus quadruples of length 25 exhaustively (row sums (a3, a1) with
a3² + 2a1² = 99: (1, 7), (7, 5), (9, 3); A1 tabulated on its autocorrelations up to shifts, all pairs (A3, A0)
probed); the first solution, in the pattern (9, 3), came after 18 minutes (one core, Python). The G-matrices of order
25 are the first found by `g_matrices(25)` (row sums 1, 1, 7, 7; 22 s). G-matrices of order 25 are known
(Georgiou–Koukouvinos, JCMCC 40 (2002), list 13 of them, as tabulated in Đoković, arXiv:0903.5053, Table 1); only the
amicable pair of order 100 is new.

## Reproduction

| | |
|---|---|
| Code | `search/sym_mm1.py` (`amicable_cores`, `g_matrices`, `skew_od_from_g`, `skew_seberry`) and `verify/verify.py` (`build_skew_amicable_plug`) |
| Command | `python3 search/sym_mm1.py propus 25` (first solution); the certificate was written from that solution and `g_matrices(25)`; `python3 verify/verify.py skew/9900` rebuilds and checks everything |
| Randomness and caps | none |

`verify/verify.py` (construction `skew-amicable-plug`) rebuilds the symmetric OD from the skew-propus quadruple and
checks it, checks that S and M are amicable with S skew, normalises them and checks every identity of the cores
(W skew, W1 = 0, (I + W)(I + W)ᵀ = 100 I − J, D symmetric, D1 = −1, D² = 100 I − J, WD + DW = 0), checks the
G-matrices (types, autocorrelation sum) and the skew-type OD(100; 1, 1, 98), assembles H and checks H Hᵀ = 9900 I,
H + Hᵀ = 2I and the SHA-256 of the canonical CSV.

## Status

To our knowledge this is the first skew-Hadamard matrix of order 9900, and the pair of order 100 the first amicable
Hadamard matrices of that order. (A symmetric Hadamard matrix of order 9900 is known; 2475 is not undecided in the
symmetric list of the census.)

**Novelty checklist** (audit of 2026-10-14 by a research agent; downloaded sources were read as text, never
executed):

1. **The Đoković–Seberry census** (PDF of 2025-05-28, md5 `79eb782c…5be`): "Finally we list the 737 undecided cases
   of odd v < 2500: … 2469, 2473, 2475, 2485, 2489, …". 2475 is in none of its skew series and not among its 41
   sporadic values.
2. **Tables of skew orders 2ᵗq**: Cati–Pasechnik, arXiv:2411.18897, Table 5 stops at n ≤ 999; Seberry–Yamada 1992,
   Table 7.1, stops at q = 999.
3. **Amicable orders**: no classical family gives 100 (p^r + 1: 99 = 3² · 11; 2(q + 1) with 2q + 1 a prime power:
   2 · 49 + 1 = 99; Wallis's 2(q + 1) needs q ≡ 5 (mod 8), and 49 ≡ 1; (p − 1)^u + 1 and products: no). The tables
   of amicable orders 2ᵗq of Seberry (1978), Seberry–Yamada (1991) and Seberry–Yamada (1992, Table 6.1) give t = 3
   for q = 25.
4. **Classical skew families** (arithmetic): Paley needs 9899 = 19 · 521; 2(q + 1) needs q = 4949 = 7² · 101;
   4(q + 1), q ≡ 9 (mod 16), needs q = 2474; Spence's q² + q + 1 skips 2475 (2451, 2551); Seberry's q(q + 1) needs
   q = 99; her 2ʳn(2ʳn − 1) needs amicable orders 50 or 25; Wallis's 1971 theorem has no solution with an integral
   auxiliary order of a skew matrix; (p − 1)^u + 1 and (|t| + 1)(q + 1) fail; products and doubling are impossible
   since 4950 is not a Hadamard order. The other splits of Seberry's theorem need OD(900; 1, 75, 824) (n = 12) or
   OD(3300; 1, 825, 2474) (n = 4); we found no source for either.
5. **SageMath** (develop branch of 2026-10-11, traced, not executed): `skew_hadamard_matrix(9900, existence=True)`
   is `Unknown`.
6. **Đoković's skew papers**, Momihara–Xiang, Cati, Seberry 2013 and the 2026 deposits (as listed for
   [skew/5700](../5700)): no mention of 2475, 9900 or amicable 100.

Not checked: Seberry–Yamada's book (Wiley 2020), Seberry's Orthogonal Designs (Springer 2017), the Handbook of
Combinatorial Designs, Koukouvinos–Stylianou (Discrete Math. 308 (2008); abstract only), MathSciNet and zbMATH. This
is a finding of nothing, not a proof that no earlier construction exists.

## Files

| File | Contents |
|---|---|
| `certificate.json` | the skew-propus quadruple ("blocks", positions of −1), its row sums, and the G-matrices ("g") |
| `sequences.txt` | A3, A1, A0, then k₀, k₁, g₂, g₃, one per line |

The matrix itself is not stored (about 245 MB as CSV); `verify/verify.py` rebuilds it in under a minute.

## Verification

```
python3 verify/verify.py skew/9900
python3 verify/verify.py skew/9900 --write && ./check_matrix < skew/9900/H9900.csv
```

SHA-256 of `H9900.csv`: `05b29efaede46727b5bbb50fb8c244200735cb5b77c32debd19a6d8a3db5913a`
