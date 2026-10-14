# Skew-Hadamard matrix of order 5700

5700 = 4 · 1425 = 76 · 75, with 1425 = 3 · 5² · 19.

## Construction

A skew-Hadamard matrix is a Hadamard matrix H with H + Hᵀ = 2I. This one comes from Seberry's theorem (Ars Combin.
1978): amicable Hadamard matrices of order n and a skew-type orthogonal design of type (1, m, mn − m − 1) in order mn
give a skew-Hadamard matrix of order mn(n − 1). Here n = 76 and m = 1. Amicable Hadamard matrices of order 76 were
not known (no classical family has this order, and Seberry's tables of amicable orders 2ᵗq give t = 3 for q = 19, see
the status below); they come from the symmetric orthogonal design of [symmetric/5700](../../symmetric/5700).

**1. Amicable Hadamard matrices of order 76.** The skew-propus quadruple of length 19 of symmetric/5700 (the first
three lines of `sequences.txt`) gives, in the Balonin array, a symmetric orthogonal design xP + yN of type
(1, 75): P a symmetric signed permutation, N symmetric with the complementary support, P² = I, N² = 75 I,
PN + NP = 0. Then S = I + PN is a skew-Hadamard matrix (PN is skew because P and N anticommute) and M = P + N a
symmetric Hadamard matrix, and since PNP = −N,

    S Mᵀ = (I + PN)(P + N) = P + N + PNP + 75 P = 76 P = (P + N)(I + NP) = M Sᵀ,

so (S, M) is an amicable pair.

**2. Amicable cores.** Let r be the first index with P[r, r] = +1 and e = S[r]. Conjugating S and M by diag(e) makes
row r of both all ones (because S[r, j] = N[r, j] = M[r, j] off the diagonal). With r moved first,
S = [[1, 1ᵀ], [−1, I + W]] and M = [[1, 1ᵀ], [1, D]], where W (order 75) is skew with W1 = 0 and
(I + W)(I + W)ᵀ = 76 I − J, D is symmetric with D1 = −1 and D² = 76 I − J, and WD + DW = 0 (from S Mᵀ = M Sᵀ).

**3. A skew-type OD(76; 1, 1, 74).** Four ±1 sequences of length 19 (the last four lines of `sequences.txt`;
"G-matrices"), k₀ and k₁ of skew type and g₂, g₃ symmetric, with row sums 1, 1, 5, 7 and
PAF_k₀ + PAF_k₁ + PAF_g₂ + PAF_g₃ = 76 δ. In the Goethals–Seidel array the circulant blocks aI + c(K₀ − I),
bI + c(K₁ − I), cG₂, cG₃ give aI + bB + cC with B and C skew, BBᵀ = I, CCᵀ = 74 I and BCᵀ + CBᵀ = 0.

**4. The plug-in.**

    H = I ⊗ (I + W) + B ⊗ J + C ⊗ D        (order 76 · 75).

H + Hᵀ = 2I because W, B and C are skew and J, D symmetric, and

    H Hᵀ = I ⊗ [(I + W)(I + W)ᵀ + J² + 74 D²] + (B + Bᵀ) ⊗ J − C ⊗ (WD + DW) + (BCᵀ + CBᵀ) ⊗ JD
         = I ⊗ [76 I − J + 75 J + 74 (76 I − J)] = 5700 I.

## Search

Nothing beyond [symmetric/5700](../../symmetric/5700) was searched for the amicable pair. The G-matrices of order 19
were found by `search/sym_mm1.py` (`g_matrices(19)`: all pairs of skew-type sequences hashed on the sum of their
autocorrelations, all symmetric pairs with row sums 5 and 7 probed; half a second). G-matrices of order 19 are
classical (Spence's series for q = 37; Georgiou–Koukouvinos, JCMCC 40 (2002), list 9 of them, as tabulated in Đoković,
arXiv:0903.5053, Table 1), so only the amicable pair is new.

**Reproduction.**

| | |
|---|---|
| Code | `search/sym_mm1.py` (`amicable_cores`, `g_matrices`, `skew_od_from_g`, `skew_seberry`) and `verify/verify.py` (`build_skew_amicable_plug`) |
| Command | none beyond symmetric/5700; `python3 verify/verify.py skew/5700` rebuilds and checks everything |
| Randomness and caps | none |

## Status

**New** as a skew-Hadamard matrix. The Đoković–Seberry census lists 1425 among its 737 undecided values, and no table
of skew orders 2ᵗq reaches it. The pair of order 76 is, to our knowledge, the first amicable pair of that order.

- The census lists 1425 among the undecided cases ("… 1417, 1423, 1425, 1427, 1429, …"). 1425 is in none of its skew
  series (Paley, Szekeres, Szekeres–Whiteman, Spence) and not among its 41 sporadic values. 737 undecided plus 513
  known is the 1250 odd v < 2500.
- The tables of skew orders 2ᵗq stop below 1425: Cati–Pasechnik Table 5 at n ≤ 999, Seberry–Yamada (1992) Table 7.1 at
  q = 999, and Seberry, Ars Combin. 6 (1978), Table 2 and Seberry–Yamada, Utilitas Math. 40 (1991) at q = 499.
- No classical family gives amicable Hadamard matrices of order 76:
  - pʳ + 1: 75 = 3 · 5².
  - 2(q + 1) with 2q + 1 a prime power: 2 · 37 + 1 = 75.
  - Wallis's 2(q + 1) with q ≡ 5 (mod 8) needs q = s² + 4, and 37 = 1² + 4 · 3² gives (|t| + 1)(q + 1) = 152 instead.
  - (p − 1)ᵘ + 1 and products do not apply.
  - The tables of amicable orders 2ᵗq of Seberry (1978, Table 1), Seberry–Yamada (1991) and Seberry–Yamada (1992,
    Table 6.1) all give t = 3 for q = 19: order 152 was known, 76 was not.
  - Two summary lines read literally would cover 76: Wallis 1971, class III, printed "2q − 1" where "2q + 1" is
    meant, and Seberry–Yamada 1992, class AV, "2(q + 1), q ≡ 5 (mod 8)", which is the |t| = 1 case of the line
    above it. The cited papers state the conditions that exclude it.
- No classical skew family reaches 5700:
  - Paley needs 5699 = 41 · 139 to be a prime power.
  - 2(q + 1) (Szekeres, Whiteman, Momihara–Xiang) needs q = 2849 = 7 · 11 · 37.
  - 4(q + 1), q ≡ 9 (mod 16), needs q = 1424.
  - Spence's q² + q + 1 never equals 1425.
  - Seberry's q(q + 1) needs q = 75. Her 2ʳn(2ʳn − 1) needs amicable orders 38 or 19.
  - Wallis's 1971 theorem has no integral solution with known orders.
  - (p − 1)ᵘ + 1 and (|t| + 1)(q + 1) fail.
  - Products and doubling are impossible, since 5700/2 = 2850 is not a Hadamard order.
  - The other splits of Seberry's theorem need OD(300; 1, 15, 284) (n = 20) or OD(1900; 1, 475, 1424) (n = 4); we
    found no source for either.
- By SageMath's source (development version of 2026-10-11), `skew_hadamard_matrix(5700, existence=True)` is
  `Unknown`. Its
  orthogonal-design branch knows only OD(28; 1, 1, 26), and its amicable test fails for 76.
- No source we searched mentions 1425, 5700 or an amicable pair of order 76.
- The Handbook of Combinatorial Designs, Koukouvinos–Stylianou (Discrete Math. 308 (2008); abstract only), and
  Spence's and Whiteman's original papers were not checked.

The sources searched are listed in [NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | the skew-propus quadruple ("blocks", positions of −1), its row sums, and the G-matrices ("g") |
| `sequences.txt` | A3, A1, A0, then k₀, k₁, g₂, g₃, one per line |
| `H5700.csv.gz` | the matrix |

## Verification

```
python3 verify/verify.py skew/5700
gzip -dc skew/5700/H5700.csv.gz | ./check_matrix
```

`verify/verify.py` (construction `skew-amicable-plug`) rebuilds the symmetric OD from the skew-propus quadruple and
checks it, checks that S and M are amicable with S skew, normalises them and checks every identity of the cores
(W skew, W1 = 0, (I + W)(I + W)ᵀ = 76 I − J, D symmetric, D1 = −1, D² = 76 I − J, WD + DW = 0), checks the
G-matrices (types, autocorrelation sum) and the skew-type OD(76; 1, 1, 74), assembles H and checks H Hᵀ = 5700 I,
H + Hᵀ = 2I and the SHA-256 of the canonical CSV.

SHA-256 of `H5700.csv`: `2bbb412f65d96bd36a0f716aff0606898065afe93b9068ad0abcdf8dfd2f57d5`
