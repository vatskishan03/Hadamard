# Skew-Hadamard matrix of order 2652

2652 = 4 · 663 = 52 · 51, with 663 = 3 · 13 · 17.

## Construction

A skew-Hadamard matrix is a Hadamard matrix H with H + Hᵀ = 2I. This one comes from Seberry's theorem (Ars Combin.
1978): amicable Hadamard matrices of order n and a skew-type orthogonal design of type (1, m, mn − m − 1) in order mn
give a skew-Hadamard matrix of order mn(n − 1). Here n = 52 and m = 1. Amicable Hadamard matrices of order 52 were
not known (no classical family has this order: 51 = 3 · 17 and 25 ≡ 1 (mod 8); Seberry's tables give t = 3 for
q = 13); they come from a skew-propus quadruple of length 13, as for [symmetric/5700](../../symmetric/5700) (which
uses length 19).

**1. Amicable Hadamard matrices of order 52.** A skew-propus quadruple of length 13 (the first three lines of
`sequences.txt`: A3 symmetric, A1, A0 of skew type, row sums 1, 5, 1, PAF_A3 + 2 PAF_A1 + PAF_A0 = 52 δ) gives, in
the Balonin array, a symmetric orthogonal design xP + yN of type (1, 51): P a symmetric signed permutation, N
symmetric with the complementary support, P² = I, N² = 51 I,
PN + NP = 0. Then S = I + PN is a skew-Hadamard matrix (PN is skew because P and N anticommute) and M = P + N a
symmetric Hadamard matrix, and since PNP = −N,

    S Mᵀ = (I + PN)(P + N) = P + N + PNP + 51 P = 52 P = (P + N)(I + NP) = M Sᵀ,

so (S, M) is an amicable pair.

**2. Amicable cores.** Let r be the first index with P[r, r] = +1 and e = S[r]. Conjugating S and M by diag(e) makes
row r of both all ones (because S[r, j] = N[r, j] = M[r, j] off the diagonal). With r moved first,
S = [[1, 1ᵀ], [−1, I + W]] and M = [[1, 1ᵀ], [1, D]], where W (order 51) is skew with W1 = 0 and
(I + W)(I + W)ᵀ = 52 I − J, D is symmetric with D1 = −1 and D² = 52 I − J, and WD + DW = 0 (from S Mᵀ = M Sᵀ).

**3. A skew-type OD(52; 1, 1, 50).** Four ±1 sequences of length 13 (the last four lines of `sequences.txt`;
"G-matrices"), k₀ and k₁ of skew type and g₂, g₃ symmetric, with row sums 1, 1, 5, 5 and
PAF_k₀ + PAF_k₁ + PAF_g₂ + PAF_g₃ = 52 δ. In the Goethals–Seidel array the circulant blocks aI + c(K₀ − I),
bI + c(K₁ − I), cG₂, cG₃ give aI + bB + cC with B and C skew, BBᵀ = I, CCᵀ = 50 I and BCᵀ + CBᵀ = 0.

**4. The plug-in.**

    H = I ⊗ (I + W) + B ⊗ J + C ⊗ D        (order 52 · 51).

H + Hᵀ = 2I because W, B and C are skew and J, D symmetric, and

    H Hᵀ = I ⊗ [(I + W)(I + W)ᵀ + J² + 50 D²] + (B + Bᵀ) ⊗ J − C ⊗ (WD + DW) + (BCᵀ + CBᵀ) ⊗ JD
         = I ⊗ [52 I − J + 51 J + 50 (52 I − J)] = 2652 I.

## Search

`search/sym_mm1.py propus 13` finds skew-propus quadruples of length 13 at once (they exist for every odd length
n ≤ 15 except n = 9, [search/NOTES.md](../../search/NOTES.md) §31); the certificate uses the first one
`skew_propus(13)` returns. The G-matrices of order 13 are the first found by `g_matrices(13)` (all pairs of skew-type
sequences hashed on the sum of their autocorrelations, all symmetric pairs with row sums 5 and 5 probed). G-matrices
of order 13 are known (Georgiou–Koukouvinos, JCMCC 40 (2002), list 8 of them, as tabulated in Đoković,
arXiv:0903.5053, Table 1); only the amicable pair of order 52 is new.

**Reproduction.**

| | |
|---|---|
| Code | `search/sym_mm1.py` (`amicable_cores`, `g_matrices`, `skew_od_from_g`, `skew_seberry`) and `verify/verify.py` (`build_skew_amicable_plug`) |
| Command | `python3 search/sym_mm1.py propus 13`; the certificate was written from `skew_propus(13)[0]` and `g_matrices(13)`; `python3 verify/verify.py skew/2652` rebuilds and checks everything |
| Randomness and caps | none |

## Status

**First explicit construction.** The Đoković–Seberry census lists 663 as undecided for skew matrices, but two older
tables imply that a skew-Hadamard matrix of order 4 · 663 = 2652 is known, and we found no construction behind either
entry. We know of no earlier skew-Hadamard matrix of order 2652, so this is the weakest of our skew novelty claims.
The pair of order 52 is, to our knowledge, the first amicable pair of that order. A symmetric Hadamard matrix of order
2652 is known: 663 is not undecided in the symmetric list of the census.

- The census lists 663 among its 737 undecided skew values ("… 659, 661, 663, 667, 669, …").
- Seberry–Yamada, "Hadamard matrices, sequences, and block designs" (Contemporary Design Theory, Wiley 1992),
  Table 7.1 (p. 500), leaves the t of q = 663 blank, and "a blank represents 2" (p. 498). The same table also leaves
  783 and 879 blank, which the census lists as undecided too.
- Cati–Pasechnik Table 5 ("Pairs n(m) where m is minimal, s.t. a skew Hadamard matrix of order 2^m n is known,
  n ≤ 999. The case m = 2 is not shown") has no entry for 663 (the row reads "655(4) 657(4) 659(…) 661(9) 667(6)
  669(3)"). It also omits 451, 783, 879 and 921, all undecided in the census. Its footnote lists the entries of
  Seberry–Yamada's book (2020, Table 9.2) with m > 2 that it dropped; 663 is not among them, which suggests (not
  checked) that the book has 663 at m = 2 as well.
- We found no construction behind these entries. Every classical family fails for 2652 by arithmetic:
  - Paley: 2651 = 11 · 241.
  - 2(q + 1): q = 1325 = 5² · 53.
  - 4(q + 1), q ≡ 9 (mod 16): q = 662.
  - Spence's q² + q + 1 skips 663.
  - Seberry's q(q + 1): q = 51. Her 2ʳn(2ʳn − 1): amicable orders 26 or 13.
  - (p − 1)ᵘ + 1, (|t| + 1)(q + 1), the small orders 4m, Wallis's 1971 theorem, and Kronecker products or doubling,
    since 1326 is not a Hadamard order.
  - Seberry's theorem itself needs amicable matrices of order 52, which her tables (1978, 1991, 1992) do not have
    (t = 3 for q = 13), or an OD(884; 1, 221, 662) for n = 4, for which we found no source.
- SageMath's source (development version of 2026-10-11) gives `Unknown`.
- The census is the most recent of the three sources and Seberry co-wrote both it and the 1992 table, so we read the
  two table entries as errors, but we could not confirm that.
- The skew tables of Seberry, Ars Combin. 6 (1978) and of Seberry–Yamada, Utilitas Math. 40 (1991) stop at q = 499.
  No source we searched mentions a skew-Hadamard matrix of order 2652 or amicable Hadamard matrices of order 52.
- The Handbook of Combinatorial Designs and Koukouvinos–Stylianou (Discrete Math. 308 (2008); abstract only) were not
  checked.

The sources searched are listed in [NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | the skew-propus quadruple ("blocks", positions of −1), its row sums, and the G-matrices ("g") |
| `sequences.txt` | A3, A1, A0, then k₀, k₁, g₂, g₃, one per line |
| `H2652.csv.gz` | the matrix |

## Verification

```
python3 verify/verify.py skew/2652
gzip -dc skew/2652/H2652.csv.gz | ./check_matrix
```

`verify/verify.py` (construction `skew-amicable-plug`) rebuilds the symmetric OD from the skew-propus quadruple and
checks it, checks that S and M are amicable with S skew, normalises them and checks every identity of the cores
(W skew, W1 = 0, (I + W)(I + W)ᵀ = 52 I − J, D symmetric, D1 = −1, D² = 52 I − J, WD + DW = 0), checks the
G-matrices (types, autocorrelation sum) and the skew-type OD(52; 1, 1, 50), assembles H and checks H Hᵀ = 2652 I,
H + Hᵀ = 2I and the SHA-256 of the canonical CSV.

SHA-256 of `H2652.csv`: `b88c557c9f6689f46a5af203b86ea4a1b2ad0f8a8c36e7db82e2417df55d8250`
