# Skew-Hadamard matrix of order 5700

5700 = 4 · 1425 = 76 · 75, with 1425 = 3 · 5² · 19.

## Construction

A skew-Hadamard matrix is a Hadamard matrix H with H + Hᵀ = 2I. This one comes from Seberry's theorem (Ars Combin.
1978): amicable Hadamard matrices of order n and a skew-type orthogonal design of type (1, m, mn − m − 1) in order mn
give a skew-Hadamard matrix of order mn(n − 1). Here n = 76 and m = 1. Amicable Hadamard matrices of order 76 were
not known (no classical family has this order); they come from the symmetric orthogonal design of
[symmetric/5700](../../symmetric/5700).

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

Nothing beyond symmetric/5700 was searched for the amicable pair. The G-matrices of order 19 were found by
`attack/sym_mm1.py` (`g_matrices(19)`: all pairs of skew-type sequences hashed on the sum of their autocorrelations,
all symmetric pairs with row sums 5 and 7 probed; half a second). G-matrices of order 19 are known (Seberry's
G-matrices; their existence for this order is in the literature), so only the amicable pair is new.

## Reproduction

| | |
|---|---|
| Code | `attack/sym_mm1.py` (`amicable_cores`, `g_matrices`, `skew_od_from_g`, `skew_seberry`) and `verify/verify.py` (`build_skew_amicable_plug`) |
| Command | none beyond symmetric/5700; `python3 verify/verify.py skew/5700` rebuilds and checks everything |
| Randomness and caps | none |

`verify/verify.py` (construction `skew-amicable-plug`) rebuilds the symmetric OD from the skew-propus quadruple and
checks it, checks that S and M are amicable with S skew, normalises them and checks every identity of the cores
(W skew, W1 = 0, (I + W)(I + W)ᵀ = 76 I − J, D symmetric, D1 = −1, D² = 76 I − J, WD + DW = 0), checks the
G-matrices (types, autocorrelation sum) and the skew-type OD(76; 1, 1, 74), assembles H and checks H Hᵀ = 5700 I,
H + Hᵀ = 2I and the SHA-256 of the canonical CSV.

## Status

The Đoković–Seberry census of skew Hadamard matrices of order 4v for odd v < 2500 (PDF of 2025-05-28) lists 1425
among its 737 undecided values (checked by the audit of symmetric/5700), and Cati–Pasechnik's Table 5 stops at 999.
The skew-specific literature audit is in progress; until it is recorded here, read this as "found by us", not as
"new".

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

SHA-256 of `H5700.csv`: `2bbb412f65d96bd36a0f716aff0606898065afe93b9068ad0abcdf8dfd2f57d5`
