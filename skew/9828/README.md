# Skew-Hadamard matrix of order 9828

9828 = 4 · 2457 = 13 · 28 · 27, with 2457 = 3³ · 7 · 13.

## Construction

A skew-Hadamard matrix is a Hadamard matrix H with H + Hᵀ = 2I. This one comes from Seberry's theorem (Ars Combin. 6
(1978), Theorem 7): an orthogonal design of type (1, m, mn − m − 1) in order mn and amicable Hadamard matrices of order
n give a skew-Hadamard matrix of order mn(n − 1). Here n = 28 (a classical amicable order, q + 1 with q = 27) and
m = 13. The new ingredient is an orthogonal design of type (1, 13, 350) in order 364.

**1. Amicable cores of order 27.** Paley's skew-Hadamard matrix S = I + C of order 28 over GF(27) (C[∞, x] = 1,
C[x, ∞] = −1, C[x, y] = χ(y − x), χ the quadratic character) and P = (−1) ⊕ (x ↦ −x) satisfy PCP = −C, since
χ(−1) = −1. So xP + yCP is a symmetric orthogonal design of type (1, 27) and (S, P + CP) an amicable pair. GF(27) is
GF(3)[x]/(x³ + 2x + 1), with elements encoded by their base-3 digits ("paley": {"p": 3, "poly": [1, 2, 0]}; the
verifier checks that the polynomial is primitive). Normalised at the fixed point of P with sign +1, the pair gives W
(order 27) skew with W1 = 0 and (I + W)(I + W)ᵀ = 28 I − J, and D symmetric with D1 = −1, D² = 28 I − J and
WD + DW = 0.

**2. A skew-type OD(364; 1, 13, 350).** Four circulant blocks of length t = 91 in the Goethals–Seidel array (codes in
`certificate.json`: 1 = a, ±2 = ±b, ±3 = ±c): A₁ = aI + cW₁ with W₁ skew, and Aᵢ = b uᵢ + c Wᵢ (i = 2, 3, 4). The b
parts lie on the subgroup 7ℤ₉₁ ≅ ℤ₁₃. With S = 7 · {7, 8, 11} (an orbit of ×3 on ℤ₁₃):

- u₂ = −δ₀ + 1_S,
- u₃ = −1_{7·{1, 3, 9}} + 1_{7·{4, 10, 12}} (skew),
- u₄ = −1_S.

Their weights are 4, 6 and 3, and their autocorrelations sum to 13δ. The c parts satisfy
Σᵢ PAF(Wᵢ)(s) = 0 for s ≠ 0 and Σᵢ Σₐ uᵢ(a)[Wᵢ(a + s) + Wᵢ(a − s)] = 0 for every s, with row sums 0, −15, −5, −10
(225 + 25 + 100 = 350). All four blocks are invariant under the multiplier group
H = ⟨79, 29⟩ = {1, 9, 16, 22, 29, 53, 74, 79, 81} of order 9, which contains no element ≡ −1 modulo 7 or 13.

**3. The plug-in.**

    H = I₃₆₄ ⊗ (I + W) + B ⊗ J + C ⊗ D        (order 364 · 27).

H + Hᵀ = 2I, the cross terms vanish as in [skew/7956](../7956), and

    H Hᵀ = I ⊗ [(I + W)(I + W)ᵀ + 13 J² + 350 D²] = I ⊗ [28 I − J + 13 · 27 J + 350 (28 I − J)] = 9828 I.

## Search

`attack/odm_engine.py 91 79,29 13` enumerates the H-invariant b parts of weight 13. These are ternary combinations of
H-orbits, with block 1 skew and autocorrelation sum 13δ, found by a meet in the middle on the PAF vectors at the
pair-orbit representatives. They give 18 templates up to multipliers, block permutations and signs, and all 18 pass
the row-sum screen. For each template, `attack/odmx.c` searches the H-invariant c parts exhaustively: per block, the
lists of c parts with the right row sum, each reduced to its vector of autocorrelations and bc terms at the
pair-orbit representatives, joined by a 4-list meet in the middle. Four templates gave 2 solutions each (8 in
total); the certificate uses the first. The whole run took about 12 seconds.

## Reproduction

| | |
|---|---|
| Code | `attack/odm_engine.py`, `attack/odmx.c` (search); `verify/verify.py` (`build_skew_amicable_plug`, "paley" over GF(pᵏ) and "od" forms) |
| Command | `gcc -O2 -o /tmp/hadamard-work/bin/odmx attack/odmx.c && python3 attack/odm_engine.py 91 79,29 13` (the first `SOL` block is the "od" of the certificate) |
| Randomness and caps | none; the search is exhaustive over its structure |

`verify/verify.py` (construction `skew-amicable-plug`) builds GF(27) from the polynomial and checks that it is
primitive. It then builds Paley's design of order 28 and checks it, checks the amicable pair and every identity of
its cores, and rebuilds the Goethals–Seidel design from the codes, checking that it is a skew-type OD(364; 1, 13, 350)
with 364 = 13 · 28. Finally it assembles H and checks H Hᵀ = 9828 I, H + Hᵀ = 2I and the SHA-256 of the canonical CSV.

## Status

To our knowledge this is the first skew-Hadamard matrix of order 9828. The literature audit is in progress. Until it
is recorded here, read this as "found by us", not as "new". What has been checked so far:

1. **The Đoković–Seberry census** (PDF of 2025-05-28): 2457 is in its list of 737 undecided skew values. It is not
   undecided in the symmetric list. Its skew series do not reach 2457:
   - Szekeres needs q = 4913 = 17³ ≡ 5 (mod 8), but 17³ ≡ 1 (mod 8).
   - Szekeres–Whiteman needs p ≡ 5 (mod 8).
   - Spence needs 1 + q + q² = 2457, but 1 + 49 + 49² = 2451.
   - 2457 is not among its 41 sporadic values.
2. **Tables of skew orders 2ᵗq**: Cati–Pasechnik (arXiv:2411.18897), Table 5, stops at n ≤ 999; Seberry–Yamada
   (1992), Table 7.1, at q = 999.
3. **Classical skew families** (arithmetic):
   - Paley needs 9827 = 31 · 317 to be a prime power.
   - 2(q + 1) needs q = 4913 = 17³ ≡ 1 (mod 8), outside the Szekeres, Szekeres–Whiteman and Spence 2(q + 1) series.
   - 4(q + 1) needs q = 2456.
   - Seberry's q(q + 1) and 2ʳn(2ʳn − 1) have no integer solution.
   - (p − 1)^u + 1 needs 9827 to be a perfect power.
   - Kronecker products and doubling are impossible, since 4914 is not a Hadamard order.
4. **SageMath**: its orthogonal-design branch has only OD(28; 1, 1, 26), so the split n = 28, m = 13 is out of its
   reach.
5. **Orthogonal designs.** The audit of [skew/7956](../7956) found no published tables of designs in orders
   4 · odd above 44. Every use of Seberry's Theorem 7 that it found has m = 1.

## Files

| File | Contents |
|---|---|
| `certificate.json` | GF(27) for the amicable pair ("paley") and the OD(364; 1, 13, 350) ("od": type and the four blocks as codes) |
| `H9828.csv.gz` | the matrix |

## Verification

```
python3 verify/verify.py skew/9828
gzip -dc skew/9828/H9828.csv.gz | ./check_matrix
```

SHA-256 of `H9828.csv`: `34ab342604885a25faa446ba1ad07a0411f924795f3f7891ac84904385554154`
