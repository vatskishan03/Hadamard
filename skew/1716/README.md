# Skew-Hadamard matrix of order 1716

1716 = 4 · 429 = 13 · 12 · 11, with 429 = 3 · 11 · 13.

## Construction

A skew-Hadamard matrix is a Hadamard matrix H with H + Hᵀ = 2I. This one comes from Seberry's theorem (Ars Combin. 6
(1978), Theorem 7): an orthogonal design of type (1, m, mn − m − 1) in order mn and amicable Hadamard matrices of
order n give a skew-Hadamard matrix of order mn(n − 1). Here n = 12 (a classical amicable order, Paley with q = 11)
and m = 13. The new ingredient is an orthogonal design of type (1, 13, 142) in order 156.

**1. Amicable cores of order 11.** Paley's skew-Hadamard matrix S = I + C of order 12 (q = 11) and
P = (−1) ⊕ (x ↦ −x) satisfy PCP = −C, so xP + yCP is a symmetric orthogonal design of type (1, 11) and (S, P + CP)
an amicable pair. Normalised at the fixed point of P with sign +1 (as in [skew/2652](../2652)), they give W
(order 11) skew with W1 = 0 and (I + W)(I + W)ᵀ = 12 I − J, and D symmetric with D1 = −1, D² = 12 I − J and
WD + DW = 0.

**2. A skew-type OD(156; 1, 13, 142).** Four circulant blocks of length t = 39 in the Goethals–Seidel array (codes in
`certificate.json`: 1 = a, ±2 = ±b, ±3 = ±c): A₁ = aI + cW₁ with W₁ skew, and Aᵢ = b uᵢ + c Wᵢ (i = 2, 3, 4). The b
parts lie on the subgroup 3ℤ₃₉ ≅ ℤ₁₃. Writing j for the position 3j:

- u₂ = −δ₀ − δ₁ − δ₃ + δ₄ + δ₇ + δ₈ − δ₉ + δ₁₀ + δ₁₁ + δ₁₂ (weight 10),
- u₃ = −δ₁ − δ₃ − δ₉ (weight 3, the orbit {1, 3, 9} of ×3),
- u₄ = 0.

Their autocorrelations sum to 13δ. The c parts satisfy Σᵢ PAF(Wᵢ)(s) = 0 for s ≠ 0 and
Σᵢ Σₐ uᵢ(a)[Wᵢ(a + s) + Wᵢ(a − s)] = 0 for every s, with row sums 0, −9, −6, −5 (81 + 36 + 25 = 142). All four blocks
are invariant under the multiplier group H = {1, 16, 22} of ℤ₃₉ (order 3), the largest admissible one, as for
[skew/7956](../7956).

**3. The plug-in.**

    H = I₁₅₆ ⊗ (I + W) + B ⊗ J + C ⊗ D        (order 156 · 11).

H + Hᵀ = 2I, the cross terms vanish as in [skew/7956](../7956), and

    H Hᵀ = I ⊗ [(I + W)(I + W)ᵀ + 13 J² + 142 D²] = I ⊗ [12 I − J + 13 · 11 J + 142 (12 I − J)] = 1716 I.

## Search

`search/odm_engine.py 39 16 13 --bsub 13` (about 3 seconds):

1. It enumerates the H-invariant b parts of weight 13 supported on the subgroup of order 13. They are ternary
   combinations of H-orbits with block 1 skew and autocorrelation sum 13δ, found by a meet in the middle. Up to
   multipliers, block permutations and signs there are 10 templates, and all 10 pass the row-sum screen.
2. For each template, `search/odmx.c` searches the H-invariant c parts exhaustively: per block, the lists of c parts
   with the right row sum, each reduced to its vector of autocorrelations and bc terms at the 7 pair-orbit
   representatives, joined by a 4-list meet in the middle.

It found 18 solutions, at most 2 per template; the certificate uses the first. Without the subgroup restriction the
template enumeration in Python is too slow at t = 39: there are 110 763 candidate b parts per block. A template with
a ×3-invariant CW(13, 9) is impossible here, for a row-sum reason given in [search/NOTES.md](../../search/NOTES.md)
§31.

**Reproduction.**

| | |
|---|---|
| Code | `search/odm_engine.py`, `search/odmx.c` (search); `verify/verify.py` (`build_skew_amicable_plug`, "paley" and "od" forms) |
| Command | `python3 search/odm_engine.py 39 16 13 --bsub 13` (the first `SOL` block is the "od" of the certificate; `search/odmx.c` is compiled into `$HADAMARD_WORK/bin` on first use, see [search/README.md](../../search/README.md)) |
| Randomness and caps | none; the search is exhaustive over its structure (`--maxsol 2` per template) |

## Status

**New** as a skew-Hadamard matrix. Cati–Pasechnik Table 5 gives 429(3): the smallest known skew order with odd part
429 was 8 · 429 = 3432. The OD(156; 1, 13, 142) is, to our knowledge, the first orthogonal design of that type.

- The Đoković–Seberry census lists 429 among its 737 undecided skew values. It is not undecided in the symmetric list.
- Seberry–Yamada (1992, Table 7.1) also gives "429 3". That entry is explained by Momihara–Xiang (arXiv:1801.08776;
  Eur. J. Combin. 76 (2019)). Their 4-block skew Hadamard difference families for q ≡ 9 (mod 16) give order 4(q + 1) =
  3432 = 8 · 429 with q = 857. Their 2-block families (order 2(q + 1)) need q ≡ 5 (mod 8), and 857 ≡ 1 (mod 8).
- SageMath's source (development version of 2026-10-11) gives Unknown. Its orthogonal-design branch tries the split n
  = 12, m = 13 and needs an OD(156; 1, 13, 142), which is not in its table.
- No classical skew family reaches 1716:
  - Paley needs 1715 = 5 · 7³ to be a prime power.
  - 2(q + 1) needs q = 857, a prime ≡ 1 (mod 8), outside the Szekeres, Szekeres–Whiteman and Momihara–Xiang 2-block
    families.
  - 4(q + 1) needs q = 428.
  - Spence's q² + q + 1 skips 429 (421, 463).
  - Seberry's q(q + 1) and 2ʳn(2ʳn − 1) have no integer solution (41 · 42 = 1722).
  - (p − 1)ᵘ + 1 needs 1715 to be a perfect power.
  - Kronecker products and doubling are impossible (858 is not a Hadamard order).
- No source we searched mentions a skew matrix of this order or the orthogonal design. The nearest use of Seberry's
  theorem in the literature is Seberry, "Some remarks on amicable orthogonal designs" (Ann. Discrete Math. 6 (1980)),
  whose amicable analogue of Theorem 7 is used only conditionally ("more knowledge of the existence of amicable
  orthogonal designs would be invaluable").
- Spence (1977), Whiteman (1971) and Ding (2008) were not read in full text.

The sources searched are listed in [NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q = 11 for the amicable pair ("paley") and the OD(156; 1, 13, 142) ("od": type and the four blocks as codes) |
| `H1716.csv.gz` | the matrix |

## Verification

```
python3 verify/verify.py skew/1716
gzip -dc skew/1716/H1716.csv.gz | ./check_matrix
```

`verify/verify.py` (construction `skew-amicable-plug`, "paley": 11) rebuilds Paley's design of order 12 and checks it.
It then checks the amicable pair and every identity of its cores, and rebuilds the Goethals–Seidel design from the
codes, checking that it is a skew-type OD(156; 1, 13, 142) with 156 = 13 · 12. Finally it assembles H and checks
H Hᵀ = 1716 I, H + Hᵀ = 2I and the SHA-256 of the canonical CSV.

SHA-256 of `H1716.csv`: `5be71948cce1d25d197252416e7ea631cbcdbef1928f167c9772077aef0e15bb`
