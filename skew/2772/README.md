# Skew-Hadamard matrix of order 2772

2772 = 4 · 693 = 21 · 12 · 11, with 693 = 3² · 7 · 11.

## Construction

A skew-Hadamard matrix is a Hadamard matrix H with H + Hᵀ = 2I. This one comes from Seberry's theorem (Ars Combin. 6
(1978), Theorem 7): an orthogonal design of type (1, m, mn − m − 1) in order mn and amicable Hadamard matrices of
order n give a skew-Hadamard matrix of order mn(n − 1). Here n = 12 (a classical amicable order, Paley with q = 11)
and m = 21. The new ingredient is an orthogonal design of type (1, 21, 230) in order 252.

**1. Amicable cores of order 11.** As in [skew/1716](../1716): from Paley's skew-Hadamard matrix of order 12 and
P = (−1) ⊕ (x ↦ −x), W (order 11) skew with W1 = 0 and (I + W)(I + W)ᵀ = 12 I − J, and D symmetric with D1 = −1,
D² = 12 I − J and WD + DW = 0.

**2. A skew-type OD(252; 1, 21, 230).** Four circulant blocks of length t = 63 in the Goethals–Seidel array (codes in
`certificate.json`: 1 = a, ±2 = ±b, ±3 = ±c): A₁ = aI + cW₁ with W₁ skew, and Aᵢ = b uᵢ + c Wᵢ (i = 2, 3, 4). The b
parts lie on the subgroup 3ℤ₆₃ ≅ ℤ₂₁. Writing j for the position 3j:

- u₂ = −δ₀ − δ₁ − δ₃ − δ₄ − δ₆ + δ₇ + δ₉ − δ₁₂ + δ₁₄ + δ₁₅ − δ₁₆ + δ₁₈ (weight 12),
- u₃ = −δ₂ + δ₇ − δ₈ − δ₉ − δ₁₁ + δ₁₄ − δ₁₅ − δ₁₈ (weight 8),
- u₄ = −δ₁₄ (weight 1).

Their autocorrelations sum to 21δ. The c parts satisfy Σᵢ PAF(Wᵢ)(s) = 0 for s ≠ 0 and
Σᵢ Σₐ uᵢ(a)[Wᵢ(a + s) + Wᵢ(a − s)] = 0 for every s, with row sums 0, −3, 5, −14 (9 + 25 + 196 = 230). The b row sums
(0, −2, −4, −1) are orthogonal to them, as Eades' condition requires. All four blocks are invariant under the
multiplier group H = {1, 4, 16} of ℤ₆₃.

**3. The plug-in.**

    H = I₂₅₂ ⊗ (I + W) + B ⊗ J + C ⊗ D        (order 252 · 11).

H + Hᵀ = 2I, the cross terms vanish as in [skew/7956](../7956), and

    H Hᵀ = I ⊗ [(I + W)(I + W)ᵀ + 21 J² + 230 D²] = I ⊗ [12 I − J + 21 · 11 J + 230 (12 I − J)] = 2772 I.

## Search

`search/odm_join.py 63 4 21 --bsub 21 --bjoin --k 2 --first` (about 5 minutes) runs in three steps:

1. **Templates.** `odm_engine.py --dry --bjoin` builds the H-invariant b parts of weight 21 on the subgroup of
   order 21 as numpy arrays and joins them with `join4x.c`, with the vector (PAF at the pair-orbit representatives,
   weight) and target (0, …, 0, 21). This was a partial sweep: join4x was asked for 200 000 of its
   1 054 944 quadruples, which gave 1848 templates up to equivalence, of which 1440 pass the row-sum screen.
2. **c parts.** For each template, `odmx` dumps the four lists of H-invariant c parts with the right row sums (each
   reduced to its vector of autocorrelations and bc terms at the pair-orbit representatives), and join4x finds the
   zero sums exactly.
3. **Decoding.** `odmx --decode` turns the list indices back into blocks.

The run found 4 solutions; the certificate uses the first.

**Reproduction.**

| | |
|---|---|
| Code | `search/odm_join.py`, `search/odm_engine.py`, `search/odmx.c`, `search/join4x.c` (search); `verify/verify.py` (`build_skew_amicable_plug`, "paley" and "od" forms) |
| Command | `python3 search/odm_join.py 63 4 21 --bsub 21 --bjoin --k 2 --first` (`search/odmx.c` and `search/join4x.c` are compiled into `$HADAMARD_WORK/bin` on first use, see [search/README.md](../../search/README.md)) |
| Randomness and caps | none; the template sweep was capped at 200 000 b quadruples (`--bk`, default) |

## Status

**New** as a skew-Hadamard matrix. Cati–Pasechnik Table 5 gives 693(4): the smallest known skew order with odd part
693 was 16 · 693. The OD(252; 1, 21, 230) is, to our knowledge, the first orthogonal design of that type.

- The Đoković–Seberry census lists 693 among its 737 undecided skew values and among its undecided symmetric values.
- Seberry–Yamada (1992, Table 7.1) also gives "693 4".
- No classical skew family reaches 2772:
  - Paley needs 2771 = 17 · 163 to be a prime power.
  - 2(q + 1) needs q = 1385 = 5 · 277, not a prime power; no two-block skew family in a group of order 1385 was
    found.
  - 4(q + 1) needs q = 692.
  - Spence's q² + q + 1 skips 693 (651, 703).
  - Seberry's q(q + 1) and 2ʳn(2ʳn − 1) have no integer solution (52 · 53 = 2756, 53 · 54 = 2862).
  - (p − 1)ᵘ + 1 needs 2771 to be a perfect power.
  - (|t| + 1)(q + 1) comes closest, with q = 461 = 19² + 4 · 5² and 6 · 462 = 2772, but needs a skew matrix of
    order 6.
  - Kronecker products and doubling are impossible (1386 is not a Hadamard order).
  - The only other splits of Seberry's theorem are n = 4 (OD(924; 1, 231, 692)) and n = 2 (degenerate).
- SageMath's orthogonal-design branch has only OD(28; 1, 1, 26), so the split n = 12, m = 21 is out of its reach.
- No source we searched mentions a skew matrix of this order or the orthogonal design. The nearest use of Seberry's
  theorem in the literature is Seberry, "Some remarks on amicable orthogonal designs" (Ann. Discrete Math. 6 (1980)),
  whose amicable analogue of Theorem 7 is used only conditionally ("more knowledge of the existence of amicable
  orthogonal designs would be invaluable").
- Spence (1977), Whiteman (1971) and Ding (2008) were not read in full text.

The sources searched are listed in [NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q = 11 for the amicable pair ("paley") and the OD(252; 1, 21, 230) ("od": type and the four blocks as codes) |
| `H2772.csv.gz` | the matrix |

## Verification

```
python3 verify/verify.py skew/2772
gzip -dc skew/2772/H2772.csv.gz | ./check_matrix
```

`verify/verify.py` (construction `skew-amicable-plug`, "paley": 11) rebuilds Paley's design of order 12 and checks it.
It then checks the amicable pair and every identity of its cores, and rebuilds the Goethals–Seidel design from the
codes, checking that it is a skew-type OD(252; 1, 21, 230) with 252 = 21 · 12. Finally it assembles H and checks
H Hᵀ = 2772 I, H + Hᵀ = 2I and the SHA-256 of the canonical CSV.

SHA-256 of `H2772.csv`: `c1bb9b575451b2e6364154f64836d2b2d2dbeb3d1855351d15669d90e1453add`
