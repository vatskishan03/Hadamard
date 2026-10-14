# Skew-Hadamard matrix of order 2772

2772 = 4 · 693 = 21 · 12 · 11, with 693 = 3² · 7 · 11.

## Construction

A skew-Hadamard matrix is a Hadamard matrix H with H + Hᵀ = 2I. This one comes from Seberry's theorem (Ars Combin. 6
(1978), Theorem 7): an orthogonal design of type (1, m, mn − m − 1) in order mn and amicable Hadamard matrices of order
n give a skew-Hadamard matrix of order mn(n − 1). Here n = 12 (a classical amicable order, Paley with q = 11) and
m = 21. The new ingredient is an orthogonal design of type (1, 21, 230) in order 252.

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

`attack/odm_join.py 63 4 21 --bsub 21 --bjoin --k 2 --first` (about 5 minutes) runs in three steps:

1. **Templates.** `odm_engine.py --dry --bjoin` builds the H-invariant b parts of weight 21 on the subgroup of
   order 21 as numpy arrays and joins them with `join4x.c`, with the vector (PAF at the pair-orbit representatives,
   weight) and target (0, …, 0, 21). This was a partial sweep: join4x was asked for 200 000 of its
   1 054 944 quadruples, which gave 1848 templates up to equivalence, of which 1440 pass the row-sum screen.
2. **c parts.** For each template, `odmx` dumps the four lists of H-invariant c parts with the right row sums (each
   reduced to its vector of autocorrelations and bc terms at the pair-orbit representatives), and join4x finds the
   zero sums exactly.
3. **Decoding.** `odmx --decode` turns the list indices back into blocks.

The run found 4 solutions; the certificate uses the first.

## Reproduction

| | |
|---|---|
| Code | `attack/odm_join.py`, `attack/odm_engine.py`, `attack/odmx.c`, `attack/join4x.c` (search); `verify/verify.py` (`build_skew_amicable_plug`, "paley" and "od" forms) |
| Command | `gcc -O2 -o /tmp/hadamard-work/bin/odmx attack/odmx.c; cc -O3 -march=native -fopenmp attack/join4x.c -o /tmp/hadamard-work/bin/join4x; python3 attack/odm_join.py 63 4 21 --bsub 21 --bjoin --k 2 --first` |
| Randomness and caps | none; the template sweep was capped at 200 000 b quadruples (`--bk`, default) |

`verify/verify.py` (construction `skew-amicable-plug`, "paley": 11) rebuilds Paley's design of order 12 and checks it.
It then checks the amicable pair and every identity of its cores, and rebuilds the Goethals–Seidel design from the
codes, checking that it is a skew-type OD(252; 1, 21, 230) with 252 = 21 · 12. Finally it assembles H and checks
H Hᵀ = 2772 I, H + Hᵀ = 2I and the SHA-256 of the canonical CSV.

## Status

To our knowledge this is the first skew-Hadamard matrix of order 2772, and the OD(252; 1, 21, 230) the first
orthogonal design of that type.

**Novelty checklist:**

1. **The Đoković–Seberry census** (PDF of 2025-05-28): 693 is in its list of 737 undecided skew values and in its
   list of undecided symmetric values.
2. **Tables of skew orders 2ᵗq.** Both give 16 · 693 as the smallest known order:
   - Cati–Pasechnik (arXiv:2411.18897), Table 5: "693(4)".
   - Seberry–Yamada (1992), Table 7.1: "693 4".
3. **Classical skew families** (arithmetic):
   - Paley needs 2771 = 17 · 163 to be a prime power.
   - 2(q + 1) needs q = 1385 = 5 · 277, not a prime power; no two-block skew family in a group of order 1385 was found.
   - 4(q + 1) needs q = 692.
   - Spence's q² + q + 1 skips 693 (651, 703).
   - Seberry's q(q + 1) and 2ʳn(2ʳn − 1) have no integer solution (52 · 53 = 2756, 53 · 54 = 2862).
   - (p − 1)^u + 1 needs 2771 to be a perfect power.
   - (|t| + 1)(q + 1) comes closest, with q = 461 = 19² + 4 · 5² and 6 · 462 = 2772, but needs a skew matrix of
     order 6.
   - Kronecker products and doubling are impossible (1386 is not a Hadamard order).
   - The only other splits of Seberry's theorem are n = 4 (OD(924; 1, 231, 692)) and n = 2 (degenerate).
4. **SageMath**: its orthogonal-design branch has only OD(28; 1, 1, 26), so the split n = 12, m = 21 is out of its
   reach.

**Searches** (audit of 2026-10-14 by a research agent; downloaded sources were read as text, never executed):
- the arXiv API (skew/Hadamard, Szekeres, supplementary difference sets, orthogonal designs, numeric queries, the latest
  papers of Đoković, Seberry, Kotsireas, Kharaghani, Pasechnik and Momihara);
- the Zenodo API, GitHub code and repository search, OEIS and web search;
- a re-grep of about 12 750 previously downloaded texts.

None mentions a skew matrix of this order or the orthogonal design. The nearest use of Seberry's theorem in the
literature is Seberry, "Some remarks on amicable orthogonal designs" (Ann. Discrete Math. 6 (1980)), whose amicable
analogue of Theorem 7 is used only conditionally ("more knowledge of the existence of amicable orthogonal designs would
be invaluable"). Not checked: the bodies of Seberry–Yamada (2020), Seberry (2017) and Geramita–Seberry (1979);
Spence (1977), Whiteman (1971) and Ding (2008) in full text; MathSciNet and zbMATH. This is a finding of nothing, not a
proof that no earlier construction exists.

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

SHA-256 of `H2772.csv`: `c1bb9b575451b2e6364154f64836d2b2d2dbeb3d1855351d15669d90e1453add`
