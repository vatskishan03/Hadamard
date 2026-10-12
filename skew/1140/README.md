# Skew-Hadamard matrix of order 1140

1140 = 4 · 285, with 285 = 3 · 5 · 19.

## Construction

A skew-Hadamard matrix is a Hadamard matrix H with H + Hᵀ = 2I.

**Difference family.** Four ±1 sequences x₁, …, x₄ on ℤ₂₈₅ satisfy two conditions:

- each is constant on the orbits of the multiplier subgroup
  H = ⟨61⟩ ≤ ℤ₂₈₅*, of order 9: {1, 16, 61, 106, 121, 196, 226, 256, 271};
- their periodic autocorrelations sum to zero at every nonzero shift.

`certificate.json` lists, for each sequence, one representative r of each orbit rH on which the sequence is −1.

**Skew first block.** The first sequence is of skew type: x₁(0) = 1 and x₁(−t) = −x₁(t) for t ≠ 0. This is possible because −1 ∉ H, so negation pairs the orbits.

**Row sums.** The row sums are 1, 17, 15, 25, and 1² + 17² + 15² + 25² = 1140.

**The matrix.** Let A, B, C, D be the circulant matrices with these first rows (`sequences.txt`) and R the back-diagonal identity. The Goethals–Seidel array

```
    [  A     BR     CR     DR  ]
H = [ -BR    A     DᵀR   -CᵀR ]
    [ -CR   -DᵀR    A     BᵀR ]
    [ -DR    CᵀR  -BᵀR    A   ]
```

is then a Hadamard matrix of order 1140. Because A + Aᵀ = 2I, it satisfies H + Hᵀ = 2I (Wallis–Whiteman).

## Search

1. **Structure: a Paley product.** 285 = 15 · 19 and H = {h : h ≡ 1 (mod 15), h a square mod 19} = ⟨61⟩, of order 9. Every block is u ⊗ 1 + w ⊗ J + z ⊗ χ over ℤ₁₅ × ℤ₁₉; see `attack/NOTES.md` §15.
2. **Reduction.** The counting identity reads 20 N_w − 4 N₋ = 60, so N_w = 3 and N₋ = 0: one w-position in each of B, C and D, translated to 0, with u = 0 there.
3. **Exact two-stage search** (`attack/paley_config.py`, general mode). Stage 1: 2,046,720 quadruples of compressions to ℤ₁₅ (6 s). Stage 2 on the first: 28 z-solutions; the first two were lifted to ℤ₂₈₅ and verified.

## Reproduction

| | |
|---|---|
| Code | commit `520f0b6` |
| Command | `python3 attack/paley_config.py 285 --maxu 6` |
| Randomness and caps | none in the search (exact enumeration); which solutions join4x prints first depends on the OpenMP schedule. |
| Result of that run | 2,046,720 y-quadruples; the first extends to 28 z-solutions; 2 verified skew-Hadamard matrices after 20 s on 4 shared cores |

The tools write hit logs, temporary files and compiled helpers to `$HADAMARD_WORK` (default `<tmp>/hadamard-work`). The certificate does not depend on the search: `verify/verify.py` rebuilds and checks the matrix from `certificate.json` alone.

## Status

To our knowledge this is the first skew-Hadamard matrix of order 1140 on record, but the order was **within reach of a
classical construction**:

- Every list we found has 1140 as unknown. Cati and Pasechnik (arXiv:2411.18897v2, Table 5) list 285(3): the smallest
  known skew order with odd part 285 was 2280 (Wallis–Whiteman 4(q + 1) with q = 569). 1140 is one of the 41 skew orders
  ≤ 1208 that their introduction counts as unknown, and it is in SageMath's `unknown_skew_hadamard` list. The
  Đoković–Seberry census (May 2025) lists 285 as undecided; Đoković's 2007–2008 lists and Koukouvinos–Stylianou (2008,
  285(3)) agree.
- **However**, Seberry's theorem (Ars Combin. 6 (1978), Thm 7; it is how the census obtains 4 · 189) gives a skew
  matrix of order mn(n − 1) from an orthogonal design OD(mn; 1, m, mn − m − 1) and amicable designs of order n. With
  m = 3, n = 20 (q = 19) it needs an OD(60; 1, 3, 56). Our audit agent found one in 0.1 s by a Goethals–Seidel search
  over ℤ₁₅ (A = a then c·(+−++−+−+−+−−+−); B = b then c·(−−−+−+−−−+−−++); C = b then c·(−−+−−+++++++++);
  D = b then c·(−−−++−−+−−−+++)), substituted a → I + P, b → J, c → (I + P)R₁ with P the Paley core of order 19, and
  checked H Hᵀ = 1140 I and H + Hᵀ = 2I; we re-ran that check. We found no record of this design or matrix.

So the honest claim is: first recorded skew-Hadamard matrix of order 1140, by a method (the ℤ₁₅ × ℤ₁₉ Paley product)
that is independent of the Seberry route; the order itself was obtainable classically.

**Checklist** (audit of 2026-10-12 by a research agent): no 2023–2026 paper or deposit we checked (arXiv; Zenodo:
Schumer, Piza, whose 1140 entry 2(569 + 1) is ordinary, Kline, Bedani; the JD-Jones-ASES repositories; Alpöge et al.)
contains it, and no infinite family reaches it: Paley I (1139 = 17 · 67); 2(q + 1) with q = 569 ≡ 1 (mod 8), which the
Szekeres, Whiteman, Spence and Momihara–Xiang families exclude; 2ᵘ(qᵉ + 1) (qᵉ = 284); Spence 4(q² + q + 1);
(p − 1)ᵘ + 1; products. Not checked: the Seberry–Yamada 2020 book and the paywalled journals of 2025–2026.

## Verification

```
python3 verify/verify.py skew/1140
gzip -dc skew/1140/H1140.csv.gz | ./check_matrix
```

SHA-256 of `H1140.csv`:
`14f790c1fe19dee7e9003337db339f1fd5a5413efa11618a2fbc7a26e553dba8`
