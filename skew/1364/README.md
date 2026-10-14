# Skew-Hadamard matrix of order 1364

1364 = 4 · 341, with 341 = 11 · 31.

## Construction

A skew-Hadamard matrix is a Hadamard matrix H with H + Hᵀ = 2I.

**Difference family.** Four ±1 sequences x₁, …, x₄ on ℤ₃₄₁ satisfy two conditions:

- each is constant on the orbits of the multiplier subgroup
  H = ⟨14⟩ ≤ ℤ₃₄₁*, of order 15: {1, 4, 14, 16, 49, 56, 64, 67, 102, 174, 196, 214, 224, 256, 268};
- their periodic autocorrelations sum to zero at every nonzero shift.

`certificate.json` lists, for each sequence, one representative r of each orbit rH on which the sequence is −1.

**Skew first block.** The first sequence is of skew type: x₁(0) = 1 and x₁(−t) = −x₁(t) for t ≠ 0. This is possible because −1 ∉ H, so negation pairs the orbits.

**Row sums.** The row sums are 1, 9, 21, 29, and 1² + 9² + 21² + 29² = 1364.

**The matrix.** Let A, B, C, D be the circulant matrices with these first rows (`sequences.txt`) and R the back-diagonal identity. The Goethals–Seidel array

```
    [  A     BR     CR     DR  ]
H = [ -BR    A     DᵀR   -CᵀR ]
    [ -CR   -DᵀR    A     BᵀR ]
    [ -DR    CᵀR  -BᵀR    A   ]
```

is then a Hadamard matrix of order 1364. Because A + Aᵀ = 2I, it satisfies H + Hᵀ = 2I (Wallis–Whiteman).

## Search

1. `search/skew_triage.py` ranks skew Goethals–Seidel multiplier structures with a Gaussian model. For H = ⟨14⟩ it predicts about 5·10⁶ solutions:
   - 25 orbits and 12 autocorrelation classes;
   - 4096 skew first blocks;
   - pools of 6.5·10⁵ to 1.3·10⁶ sequences for the other three blocks, row sums (9, 21, 29).
2. `search/skew_ktree.py` (four-list k-tree matcher `search/match4.c`) enumerates the pools and matches the four lists.
3. It found 64 solutions for the most likely partial-sum value, in about a second.

## Reproduction

| | |
|---|---|
| Code | commit `90ab79d` of this branch (the options used behave the same at later commits) |
| Command | `python3 search/skew_ktree.py 341 14 --maxtau 400 --sols 2` |
| Randomness and caps | complete pools (no sampling, no seed); clamp S = 4; top 6; at most 400 values of τ |
| Result of that run | 64 solutions at τ #0 |

The tools write their hit logs, temporary list files and compiled helpers to `$HADAMARD_WORK` (default
`<tmp>/hadamard-work`), never into the repository. `search/skew_certify.py` (or `skew_group.py certify` for the
non-cyclic orders) turns a logged hit into the files in this directory. The certificate does not depend on
the search: `verify/verify.py` rebuilds and checks the matrix from `certificate.json` alone.

## Status

To our knowledge this is the first skew-Hadamard matrix of order 1364.

- Cati and Pasechnik, *A database of constructions of Hadamard matrices*, arXiv:2411.18897v2 (Aug. 2025), Table 5, lists 341(4). So the smallest known skew order of the form 2ᵐ · 341 was 5456, and 1364 was unknown.
- A literature search in October 2026 found no later construction.

The doubling construction (Seberry) then also gives skew-Hadamard matrices of orders 2728.

**Novelty checklist** (audit of 2026-10-12, by independent agents that also re-verified the matrix):

1. Table 5 of Cati–Pasechnik lists 341(4): the smallest known skew order with odd part 341 was 2⁴ · 341 = 5456. This order is above 1208, the range they checked against implemented constructions, so the entry was not taken on trust (items 2 and 3).
2. Đoković's papers (1992–2024, alone and with Kotsireas, Golubitsky or Balonin), Koukouvinos–Stylianou
   (Discrete Math. 2008, Tables 1 and 3) and the Đoković–Seberry census of skew orders 4v for odd v < 2500
   (May 2025) do not construct it. The census lists 341 among its 737 undecided values.
3. The 2026 public deposits (Schumer; Alpöge et al.; JD-Jones-ASES; Schneider; Karoui; Piza) do not contain
   it. SageMath's `unknown_skew_hadamard` list stops at 1196, below this order.
4. No infinite family reaches it: Paley I, Szekeres 2(q + 1), Wallis–Whiteman, Spence, Momihara–Xiang,
   the product theorem n(pᵏ + 1) and amicable products were checked.

This is a finding of nothing, not a proof that no earlier construction exists.

## Verification

```
python3 verify/verify.py skew/1364
gzip -dc skew/1364/H1364.csv.gz | ./check_matrix
```

SHA-256 of `H1364.csv`:
`6f7bdf523e0e79ac996233df286313b5105331522292a050c0d2c238fd42132b`
