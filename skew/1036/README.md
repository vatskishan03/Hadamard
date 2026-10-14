# Skew-Hadamard matrix of order 1036

1036 = 4 · 259, with 259 = 7 · 37.

## Construction

A skew-Hadamard matrix is a Hadamard matrix H with H + Hᵀ = 2I.

**Difference family.** Four ±1 sequences x₁, …, x₄ on ℤ₂₅₉ satisfy two conditions:

- each is constant on the orbits of the multiplier subgroup
  H = ⟨9⟩ ≤ ℤ₂₅₉*, of order 9: {1, 9, 16, 81, 86, 144, 211, 232, 256};
- their periodic autocorrelations sum to zero at every nonzero shift.

`certificate.json` lists, for each sequence, one representative r of each orbit rH on which the sequence is −1.

**Skew first block.** The first sequence is of skew type: x₁(0) = 1 and x₁(−t) = −x₁(t) for t ≠ 0. This is possible because −1 ∉ H, so negation pairs the orbits.

**Row sums.** The row sums are 1, 7, 31, 5, and 1² + 7² + 31² + 5² = 1036.

**The matrix.** Let A, B, C, D be the circulant matrices with these first rows (`sequences.txt`) and R the back-diagonal identity. The Goethals–Seidel array

```
    [  A     BR     CR     DR  ]
H = [ -BR    A     DᵀR   -CᵀR ]
    [ -CR   -DᵀR    A     BᵀR ]
    [ -DR    CᵀR  -BᵀR    A   ]
```

is then a Hadamard matrix of order 1036. Because A + Aᵀ = 2I, it satisfies H + Hᵀ = 2I (Wallis–Whiteman).

## Search

1. `search/skew_rank.py` ranks skew Goethals–Seidel multiplier structures with a Gaussian model, counting only the solutions that survive when each list is cut to 4·10⁶ sequences. For H = ⟨9⟩ it predicts about 10⁹ solutions in all and 10⁵·⁸ in the cut lists:
   - 31 orbits and 15 autocorrelation classes;
   - 32768 skew first blocks;
   - pools of 3·10⁷ to 4·10⁷ sequences for the other three blocks.
2. `search/skew_ktree.py` (four-list k-tree matcher `search/match4.c`) matched uniform samples of 3.9·10⁶ sequences from each pool, row sums (7, 31, 5).
3. It found 2 solutions for the most likely partial-sum value.

## Reproduction

| | |
|---|---|
| Code | commit `f6f42b2` of this branch (the options used behave the same at later commits) |
| Command | `python3 search/skew_ktree.py 259 9 --maxtau 300 --sols 2 --maxpool 4e6 --seed 1` |
| Randomness and caps | uniform samples of at most 4·10⁶ distinct sequences per pool (seed 1); clamp S = 6; top 6; at most 300 values of τ |
| Result of that run | pattern (5, 7, 31): 2 solutions at τ #0 |

The tools write their hit logs, temporary list files and compiled helpers to `$HADAMARD_WORK` (default
`<tmp>/hadamard-work`), never into the repository. `search/skew_certify.py` (or `skew_group.py certify` for the
non-cyclic orders) turns a logged hit into the files in this directory. The certificate does not depend on
the search: `verify/verify.py` rebuilds and checks the matrix from `certificate.json` alone.

## Status

To our knowledge this is the first skew-Hadamard matrix of order 1036.

- Cati and Pasechnik, *A database of constructions of Hadamard matrices*, arXiv:2411.18897v2 (Aug. 2025), Table 5, lists 259(4). So the smallest known skew order of the form 2ᵐ · 259 was 4144, and 1036 was unknown.
- A literature search in October 2026 found no later construction.

The doubling construction (Seberry) then also gives skew-Hadamard matrices of orders 2072.

**Novelty checklist** (audit of 2026-10-12, by independent agents that also re-verified the matrix):

1. Table 5 of Cati–Pasechnik lists 259(4): the smallest known skew order with odd part 259 was 2⁴ · 259 = 4144. This order is inside the range (≤ 1208) that they checked against implemented constructions.
2. Đoković's papers (1992–2024, alone and with Kotsireas, Golubitsky or Balonin), Koukouvinos–Stylianou
   (Discrete Math. 2008, Tables 1 and 3) and the Đoković–Seberry census of skew orders 4v for odd v < 2500
   (May 2025) do not construct it. The census lists 259 among its 737 undecided values.
3. The 2026 public deposits (Schumer; Alpöge et al.; JD-Jones-ASES; Schneider; Karoui; Piza) do not contain
   it. SageMath's `unknown_skew_hadamard` list (develop branch, fetched 2026-10-12) contains 1036.
4. No infinite family reaches it: Paley I, Szekeres 2(q + 1), Wallis–Whiteman, Spence, Momihara–Xiang,
   the product theorem n(pᵏ + 1) and amicable products were checked.

This is a finding of nothing, not a proof that no earlier construction exists.

## Verification

```
python3 verify/verify.py skew/1036
gzip -dc skew/1036/H1036.csv.gz | ./check_matrix
```

SHA-256 of `H1036.csv`:
`ab03b54030c8a75204713786ec15f3b4a05a81af66c6b01d1dcd8ba515c79169`
