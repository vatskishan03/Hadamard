# Skew-Hadamard matrix of order 1612

1612 = 4 · 403, with 403 = 13 · 31.

## Construction

A skew-Hadamard matrix is a Hadamard matrix H with H + Hᵀ = 2I.

**Difference family.** Four ±1 sequences x₁, …, x₄ on ℤ₄₀₃ satisfy two conditions:

- each is constant on the orbits of the multiplier subgroup
  H = ⟨100⟩ ≤ ℤ₄₀₃*, of order 15: {1, 66, 87, 100, 133, 152, 157, 237, 287, 289, 315, 326, 328, 360, 386};
- their periodic autocorrelations sum to zero at every nonzero shift.

`certificate.json` lists, for each sequence, one representative r of each orbit rH on which the sequence is −1.

**Skew first block.** The first sequence is of skew type: x₁(0) = 1 and x₁(−t) = −x₁(t) for t ≠ 0. This is possible because −1 ∉ H, so negation pairs the orbits.

**Row sums.** The row sums are 1, 5, 25, 31, and 1² + 5² + 25² + 31² = 1612.

**The matrix.** Let A, B, C, D be the circulant matrices with these first rows (`sequences.txt`) and R the back-diagonal identity. The Goethals–Seidel array

```
    [  A     BR     CR     DR  ]
H = [ -BR    A     DᵀR   -CᵀR ]
    [ -CR   -DᵀR    A     BᵀR ]
    [ -DR    CᵀR  -BᵀR    A   ]
```

is then a Hadamard matrix of order 1612. Because A + Aᵀ = 2I, it satisfies H + Hᵀ = 2I (Wallis–Whiteman).

## Search

1. `attack/skew_rank.py` ranks skew Goethals–Seidel multiplier structures with a Gaussian model, counting only the solutions that survive when each list is cut to a sample. For H = ⟨100⟩ (order 15) it predicts about 10⁶ solutions in all and 10⁴·⁵ in the cut lists:
   - 31 orbits and 15 autocorrelation classes;
   - 32768 skew first blocks;
   - pools of 4·10⁷ to 6·10⁷ sequences for the other three blocks, row sums (5, 25, 31).
2. `attack/skew_compress.py` confirms that the compressions to ℤ₁₃ and ℤ₃₁ admit solutions.
3. `attack/skew_ktree.py` (four-list k-tree matcher `attack/match4.c`) matched uniform samples of about 1.05·10⁷ sequences from each pool.
4. It found 4 solutions among the 154 most likely partial-sum values. The first pattern tried, (5, 19, 35), gave none in 400.

## Reproduction

| | |
|---|---|
| Code | commit `6e27e84` of this branch (the options used behave the same at later commits) |
| Command | `python3 attack/skew_ktree.py 403 100 --maxpool 1.2e7 --maxtau 400 --top 6 --sols 2 --seed 1` |
| Randomness and caps | uniform samples of at most 1.2·10⁷ distinct sequences per pool (seed 1); clamp S = 6; top 6; at most 400 values of τ per pattern |
| Result of that run | pattern (5, 19, 35): none in 400 values of τ; pattern (5, 25, 31): 4 solutions at τ #153 |

The tools write their hit logs, temporary list files and compiled helpers to `$HADAMARD_WORK` (default
`<tmp>/hadamard-work`), never into the repository. `attack/skew_certify.py` (or `skew_group.py certify` for the
non-cyclic orders) turns a logged hit into the files in this directory. The certificate does not depend on
the search: `verify/verify.py` rebuilds and checks the matrix from `certificate.json` alone.

## Status

To our knowledge this is the first skew-Hadamard matrix of order 1612.

- Cati and Pasechnik, *A database of constructions of Hadamard matrices*, arXiv:2411.18897v2 (Aug. 2025), Table 5, lists 403(5). So the smallest known skew order of the form 2ᵐ · 403 was 12896, and 1612 was unknown.
- A literature search in October 2026 found no later construction.

The doubling construction (Seberry) then also gives skew-Hadamard matrices of orders 3224, 6448.

**Novelty checklist** (audit of 2026-10-12, by independent agents that also re-verified the matrix):

1. Table 5 of Cati–Pasechnik lists 403(5): the smallest known skew order with odd part 403 was 2⁵ · 403 = 12896. This order is above 1208, the range they checked against implemented constructions, so the entry was not taken on trust (items 2 and 3).
2. Đoković's papers (1992–2024, alone and with Kotsireas, Golubitsky or Balonin), Koukouvinos–Stylianou
   (Discrete Math. 2008, Tables 1 and 3) and the Đoković–Seberry census of skew orders 4v for odd v < 2500
   (May 2025) do not construct it. The census lists 403 among its 737 undecided values.
3. The 2026 public deposits (Schumer; Alpöge et al.; JD-Jones-ASES; Schneider; Karoui; Piza) do not contain
   it. SageMath's `unknown_skew_hadamard` list stops at 1196, below this order.
4. No infinite family reaches it: Paley I, Szekeres 2(q + 1), Wallis–Whiteman, Spence, Momihara–Xiang,
   the product theorem n(pᵏ + 1) and amicable products were checked.

This is a finding of nothing, not a proof that no earlier construction exists.

## Verification

```
python3 verify/verify.py skew/1612
gzip -dc skew/1612/H1612.csv.gz | ./check_matrix
```

SHA-256 of `H1612.csv`:
`34b140336882e45a3af6f21c499d2242280292d973a831150c1489974868cb6a`
