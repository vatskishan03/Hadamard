# Skew-Hadamard matrix of order 1148

1148 = 4 · 287, with 287 = 7 · 41.

## Construction

A skew-Hadamard matrix is a Hadamard matrix H with H + Hᵀ = 2I.

**Difference family.** Four ±1 sequences x₁, …, x₄ on ℤ₂₈₇ satisfy two conditions:

- each is constant on the orbits of the multiplier subgroup
  H = ⟨16⟩ ≤ ℤ₂₈₇*, of order 15: {1, 16, 18, 37, 51, 57, 78, 92, 100, 141, 165, 221, 242, 247, 256};
- their periodic autocorrelations sum to zero at every nonzero shift.

`certificate.json` lists, for each sequence, one representative r of each orbit rH on which the sequence is −1.

**Skew first block.** The first sequence is of skew type: x₁(0) = 1 and x₁(−t) = −x₁(t) for t ≠ 0. This is possible because −1 ∉ H, so negation pairs the orbits.

**Row sums.** The row sums are 1, 3, 33, 7, and 1² + 3² + 33² + 7² = 1148.

**The matrix.** Let A, B, C, D be the circulant matrices with these first rows (`sequences.txt`) and R the back-diagonal identity. The Goethals–Seidel array

```
    [  A     BR     CR     DR  ]
H = [ -BR    A     DᵀR   -CᵀR ]
    [ -CR   -DᵀR    A     BᵀR ]
    [ -DR    CᵀR  -BᵀR    A   ]
```

is then a Hadamard matrix of order 1148. Because A + Aᵀ = 2I, it satisfies H + Hᵀ = 2I (Wallis–Whiteman).

## Search

1. `attack/skew_triage.py` ranks skew Goethals–Seidel multiplier structures with a Gaussian model. For H = ⟨16⟩ it predicts a few hundred solutions:
   - 27 orbits and 13 autocorrelation classes;
   - 8192 skew first blocks;
   - pools of 8.8·10⁵ to 1.1·10⁶ sequences for the other three blocks, row sums (3, 33, 7).
2. `attack/skew_ktree.py` (four-list k-tree matcher `attack/match4.c`) enumerates the pools and matches the four lists.
3. It found 16 solutions among the 107 most likely partial-sum values, in 51 s.

## Reproduction

| | |
|---|---|
| Code | commit `90ab79d` of this branch (the options used behave the same at later commits) |
| Command | `python3 attack/skew_ktree.py 287 16 --maxtau 400 --sols 2` |
| Randomness and caps | complete pools (no sampling, no seed); clamp chosen automatically (S = 4); top 6 values per class; at most 400 values of τ |
| Result of that run | pattern (3, 7, 33): 16 solutions, at τ #106 |

The tools write their hit logs, temporary list files and compiled helpers to `$HADAMARD_WORK` (default
`<tmp>/hadamard-work`), never into the repository. `attack/skew_certify.py` (or `skew_group.py certify` for the
non-cyclic orders) turns a logged hit into the files in this directory. The certificate does not depend on
the search: `verify/verify.py` rebuilds and checks the matrix from `certificate.json` alone.

## Status

To our knowledge this is the first skew-Hadamard matrix of order 1148.

- Cati and Pasechnik, *A database of constructions of Hadamard matrices*, arXiv:2411.18897v2 (Aug. 2025), Table 5, lists 287(4). So the smallest known skew order of the form 2ᵐ · 287 was 4592, and 1148 was unknown.
- A literature search in October 2026 found no later construction.

The doubling construction (Seberry) then also gives skew-Hadamard matrices of orders 2296.

**Novelty checklist** (audit of 2026-10-12, by independent agents that also re-verified the matrix):

1. Table 5 of Cati–Pasechnik lists 287(4): the smallest known skew order with odd part 287 was 2⁴ · 287 = 4592. This order is inside the range (≤ 1208) that they checked against implemented constructions.
2. Đoković's papers (1992–2024, alone and with Kotsireas, Golubitsky or Balonin), Koukouvinos–Stylianou
   (Discrete Math. 2008, Tables 1 and 3) and the Đoković–Seberry census of skew orders 4v for odd v < 2500
   (May 2025) do not construct it. The census lists 287 among its 737 undecided values.
3. The 2026 public deposits (Schumer; Alpöge et al.; JD-Jones-ASES; Schneider; Karoui; Piza) do not contain
   it. SageMath's `unknown_skew_hadamard` list (develop branch, fetched 2026-10-12) contains 1148.
4. No infinite family reaches it: Paley I, Szekeres 2(q + 1), Wallis–Whiteman, Spence, Momihara–Xiang,
   the product theorem n(pᵏ + 1) and amicable products were checked.

This is a finding of nothing, not a proof that no earlier construction exists.

## Verification

```
python3 verify/verify.py skew/1148
gzip -dc skew/1148/H1148.csv.gz | ./check_matrix
```

SHA-256 of `H1148.csv`:
`9528531de381c5d5044b38503b855aa0a5622efda205da58f330e57dcb7701b9`
