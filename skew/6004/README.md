# Skew-Hadamard matrix of order 6004

6004 = 4 · 1501, with 1501 = 19 · 79.

## Construction

A skew-Hadamard matrix is a Hadamard matrix H with H + Hᵀ = 2I.

**Difference family.** Four ±1 sequences x₁, …, x₄ on ℤ₁₅₀₁ satisfy two conditions:

- each is constant on the orbits of the multiplier subgroup
  H = ⟨799⟩ ≤ ℤ₁₅₀₁*, of order 39: {1, 20, 115, 134, 210, 248, 286, 324, 362, 381, 400, 457, 476, 495, 514, 571, 704, 742, 761, 799, 894, 913, 970, 1046, 1065, 1103, 1122, 1179, 1198, 1217, 1236, 1274, 1331, 1369, 1388, 1407, 1426, 1445, 1464};
- their periodic autocorrelations sum to zero at every nonzero shift.

`certificate.json` lists, for each sequence, one representative r of each orbit rH on which the sequence is −1.

**Skew first block.** The first sequence is of skew type: x₁(0) = 1 and x₁(−t) = −x₁(t) for t ≠ 0. This is possible because −1 ∉ H, so negation pairs the orbits.

**Row sums.** The row sums are 1, 77, 7, −5, and 1² + 77² + 7² + (−5)² = 6004.

**The matrix.** Let A, B, C, D be the circulant matrices with these first rows (`sequences.txt`) and R the back-diagonal identity. The Goethals–Seidel array

```
    [  A     BR     CR     DR  ]
H = [ -BR    A     DᵀR   -CᵀR ]
    [ -CR   -DᵀR    A     BᵀR ]
    [ -DR    CᵀR  -BᵀR    A   ]
```

is then a Hadamard matrix of order 6004. Because A + Aᵀ = 2I, it satisfies H + Hᵀ = 2I (Wallis–Whiteman).

## Search

1. **Structure: a Paley product.** 1501 = 19 · 79 and H = {h : h ≡ 1 (mod 19), h a square mod 79} = ⟨799⟩, of order 39. Every block is u ⊗ 1 + w ⊗ J + z ⊗ χ over ℤ₁₉ × ℤ₇₉ (χ the Legendre sequence of ℤ₇₉, J the all-ones sequence, u, w, z sequences on ℤ₁₉); see `attack/NOTES.md` §15.
2. **Reduction.** The compressions to ℤ₁₉ satisfy Σ_k Σ_a y_k(a)² = 4 · 19 · 79, which reads 80 N_w − 4 N₋ = 76: exactly one position carries w ≠ 0 (N_w = 1, N₋ = 1, q = 4m + 3). Because 4m < 2q, the four u-sequences and the four z-sequences must each have autocorrelations summing to zero, the u-sequence of the block with the w-position is skew-type away from 0, and the cross-correlations Σ_k (z_k ⋆ u_k − u_k ⋆ z_k) vanish. The row-sum identity then needs 2m − 1 = 37 to be a sum of two squares (37 = 1 + 36); the same test rules out the census values 473, 517, 1139 and 1207 in this structure.
3. **Exact two-stage search** (`attack/paley_config.py`). Stage 1 lists all 95,552,064 u-quadruples (26 s). Stage 2 joins, for the first of them, every z-quadruple on the 18 remaining classes (`attack/join4x.c`, about 1.4·10¹¹ probe pairs): 24 solutions. The first two were lifted to ℤ₁₅₀₁ and verified.

## Reproduction

| | |
|---|---|
| Code | commit `520f0b6` |
| Command | `python3 attack/paley_config.py 1501 --maxu 3` |
| Randomness and caps | none in the search (exact enumeration); which solutions join4x prints first depends on the OpenMP schedule |
| Result of that run | 95,552,064 u-quadruples; the first extends to 24 z-solutions; 2 verified skew-Hadamard matrices after 35 minutes on 4 heavily shared cores |

The tools write hit logs, temporary files and compiled helpers to `$HADAMARD_WORK` (default `<tmp>/hadamard-work`). The certificate does not depend on the search: `verify/verify.py` rebuilds and checks the matrix from `certificate.json` alone.

## Status

To our knowledge this is the first skew-Hadamard matrix of order 6004.

- The Đoković–Seberry census of skew orders 4v for odd v < 2500 (May 2025) lists 1501 among its 737 undecided values.
- Every other skew table we found stops below 1501 (Cati–Pasechnik Table 5 at 999).
- Before this, the smallest known skew order with odd part 1501 was 2³ · 1501 = 12008: Paley I applies because
  12007 is prime and ≡ 3 (mod 4), and Wallis–Whiteman 4(q + 1) applies with q = 3001.
- A *symmetric* Hadamard matrix of order 6004 is classical (Paley II with q = 3001), which is why 1501 appears in some
  tables of ordinary and Williamson-type orders. It is not skew.
- A literature search in October 2026 found no later construction.

Doubling also gives every order 2ᵗ · 6004, but only 6004 itself is new: 12008 and its doublings were already known.

**Novelty checklist** (audit of 2026-10-13 by a research agent):

1. The Đoković–Seberry census (PDF of 2025-05-28) lists 1501 among its 737 undecided odd v < 2500, and not in its
   Paley, Szekeres, Szekeres–Whiteman or Spence series or among its 41 sporadic values. Its undecided and known lists
   together cover all 1250 odd v < 2500, with no overlap.
2. Cati–Pasechnik Table 5 (arXiv:2411.18897, v1 and v2) stops at n = 999; neither version contains 1501 or 6004.
   Seberry–Yamada (1992) Table 7.1 stops at 999, Koukouvinos–Stylianou (2008) at 493, Seberry–Yamada (Utilitas Math.
   1991) at 499.
3. SageMath (develop branch, 2026-10-11): `skew_hadamard_matrix(6004, existence=True)` returns Unknown. Paley I fails
   (6003 = 3² · 23 · 29), the Spence 2(q + 1) branch needs q ≡ 5 (mod 8) (q = 3001 ≡ 1), there is no stored skew
   supplementary difference set or good-matrix data for v = 1501, and 6004 ≡ 4 (mod 8) skips the doubling branch.
4. The 2023–2026 papers and deposits we checked do not contain it:
   - arXiv (26 full texts searched for 6004, 1501, 3001 and 12008);
   - Zenodo (165 Hadamard records; files read for Schumer, Piza, Kline and Bedani);
   - the JD-Jones-ASES repositories and jeff-kline/order18-jacobi-triples;
   - Alpöge et al. (12 August 2026) and Schneider (order 2060).
5. No infinite family reaches it:
   - Paley I (6003 = 3² · 23 · 29);
   - Szekeres, Spence and Szekeres–Whiteman 2(q + 1) with q = 3001 ≡ 1 (mod 8), and the Momihara–Xiang two-block
     family (needs q ≡ 5 (mod 8); q = 3001 ≡ 9 (mod 16) gives 4(q + 1) = 12008 instead);
   - Wallis–Whiteman and Momihara–Xiang 2ᵘ(qᵉ + 1), with qᵉ = 1500;
   - Spence 4(q² + q + 1) (6001 is not a square) and Spence G-matrices;
   - (p − 1)ᵘ + 1, (|t| + 1)(q + 1), mn(n − 1) and 2q(q − 1);
   - products 2ᵗ ∏(qᵢ + 1), skew × amicable products (6004 has only two factors of 2);
   - good matrices (known only for small orders) and best matrices (6001 is not a square).

Not checked: Seberry–Yamada (2020) Table 9.2 directly (described by Cati–Pasechnik as covering n ≤ 999), the Handbook
of Combinatorial Designs Table 1.51, paywalled journals of 2025–2026, and posts on social media. This is a finding of
nothing, not a proof that no earlier construction exists.

## Verification

```
python3 verify/verify.py skew/6004
gzip -dc skew/6004/H6004.csv.gz | ./check_matrix
```

SHA-256 of `H6004.csv`:
`f08cd0a71e4f45119f3c0ebbb63360e93c0a02a8d9574e83d484c996392b4ae8`
