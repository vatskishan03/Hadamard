# Skew-Hadamard matrix of order 1804

1804 = 4 · 451, with 451 = 11 · 41.

## Construction

A skew-Hadamard matrix is a Hadamard matrix H with H + Hᵀ = 2I.

**Difference family.** Four ±1 sequences x₁, …, x₄ on ℤ₄₅₁ satisfy two conditions:

- each is constant on the orbits of the multiplier subgroup
  H = ⟨16, 42⟩ ≤ ℤ₄₅₁*, of order 25: {1, 16, 37, 42, 59, 78, 92, 100, 119, 124, 133, 141, 174, 180, 201, 221, 223, 247, 256, 262, 324, 344, 346, 379, 411};
- their periodic autocorrelations sum to zero at every nonzero shift.

`certificate.json` lists, for each sequence, one representative r of each orbit rH on which the sequence is −1.

**Skew first block.** The first sequence is of skew type: x₁(0) = 1 and x₁(−t) = −x₁(t) for t ≠ 0. This is possible because −1 ∉ H, so negation pairs the orbits.

**Row sums.** The row sums are 1, 1, 11, 41, and 1² + 1² + 11² + 41² = 1804.

**The matrix.** Let A, B, C, D be the circulant matrices with these first rows (`sequences.txt`) and R the back-diagonal identity. The Goethals–Seidel array

```
    [  A     BR     CR     DR  ]
H = [ -BR    A     DᵀR   -CᵀR ]
    [ -CR   -DᵀR    A     BᵀR ]
    [ -DR    CᵀR  -BᵀR    A   ]
```

is then a Hadamard matrix of order 1804. Because A + Aᵀ = 2I, it satisfies H + Hᵀ = 2I (Wallis–Whiteman).

## Search

1. **Structure.** H = ⟨16, 42⟩ ≤ ℤ₄₅₁* ≅ ℤ₁₀ × ℤ₄₀ is non-cyclic of order 25: the order-5 subgroup of ℤ₁₁* times the order-5 subgroup of ℤ₄₁*. It has 27 orbits and 13 autocorrelation classes; the exact compression screen (`attack/census_screen.py`) keeps the rows (1, 11, 41) alive.
2. **Exact compression chain** (`attack/chain.py`). A dynamic program lists every joint compression of each block to ℤ₁₁; the four lists are joined exactly (60 level-1 branches). In each branch, every sign vector with the branch's compressions is listed by meet in the middle over the orbits, and the four lists are joined on all 13 classes by the exhaustive hash-probe joiner `attack/join4x.c`.
3. The first four branches held no solution; the fifth held 16, and the first two rebuilt ones verified. The same exhaustive search **rules out** the other surviving pattern of this structure, rows (1, 29, 31) (20 branches, 12 minutes), and the structure 623/⟨2⟩ with rows (1, 5, 21, 45). The Gaussian model had given these structures 10^1.5, 10^1.1 and 10^1.6 solutions.

## Reproduction

| | |
|---|---|
| Code | commit `f2e5a51` |
| Command | `python3 attack/chain.py 451 16,42 1,11,41 --budget 30000` |
| Randomness and caps | none (exact search); default caps --maxfibre 3e6, --maxstates 2e7 |
| Result of that run | found in level-1 branch 5 of 60 (16 solutions counted) after 241 s on 4 shared cores |

The tools write hit logs, temporary files and compiled helpers to `$HADAMARD_WORK` (default `<tmp>/hadamard-work`). The certificate does not depend on the search: `verify/verify.py` rebuilds and checks the matrix from `certificate.json` alone.

## Status

To our knowledge this is the first skew-Hadamard matrix of order 1804.

- The Đoković–Seberry census of skew orders 4v for odd v < 2500 (May 2025) lists 451 among its 737 undecided values.
- Seberry–Yamada (1992), Table 7.1, and Koukouvinos–Stylianou (Discrete Math. 2008), Table 1, both list 451(3): the
  smallest known skew order with odd part 451 was 2³ · 451 = 3608, from Paley I (3607 is prime).
- **Caveat.** Cati–Pasechnik's Table 5 (arXiv:2411.18897, v1 and v2) has no entry for 451. By that table's convention a
  missing value means m = 2, i.e. that order 1804 is known. We found nothing behind this: the paper gives no construction
  or citation for it, 1804 lies outside their verified range (orders ≤ 1208), the census lists 451 as undecided, and the
  older tables give 451(3). We read the omission as a transcription gap, not as a known matrix.
- A literature search in October 2026 found no later construction.

## Verification

```
python3 verify/verify.py skew/1804
gzip -dc skew/1804/H1804.csv.gz | ./check_matrix
```

SHA-256 of `H1804.csv`:
`1f9121d6103043e3e6dcfc406f305d4fb4418d1ee8a3e1ddf5224c8f275e54ce`

**Novelty checklist** (audit of 2026-10-12 by a research agent):

1. Census (May 2025): 451 is undecided, and it is in none of the census's Paley, Szekeres, Szekeres–Whiteman or Spence
   series or its 41 sporadic values. Seberry–Yamada 1992 and Koukouvinos–Stylianou 2008: 451(3). Cati–Pasechnik
   Table 5: no entry (see the caveat above).
2. SageMath (develop, 2026-10-11) has no route to skew order 1804: no Goethals–Seidel or difference-family data for
   v = 451, no Paley or Spence parameters, and 8 ∤ 1804 rules out doubling.
3. The 2023–2026 papers and deposits we checked do not contain it:
   - arXiv;
   - Zenodo: Schumer, Piza (whose 1804 entry, "CW(11, 41)", is an ordinary Cooper–Wallis matrix), Kline and Bedani;
   - the JD-Jones-ASES repositories;
   - Alpöge et al. and Schneider.
4. No infinite family reaches it:
   - Paley I (1803 = 3 · 601);
   - Szekeres and Momihara–Xiang 2(q + 1), with q = 901 = 17 · 53;
   - Wallis–Whiteman 4(q + 1), with q = 450;
   - Spence 4(q² + q + 1), (p − 1)ᵘ + 1, (|t| + 1)(q + 1) and mn(n − 1);
   - products and skew × amicable products (every split leaves the odd cofactor 11, 41 or 451).

Not checked: the Seberry–Yamada 2020 book (Table 9.2) and the paywalled journals of 2025–2026. This is a finding of
nothing, not a proof that no earlier construction exists.
