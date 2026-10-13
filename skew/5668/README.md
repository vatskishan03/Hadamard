# Skew-Hadamard matrix of order 5668

5668 = 4 · 1417, with 1417 = 13 · 109.

## Construction

A skew-Hadamard matrix is a Hadamard matrix H with H + Hᵀ = 2I.

**Difference family.** Four ±1 sequences x₁, …, x₄ on ℤ₁₄₁₇ satisfy two conditions:

- each is constant on the orbits of the multiplier subgroup
  H = ⟨1405, 328⟩ ≤ ℤ₁₄₁₇*, of order 81: H = {h : h mod 13 ∈ {1, 3, 9}, h mod 109 a fourth power};
- their periodic autocorrelations sum to zero at every nonzero shift.

`certificate.json` lists, for each sequence, one representative r of each orbit rH on which the sequence is −1.

**Skew first block.** The first sequence is of skew type: x₁(0) = 1 and x₁(−t) = −x₁(t) for t ≠ 0. This is possible because −1 ∉ H (−1 mod 13 = 12 ∉ {1, 3, 9}), so negation pairs the orbits.

**Row sums.** The row sums are 1, 7, −53, −53, and 1² + 7² + (−53)² + (−53)² = 5668.

**The matrix.** Let A, B, C, D be the circulant matrices with these first rows (`sequences.txt`) and R the back-diagonal identity. The Goethals–Seidel array

```
    [  A     BR     CR     DR  ]
H = [ -BR    A     DᵀR   -CᵀR ]
    [ -CR   -DᵀR    A     BᵀR ]
    [ -DR    CᵀR  -BᵀR    A   ]
```

is then a Hadamard matrix of order 5668. Because A + Aᵀ = 2I, it satisfies H + Hᵀ = 2I (Wallis–Whiteman).

## Search

1. **Structure: a cyclotomic product with a multiplier on ℤ₁₃** (`attack/NOTES.md` §16). Over ℤ₁₃ × GF(109),
   H acts on GF(109) as the fourth powers C₀ (index e = 4, F = 27) and on ℤ₁₃ as K₁₃ = ⟨3⟩ = {1, 3, 9}. A block is
   given at each position a ∈ ℤ₁₃ by its value f(a) at (a, 0) and its values on the four classes (a, C_j); the
   positions fall into the K₁₃-orbits {0}, {1, 3, 9}, {2, 5, 6}, {4, 10, 12} and {7, 8, 11}. There are 25 H-orbits
   per block.
2. **Stage 1: the compressions to ℤ₁₃.** With s(a) the sum of the four class values at a, the compression is
   y(a) = f(a) + 27 s(a), and Σ_k Σ_a y_k(a)² = 4 · 1417 becomes Σ over the 52 positions of (27 s² + 2 f s) = 208.
   A position with s = ±2 costs at least 104, and K₁₃ forces a spike off 0 to come with its two images, so there
   are exactly two spikes, both at position 0, each with y(0) = ∓53. This leaves 12 stage-1 classes up to
   equivalence (row sums (1, ±7, −53, −53)).
3. **Exact join** (`attack/cyclo_product.py` with `attack/join3.c`). For each class the four lists have 144 (skew
   block), 7,776 and 5,184 (twice) members, grouped by the autocorrelations of their compressions to GF(109) (level B),
   then joined on the remaining shifts (level C). The complete join of all 12 classes finds 8 solutions in each of the
   classes 2, 3, 6 and 11 and none in the others, in about 1 s per class.
4. **How it was found.** `attack/skew_kmult_sweep.py` runs stage 1 for every census value v = m q (q prime), every
   even e with (q − 1)/e odd and every cyclic K_m ≤ ℤ_m* without −1 whose structure has at most 22 position orbits
   and at most 2·10⁷ members per block: 210 structures, 8 with stage-1 classes, all 8 joined. 1417 was the only hit.

## Reproduction

| | |
|---|---|
| Code | commit `2281188` (`attack/cyclo_product.py`); the sweep script is `attack/skew_kmult_sweep.py` (commit `c789279`) |
| Command | `python3 attack/cyclo_product.py 1417 --q 109 --e 4 --kgen 3 --config 2 --K 1 --sub 8000 --threads 2` |
| Randomness and caps | none (exact enumeration; `--sub 8000` keeps every level-B group here, since no list has more than 246 groups); which solution join3 reports first depends on the thread schedule. `--K 1000` instead of `--K 1` counts all solutions of the class |
| Result of that run | 12 stage-1 classes; class 2 has lists of 144, 7,776, 5,184 and 5,184 members and gives the solution of the certificate; under 1 s on 2 cores |

The tools write hit logs, temporary files and compiled helpers to `$HADAMARD_WORK` (default `<tmp>/hadamard-work`). The certificate does not depend on the search: `verify/verify.py` rebuilds and checks the matrix from `certificate.json` alone.

## Status

To our knowledge this is the first skew-Hadamard matrix of order 5668.

- The Đoković–Seberry census of skew orders 4v for odd v < 2500 (May 2025) lists 1417 among its 737 undecided values.
- Every other skew table we found stops below 1417 (Cati–Pasechnik Table 5 at 999).
- Before this, the smallest known skew order with odd part 1417 was 2⁵ · 1417 = 45344 (Paley I: 45343 is prime and
  ≡ 3 (mod 4)). We found no construction of 11336 or 22672 either.
- A *symmetric* Hadamard matrix of order 5668 is classical (Paley II with q = 2833), so 5668 was never an open order
  for ordinary Hadamard matrices. It is not skew.
- A literature search in October 2026 found no later construction.

Doubling (Seberry) then also gives skew-Hadamard matrices of orders 11336 and 22672, which are new as far as we
know; 45344 and its doublings were already known.

**Novelty checklist** (audit of 2026-10-13 by a research agent):

1. The Đoković–Seberry census (PDF of 2025-05-28) lists 1417 among its 737 undecided odd v < 2500, and not in its
   Paley, Szekeres, Szekeres–Whiteman or Spence series or among its 41 sporadic values. Its undecided and known lists
   together cover all 1250 odd v < 2500, with no overlap.
2. Cati–Pasechnik Table 5 and Seberry–Yamada (1992) Table 7.1 stop at 999, Seberry–Yamada (Utilitas Math. 1991) at
   499 and Koukouvinos–Stylianou (2008) at 493.
3. SageMath (passagemath 10.8.13; develop branch, 2026-10-11): `skew_hadamard_matrix(n, existence=True)` returns
   Unknown for n = 5668, 11336 and 22672.
4. The 2023–2026 papers and deposits we checked do not contain it:
   - arXiv (e-print sources of 587 Hadamard-related papers, 2000 to 2026-10-08);
   - Zenodo (about 1,500 records; files read for Schumer, Piza, Kline and Bedani);
   - 13 GitHub repositories, among them the JD-Jones-ASES repositories and jeff-kline/order18-jacobi-triples;
   - Alpöge et al. (12 August 2026) and Schneider (order 2060);
   - Balonin's catalogue, the Palomar registry and the OEIS.
5. No infinite family reaches it:
   - Paley I (5667 = 3 · 1889);
   - Szekeres, Whiteman, Spence and Momihara–Xiang 2(q + 1), with q = 2833 ≡ 1 (mod 8);
   - Wallis–Whiteman and Momihara–Xiang 2ᵘ(qᵉ + 1), with qᵉ = 1416 (and q = 2833 ≡ 1 (mod 16) for 11336);
   - Spence 4(q² + q + 1) (5665 is not a square);
   - (p − 1)ᵘ + 1, (|t| + 1)(q + 1), 2ᵗq and mn(n − 1);
   - skew × amicable products (5668 has only two factors of 2), good matrices and best matrices.

   We found no product theorem that would combine skew matrices of orders 52 and 436 into one of order 11336:
   Seberry (1978) and Seberry–Yamada (1992) state that Kronecker products of skew matrices are not skew, and their
   only product is skew × amicable.

Not checked: Seberry–Yamada (2020), the Handbook of Combinatorial Designs skew table (v < 200), Koukouvinos–Stylianou
(2008) in full, and paywalled journals of 2025–2026. This is a finding of nothing, not a proof that no earlier
construction exists.

## Verification

```
python3 verify/verify.py skew/5668
gzip -dc skew/5668/H5668.csv.gz | ./check_matrix
```

SHA-256 of `H5668.csv`:
`0445a99fdbbee9be8e319ca2b4a277099110bb47073758f1f0ff0866ab747f0c`
