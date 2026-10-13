# Skew-Hadamard matrix of order 1044

1044 = 4 · 261, with 261 = 9 · 29.

## Construction

A skew-Hadamard matrix is a Hadamard matrix H with H + Hᵀ = 2I.

**Difference family.** Four ±1 sequences x₁, …, x₄ on ℤ₂₆₁ satisfy two conditions:

- each is constant on the orbits of the multiplier subgroup H = ⟨190⟩ = {1, 82, 136, 181, 190, 199, 226} ≤ ℤ₂₆₁*, of order 7:
  H = {h : h ≡ 1 (mod 9), h mod 29 a fourth power};
- their periodic autocorrelations sum to zero at every nonzero shift.

`certificate.json` lists, for each sequence, one representative r of each orbit rH on which the sequence is −1.

**Skew first block.** The first sequence is of skew type: x₁(0) = 1 and x₁(−t) = −x₁(t) for t ≠ 0. This is
possible because −1 ∉ H, so negation pairs the orbits.

**Row sums.** The row sums are 1, 9, −1, 31, and 1² + 9² + 1² + 31² = 1044.

**The matrix.** Let A, B, C, D be the circulant matrices with these first rows (`sequences.txt`) and R the
back-diagonal identity. The Goethals–Seidel array

```
    [  A     BR     CR     DR  ]
H = [ -BR    A     DᵀR   -CᵀR ]
    [ -CR   -DᵀR    A     BᵀR ]
    [ -DR    CᵀR  -BᵀR    A   ]
```

is then a Hadamard matrix of order 1044. Because A + Aᵀ = 2I, it satisfies H + Hᵀ = 2I (Wallis–Whiteman).

## Search

1. **Structure: a cyclotomic product** (`attack/NOTES.md` §16). Over ℤ₉ × GF(29), H acts on GF(29) as the classes of index e = 4
   (class size F = 7) and trivially on ℤ₉. A block is given at each position a ∈ ℤ₉ by its value f(a) at (a, 0) and its values on the four classes (a, C_j): 45 H-orbits per block.
2. **Stage 1: the compressions to ℤ₉.** The compression is y(a) = f(a) + 7 s(a), with s(a) the sum of the four class values at a, and Σ_k Σ_a y_k(a)² = 4 · 261 becomes Σ over the 36 positions of (7 s² + 2 f s) = 144. Stage 1 has about 3.3 · 10⁶ compressions per block, too many to list all quadruples, so `cp_lowe.stage1_sample` hashes block 3's compressions by cost and autocorrelations and completes random triples of the other blocks (seed 1). Its first quadruple is y = (1, −1, 1, 1, −1, 1, −1, −1, 1), (1, −1, −1, −1, 15, −1, −1, −1, −1), (1, −1, 1, 1, −15, −1, −1, 1, 13), (−1, 1, 1, 15, 1, 1, 1, −1, 13).
3. **Lists.** The blocks of this class have 5,184 (skew block), 6,718,464, 4,478,976 and 4,478,976 members: small enough to keep whole, but the structure had been skipped because its stage 1 could not be listed. The Gaussian estimate of the class (18 level-B/level-C coordinates, rank 17) is 10^3.62 solutions.
4. **Low-energy k-tree** (`attack/cp_lowe.py`, `attack/cp_enum.c`, `attack/kt4.c`; `attack/NOTES.md` §20). `cp_enum.c` enumerates the blocks (all members kept at `--L 1e7`) and `kt4.c` joins them over 3000 split values of five merge classes, ranked by expected solutions per unit of work. The first two solutions appeared within the first 65 split values.

## Reproduction

| | |
|---|---|
| Code | commit `15deec6` (`attack/cp_lowe.py`, `attack/cp_enum.c`, `attack/kt4.c`) |
| Command | `python3 attack/cp_lowe.py 261 --q 29 --e 4 --sample-stage1 5 --config 0 --L 1e7 --nsig 3000 --threads 3 --K 20` |
| Randomness and caps | `--seed 1` (default; numpy `default_rng`) draws the sampled stage-1 quadruples, the 2·10⁵ samples per block for the energy statistics and the 4·10⁶ pairs per side that rank the split values; `--L 1e7` keeps every member here |
| Result of that run | first verified solution after 215 s on 3 threads (enumeration 115 s); the certificate is the first hit recorded |

The tools write hit logs, temporary files and compiled helpers to `$HADAMARD_WORK` (default `<tmp>/hadamard-work`).
The certificate does not depend on the search: `verify/verify.py` rebuilds and checks the matrix from
`certificate.json` alone.

## Status

To our knowledge this is the first skew-Hadamard matrix of order 1044.

- The Đoković–Seberry census of skew orders 4v for odd v < 2500 (PDF of 2025-05-28) lists 261 among its 737 undecided
  values.
- Cati–Pasechnik Table 5 (arXiv:2411.18897, v1 and v2) gives 261(3): the smallest known skew order with odd part 261 was
  8 · 261 = 2088, from Paley I (2087 is prime).
- An *ordinary* (and a symmetric) Hadamard matrix of order 1044 is classical (Paley II with q = 521), so 1044 was never
  an open Hadamard order. Only the skew property is new.

**Novelty checklist** (audit of 2026-10-13 by a research agent):

1. **The Đoković–Seberry census** (PDF of 2025-05-28, 11 pages)
   - It lists 261 among the undecided v ("… 259, 261, 265, …").
   - 261 is in none of its Paley-3, Szekeres, Szekeres–Whiteman or Spence series, nor among its 41 sporadic values.
   - No later version was found (arXiv, Crossref, the workshop pages).
2. **Cati–Pasechnik** (arXiv:2411.18897, v1 of 2024-11-28 and v2 of 2025-08-30)
   - Table 5 gives 261(3), in plain font. Their Table 3 (skew constructions) is empty at 261.
   - 1044 is among the 41 unknown skew orders up to 1208 named in their introduction.
3. **Older tables and Đoković's papers.** Seberry (Wallis), Ars Combin. 6 (1978), Table 2: "261 3". Seberry–Yamada
   1991 and 1992: "261 3" (skew and amicable tables). Koukouvinos–Stylianou 2008, Table 1: 261(3). Đoković's lists of
   unknown skew orders (Australas. J. Combin. 10 (1994); arXiv:0704.0640; arXiv:0706.1973) include 261. Thirteen later
   papers of Đoković, Kotsireas, Bright et al. and others do not construct it.
4. **SageMath** (develop 0835358, stable 10.10, passagemath main): the module lists 1044 in `unknown_skew_hadamard`, and
   no branch of `skew_hadamard_matrix` applies.
5. **2026 deposits**: none has a skew matrix of order 1044.
   - Checked: Alpöge et al., Schneider, JD-Jones-ASES, Schumer, Piza (its "1044 = 2(521 + 1)" is Paley II), Kline,
     Bedani, Karoui.
   - Zenodo metadata and GitHub code searches found only copies of Sage's file, where 1044 is listed as unknown.
6. **arXiv**: 62 full texts of 2024–2026 (Hadamard matrices, skew, Goethals–Seidel, Legendre pairs, Williamson type,
   Hadamard difference sets) do not contain 1044.
7. **No classical construction gives skew 1044**:
   - Paley I: 1043 = 7 · 149. Szekeres 2(q + 1), Spence's G-matrices and Momihara–Xiang need q = 521 ≡ 5 (mod 8), but
     521 ≡ 1.
   - Wallis–Whiteman 4(q + 1) needs q = 260; q = 521 ≡ 9 (mod 16) gives 2088.
   - Spence 4(q² + q + 1), the Seberry–Yamada families SI–SXI, Williamson n(pᵏ + 1), doubling (8 ∤ 1044) and amicable
     products fail by arithmetic. Seberry's mn(n − 1) needs an OD(348; 1, 87, 260), which no source gives.
   - Good, G-, best or skew Williamson-type matrices of order 261 are not recorded. The only groups of order 261 are
     ℤ₂₆₁ and ℤ₃ × ℤ₈₇.

Not checked: the Handbook (2007) Table 1.51, the Seberry–Yamada 2020 book (Table 9.2), Seberry 2017 (Table 4.5) and
paywalled 2025–2026 journals. This is a finding of nothing, not a proof that no earlier construction exists.

Doubling (Seberry) gives skew orders 2088 (already known, see above), 4176, …

## Files

| File | Contents |
|---|---|
| `certificate.json` | the four sequences as orbit representatives, H, the row sums |
| `sequences.txt` | the four sequences, one per line |
| `H1044.csv.gz` | the matrix |

## Verification

```
python3 verify/verify.py skew/1044
gzip -dc skew/1044/H1044.csv.gz | ./check_matrix
```

SHA-256 of `H1044.csv`: `fe17f4ecc085fe9aca23f11f18bff5eae446318455bd2802648aa91e788c2f0e`
