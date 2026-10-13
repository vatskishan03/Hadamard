# Skew-Hadamard matrix of order 1204

1204 = 4 · 301, with 301 = 7 · 43.

## Construction

A skew-Hadamard matrix is a Hadamard matrix H with H + Hᵀ = 2I.

**Difference family.** Four ±1 sequences x₁, …, x₄ on ℤ₃₀₁ satisfy two conditions:

- each is constant on the orbits of the multiplier subgroup H = ⟨127⟩ = {1, 64, 78, 127, 176, 183, 274} ≤ ℤ₃₀₁*, of order 7:
  H = {h : h ≡ 1 (mod 7), h mod 43 a sixth power};
- their periodic autocorrelations sum to zero at every nonzero shift.

`certificate.json` lists, for each sequence, one representative r of each orbit rH on which the sequence is −1.

**Skew first block.** The first sequence is of skew type: x₁(0) = 1 and x₁(−t) = −x₁(t) for t ≠ 0. This is
possible because −1 ∉ H, so negation pairs the orbits.

**Row sums.** The row sums are 1, −11, −11, −31, and 1² + 11² + 11² + 31² = 1204.

**The matrix.** Let A, B, C, D be the circulant matrices with these first rows (`sequences.txt`) and R the
back-diagonal identity. The Goethals–Seidel array

```
    [  A     BR     CR     DR  ]
H = [ -BR    A     DᵀR   -CᵀR ]
    [ -CR   -DᵀR    A     BᵀR ]
    [ -DR    CᵀR  -BᵀR    A   ]
```

is then a Hadamard matrix of order 1204. Because A + Aᵀ = 2I, it satisfies H + Hᵀ = 2I (Wallis–Whiteman).

## Search

1. **Structure: a cyclotomic product** (`attack/NOTES.md` §16). Over ℤ₇ × GF(43), H acts on GF(43) as the classes of index e = 6
   (class size F = 7) and trivially on ℤ₇. A block is given at each position a ∈ ℤ₇ by its value f(a) at (a, 0) and its values on the six classes (a, C_j): 49 H-orbits per block.
2. **Stage 1: the compressions to ℤ₇.** The compression is y(a) = f(a) + 7 s(a), and Σ_k Σ_a y_k(a)² = 4 · 301 becomes Σ over the 28 positions of (7 s² + 2 f s) = 168. Stage 1 is sampled by `cp_lowe.stage1_sample` (seed 1); its first quadruple is y = (1, 1, 1, −1, 1, −1, −1), (−13, −1, −1, 1, 1, 1, 1), (−1, 1, −13, 1, −15, 1, 15), (−1, −13, −1, −1, −15, −1, 1).
3. **Lists.** The blocks of this class have 64,000 (skew block), 960,000,000, 540,000,000 and 720,000,000 members. The Gaussian estimate of the class (21 coordinates, rank 20) is 10^6.03 solutions.
4. **Low-energy k-tree** (`attack/cp_lowe.py`, `attack/cp_enum.c`, `attack/kt4.c`; `attack/NOTES.md` §20). `cp_enum.c` enumerates every member of each non-skew block and keeps those whose autocorrelation vector is close to its conditional mean given that the four vectors sum to zero: 9,912,243, 10,109,883 and 9,927,286 members, with a Gaussian estimate of 10^2.52 solutions inside the lists. `kt4.c` joins them over split values of five merge classes ranked by expected solutions per unit of work; the first solution came at the 24th split value. The structure was selected by the estimate-then-search loop `attack/cp_hunt.py` (§20).

## Reproduction

| | |
|---|---|
| Code | commit `15deec6` (`attack/cp_lowe.py`, `attack/cp_enum.c`, `attack/kt4.c`) |
| Command | `python3 attack/cp_lowe.py 301 --q 43 --e 6 --kgen 1 --twist 0 --sample-stage1 3 --config 0 --L 1e7 --nsig 3000 --threads 3 --K 3` |
| Randomness and caps | `--seed 1` (default) draws the sampled stage-1 quadruples, the 2·10⁵ samples per block for the energy statistics and the 4·10⁶ pairs per side that rank the split values; `--L 1e7` members kept per non-skew block |
| Result of that run | first verified solution after about 6 minutes on 3 threads (enumeration 121 s, split values 1–24 in 172 s); the certificate is the first hit recorded |

The tools write hit logs, temporary files and compiled helpers to `$HADAMARD_WORK` (default `<tmp>/hadamard-work`).
The certificate does not depend on the search: `verify/verify.py` rebuilds and checks the matrix from
`certificate.json` alone.

## Status

To our knowledge this is the first skew-Hadamard matrix of order 1204.

- The Đoković–Seberry census of skew orders 4v for odd v < 2500 (PDF of 2025-05-28) lists 301 among its 737 undecided
  values.
- Cati–Pasechnik Table 5 (arXiv:2411.18897, v1 and v2) gives 301(3): the smallest known skew order with odd part 301 was
  8 · 301 = 2408, from Wallis–Whiteman with q = 601 ≡ 9 (mod 16).
- An *ordinary* Hadamard matrix of order 1204 is classical (Paley II with q = 601; the census also lists 301 in its
  symmetric Paley series). Only the skew property is new.

**Novelty checklist** (audit of 2026-10-13 by a research agent):

1. **The Đoković–Seberry census** (PDF of 2025-05-28, read through a cached copy)
   - It lists 301 among the undecided v ("… 299, 301, 303, …").
   - 301 is in none of its Paley-3, Szekeres, Szekeres–Whiteman or Spence series, nor among its 41 sporadic values.
   - No later version was found.
2. **Cati–Pasechnik** (arXiv:2411.18897, v1 of 2024-11-28 and v2 of 2025-08-30)
   - Table 5 gives 301(3), in plain font. Their Table 3 (skew constructions, n ≤ 1200) ends at n = 299.
   - They state that the first skew order missing from Sage is 1252, and Sage has no 1204.
3. **Surveys and Đoković's papers.** Seberry–Yamada 1992, Table 7.1: "301 3". Koukouvinos–Stylianou 2008, Table 1:
   301(3). The full texts of twelve Đoković and Đoković–Kotsireas papers (0704.0640 to 2508.17141) contain no skew 1204
   or 2408 with odd part 301.
4. **SageMath** (develop, 2026-10-11; passagemath 10.8.13): `skew_hadamard_matrix(1204, existence=True)` returns
   Unknown. Its list of unknown skew orders stops at 1200, so it neither lists nor excludes 1204.
5. **2026 deposits**: none has a skew matrix of order 1204.
   - Checked: Alpöge et al., Schneider, JD-Jones-ASES, Schumer, Piza (its "1204 = 2(601 + 1)" is Paley II), Kline,
     Bedani, Karoui.
   - Zenodo and GitHub code searches for skew 1204 and for 301 found nothing.
6. **arXiv**: 22 Hadamard-related papers of 2024–2026 were searched in full text for 1204, 301 and 2408; the only "301"
   is a block size in arXiv:2404.14375.
7. **No classical construction gives skew 1204**:
   - Paley I: 1203 = 3 · 401. Szekeres 2(q + 1) and Spence's G-matrices need q = 601 ≡ 5 (mod 8), but 601 ≡ 1.
   - Wallis–Whiteman 4(q + 1) needs q = 300; q = 601 gives 2408.
   - Products and doubling need 8 | 1204. Spence 4(q² + q + 1), (p − 1)^u + 1, (|t| + 1)(q + 1), 2ᵗ q and
     Seberry's mn(n − 1) all fail by arithmetic.
   - Good matrices are not known for order 301.

Not checked: the Handbook (2007) Table 1.51, the Seberry–Yamada 2020 book (Table 9.2), Seberry 2017 (Table 4.5) and
paywalled 2025–2026 journals. This is a finding of nothing, not a proof that no earlier construction exists.

Doubling (Seberry) gives skew orders 2408 (already known, see above), 4816, …

## Files

| File | Contents |
|---|---|
| `certificate.json` | the four sequences as orbit representatives, H, the row sums |
| `sequences.txt` | the four sequences, one per line |
| `H1204.csv.gz` | the matrix |

## Verification

```
python3 verify/verify.py skew/1204
gzip -dc skew/1204/H1204.csv.gz | ./check_matrix
```

SHA-256 of `H1204.csv`: `076a687c355073792481e4f498f79d6080e0c46abd492c39b6cb9cb4bbcaa61a`
