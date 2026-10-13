# Skew-Hadamard matrix of order 820

820 = 4 · 205, with 205 = 5 · 41.

## Construction

A skew-Hadamard matrix is a Hadamard matrix H with H + Hᵀ = 2I.

**Difference family.** Four ±1 sequences x₁, …, x₄ on ℤ₂₀₅ satisfy two conditions:

- each is constant on the orbits of the multiplier subgroup H = ⟨51⟩ = {1, 16, 51, 141, 201} ≤ ℤ₂₀₅*, of order 5:
  H = {h : h ≡ 1 (mod 5), h mod 41 an eighth power};
- their periodic autocorrelations sum to zero at every nonzero shift.

`certificate.json` lists, for each sequence, one representative r of each orbit rH on which the sequence is −1.

**Skew first block.** The first sequence is of skew type: x₁(0) = 1 and x₁(−t) = −x₁(t) for t ≠ 0. This is
possible because −1 ∉ H, so negation pairs the orbits.

**Row sums.** The row sums are 1, 1, −17, −23, and 1² + 1² + 17² + 23² = 820.

**The matrix.** Let A, B, C, D be the circulant matrices with these first rows (`sequences.txt`) and R the
back-diagonal identity. The Goethals–Seidel array

```
    [  A     BR     CR     DR  ]
H = [ -BR    A     DᵀR   -CᵀR ]
    [ -CR   -DᵀR    A     BᵀR ]
    [ -DR    CᵀR  -BᵀR    A   ]
```

is then a Hadamard matrix of order 820. Because A + Aᵀ = 2I, it satisfies H + Hᵀ = 2I (Wallis–Whiteman).

## Search

1. **Structure: a cyclotomic product** (`attack/NOTES.md` §16). Over ℤ₅ × GF(41), H acts on GF(41) as the eighth
   powers C₀ (index e = 8, class size F = 5) and trivially on ℤ₅. A block is given at each position a ∈ ℤ₅ by its value
   f(a) at (a, 0) and its values on the eight classes (a, C_j): 45 H-orbits per block.
2. **Stage 1: the compressions to ℤ₅.** With s(a) the sum of the eight class values at a, the compression is
   y(a) = f(a) + 5 s(a), and Σ_k Σ_a y_k(a)² = 4 · 205 becomes Σ over the 20 positions of (5 s² + 2 f s) = 160.
   A spike (s ≠ 0) costs 16 to 88, so up to ten spikes are possible. There are 1084 stage-1 classes up to negation,
   translation and reversal of the non-skew blocks. The search uses the first:
   y = (1, 1, 1, −1, −1), (−1, 1, 1, −1, 1), (−9, 1, −9, 1, −1), (−21, −11, 1, −1, 9).
3. **Lists.** The blocks of this class have 78,400 (skew block), 1,680,700,000, 1,075,648,000 and 430,259,200
   members. Earlier sweeps skipped the structure because these lists are far too large to store. The Gaussian
   estimate of the class (in the 20 level-B/level-C coordinates, rank 19) is 10^8.47 solutions.
4. **Low-energy k-tree** (`attack/cp_lowe.py`, `attack/cp_enum.c`, `attack/kt4.c`; `attack/NOTES.md` §20).
   - `cp_enum.c` enumerates every member of each non-skew block. It keeps the members whose autocorrelation vector
     is close to its conditional mean given that the four vectors sum to zero (whitened energy): 5,156,254,
     4,774,328 and 4,908,954 members. The Gaussian estimate inside these lists is 10^3.84.
   - `kt4.c` joins the lists over 3000 split values of five merge classes, ranked by expected solutions per unit
     of work. The first split value gave 6 solutions; the run stopped at 20 (all verified).

## Reproduction

| | |
|---|---|
| Code | commit `15deec6` (`attack/cp_lowe.py`, `attack/cp_enum.c`, `attack/kt4.c`); the original run used the uncommitted code that became commit `669dc2f` (same algorithm and seed) |
| Command | `python3 attack/cp_lowe.py 205 --q 41 --e 8 --config 0 --L 5e6 --nsig 3000 --threads 3 --K 20` |
| Randomness and caps | `--seed 1` (default; numpy `default_rng`) draws the 2·10⁵ samples per block that set the energy thresholds and the 4·10⁶ pairs per side that rank the split values; `--L 5e6` members kept per non-skew block; `--K 20` stops kt4 after 20 solutions |
| Result of that run | 20 verified solutions in 268 s on 3 threads (block enumeration 182 s, join 40 s); the certificate is the first hit recorded. A rerun with commit `15deec6` keeps the same 5,156,254 / 4,774,328 / 4,908,954 members and again gives 20 verified solutions (597 s on 3 shared threads), the certificate's among them; kt4's order of output depends on the thread schedule |

The tools write hit logs, temporary files and compiled helpers to `$HADAMARD_WORK` (default `<tmp>/hadamard-work`).
The certificate does not depend on the search: `verify/verify.py` rebuilds and checks the matrix from
`certificate.json` alone.

## Status

To our knowledge this is the first skew-Hadamard matrix of order 820.

- The Đoković–Seberry census of skew orders 4v for odd v < 2500 (PDF of 2025-05-28) lists 205 among its 737 undecided
  values.
- Cati–Pasechnik Table 5 (arXiv:2411.18897, v1 and v2, identical tables) gives 205(3): the smallest known skew order
  with odd part 205 was 8 · 205 = 1640, from Wallis–Whiteman with q = 409 ≡ 9 (mod 16).
- An *ordinary* (and a symmetric) Hadamard matrix of order 820 is classical (Paley II with q = 409), so 820 was never
  an open Hadamard order. Only the skew property is new.

**Novelty checklist** (audit of 2026-10-13 by a research agent):

1. **The Đoković–Seberry census** (PDF of 2025-05-28; HADAMARD 2025 workshop, Sevilla)
   - It lists 205 among the undecided v.
   - 205 is in none of its Paley-3, Szekeres, Szekeres–Whiteman or Spence series, nor among its 41 sporadic values.
   - No later version was found, on arXiv or on the workshop pages.
2. **Cati–Pasechnik** (arXiv:2411.18897, v1 of 2024-11-28 and v2 of 2025-08-30)
   - Table 5 gives 205(3), in plain font. Table 3 (skew constructions) is blank at 205.
   - A commented-out list of unknown skew orders in the TeX source includes 820.
3. **Surveys and Đoković's papers.** Every one whose range covers 205 has it as unknown or with exponent 3:
   - Seberry–Yamada 1992, Table 7.1: "205 3".
   - Koukouvinos–Stylianou 2008, Table 1: 205(3).
   - Đoković's lists of unknown skew orders: Australas. J. Combin. 10 (1994); arXiv:0704.0640; arXiv:0706.1973;
     Đoković–Kotsireas arXiv:1511.08734 ("201, 205, 209, …").
   - arXiv:1301.3671, 2301.02751, 2404.14375, 2508.17141, 1802.00556, 1811.05094 and 1907.04987 do not construct it.
4. **SageMath** develop (2026-10-13)
   - The module lists 820 in `unknown_skew_hadamard`, and `skew_hadamard_matrix(820, existence=True)` returns Unknown
     (passagemath 10.8.13, same code).
   - It also returns Unknown for 1640.
5. **2026 deposits**: none has a skew matrix of order 820.
   - Checked: Alpöge et al., Schneider, JD-Jones-ASES, Schumer (skew 3676), Piza, Kline, Bedani, Karoui, Blackthorn.
   - Blackthorn's Legendre pairs of length 205 give an ordinary matrix of order 412.
   - A sweep of all 910 Zenodo records matching "hadamard" and GitHub code searches found nothing.
6. **arXiv**: the sources of 148 Hadamard-related papers (2024-01 to 2026-10-08) and the math.CO listings to
   2026-10-12 use 820 or 205 in a skew sense only in Cati–Pasechnik.
7. **No classical construction gives skew 820**:
   - Paley I: 819 is not a prime power. No product 2ᵗ ∏(qᵢ + 1), since 16 ∤ 820.
   - Szekeres 2(q + 1): q = 409 ≡ 1 (mod 8). Wallis–Whiteman 4(q + 1) needs q = 204, which is not a prime power;
     q = 409 gives 1640.
   - Spence and Whiteman families; Seberry–Yamada SII, SVIII and SX; skew × amicable products.
   - Seberry's mn(n − 1) would need an OD(820; 1, 410, 409).
   - Good matrices and skew Williamson-type matrices of order 205 are not recorded.

Not checked: the Handbook (2007) skew table (paywalled), the Seberry–Yamada 2020 book, and paywalled 2025–2026
journals beyond their metadata. This is a finding of nothing, not a proof that no earlier construction exists.

Doubling (Seberry) gives skew orders 1640 (already known, see above), 3280, …

## Files

| File | Contents |
|---|---|
| `certificate.json` | the four sequences as orbit representatives, H, the row sums |
| `sequences.txt` | the four sequences, one per line |
| `H820.csv.gz` | the matrix |

## Verification

```
python3 verify/verify.py skew/820
gzip -dc skew/820/H820.csv.gz | ./check_matrix
```

SHA-256 of `H820.csv`: `8bb20ef554a33c7f179843831dc125384c38a08e33fc69c2e2bb8de4f2e81774`
