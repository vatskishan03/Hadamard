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

1. **Structure: a cyclotomic product** ([search/NOTES.md](../../search/NOTES.md) §16). Over ℤ₅ × GF(41), H acts on
   GF(41) as the eighth powers C₀ (index e = 8, class size F = 5) and trivially on ℤ₅. A block is given at each
   position a ∈ ℤ₅ by its value f(a) at (a, 0) and its values on the eight classes (a, C_j): 45 H-orbits per block.
2. **Stage 1: the compressions to ℤ₅.** With s(a) the sum of the eight class values at a, the compression is
   y(a) = f(a) + 5 s(a), and Σ_k Σ_a y_k(a)² = 4 · 205 becomes Σ over the 20 positions of (5 s² + 2 f s) = 160.
   A spike (s ≠ 0) costs 16 to 88, so up to ten spikes are possible. There are 1084 stage-1 classes up to negation,
   translation and reversal of the non-skew blocks. The search uses the first:
   y = (1, 1, 1, −1, −1), (−1, 1, 1, −1, 1), (−9, 1, −9, 1, −1), (−21, −11, 1, −1, 9).
3. **Lists.** The blocks of this class have 78,400 (skew block), 1,680,700,000, 1,075,648,000 and 430,259,200
   members. Earlier sweeps skipped the structure because these lists are far too large to store. The Gaussian
   estimate of the class (in the 20 level-B/level-C coordinates, rank 19) is 10^8.47 solutions.
4. **Low-energy k-tree** (`search/cp_lowe.py`, `search/cp_enum.c`, `search/kt4.c`;
   [search/NOTES.md](../../search/NOTES.md) §20).
   - `cp_enum.c` enumerates every member of each non-skew block. It keeps the members whose autocorrelation vector
     is close to its conditional mean given that the four vectors sum to zero (whitened energy): 5,156,254,
     4,774,328 and 4,908,954 members. The Gaussian estimate inside these lists is 10^3.84.
   - `kt4.c` joins the lists over 3000 split values of five merge classes, ranked by expected solutions per unit
     of work. The first split value gave 6 solutions; the run stopped at 20 (all verified).

**Reproduction.**

| | |
|---|---|
| Code | commit `15deec6` (`search/cp_lowe.py`, `search/cp_enum.c`, `search/kt4.c`); the original run used an earlier state of the same code, committed as `669dc2f` (same algorithm and seed) |
| Command | `python3 search/cp_lowe.py 205 --q 41 --e 8 --config 0 --L 5e6 --nsig 3000 --threads 3 --K 20` |
| Randomness and caps | `--seed 1` (default; numpy `default_rng`) draws the 2·10⁵ samples per block that set the energy thresholds and the 4·10⁶ pairs per side that rank the split values; `--L 5e6` members kept per non-skew block; `--K 20` stops kt4 after 20 solutions |
| Result | 20 verified solutions in 268 s on 3 threads (block enumeration 182 s, join 40 s); the certificate is the first hit recorded. A rerun at commit `15deec6` keeps the same 5,156,254 / 4,774,328 / 4,908,954 members and again gives 20 verified solutions, the certificate's among them; the order of kt4's output depends on the thread schedule |

## Status

**New** as a skew-Hadamard matrix. An ordinary (and a symmetric) Hadamard matrix of order 820 is classical (Paley II
with q = 409); the skew property is what is new.

- The Đoković–Seberry census lists 205 among its 737 undecided v. 205 is in none of its Paley, Szekeres,
  Szekeres–Whiteman or Spence series, nor among its 41 sporadic values.
- Cati–Pasechnik Table 5 gives 205(3): the smallest known skew order with odd part 205 was 8 · 205 = 1640, from
  Wallis–Whiteman with q = 409 ≡ 9 (mod 16). Their Table 3 (skew constructions) is blank at 205. Seberry–Yamada (1992,
  Table 7.1) and Koukouvinos–Stylianou (2008, Table 1) also give 205(3), and Đoković and Kotsireas list 205 as open
  (arXiv:1511.08734).
- SageMath lists 820 in `unknown_skew_hadamard` and returns Unknown for 820 and for 1640.
- No classical construction gives a skew matrix of order 820:
  - Paley I: 819 is not a prime power. No product 2ᵗ ∏(qᵢ + 1) applies, since 16 ∤ 820.
  - Szekeres 2(q + 1) needs q ≡ 5 (mod 8), and q = 409 ≡ 1 (mod 8). Wallis–Whiteman 4(q + 1) needs q = 204, which is
    not a prime power.
  - The Spence and Whiteman families, Seberry–Yamada's SII, SVIII and SX, and skew × amicable products do not apply.
  - Seberry's mn(n − 1) construction would need an OD(820; 1, 410, 409).
  - Good matrices and skew Williamson-type matrices of order 205 are not recorded.

The sources searched are listed in [NOVELTY.md](../../NOVELTY.md). Doubling (Seberry) gives skew-Hadamard matrices of
orders 1640 (already known), 3280, …

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
