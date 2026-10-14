# Skew-Hadamard matrix of order 1204

1204 = 4 · 301, with 301 = 7 · 43.

## Construction

A skew-Hadamard matrix is a Hadamard matrix H with H + Hᵀ = 2I.

**Difference family.** Four ±1 sequences x₁, …, x₄ on ℤ₃₀₁ satisfy two conditions:

- each is constant on the orbits of the multiplier subgroup H = ⟨127⟩ = {1, 64, 78, 127, 176, 183, 274} ≤ ℤ₃₀₁*, of
  order 7: H = {h : h ≡ 1 (mod 7), h mod 43 a sixth power};
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

1. **Structure: a cyclotomic product** ([search/NOTES.md](../../search/NOTES.md) §16). Over ℤ₇ × GF(43), H acts on
   GF(43) as the classes of index e = 6 (class size F = 7) and trivially on ℤ₇. A block is given at each position
   a ∈ ℤ₇ by its value f(a) at (a, 0) and its values on the six classes (a, C_j): 49 H-orbits per block.
2. **Stage 1: the compressions to ℤ₇.** The compression is y(a) = f(a) + 7 s(a), and Σ_k Σ_a y_k(a)² = 4 · 301 becomes
   Σ over the 28 positions of (7 s² + 2 f s) = 168. Stage 1 is sampled by `cp_lowe.stage1_sample` (seed 1); its first
   quadruple is y = (1, 1, 1, −1, 1, −1, −1), (−13, −1, −1, 1, 1, 1, 1), (−1, 1, −13, 1, −15, 1, 15),
   (−1, −13, −1, −1, −15, −1, 1).
3. **Lists.** The blocks of this class have 64,000 (skew block), 960,000,000, 540,000,000 and 720,000,000 members. The
   Gaussian estimate of the class (21 coordinates, rank 20) is 10^6.03 solutions.
4. **Low-energy k-tree** (`search/cp_lowe.py`, `search/cp_enum.c`, `search/kt4.c`;
   [search/NOTES.md](../../search/NOTES.md) §20). `cp_enum.c` enumerates every member of each non-skew block and keeps
   those whose autocorrelation vector is close to its conditional mean given that the four vectors sum to zero:
   9,912,243, 10,109,883 and 9,927,286 members, with a Gaussian estimate of 10^2.52 solutions inside the lists.
   `kt4.c` joins them over split values of five merge classes ranked by expected solutions per unit of work; the first
   solution came at the 24th split value. The structure was selected by the estimate-then-search loop
   `search/cp_hunt.py` (§20).

**Reproduction.**

| | |
|---|---|
| Code | commit `ed9f756` (`search/cp_lowe.py`, `search/cp_enum.c`, `search/kt4.c`) |
| Command | `python3 search/cp_lowe.py 301 --q 43 --e 6 --kgen 1 --twist 0 --sample-stage1 3 --config 0 --L 1e7 --nsig 3000 --threads 3 --K 3` |
| Randomness and caps | `--seed 1` (default) draws the sampled stage-1 quadruples, the 2·10⁵ samples per block for the energy statistics and the 4·10⁶ pairs per side that rank the split values; `--L 1e7` members kept per non-skew block |
| Result | first verified solution after about 6 minutes on 3 threads (enumeration 121 s, split values 1–24 in 172 s); the certificate is the first hit recorded |

## Status

**New** as a skew-Hadamard matrix. An ordinary Hadamard matrix of order 1204 is classical (Paley II with q = 601; the
census also lists 301 in its symmetric Paley series); the skew property is what is new.

- The Đoković–Seberry census lists 301 among its 737 undecided v. 301 is in none of its Paley, Szekeres,
  Szekeres–Whiteman or Spence series, nor among its 41 sporadic values.
- Cati–Pasechnik Table 5 gives 301(3): the smallest known skew order with odd part 301 was 8 · 301 = 2408, from
  Wallis–Whiteman with q = 601 ≡ 9 (mod 16). Their Table 3 (skew constructions, n ≤ 1200) ends at n = 299. They state
  that the first skew order missing from Sage is 1252, and Sage has no 1204.
- Seberry–Yamada (1992, Table 7.1) gives "301 3" and Koukouvinos–Stylianou (2008, Table 1) gives 301(3). The papers of
  Đoković and of Đoković and Kotsireas contain no skew matrix of order 1204 or 2408.
- In SageMath, `skew_hadamard_matrix(1204, existence=True)` returns Unknown. Its list of unknown skew orders stops at
  1200, so it neither lists nor excludes 1204.
- Piza's deposit lists "1204 = 2(601 + 1)", which is Paley II. arXiv:2404.14375 contains 301 only as a block size.
- No classical construction gives skew 1204:
  - Paley I: 1203 = 3 · 401. Szekeres 2(q + 1) and Spence's G-matrices need q = 601 ≡ 5 (mod 8), but 601 ≡ 1 (mod 8).
  - Wallis–Whiteman 4(q + 1) needs q = 300; q = 601 gives 2408.
  - Products and doubling need 8 | 1204. Spence 4(q² + q + 1), (p − 1)ᵘ + 1, (|t| + 1)(q + 1), 2ᵗ q and Seberry's
    mn(n − 1) all fail by arithmetic.
  - Good matrices are not known for order 301.

The sources searched are listed in [NOVELTY.md](../../NOVELTY.md). Doubling (Seberry) gives skew-Hadamard matrices of
orders 2408 (already known), 4816, …

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
