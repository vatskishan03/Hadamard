# Skew-Hadamard matrix of order 1276

1276 = 4 · 319, with 319 = 11 · 29.

## Construction

A skew-Hadamard matrix is a Hadamard matrix H with H + Hᵀ = 2I.

**Difference family.** Four ±1 sequences x₁, …, x₄ on ℤ₃₁₉ satisfy two conditions:

- each is constant on the orbits of the multiplier subgroup H = ⟨45⟩ = {1, 23, 45, 78, 111, 199, 210} ≤ ℤ₃₁₉*, of
  order 7: H = {h : h ≡ 1 (mod 11), h mod 29 a fourth power};
- their periodic autocorrelations sum to zero at every nonzero shift.

`certificate.json` lists, for each sequence, one representative r of each orbit rH on which the sequence is −1.

**Skew first block.** The first sequence is of skew type: x₁(0) = 1 and x₁(−t) = −x₁(t) for t ≠ 0. This is
possible because −1 ∉ H (−1 is a square but not a fourth power mod 29), so negation pairs the orbits.

**Row sums.** The row sums are 1, 25, 19, 17, and 1² + 25² + 19² + 17² = 1276.

**The matrix.** Let A, B, C, D be the circulant matrices with these first rows (`sequences.txt`) and R the
back-diagonal identity. The Goethals–Seidel array

```
    [  A     BR     CR     DR  ]
H = [ -BR    A     DᵀR   -CᵀR ]
    [ -CR   -DᵀR    A     BᵀR ]
    [ -DR    CᵀR  -BᵀR    A   ]
```

is then a Hadamard matrix of order 1276. Because A + Aᵀ = 2I, it satisfies H + Hᵀ = 2I (Wallis–Whiteman).

## Search

1. **Structure: a cyclotomic product** ([search/NOTES.md](../../search/NOTES.md) §16). Over ℤ₁₁ × GF(29), H acts on
   GF(29) as the classes of index e = 4 (class size F = 7) and trivially on ℤ₁₁. A block is given at each position
   a ∈ ℤ₁₁ by its value f(a) at (a, 0) and its values on the four classes (a, C_j): 55 H-orbits per block.
2. **Stage 1: the compressions to ℤ₁₁.** The compression is y(a) = f(a) + 7 s(a), and Σ_k Σ_a y_k(a)² = 4 · 319
   becomes Σ over the 44 positions of (7 s² + 2 f s) = 160. Stage 1 is sampled by `cp_lowe.stage1_sample` (seed 1);
   its first quadruple is y = (1, −1, −1, 1, −1, −1, 1, 1, −1, 1, 1), (−1, −1, −1, 1, 1, 29, −1, −1, −1, 1, −1),
   (1, 1, 1, 1, −1, 1, 1, 1, 13, 1, −1), (1, 1, 1, 1, 1, −1, 15, −1, −1, 1, −1).
3. **Lists.** The blocks of this class have 31,104 (skew block), 60,466,176, 241,864,704 and 241,864,704 members.
   The Gaussian estimate of the class (22 coordinates, rank 21) is 10^3.19 solutions.
4. **Low-energy k-tree** (`search/cp_lowe.py`, `search/cp_enum.c`, `search/kt4.c`;
   [search/NOTES.md](../../search/NOTES.md) §20–22). The skew block and the second block are kept whole; `cp_enum.c`
   keeps the 9,966,434 and 10,012,823 members of the other two blocks whose autocorrelation vectors are closest to
   their conditional means given that the four vectors sum to zero. The Gaussian estimate inside these lists is
   10^1.64 solutions. `kt4.c` joins the lists over split values of the merge classes 0, 1, 3, 4, 5, ranked by expected
   solutions per unit of work (6000 values holding a predicted share of 0.09 of the solutions in the lists).
5. **The run.** A first run (old join) searched 267 split values without a solution and was stopped to switch to the
   radix-partitioned join. The second run rebuilt the lists and skipped the 131 values it had in common with the first
   (`--skip-sigmas`); the solution came at its 321st split value, σ = (−6, 18, 6, −2, −6), after 3203 s of joins on 2
   threads (about 69 minutes including the enumeration).

**Reproduction.**

| | |
|---|---|
| Code | commit `45afdab` (`search/cp_lowe.py`, `search/cp_enum.c`, `search/kt4.c`) |
| Command | `python3 search/cp_lowe.py 319 --q 29 --e 4 --sample-stage1 3 --config 0 --L 6.5e7,1e7,1e7 --nsig 6000 --threads 2 --K 3 --stream --samples 1000000 --S 0,1,3,4,5` (the run that found it also passed `--skip-sigmas` with the 267 split values of the first run) |
| Randomness and caps | `--seed 1` (default) draws the sampled stage-1 quadruples, the 10⁶ samples per block for the energy statistics and the 4·10⁶ pairs per side that rank the split values; `--L 6.5e7,1e7,1e7` members kept per non-skew block (the first non-skew block has fewer, so it is kept whole) |
| Determinism | the lists are the same sets on every run, but their order follows the threaded enumeration, and the samples that rank the split values are drawn by position; the ranking can therefore differ between runs (131 of the first 267 values recurred) |
| Result | one verified solution at the 321st split value; the certificate is that hit |

## Status

**New** as a skew-Hadamard matrix. An ordinary Hadamard matrix of order 1276 is known (Cati–Pasechnik Table 4 has no
entry for 319; Piza's table gives it by Cooper–Wallis from orders 11 and 29); the skew property is what is new.

- The Đoković–Seberry census lists 319 among its 737 undecided v. 319 is in none of its Paley, Szekeres,
  Szekeres–Whiteman or Spence series, nor among its 41 sporadic values.
- Cati–Pasechnik Table 5 gives 319(3): the smallest known skew order with odd part 319 was 8 · 319 = 2552, from
  Paley I (2551 is a prime ≡ 3 (mod 4)). Their Table 3 (skew constructions) ends at n = 299. They state that the first
  skew order missing from Sage is 1252.
- Seberry–Yamada (1992, Table 7.1, p. 499) gives "319 3". The papers of Đoković and his coauthors contain neither 1276
  nor 2552.
- SageMath's source contains neither 1276 nor 319, and every branch of `skew_hadamard_matrix(1276, existence=True)`
  fails, so it returns Unknown. Its list of unknown skew orders stops at 1196.
- No classical construction gives skew 1276:
  - Paley I: 1275 = 3 · 5² · 17. Szekeres, Szekeres–Whiteman, Spence's G-matrices and Momihara–Xiang 2(q + 1):
    q = 637 = 7² · 13. Wallis–Whiteman and Momihara–Xiang 4(q + 1): q = 318.
  - Spence 4(q² + q + 1): q² + q + 1 = 319 has no integer root. The product theorem n(pᵏ + 1) would need n = 29 or
    319, which are odd. Amicable products need 16 | 1276.
  - Seberry's mn(n − 1) works only with n = 2, which needs an OD(1276; 1, 638, 637), i.e. a skew matrix of order 1276
    already. Doubling needs 8 | 1276.

Not checked: Koukouvinos–Stylianou (2008, Table 1).

The sources searched are listed in [NOVELTY.md](../../NOVELTY.md). Doubling (Seberry) gives skew-Hadamard matrices of
orders 2552 (already known), 5104, …

## Files

| File | Contents |
|---|---|
| `certificate.json` | the four sequences as orbit representatives, H, the row sums |
| `sequences.txt` | the four sequences, one per line |
| `H1276.csv.gz` | the matrix |

## Verification

```
python3 verify/verify.py skew/1276
gzip -dc skew/1276/H1276.csv.gz | ./check_matrix
```

SHA-256 of `H1276.csv`: `4e33b918c799f84f0eac51e00ec1c46ac9c4b87d75eae73b0cb01d061e0cbd3f`
