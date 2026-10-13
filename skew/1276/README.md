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

1. **Structure: a cyclotomic product** (`attack/NOTES.md` §16). Over ℤ₁₁ × GF(29), H acts on GF(29) as the classes
   of index e = 4 (class size F = 7) and trivially on ℤ₁₁. A block is given at each position a ∈ ℤ₁₁ by its value
   f(a) at (a, 0) and its values on the four classes (a, C_j): 55 H-orbits per block.
2. **Stage 1: the compressions to ℤ₁₁.** The compression is y(a) = f(a) + 7 s(a), and Σ_k Σ_a y_k(a)² = 4 · 319
   becomes Σ over the 44 positions of (7 s² + 2 f s) = 160. Stage 1 is sampled by `cp_lowe.stage1_sample` (seed 1);
   its first quadruple is y = (1, −1, −1, 1, −1, −1, 1, 1, −1, 1, 1), (−1, −1, −1, 1, 1, 29, −1, −1, −1, 1, −1),
   (1, 1, 1, 1, −1, 1, 1, 1, 13, 1, −1), (1, 1, 1, 1, 1, −1, 15, −1, −1, 1, −1).
3. **Lists.** The blocks of this class have 31,104 (skew block), 60,466,176, 241,864,704 and 241,864,704 members.
   The Gaussian estimate of the class (22 coordinates, rank 21) is 10^3.19 solutions.
4. **Low-energy k-tree** (`attack/cp_lowe.py`, `attack/cp_enum.c`, `attack/kt4.c`; `attack/NOTES.md` §20–22).
   The skew block and the second block are kept whole; `cp_enum.c` keeps the 9,966,434 and 10,012,823 members of the
   other two blocks whose autocorrelation vectors are closest to their conditional means given that the four vectors
   sum to zero. The Gaussian estimate inside these lists is 10^1.64 solutions. `kt4.c` joins the lists over split
   values of the merge classes 0, 1, 3, 4, 5, ranked by expected solutions per unit of work (6000 values holding a
   predicted share of 0.09 of the solutions in the lists).
5. **The run.** A first run (old join) searched 267 split values without a solution and was stopped to switch to the
   radix-partitioned join. The second run rebuilt the lists and skipped the 131 values it had in common with the first
   (`--skip-sigmas`); the solution came at its 321st split value, σ = (−6, 18, 6, −2, −6), after 3203 s of joins on 2
   threads (about 69 minutes including the enumeration).

## Reproduction

| | |
|---|---|
| Code | commit `c29fd15` (`attack/cp_lowe.py`, `attack/cp_enum.c`, `attack/kt4.c`) |
| Command | `python3 attack/cp_lowe.py 319 --q 29 --e 4 --sample-stage1 3 --config 0 --L 6.5e7,1e7,1e7 --nsig 6000 --threads 2 --K 3 --stream --samples 1000000 --S 0,1,3,4,5` (the run that found it also passed `--skip-sigmas` with the 267 split values of the first run) |
| Randomness and caps | `--seed 1` (default) draws the sampled stage-1 quadruples, the 10⁶ samples per block for the energy statistics and the 4·10⁶ pairs per side that rank the split values; `--L 6.5e7,1e7,1e7` members kept per non-skew block (the first non-skew block has fewer, so it is kept whole) |
| Determinism | the lists are the same sets on every run, but their order follows the threaded enumeration, and the samples that rank the split values are drawn by position; the ranking can therefore differ between runs (131 of the first 267 values recurred) |
| Result of that run | one verified solution at the 321st split value; the certificate is that hit |

The tools write hit logs, temporary files and compiled helpers to `$HADAMARD_WORK` (default `<tmp>/hadamard-work`).
The certificate does not depend on the search: `verify/verify.py` rebuilds and checks the matrix from
`certificate.json` alone.

## Status

To our knowledge this is the first skew-Hadamard matrix of order 1276.

- The Đoković–Seberry census of skew orders 4v for odd v < 2500 (PDF of 2025-05-28) lists 319 among its undecided
  values.
- Cati–Pasechnik Table 5 (arXiv:2411.18897) gives 319(3): the smallest known skew order with odd part 319 was
  8 · 319 = 2552.
- An *ordinary* Hadamard matrix of order 1276 is known (319 is not among the unresolved n of Cati–Pasechnik Table 4).
  Only the skew property is new.

**Novelty checklist:** pending (a literature audit in the format of [skew/1204](../1204/README.md) is being run).

Doubling (Seberry) gives skew orders 2552 (already known, see above), 5104, …

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
