# Order 3732

3732 = 4 · 933, with 933 = 3 · 311.

The matrix is a **skew**-Hadamard matrix: H Hᵀ = 3732 I and H + Hᵀ = 2I.

## Construction

**Difference family.** Four ±1 sequences x₁, …, x₄ on ℤ₉₃₃ satisfy two conditions:

- each is constant on the orbits of the multiplier subgroup H = ⟨7⟩ ≤ ℤ₉₃₃*, of order 31:
  {1, 7, 13, 49, 91, 121, 169, 250, 265, 331, 343, 358, 394, 400, 424, 451, 457, 535, 571, 637, 640, 646, 727, 748,
  790, 817, 847, 856, 865, 892, 922};
- their periodic autocorrelations sum to zero at every nonzero shift.

`certificate.json` lists, for each sequence, one representative r of each orbit rH on which the sequence is −1.

**Skew first block.** x₁(0) = 1 and x₁(−t) = −x₁(t) for t ≠ 0. This is possible because −1 ∉ H.

**Row sums.** 1, 61, 3, 1, and 1² + 61² + 3² + 1² = 3732.

**The matrix.** Let A, B, C, D be the circulant matrices with these first rows (`sequences.txt`) and R the
back-diagonal identity. The Goethals–Seidel array

```
    [  A     BR     CR     DR  ]
H = [ -BR    A     DᵀR   -CᵀR ]
    [ -CR   -DᵀR    A     BᵀR ]
    [ -DR    CᵀR  -BᵀR    A   ]
```

is a Hadamard matrix of order 3732, and H + Hᵀ = 2I because A + Aᵀ = 2I.

## The structure: a cyclotomic product of index 10

ℤ₉₃₃ ≅ ℤ₃ × GF(311), and H = {h : h ≡ 1 (mod 3), h mod 311 a 10th power}. Let C₀, …, C₉ be the cyclotomic classes
of index 10 in GF(311)* (each of size 31; −1 lies in C₅). An H-invariant sequence is given, at each a ∈ ℤ₃, by its
value f(a) at (a, 0) and its values x₀(a), …, x₉(a) on (a, C₀), …, (a, C₉): 33 orbits. For index 2 this is the
Paley product used for the skew orders 1012, 1860, 3612 and 6004; here the index is 10.

**Compression identity.** The compression of a block to ℤ₃ is y(a) = f(a) + 31 s(a), with s(a) = Σⱼ xⱼ(a). The
compressions of a Goethals–Seidel quadruple satisfy Σₖ Σₐ yₖ(a)² = 4 · 933, which is equivalent to

  Σ over the 12 positions of (31 s² + 2 f s) = 120.

A position with s ≠ 0 costs at least 31 · 4 − 4 = 120, so **exactly one** position carries s = ±2, with f = −sign(s)
and y = ±61; everywhere else s = 0 and y = f = ±1. The ℤ₃-condition Σₖ PAF_{yₖ}(1) = 0 then leaves one quadruple of
compressions up to symmetry: (1, 1, −1) for the skew block, (61, 1, −1), (1, 1, 1) and (1, 1, −1). This is why the
row sums are (1, 61, 3, 1).

**Three levels.** With the compressions fixed, each block is a choice of class patterns at its three positions:
8,064 choices for the skew block, 13,335,840 for the block with the spike and 16,003,008 for each of the other two.
The remaining conditions split into

- level B: the compression to GF(311) (summing over ℤ₃), whose autocorrelations on the classes C₀, …, C₄ must sum to
  zero (5 equations);
- level C: the autocorrelations at the shifts (1, Cⱼ), j = 0, …, 9 (10 equations).

Together with the ℤ₃-level these are equivalent to the full condition.

## Search

`attack/cyclo_product.py` lists the four blocks, groups each list by its level-B vector (512, 55,915, 57,892 and
57,917 groups), and `attack/join3.c` joins them exactly: group quadruples whose level-B vectors sum to zero are found
by hashing, and inside each one the members are joined on the level-C vectors. Every hit is rebuilt on ℤ₉₃₃ and its
autocorrelations are checked before it is recorded.

The search is exact, but it was stopped early. After 6,992 of the roughly 6 · 10⁷ matching group quadruples it had
found 24 solutions. Extrapolated, the structure holds of the order of 10⁵ skew solutions. A multivariate Gaussian
model over the 14 independent equations, conditioned on the forced compressions, predicts 10^4.6. The lottery ranking
of this subgroup in `attack/NOTES.md` §5 (log₁₀ E = 0.3, about 10¹⁷ pair operations) was too pessimistic by about five
orders of magnitude, in solutions and in work.

## Reproduction

| | |
|---|---|
| Code | the commit that adds this directory (`attack/cyclo_product.py`, `attack/join3.c`) |
| Command | `python3 attack/cyclo_product.py 933 --q 311 --e 10 --config 1 --K 3 --threads 3` |
| Randomness and caps | none in the search (exact enumeration); which solutions appear first depends on the OpenMP schedule; `--K 3` stops after three printed solutions |
| Result of that run | 24 solutions in 6,992 group quadruples, 267 s on 3 threads (lists 55 s, level vectors 100 s, join 110 s) |

The matrix in this directory comes from the first solution of an earlier run of the same search (a scratchpad
prototype of the same code). The tools write hit logs, temporary files and compiled helpers to `$HADAMARD_WORK`
(default `<tmp>/hadamard-work`). The certificate does not depend on the search: `verify/verify.py` rebuilds and checks
the matrix from `certificate.json` alone.

## Status

To our knowledge no Hadamard matrix of order 3732 was known before this one. 3732 is among the unresolved orders in
[FRONTIER.md](../FRONTIER.md) (Cati–Pasechnik Table 4 together with the 2026 constructions), and 933 is undecided in
the Đoković–Seberry census of skew orders 4v (May 2025). A full literature audit is in progress; this section will
carry its checklist.

## Verification

```
python3 verify/verify.py 3732
gzip -dc 3732/H3732.csv.gz | ./check_matrix
```

SHA-256 of `H3732.csv`:
`fd0bf96e4fe3064d96326a339f7983344983f8bf08c55c8c888a10bc41142c2c`

## Files

| File | Contents |
|---|---|
| `certificate.json` | group, multiplier, orbit representatives of the −1 entries, row sums, SHA-256 |
| `sequences.txt` | x₁, …, x₄ as strings of + and − |
| `H3732.csv.gz` | the matrix |
