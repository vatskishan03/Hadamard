# Hadamard matrix of order 3732

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

### The structure: a cyclotomic product of index 10

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

`search/cyclo_product.py` lists the four blocks, groups each list by its level-B vector (512, 55,915, 57,892 and
57,917 groups), and `search/join3.c` joins them exactly: group quadruples whose level-B vectors sum to zero are found
by hashing, and inside each one the members are joined on the level-C vectors. Every hit is rebuilt on ℤ₉₃₃ and its
autocorrelations are checked before it is recorded.

The search is exact, but it was stopped early. After 6,978 of the roughly 6 · 10⁷ matching group quadruples it had
found 24 solutions. Extrapolated, the structure holds of the order of 10⁵ skew solutions. A multivariate Gaussian
model over the 14 independent equations, conditioned on the forced compressions, predicts 10^4.6. The lottery ranking
of this subgroup in [search/NOTES.md](../../search/NOTES.md) §5 (log₁₀ E = 0.3, about 10¹⁷ pair operations) was too
pessimistic by about five orders of magnitude, in solutions and in work.

**Reproduction.**

| | |
|---|---|
| Code | commit `3d9645a` (`search/cyclo_product.py`, `search/join3.c`) |
| Command | `python3 search/cyclo_product.py 933 --q 311 --e 10 --config 0 --K 3 --threads 3` |
| Randomness and caps | none in the search (exact enumeration); which solutions appear first depends on the OpenMP schedule; `--K 3` stops after three printed solutions |
| Result | 24 solutions in 6,978 group quadruples, 306 s on 3 threads (lists 54 s, level vectors 125 s, join 127 s) |

The matrix in this directory comes from the first solution of an earlier run of the same search (a prototype of the
same code).

## Status

**Concurrent** as a Hadamard matrix, **New** as a skew-Hadamard matrix. Cati–Pasechnik Table 4 lists 933(4): the
smallest known Hadamard order of the form 2ᵐ · 933 was 14928, so 3732 was unknown. Two other groups constructed
Hadamard matrices of order 3732 independently: Egan, Kotsireas, London and Ó Catháin (arXiv:2610.27175,
19 October 2026), by a Goethals–Seidel array with a multiplier group of order 31 as here, and Cooley's repository
[Middletownbooks/Hadamard](https://github.com/Middletownbooks/Hadamard), committed on 14 October 2026
(Goethals–Seidel over ℤ₃ × ℤ₃₁₁ with the cyclotomic classes of index 10). Ours was committed on 13 October 2026;
we claim no priority for the order. The matrix here is skew, and neither of the others is stated to be skew; no
skew-Hadamard matrix of order 3732 was known.

- Seberry–Yamada (1992, Table A.2: "933 4 d1", Agaian multiplication) and the Handbook of Combinatorial Designs (2007,
  Table 1.53: 933 → 4 and 311 → 3) give the same entry. Seberry–Yamada Table A.1 has no Williamson or Williamson-type
  matrices of order 311 or 933.
- Cati–Pasechnik mark no update at this entry, and 933 is not among the entries dropped from Seberry–Yamada (2020),
  Table A.17.
- For skew matrices, Cati–Pasechnik Table 5 lists 933(10). They fill this entry in via Paley's construction, where
  Seberry–Yamada (2020) lists m as unknown: 955391 is prime, and 2ᵐ · 933 − 1 is composite for m = 2, …, 9. The
  smallest known skew order of the form 2ᵐ · 933 was 955392. The Đoković–Seberry census (May 2025) lists 933 among its
  737 undecided values, and Seberry–Yamada (1992, Table 7.1) has 933 as unknown.
- SageMath returns Unknown for 3732, as a Hadamard matrix and as a skew-Hadamard matrix.
- The 2026 matrix of order 1244 over ℤ₃₁₁ is a plain Goethals–Seidel matrix, not Williamson-type, so no Cooper–Wallis
  or tripling route takes it to 3732.
- No recursive construction reaches 3732:
  - Paley I and II: 3731 = 7 · 13 · 41 and 1865 = 5 · 373 are not prime powers.
  - Miyamoto: 933 is not a prime power, and a symmetric conference matrix of order 934 cannot exist, since 933 is not
    a sum of two squares.
  - Yamada and Kiyasu (q = 931 = 7² · 19), Mathon, Turyn, Whiteman and Spence.
  - T-matrices and Williamson-type matrices of order 311 or 933 are unknown.
  - The Wallis–Whiteman and Spence skew families.
  - Product theorems: 3732 ≡ 4 (mod 8).

**Consequences.**

- Doubling gives skew-Hadamard matrices of all orders 2ᵐ · 933 with m ≥ 2, of which m = 2, …, 9 were unknown.
- 3733 is a prime ≡ 1 (mod 4), so Miyamoto's theorem gives a Hadamard matrix of order 4 · 3733 = 14932; see
  [14932](../14932).

The sources searched are listed in [NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | group, multiplier, orbit representatives of the −1 entries, row sums, SHA-256 |
| `sequences.txt` | x₁, …, x₄ as strings of + and − |
| `H3732.csv.gz` | the matrix |

## Verification

```
python3 verify/verify.py 3732
gzip -dc hadamard/3732/H3732.csv.gz | ./check_matrix
```

SHA-256 of `H3732.csv`: `fd0bf96e4fe3064d96326a339f7983344983f8bf08c55c8c888a10bc41142c2c`
