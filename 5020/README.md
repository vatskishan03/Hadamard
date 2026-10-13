# Order 5020

5020 = 4 · 1255, with 1255 = 5 · 251.

## Construction

**Difference family.** Four ±1 sequences x₁, …, x₄ on ℤ₁₂₅₅ satisfy two conditions:

- each is constant on the orbits of the multiplier subgroup H = ⟨276, 754⟩ ≤ ℤ₁₂₅₅*, of order 50:
  {1, 4, 16, 51, 64, 69, 91, 94, 149, 201, 204, 211, 219, 241, 249, 256, 271, 276, 314, 331, 351, 364, 374, 376, 494,
  506, 566, 571, 596, 651, 706, 721, 751, 754, 769, 804, 816, 844, 866, 876, 954, 964, 994, 996, 1009, 1024, 1029, 1084,
  1104, 1129};
- their periodic autocorrelations sum to zero at every nonzero shift.

`certificate.json` lists, for each sequence, one representative r of each orbit rH on which the sequence is −1.

**Row sums.** 1, −1, −53, −47, and 1² + 1² + 53² + 47² = 5020.

**The matrix.** With A, B, C, D the circulant matrices with these first rows (`sequences.txt`) and R the
back-diagonal identity, the Goethals–Seidel array

```
    [  A     BR     CR     DR  ]
H = [ -BR    A     DᵀR   -CᵀR ]
    [ -CR   -DᵀR    A     BᵀR ]
    [ -DR    CᵀR  -BᵀR    A   ]
```

is a Hadamard matrix of order 5020.

## The structure: a cyclotomic product with a multiplier on ℤ₅

ℤ₁₂₅₅ ≅ ℤ₅ × GF(251), and H = {h : h ≡ ±1 (mod 5), h mod 251 a 10th power}: it acts on ℤ₅ by ±1 (orbits {0},
{±1}, {±2}) and on GF(251) by the 10th powers, of order 25 (the cyclotomic classes C₀, …, C₉ of index 10). An
H-invariant sequence is given, on each of the three ℤ₅-orbits, by its value f at (a, 0) and its values x₀, …, x₉ on
(a, C₀), …, (a, C₉): 33 orbits. This is the structure of [order 3732](../3732) with one more multiplier.

**Compression identity.** The compression to ℤ₅ is y(a) = f(a) + 25 s(a), s(a) = Σⱼ xⱼ(a), and the four compressions
satisfy Σₖ Σₐ yₖ(a)² = 4 · 1255, that is

  Σ over the 20 positions of (25 s² + 2 f s) = 200.

A position with s = ±2 costs 96 (f = −sign s) or 104 (f = sign s), twice that on the orbits {±1}, {±2}. So exactly two
positions carry s ≠ 0: position 0 of two different blocks, with y = −49 and −51. Up to symmetry this leaves one
quadruple of compressions, with row sums (1, −1, −53, −47).

## Search

`attack/cyclo_product.py` (plain mode, `--kgen 4`) lists the blocks: 16,003,008 members for each of the two blocks
without a spike and 13,335,840 for each spike block. It groups them by the autocorrelations of their compression to
GF(251) (level B, 5 equations), and `attack/join3.c` joins group quadruples exactly on the remaining 20 equations.
The level-B groups number about 10⁶ per list, too many for a complete join, so the run kept the 8,000 largest groups
of each list (7.6 % and 9.0 % of the members). It found 32 solutions after 17,273 group quadruples. All were rebuilt
on ℤ₁₂₅₅ and checked.

## Reproduction

| | |
|---|---|
| Code | commit `81b18b6` (`attack/cyclo_product.py`, `attack/join3.c`) |
| Command | `python3 attack/cyclo_product.py 1255 --q 251 --e 10 --plain --kgen 4 --config 0 --K 3 --threads 3 --sub 8000` |
| Randomness and caps | none in the search; `--sub 8000` keeps the 8,000 largest level-B groups of each list (ties broken by a stable sort); which solutions appear first depends on the OpenMP schedule |
| Result of that run | 32 solutions in 17,273 group quadruples, 349 s on 3 threads (lists 19 s, levels and grouping 280 s, join 48 s) |

The certificate is the first recorded solution (`attack/gs_certify.py 1255`). The tools write hit logs, temporary
files and compiled helpers to `$HADAMARD_WORK` (default `<tmp>/hadamard-work`). `verify/verify.py` rebuilds and
checks the matrix from `certificate.json` alone.

## Status

To our knowledge no Hadamard matrix of order 5020 was known before this one. 5020 is among the unresolved orders in
[FRONTIER.md](../FRONTIER.md) (Cati–Pasechnik Table 4 together with the 2026 constructions). A full literature audit
is in progress; this section will carry its checklist.

5021 is a prime ≡ 1 (mod 4), so Miyamoto's theorem gives a Hadamard matrix of order 4 · 5021 = 20084 from this one.

## Verification

```
python3 verify/verify.py 5020
gzip -dc 5020/H5020.csv.gz | ./check_matrix
```

SHA-256 of `H5020.csv`:
`71d4cb944767f77cf49fae68b9d164bd722eb36786df41a642d7c8a97428e8f7`

## Files

| File | Contents |
|---|---|
| `certificate.json` | group, multipliers, orbit representatives of the −1 entries, row sums, SHA-256 |
| `sequences.txt` | x₁, …, x₄ as strings of + and − |
| `H5020.csv.gz` | the matrix |
