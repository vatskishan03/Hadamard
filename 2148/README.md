# Order 2148

2148 = 4 · 537, with 537 = 3 · 179.

## Construction

**Bordered Goethals–Seidel array with a border of size 12** (`attack/NOTES.md` §18). Let G = ℤ₅₃₄ ≅ ℤ₃ × ℤ₁₇₈
(the certificate uses the coordinates of ℤ₃ × ℤ₁₇₈) and K = 0 × ℤ₁₇₈, the subgroup of index 3. Four ±1 sequences
x₁, …, x₄ on G satisfy

    Σᵢ PAFᵢ(g) = −12  for g ∈ K \ {0},      Σᵢ PAFᵢ(g) = 0  for g ∉ K.

Their row sums are 2, −2, 2, 0. The Goethals–Seidel array developed over G (order 4 · 534 = 2136) is completed by
12 border rows and columns that are constant on the cosets of K:

    H = [ X  R ]      R[x, (i, g)] = R̂[x, i, g + K],    C[(i, g), x] = Ĉ[i, g + K, x],
        [ C  M ]

with X a Hadamard matrix of order 12 and R̂, Ĉ ±1 matrices of sizes 12 × 12 that depend only on the compressions of
the sequences to G/K ≅ ℤ₃. This is the case Q = ℤ₃, T = all characters of Q of the construction in §18: the border
absorbs the three characters of G that are trivial on K. `verify/verify.py` (`build_bordered_gs`) rebuilds the
matrix from `certificate.json`.

**Symmetry.** Over ℤ₅₃₄ ≅ ℤ₆ × GF(89), the sequences are invariant under the multipliers h ≡ 1 (mod 6) with h mod 89
an eighth power, a subgroup of order 11 of ℤ₅₃₄*: {1, 67, 91, 97, 121, 217, 223, 271, 283, 331, 523}. On each fibre
{a} × GF(89) a sequence is its value at (a, 0) and one value on each of the eight cyclotomic classes of index 8
(54 free signs per sequence). In the certificate's coordinates this is the symmetry (c₀, c₁) ↦ (c₀, 153 c₁).

## Search

1. **Stage 1: the compressions to ℤ₆.** With F = 11 (class size) the compression of a sequence to ℤ₆ is
   y(a) = f(a) + 11 s(a), where f(a) is the value at (a, 0) and s(a) the sum of the eight class values. The bordered
   conditions on the characters of ℤ₆ force exactly two "spikes" (s = ±2), both in one sequence at the positions a
   and a + 3. One class remains up to symmetry: y = (1, 1, −1, 1, −1, 1), (1, −1, −1, −1, 1, −1),
   (1, 1, 1, 1, −1, −1), (−23, 1, 1, 23, −1, −1) (`attack/bgs.py`). This leaves about 1.2 · 10¹¹ sequences per
   block, too many for the joins used for the orders 3732 and 5020.
2. **The ℤ₂-split.** Pair the positions a and a + 3 of ℤ₆: u = (x(a) + x(a + 3))/2 and w = (x(a) − x(a + 3))/2 are
   0/±1 functions on ℤ₃ × GF(89) with complementary supports. The conditions split into one on the u's
   (Σ PAF = −6 on 0 × GF(89) \ {0}, 0 elsewhere; 14 classes, 12 independent) and one on the w's (Σ PAF = 0 off 0;
   11 independent classes), coupled only through the supports. Given the u's, the w's of each block are the sign
   choices on the zero cells of u with fixed sums, about 200 per block. The support condition forces exactly 48 of
   the 96 class cells to be nonzero in the u's; the search takes 12 in each block.
3. **U-solutions by a four-list k-tree** (`attack/kt4.c`). Each list holds 10⁷ random u-vectors of one block (out of
   1.55 · 10⁸ with 12 nonzero class cells), keyed by the 14 class values. Merging on five classes at two split
   values gave 284,292 U-solutions in 17 s.
4. **W-completion.** For each U-solution the four w-lists are joined exactly on the w-conditions
   (`UW.wjoin_fast` in `attack/bgs_uw.py`). One of the 284,292 U-solutions completed. The rate, 3.5 · 10⁻⁶, matches
   a Gaussian estimate of 10^−5.4 per U-solution (with the factor 2¹¹ of the parity lattice of the w-classes).
5. **Check.** The completed sequences satisfy the bordered conditions on ℤ₅₃₄; `attack/border3.c` finds the border for
   their compressions to ℤ₃, and the assembled matrix satisfies H Hᵀ = 2148 I.

The same pipeline, run on the analogue n = 177 (G = ℤ₆ × GF(29), e = 4), finds Hadamard matrices of
order 708 in seconds (171 from 2,000 U-solutions in 4 s); there the structure holds 26,095,424 solutions in its single stage-1 class (counted exactly
with `attack/join3.c`).

## Reproduction

| | |
|---|---|
| Code | `attack/bgs_uw.py`, `attack/kt4.c` (this commit) |
| Command | `python3 attack/bgs_uw.py search 89 8 --L 10000000 --sigmas 2 --maxu 300000 --threads 4 --seed 2` |
| Randomness and caps | `--seed 2` (numpy `default_rng`) draws the four lists of 10⁷ u-vectors; the k-tree uses the two split values nearest the centre; all U-solutions found (below the cap `--maxu`) are completed, so the set of completions does not depend on the thread schedule |
| Result of that run | 284,292 U-solutions, one completion, verified; 197 s on 4 shared cores (2026-10-13) |

The tools write hit logs, temporary files and compiled helpers to `$HADAMARD_WORK` (default `<tmp>/hadamard-work`).
The certificate does not depend on the search.

## Status

Novelty audit in progress. Before this construction the order was unresolved as far as we know: 537(3) is
unresolved in Cati–Pasechnik Table 4 and not covered by Đoković's updates, and after 2060 and 2092 it was the
smallest order we knew to be unresolved.

## Files

| File | Contents |
|---|---|
| `certificate.json` | the four sequences (over ℤ₃ × ℤ₁₇₈), the symmetry, the targets and the border X, R̂, Ĉ |
| `sequences.txt` | the four sequences, one per line |
| `H2148.csv.gz` | the matrix |

## Verification

```
python3 verify/verify.py 2148
gzip -dc 2148/H2148.csv.gz | ./check_matrix
```

SHA-256 of `H2148.csv`:
`c52fc913053fb8f3ef9f9a7cec1b22500eccc00ed000957166724026d8b940f8`
