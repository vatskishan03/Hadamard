# Skew-Hadamard matrix of order 1540

1540 = 4 · 385, with 385 = 5 · 7 · 11.

## Construction

A skew-Hadamard matrix is a Hadamard matrix H with H + Hᵀ = 2I.

**Difference family.** Four ±1 sequences x₁, …, x₄ on ℤ₃₈₅ satisfy two conditions:

- each is constant on the orbits of the multiplier subgroup H = ⟨246, 221⟩ ≤ ℤ₃₈₅*, of order 15. In ℤ₃₈₅ = ℤ₃₅ × ℤ₁₁,
  246 = (1, 4) and 221 = (11, 1), so H is ⟨11⟩ = {1, 11, 16} on ℤ₃₅ times the squares of F₁₁*;
- their periodic autocorrelations sum to zero at every nonzero shift.

`certificate.json` lists, for each sequence, the least element r of each orbit rH on which the sequence is −1.

**Skew first block.** The first sequence is of skew type: x₁(0) = 1 and x₁(−t) = −x₁(t) for t ≠ 0. This is possible
because −1 ∉ H (H has odd order), so negation pairs the orbits.

**Row sums.** The row sums are 1, −7, 31, 23, and 1² + 7² + 31² + 23² = 1540.

**The matrix.** Let A, B, C, D be the circulant matrices with these first rows (`sequences.txt`) and R the
back-diagonal identity. The Goethals–Seidel array

```
    [  A     BR     CR     DR  ]
H = [ -BR    A     DᵀR   -CᵀR ]
    [ -CR   -DᵀR    A     BᵀR ]
    [ -DR    CᵀR  -BᵀR    A   ]
```

gives a Hadamard matrix of order 1540, and H + Hᵀ = 2I because A + Aᵀ = 2I.

## Search

1. **Structure: a Paley product with a multiplier on K** ([search/NOTES.md](../../search/NOTES.md) §15, §87–88,
   §91). Over K × F₁₁ with K = ℤ₃₅, each block is u ⊗ 1 + w ⊗ J + z ⊗ χ with χ the Legendre sequence of F₁₁; ⟨11⟩ has
   15 orbits on ℤ₃₅ (the five multiples of 7 are fixed), so a block has 45 orbits. The counting identity of §15 reads
   12 N_w − 4 N₋ = 140.
2. **The solution.** The w-positions are {5, 10, 15, 20, 25, 30} in the first block, 0 in the second, 21 in the third
   and 0, 15, 25, 30 in the fourth (N_w = 12), the position 0 of the fourth of the second type (N₋ = 1): 144 − 4 = 140.
   The compressions to ℤ₃₅ take the values ±1, ±9 and ±11.
3. **Engine** (`search/skewplug`, `exh38.py`). The sampler `prod_dedup.py` met each class of level-1 instances of this
   shape about 2300 times, so it found few new ones. `exh38.py` instead generates, for each block, defect weight and
   row sum, the list of blocks until it covers that block space. It keeps one block per orbit of the exact block
   symmetries (reflection, negation, translation by the points of ℤ₃₅ fixed by ⟨11⟩) and reduces the skew block by the
   units of ℤ₃₅. It then joins every row-sum pattern of a slot exactly (`gl1.c`) and sends the new classes to the exact
   level 2. It runs about 100 classes per second against 2 to 3 for the sampler. The slot (6, 1, 1, 4) / (0, 0, 0, 1)
   is the 7th of the 2747 row-sum-feasible slots of the shape in the order of the estimated completions per class
   (logE −1.54, logρ −7.49; heuristic), below the reachability threshold of §90. 2 of its 100 classes complete; the
   second gives a different matrix with the same row sums.

**Reproduction.**

| | |
|---|---|
| Code | `search/skewplug` (`exh38.py`; build commands in its README) |
| Command | copy `data/repro385/est_hit_slot.jsonl` to the directory, then `PROD_CONTINUE=1 python3 -u exh38.py hit385 600 "[35]" 11 "[[11]]" est_hit_slot.jsonl 1 file skf` and `python3 cyclic.py out/skew_1540_hit385_s0_<time>_<i>.json OUT.json` |
| Randomness and caps | none: the block lists come from seeds fixed by the block data, and the join and level 2 are deterministic |
| Result | two `HIT` lines with the group-form SHA-256 `9ea5b315…` and `1b4d57ec…` (over ℤ₃₅ × ℤ₁₁), in about 14 s; `cyclic.py` turns the first into `certificate.json` without its `sha256` field |

## Status

**New** as a skew-Hadamard matrix; it settles v = 385 in the Đoković–Seberry census. An ordinary Hadamard matrix of
order 1540 is known (385 is not in Cati–Pasechnik Table 4); the skew property is what is new.

- The Đoković–Seberry census lists 385 among its undecided v.
- Cati–Pasechnik Table 5 gives 385(3): the smallest known skew order with odd part 385 was 8 · 385 = 3080.
- SageMath's `skew_hadamard_matrix` (development source of 21 October 2026, read, not run) has no construction of
  order 1540, and 385 is in none of its difference-set tables.
- We found no classical construction:
  - Paley I needs 1539 = 3⁴ · 19.
  - The two-block families of order 2(q + 1) need q = 769 ≡ 5 (mod 8) (Szekeres), a power q = pᵗ with p ≡ 5 (mod 8)
    (Whiteman), and Momihara–Xiang's difference families give the order 2ᵘ⁻¹(qᵉ + 1) with q ≡ 2ᵘ + 1 (mod 2ᵘ⁺¹),
    which for 1540 would need u = 2, so qᵉ = 769 ≡ 5 (mod 8), or u = 3 and qᵉ = 384. 769 ≡ 1 (mod 8).
  - Wallis–Whiteman 4(q + 1) needs 384; Spence 4(q² + q + 1) needs q² + q + 1 = 385; Turyn (p − 1)ᵘ + 1 needs a
    perfect power 1539; the family 3(q + 1) of [PALEY_SKEW.md](../PALEY_SKEW.md) needs 3 | 1540.
  - Skew × amicable products: in 1540 = 2 · 770 = 4 · 385 = 20 · 77 = 28 · 55 = 44 · 35 = 140 · 11 one factor is not a
    Hadamard order.
  - Good matrices and skew Williamson-type matrices of order 385 are not recorded.

The sources searched are listed in [NOVELTY.md](../../NOVELTY.md). Doubling gives skew-Hadamard matrices of orders
3080 (Table 5's value), 6160, …

## Files

| File | Contents |
|---|---|
| `certificate.json` | the multipliers, the orbit representatives of the −1 positions of each sequence, the row sums |
| `sequences.txt` | the four sequences, one per line |
| `H1540.csv.gz` | the matrix |

## Verification

```
python3 verify/verify.py skew/1540
gzip -dc skew/1540/H1540.csv.gz | ./check_matrix
```

SHA-256 of `H1540.csv`: `f98a3775459959ba41f39b69220e52ed166feb71020107d8d11c2e08eaab1b5d`
