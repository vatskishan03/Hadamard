# Skew-Hadamard matrix of order 3124

3124 = 4 · 781, with 781 = 11 · 71.

## Construction

A skew-Hadamard matrix is a Hadamard matrix H with H + Hᵀ = 2I.

**Difference family.** Four ±1 sequences x₁, …, x₄ on ℤ₇₈₁ satisfy two conditions:

- each is constant on the orbits of the multiplier subgroup H = ⟨356, 309⟩ ≤ ℤ₇₈₁*, of order 25. In ℤ₇₈₁ = ℤ₇₁ × ℤ₁₁,
  H is ⟨25⟩ (order 5) on ℤ₇₁ times the squares of F₁₁*;
- their periodic autocorrelations sum to zero at every nonzero shift.

`certificate.json` lists, for each sequence, the least element r of each orbit rH on which the sequence is −1.

**Skew first block.** The first sequence is of skew type: x₁(0) = 1 and x₁(−t) = −x₁(t) for t ≠ 0. This is possible
because −1 ∉ H (H has odd order), so negation pairs the orbits.

**Row sums.** The row sums are 1, −51, 9, 21, and 1² + 51² + 9² + 21² = 3124.

**The matrix.** Let A, B, C, D be the circulant matrices with these first rows (`sequences.txt`) and R the
back-diagonal identity. The Goethals–Seidel array

```
    [  A     BR     CR     DR  ]
H = [ -BR    A     DᵀR   -CᵀR ]
    [ -CR   -DᵀR    A     BᵀR ]
    [ -DR    CᵀR  -BᵀR    A   ]
```

gives a Hadamard matrix of order 3124, and H + Hᵀ = 2I because A + Aᵀ = 2I.

## Search

1. **Structure: a Paley product with a multiplier on K** ([search/NOTES.md](../../search/NOTES.md) §15, §87–88,
   §91, §93). Over K × F₁₁ with K = ℤ₇₁, each block is u ⊗ 1 + w ⊗ J + z ⊗ χ with χ the Legendre sequence of F₁₁;
   ⟨25⟩ has 15 orbits on ℤ₇₁ (0 and 14 orbits of size 5), so a block has 45 orbits. The counting identity of §15 reads
   12 N_w − 4 N₋ = 284.
2. **The solution.** The first block has no w-position; the second has 0, 1, 5, 25, 54, 57; the third 1, 5, 7, 23,
   25, 33, 35, 44, 54, 57 and the fourth 0, 6, 8, 11, 26, 30, 40, 55, 58, 59, 62, five of each of the last two of the
   second type (N_w = 27, N₋ = 10): 324 − 40 = 284. The compressions to ℤ₇₁ take the values ±1, ±9 and ±11. The slot
   (0, 6, 10, 11) / (0, 0, 5, 5) has the lowest row-sum-feasible energy of the shape, E = 27.
3. **Engine** (`search/skewplug`, `exh38.py`, NOTES §91). The sweep of round 40 enumerated the lowest-energy slots of
   every open real skew census shape. This slot is heavy (block specifications with up to 6 · 10⁶ raw
   configurations): its 1841 classes of level-1 instances give one completion.

**Reproduction.**

| | |
|---|---|
| Code | `search/skewplug` (`exh38.py`, `cyclic.py`; build commands in its README) |
| Command | copy `data/repro781/est_hit781.jsonl` to the directory, then `EX_COV=0.99 EX_PERSIST_MB=0 EX_MAXDRAW=6000000 EX_JOINCAP=2e9 EX_SLOT_SECS=600 PROD_CONTINUE=1 python3 -u exh38.py hit781 1500 "[71]" 11 "[[25]]" est_hit781.jsonl 1 file sl2` and `python3 cyclic.py out/skew_3124_hit781_s0_<time>_<i>.json OUT.json` |
| Randomness and caps | none: the block lists come from seeds fixed by the block data (up to 6 · 10⁶ draws per specification), and the join and level 2 are deterministic |
| Result | one `HIT` line, group-form SHA-256 `0e0e2e2f…` (over ℤ₇₁ × ℤ₁₁), 1841 classes, in 1.5 to 3 minutes; `cyclic.py` turns it into `certificate.json` without its `sha256` field |

## Status

**New** as a skew-Hadamard matrix; it settles v = 781 in the Đoković–Seberry census. An ordinary Hadamard matrix of
order 3124 is known (781 is not in Cati–Pasechnik Table 4); the skew property is what is new.

- The Đoković–Seberry census lists 781 among its undecided v.
- Cati–Pasechnik Table 5 gives 781(3): the smallest known skew order with odd part 781 was 8 · 781 = 6248.
- We found no classical construction:
  - Paley I needs 3123 = 3² · 347.
  - The two-block families of order 2(q + 1) need q = 1561 = 7 · 223; Momihara–Xiang's order 2ᵘ⁻¹(qᵉ + 1) would need
    qᵉ = 1561 (u = 2) or 780 (u = 3).
  - Wallis–Whiteman 4(q + 1) needs 780; Spence 4(q² + q + 1) needs q² + q + 1 = 781; Turyn (p − 1)ᵘ + 1 needs a
    perfect power 3123; the family 3(q + 1) of [PALEY_SKEW.md](../PALEY_SKEW.md) needs 3 | 3124.
  - Skew × amicable products: in 3124 = 2 · 1562 = 4 · 781 = 44 · 71 = 284 · 11 one factor is not a Hadamard order.
  - Good matrices and skew Williamson-type matrices of order 781 are not recorded.

The sources searched are listed in [NOVELTY.md](../../NOVELTY.md). Doubling gives skew-Hadamard matrices of orders
6248 (Table 5's value), 12496, …

## Files

| File | Contents |
|---|---|
| `certificate.json` | the multipliers, the orbit representatives of the −1 positions of each sequence, the row sums |
| `sequences.txt` | the four sequences, one per line |
| `H3124.csv.gz` | the matrix |

## Verification

```
python3 verify/verify.py skew/3124
gzip -dc skew/3124/H3124.csv.gz | ./check_matrix
```

SHA-256 of `H3124.csv`: `30dc5363c61e52c101c25f7df36100d0b3fd534362942613c5b45fded4856fad`
