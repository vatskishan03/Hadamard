# Skew-Hadamard matrix of order 1892

1892 = 4 · 473, with 473 = 11 · 43.

## Construction

A skew-Hadamard matrix is a Hadamard matrix H with H + Hᵀ = 2I.

**Difference family.** Four ±1 sequences x₁, …, x₄ on ℤ₄₇₃ satisfy two conditions:

- each is constant on the orbits of the multiplier subgroup H = ⟨345, 221⟩ ≤ ℤ₄₇₃*, of order 15:
  345 ≡ 1 (mod 43) and is a square mod 11 (of order 5), 221 ≡ 1 (mod 11) and ≡ 6 (mod 43) (of order 3);
- their periodic autocorrelations sum to zero at every nonzero shift.

`certificate.json` lists, for each sequence, the least element r of each orbit rH on which the sequence is −1.

**Skew first block.** The first sequence is of skew type: x₁(0) = 1 and x₁(−t) = −x₁(t) for t ≠ 0. This is
possible because −1 ∉ H, so negation pairs the orbits.

**Row sums.** The row sums are 1, −9, 17, −39, and 1² + 9² + 17² + 39² = 1892.

**The matrix.** Let A, B, C, D be the circulant matrices with these first rows (`sequences.txt`) and R the
back-diagonal identity. The Goethals–Seidel array

```
    [  A     BR     CR     DR  ]
H = [ -BR    A     DᵀR   -CᵀR ]
    [ -CR   -DᵀR    A     BᵀR ]
    [ -DR    CᵀR  -BᵀR    A   ]
```

is then a Hadamard matrix of order 1892. Because A + Aᵀ = 2I, it satisfies H + Hᵀ = 2I (Wallis–Whiteman).

## Search

1. **Structure: a Paley product with a multiplier on K** ([search/NOTES.md](../../search/NOTES.md) §15, §87, §88).
   Over ℤ₄₇₃ = ℤ₄₃ × F₁₁, the squares of F₁₁* act on F₁₁ and 6 (of order 3) acts on ℤ₄₃: 15 orbits of ℤ₄₃ and
   3 · 15 = 45 orbits per block. Each block is u ⊗ 1 + w ⊗ J + z ⊗ χ with χ the Legendre sequence of F₁₁, and the
   counting identity of §15 reads 12 N_w − 4 N₋ = 172.
2. **The solution.** The first block has no w-position; the others have 4, 7 and 10 (orbits of ⟨6⟩ together with 0).
   All 21 are of type u = −2w except one in the third block: N_w = 21, N₋ = 20. The compressions to ℤ₄₃ take the
   values ±1, ±9 and 11, and their autocorrelations sum to 1892 δ.
3. **Two-level search** (`search/skewplug`, `prod_ks.py`): rounds draw an allocation of w-positions, weighted by a
   Gaussian estimate, and an exact row-sum pattern. Level 1 joins the complete list of skew first blocks against
   sampled lists of the other three; level 2 fixes the signs of z by a meet in the middle. The first hit came in the
   tenth round (allocation (0, 4, 7, 10)), after about 5 minutes on one core.
4. **Rate.** Four more matrices of this order came from the same loop with other seeds (allocations (6, 1, 3, 10),
   (0, 3, 7, 7) twice and (6, 1, 7, 7)), about one per 2 · 10⁴ to 3 · 10⁴ level-2 runs. The estimators of this
   repository put the class 3 to 4 decades lower (HEURISTIC).

**Reproduction.**

| | |
|---|---|
| Code | `search/skewplug` (build commands in its README) |
| Command | `python3 -u prod_ks.py k473a 1 9000 "[43]" 11 "[6]" data/cost473.jsonl 40 10000 60000 60000 sl2`, then `python3 cyclic.py out/skew_1892_k473a.json out/skew_1892_cyclic.json` |
| Randomness and caps | rounds are seeded by the seed base and the round index; slots (19 allocations) are drawn with the efficiency weights of `data/cost473.jsonl`; lists of 10000, 60000 and 60000 sampled members; every level-1 instance goes to level 2 |
| Result | the tenth round hits (6390 level-1 instances, 6058 level-2 runs in all, about 5 minutes on one core), and `cyclic.py` gives this certificate without its `sha256` field |

## Status

**New** as a skew-Hadamard matrix. An ordinary Hadamard matrix of order 1892 is known (473 is not in Cati–Pasechnik
Table 4); the skew property is what is new.

- The Đoković–Seberry census lists 473 among its undecided v.
- Cati–Pasechnik Table 5 gives 473(4): the smallest known skew order with odd part 473 was 16 · 473 = 7568. The
  doubling of this matrix, of order 3784, is new as well.
- [search/NOTES.md](../../search/NOTES.md) §15 showed that the other Paley product, over ℤ₁₁ × F₄₃, has no solution
  (2m − 1 = 21 is not a sum of two squares); this one is over ℤ₄₃ × F₁₁, with a multiplier on ℤ₄₃.
- We found no classical construction:
  - Paley I and the products 2ᵗ ∏(qᵢ + 1) (qᵢ ≡ 3 mod 4): the only candidate is q = 1891 = 31 · 61, not a prime power.
  - Szekeres 2(q + 1) needs q = 945, and Wallis–Whiteman 4(q + 1) needs q = 472; neither is a prime power.
  - Skew × amicable products: in 1892 = 2 · 946 = 4 · 473 = 11 · 172 = 22 · 86 = 43 · 44 one factor is not a
    Hadamard order.
  - Good matrices and skew Williamson-type matrices of order 473 are not recorded.

The sources searched are listed in [NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | the four sequences as orbit representatives, H, the row sums |
| `sequences.txt` | the four sequences, one per line |
| `H1892.csv.gz` | the matrix |

## Verification

```
python3 verify/verify.py skew/1892
gzip -dc skew/1892/H1892.csv.gz | ./check_matrix
```

SHA-256 of `H1892.csv`: `8ab883d6a7268bc90086c49403091fc9b96d498d4782adca97b62be521d62572`
