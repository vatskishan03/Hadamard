# Skew-Hadamard matrix of order 476

476 = 4 · 119, with 119 = 7 · 17.

## Construction

A skew-Hadamard matrix is a Hadamard matrix H with H + Hᵀ = 2I.

**Difference family.** Four ±1 sequences x₁, …, x₄ on ℤ₁₁₉ satisfy two conditions:

- each is constant on the orbits of the multiplier subgroup H = ⟨86⟩ = {1, 18, 86} ≤ ℤ₁₁₉*, of order 3:
  H = {h : h ≡ 1 (mod 17), h mod 7 a square};
- their periodic autocorrelations sum to zero at every nonzero shift.

`certificate.json` lists, for each sequence, the least element r of each orbit rH on which the sequence is −1.

**Skew first block.** The first sequence is of skew type: x₁(0) = 1 and x₁(−t) = −x₁(t) for t ≠ 0. This is
possible because −1 ∉ H, so negation pairs the orbits.

**Row sums.** The row sums are 1, −5, 15, −15, and 1² + 5² + 15² + 15² = 476.

**The matrix.** Let A, B, C, D be the circulant matrices with these first rows (`sequences.txt`) and R the
back-diagonal identity. The Goethals–Seidel array

```
    [  A     BR     CR     DR  ]
H = [ -BR    A     DᵀR   -CᵀR ]
    [ -CR   -DᵀR    A     BᵀR ]
    [ -DR    CᵀR  -BᵀR    A   ]
```

is then a Hadamard matrix of order 476. Because A + Aᵀ = 2I, it satisfies H + Hᵀ = 2I (Wallis–Whiteman).

## Search

1. **Structure: a Paley product** ([search/NOTES.md](../../search/NOTES.md) §15). Over ℤ₁₁₉ = ℤ₁₇ × F₇, H acts on
   F₇ as the squares and trivially on ℤ₁₇: 51 orbits per block. With χ the Legendre sequence of F₇, every block is
   u ⊗ 1 + w ⊗ J + z ⊗ χ with u, w, z on ℤ₁₇. The counting identity (q + 1) N_w − 4 N₋ = 4m of §15 reads
   2 N_w − N₋ = 17, so 9 to 17 of the 68 positions are w-positions. The search of §15 handled at most five, and this
   structure had never been searched.
2. **The skew block.** In the first block w is odd and z even under a ↦ −a, so it has no w-position at 0 and its
   w-positions come in pairs ±a.
3. **The solution.** The w-positions are a = 7, 10 (first block), 0, 13 (second), 8, 16 (third) and 7, 12, 15
   (fourth): N_w = 9, and N₋ = 1 (at a = 8 of the third block). The compressions to ℤ₁₇ (sums over the fibres of
   F₇) take the values ±1, ±7 and ±5, and their autocorrelations sum to 476 δ.
4. **Two-level search** (`search/skewplug`, [search/NOTES.md](../../search/NOTES.md) §87). Level 1 joins the complete
   list of skew first blocks against sampled lists of the other three blocks (one allocation of w-positions and one
   exact row-sum pattern per round) on the compressed equation over ℤ₁₇. Level 2 fixes the signs of z off the
   w-positions by a meet in the middle. The first round gave the matrix. Solutions are plentiful in this structure:
   three more were found in other allocations of w-positions and a fourth by an independent judge with its own code,
   and the estimated total is 10¹⁹ to 10²⁰ ordered quadruples (HEURISTIC).

**Reproduction.**

| | |
|---|---|
| Code | `search/skewplug` (build commands in its README) |
| Command | `python3 -u sk_worker.py w119a 1 900 "[17]" 7 "[[2,2,2,3]]" 2000 20000 20000 30 1.0`, then `python3 cyclic.py out/skew_476_w119a.json out/skew_476_cyclic.json` |
| Randomness and caps | rounds are seeded by the seed base 1 and the round index; lists of 2000, 20000 and 20000 sampled members for the second, third and fourth blocks; 30 level-1 instances per round go to level 2 |
| Result | the first round hits, in 19 s on one core; the output over ℤ₁₇ × ℤ₇ (SHA-256 of its matrix `1f4f883c…`) becomes this certificate under x ↦ (x mod 17, x mod 7) |

## Status

**New** as a skew-Hadamard matrix. An ordinary Hadamard matrix of order 476 is classical (Cooper–Wallis with
T-matrices of order 7); the skew property is what is new.

- The Đoković–Seberry census lists 119 among its 737 undecided v; it is the fourth smallest, after 89, 101 and 107.
- Cati–Pasechnik Table 5 gives 119(4): the smallest known skew order with odd part 119 was 16 · 119 = 1904. So the
  doubling of this matrix, a skew-Hadamard matrix of order 952, is new as well.
- SageMath lists 476 in `unknown_skew_hadamard`. The survey of Egan, Kotsireas, London and Ó Catháin
  (arXiv:2610.27175, 19 October 2026) gives 356 as the smallest open skew order and no construction of 476.
- No classical construction gives a skew matrix of order 476:
  - Paley I and the products 2ᵗ ∏(qᵢ + 1) (qᵢ ≡ 3 mod 4): every factor qᵢ + 1 is divisible by 4 and 476 = 4 · 119,
    so the only candidate is q = 475 = 5² · 19, which is not a prime power.
  - Szekeres 2(q + 1) needs q = 237 = 3 · 79, and Wallis–Whiteman 4(q + 1) needs q = 118; neither is a prime power.
  - Skew × amicable products: 476 is not a product of two Hadamard orders greater than 1 (2 · 238, 4 · 119 and 28 · 17
    each have a factor that is not one).
  - A cyclic low-energy search over ℤ₁₁₉ ([search/NOTES.md](../../search/NOTES.md) §10) found nothing.

The sources searched are listed in [NOVELTY.md](../../NOVELTY.md). Doubling (Seberry) gives skew-Hadamard matrices of
orders 952, 1904 (already known), …; the doubling of order 952 gives the symmetric Hadamard matrix of order 7612 by
the symmetric Yamada–Kiyasu theorem ([symmetric/7612](../../symmetric/7612)).

## Files

| File | Contents |
|---|---|
| `certificate.json` | the four sequences as orbit representatives, H, the row sums |
| `sequences.txt` | the four sequences, one per line |
| `H476.csv.gz` | the matrix |

## Verification

```
python3 verify/verify.py skew/476
gzip -dc skew/476/H476.csv.gz | ./check_matrix
```

SHA-256 of `H476.csv`: `767afb6f57c5f3e27dd251890cca26ac8b603a3db6de6d3aa125ff0fbc136443`
