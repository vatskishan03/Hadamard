# Skew-Hadamard matrix of order 2452

2452 = 4 · 613 (prime).

## Construction

A skew-Hadamard matrix is a Hadamard matrix H with H + Hᵀ = 2I.

**Difference family.** Four ±1 sequences x₁, …, x₄ on ℤ₆₁₃ satisfy two conditions:

- each is constant on the orbits of the multiplier subgroup H = ⟨37⟩ ≤ ℤ₆₁₃*, of order 17: {1, 37, 116, 143, 171, 190,
  197, 198, 220, 227, 287, 387, 430, 546, 583, 585, 586};
- their periodic autocorrelations sum to zero at every nonzero shift.

`certificate.json` lists, for each sequence, one representative r of each orbit rH on which the sequence is −1.

**Skew first block.** The first sequence is of skew type: x₁(0) = 1 and x₁(−t) = −x₁(t) for t ≠ 0. This is possible
because −1 ∉ H, so negation pairs the orbits.

**Row sums.** The row sums are 1, 1, 35, 35, and 1² + 1² + 35² + 35² = 2452.

**The matrix.** Let A, B, C, D be the circulant matrices with these first rows (`sequences.txt`) and R the
back-diagonal identity. The Goethals–Seidel array

```
    [  A     BR     CR     DR  ]
H = [ -BR    A     DᵀR   -CᵀR ]
    [ -CR   -DᵀR    A     BᵀR ]
    [ -DR    CᵀR  -BᵀR    A   ]
```

is then a Hadamard matrix of order 2452. Because A + Aᵀ = 2I, it satisfies H + Hᵀ = 2I (Wallis–Whiteman).

## Search

1. **Ranking.** `search/skew_rank.py` over the census primes below 1000 put 613 first: H = ⟨37⟩ (order 17, 37 orbits,
   18 autocorrelation classes), rows (1, 1, 35, 35), about 10^8.3 solutions in the structure (v is prime, where the
   Gaussian model has been reliable).
2. **Lists of unequal size** (`search/skew_lowE.py`):
   - all 262,144 skew first blocks;
   - 1,290,242 low-energy candidates with row sum 1, from `search/skew_gen.c`: local minima of the distance of the
     autocorrelations from their mean, with three descent steps;
   - 7,075,788 low-energy candidates with row sum 35, one list shared by the roles C and D.

   The script's Gaussian estimate, from the lists' own class statistics, was 10^1.65 solutions inside the lists.
3. **Four-list k-tree matcher** (`search/match4.c`). Pairs (A, D) and (B, C) are formed with the partial sums on the
   first 6 classes fixed to τ, taken in decreasing order of probability, and joined on the other 12. It stops after 4
   solutions (`-k 4`), and it reached them after 12,792 s; all 4 were verified. The first 40,000 values of τ carry
   only about 4% of the τ-weight, which predicts about 2 solutions, so the yield is in line with the estimate.

**Reproduction.**

| | |
|---|---|
| Code | commit `9ef2547` |
| Command | `python3 search/skew_lowE.py 613 37 1,35,35 --nbc 1.5e6 --nd 1.6e7 --steps 3 --maxtau 40000 --top 10 --seed 23` |
| Randomness and caps | seed 23; the candidate lists also depend on the OpenMP thread schedule of skew_gen, so a rerun gives different lists with the same statistics |
| Result | 4 verified solutions after 12,792 s of matching on 4 cores (26 minutes of list generation before that) |

## Status

**New** as a skew-Hadamard matrix. A symmetric Hadamard matrix of order 2452 is known (Turyn); no skew one was.

- Cati–Pasechnik Table 5 lists 613(3) as a value carried over from Seberry–Yamada. Paley I gives 4904, since 4903 is
  prime, so the smallest known skew order with odd part 613 was 2³ · 613 = 4904. The order is above 1208, the range
  they checked against implemented constructions.
- The Đoković–Seberry census lists 613 among its 737 undecided v. 613 is in none of its Paley, Szekeres,
  Szekeres–Whiteman or Spence series nor among its sporadic values.
- SageMath's `unknown_skew_hadamard` list stops at 1196, below this order.
- No infinite family reaches it:
  - Paley I (2451 = 3 · 19 · 43).
  - Szekeres, Whiteman and Momihara–Xiang 2(q + 1), with q = 1225 = 35².
  - Wallis–Whiteman 4(q + 1), with q = 612.
  - (|t| + 1)(q + 1).
  - Spence 4(q² + q + 1).
  - (p − 1)ᵘ + 1.
  - Products and amicable products.
  - Đoković's Szekeres + D-optimal route, which would need a D-optimal family of length 613; none is known.

The sources searched are listed in [NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | the four sequences as orbit representatives, H, the row sums |
| `sequences.txt` | the four sequences, one per line |
| `H2452.csv.gz` | the matrix |

## Verification

```
python3 verify/verify.py skew/2452
gzip -dc skew/2452/H2452.csv.gz | ./check_matrix
```

SHA-256 of `H2452.csv`: `1778d833c89379a0ad13c21847fffa98e6a611c2d6861b1779af34fcf55f0415`
