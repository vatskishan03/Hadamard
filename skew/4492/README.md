# Skew-Hadamard matrix of order 4492

4492 = 4 · 1123 (prime).

## Construction

A skew-Hadamard matrix is a Hadamard matrix H with H + Hᵀ = 2I.

**Difference family.** Four ±1 sequences x₁, …, x₄ on ℤ₁₁₂₃ satisfy two conditions:

- each is constant on the orbits of the multiplier subgroup H = ⟨691⟩ ≤ ℤ₁₁₂₃*, of order 33, which is the group of 34th
  powers modulo 1123: {1, 7, 33, 49, 60, 89, 155, 169, 206, 231, 316, 319, 321, 343, 345, 366, 384, 420, 442, 486, 494,
  580, 623, 691, 694, 848, 857, 885, 992, 1032, 1085, 1089, 1110};
- their periodic autocorrelations sum to zero at every nonzero shift.

`certificate.json` lists, for each sequence, one representative r of each orbit rH on which the sequence is −1. The
four sequences are −1 on 17, 17, 17 and 16 of the 34 nonzero orbits (561, 561, 561 and 528 positions), a 4-{1123; 561,
561, 561, 528; 1088} supplementary difference set.

**Skew first block.** The first sequence is of skew type: x₁(0) = 1 and x₁(−t) = −x₁(t) for t ≠ 0. This is possible
because −1 ∉ H (H has odd order), so negation pairs the 34 orbits into 17 pairs and x₁ is −1 on exactly one orbit of
each pair.

**Row sums.** The row sums are 1, 1, 1, 67, and 1² + 1² + 1² + 67² = 4492. A sequence constant on the orbits of H has
row sum ±1 + 33 · (an even number), that is ≡ ±1 (mod 66), and 4492 = a² + b² + c² + d² has only one solution in such
numbers up to order and sign.

**The matrix.** Let A, B, C, D be the circulant matrices with these first rows (`sequences.txt`) and R the
back-diagonal identity. The Goethals–Seidel array

```
    [  A     BR     CR     DR  ]
H = [ -BR    A     DᵀR   -CᵀR ]
    [ -CR   -DᵀR    A     BᵀR ]
    [ -DR    CᵀR  -BᵀR    A   ]
```

is then a Hadamard matrix of order 4492. Because A + Aᵀ = 2I, it satisfies H + Hᵀ = 2I (Wallis–Whiteman).

## Search

1. **Structure.** Over the prime field ℤ₁₁₂₃ each block is a choice of sign on the 34 orbits of H and on the point 0.
   The row sums (1, 1, 1, 67) are forced, which leaves one class of block shapes:
   - 2¹⁷ = 131,072 skew first blocks;
   - C(34, 17) = 2,333,606,220 choices each for B and C;
   - C(34, 16) = 2,203,961,430 choices for D.

   The complementarity condition is Σₖ Pₖ(s) = 0 on the 17 orbits of ⟨H, −1⟩ on the nonzero shifts, of rank 16 once the
   row sums are fixed. The Gaussian estimate of the number of solutions of the class is 10^6.63.
2. **Low-energy k-tree** (`search/cp_lowe.py`, `search/cp_enum.c`, `search/kt4.c`;
   [search/NOTES.md](../../search/NOTES.md) §20–22, §110).
   - `cp_enum.c` enumerates the orbit-sign vectors of B, C and D. It keeps the members whose autocorrelation vector is
     closest to its conditional mean given that the four vectors sum to zero (whitened energy): 97,499,345, 30,445,271
     and 30,439,043 members. All 131,072 skew blocks are kept. The Gaussian estimate inside these lists is 10^3.26
     solutions; it is a heuristic, not a rate.
   - `kt4.c` joins the four lists over 10⁵ split values of five merge classes (1, 5, 10, 11, 12), ranked by expected
     solutions per unit of work, with about 4 · 10⁵ pairs on the AB side and on average 2.8 · 10⁸ (7 · 10⁷ to 5 · 10⁸
     for 80 % of them) on the CD side of each. The first 818 logged split values gave nothing. The run stopped at its
     limit of 3 solutions, all verified.
3. **Earlier attempt.** A run with smaller lists (3 · 10⁷, 10⁷ and 10⁷ members, 6000 split values, 4932 s on 3 threads)
   found no solution. Lists about three times larger, as above, produced three solutions (the run stopped at its
   limit K = 3).

**Reproduction.**

| | |
|---|---|
| Code | commit `919e73d` (`search/cp_lowe.py`, `search/cp_enum.c`, `search/kt4.c`; `search/skew_certify.py` for the certificate) |
| Command | `python3 search/cp_lowe.py 1123 --q 1123 --e 34 --config 0 --L 1e8,3e7,3e7 --nsig 100000 --threads 8 --K 3 --seed 1 --groups 4 --samples 200000 --mcpairs 10000000 --maxpairs 300000000 --stream` |
| Randomness and caps | `--seed 1` (numpy `default_rng`) draws the 2 · 10⁵ samples per block that set the energy thresholds and the 10⁷ pairs per side that rank the split values; `--L 1e8,3e7,3e7` the target number of members kept per non-skew block (the threshold is a quantile of the sampled energies, so the lists hold 97,499,345, 30,445,271 and 30,439,043); `--nsig 100000` split values; `--K 3` stops kt4 after 3 solutions |
| Result | 3 verified solutions after 1230 s on 8 threads (block enumeration and list building 159 s, then about 1.2 s per split value); about 11 GB of memory. The order of kt4's output depends on the thread schedule; the certificate is the third hit logged. |
| Certificate | `search/skew_certify.py` turns a logged hit into the files in this directory |

## Status

**New** as a skew-Hadamard matrix. An ordinary Hadamard matrix of order 4492 is known (see below); the skew property is
what is new.

- The Đoković–Seberry census lists 1123 among its 737 undecided odd v < 2500 for skew matrices, and not in its Paley,
  Szekeres, Szekeres–Whiteman or Spence series or among its sporadic values.
- Every other skew table we found stops below 1123: Cati–Pasechnik Table 5 and Seberry–Yamada (1992) Table 7.1 at 999,
  Seberry–Yamada (Utilitas Math. 1991) at 499, and Koukouvinos–Stylianou (2008) at 493.
- Before this, the smallest skew order with odd part 1123 that we found from a published construction was
  2⁷ · 1123 = 143744 (Paley I: 143743 is prime and ≡ 3 (mod 4)). For N = 4492, 8984, 17968, 35936 and 71872 none of
  the numbers N / 2ᵘ − 1 (u = 0, 1, …) is a prime power, so the Paley, Szekeres and Wallis–Whiteman shapes do not
  apply.
- The ordinary matrix of order 4492 is a Goethals–Seidel matrix with the same multiplier group of order 33, constructed
  independently here (4-{1123; 561, 561, 561, 528; 1088}, no skew block) and by Egan, Kotsireas, London and Ó Catháin
  (arXiv:2610.27175, 19 October 2026). We found no statement in the survey that its matrix is skew, and its list of
  known skew orders (Paley's q + 1, doubling, Spence's family, 4n for odd n < 250, 276 and 292, with 356 the smallest
  open order) does not reach 4492. The Hadamard Matrix Explorer builds 4492 as `G(A(1123))`, from supplementary
  difference sets, and not from its skew-type sequences `SK`. Decoded from the Explorer's image on 2 November 2026,
  their matrix is a Goethals–Seidel matrix that is not skew ([NOVELTY.md](../../NOVELTY.md)).
- No classical construction gives a skew matrix of order 4492:
  - Paley I: 4491 = 3² · 499 is not a prime power.
  - Szekeres, Whiteman, Spence and Momihara–Xiang 2(qᵉ + 1): 2245 = 5 · 449 is not a prime power.
  - Wallis–Whiteman and Momihara–Xiang 4(qᵉ + 1): 1122 = 2 · 3 · 11 · 17 is not a prime power.
  - Spence 4(q² + q + 1): 1123 = 33² + 33 + 1, but 33 is not a prime power.
  - Turyn (h − 1)ᵘ + 1 from a skew matrix of order h, u odd: 4491 = 3² · 499 is not a perfect power, so only u = 1 (the
    order itself) remains.
  - Doubling needs a skew matrix of order 2246, which is not a multiple of 4.
  - Skew × amicable products need 4492 = mn with m and n Hadamard orders, and every nontrivial factorisation has a
    factor 1123 or 2246.
  - Seberry's mn(n − 1) is possible only with n = 2, which needs an OD(4492; 1, 2246, 2245).
  - The 3(q + 1) theorem needs 3 | 4492. Blatt–Szekeres needs complementary difference sets in a group of order
    2 · 1122 + 1 = 2245 = 5 · 449; the known families need 2245 or 2 · 2245 + 1 = 4491 = 3² · 499 to be a prime
    power, and neither is. Williamson's h(pᵏ + 1), with pᵏ + 1 ≡ 0 (mod 4), needs pᵏ + 1 = 4, that is a skew matrix of
    the odd order h = 1123, which cannot exist.
  - Good matrices and skew Williamson-type matrices of order 1123 are not recorded.

The sources searched are listed in [NOVELTY.md](../../NOVELTY.md). Doubling (Seberry) gives skew-Hadamard matrices of
orders 8984, 17968, 35936 and 71872, which are new as far as we know; 143744 was already known.

## Files

| File | Contents |
|---|---|
| `certificate.json` | the four sequences as orbit representatives, H, the row sums |
| `sequences.txt` | the four sequences, one per line |
| `H4492.csv.gz` | the matrix |

## Verification

```
python3 verify/verify.py skew/4492
gzip -dc skew/4492/H4492.csv.gz | ./check_matrix
```

SHA-256 of `H4492.csv`: `90d6bc84e8534dcdda3e75b8e62434ef480cde10a63209f0a2418eb9ceeee80c`
