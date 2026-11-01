# Skew-Hadamard matrix of order 2212

2212 = 4 · 553, with 553 = 7 · 79.

## Construction

A skew-Hadamard matrix is a Hadamard matrix H with H + Hᵀ = 2I.

**Difference family.** Four ±1 sequences x₁, …, x₄ on ℤ₅₅₃ ≅ ℤ₇ × ℤ₇₉ satisfy two conditions:

- each is constant on the orbits of the multiplier subgroup H = ⟨176⟩ ≤ ℤ₅₅₃*, of order 13: {1, 8, 22, 64, 141, 176,
  204, 225, 302, 337, 484, 512, 526} = {h : h ≡ 1 (mod 7), h mod 79 a sixth power};
- their periodic autocorrelations sum to zero at every nonzero shift.

H has 49 orbits on ℤ₅₅₃: the 7 fixed points 0, 79, 158, …, 474 (the elements (a, 0) of ℤ₇ × ℤ₇₉) and 42 orbits of size
13 (for each a ∈ ℤ₇, the six sets {a} × (a coset of the sixth powers in ℤ₇₉*)). `certificate.json` lists, for each
sequence, one representative r of each orbit rH on which the sequence is −1. The four sequences are −1 on 24, 21, 25
and 25 of the 49 orbits (276, 261, 289 and 289 positions), a 4-{553; 276, 261, 289, 289; 562} supplementary difference
set.

**Skew first block.** The first sequence is of skew type: x₁(0) = 1 and x₁(−t) = −x₁(t) for t ≠ 0. This is possible
because −1 ∉ H (H has odd order): −1 = 3³⁹ modulo 79 is not a sixth power, so negation maps the orbit of (a, b) to a
different orbit, and pairs the 48 nonzero orbits into 24 pairs (3 pairs of fixed points, 21 pairs of size 13). x₁ is −1
on exactly one orbit of each pair.

**Row sums.** The row sums are 1, 31, −25, −25, and 1² + 31² + 25² + 25² = 2212. A sequence constant on the orbits of H
has row sum u + 26w, where u is the sum over the 7 fixed points (odd, |u| ≤ 7) and w an integer. A skew first block has
row sum 1, and 2212 = 1 + b² + c² + d² has three solutions in numbers of this form, up to order and sign: (b, c, d) =
(1, 1, 47), (23, 29, 29) and (25, 25, 31). The compression to ℤ₇ (below) admits only the last two.

**The matrix.** Let A, B, C, D be the circulant matrices with these first rows (`sequences.txt`) and R the
back-diagonal identity. The Goethals–Seidel array

```
    [  A     BR     CR     DR  ]
H = [ -BR    A     DᵀR   -CᵀR ]
    [ -CR   -DᵀR    A     BᵀR ]
    [ -DR    CᵀR  -BᵀR    A   ]
```

is then a Hadamard matrix of order 2212. Because A + Aᵀ = 2I, it satisfies H + Hᵀ = 2I (Wallis–Whiteman).

## Search

1. **Structure: a cyclotomic product** ([search/NOTES.md](../../search/NOTES.md) §16, §20). Over ℤ₇ × GF(79), H acts
   trivially on ℤ₇ and on GF(79) as the sixth powers, so the classes of index e = 6 (class size F = 13) are fixed. A
   block is given at each position a ∈ ℤ₇ by its value f(a) at (a, 0) and its values on the six classes (a, C_j): 49
   H-orbits per block.
2. **Stage 1: the compressions to ℤ₇.** The compression is y(a) = f(a) + 13 s(a), with s(a) the sum of the six class
   values at a, and Σ_k Σ_a y_k(a)² = 4 · 553 becomes Σ over the 28 positions of (13 s² + 2 f s) = 168. A position with
   s ≠ 0 costs 48 or 56 for s = ±2, and at least 200 for |s| ≥ 4, so 168 = 3 · 56: exactly three spikes y = ±27 fit.
   Stage 1 lists 36 classes (the three non-skew blocks are taken up to negation, translation, reversal and order, the
   skew block is not reduced): 18 with row-sum moduli (1, 31, 25, 25) and 18 with (1, 29, 29, 23), each with one spike
   in every non-skew block and none in the skew block. The class searched has
   y = (1, 1, −1, 1, −1, 1, −1), (1, 27, −1, 1, 1, 1, 1), (−27, 1, 1, 1, 1, −1, −1), (−27, 1, 1, −1, 1, 1, −1), one
   spike in each non-skew block, and row sums (1, 31, −25, −25); the compressions of the sequences in `sequences.txt`
   are exactly these vectors. The Gaussian estimate of its number of solutions is 10^4.06. A recomputation of stage 1
   on another machine (arm64 macOS, the same numpy 1.26.4) gives the same 36 classes in the same order with the same
   skew blocks, but may print another representative of a non-skew block; there `--config 10` shows
   (1, 1, 27, −1, 1, 1, 1) for the second block, a translate of the vector above. The vectors above are the ones
   searched.
3. **Lists.** The blocks of this class have 64,000 (skew block: 8 · 20³) and 960,000,000 (15 · 20⁶, each) members. The
   key of a member is 21-dimensional (3 coordinates for the compression to GF(79) and 18 for the autocorrelations at the
   shifts (b, C_c), 1 ≤ b ≤ 3), of rank 20.
4. **Low-energy k-tree** (`search/cp_lowe.py`, `search/cp_enum.c`, `search/kt4.c`;
   [search/NOTES.md](../../search/NOTES.md) §20, §110).
   - `cp_enum.c` enumerates every member of each non-skew block. It keeps the members whose autocorrelation vector is
     closest to its conditional mean given that the four vectors sum to zero (whitened energy): 200,347,074,
     9,911,542 and 10,163,372 members. All 64,000 skew blocks are kept. The Gaussian estimate inside these lists is
     10^1.28 solutions; it is a heuristic, not a rate.
   - `kt4.c` joins the four lists over 10⁵ split values of five merge classes (0, 1, 2, 16, 18), chosen from 1,223,596
     sampled values and ranked by expected solutions per unit of work (they hold 0.335 of the expected solutions), with
     on average 3.5 · 10⁶ pairs on the AB side and 8 · 10⁷ on the CD side of each. The first 935 split values gave
     nothing. The 936th gave one solution, verified, and the run stopped at its limit of 1 solution after 3.3 · 10⁹ AB
     pairs and 7.6 · 10¹⁰ CD pairs.
5. **Earlier attempt.** A run on the same class with smaller lists (2 · 10⁷, 10⁶ and 10⁶ members, Gaussian estimate
   10^−1.20 inside them, 2000 split values, 1 thread, 234 s) found no solution. Lists ten times larger, as above,
   produced one.

**Reproduction.**

| | |
|---|---|
| Code | commit `919e73d` (`search/cp_lowe.py`, `search/cp_enum.c`, `search/kt4.c`; `search/skew_certify.py` for the certificate) |
| Command | `python3 search/cp_lowe.py 553 --q 79 --e 6 --config 10 --L 2e8,1e7,1e7 --nsig 100000 --threads 12 --K 1 --seed 1 --groups 4 --samples 500000 --mcpairs 10000000 --maxpairs 300000000` |
| Randomness and caps | `--seed 1` (numpy `default_rng`) draws the 5 · 10⁵ samples per block that set the energy thresholds and the 10⁷ pairs per side that rank the split values; `--L 2e8,1e7,1e7` the target number of members kept per non-skew block (the threshold is a quantile of the sampled energies, so the lists hold 200,347,074, 9,911,542 and 10,163,372); `--nsig 100000` split values; `--K 1` stops kt4 at the first solution; `--kgen 1 --twist 0` (the defaults) leave ℤ₇ unmultiplied and untwisted; the same command computes the 36 stage-1 classes, and `--config 10` is the eleventh in their order |
| Result | one verified solution after 635 s on 12 threads, about 18 GB of memory. Elapsed times from the start: stage 1 21 s, the four lists 44 s, the ranking of the split values 72 s, kt4 started at 81 s and needed about 65 s more to map and index the lists; its 936 split values then took 489 s (about 0.5 s each). |
| Certificate | `search/skew_certify.py` turns the logged hit into the files in this directory |

## Status

**New** as a skew-Hadamard matrix. An ordinary Hadamard matrix of order 2212 is classical (Seberry 1973, Hadamard
matrices of orders 14(q + 1) for q ≡ 1 (mod 4) a prime power, here q = 157; the abstract lists 2212; the Hadamard
Matrix Explorer builds it as `G(T(7), W(79))`); the skew property is what is new.

- The Đoković–Seberry census lists 553 among its 737 undecided odd v < 2500 for skew matrices, and not in its Paley,
  Szekeres, Szekeres–Whiteman or Spence series or among its sporadic values.
- Cati–Pasechnik Table 5 (v2, 30 August 2025; entries n(m) give the smallest m with a skew matrix of order 2ᵐn known,
  and m = 2 is not shown) lists 553(3): the smallest known skew order with odd part 553 was 8 · 553 = 4424, from Paley I
  (4423 is prime and ≡ 3 (mod 4)). Koukouvinos–Stylianou (2008), Table 1 stops at 493, before 553.
- The survey of Egan, Kotsireas, London and Ó Catháin (arXiv:2610.27175, 19 October 2026) lists in its Table 3 as known
  skew orders Paley's q + 1, doubling, Spence's family, 4n for odd n < 250, 276 and 292, with 356 the smallest open
  order; neither 2212 nor 553 occurs anywhere in the paper. SageMath's `skew_hadamard_matrix` implements the known
  orders up to 1200 and only some above (Paley I, Spence 1975 and 1977, complementary difference sets, orthogonal
  designs, doubling, Williamson, stored small cases); none of these gives 2212 (see below).
- 553 = 23² + 23 + 1 has the shape of Spence's family 4(q² + q + 1) with q = 23. Spence's theorem (Canad. J. Math. 27
  (1975) 555–560) needs 1 + q + q² a prime ≡ 3, 5 or 7 (mod 8) or 3 + 2q + 2q² a prime power; here 553 = 7 · 79 is not a
  prime and 3 + 2q + 2q² = 1107 = 3³ · 41 is not a prime power. Its general form (Theorem 3 of the paper) needs two
  complementary difference sets in ℤ₅₅₃ together with a cyclic projective plane of order q² = 529, which exists because
  q = 23 is a prime power; no such sets are known: searches for them with multiplier groups of orders 39 and 117 are
  empty. The paper's own list of the q it covers (2, 3, 4, 5, 13, 16, 17, 25, 27, 31, …, orders 28, 52, 84, 124, 732,
  1092, 1228, 2604, 3028, 3972, …) omits 23.
- No classical construction gives a skew matrix of order 2212:
  - Paley I: 2211 = 3 · 11 · 67 is not a prime power.
  - Szekeres, Whiteman, Spence and Momihara–Xiang 2(qᵉ + 1): 1105 = 5 · 13 · 17 is not a prime power.
  - Wallis–Whiteman and Momihara–Xiang 4(qᵉ + 1): 552 = 2³ · 3 · 23 is not a prime power.
  - Turyn (h − 1)ᵘ + 1 from a skew matrix of order h, u odd: 2211 is not a perfect power, so only u = 1 (the order
    itself) remains.
  - Doubling needs a skew matrix of order 1106, which is not a multiple of 4.
  - Skew × amicable products need 2212 = mn with m and n Hadamard orders, and the only such factorisation is
    1 · 2212 (the others have a factor 553, 1106, 7, 14, 79 or 158). Williamson's h(pᵏ + 1), with pᵏ + 1 ≡ 0 (mod 4),
    needs pᵏ + 1 = 4 or 28, that is a skew matrix of the odd order h = 553 or 79, which cannot exist.
  - Seberry's mn(n − 1) is possible only with n = 2, which needs an OD(2212; 1, 1106, 1105).
  - The 3(q + 1) theorem needs 3 | 2212. Blatt–Szekeres needs complementary difference sets in a group of order
    2 · 552 + 1 = 1105 = 5 · 13 · 17; the known families need 1105 or 2 · 1105 + 1 = 2211 to be a prime power, and
    neither is.
  - Good matrices and skew Williamson-type matrices of order 553 are not recorded.

The sources searched are listed in [NOVELTY.md](../../NOVELTY.md). Doubling (Seberry) gives skew-Hadamard matrices of
orders 4424, 8848, 17696, …, all already known: Paley I gives 4424, and doubling the rest.

## Files

| File | Contents |
|---|---|
| `certificate.json` | the four sequences as orbit representatives, H, the row sums |
| `sequences.txt` | the four sequences, one per line |
| `H2212.csv.gz` | the matrix |

## Verification

```
python3 verify/verify.py skew/2212
gzip -dc skew/2212/H2212.csv.gz | ./check_matrix
```

SHA-256 of `H2212.csv`: `edae065395b6ed1052b13f0ccac80fecc1e91e8af5185a2ca88d8ec6f4d5e02b`
