# Skew-Hadamard matrix of order 4712

4712 = 8 · 589, with 589 = 19 · 31.

## Construction

A skew-Hadamard matrix is a Hadamard matrix H with H + Hᵀ = 2I.

**Sequences.** Four ±1 sequences of length 1178 = 2 · 589 (`sequences.txt`) whose negaperiodic autocorrelations sum to
4712 at shift 0 and to 0 elsewhere. The first is palindromic with first entry 1, so its negacyclic matrix A satisfies
A + Aᵀ = 2I. The sums of the four sequences are −10, 34, −46, 12 (negacyclic blocks have no constant row sum). The
negacyclic Goethals–Seidel array of [NEGACYCLIC.md](../NEGACYCLIC.md) then gives the matrix.

**Quaternary form.** The four sequences are u_X on ℤ₅₈₉ with values in {1, i, −1, −i}, constant on the orbits of
H = ⟨311, 381⟩ ≤ ℤ₅₈₉*, of order 45 (21 orbits), with Σ u_X u_X* = 2356 δ, u_A(−g) = i · conj(u_A(g)) for g ≠ 0 and
u_A(0) = 1. Over ℤ₅₈₉ = ℤ₁₉ × F₃₁, H is ⟨7⟩ (order 3) on ℤ₁₉ times the squares of F₃₁* (order 15): 311 = (7, 1) and
381 = (1, 9) mod (19, 31). The block sums are 19 + 18i, 31, −5 and −26 + 3i (685 + 961 + 25 + 685 = 2356).
`certificate.json` records the orbit values as exponents of i (`quaternary_orbit_exponents`).

## Search

1. **Structure: a Paley plug-in in the negacyclic ring** ([search/NOTES.md](../../search/NOTES.md) §89), over
   ℤ₁₉ × F₃₁ with 7 orbits of ⟨7⟩ on ℤ₁₉ (21 orbits per block).
2. **The solution.** All defects sit at k = 0: the half defect of A (weight (q + 1)/2 = 16) and full defects in B
   (q + 1 = 32) and D (q − 3 = 28), none in C: 16 + 32 + 28 = 76 = 4 · 19.
3. **Two-level search** (`search/negaskew/qp`). Level 1 is exhaustive: 256 skew first blocks, 36,863 other blocks,
   261,984 level-1 instances (23 s). Level 2 on a 1% sample of them (2540 instances, seed 3) completed 84, with 168
   solutions: the class is rich. This solution is the first one found in that sample, after 2.4 s.

**Reproduction.**

| | |
|---|---|
| Code | `search/negaskew` (`qp.c`, `qshape.py`, `plugcert.py`; build commands in its README) |
| Command | `python3 qshape.py '[19]' 31 '[7]' shapes/p589b.txt`, `./qp shapes/p589b.txt l2c=1 l2p=0.01 first=1 maxprint=1 seed=3`, then `python3 plugcert.py 19 31 "[7]" "0031210\|1123120\|0301302 0010322\|0031232\|0213010 0223120\|2100213\|0322031 0311200\|2013232\|2231010"` |
| Randomness and caps | `seed=3` draws the 1% of level-1 instances (`l2p=0.01`) that go to level 2; `first=1` stops at the first completed instance |
| Result | the solution above in 2.4 s; `plugcert.py` writes `out/cert_4712.json`, identical to `certificate.json` |

## Status

**New.** An ordinary Hadamard matrix of order 4712 is known (589 is not in Cati–Pasechnik Table 4); the skew property
is what is new.

- Cati–Pasechnik Table 5 gives 589(5), and [skew/9424](../9424) of this repository (Wallis's product, 16 · 589) made it
  589(4); this matrix makes it 589(3). 589 is undecided in the Đoković–Seberry census (order 4 · 589), which this
  matrix does not settle.
- SageMath's `skew_hadamard_matrix` (development source of 21 October 2026, read, not run) has no construction of
  order 4712.
- We found no classical construction:
  - Paley I needs 4711 = 7 · 673; Szekeres and Momihara–Xiang 2(qᵉ + 1) need 2355 = 3 · 5 · 157; Wallis–Whiteman and
    Momihara–Xiang 2ᵘ(qᵉ + 1) need 1177 = 11 · 107; none is a prime power.
  - Spence 4(q² + q + 1) would need q² + q + 1 = 1178, and Turyn (p − 1)ᵘ + 1 a perfect power 4711: neither exists.
  - Seberry's mn(n − 1) construction only fits n = 2, which needs an OD(4712; 1, 2356, 2355).
  - Doubling needs a skew-Hadamard matrix of order 2356 (census-undecided). Skew × amicable products have a factor
    that is not a Hadamard order (4 · 1178, 8 · 589, 76 · 62, 124 · 38, 152 · 31, 248 · 19), except 2 · 2356, which is
    doubling.

The sources searched are listed in [NOVELTY.md](../../NOVELTY.md). Doubling gives the skew orders 9424 (already in
this repository), 18848, …

## Files

| File | Contents |
|---|---|
| `certificate.json` | the four ±1 sequences, the quaternary orbit exponents and H, the sequence sums |
| `sequences.txt` | the four sequences, one per line |
| `H4712.csv.gz` | the matrix |

## Verification

```
python3 verify/verify.py skew/4712
gzip -dc skew/4712/H4712.csv.gz | ./check_matrix
```

SHA-256 of `H4712.csv`: `84f9d41f13e6fdfc302d6725162c0e6a9806f07f2d33e5de550538bbf4a496a7`
