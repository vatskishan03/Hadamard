# Skew-Hadamard matrix of order 2392

2392 = 8 · 299, with 299 = 13 · 23.

## Construction

A skew-Hadamard matrix is a Hadamard matrix H with H + Hᵀ = 2I.

**Sequences.** Four ±1 sequences of length 598 = 2 · 299 (`sequences.txt`) whose negaperiodic autocorrelations sum to
2392 at shift 0 and to 0 elsewhere. The first is palindromic with first entry 1, so its negacyclic matrix A satisfies
A + Aᵀ = 2I. The sums of the four sequences are −60, −8, 12, 68 (negacyclic blocks have no constant row sum). The
negacyclic Goethals–Seidel array of [NEGACYCLIC.md](../NEGACYCLIC.md) then gives the matrix.

**Quaternary form.** The four sequences are u_X on ℤ₂₉₉ with values in {1, i, −1, −i}, constant on the orbits of
H = ⟨185, 209⟩ ≤ ℤ₂₉₉*, of order 33 (15 orbits), with Σ u_X u_X* = 1196 δ, u_A(−g) = i · conj(u_A(g)) for g ≠ 0 and
u_A(0) = 1. Over ℤ₂₉₉ = ℤ₁₃ × F₂₃, H is ⟨3⟩ (order 3) on ℤ₁₃ times the squares of F₂₃* (order 11): 209 = (1, 2) and
185 = (3, 1) mod (13, 23). The block sums are −10 − 11i, 4 − 3i, −24 − 3i and −2 − 19i (221 + 25 + 585 + 365 = 1196).
`certificate.json` records the orbit values as exponents of i (`quaternary_orbit_exponents`).

## Search

1. **Structure: a Paley plug-in in the negacyclic ring** ([search/NOTES.md](../../search/NOTES.md) §89). Each block
   is s ⊗ δ + t ⊗ Q + e ⊗ J over ℤ₁₃ × F₂₃ with Gaussian-integer s, t, e: at each position k of ℤ₁₃ it takes the value
   x₀ at (k, 0), x₁ on (k, squares) and x₂ on (k, non-squares). The compressions w = x₀ + 22 e to ℤ₁₃ satisfy
   Σ_X w_X w_X* = 1196 δ (level 1), and t follows by a meet in the middle (level 2).
2. **The solution.** All defects (positions with x₂ ≠ −x₁) sit at k = 0: a half defect in A (x₂ = −i x₁, the one the
   skew condition forces), full defects in C and D (x₂ = x₁), none in B; their weights 10 + 20 + 22 = 52 = 4 · 13 fill
   the defect identity.
3. **Exhaustive count.** `search/negaskew/qp` lists all 64 skew first blocks and 2303 other blocks, joins them into
   36,000 level-1 instances and completes 4800 of them: 13,056 solutions with B, C, D normalised to u(0) = 1
   (835,584 ordered quadruples), in 11 s on one core. This one is the 7361st printed.
4. **Independent checks.** A second program (a K-level pattern join written separately) found a different solution
   of the same class, and an independent reimplementation of both levels gave the same 13,056 solutions as a set.

**Reproduction.**

| | |
|---|---|
| Code | `search/negaskew` (`qp.c`, `qshape.py`, `plugcert.py`; build commands in its README) |
| Command | `python3 qshape.py '[13]' 23 '[3]' shapes/p299.txt`, `./qp shapes/p299.txt l2c=1 maxprint=20000`, then `python3 plugcert.py 13 23 "[3]" "03120\|31320\|23102 00313\|10322\|32100 01233\|20130\|22312 01022\|30132\|32310"` |
| Randomness and caps | none: both levels are exhaustive and deterministic |
| Result | 13,056 solutions; `plugcert.py` writes `out/cert_2392.json`, identical to `certificate.json` |

## Status

**New.** An ordinary Hadamard matrix of order 2392 is known (299 is not in Cati–Pasechnik Table 4); the skew property
is what is new.

- Cati–Pasechnik Table 5 gives 299(4): the smallest known skew order with odd part 299 was 16 · 299 = 4784. 299 is
  undecided in the Đoković–Seberry census (order 4 · 299), which this matrix does not settle.
- SageMath's `skew_hadamard_matrix` (development source of 21 October 2026, read, not run) has no construction of
  order 2392, and lists 1196 in `unknown_skew_hadamard`.
- We found no classical construction:
  - Paley I needs 2391 = 3 · 797; Szekeres and Momihara–Xiang 2(qᵉ + 1) need 1195 = 5 · 239; Wallis–Whiteman and
    Momihara–Xiang 2ᵘ(qᵉ + 1) need 597 = 3 · 199; none is a prime power.
  - Spence 4(q² + q + 1) would need q² + q + 1 = 598, and Turyn (p − 1)ᵘ + 1 a perfect power 2391: neither exists.
  - Seberry's mn(n − 1) construction only fits n = 2, which needs an OD(2392; 1, 1196, 1195).
  - Doubling needs a skew-Hadamard matrix of order 1196 (census-undecided). Skew × amicable products have a factor
    that is not a Hadamard order (4 · 598, 8 · 299, 52 · 46, 92 · 26, 104 · 23, 184 · 13), except 2 · 1196, which
    is doubling.

The sources searched are listed in [NOVELTY.md](../../NOVELTY.md). Doubling gives the skew orders 4784 (Table 5's
value), 9568, …

## Files

| File | Contents |
|---|---|
| `certificate.json` | the four ±1 sequences, the quaternary orbit exponents and H, the sequence sums |
| `sequences.txt` | the four sequences, one per line |
| `H2392.csv.gz` | the matrix |

## Verification

```
python3 verify/verify.py skew/2392
gzip -dc skew/2392/H2392.csv.gz | ./check_matrix
```

SHA-256 of `H2392.csv`: `72e8e896f9f6c9bce173e2334d4c13c8d9ead49314a3c0f2c4a299d07411a376`
