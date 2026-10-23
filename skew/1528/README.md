# Skew-Hadamard matrix of order 1528

1528 = 8 · 191.

## Construction

A negacyclic Goethals–Seidel array ([NEGACYCLIC.md](../NEGACYCLIC.md)): four ±1 sequences of length 382 = 2 · 191
(`sequences.txt`) whose negaperiodic autocorrelations sum to 1528 at shift 0 and to 0 elsewhere, the first palindromic
with first entry 1 (so its negacyclic matrix A satisfies A + Aᵀ = 2I). The sums of the four sequences are 6, −6, −22,
18 (negacyclic blocks have no constant row sum).

In quaternary form the four sequences are u_X on ℤ₁₉₁ with values in {1, i, −1, −i}, constant on the orbits of H =
⟨180⟩ ≤ ℤ₁₉₁*, of order 19 (11 orbits), with Σ u_X u_X* = 764 δ and u_A(−g) = i · conj(u_A(g)) for g ≠ 0, u_A(0) = −i; the squared
moduli of their sums are 761, 1, 1 and 1 (761 = 19² + 20²). `certificate.json` records the orbit values as exponents
of i (`quaternary_orbit_exponents`).

## Search

`search/negaskew`: the exact orbit-level join `qskew` over the shapes of [NEGACYCLIC.md](../NEGACYCLIC.md#search)
found the solution in 4 s on one core.

**Reproduction.**

| | |
|---|---|
| Code | `search/negaskew` (build commands in its README) |
| Command | `python3 -c "from qdrv import *; print(run(191, [180], 1, 1, 1)[1])"`, then `python3 certify.py 191 "[180]" "32032213013 03001230212 02232120000 03303111220"` |
| Randomness and caps | none: the join is exhaustive and deterministic, and stops at the first solution |
| Result | the solution above; `certify.py` writes `out/cert_1528.json`, which is `certificate.json` without its `sha256` field |

## Status

**New.** An ordinary Hadamard matrix of order 1528 is known (191 is not in Cati–Pasechnik Table 4); the skew property
is what is new.

- Cati–Pasechnik Table 5 gives 191(226): the smallest known skew order with odd part 191 was 2²²⁶ · 191, from Paley's
  construction with the prime 2²²⁶ · 191 − 1. 191 is undecided in the Đoković–Seberry census (order 4 · 191), which
  this matrix does not settle.
- SageMath's `skew_hadamard_matrix` (development source of 21 October 2026, read, not run) has no construction of
  order 1528, and lists 764 in `unknown_skew_hadamard`.
- We found no classical construction:
  - Paley I needs 1527 = 3 · 509; Szekeres and Momihara–Xiang 2(qᵉ + 1) need 763 = 7 · 109; Wallis–Whiteman and
    Momihara–Xiang 2ᵘ(qᵉ + 1) need 381 = 3 · 127 or 190; none is a prime power.
  - Spence 4(q² + q + 1) would need q² + q + 1 = 382, and Turyn (p − 1)ᵘ + 1 a perfect power 1527: neither exists.
  - Seberry's mn(n − 1) construction only fits n = 2, which needs an OD(1528; 1, 764, 763).
  - Doubling needs a skew-Hadamard matrix of order 764 (census-undecided). Skew × amicable products have a factor
    that is not a Hadamard order (4 · 382, 8 · 191), except 2 · 764, which is doubling.

See [NEGACYCLIC.md](../NEGACYCLIC.md#status) for the literature on the array. The sources searched are listed in
[NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | the four ±1 sequences, the quaternary orbit exponents and H, the sequence sums |
| `sequences.txt` | the four sequences, one per line |
| `H1528.csv.gz` | the matrix |

## Verification

```
python3 verify/verify.py skew/1528
gzip -dc skew/1528/H1528.csv.gz | ./check_matrix
```

SHA-256 of `H1528.csv`: `680854e328d9448ce6165e7cc058cb6145ddeaaa9b8da3a82953fb8a663a2729`
