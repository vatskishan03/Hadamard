# Skew-Hadamard matrix of order 2216

2216 = 8 · 277.

## Construction

A negacyclic Goethals–Seidel array ([NEGACYCLIC.md](../NEGACYCLIC.md)): four ±1 sequences of length 554 = 2 · 277
(`sequences.txt`) whose negaperiodic autocorrelations sum to 2216 at shift 0 and to 0 elsewhere, the first palindromic
with first entry 1 (so its negacyclic matrix A satisfies A + Aᵀ = 2I). The sums of the four sequences are −28, −34,
−14, 32 (negacyclic blocks have no constant row sum).

In quaternary form the four sequences are u_X on ℤ₂₇₇ with values in {1, i, −1, −i}, constant on the orbits of H =
⟨273⟩ ≤ ℤ₂₇₇*, of order 23 (13 orbits), with Σ u_X u_X* = 1108 δ and u_A(−g) = i · conj(u_A(g)) for g ≠ 0, u_A(0) = −i; the squared
moduli of their sums are 1, 1, 1 and 1105 (1105 = 23² + 24²). `certificate.json` records the orbit values as exponents
of i (`quaternary_orbit_exponents`).

## Search

`search/negaskew`: the exact orbit-level join `qskew` over the shapes of [NEGACYCLIC.md](../NEGACYCLIC.md#search)
found the solution in about 3 to 4 minutes on one core.

**Reproduction.**

| | |
|---|---|
| Code | `search/negaskew` (build commands in its README) |
| Command | `python3 -c "from qdrv import *; print(run(277, [273], 1, 1, 1)[1])"`, then `python3 certify.py 277 "[273]" "3313220001321 0301322112300 0023212321000 0033103132201"` |
| Randomness and caps | none: the join is exhaustive and deterministic, and stops at the first solution |
| Result | the solution above; `certify.py` writes `out/cert_2216.json`, which is `certificate.json` without its `sha256` field |

## Status

**New.** An ordinary Hadamard matrix of order 2216 is known (277 is not in Cati–Pasechnik Table 4); the skew property
is what is new.

- Cati–Pasechnik Table 5 gives 277(5): the smallest known skew order with odd part 277 was 32 · 277 = 8864. 277 is
  undecided in the Đoković–Seberry census (order 4 · 277), which this matrix does not settle.
- SageMath's `skew_hadamard_matrix` (development source of 21 October 2026, read, not run) has no construction of
  order 2216, and lists 1108 in `unknown_skew_hadamard`.
- We found no classical construction:
  - Paley I needs 2215 = 5 · 443; Szekeres and Momihara–Xiang 2(qᵉ + 1) need 1107 = 3³ · 41; Wallis–Whiteman and
    Momihara–Xiang 2ᵘ(qᵉ + 1) need 553 = 7 · 79 or 276; none is a prime power.
  - Spence 4(q² + q + 1) would need q² + q + 1 = 554, and Turyn (p − 1)ᵘ + 1 a perfect power 2215: neither exists.
  - Seberry's mn(n − 1) construction only fits n = 2, which needs an OD(2216; 1, 1108, 1107).
  - Doubling needs a skew-Hadamard matrix of order 1108 (census-undecided). Skew × amicable products have a factor
    that is not a Hadamard order (4 · 554, 8 · 277), except 2 · 1108, which is doubling.

See [NEGACYCLIC.md](../NEGACYCLIC.md#status) for the literature on the array. The sources searched are listed in
[NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | the four ±1 sequences, the quaternary orbit exponents and H, the sequence sums |
| `sequences.txt` | the four sequences, one per line |
| `H2216.csv.gz` | the matrix |

## Verification

```
python3 verify/verify.py skew/2216
gzip -dc skew/2216/H2216.csv.gz | ./check_matrix
```

SHA-256 of `H2216.csv`: `c15040b04ae41217c0ff0b4e6488e5114113b7f2a8fb7334d4cc21d9ade97929`
