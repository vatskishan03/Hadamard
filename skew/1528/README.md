# Skew-Hadamard matrix of order 1528

1528 = 8 · 191.

## Construction

A negacyclic Goethals–Seidel array ([NEGACYCLIC.md](../NEGACYCLIC.md)): four ±1 sequences of length 382 = 2 · 191 (`sequences.txt`)
whose negaperiodic autocorrelations sum to 1528 at shift 0 and to 0 elsewhere, the first palindromic with first entry 1 (so
its negacyclic matrix A satisfies A + Aᵀ = 2I). The sums of the four sequences are 6, −6, −22, 18 (negacyclic
blocks have no constant row sum).

In quaternary form the four sequences are u_X on ℤ_191 with values in {1, i, −1, −i}, constant on the orbits of
H = ⟨180⟩ ≤ ℤ₁₉₁*, of order 19 (11 orbits), with Σ u_X u_X* = 764 δ and u_A(−g) = i · conj(u_A(g)); the squared moduli of their sums
are 761, 1, 1 and 1 (761 = 19² + 20²). `certificate.json` records the orbit values as exponents of i (`quaternary_orbit_exponents`).

## Search

`search/negaskew`: the exact orbit-level join `qskew` over the shapes of [NEGACYCLIC.md](../NEGACYCLIC.md#search)
found the solution in 4 s on one core.

**Reproduction.**

| | |
|---|---|
| Code | `search/negaskew` (build commands in its README) |
| Command | `python3 -c "from qdrv import *; print(run(191, [180], 1, 1, 1)[1])"`, then `python3 certify.py 191 "[180]" "32032213013 03001230212 02232120000 03303111220"` |
| Randomness and caps | none: the join is exhaustive and deterministic, and stops at the first solution |
| Result | the solution above; `certify.py` writes `out/cert_1528.json` with the sequences of `certificate.json` |

## Status

**New.** Cati–Pasechnik Table 5 gives 191(226): the smallest known skew order with odd part 191 was 2²²⁶ · 191, from Paley's construction with the prime 2²²⁶ · 191 − 1. 191 is undecided in the Đoković–Seberry census (order 4 · 191), which
this matrix does not settle. See [NEGACYCLIC.md](../NEGACYCLIC.md#status) for the literature. The sources searched are
listed in [NOVELTY.md](../../NOVELTY.md).

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
