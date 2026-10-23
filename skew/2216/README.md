# Skew-Hadamard matrix of order 2216

2216 = 8 · 277.

## Construction

A negacyclic Goethals–Seidel array ([NEGACYCLIC.md](../NEGACYCLIC.md)): four ±1 sequences of length 554 = 2 · 277 (`sequences.txt`)
whose negaperiodic autocorrelations sum to 2216 at shift 0 and to 0 elsewhere, the first palindromic with first entry 1 (so
its negacyclic matrix A satisfies A + Aᵀ = 2I). The sums of the four sequences are −28, −34, −14, 32 (negacyclic
blocks have no constant row sum).

In quaternary form the four sequences are u_X on ℤ_277 with values in {1, i, −1, −i}, constant on the orbits of
H = ⟨273⟩ ≤ ℤ₂₇₇*, of order 23 (13 orbits), with Σ u_X u_X* = 1108 δ and u_A(−g) = i · conj(u_A(g)); the squared moduli of their sums
are 1, 1, 1 and 1105 (1105 = 23² + 24²). `certificate.json` records the orbit values as exponents of i (`quaternary_orbit_exponents`).

## Search

`search/negaskew`: the exact orbit-level join `qskew` over the shapes of [NEGACYCLIC.md](../NEGACYCLIC.md#search)
found the solution in 3 minutes on one core.

**Reproduction.**

| | |
|---|---|
| Code | `search/negaskew` (build commands in its README) |
| Command | `python3 -c "from qdrv import *; print(run(277, [273], 1, 1, 1)[1])"`, then `python3 certify.py 277 "[273]" "3313220001321 0301322112300 0023212321000 0033103132201"` |
| Randomness and caps | none: the join is exhaustive and deterministic, and stops at the first solution |
| Result | the solution above; `certify.py` writes `out/cert_2216.json` with the sequences of `certificate.json` |

## Status

**New.** Cati–Pasechnik Table 5 gives 277(5): the smallest known skew order with odd part 277 was 32 · 277 = 8864. 277 is undecided in the Đoković–Seberry census (order 4 · 277), which
this matrix does not settle. See [NEGACYCLIC.md](../NEGACYCLIC.md#status) for the literature. The sources searched are
listed in [NOVELTY.md](../../NOVELTY.md).

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
