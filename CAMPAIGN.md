# Branch `claude/exciting-bell-jonrin`: results of the October 14 campaign

This branch extends the repository (based on `claude/hadamard-progress-review-t7i8ji` at 35ca068) with the results
below. Every object has a certificate that the repository's checkers rebuild and verify exactly. Novelty statements
and their caveats are in each linked README.

## Real Hadamard matrices (Cooper–Wallis products, orders above 4 · 2999)

| Order | Construction | Novelty |
|---:|---|---|
| [12596](12596) = 4 · 47 · 67 | T-matrices of order 67 × bordered Williamson-type matrices of order 47 | probably new (Handbook Table 1.53 audit, NOTES.md §7) |
| [19028](19028) = 4 · 67 · 71 | T-sequences of length 71 (Turyn-type TT(24)) × near-Williamson WT(67) | probably new (same audit) |
| [19364](19364) = 4 · 47 · 103 | T-matrices of order 103 × bordered Williamson-type matrices of order 47 | probably new (same audit) |
| [34196](34196) = 4 · 83 · 103 | T-matrices of order 103 × bordered Williamson-type matrices of order 83 | probably new (same audit) |

`python3 verify/verify.py <order>` rebuilds each one and checks H Hᵀ = nI.

## Williamson-type and T-matrices

- **Near-Williamson matrices of order 111** ([near-williamson/NW111.txt](near-williamson/NW111.txt)). To our
  knowledge they are the first Williamson-type matrices of that order: Seberry–Yamada 1992 Table A.1 has t = 3 at
  111, and 111 is not among the 2026 near-Williamson orders. They come from the order-3 twist of the 2026
  near-Williamson paper, searched exactly by `attack/wt_twist.py` and `attack/twist_join.c`.
- **Bordered Williamson-type matrices** ([williamson-type](williamson-type)) of orders 47 and 83. These orders were
  known; the construction is new to the repository, with its proof in the README.
- **T-matrices of order 131** ([t-matrices/T131.txt](t-matrices/T131.txt)). No earlier construction found: 131 is not
  in the T-matrix orders of Seberry–Yamada 1992 (§3). The branch also stores T(67), T(71) and T(101) (classical
  orders, our instances).

## Complex (quaternary) Hadamard matrices of orders 2N, N odd

| Orders | Construction |
|---|---|
| 142, 154, 238, 262, 286, 302, 310, 326, 370, 494, 534 | Szöllősi's two-symmetric array (arXiv:2603.09572, Theorem 4) with multiplier groups; exact joins, sampled lists or the exact compression chain (310) |
| 222 | from the near-Williamson matrices of order 111 |

The details are in [complex-hadamard](complex-hadamard), and `python3 verify/verify_complex.py` checks K K* = 2N I
exactly. Every order is blank in Table 11.2 of Seberry–Yamada (1992), the table of known complex Hadamard orders
that the 2026 papers update, and none comes from the later work we checked (§ Status of that README). Order 142 is
the smallest complex Hadamard order left open by the 2026 papers. After 310 the smallest open order in our accounting is 214 (N = 107, a safe prime with no
usable multiplier structure), then 322 (N = 161).

## Negative and structural results (NOTES.md §7–8)

- **Order-3 twist.** It cannot reach N ≡ 2 (mod 3): reduction to ℤ₂ forces N = x² + 3y². That excludes 71, 77,
  107, 191 (573), 263 (789), …
- **Quaternary Legendre pairs of lengths 70 and 76.** Any such pair needs a Gray pair without a common centre of
  symmetry. Centres at −1/2 are ruled out by ℓ + 1 = 71, 77 not being sums of two squares. Binary-A structures with
  multipliers at the open lengths all predict at most 10^−2.2 solutions; the exhaustive searches found none.
- **Twisted and sampled near-Williamson structures** for the other Williamson-type-unknown orders are infeasible
  (lists of 10^14 and more) or empty.
- **Calibration.** For prime N the Gaussian estimate of solutions inside sampled lists is accurate: 1445 predicted
  and 1430 found at N = 71, 10^4.0 and 12618 at 131, 178 and 208 at 151. For composite N with Paley-product groups
  it underestimates by up to 10⁴: at 77 it predicts 10^−1.7 and the join finds 2352.

## Tools added

`attack/bnw.py`, `attack/bnw_enum.c`, `attack/bnw_certify.py` (bordered Williamson-type);
`attack/wt_twist.py`, `attack/twist_join.c` (order-3 twist, exact streaming join); `attack/tmat_mult.py`
(T-matrices by multipliers); `attack/ssxx.py`, `attack/ssxx_sweep.py`, `attack/ssxx_sample.py`,
`attack/ssxx_group.py`, `attack/ssxx_certify.py`, `attack/ssxx_table.py` (Szöllősi array: cyclic, sampled,
non-cyclic); `attack/qlp_bin.py` (binary-A quaternary Legendre pairs); `attack/lp_mult.py` (Legendre pairs with
multipliers); `verify/verify_complex.py`; data tables `attack/data/wt_known_1992.txt` and
`attack/data/complex_known_1992.txt` (parsed from Seberry–Yamada 1992).
