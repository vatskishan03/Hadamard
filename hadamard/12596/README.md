# Hadamard matrix of order 12596

12596 = 4 · 3149, with 3149 = 47 · 67.

## Construction

Cooper–Wallis: T-matrices of order 67 ([`t-matrices/T67.txt`](../../t-matrices/T67.txt)) give an orthogonal design
OD(268; 67, 67, 67, 67) through the Goethals–Seidel array; substituting four Williamson-type matrices of order 47 for
its variables gives a Hadamard matrix of order 4 · 67 · 47 = 12596 (J. Cooper and J. Wallis, Bull. Austral. Math.
Soc. 7 (1972)).

The T-matrices were found by `search/tmat_mult.py 67 29` (multiplier group ⟨29⟩ of order 3, row sums −3, 0, 7, 3).
T-matrices of this order were known before, from T-sequences of length 67.

The Williamson-type matrices are bordered ([`williamson-type/WT47.json`](../../williamson-type/WT47.json)): four ±1
sequences on ℤ₄₆, the last three symmetric, with PAF₁ + … + PAF₄ = −4 off 0 and row sums 2, 0, 0, 0, all constant on
the orbits of ⟨45⟩ = {±1} (so all four are symmetric). `build_cooper_wallis_bordered` in `verify/verify.py` rebuilds
W₁, …, W₄, checks WₖWₗᵀ = WₗWₖᵀ and Σ WₖWₖᵀ = 188 I, assembles the matrix in int8 and checks H Hᵀ = 12596 I.

## Search

`python3 search/bnw.py search 47 45 --pat A --sub 400000 --nsig 20000`: the four blocks have 1.35–1.41 · 10⁶ members
each (complete lists); uniform samples of 4 · 10⁵ were matched by the k-tree (`search/kt4.c`), which returned 4
solutions in 288 s. The T-matrices of order 67 came from the four-list split of `search/tmat_mult.py` at its first
random support.

## Status

**New, beyond the tables.** The order is listed as a candidate in [NOTES.md](../../NOTES.md) §7. It lies above the
range of Cati–Pasechnik Table 4, which stops at n = 2999.

- Paley I and II need 12595 = 5 · 11 · 229 and 6297 = 3 · 2099 to be prime powers.
- 3149 is not a prime power (Miyamoto). 3147 = 3 · 1049 is not one either (Yamada, Kiyasu). 3149 is not of Spence's
  form q² + q + 1.
- T-sequences of lengths 47 and 67 are classical (base sequences), but Williamson matrices of order 47 do not exist
  (Đoković 1993). Williamson-type matrices of orders 47 and 67 have been known only since the near-Williamson
  matrices of 2026 (arXiv:2605.08661, every odd order up to 69). So both factorisations need a 2026 ingredient, and
  what is new here is the combination.

The sources searched are listed in [NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | T-matrices of order 67, the four sequences on ℤ₄₆, the multiplier, SHA-256 |

The matrix is not stored.

## Verification

```
python3 verify/verify.py 12596
```

This rebuilds the matrix and checks it.

SHA-256 of `H12596.csv`: `406b4241ac45a64ac6f8eaf6c768dac4e76a1288c966a2bb6d04eea9dc06aeff`
