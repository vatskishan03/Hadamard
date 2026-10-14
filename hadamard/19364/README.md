# Order 19364

19364 = 4 · 4841, with 4841 = 47 · 103.

## Construction

Cooper–Wallis: T-matrices of order 103 ([`t-matrices/T103.txt`](../../t-matrices/T103.txt)) give an orthogonal design OD(412; 103, 103, 103, 103) through the
Goethals–Seidel array; substituting four Williamson-type matrices of order 47 for its variables gives a Hadamard matrix
of order 4 · 103 · 47 = 19364 (J. Cooper and J. Wallis, Bull. Austral. Math. Soc. 7 (1972)).

The Williamson-type matrices are bordered ([`williamson-type/WT47.json`](../../williamson-type/WT47.json)): four ±1
sequences on ℤ₄₆, the last three symmetric, with PAF₁ + … + PAF₄ = −4 off 0 and row sums 2, 0, 0, 0, all constant on
the orbits of ⟨45⟩ = {±1} (so all four are symmetric). `build_cooper_wallis_bordered` in
[`verify/verify.py`](../../verify/verify.py) rebuilds W₁, …, W₄, checks WₖWₗᵀ = WₗWₖᵀ and Σ WₖWₖᵀ = 188 I, assembles the
matrix in int8 and checks H Hᵀ = 19364 I.

## Search

`python3 search/bnw.py search 47 45 --pat A --sub 400000 --nsig 20000`: the four blocks have 1.35–1.41 · 10⁶ members
each (complete lists); uniform samples of 4 · 10⁵ were matched by the k-tree (`search/kt4.c`), which returned 4
solutions in 288 s. 

## Status

Probably new. The order was listed as a further candidate by this repository's audit of Handbook Table 1.53 in
[NOTES.md](../../NOTES.md) §7. Checks made here: Paley I and II need 19363 = 17² · 67 and 9681 = 3 · 7 · 461 to be prime powers; 4841 is not a
prime power (Miyamoto); 4839 = 3 · 1613 is not one either (Yamada, Kiyasu); 4841 is not of Spence's form q² + q + 1. A product route needs an object of order 103 (T-matrices of order 103 are from 2026) or T-matrices of order 47 with Williamson-type matrices of order 103 (known only from the near-Williamson matrices of 2026).
Williamson-type matrices of order 47 were known before (near-Williamson matrices exist for all odd orders up to 63,
arXiv:2605.08661, 2026), so the new step is only the combination. The caveats of the other orders above 4 · 2999
apply: Table 1.53 was read from an unofficial copy, and the claim rests on a search that found nothing.

## Files

| File | Contents |
|---|---|
| `certificate.json` | T-matrices of order 103, the four sequences on ℤ₄₆, the multiplier, SHA-256 |

The matrix is not stored.

```
python3 verify/verify.py 19364
```

rebuilds and checks it.

SHA-256 of `H19364.csv`:
`64e2a7a46411f36c225985285c50415cc153fa757b35157177f397d746abd67f`
