# Hadamard matrix of order 20084

20084 = 4 · 5021, with 5021 prime and 5021 ≡ 1 (mod 4).

## Construction

Miyamoto's theorem: if q ≡ 1 (mod 4) is a prime power and a Hadamard matrix of order q − 1 exists, then a Hadamard
matrix of order 4q exists (M. Miyamoto, A construction of Hadamard matrices, J. Combin. Theory Ser. A 57 (1991)
86–108).

Here q = 5021. The inputs are

- Paley's symmetric conference matrix of order 5022: the matrix χ(y − x), x, y ∈ ℤ₅₀₂₁, with χ the quadratic
  character, bordered by a row and a column of ones;
- the Hadamard matrix of order 5020 = q − 1 in [5020](../5020), a Goethals–Seidel matrix over ℤ₁₂₅₅.

The blocks are arranged as in SageMath's `hadamard_matrix_miyamoto_construction`; `build_miyamoto` in
`verify/verify.py` spells out the layout. `certificate.json` contains q, the choice of conference matrix, the complete
certificate of the order-5020 matrix and the SHA-256 of the result.

## Status

**First explicit construction.** Handbook Table 1.53 (2007) marks n = 5021 with "." (t = 2), as already known, but we
found no construction behind that entry. The order follows from Miyamoto's theorem only now that a Hadamard matrix of
order 5020 exists (Cati–Pasechnik Table 4 lists 1255(3)).

- Paley I and II fail (20083 = 7 · 19 · 151, 10041 = 3 · 3347). A conference matrix of order 10042 cannot exist, since
  10041 is not a sum of two squares. Miyamoto needed H(5020), which the same table lists as unknown (1255 → 3).
  Yamada and Yamada–Kiyasu at q + 2 = 5021 need q = 5019 = 3 · 7 · 239 to be a prime power. No Turyn, Spence,
  T-matrix or Williamson-type construction of the needed orders is known.
- For 5021 the "." may come from 2q − 3 = 10039 being a prime ≡ 3 (mod 4): Paley I then gives H(2(q − 1)) = H(10040),
  not H(q − 1).
- Cati–Pasechnik and Seberry–Yamada (1992) stop at n ≤ 2999. SageMath returns Unknown for 20084 (and for 40168, 80336
  and 160672). No paper, 2026 deposit or Zenodo record we checked mentions 20084 or 5021.

The sources searched are listed in [NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, conference matrix, certificate of the order-5020 matrix, SHA-256 |

The matrix is not stored.

## Verification

```
python3 verify/verify.py 20084
```

This rebuilds the matrix, checks H Hᵀ = 20084 I and compares the SHA-256.

SHA-256 of `H20084.csv`: `1878bd6ce19fe688dff004275db89f7816a80a62a56716193f6891b4b7e0ae97`
