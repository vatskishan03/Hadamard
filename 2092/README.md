# Order 2092

2092 = 4 · 523, with 523 prime.

## Construction

A supplementary difference set with parameters

```
4-{523; 261, 252, 252, 243; 485}
```

in ℤ₅₂₃. Each block is a union of cosets rH of the multiplier subgroup

```
H = ⟨468⟩ = {1, 19, 60, 94, 217, 361, 410, 462, 468} ≤ ℤ₅₂₃*
```

and `certificate.json` lists one representative r for each coset in each block.
Setting xᵢ(t) = −1 when t lies in block i and +1 otherwise gives four sequences
with row sums 1, 19, 19, 37 (1² + 19² + 19² + 37² = 2092) whose periodic
autocorrelations sum to zero at every nonzero shift. With A, B, C, D the
circulant matrices having these first rows (`sequences.txt`) and R the
back-diagonal identity matrix, the Goethals–Seidel array

```
    [  A     BR     CR     DR  ]
H = [ -BR    A     DᵀR   -CᵀR ]
    [ -CR   -DᵀR    A     BᵀR ]
    [ -DR    CᵀR  -BᵀR    A   ]
```

is a Hadamard matrix of order 2092.

## Status

New. At the time of construction (October 11, 2026), 2092 was the smallest
order without a published Hadamard matrix. The entry 523(3) is unresolved in
Cati–Pasechnik Table 4, and 523 does not satisfy the hypotheses of the Paley,
Miyamoto, Yamada, Mathon or Spence constructions. The closest of these,
Yamada's theorem with q = 521 ≡ 1 (mod 8), would require a symmetric conference
matrix of order 262, which is not known.

## Files

| File | Contents |
|---|---|
| `certificate.json` | v, multiplier, coset representatives, row sums, SHA-256 |
| `sequences.txt` | first rows of A, B, C, D (`+` = 1, `-` = −1) |
| `H2092.csv.gz` | the matrix |

SHA-256 of `H2092.csv`:
`33466cc80ff769c4682a30121fb736f0cc6075fbe6fcbf5b8f694e475e404655`
