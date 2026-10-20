# Hadamard matrix of order 3532

3532 = 4 · 883, with 883 = 2 · 21² + 1 prime.

## Construction

A supplementary difference set with parameters

```
4-{883; 462, 441, 441, 420; 881}
```

in ℤ₈₈₃, each block a union of cosets of the multiplier subgroup H = ⟨86⟩ of order 21. The coset representatives are in
`certificate.json`. The four ±1 sequences (`sequences.txt`) have row sums −41, 1, 1, 43, and 41² + 2 + 43² = 3532. They
are placed in the Goethals–Seidel array as described for [order 2092](../2092).

## Status

**Known.** Đoković (arXiv:1008.2043, Proposition 3.1) obtained a Hadamard matrix of order 4 · 883 from Yamada's theorem
with q = 881, using Mathon's symmetric conference matrix of order 442. Cati–Pasechnik Table 4 still lists 883(3) as
unresolved. The matrix here is an explicit Goethals–Seidel matrix of this order, found by the search of
[METHOD.md](../../METHOD.md) independently of that construction. The sources searched are listed in
[NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | v, multiplier, coset representatives, row sums, SHA-256 |
| `sequences.txt` | first rows of A, B, C, D |
| `H3532.csv.gz` | the matrix |

## Verification

```
python3 verify/verify.py 3532
gzip -dc hadamard/3532/H3532.csv.gz | ./check_matrix
```

SHA-256 of `H3532.csv`: `01a1b20a5c432463fa7b7b63a8294cdf03f75aa4bba2e16fac814bdb044845de`