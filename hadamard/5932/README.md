# Hadamard matrix of order 5932

5932 = 4 · 1483, with 1483 = 38² + 38 + 1 prime.

## Construction

A supplementary difference set with parameters

```
4-{1483; 741, 741, 741, 780; 1520}
```

in ℤ₁₄₈₃, each block a union of cosets of the multiplier subgroup H = ⟨28⟩ of order 39. The coset representatives are in
`certificate.json`. The four ±1 sequences (`sequences.txt`) have row sums 1, 1, 1, −77, and 3 + 77² = 5932. They are
placed in the Goethals–Seidel array exactly as described for [order 2092](../2092).

Here 1483 − 1 = 38 · 39, so sequences constant on the cosets of H have row sum ≡ ±1 (mod 78), which admits (1, 1, 1,
−77). This is the same family as [order 4492](../4492), with the subgroup of order m + 1 instead of m.

## Status

**Concurrent.** The entry 1483(3) is unresolved in Cati–Pasechnik Table 4. Egan, Kotsireas, London and Ó Catháin
(arXiv:2610.27175, 19 October 2026) constructed a Hadamard matrix of order 5932 independently, by a Goethals–Seidel
array with a multiplier group of order 39 as here; ours was committed on 11 October 2026. We claim no priority. None of
the standard recursive constructions applies:

- Paley: 5931 = 3² · 659 and 2965 = 5 · 593 are not prime powers.
- Miyamoto: 1483 ≡ 3 (mod 4).
- Yamada with q = 1481 ≡ 1 (mod 8) would need a symmetric conference matrix of order 742. None exists, since
  741 = 3 · 13 · 19 is not a sum of two squares.
- Mathon: 2965 is not of the form q²(q + 2).
- Spence's construction of order 4(q² + q + 1) needs q to be a prime power, and q = 38 is not.

The sources searched are listed in [NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | v, multiplier, coset representatives, row sums, SHA-256 |
| `sequences.txt` | first rows of A, B, C, D |
| `H5932.csv.gz` | the matrix |

## Verification

```
python3 verify/verify.py 5932
gzip -dc hadamard/5932/H5932.csv.gz | ./check_matrix
```

SHA-256 of `H5932.csv`: `a30eaf309bf1c842f9da0b12f020a497e409f7807f03c4b3a043322cf291f940`
