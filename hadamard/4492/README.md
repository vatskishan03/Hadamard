# Hadamard matrix of order 4492

4492 = 4 · 1123, with 1123 = 33² + 33 + 1 prime.

## Construction

A supplementary difference set with parameters

```
4-{1123; 561, 561, 561, 528; 1088}
```

in ℤ₁₁₂₃, each block a union of cosets of the multiplier subgroup H = ⟨60⟩ of order 33. The coset representatives are in
`certificate.json`. The four ±1 sequences (`sequences.txt`) have row sums 1, 1, 1, 67, and 3 + 67² = 4492. They are
placed in the Goethals–Seidel array exactly as described for [order 2092](../2092).

Because 1123 − 1 = 33 · 34, every sequence that is constant on the cosets of H has row sum ≡ ±1 (mod 66), and the
pattern (1, 1, 1, 67) is admissible.

## Search

The search space is small, 2³⁵ sequences, and rich in solutions. The search found one within seconds of matching.

## Status

**New.** The entry 1123(3) is unresolved in Cati–Pasechnik Table 4. None of the standard recursive constructions
applies:

- Paley: 4491 = 3² · 499 and 2245 = 5 · 449 are not prime powers.
- Miyamoto and Yamada: 1123 ≡ 3 (mod 4), and q = 1121 = 19 · 59 is not a prime power.
- Mathon: 2245 is not of the form q²(q + 2).
- Spence's construction of order 4(q² + q + 1) needs q to be a prime power, and q = 33 is not.

The sources searched are listed in [NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | v, multiplier, coset representatives, row sums, SHA-256 |
| `sequences.txt` | first rows of A, B, C, D |
| `H4492.csv.gz` | the matrix |

## Verification

```
python3 verify/verify.py 4492
gzip -dc hadamard/4492/H4492.csv.gz | ./check_matrix
```

SHA-256 of `H4492.csv`: `fe6d71ffa75342d420ff0e8811a911ab8d98e6503ad53d4e12c700cc94fe7add`