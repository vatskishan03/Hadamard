# Order 2284

2284 = 4 · 571, with 571 prime.

## Construction

A supplementary difference set with parameters

```
4-{571; 270, 275, 275, 275; 524}
```

in ℤ₅₇₁. Each block is a union of cosets of the multiplier subgroup H = ⟨481⟩ of
order 5. The last three blocks are images of one another under the multiplier
μ = 109:

```
X₂(t) = X₁(μt),   X₃(t) = X₁(μ²t).
```

μ has order 3 modulo ⟨H, −1⟩, and the first block is also invariant under
K = ⟨H, μ⟩, of order 15. This is the structure Đoković introduced
(arXiv:2404.14375) and Schumer used for [order 3676](../3676). The row sums are
31, 21, 21, 21, and 31² + 3 · 21² = 2284. The four ±1 sequences
(`sequences.txt`) go into the Goethals–Seidel array as described for
[order 2092](../2092).

## Search

Written on the 19 classes of ⟨K, −1⟩, the condition leaves 38 free signs for the
first block and 114 for the second. About 10¹⁷ solutions are expected, but the
plain structure over ℤ₅₇₁ needs at least 114 cosets per block and is far out of
reach.

1. We enumerated all C(38, 18) ≈ 3.4 · 10¹⁰ choices of the first block in 171 s.
2. We kept the 3.1 · 10⁹ whose class values lie closest to the centre of their
   distribution, in a 39 GB hash table.
3. A tabu walk over the second block evaluated all 3245 swap neighbours of every
   point it visited. Each neighbour close enough to the centre was looked up in
   the table, about 10⁸ lookups per second on 14 cores.

The Gaussian model predicted 3.5 solutions per hour. The first one came after
33 minutes of walking.

## Status

New. The entry 571(3) is unresolved in Cati–Pasechnik Table 4 and does not
appear in Đoković's updates (arXiv:0912.5091, arXiv:1008.2043). None of the
standard recursive constructions apply:

- Paley: 2283 = 3 · 761 and 1141 = 7 · 163 are not prime powers.
- Miyamoto: 571 ≡ 3 (mod 4).
- Yamada: 569 is a prime ≡ 1 (mod 8), but the symmetric conference matrix of
  order 286 it needs cannot exist, because 285 is not a sum of two squares.
- Mathon: 1141 is not of the form q²(q + 2).

## Files

| File | Contents |
|---|---|
| `certificate.json` | v, multipliers, coset representatives, row sums, SHA-256 |
| `sequences.txt` | first rows of A, B, C, D |
| `H2284.csv.gz` | the matrix |

SHA-256 of `H2284.csv`:
`5a1a637a71cabb52b20d0a2d7d39b36647c1a98af04c3a27155afb32f9d86dbc`
