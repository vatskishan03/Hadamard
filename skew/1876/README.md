# Skew-Hadamard matrix of order 1876

1876 = 4 · 469, with 469 = 7 · 67.

## Construction

A skew-Hadamard matrix is a Hadamard matrix H with H + Hᵀ = 2I.

**Difference family.** Four ±1 sequences x₁, …, x₄ on ℤ₄₆₉ satisfy two conditions:

- each is constant on the orbits of the multiplier subgroup H = ⟨9⟩ ≤ ℤ₄₆₉*, of order 33: {1, 9, 15, 22, 25, 64, 81,
  92, 107, 135, 148, 149, 156, 158, 193, 198, 225, 226, 260, 263, 277, 282, 330, 344, 359, 375, 394, 403, 417, 424,
  442, 464, 466};
- their periodic autocorrelations sum to zero at every nonzero shift.

`certificate.json` lists, for each sequence, one representative r of each orbit rH on which the sequence is −1.

**Skew first block.** The first sequence is of skew type: x₁(0) = 1 and x₁(−t) = −x₁(t) for t ≠ 0. This is possible
because −1 ∉ H, so negation pairs the orbits.

**Row sums.** The row sums are 1, 1, 5, 43, and 1² + 1² + 5² + 43² = 1876.

**The matrix.** Let A, B, C, D be the circulant matrices with these first rows (`sequences.txt`) and R the
back-diagonal identity. The Goethals–Seidel array

```
    [  A     BR     CR     DR  ]
H = [ -BR    A     DᵀR   -CᵀR ]
    [ -CR   -DᵀR    A     BᵀR ]
    [ -DR    CᵀR  -BᵀR    A   ]
```

is then a Hadamard matrix of order 1876. Because A + Aᵀ = 2I, it satisfies H + Hᵀ = 2I (Wallis–Whiteman).

## Search

1. `search/skew_triage.py` ranks skew Goethals–Seidel multiplier structures with a Gaussian model. For H = ⟨9⟩ it
   predicts about 0.2 solutions. The search was run anyway because the structure is small:
   - 21 orbits and 10 autocorrelation classes;
   - 1024 skew first blocks;
   - pools of 2·10⁴ to 4·10⁴ sequences for the other three blocks.
2. `search/skew_exhaustive.py` matches the four lists exactly (`search/find4b.c`).
3. The only admissible row-sum pattern, (1; 1, 5, 43), gives 192 solutions in about a minute.

**Reproduction.**

| | |
|---|---|
| Code | commit `f6f42b2` (the options used behave the same at later commits) |
| Command | `python3 search/skew_exhaustive.py 469 9` |
| Randomness and caps | exact (no sampling, no seed); pair-work cap 4·10¹⁰ |
| Result | pattern (1; 1, 5, 43): 192 solutions |
| Certificate | `search/skew_certify.py` turns a logged hit into the files in this directory |

## Status

**New** as a skew-Hadamard matrix.

- Cati–Pasechnik Table 5 lists 469(3): the smallest known skew order with odd part 469 was 2³ · 469 = 3752. The order
  is above 1208, the range they checked against implemented constructions.
- The Đoković–Seberry census lists 469 among its 737 undecided v.
- SageMath's `unknown_skew_hadamard` list stops at 1196, below this order.
- No infinite family reaches it: Paley I, Szekeres 2(q + 1), Wallis–Whiteman, Spence, Momihara–Xiang, the product
  theorem n(pᵏ + 1) and amicable products do not apply.

The sources searched are listed in [NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | the four sequences as orbit representatives, H, the row sums |
| `sequences.txt` | the four sequences, one per line |
| `H1876.csv.gz` | the matrix |

## Verification

```
python3 verify/verify.py skew/1876
gzip -dc skew/1876/H1876.csv.gz | ./check_matrix
```

SHA-256 of `H1876.csv`: `15823c5abedf92aeea0792bc9af5dcf46032d2dd8f5957f60a9c697cee6092b5`
