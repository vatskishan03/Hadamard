# Skew-Hadamard matrix of order 3556

3556 = 4 · 889, with 889 = 7 · 127.

## Construction

A skew-Hadamard matrix is a Hadamard matrix H with H + Hᵀ = 2I.

**Difference family.** Four ±1 sequences x₁, …, x₄ on ℤ₈₈₉ satisfy two conditions:

- each is constant on the orbits of the multiplier subgroup H = ⟨2, 25⟩ ≤ ℤ₈₈₉*, of order 63;
- their periodic autocorrelations sum to zero at every nonzero shift.

`certificate.json` lists, for each sequence, one representative r of each orbit rH on which the sequence is −1.

**Skew first block.** The first sequence is of skew type: x₁(0) = 1 and x₁(−t) = −x₁(t) for t ≠ 0. This is possible
because −1 ∉ H, so negation pairs the orbits.

**Row sums.** The row sums are 1, 5, 41, 43, and 1² + 5² + 41² + 43² = 3556.

**The matrix.** Let A, B, C, D be the circulant matrices with these first rows (`sequences.txt`) and R the
back-diagonal identity. The Goethals–Seidel array

```
    [  A     BR     CR     DR  ]
H = [ -BR    A     DᵀR   -CᵀR ]
    [ -CR   -DᵀR    A     BᵀR ]
    [ -DR    CᵀR  -BᵀR    A   ]
```

is then a Hadamard matrix of order 3556. Because A + Aᵀ = 2I, it satisfies H + Hᵀ = 2I (Wallis–Whiteman).

## Search

1. `search/skew_sweep.py` runs the exact matcher on every skew multiplier structure with at most 33 orbits and pair
   work at most 1.5·10⁹, whatever the model predicts. Here:
   - H = ⟨2, 25⟩ of order 63, with 21 orbits and 10 autocorrelation classes;
   - 1024 skew first blocks;
   - pools of 20064, 37224 and 37224 sequences for row sums 5, 41 and 43.
2. `search/skew_exhaustive.py` matches the four lists exactly (`search/find4b.c`).
3. It found the solutions on the row-sum pattern (1; 5, 41, 43), in about five minutes.

**Reproduction.**

| | |
|---|---|
| Code | commit `03c9a52` (the options used behave the same at later commits) |
| Command | `python3 search/skew_sweep.py 33 1.5e9 889   (calls skew_exhaustive.run(889, [2, 25], cap=1.5e9))` |
| Randomness and caps | exact (no sampling, no seed); pair-work cap 1.5·10⁹ |
| Result | pattern (1; 5, 41, 43): solutions found (2 verified) |
| Certificate | `search/skew_certify.py` turns a logged hit into the files in this directory |

## Status

**New** as a skew-Hadamard matrix.

- Cati–Pasechnik Table 5 lists 889(4): the smallest known skew order with odd part 889 was 2⁴ · 889 = 14224. The order
  is above 1208, the range they checked against implemented constructions.
- The Đoković–Seberry census lists 889 among its 737 undecided v.
- SageMath's `unknown_skew_hadamard` list stops at 1196, below this order.
- No infinite family reaches it: Paley I, Szekeres 2(q + 1), Wallis–Whiteman, Spence, Momihara–Xiang, the product
  theorem n(pᵏ + 1) and amicable products do not apply.

The sources searched are listed in [NOVELTY.md](../../NOVELTY.md). Doubling (Seberry) also gives a skew-Hadamard
matrix of order 7112.

## Files

| File | Contents |
|---|---|
| `certificate.json` | the four sequences as orbit representatives, H, the row sums |
| `sequences.txt` | the four sequences, one per line |
| `H3556.csv.gz` | the matrix |

## Verification

```
python3 verify/verify.py skew/3556
gzip -dc skew/3556/H3556.csv.gz | ./check_matrix
```

SHA-256 of `H3556.csv`: `a7342b20c2d9a65bc38e2ff4ed4e36a3e90fd9856eec207e7c96a30e9dad34be`
