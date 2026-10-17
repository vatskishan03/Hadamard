# Skew-Hadamard matrix of order 1036

1036 = 4 · 259, with 259 = 7 · 37.

## Construction

A skew-Hadamard matrix is a Hadamard matrix H with H + Hᵀ = 2I.

**Difference family.** Four ±1 sequences x₁, …, x₄ on ℤ₂₅₉ satisfy two conditions:

- each is constant on the orbits of the multiplier subgroup H = ⟨9⟩ ≤ ℤ₂₅₉*, of order 9: {1, 9, 16, 81, 86, 144, 211,
  232, 256};
- their periodic autocorrelations sum to zero at every nonzero shift.

`certificate.json` lists, for each sequence, one representative r of each orbit rH on which the sequence is −1.

**Skew first block.** The first sequence is of skew type: x₁(0) = 1 and x₁(−t) = −x₁(t) for t ≠ 0. This is possible
because −1 ∉ H, so negation pairs the orbits.

**Row sums.** The row sums are 1, 7, 31, 5, and 1² + 7² + 31² + 5² = 1036.

**The matrix.** Let A, B, C, D be the circulant matrices with these first rows (`sequences.txt`) and R the back-diagonal
identity. The Goethals–Seidel array

```
    [  A     BR     CR     DR  ]
H = [ -BR    A     DᵀR   -CᵀR ]
    [ -CR   -DᵀR    A     BᵀR ]
    [ -DR    CᵀR  -BᵀR    A   ]
```

is then a Hadamard matrix of order 1036. Because A + Aᵀ = 2I, it satisfies H + Hᵀ = 2I (Wallis–Whiteman).

## Search

1. `search/skew_rank.py` ranks skew Goethals–Seidel multiplier structures with a Gaussian model, counting only the
   solutions that survive when each list is cut to 4·10⁶ sequences. For H = ⟨9⟩ it predicts about 10⁹ solutions in all
   and 10⁵·⁸ in the cut lists:
   - 31 orbits and 15 autocorrelation classes;
   - 32768 skew first blocks;
   - pools of 3·10⁷ to 4·10⁷ sequences for the other three blocks.
2. `search/skew_ktree.py` (four-list k-tree matcher `search/match4.c`) matched uniform samples of 3.9·10⁶ sequences from
   each pool, row sums (7, 31, 5).
3. It found 2 solutions for the most likely partial-sum value.

**Reproduction.**

| | |
|---|---|
| Code | commit `42b3ba6` (the options used behave the same at later commits) |
| Command | `python3 search/skew_ktree.py 259 9 --maxtau 300 --sols 2 --maxpool 4e6 --seed 1` |
| Randomness and caps | uniform samples of at most 4·10⁶ distinct sequences per pool (seed 1); clamp S = 6; top 6; at most 300 values of τ |
| Result | pattern (5, 7, 31): 2 solutions at τ #0 |

## Status

**New** as a skew-Hadamard matrix. Cati–Pasechnik Table 5 lists 259(4): the smallest known skew order 2ᵐ · 259 was 4144,
so 1036 was unknown.

- 1036 is inside the range (≤ 1208) that Cati and Pasechnik checked against implemented constructions.
- The Đoković–Seberry census lists 259 among its 737 undecided v. Koukouvinos–Stylianou (2008, Tables 1 and 3) and
  Đoković's papers of 1992–2024 do not construct 1036.
- SageMath lists 1036 in `unknown_skew_hadamard`.
- No infinite family reaches it: we checked Paley I, Szekeres 2(q + 1), Wallis–Whiteman, Spence, Momihara–Xiang, the
  product theorem n(pᵏ + 1) and amicable products.

The sources searched are listed in [NOVELTY.md](../../NOVELTY.md). Doubling (Seberry) gives a skew-Hadamard matrix of
order 2072.

## Files

| File | Contents |
|---|---|
| `certificate.json` | the four sequences as orbit representatives, H, the row sums |
| `sequences.txt` | the four sequences, one per line |
| `H1036.csv.gz` | the matrix |

## Verification

```
python3 verify/verify.py skew/1036
gzip -dc skew/1036/H1036.csv.gz | ./check_matrix
```

SHA-256 of `H1036.csv`: `ab03b54030c8a75204713786ec15f3b4a05a81af66c6b01d1dcd8ba515c79169`