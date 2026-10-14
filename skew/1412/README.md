# Skew-Hadamard matrix of order 1412

1412 = 4 · 353 (prime).

## Construction

A skew-Hadamard matrix is a Hadamard matrix H with H + Hᵀ = 2I.

**Difference family.** Four ±1 sequences x₁, …, x₄ on ℤ₃₅₃ satisfy two conditions:

- each is constant on the orbits of the multiplier subgroup
  H = ⟨22⟩ ≤ ℤ₃₅₃*, of order 11: {1, 22, 58, 131, 140, 185, 187, 217, 231, 256, 337};
- their periodic autocorrelations sum to zero at every nonzero shift.

`certificate.json` lists, for each sequence, one representative r of each orbit rH on which the sequence is −1.

**Skew first block.** The first sequence is of skew type: x₁(0) = 1 and x₁(−t) = −x₁(t) for t ≠ 0. This is possible
because −1 ∉ H, so negation pairs the orbits.

**Row sums.** The row sums are 1, 21, 23, 21, and 1² + 21² + 23² + 21² = 1412.

**The matrix.** Let A, B, C, D be the circulant matrices with these first rows (`sequences.txt`) and R the
back-diagonal identity. The Goethals–Seidel array

```
    [  A     BR     CR     DR  ]
H = [ -BR    A     DᵀR   -CᵀR ]
    [ -CR   -DᵀR    A     BᵀR ]
    [ -DR    CᵀR  -BᵀR    A   ]
```

is then a Hadamard matrix of order 1412. Because A + Aᵀ = 2I, it satisfies H + Hᵀ = 2I (Wallis–Whiteman).

## Search

1. `search/skew_rank.py` ranks skew Goethals–Seidel multiplier structures with a Gaussian model, counting only the
   solutions that survive when each list is cut to a sample. For H = ⟨22⟩ (order 11) the model predicts 10⁹ solutions
   in all:
   - 33 orbits and 16 autocorrelation classes;
   - 65536 skew first blocks;
   - pools of 5.7·10⁸ sequences for each of the other three blocks, row sums (21, 21, 23).
2. `search/skew_ktree.py` (four-list k-tree matcher `search/match4.c`) matched uniform samples of 1.2·10⁷ sequences
   from each pool.
3. It found 2 solutions among the 33 most likely partial-sum values.

**Reproduction.**

| | |
|---|---|
| Code | commit `03c9a52` (the options used behave the same at later commits) |
| Command | `python3 search/skew_ktree.py 353 22 --maxpool 1.2e7 --maxtau 400 --top 6 --sols 2 --seed 1` |
| Randomness and caps | uniform samples of at most 1.2·10⁷ distinct sequences per pool (seed 1); clamp S = 6; top 6; at most 400 values of τ |
| Result | 2 solutions at τ #32 |

## Status

**New** as a skew-Hadamard matrix. Cati–Pasechnik Table 5 lists 353(4): the smallest known skew order of the form
2ᵐ · 353 was 2⁴ · 353 = 5648, so 1412 was unknown. The order is above 1208, the range Cati and Pasechnik checked
against implemented constructions.

- Đoković's papers of 1992–2024 (alone and with Kotsireas, Golubitsky or Balonin), Koukouvinos–Stylianou (Discrete
  Math. 2008, Tables 1 and 3) and the Đoković–Seberry census do not construct it. The census lists 353 among its 737
  undecided values.
- SageMath's `unknown_skew_hadamard` list stops at 1196, below this order.
- No infinite family we know of reaches it. We checked Paley I, Szekeres 2(q + 1), Wallis–Whiteman, Spence,
  Momihara–Xiang, the product theorem n(pᵏ + 1) and amicable products.

The sources searched are listed in [NOVELTY.md](../../NOVELTY.md). Doubling (Seberry) gives a skew-Hadamard matrix of
order 2824.

## Files

| File | Contents |
|---|---|
| `certificate.json` | the four sequences as orbit representatives, H, the row sums |
| `sequences.txt` | the four sequences, one per line |
| `H1412.csv.gz` | the matrix |

## Verification

```
python3 verify/verify.py skew/1412
gzip -dc skew/1412/H1412.csv.gz | ./check_matrix
```

SHA-256 of `H1412.csv`: `a3298cc6fed6b9c74ce477f5332a5a051124078423b24b052a589e771d7d253e`
