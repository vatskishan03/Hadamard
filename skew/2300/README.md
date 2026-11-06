# Skew-Hadamard matrix of order 2300

2300 = 4 · 575, with 575 = 5² · 23.

## Construction

A skew-Hadamard matrix is a Hadamard matrix H with H + Hᵀ = 2I.

**Difference family.** Four ±1 sequences x₁, …, x₄ on ℤ₅₇₅ satisfy two conditions:

- each is constant on the orbits of the multiplier subgroup H = ⟨301⟩ ≤ ℤ₅₇₅*, of order 11: {1, 26, 101, 151, 301, 326,
  351, 376, 426, 476, 501}. Under ℤ₅₇₅ ≅ ℤ₂₅ × ℤ₂₃, H fixes the factor ℤ₂₅ (301 ≡ 1 mod 25) and acts on ℤ₂₃ as the
  group of non-zero squares (301 ≡ 2 mod 23): the Paley product of [search/NOTES.md](../../search/NOTES.md) §10–16 with
  m = 25 and q = 23;
- their periodic autocorrelations sum to zero at every nonzero shift.

H has 74 orbits on ℤ₅₇₅ ∖ {0}: the 24 points (a, 0) with a ≠ 0, and for each a ∈ ℤ₂₅ the squares and the non-squares
of the fibre {a} × ℤ₂₃*. `certificate.json` lists, for each sequence, one representative r of each orbit rH on which
the sequence is −1 (0 is listed when the sequence is −1 at 0).

**Skew first block.** The first sequence is of skew type: x₁(0) = 1 and x₁(−t) = −x₁(t) for t ≠ 0. This is possible
because −1 ∉ H (H has odd order), so negation pairs the orbits.

**Paley-product form.** On each fibre {a} × ℤ₂₃ a block is u(a)δ + w(a)J + z(a)χ, with χ the quadratic character of
𝔽₂₃. In the four blocks w ≠ 0 at N_w = 5 positions (0, 2, 2 and 1 of them), and u = −2w at all five, as the counting
identity (q + 1)N_w − 4N₋ = 4m, that is 24 · 5 − 4 · 5 = 100, requires.

**Row sums.** The row sums are 1, 43, −3, 21, and 1² + 43² + 3² + 21² = 2300.

**The matrix.** Let A, B, C, D be the circulant matrices with these first rows (`sequences.txt`) and R the
back-diagonal identity. The Goethals–Seidel array

```
    [  A     BR     CR     DR  ]
H = [ -BR    A     DᵀR   -CᵀR ]
    [ -CR   -DᵀR    A     BᵀR ]
    [ -DR    CᵀR  -BᵀR    A   ]
```

is then a Hadamard matrix of order 2300. Because A + Aᵀ = 2I, it satisfies H + Hᵀ = 2I (Wallis–Whiteman).

## Search

Low-energy k-tree over the Paley-product structure (`search/cp_lowe.py`; [search/NOTES.md](../../search/NOTES.md) §20–22
and §115), run by the census scheduler on a 32-core machine:

| | |
|---|---|
| Structure | cyclotomic product with q = 23, e = 2 (the Paley product), m = 25; layout `config 5`, with N_w = 5 |
| Lists | 4.2 · 10⁶, 3.2 · 10⁶ and 3.2 · 10⁶ members for the three non-skew blocks, seed 1 |
| Run | 2 November 2026, 13.6 h wall on 8 threads; one solution, verified exactly |

## Status

**New** as a skew-Hadamard matrix.

- The Đoković–Seberry census lists 575 among its 737 undecided odd v < 2500 for skew matrices.
- Cati–Pasechnik Table 5 gives 575(4): the smallest skew-Hadamard order known with odd part 575 was 16 · 575 = 9200
  (Paley I, 9199 prime). The matrix here lowers the exponent from 4 to 2.

The sources searched are listed in [NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | the four sequences as orbit representatives, H, the row sums |
| `sequences.txt` | the four sequences, one per line |
| `H2300.csv.gz` | the matrix |

## Verification

```
python3 verify/verify.py skew/2300
gzip -dc skew/2300/H2300.csv.gz | ./check_matrix
```

SHA-256 of `H2300.csv`: `341c306d57677f5265da26c8e6d92bc09467d6494bc863e28ef562e5ec93605d`
