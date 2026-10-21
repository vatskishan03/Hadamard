# Symmetric Hadamard matrix of order 532

532 = 4 · 133, with 133 = 7 · 19.

## Construction

A symmetric Hadamard matrix is a Hadamard matrix H with H = Hᵀ.

**Sequences.** Three ±1 sequences A, B, D on ℤ₁₃₃ satisfy

    PAF_A(s) + 2 PAF_B(s) + PAF_D(s) = 0   for s ≠ 0,

with A symmetric (A(−s) = A(s)). Row sums 11, 5, 19: 11² + 2 · 5² + 19² = 532. All three sequences are constant on the
orbits of the multiplier subgroup H = ⟨26⟩ = {1, 11, 20, 26, 87, 121} ≤ ℤ₁₃₃*, of order 6 (26 ≡ 5 (mod 7) has order 6,
26 ≡ 7 (mod 19) has order 3). Since −1 ∉ H, B and D are not forced to be symmetric; A is constant on the orbits of
⟨H, −1⟩. On ℤ₁₃₃ = ℤ₇ × ℤ₁₉, H has 26 orbits ({0}, ℤ₇* × {0}, six orbits of size 3 in {0} × ℤ₁₉*, eighteen of size 6)
and ⟨H, −1⟩ has 14. `certificate.json` lists, for each sequence, one representative r of each orbit rH on which the
sequence is −1.

**The matrix.** With A, B, D the circulant matrices of the sequences and R the back-diagonal identity,

```
    [ A     B     B     D  ]
H = [ Bᵀ    DR   −A    −BR ]
    [ Bᵀ   −A    −DR    BR ]
    [ Dᵀ   −BR    BR   −A  ]
```

is a Hadamard matrix of order 532. It is symmetric because A is, and because XR is symmetric for every circulant X.
This is the propus layout (A, B, B, D) with one symmetric block (Seberry and Balonin, Australas. J. Combin. 69 (2017);
Balonin and Đoković's propus difference families, arXiv:1803.08787); the array is derived and checked in
[search/NOTES.md](../../search/NOTES.md) §17, where it also reproduces the published order-116 family.

## Search

`search/propus.py` ([search/NOTES.md](../../search/NOTES.md) §17) searches multiplier structures exhaustively. For
every subgroup H of ℤ_v* (cyclic subgroups and products of two) whose structures have at most 30 orbits, smallest
first, and for each of three variants (A and D symmetric with B symmetric or not; only A symmetric), it lists for
every row-sum pattern r_A² + 2 r_B² + r_D² = 4v the complete pools of A, B and D and joins them exactly on the
autocorrelation classes with `search/join4x.c`. A hit is rebuilt as a matrix and checked (H = Hᵀ, H Hᵀ = 4v I).

For v = 133 the hit came in the structure H = ⟨26⟩, only A symmetric, row pattern (11, 5, 19): pools of 966 (A),
2,855,320 (B) and 2,368,236 (D) sequences, 13 autocorrelation classes, 72 solutions in 175 s. Every earlier structure
of the sweep (fewer orbits, and the other variants of this H) had no solution. The structure has 26 orbits for B and
D, beyond the cap of 24 of the first symmetric sweep (§17).

**Reproduction.**

| | |
|---|---|
| Code | commit `75a2fce` (`search/propus.py`, `search/join4x.c`) |
| Command | `python3 search/propus.py 133 --sweep --mino 1 --maxo 30 --maxpool 5e7` |
| Randomness and caps | none: the search is exhaustive and deterministic; `--maxpool 5e7` skips row patterns with a larger pool (none in this structure), and the sweep stops at the first structure with a solution |
| Result | 72 solutions of the row pattern (11, 5, 19); the certificate is the first (`search/sym_certify.py 133`) |

## Status

**New** as a symmetric Hadamard matrix. Ordinary and skew Hadamard matrices of order 532 are classical
(Cati–Pasechnik Table 2: "133 CW(7)"; Table 3: "133 SDS"); the symmetry is what is new.

- The Đoković–Seberry census lists 133 among its 595 undecided v for symmetric matrices. 133 is in none of its eight
  symmetric series nor among its 12 sporadic values.
- The propus literature lists 133 as open: Seberry–Balonin (arXiv:1512.01732; Australas. J. Combin. 69 (2017)),
  Balonin–Đoković (arXiv:1803.08787), Balonin–Seberry (Math. Scientist 43 (2018) 137–143, Table 2: "133 —", no
  symmetric Hadamard matrix known) and Đoković (arXiv:2101.05429v2). The Handbook (2007, Table V.1.52) lists 133(3)
  among the odd n < 200 for which no symmetric Hadamard matrix of order 4n is known.
- SageMath has no general constructor for symmetric matrices, and its symmetric routes (regular symmetric matrices
  with constant diagonal, conference matrices of order 266, amicable pairs with 531 a prime power) do not apply.
- Kharaghani, Mohammadian, Tayfeh-Rezaie and Zaitsev (arXiv:2605.08661v2) give near-Williamson matrices of order 133,
  hence an ordinary matrix of order 532. No two of their blocks share an autocorrelation, so they do not fill the
  propus array.
- No classical route gives a symmetric matrix of order 532:
  - Paley II needs q = 265 = 5 · 53 to be a prime power.
  - A conference matrix of order 266 is unknown: it needs an SRG(265, 132, 65, 66), open in Brouwer's table.
  - Kronecker products and doubling would need Hadamard matrices of order 266 or 133, which cannot exist.
  - The families of Seberry–Whiteman, Turyn, Mathon, Muzychuk–Xiang, Seberry–Balonin (q = 131 ≡ 3 (mod 4)) and
    Đoković's q(q + 1) fail by arithmetic.
  - Miyamoto would need a conference matrix of order 134, and 133 is not a sum of two squares.
  - Two-circulant constructions need 532 = a² + b², which is impossible.

Doubling ([[H, H], [H, −H]]) gives symmetric Hadamard matrices of orders 1064, 2128, …

The sources searched are listed in [NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | H, the row sums and the three sequences as orbit representatives (`propus`) |
| `sequences.txt` | A, B and D, one per line |
| `H532.csv.gz` | the matrix |

## Verification

```
python3 verify/verify.py symmetric/532
cc -O2 verify/check_matrix.c -o check_matrix
gzip -dc symmetric/532/H532.csv.gz | ./check_matrix
```

`verify/verify.py` rebuilds the matrix from `certificate.json` (construction `propus`) and checks the symmetry of A,
the autocorrelation identity, H Hᵀ = 532 I and H = Hᵀ.

SHA-256 of `H532.csv`: `b071e7f8fab6bf616e81acf98f23cf481d8380a6779b89b4ca27a71d9950ad58`
