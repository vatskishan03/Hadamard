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
`attack/NOTES.md` §17, where it also reproduces the published order-116 family.

## Search

`attack/propus.py` (`attack/NOTES.md` §17) searches multiplier structures exhaustively. For every subgroup H of ℤ_v*
(cyclic subgroups and products of two) whose structures have at most 30 orbits, smallest first, and for each of three
variants (A and D symmetric with B symmetric or not; only A symmetric), it lists for every row-sum pattern
r_A² + 2 r_B² + r_D² = 4v the complete pools of A, B and D and joins them exactly on the autocorrelation classes with
`attack/join4x.c`. A hit is rebuilt as a matrix and checked (H = Hᵀ, H Hᵀ = 4v I).

For v = 133 the hit came in the structure H = ⟨26⟩, only A symmetric, row pattern (11, 5, 19): pools of 966 (A),
2,855,320 (B) and 2,368,236 (D) sequences, 13 autocorrelation classes, 72 solutions in 175 s. Every earlier structure
of the sweep (fewer orbits, and the other variants of this H) had no solution. The structure has 26 orbits for B and
D, beyond the cap of 24 of the first symmetric sweep (`attack/NOTES.md` §17).

## Reproduction

| | |
|---|---|
| Code | commit `2281188` (`attack/propus.py`, `attack/join4x.c`) |
| Command | `python3 attack/propus.py 133 --sweep --mino 1 --maxo 30 --maxpool 5e7` |
| Randomness and caps | none: the search is exhaustive and deterministic; `--maxpool 5e7` skips row patterns with a larger pool (none in this structure), and the sweep stops at the first structure with a solution |
| Result of that run | 72 solutions of the row pattern (11, 5, 19); the certificate is the first (`attack/sym_certify.py 133`) |

The tools write hit logs, temporary files and compiled helpers to `$HADAMARD_WORK` (default `<tmp>/hadamard-work`).
The certificate does not depend on the search: `verify/verify.py` rebuilds and checks the matrix from
`certificate.json` alone (construction `propus`: the symmetry of A, the autocorrelation identity, H Hᵀ = 532 I and
H = Hᵀ).

## Status

To our knowledge this is the first symmetric Hadamard matrix of order 532.

- The Đoković–Seberry census of symmetric and skew Hadamard matrices of order 4v for odd v < 2500 (PDF of
  2025-05-28) lists 133 among its 595 undecided values for symmetric matrices.
- Ordinary and skew Hadamard matrices of order 532 are classical (Cati–Pasechnik Table 2: "133 CW(7)"; Table 3:
  "133 SDS"); only the symmetry is new.

**Novelty checklist** (audit of 2026-10-13 by a research agent):

1. **The Đoković–Seberry census** (PDF of 2025-05-28, downloaded again, unchanged): "the 595 odd values of v < 2500
   for which the existence of symmetric Hadamard matrices of order 4v is undecided: 65, 89, 101, 107, 119, 133, …".
   133 is in none of its eight symmetric series nor among its 12 sporadic values. No later version was found.
2. **The propus literature** lists 133 as open and contains no 532: Seberry–Balonin (arXiv:1512.01732; Australas. J.
   Combin. 69 (2017)): unresolved n < 200 include 133; Balonin–Đoković (arXiv:1803.08787): "all odd v < 200 except
   for 59, 65, 81, 89, 93, 101, 107, 119, 127, 133, …"; Đoković (arXiv:2101.05429v2): "reduces to v = 65, 89, 93,
   101, 107, 119, 133, …". Also read: arXiv:1503.04226, 1708.05098, 1710.03037, 1802.00556, 1801.07627, 1508.00640,
   2301.02751, 2404.14375, 2508.17141 (which settles 4 · 189 and 4 · 473) and Abuzin et al. (2019).
3. **Cati–Pasechnik** (v2 of 2025-08-30) have no symmetric tables ("a similar programme for symmetric Hadamard
   matrices … is planned"). **SageMath** (develop of 2026-10-11, read from the source) has no general symmetric
   constructor; its symmetric routes (regular symmetric matrices with constant diagonal, conference matrices of order
   266, amicable pairs with 531 a prime power) do not apply.
4. **Handbook of Combinatorial Designs (2007), Table V.1.52** (web copy): "Odd n < 200 for which a symmetric Hadamard
   matrix of order 4n is not known … 133(3)".
5. **2026 deposits**: none contains a symmetric matrix of order 532 (Alpöge et al., Schneider, JD-Jones-ASES,
   Schumer, Kline, Bedani, Karoui; Piza gives "532 = CW(7, 19)", an ordinary matrix). Kharaghani, Mohammadian,
   Tayfeh-Rezaie and Zaitsev (arXiv:2605.08661v2, September 2026) publish near-Williamson matrices of order 133, which
   give an ordinary matrix of order 532; no two of their four blocks share an autocorrelation, so they do not fill
   the propus array. Searches of Zenodo, GitHub and arXiv found nothing else.
6. **No classical route gives a symmetric matrix of order 532**: Paley II needs q = 265 = 5 · 53; a conference matrix
   of order 266 is unknown (it needs an SRG(265, 132, 65, 66), marked "?" in Brouwer's table); Kronecker products and
   doubling would need Hadamard matrices of order 266 or 133, which cannot exist; the families of Seberry–Whiteman,
   Turyn, Mathon, Muzychuk–Xiang, Seberry–Balonin (q = 131 ≡ 3 (mod 4)) and Đoković's q(q + 1) fail by arithmetic;
   Miyamoto would need a conference matrix of order 134 (133 is not a sum of two squares); two-circulant
   constructions need 532 = a² + b², which is impossible.

Not checked: the Seberry–Yamada 2020 book, a print copy of the Handbook, MathSciNet and zbMATH, and paywalled
2025–2026 journals. This is a finding of nothing, not a proof that no earlier construction exists.

Doubling ([[H, H], [H, −H]]) gives symmetric Hadamard matrices of orders 1064, 2128, …

## Files

| File | Contents |
|---|---|
| `certificate.json` | H, the row sums and the three sequences as orbit representatives (`propus`) |
| `sequences.txt` | A, B and D, one per line |
| `H532.csv.gz` | the matrix |

## Verification

```
python3 verify/verify.py symmetric/532
gzip -dc symmetric/532/H532.csv.gz | ./check_matrix
```

SHA-256 of `H532.csv`: `b071e7f8fab6bf616e81acf98f23cf481d8380a6779b89b4ca27a71d9950ad58`
