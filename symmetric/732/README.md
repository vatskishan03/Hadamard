# Symmetric Hadamard matrix of order 732

732 = 4 · 183, with 183 = 3 · 61.

## Construction

A symmetric Hadamard matrix is a Hadamard matrix H with H = Hᵀ.

**Sequences.** Three ±1 sequences A, B, D on ℤ₁₈₃ satisfy

    PAF_A(s) + 2 PAF_B(s) + PAF_D(s) = 0   for s ≠ 0,

with A symmetric (A(−s) = A(s)). Row sums 21, 1, 17: 21² + 2 · 1² + 17² = 732. All three sequences are constant on
the orbits of the multiplier subgroup H = ⟨34⟩ = {1, 34, 58, 70, 142} ≤ ℤ₁₈₃*, of order 5 (34 ≡ 1 (mod 3); 34 has
order 5 modulo 61). On ℤ₁₈₃ = ℤ₃ × ℤ₆₁, H fixes the ℤ₃ coordinate and acts on ℤ₆₁* as its subgroup of order 5 (the
twelfth powers), whose 12 cosets are the cyclotomic classes of index 12: H has 3 · (1 + 12) = 39 orbits, a cyclotomic
product ℤ₃ × GF(61) of index 12. Since −1 ∉ H, B and D are not forced to be symmetric; A is constant on the 20 orbits
of ⟨H, −1⟩. `certificate.json` lists, for each sequence, one representative r of each orbit rH on which the sequence
is −1.

**The matrix.** With A, B, D the circulant matrices of the sequences and R the back-diagonal identity,

```
    [ A     B     B     D  ]
H = [ Bᵀ    DR   −A    −BR ]
    [ Bᵀ   −A    −DR    BR ]
    [ Dᵀ   −BR    BR   −A  ]
```

is a Hadamard matrix of order 732. It is symmetric because A is, and because XR is symmetric for every circulant X.
This is the propus layout (A, B, B, D) with one symmetric block, as for [532](../532) and [876](../876) (Seberry and
Balonin, Australas. J. Combin. 69 (2017); Balonin and Đoković's propus difference families, arXiv:1803.08787); the
array is derived and checked in `search/NOTES.md` §17.

## Search

`search/prow.c` (`search/NOTES.md` §30) is a dedicated search for this kind of structure. In coordinates
(x, y) ∈ ℤ₃ × GF(61) an H-invariant sequence is three rows x₀, x₁, x₂ over GF(61), each constant on {0} and on the
12 cyclotomic classes (13 bits per row); the 19 autocorrelation classes are sums of row autocorrelations
(x = 0) and row cross-correlations (x ≠ 0), read from tables. The program lists the ℤ₃ compressions (row sums of the
three rows) compatible with the compressed identity, with the row sums of B and of D sorted (translations of ℤ₃ and
reversal permute the rows) and A canonical under the multipliers of GF(61), cuts B and D by power-spectrum filters
(|B̂(χ)|² ≤ 2v, |D̂(χ)|² ≤ 4v, Â(χ)² ≤ 4v), and joins exactly: D is hashed in parts by an additive key of its
x = 0 classes and every part is probed by every pair (A, B). Every hit is rebuilt on ℤ₁₈₃ and checked before it is
printed.

For 183/⟨34⟩ there are 8 compatible ℤ₃ triples. The fifth, row pattern (21, 1, 17) with row sums by rows
A (−1, 11, 11), B (−9, 1, 9), D (−1, 9, 9), has 464 sequences A, 3.8 · 10⁶ B and 4.3 · 10⁸ D after the filters
(1.75 · 10⁹ pairs (A, B)); it holds two solutions, found in the third of its 29 D parts. The run used two processes
(triples 0, 2, 4, 6 and 1, 3, 5, 7) and reached the hit after about 1.5 hours; the remaining triples were not
finished. The Gaussian estimate of `search/propus_rank.py` for the structure was 0.69 expected solution orbits
(0.06 for this row pattern; a solution orbit under the multipliers, the translations by the fixed points of H and
reversal has 432 elements).

## Reproduction

| | |
|---|---|
| Code | commit `7a9f079` (`search/prow.c`, the version that ran); reproduced with commit `a0fc82d` (`prow.c` as of `b407e96`, which changed only the part sizes of the join) |
| Command | `cc -O3 -march=native -o prow search/prow.c -lm`, then `./prow 61 12 all run 0 2 15000000` (as run, triples 0, 2, 4, 6; `ulimit -v` 6 GB); the triple alone: `./prow 61 12 all run 4 8 15000000` |
| Randomness and caps | none: the search is exhaustive and deterministic; `15000000` caps the D entries per part (it only sets the number of parts) |
| Result | two solutions (same A and B) in the triple (21, 1, 17); the rerun of that triple alone with the later code prints the same two lines and finishes the triple (29 parts, 1.75 · 10⁹ pairs) with no other; `search/hit2rec.py LOG 183 34` records them, `search/sym_certify.py 183` writes the certificate from the first |

The tools write hit logs, temporary files and compiled helpers to `$HADAMARD_WORK` (default `<tmp>/hadamard-work`).
The certificate does not depend on the search: `verify/verify.py` rebuilds and checks the matrix from
`certificate.json` alone (construction `propus`: the symmetry of A, the autocorrelation identity, H Hᵀ = 732 I and
H = Hᵀ).

## Status

To our knowledge this is the first symmetric Hadamard matrix of order 732.

- The Đoković–Seberry census of symmetric and skew Hadamard matrices of order 4v for odd v < 2500 (PDF of
  2025-05-28) lists 183 among its 595 undecided values for symmetric matrices.
- Ordinary and skew Hadamard matrices of order 732 are classical (Cati–Pasechnik Table 2: "183 CW(3)"; Table 3:
  "183 Spence(13)"); only the symmetry is new.

**Novelty checklist** (audit of 2026-10-14 by a research agent; downloaded sources were read as text, never
executed):

1. **The Đoković–Seberry census** (PDF of 2025-05-28, SHA-256 `4b61953c…f04a`, unchanged): "We list below the 595
   odd values of v < 2500 for which the existence of symmetric Hadamard matrices of order 4v is undecided: 65, 89,
   101, 107, 119, 133, 149, 153, 163, 167, 179, 183, …". 183 is in none of its symmetric series (Paley, XXSW,
   Seberry–Balonin, Turyn C, Turyn, Mathon, Seberry–Whiteman, Muzychuk–Xiang) nor among its 12 sporadic values. No
   later version was found.
2. **Every list of open orders we could read has 183**: Seberry–Balonin (arXiv:1512.01732; Australas. J. Combin.
   69 (2017)), "… 179, 183, 185, 189, 191, 193, 197"; Balonin–Đoković (arXiv:1803.08787, 2018), "… 179, 183, 189,
   191, 193"; Balonin–Seberry, Math. Scientist 43 (2018) 137–143, Table 2, "183 —" ("a long dash denotes that no
   symmetric Hadamard matrix is known"); Đoković (arXiv:2101.05429v2, 2022), "… 179, 183, 189, 193"; the Handbook
   of Combinatorial Designs (2007) Table V.1.52 as quoted there. Also read: arXiv:1503.04226, 1708.05098,
   1710.03037, 1802.00556, 2301.02751, 2404.14375, 2508.17141 (order q(q + 1); 732 is not of this form) and
   Abuzin et al. (2019).
3. **Cati–Pasechnik** (arXiv:2411.18897v2) have no symmetric tables. **SageMath** (develop, read 2026-10-14, not
   executed) has no general symmetric constructor; its symmetric routes (Paley II with q = 365, a conference matrix
   of order 366, constant-diagonal regular matrices, amicable pairs with 731 a prime power) do not apply.
4. **2026 deposits**: none contains a symmetric matrix of order 732 (Alpöge et al., Schneider, JD-Jones-ASES,
   Schumer, Kline, Bedani, Karoui; Piza lists "732 = CW(3, 61)", an ordinary matrix). Kline's unrefereed nonexistence
   claim for Turyn pairs with q = 365 concerns quadruples with A − D = 2I and B = C symmetric, which this one is not.
5. **Classical routes.** Paley II needs q = 365 = 5 · 73; XXSW needs 731 = 17 · 43; Mathon, Seberry–Whiteman,
   Turyn, Muzychuk–Xiang and Seberry–Balonin fail by arithmetic; Kronecker products, doubling, constant-diagonal and
   two-circulant constructions are impossible for 732; q(q + 1) and m(m − 1) have no integer solution. The one route
   not ruled out is a symmetric conference matrix of order 366 (365 = 13² + 14²), which is not known: Brouwer's
   table lists srg(365, 182, 90, 91) with "?".

Not checked: a print copy of the Handbook, the Seberry–Yamada 2020 book, MathSciNet and zbMATH, and paywalled
2025–2026 journals. This is a finding of nothing, not a proof that no earlier construction exists.

Doubling ([[H, H], [H, −H]]) gives symmetric Hadamard matrices of orders 1464, 2928, …

## Files

| File | Contents |
|---|---|
| `certificate.json` | H, the row sums and the three sequences as orbit representatives (`propus`) |
| `sequences.txt` | A, B and D, one per line |
| `H732.csv.gz` | the matrix |

## Verification

```
python3 verify/verify.py symmetric/732
gzip -dc symmetric/732/H732.csv.gz | ./check_matrix
```

SHA-256 of `H732.csv`: `7900ec94d9a949e986e2c8676641b93f1a4b04de8984f7b4c91b99e0aaf78abf`
