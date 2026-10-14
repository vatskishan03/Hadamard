# Symmetric Hadamard matrix of order 876

876 = 4 · 219, with 219 = 3 · 73.

## Construction

A symmetric Hadamard matrix is a Hadamard matrix H with H = Hᵀ.

**Sequences.** Three ±1 sequences A, B, D on ℤ₂₁₉ satisfy

    PAF_A(s) + 2 PAF_B(s) + PAF_D(s) = 0   for s ≠ 0,

with A symmetric (A(−s) = A(s)). Row sums 3, 17, 17: 3² + 2 · 17² + 17² = 876. All three sequences are constant on the
orbits of the multiplier subgroup H = ⟨4⟩ ≤ ℤ₂₁₉*, of order 9 (4 ≡ 1 (mod 3); 4 has order 9 modulo 73). On
ℤ₂₁₉ = ℤ₃ × ℤ₇₃, H fixes the ℤ₃ coordinate and acts on ℤ₇₃* as its subgroup of order 9 (the eighth powers), whose
8 cosets are the cyclotomic classes of index 8: H has 3 · (1 + 8) = 27 orbits, a cyclotomic product ℤ₃ × GF(73) of
index 8. Since −1 ∉ H, B and D are not forced to be symmetric; A is constant on the 14 orbits
of ⟨H, −1⟩. `certificate.json` lists, for each sequence, one representative r of each orbit rH on which the sequence
is −1.

**The matrix.** With A, B, D the circulant matrices of the sequences and R the back-diagonal identity,

```
    [ A     B     B     D  ]
H = [ Bᵀ    DR   −A    −BR ]
    [ Bᵀ   −A    −DR    BR ]
    [ Dᵀ   −BR    BR   −A  ]
```

is a Hadamard matrix of order 876. It is symmetric because A is, and because XR is symmetric for every circulant X.
This is the propus layout (A, B, B, D) with one symmetric block, as for [532](../532) (Seberry and Balonin, Australas.
J. Combin. 69 (2017); Balonin and Đoković's propus difference families, arXiv:1803.08787); the array is derived and
checked in `attack/NOTES.md` §17.

## Search

`attack/propus_chain.py` (`attack/NOTES.md` §25) searches a multiplier structure exhaustively through compressions.
For d | v the compressions of A, B, D to ℤ_d satisfy the same identity on ℤ_d and are constant on the H-orbits of ℤ_d.
Stage 1 lists every compression each block attains (a dynamic program over the orbits) and joins the three lists
exactly; stage 2 lists, for each compatible triple of compressions, all members of the three blocks (meet in the
middle) and joins them on the 13 autocorrelation classes with `attack/join4x.c`. Every hit is rebuilt as a matrix and
checked (H = Hᵀ, H Hᵀ = 4v I).

The run was the sweep of the composite census values v ≤ 1000 over every multiplier structure with at most 36 orbits
(smallest first). For v = 219 the first structure with an admissible row pattern, ⟨2⟩ (order 18, 18 orbits), has none;
in the second, ⟨4⟩, the compressions to ℤ₇₃ give 768 compatible triples (pair work 3.4 · 10⁷, against 6.9 · 10⁹ for
the direct join of the pools of 10^3.0, 10^6.9 and 10^6.9 sequences), and the 573rd triple (lists of 27, 6561 and 243
members) holds 18 solutions. Over all 768 triples the structure and row pattern have 288 solutions, in 16 triples:
the 16 cosets of H in ℤ₂₁₉* permute them, so there are about 18 up to multipliers.

The Gaussian estimate of `attack/propus_rank.py` for this structure and row pattern is 10^0.6 solutions (13 classes,
full covariance, lattice of the class values), against 288 found. (An earlier version of the estimator dropped the
lattice factor and gave 10^−8.1; see `attack/NOTES.md` §29.)

## Reproduction

| | |
|---|---|
| Code | commit `13d8d89` (`attack/propus_chain.py`, `attack/join4x.c`) |
| Command | `python3 attack/propus_chain.py 219 --sweep --mino 1 --maxo 36 --maxwork 2e10` (as run, in `chain_queue.sh` under `ulimit -v 6000000`); the structure alone: `python3 attack/propus_chain.py 219 4 --rows 3,17,17` |
| Randomness and caps | none: the hash weights are fixed (seeds 7 and 11), the search is exhaustive; `--maxwork` (not reached) caps the pair work, and the sweep stops at the first structure with a solution |
| Result | 18 solutions in the 573rd compatible triple after 4 minutes (135 s when rerun alone); the certificate is the first (`attack/sym_certify.py 219`); the rerun gives the same sequences |

The tools write hit logs, temporary files and compiled helpers to `$HADAMARD_WORK` (default `<tmp>/hadamard-work`).
The certificate does not depend on the search: `verify/verify.py` rebuilds and checks the matrix from
`certificate.json` alone (construction `propus`: the symmetry of A, the autocorrelation identity, H Hᵀ = 876 I and
H = Hᵀ).

## Status

To our knowledge this is the first symmetric Hadamard matrix of order 876.

- The Đoković–Seberry census of symmetric and skew Hadamard matrices of order 4v for odd v < 2500 (PDF of
  2025-05-28) lists 219 among its 595 undecided values for symmetric matrices.
- Ordinary and skew Hadamard matrices of order 876 are classical (Cati–Pasechnik Tables 2 and 3: "219 SDS"; skew
  since Đoković 1992); only the symmetry is new.

**Novelty checklist** (audit of 2026-10-14 by a research agent; the downloaded sources were read, never executed):

1. **The Đoković–Seberry census** ("Census of symmetric and skew Hadamard matrices of order 4v for odd v < 2500",
   invited talk at HADAMARD 2025, Sevilla; PDF of 2025-05-28, 11 pages, SHA-256
   `4b61953c0754f7e2da48a136ee8d3db4512c9dfbd96a76abebf2689194d0f04a`): "We list below the 595 odd values of
   v < 2500 for which the existence of symmetric Hadamard matrices of order 4v is undecided: 65, 89, 101, 107, 119,
   133, 149, 153, 163, 167, 179, 183, 189, 193, 209, 213, 219, …". A mechanical parse of its series (Paley I, XXSW,
   Seberry–Balonin, the Turyn C-series, Turyn, Mathon, Seberry–Whiteman) and its 12 sporadic values finds 219 in
   none of them; its 655 known and 595 undecided values are exactly the 1250 odd v < 2500. The census warns that
   "it is likely we missed some orders". No later version was found.
2. **The propus and symmetric literature** contains neither 219 nor 876 (full texts searched): Seberry–Balonin
   (arXiv:1512.01732; Australas. J. Combin. 69 (2017)) covers n < 200, and its Corollary 1 would need q = 217, not a
   prime power; Balonin–Đoković (arXiv:1803.08787) and Đoković (arXiv:2101.05429v2, 2301.02751v3) list the open
   v < 200; Đoković (arXiv:2508.17141, 2025) builds symmetric matrices of order q(q + 1) for q ≡ 3 (mod 8), and
   q(q + 1) = 876 has no solution. Also read: arXiv:1503.04226, 1708.05098, 1710.03037, 2404.14375, 1802.00556,
   1801.07627, 1301.3671, 1511.08734, 2512.16653, 2608.06920, 2206.11625, 1212.3021 and Abuzin et al. (2019);
   Balonin's online propus catalogue stops near order 400.
3. **Cati–Pasechnik** (arXiv:2411.18897v2, 2025-08-30) have no symmetric tables ("a similar programme for symmetric
   Hadamard matrices … is planned"). **SageMath** (develop, read 2026-10-14) has no general symmetric constructor; its
   symmetric routes (Paley II, symmetric conference matrices, regular symmetric matrices with constant diagonal,
   amicable pairs) need q = 437, an SRG(437, 218, 108, 109), a square order or a prime power 875, and all fail.
4. **Handbook of Combinatorial Designs (2007), Table V.1.52** covers odd n < 200 only (seen through secondary
   quotes), so it says nothing about 219.
5. **2026 deposits**: none contains a symmetric matrix of order 876. Piza (doi:10.5281/zenodo.23244254) lists
   "876 = CW(73, 3)", an ordinary Cooper–Wallis matrix; the orders of Alpöge et al. are 668, 716, 892, 1132, …
   (ordinary); Kharaghani, Mohammadian, Tayfeh-Rezaie and Zaitsev (arXiv:2605.08661v2) give near-Williamson matrices
   of orders up to 181 and no symmetric claims; Schneider, Schumer and JD-Jones-ASES have no 876. Web, Zenodo and
   GitHub searches found nothing else.
6. **No classical route gives a symmetric matrix of order 876.** Paley II needs q = 437 = 19 · 23. A symmetric
   conference matrix of order 438 cannot exist, since 437 is not a sum of two squares (Belevitch; Brouwer's table
   marks srg(437, 218, 108, 109) as impossible), which rules out every conference-matrix route at once (Paley,
   Turyn, Mathon, Seberry–Whiteman). Kronecker products and doubling would need Hadamard matrices of order 438, 219
   or 73. Muzychuk–Xiang and the Bush-type constructions need square orders; XXSW needs q = 875 = 5³ · 7;
   Seberry–Balonin's Lemma 2 needs 437 = a² + b²; two-circulant constructions need 876 = a² + b²; Miyamoto's route
   needs v ≡ 1 (mod 4).

Not checked: **Balonin and Seberry, "Computer construction conjecture for symmetric Hadamard matrices", Math.
Scientist 43(2) (2018) 137–143**, whose abstract mentions tables of symmetric Hadamard matrices of order 4n for
n < 400 (the full text was not reachable; Seberry's 2025 census, which would include it, lists 219 as undecided);
the Seberry–Yamada 2020 book; a print copy of the Handbook; MathSciNet and zbMATH; paywalled 2025–2026 journals.
This is a finding of nothing, not a proof that no earlier construction exists.

Doubling ([[H, H], [H, −H]]) gives symmetric Hadamard matrices of orders 1752, 3504, …

## Files

| File | Contents |
|---|---|
| `certificate.json` | H, the row sums and the three sequences as orbit representatives (`propus`) |
| `sequences.txt` | A, B and D, one per line |
| `H876.csv.gz` | the matrix |

## Verification

```
python3 verify/verify.py symmetric/876
gzip -dc symmetric/876/H876.csv.gz | ./check_matrix
```

SHA-256 of `H876.csv`: `e09d8d2c93e49a361fc3e3c2188f409da338c697f64e03bdcc967be8198cbaf9`
