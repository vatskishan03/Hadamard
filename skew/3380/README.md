# Skew-Hadamard matrix of order 3380

3380 = 4 · 845, with 845 = 5 · 13².

## Construction

A skew-Hadamard matrix is a Hadamard matrix H with H + Hᵀ = 2I.

**Group.** The difference family lives on G = ℤ₅ × GF(169)⁺ = ℤ₅ × ℤ₁₃ × ℤ₁₃, a non-cyclic group of order 845.
Elements are triples (a, c₀, c₁) in mixed-radix order (last coordinate fastest); (c₀, c₁) are the coordinates of
c₀ + c₁w in GF(169) = F₁₃[x]/(x² + x + 2), w the class of x (a primitive element; `attack/cp_field.py`).

**Symmetry: a cyclotomic product of index 8.** All four sequences are invariant under H = ⟨w⁸⟩, the subgroup of
order 21 of GF(169)*, acting by multiplication on the GF(169) factor and trivially on ℤ₅. On coordinates its
generator is the matrix [[1, 0, 0], [0, 12, 7], [0, 3, 9]] stored in `certificate.json` under `automorphisms`.
Since 168/8 = 21 is odd, −1 = w⁸⁴ is not in H (it lies in the class w⁴H), so x ↦ −x pairs the orbits and a skew first
block is possible. A block is given at each position a ∈ ℤ₅ by its value at (a, 0) and its values on the eight
classes wʲH: 45 H-orbits per block.

**Difference family.** x₁, …, x₄ are ±1 functions on G whose periodic autocorrelations sum to zero at every
nonzero shift. The first is of skew type: x₁(0) = 1 and x₁(−g) = −x₁(g). The row sums are 1, 43, −3, −39, and
1² + 43² + 3² + 39² = 3380.

**The matrix.** With group-developed blocks X[g, h] = x(h − g) and R[g, h] = 1 when g + h = 0, the
Goethals–Seidel array

```
    [  A     BR     CR     DR  ]
H = [ -BR    A     DᵀR   -CᵀR ]
    [ -CR   -DᵀR    A     BᵀR ]
    [ -DR    CᵀR  -BᵀR    A   ]
```

is a Hadamard matrix of order 3380 with H + Hᵀ = 2I.

## Search

1. **Structure** (`attack/NOTES.md` §16, §21). The cyclotomic product ℤ₅ × GF(169) with classes of index e = 8
   (class size F = 21) and no multiplier on ℤ₅, found by `attack/cp_hunt.py --fibre-only --estimate-only` over the
   census values with a prime-power fibre (§22).
2. **Stage 1: the compressions to ℤ₅.** The compression is y(a) = f(a) + 21 s(a), and Σ_k Σ_a y_k(a)² = 4 · 845
   becomes Σ over the 20 positions of (21 s² + 2 f s) = 160. Stage 1 is sampled by `cp_lowe.stage1_sample` (seed 1);
   the second of its three quadruples (class 1) is y = (1, 1, −1, 1, −1), (−1, 41, 1, 1, 1), (−1, −1, −1, 1, −1),
   (1, 1, −1, −41, 1).
3. **Lists.** The blocks of this class have 78,400 (skew block), 1,344,560,000, 1,680,700,000 and 1,344,560,000
   members. The Gaussian estimate of the class (20 coordinates, rank 19) is 10^4.54 solutions.
4. **Low-energy k-tree** (`attack/cp_lowe.py --field 13,2 --stream`, `attack/cp_enum.c`, `attack/kt4.c`). `cp_enum.c`
   keeps the 50,002,409, 25,105,989 and 25,135,175 members of the non-skew blocks whose autocorrelation vectors are
   closest to their conditional means given that the four vectors sum to zero; the Gaussian estimate inside these
   lists is 10^1.48. `kt4.c` joins them over 3000 split values of the merge classes 0, 1, 2, 3, 12 (predicted share of
   the solutions in the lists: 0.09). The 26th split value holds four solutions; the run stopped after recording three
   (each checked for the autocorrelation identity, the skew first block and the invariance under H).

## Reproduction

| | |
|---|---|
| Code | commit `c29fd15` (`attack/cp_lowe.py`, `attack/cp_field.py`, `attack/cp_enum.c`, `attack/kt4.c`) |
| Command | `python3 attack/cp_lowe.py 845 --q 169 --field 13,2 --e 8 --sample-stage1 3 --config 1 --L 5e7,2.5e7,2.5e7 --nsig 3000 --threads 3 --K 3 --stream --samples 1000000` |
| Randomness and caps | `--seed 1` (default) draws the sampled stage-1 quadruples, the 10⁶ samples per block for the energy statistics and the 4·10⁶ pairs per side that rank the split values; `--L 5e7,2.5e7,2.5e7` members kept per non-skew block; `--K 3` stops after three solutions |
| Determinism | the lists are the same sets on every run, but their order follows the threaded enumeration and the split values are ranked from samples drawn by position, so the ranking can differ between runs |
| Result of that run | three verified solutions after 732 s on 3 shared cores (enumeration 485 s, 26 split values in 240 s); the certificate is the first (`attack/skew_certify_group.py 845 0`) |

The tools write hit logs, temporary files and compiled helpers to `$HADAMARD_WORK` (default `<tmp>/hadamard-work`).
The certificate does not depend on the search: `verify/verify.py` rebuilds and checks the matrix from
`certificate.json` alone, including the invariance under the stated automorphism.

## Status

To our knowledge this is the first skew-Hadamard matrix of order 3380, and its doubling 6760 is new too.

- The Đoković–Seberry census of skew orders 4v for odd v < 2500 (PDF of 2025-05-28) lists 845 among its 737 undecided
  values.
- Cati–Pasechnik Table 5 (arXiv:2411.18897, v1 and v2) gives 845(4), one of their own updates (Williamson's product
  theorem): the smallest known skew order with odd part 845 was 16 · 845 = 13520 = 676 · (19 + 1).
- An *ordinary* Hadamard matrix of order 3380 is known (Cati–Pasechnik Table 4 has no entry for 845). Only the skew
  property is new.

**Novelty checklist** (audit of 2026-10-13 by a research agent):

1. **The Đoković–Seberry census** (PDF of 2025-05-28, downloaded again; byte-identical to the cached copy)
   - It lists 845 among the undecided v ("… 841, 845, 849 …").
   - 845 is in none of its series (Paley-3 would need 3379 = 31 · 109 a prime power, Szekeres 1689 = 3 · 563) nor
     among its 41 sporadic values. No later version was found.
2. **Cati–Pasechnik** (v1 and v2): Table 5 gives 845(4), marked as their update by Williamson's construction (their
   Theorem 9). Table 3 ends at n = 299. The first skew order missing from Sage is 1252.
3. **Surveys and papers.** Seberry–Yamada 1992, Table 7.1 (p. 500): "845 6". The e-prints of 40 Đoković and
   Đoković–Kotsireas papers (including those on abelian groups, arXiv:1801.07627 and 1802.00556) and of Momihara–Xiang
   contain none of 3380, 6760, 13520 or 845.
4. **SageMath** (develop of 2026-10-11, master and passagemath; read from the source, not run):
   `skew_hadamard_matrix(n, existence=True)` gives Unknown for 3380 and 6760 and True for 13520 (Williamson, 19 and
   676).
5. **2026 deposits**: none mentions 3380 or 6760 (Alpöge et al., Schneider, JD-Jones-ASES, Schumer, Piza, Kline,
   Bedani, Karoui). Zenodo searches for skew 3380 and for Hadamard 3380, 6760 and 13520 return nothing; GitHub code
   search gives only false positives.
6. **arXiv**: the e-prints of 157 Hadamard-related papers of 2024–2026 and 94 older ones; the only match is Cati–
   Pasechnik's Table 5 entry.
7. **No classical construction gives skew 3380 or 6760**:
   - Paley I: 3379 and 6759 are not prime powers. Szekeres, Spence, Wallis–Whiteman and Momihara–Xiang (2(qᵉ + 1) and
     4(qᵉ + 1)): 1689, 3379 and 844 are not prime powers. Spence 4(q² + q + 1): no q.
   - Williamson's product n(pᵏ + 1): the cofactors 845, 169, 1690, 338 and 65 are not skew orders.
   - Seberry's mn(n − 1) works only with n = 2, which needs an OD(3380; 1, 1690, 1689) or OD(6760; 1, 3380, 3379); no
     source gives either.
   - Skew times amicable (order hg): every split needs a skew matrix of order 3380 itself. A halved product rule
     (order hg/2) would give 6760 from 676 and 20, but Seberry–Yamada state only hg, and Cati–Pasechnik's 177(4) and
     253(4) would read 3 under such a rule.
   - Doubling: 8 ∤ 3380, and 6760 needs 3380.

Not checked: Koukouvinos–Stylianou 2008 (Table 1), the Seberry–Yamada 2020 book (Table 9.2), Seberry 2017, the Handbook
(2007) and paywalled 2025–2026 journals. This is a finding of nothing, not a proof that no earlier construction exists.

Doubling (Seberry) gives skew orders 6760 (new, see above), 13520 (already known), …

## Files

| File | Contents |
|---|---|
| `certificate.json` | group, automorphism, row sums and the four sequences (`goethals-seidel-group`) |
| `sequences.txt` | the four sequences, one per line, elements in mixed-radix order |
| `H3380.csv.gz` | the matrix |

## Verification

```
python3 verify/verify.py skew/3380
gzip -dc skew/3380/H3380.csv.gz | ./check_matrix
```

SHA-256 of `H3380.csv`: `ef859d237d18bc2ce5effc3fdbda3738056b82cb4dada1adb4977be0825871ce`
