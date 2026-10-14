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

The Gaussian estimate of `attack/propus_rank.py` for this structure and row pattern is 10^−8.1 solutions
(13 classes, full covariance), against 288 found: a miss by a factor of about 10^10.6, as for 532 (10^9.6).

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

The Đoković–Seberry census (PDF of 2025-05-28, `attack/NOTES.md` §17) lists 219 among its 595 undecided values of
v for symmetric Hadamard matrices of order 4v. Ordinary and skew Hadamard matrices of order 876 are classical; only
the symmetry is at issue. The novelty audit (literature, databases, classical routes, 2026 deposits) is in
progress; until it is recorded here, treat the order as "not found in the census", not as new.

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
