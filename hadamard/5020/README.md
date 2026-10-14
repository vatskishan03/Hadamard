# Hadamard matrix of order 5020

5020 = 4 · 1255, with 1255 = 5 · 251.

## Construction

**Difference family.** Four ±1 sequences x₁, …, x₄ on ℤ₁₂₅₅ satisfy two conditions:

- each is constant on the orbits of the multiplier subgroup H = ⟨276, 754⟩ ≤ ℤ₁₂₅₅*, of order 50:
  {1, 4, 16, 51, 64, 69, 91, 94, 149, 201, 204, 211, 219, 241, 249, 256, 271, 276, 314, 331, 351, 364, 374, 376, 494,
  506, 566, 571, 596, 651, 706, 721, 751, 754, 769, 804, 816, 844, 866, 876, 954, 964, 994, 996, 1009, 1024, 1029,
  1084, 1104, 1129};
- their periodic autocorrelations sum to zero at every nonzero shift.

`certificate.json` lists, for each sequence, one representative r of each orbit rH on which the sequence is −1.

**Row sums.** 1, −1, −53, −47, and 1² + 1² + 53² + 47² = 5020.

**The matrix.** With A, B, C, D the circulant matrices with these first rows (`sequences.txt`) and R the
back-diagonal identity, the Goethals–Seidel array

```
    [  A     BR     CR     DR  ]
H = [ -BR    A     DᵀR   -CᵀR ]
    [ -CR   -DᵀR    A     BᵀR ]
    [ -DR    CᵀR  -BᵀR    A   ]
```

is a Hadamard matrix of order 5020.

## The structure: a cyclotomic product with a multiplier on ℤ₅

ℤ₁₂₅₅ ≅ ℤ₅ × GF(251), and H = {h : h ≡ ±1 (mod 5), h mod 251 a 10th power}: it acts on ℤ₅ by ±1 (orbits {0},
{±1}, {±2}) and on GF(251) by the 10th powers, of order 25 (the cyclotomic classes C₀, …, C₉ of index 10). An
H-invariant sequence is given, on each of the three ℤ₅-orbits, by its value f at (a, 0) and its values x₀, …, x₉ on
(a, C₀), …, (a, C₉): 33 orbits. This is the structure of [order 3732](../3732) with one more multiplier.

**Compression identity.** The compression to ℤ₅ is y(a) = f(a) + 25 s(a), s(a) = Σⱼ xⱼ(a), and the four compressions
satisfy Σₖ Σₐ yₖ(a)² = 4 · 1255, that is

  Σ over the 20 positions of (25 s² + 2 f s) = 200.

A position with s = ±2 costs 96 (f = −sign s) or 104 (f = sign s), twice that on the orbits {±1}, {±2}. So exactly two
positions carry s ≠ 0: position 0 of two different blocks, with y = −49 and −51. Up to symmetry this leaves one
quadruple of compressions, with row sums (1, −1, −53, −47).

## Search

`search/cyclo_product.py` (plain mode, `--kgen 4`) lists the blocks: 16,003,008 members for each of the two blocks
without a spike and 13,335,840 for each spike block. It groups them by the autocorrelations of their compression to
GF(251) (level B, 5 equations), and `search/join3.c` joins group quadruples exactly on the remaining 20 equations.
The level-B groups number about 10⁶ per list, too many for a complete join, so the run kept the 8,000 largest groups
of each list (7.6 % and 9.0 % of the members). It found 32 solutions after 17,273 group quadruples. All were rebuilt
on ℤ₁₂₅₅ and checked.

**Reproduction.**

| | |
|---|---|
| Code | commit `eddf924` (`search/cyclo_product.py`, `search/join3.c`) |
| Command | `python3 search/cyclo_product.py 1255 --q 251 --e 10 --plain --kgen 4 --config 0 --K 3 --threads 3 --sub 8000` |
| Randomness and caps | none in the search; `--sub 8000` keeps the 8,000 largest level-B groups of each list (ties broken by a stable sort); which solutions appear first depends on the OpenMP schedule |
| Result | 32 solutions in 17,273 group quadruples, 349 s on 3 threads (lists 19 s, levels and grouping 280 s, join 48 s) |

The certificate is the first recorded solution (`search/gs_certify.py 1255`).

## Status

**New.** Cati–Pasechnik Table 4 lists 1255(3): the smallest known Hadamard order of the form 2ᵐ · 1255 was 10040
(Paley I, 10039 being a prime ≡ 3 (mod 4)), so 5020 was unknown. The matrix is neither symmetric nor skew; the
Đoković–Seberry census (May 2025) lists 1255 as undecided for both, at order 4 · 1255.

- Seberry–Yamada (1992, Table A.2: "1255 3 a1", Paley I at 10040) and the Handbook of Combinatorial Designs (2007,
  Table 1.53: 1255 → 3 and 251 → 3) give the same entry, and [FRONTIER.md](../../FRONTIER.md) listed 5020 as
  unresolved. Seberry–Yamada Table A.1 has no Williamson or Williamson-type matrices of order 1255, and only
  8-Williamson matrices of order 251 ("251 3 w#q").
- Đoković–Golubitsky–Kotsireas (arXiv:1301.3671) settled 4 · 251 = 1004, not 4 · 1255.
- SageMath returns Unknown for 5020, as a Hadamard matrix and as a skew-Hadamard matrix.
- No recursive construction reaches 5020:
  - Paley I and II: 5019 = 3 · 7 · 239 and 2509 = 13 · 193 are not prime powers.
  - Miyamoto: 1255 is neither a prime power nor ≡ 1 (mod 4).
  - Yamada and Kiyasu (q = 1253 = 7 · 179), Mathon, Turyn, Whiteman and Spence.
  - Cooper–Wallis with T(5) and Williamson-type matrices of order 251, with T(251) or with T(1255): none is known.
  - A symmetric conference matrix of order 2510 is not known, though not excluded, since 2509 = 50² + 3².
  - Product theorems: 5020 ≡ 4 (mod 8).
- We could not read the lists in §7 of Miyamoto (J. Combin. Theory Ser. A 57, 1991); Cati–Pasechnik could not verify
  them either. By Corollary 4 of that paper, 1255 would need a W(250, 249) or a W(1254, 1253), and neither can exist,
  since 249 and 1253 are not sums of two squares.

The sources searched are listed in [NOVELTY.md](../../NOVELTY.md). Doubling gives 10040, which was already known. 5021
is a prime ≡ 1 (mod 4), so Miyamoto's theorem gives a Hadamard matrix of order 4 · 5021 = 20084; see
[20084](../20084).

## Files

| File | Contents |
|---|---|
| `certificate.json` | group, multipliers, orbit representatives of the −1 entries, row sums, SHA-256 |
| `sequences.txt` | x₁, …, x₄ as strings of + and − |
| `H5020.csv.gz` | the matrix |

## Verification

```
python3 verify/verify.py 5020
gzip -dc hadamard/5020/H5020.csv.gz | ./check_matrix
```

SHA-256 of `H5020.csv`: `71d4cb944767f77cf49fae68b9d164bd722eb36786df41a642d7c8a97428e8f7`
