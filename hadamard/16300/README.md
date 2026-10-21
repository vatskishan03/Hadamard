# Hadamard matrix of order 16300

16300 = 4 · 4075, with 4075 = 5² · 163.

## Construction

Paley's skew-Hadamard matrix of order 3260 (q = 3259, a prime ≡ 3 mod 8) is a 4 × 4 array of circulant blocks of
order s = 815 once its points are listed along the four orbits of a torus element of PSL(2, 3259) and rescaled by
signs. Multiplying it by T-matrices of order 5 (Seberry–Yamada 1990, Theorem 4) gives four ±1 sequences on the group
ℤ₅ × ℤ₈₁₅ of order 4075 with Σ PAF = 0 off 0, and the Goethals–Seidel array gives the Hadamard matrix. The
construction and its proof are in [PALEY_BLOCKS.md](../PALEY_BLOCKS.md).

- Torus element: x ↦ −1/(x + τ) with τ = 4; sign normalisation c = +1.
- T-matrices of order 5 (from base sequences BS(3, 2), built by the program): first rows `++000`, `00+00`, `000+-`,
  `00000`.
- Row sums of the four sequences: −51, 107, −15, −45 ((−51)² + 107² + (−15)² + (−45)² = 16300).
- Since gcd(s, t) = 5, the sequences live on ℤ_5 × ℤ_815, not on ℤ_4075. `sequences.txt` lists the entry at (i, a) in
  position i · 815 + a, and the Goethals–Seidel array uses the group-developed matrices and the reversal x ↦ −x of
  the group.

## Search

None: the matrix is determined by q and the T-matrices. `python3 search/paley_blocks.py certify 3259 5` writes the
certificate.

## Status

**New, beyond the tables, but obtainable classically** ([PALEY_BLOCKS.md](../PALEY_BLOCKS.md#status)). In Table 1.53
of the Handbook of Combinatorial Designs (2nd ed., 2007), n = 4075 has t = 3: the smallest known Hadamard matrix with
odd part 4075 had order 2³ · 4075. The order lies above the range of Cati–Pasechnik Table 4, which stops at n = 2999.

- The entry does not appear in Đoković's updates (arXiv:0912.5091, arXiv:1008.2043).
- It follows from two published theorems that we found combined nowhere: Paley's matrix of order 3260 is equivalent
  to an array of sixteen circulant blocks of order 815 (Schmidt 1999, after Ito; Yamada 1991), and Seberry and Yamada
  (1990, Theorem 4) multiply such an array by T-matrices of order 5.
- Paley's constructions do not apply: 4n − 1 = 16299 = 3² · 1811 and 2n − 1 = 8149 = 29 · 281 are not prime powers.
- Miyamoto's theorem does not apply: it needs a symmetric conference matrix of order n + 1, which cannot exist for n ≡ 3
  (mod 4). Yamada's theorem with q = n − 2 = 4073 needs a symmetric conference matrix of order 2038; 2037 = 3 · 7 · 97
  is not a sum of two squares, so no conference matrix of order 2038 exists.
- SageMath's development version (17 October 2026) does not construct it: `hadamard_matrix(16300, existence=True)`
  returns Unknown.

Since 29 August 2026 the order also follows by Cooper–Wallis from T-matrices of order 163 (JD-Jones-ASES/Hadamard-T)
and Williamson-type matrices of order 25 (Seberry–Yamada 1992, Table A.1), but we found no construction or statement
of it.

We found no construction or statement of order 16300. The sources searched are listed in
[NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, τ, power, c, the T-matrices, the four sequences and the SHA-256 of the matrix |
| `sequences.txt` | the four sequences of length 4075 |

The matrix itself is not stored; `verify/verify.py` rebuilds it.

## Verification

```
python3 verify/verify.py 16300
```

SHA-256 of `H16300.csv`: `c1e65de0f7ab324d11c506bc76742fc0a63ba3e71dd83ce4e26f61a818612e67`
