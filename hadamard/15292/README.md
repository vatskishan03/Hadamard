# Hadamard matrix of order 15292

15292 = 4 · 3823, with 3823 a prime.

## Construction

Yamada's theorem: if q ≡ 5 (mod 8) is a prime power and a skew-Hadamard matrix of order (q + 3)/2 exists, then a
Hadamard matrix of order 4(q + 2) exists (M. Yamada, Some new series of Hadamard matrices, J. Austral. Math. Soc.
Ser. A 46 (1989) 371–383, Theorem 4; Yamada credits it to Z. Kiyasu; Đoković, arXiv:1008.2043, Theorem 2.3). A
skew-Hadamard matrix H satisfies H + Hᵀ = 2I.

Here q = 3821, a prime ≡ 5 (mod 8), and the inputs are

- the circulants α and β of order 1911 from GF(3821²) = GF(3821)[w]/(w² − 2), with first rows ψ(Tr ξ⁴ᵐ) and
  ψ(Tr ξ⁴ᵐ⁺¹⁹¹¹) for the primitive element ξ = 8 + w and the quadratic character ψ; they satisfy
  ααᵀ + ββᵀ = 3821 I (Yamada's Lemma 3);
- a skew-Hadamard matrix Q of order 1912, the double [[S, S], [−Sᵀ, Sᵀ]] of a skew-Hadamard matrix S of order
  956 = 4 · 239.

S is the skew-Hadamard matrix of Đoković and Kotsireas (arXiv:1511.08734, §3, first family): a Goethals–Seidel matrix
over ℤ₂₃₉ whose first block has −1 on the nonzero squares (Paley–Todd, skew type), and whose other three blocks are
unions of orbits of H = ⟨10⟩ = {1, 10, 24, 44, 98, 100, 201} with their orbit representatives J₁, K₁, L₁, a
(239; 119, 112, 106; 158) difference family. The row sums are 1, 1, 15, 27, with 1 + 1 + 15² + 27² = 956. The data are
in [`search/data/skew956_dk.json`](../../search/data/skew956_dk.json).

Q is normalised to [[1, eᵀ], [−e, S′ + I]], and the four matrices of order 3822

```
A = P ⊗ β + J₂ ⊗ I,   B = P ⊗ α,   C = D = J₂ ⊗ S′ + P ⊗ I,   P = [[1, −1], [−1, 1]]
```

go into the bordered array of Yamada's Theorem 1. `build_yamada` in `verify/verify.py` gives the layout.

`certificate.json` contains q, the field data, the complete certificate of Q (which contains that of S) and the
SHA-256 of the result.

## Search

None. `python3 search/make_yamada.py search/data/skew956_dk.json 1 hadamard/15292` doubles S, applies Yamada's
construction, checks H Hᵀ = 15292 I and writes the certificate. (Our own search for the Paley-plus-three structure on
ℤ₂₃₉ under the subgroup of order 17, `python3 search/skew_wallis.py paley3 239 17`, is exhaustive and finds nothing,
so we used Đoković and Kotsireas's family under the subgroup of order 7.)

## Status

**New, beyond the tables, but obtainable classically.** In Table 1.53 of the Handbook of Combinatorial Designs (2nd
ed., 2007), n = 3823 has t = 3: the smallest known Hadamard matrix with odd part 3823 had order 8 · 3823. The order
lies above the range of Cati–Pasechnik Table 4, which stops at n = 2999.

- Every ingredient was published before 2026: Yamada's theorem (1989), the skew-Hadamard matrix of order 956
  (Đoković and Kotsireas, arXiv:1511.08734; the Đoković–Seberry census lists 239 among the decided skew values,
  citing them) and doubling. We found no construction or statement of order 15292.
- The entry does not appear in Đoković's updates (arXiv:0912.5091, arXiv:1008.2043).
- Among the Handbook entries with t ≥ 3 that Đoković's updates leave, it is the only one that the closure of the
  published theorems reaches without the T-matrices of order 131 ([NOTES.md](../../NOTES.md) §28).
- Paley I and II need 15291 = 3² · 1699 and 7645 = 5 · 11 · 139 to be prime powers. Miyamoto's theorem needs a
  symmetric conference matrix of order n + 1 = 3824, which cannot exist since 3824 ≡ 0 (mod 4). Since 3823 is a
  prime, no product theorem applies.
- SageMath's development version (17 October 2026) does not construct it: `hadamard_matrix(15292, existence=True)`
  returns Unknown.

The sources searched are listed in [NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, field data, certificate of Q, SHA-256 |

The matrix is not stored. Its CSV encoding takes about 0.6 GB.

## Verification

```
python3 verify/verify.py --write 15292
```

SHA-256 of `H15292.csv`: `f2bebcab7d45e8a263293e50aeb468926f974bdec7b3a650d92c6b89ce3e3b3d`
