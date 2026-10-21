# Hadamard matrix of order 24124

24124 = 4 · 6031, with 6031 = 37 · 163.

## Construction

Yamada's theorem: if q ≡ 5 (mod 8) is a prime power and a skew-Hadamard matrix of order (q + 3)/2 exists, then a
Hadamard matrix of order 4(q + 2) exists (M. Yamada, Some new series of Hadamard matrices, J. Austral. Math. Soc.
Ser. A 46 (1989) 371–383, Theorem 4; Yamada credits it to Z. Kiyasu). A skew-Hadamard matrix H satisfies H + Hᵀ = 2I.

Here q = 6029, a prime, and the inputs are

- the circulants α and β of order 3015 from GF(6029²) = GF(6029)[w]/(w² − 2), with first rows ψ(Tr ξ⁴ᵐ) and
  ψ(Tr ξ⁴ᵐ⁺³⁰¹⁵) for the primitive element ξ = 1 + 4w and the quadratic character ψ; they satisfy
  ααᵀ + ββᵀ = 6029 I (Yamada's Lemma 3);
- a skew-Hadamard matrix Q of order 3016, the double [[S, S], [−Sᵀ, Sᵀ]] of a skew-Hadamard matrix S of order 1508.

S is a Goethals–Seidel matrix over ℤ₃₇₇ (377 = 13 · 29). Its four blocks are constant on the cosets of the multiplier
subgroup ⟨16⟩ of order 21, the first block is of skew type (x(0) = 1, x(−t) = −x(t)), and the row sums are
1, 27, 27, 7, with 1 + 27² + 27² + 7² = 1508. It came out of our skew-Hadamard searches.

Q is normalised to [[1, eᵀ], [−e, S′ + I]], and the four matrices of order 6030

```
A = P ⊗ β + J₂ ⊗ I,   B = P ⊗ α,   C = D = J₂ ⊗ S′ + P ⊗ I,   P = [[1, −1], [−1, 1]]
```

go into the bordered array of Yamada's Theorem 1. `build_yamada` in `verify/verify.py` gives the layout. Before
building the order-24124 matrix we checked the same code on six known cases, orders 60 to 700.

`certificate.json` contains q, the field data, the complete certificate of Q (which contains that of S) and the
SHA-256 of the result.

## Status

**New, beyond the tables.** In Table 1.53 of the Handbook of Combinatorial Designs (2nd ed., 2007), n = 6031 has
t = 3: the smallest known Hadamard matrix with odd part 6031 had order 8 · 6031. The order lies above the range of
Cati–Pasechnik Table 4, which stops at n = 2999.

- The entry does not appear in Đoković's updates (arXiv:0912.5091, arXiv:1008.2043).
- The order has followed from published objects since 29 August 2026 (the second route below), but we found no
  construction or statement of it.

The skew input, the double of [skew 1508](../../skew/1508), is new as far as we know. Cati–Pasechnik Table 5 lists
377(5), so the smallest known skew-Hadamard matrix with odd part 377 had order 2⁵ · 377. Their table is only
implemented and checked up to order 1208, so this is less certain than the main claim.

A second route reaches the same order: Cooper–Wallis with T-matrices of order 163 (ours are in
[`t-matrices/T163.txt`](../../t-matrices/T163.txt)) and Williamson matrices of order 37, which exist (Holzmann,
Kharaghani, Tayfeh-Rezaie, Des. Codes Cryptogr. 46 (2008) 343–352). JD-Jones-ASES/Hadamard-T published T-matrices of
order 163 on 29 August 2026, together with their products for Williamson orders 11 and 13 (orders 7172 and 8476).
T-matrices of order 163 may appear earlier in Lou, Xia, Xia, *On T-matrices of order 6m + 1* (2019/2022), which we
could not read.

The sources searched are listed in [NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, field data, certificate of Q, SHA-256 |

The matrix is not stored. Its CSV encoding takes 1.45 GB (113 MB gzipped).

## Verification

```
python3 verify/verify.py --write 24124
```

This rebuilds and checks the matrix and writes `hadamard/24124/H24124.csv`. On four cores it takes 45 seconds and
1.6 GB of memory.

SHA-256 of `H24124.csv`: `4e3b57c4c1c3d1a8661a15899a687318aa3632c3672cee0cd9f7257c8330e3f7`
