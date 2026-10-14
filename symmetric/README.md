# Symmetric Hadamard matrices

A Hadamard matrix H is symmetric if H = Hᵀ. The Đoković–Seberry census (2025-05-28) lists 595 odd v < 2500 for which
no symmetric Hadamard matrix of order 4v was known; three of them, 189, 473 and 1139, were settled by Đoković in
August 2025 (arXiv:2508.17141). Each directory here holds one order with its certificate, matrix and README. Status
labels are defined in [NOVELTY.md](../NOVELTY.md).

| Order | v | Construction | Census | Status |
|------:|---|---|---|---|
| [532](532) | 133 = 7 · 19 | propus quadruple over ℤ₁₃₃, H = ⟨26⟩ of order 6; row sums 11, 5, 19 | undecided | New |
| [732](732) | 183 = 3 · 61 | propus quadruple over ℤ₁₈₃, H = ⟨34⟩ of order 5 (cyclotomic product ℤ₃ × GF(61) of index 12); row sums 21, 1, 17 | undecided | New |
| [876](876) | 219 = 3 · 73 | propus quadruple over ℤ₂₁₉, H = ⟨4⟩ of order 9 (cyclotomic product ℤ₃ × GF(73) of index 8); row sums 3, 17, 17 | undecided | New |
| [4828](4828) | 1207 = 17 · 71 | m(m + 3) with m = 68: Paley's symmetric OD(68; 1, 67), x ↦ J − 2I, y ↦ the core of a symmetric Hadamard matrix of order 72 | undecided | New |
| [5700](5700) | 1425 = 19 · 75 | m(m − 1) with m = 76: a skew-propus quadruple of length 19 (row sums 5, 5, 1, invariant under the multiplier 7) gives a symmetric OD(76; 1, 75); x ↦ J, y ↦ the core of P + N | undecided | New |

## Constructions

**Propus quadruples.** Circulant ±1 matrices (A, B, B, D) of order v with A symmetric and
A Aᵀ + 2 B Bᵀ + D Dᵀ = 4v I fill a symmetric array of order 4v (Seberry and Balonin, Australas. J. Combin. 69 (2017);
Balonin and Đoković, arXiv:1803.08787). The quadruples here are invariant under a multiplier group H and were found by
exact joins (`search/propus.py`, [search/NOTES.md](../search/NOTES.md) §17), through compressions for 876
(`search/propus_chain.py`, §25 and §28), and by the row search `search/prow.c` for 732 (§30 and §31).

**Order m(m − 1).** A skew-Hadamard matrix S = I + C of order m with an involutive signed permutation P such that
PCP = −C is the same as a symmetric orthogonal design xP + yN of type (1, m − 1), with N = CP. Replacing x by J and y
by the core D of the symmetric Hadamard matrix P + N gives the symmetric Hadamard matrix P ⊗ J + N ⊗ D of order
m(m − 1) (§31, `search/sym_mm1.py`). Đoković's theorem (arXiv:2508.17141) is the case m − 1 = q ≡ 3 (mod 8) a prime
power. The design can also come from a skew-propus quadruple of any odd length n (m = 4n), found here by exhaustive
search, or from a Legendre pair with one skew and one symmetric block.

**Order m(m + 3).** Replacing x by J − 2I instead of J, with the core of a symmetric Hadamard matrix of order m + 4,
gives order m(m + 3) (§31). This is the symmetric form of Wallis's construction of Hadamard and Williamson-type
matrices of order m(m + 3) (1970, 1975).

`python3 verify/verify.py symmetric/<order>` rebuilds a matrix from its certificate and checks H Hᵀ = nI and H = Hᵀ.
