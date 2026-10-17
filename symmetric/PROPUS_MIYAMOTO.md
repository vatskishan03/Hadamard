# Symmetric Hadamard matrices of order 4q from Miyamoto's lemma

**Theorem.** Let q ≡ 1 (mod 8) be a prime power. If a symmetric Hadamard matrix of order (q − 1)/2 exists, then a
symmetric Hadamard matrix of order 4q exists.

Both halves of the proof are published: Seberry and Yamada's form of Miyamoto's lemma, which gives Williamson-type
matrices W₁, W₂, W₃, W₄ of order q with W₃ = W₄ (J. Combin. Math. Combin. Comput. 7 (1990), Lemma 26 and Corollary
29(ii)), and Seberry and Balonin's propus array, into which such matrices fit (Australas. J. Combin. 69 (2017), Theorem
3; Lemma 3 of arXiv:1512.01732). Seberry and Balonin apply it in their Corollary 1 (v = q + 2) only, and write that
"appropriate Williamson type matrices may also be used to give propus-Hadamard matrices but do not pursue this avenue in
this paper". We found the theorem stated nowhere. It settles 48 of the 595 values v that the Đoković–Seberry census
(2025-05-28) lists as undecided for symmetric Hadamard matrices of order 4v: 44 with inputs of Paley and Sylvester,
three with published propus matrices of orders 116, 236 and 604, and v = 1753 with our symmetric matrix of order 876
([876](876)).

H₂ = [[1, 1], [1, −1]]. Paley I is the symmetric member P + N of Paley's amicable pair of order q′ + 1 (q′ ≡ 3 mod 4),
Paley II is Paley's symmetric matrix of order 2(q′ + 1) (q′ ≡ 1 mod 4). The propus matrices are those of Di Matteo,
Đoković and Kotsireas (order 116 = 4 · 29, Spec. Matrices 3 (2015), arXiv:1503.04226v2, §4, the first of fifteen
SDSs), Đoković (236 = 4 · 59, arXiv:2101.05429v2, §4, the first PDF for v = 59; the census credits v = 59 to Abuzin,
Balonin, Đoković and Kotsireas, 2019) and Balonin and Đoković (604 = 4 · 151, Inf. Control Syst. 2018 no. 4,
arXiv:1803.08787, Table 1); the certificates list their blocks, transcribed from the LaTeX sources, and `verify.py`
rebuilds them in the propus array of [532](532) and checks them.

| Order | q | (q − 1)/2 | M | Status |
|------:|---:|---:|---|---|
| [356](356) | 89 | 44 | Paley I, order 44 (q′ = 43) | New, but obtainable classically |
| [772](772) | 193 | 96 | H₂ ⊗ Paley I, order 48 (q′ = 47) | New, but obtainable classically |
| [932](932) | 233 | 116 | propus, order 116 (Di Matteo–Đoković–Kotsireas) | New, but obtainable classically |
| [964](964) | 241 | 120 | H₂ ⊗ Paley I, order 60 (q′ = 59) | New, but obtainable classically |
| [1028](1028) | 257 | 128 | Paley I, order 128 (q′ = 127) | New, but obtainable classically |
| [1412](1412) | 353 | 176 | H₂^⊗2 ⊗ Paley I, order 44 (q′ = 43) | New, but obtainable classically |
| [1444](1444) | 361 = 19² | 180 | Paley I, order 180 (q′ = 179) | New, but obtainable classically |
| [1604](1604) | 401 | 200 | Paley I, order 200 (q′ = 199) | New, but obtainable classically |
| [1636](1636) | 409 | 204 | Paley II, order 204 (q′ = 101) | New, but obtainable classically |
| [1732](1732) | 433 | 216 | H₂ ⊗ Paley I, order 108 (q′ = 107) | New, but obtainable classically |
| [1796](1796) | 449 | 224 | Paley I, order 224 (q′ = 223) | New, but obtainable classically |
| [1828](1828) | 457 | 228 | Paley I, order 228 (q′ = 227) | New, but obtainable classically |
| [2116](2116) | 529 = 23² | 264 | Paley I, order 264 (q′ = 263) | New, but obtainable classically |
| [2276](2276) | 569 | 284 | Paley I, order 284 (q′ = 283) | New, but obtainable classically |
| [2564](2564) | 641 | 320 | H₂^⊗2 ⊗ Paley I, order 80 (q′ = 79) | New, but obtainable classically |
| [2692](2692) | 673 | 336 | H₂ ⊗ Paley I, order 168 (q′ = 167) | New, but obtainable classically |
| [3076](3076) | 769 | 384 | Paley I, order 384 (q′ = 383) | New, but obtainable classically |
| [3524](3524) | 881 | 440 | Paley I, order 440 (q′ = 439) | New, but obtainable classically |
| [3716](3716) | 929 | 464 | Paley I, order 464 (q′ = 463) | New, but obtainable classically |
| [3844](3844) | 961 = 31² | 480 | Paley I, order 480 (q′ = 479) | New, but obtainable classically |
| [4132](4132) | 1033 | 516 | Paley II, order 516 (q′ = 257) | New, but obtainable classically |
| [4196](4196) | 1049 | 524 | Paley I, order 524 (q′ = 523) | New, but obtainable classically |
| [4388](4388) | 1097 | 548 | Paley I, order 548 (q′ = 547) | New, but obtainable classically |
| [4516](4516) | 1129 | 564 | Paley I, order 564 (q′ = 563) | New, but obtainable classically |
| [4612](4612) | 1153 | 576 | H₂^⊗3 ⊗ Paley I, order 72 (q′ = 71) | New, but obtainable classically |
| [4868](4868) | 1217 | 608 | Paley I, order 608 (q′ = 607) | New, but obtainable classically |
| [4996](4996) | 1249 | 624 | H₂ ⊗ Paley I, order 312 (q′ = 311) | New, but obtainable classically |
| [5156](5156) | 1289 | 644 | Paley I, order 644 (q′ = 643) | New, but obtainable classically |
| [5284](5284) | 1321 | 660 | Paley I, order 660 (q′ = 659) | New, but obtainable classically |
| [5476](5476) | 1369 = 37² | 684 | Paley I, order 684 (q′ = 683) | New, but obtainable classically |
| [5636](5636) | 1409 | 704 | H₂^⊗4 ⊗ Paley I, order 44 (q′ = 43) | New, but obtainable classically |
| [5956](5956) | 1489 | 744 | Paley I, order 744 (q′ = 743) | New, but obtainable classically |
| [6404](6404) | 1601 | 800 | H₂^⊗2 ⊗ Paley I, order 200 (q′ = 199) | New, but obtainable classically |
| [6788](6788) | 1697 | 848 | H₂^⊗2 ⊗ Paley I, order 212 (q′ = 211) | New, but obtainable classically |
| [7012](7012) | 1753 | 876 | our symmetric matrix of order 876 | New |
| [7108](7108) | 1777 | 888 | Paley I, order 888 (q′ = 887) | New, but obtainable classically |
| [7204](7204) | 1801 | 900 | Paley II, order 900 (q′ = 449) | New, but obtainable classically |
| [7492](7492) | 1873 | 936 | H₂ ⊗ Paley I, order 468 (q′ = 467) | New, but obtainable classically |
| [7556](7556) | 1889 | 944 | H₂^⊗2 ⊗ propus, order 236 (Đoković) | New, but obtainable classically |
| [8068](8068) | 2017 | 1008 | H₂ ⊗ Paley I, order 504 (q′ = 503) | New, but obtainable classically |
| [8324](8324) | 2081 | 1040 | Paley I, order 1040 (q′ = 1039) | New, but obtainable classically |
| [8452](8452) | 2113 | 1056 | H₂^⊗2 ⊗ Paley I, order 264 (q′ = 263) | New, but obtainable classically |
| [8516](8516) | 2129 | 1064 | Paley I, order 1064 (q′ = 1063) | New, but obtainable classically |
| [8644](8644) | 2161 | 1080 | H₂ ⊗ Paley II, order 540 (q′ = 269) | New, but obtainable classically |
| [8836](8836) | 2209 = 47² | 1104 | Paley I, order 1104 (q′ = 1103) | New, but obtainable classically |
| [9508](9508) | 2377 | 1188 | Paley I, order 1188 (q′ = 1187) | New, but obtainable classically |
| [9668](9668) | 2417 | 1208 | H₂ ⊗ propus, order 604 (Balonin–Đoković) | New, but obtainable classically |
| [9892](9892) | 2473 | 1236 | Paley II, order 1236 (q′ = 617) | New, but obtainable classically |

## Proof

Let n = (q − 1)/2, g a primitive element of GF(q), χ the quadratic character, s_i = g²ⁱ (the squares) and
t_j = g¹⁻²ʲ (the non-squares), i, j ∈ ℤ_n, and M a symmetric Hadamard matrix of order n. Put

```
U₁ = −[χ(s_i − s_j)],  U₂ = [χ(s_i − t_j)],  U₃ = U₄ = 0,      V₁ = I,  V₂ = 0,  V₃ = V₄ = M.
```

These satisfy the conditions of Seberry and Balonin's Theorem 3:

- (i) U₁ is a symmetric circulant, since χ(s_i − s_j) = χ(g²⁽ⁱ⁻ʲ⁾ − 1) and χ(−1) = 1. U₂ is back-circulant, since
  χ(s_i − t_j) = χ(g⁻²ʲ) χ(g²⁽ⁱ⁺ʲ⁾ − g) depends on i + j only. A product of a symmetric circulant and a back-circulant
  is back-circulant, hence symmetric, so U₁U₂ = (U₁U₂)ᵀ = U₂U₁: the U_k are pairwise amicable.
- (ii) I, 0, M, M are pairwise amicable because M is symmetric.
- (iii) U_k ± V_k is a ±1 matrix: U₁ has zero diagonal, s_i ≠ t_j, and U₃ = U₄ = 0.
- (iv) U₁ has row sum 1 and U₂ row sum 0, by the cyclotomic numbers of order 2 ((0, 0) = (q − 5)/4,
  (0, 1) = (q − 1)/4).
- (v) U₁U₁ᵀ + U₂U₂ᵀ = qI − 2J (Seberry–Yamada, from Paley's core), and V₁V₁ᵀ + V₃V₃ᵀ + V₄V₄ᵀ = I + 2nI = qI.

With S_k = U_k ⊗ J₂ + V_k ⊗ [[1, −1], [−1, 1]], W₁ = [[1, −e], [−eᵀ, S₁]] and W_k = [[1, e], [eᵀ, S_k]] (k = 2, 3,
4), the theorem gives symmetric, pairwise commuting ±1 matrices of order q with W₁² + W₂² + W₃² + W₄² = 4qI, and
W₃ = W₄ because U₃ = U₄ and V₃ = V₄. These are Seberry and Yamada's matrices of Corollary 29(ii). Then, with A = W₁,
B = W₃ and D = W₂,

```
H = [[A,  B,  B,  D],
     [B,  D, −A, −B],
     [B, −A, −D,  B],
     [D, −B,  B, −A]]
```

is symmetric, its diagonal blocks of H Hᵀ = H² are A² + 2B² + D² = 4qI, and each off-diagonal block is a sum of
commutators, e.g. AB + BD − BA − DB = 0. So H is a symmetric Hadamard matrix of order 4q. The condition q ≡ 1 (mod 8)
is only that n = (q − 1)/2 is the order of a Hadamard matrix (q = 5 gives order 20).

Two of the four matrices must be equal: no symmetric orthogonal design of order 4 has four distinct variables, since
at most three real symmetric 4 × 4 matrices square to I and anticommute pairwise (the Clifford algebra with four such
generators is M₂(ℍ), whose modules have dimension divisible by 8). So the conference-matrix form of Corollary 29(ii)
(V₃, V₄ = C ± I) and its repeated-block form (V₃ ≠ V₄, [NOTES.md](../NOTES.md) §7) do not fit the propus array.
Corollaries 27(i) and 31 also give a repeated block; every v < 2500 that they reach is decided in the census.

## Reach

Of the 595 undecided v of the census, 58 are prime powers q ≡ 1 (mod 8). A symmetric Hadamard matrix of order
(q − 1)/2 is

- built from Paley's matrices and H₂ for 44 of them;
- a published propus matrix (times H₂^⊗a) for q = 233, 1889 and 2417 ((q − 1)/2 = 4 · 29, 16 · 59 and 8 · 151; the
  census lists 29, 59 and 151 among its 12 sporadic values);
- our matrix of order 876 for q = 1753;
- unknown for the other ten, q = 809, 857, 953, 1193, 1433, 1913, 1993, 2153, 2393 and 2441: (q − 1)/2 = 4m with
  m = 101, 107, 119, 149, 179, 239, 249, 269, 299 and 305, all undecided in the census (and a Kronecker product cannot
  have order 4m with m odd).

The theorem feeds itself: a new symmetric matrix of order 4m gives q = 8m + 1 (or 2ᵃ · 8m + 1). Below 2500 this
happens only for m = 219 (q = 1753); 8 · 89 + 1, 8 · 133 + 1, 8 · 183 + 1, 8 · 193 + 1, 8 · 241 + 1 and 8 · 257 + 1
are not prime powers, and 16 · 133 + 1 = 2129 already has a Paley input.

## Status

- 47 orders are **New, but obtainable classically**: Seberry–Yamada's lemma and corollary (1990), the propus array
  (Seberry–Balonin 2015/2017) and the inputs (Paley's matrices, the propus matrices of 2015–2019) were all published,
  and only their combination is new. 7012 is **New**: its input, our symmetric matrix of order 876, is from 2026.
- The census lists all 48 values as undecided. Its Seberry–Balonin series is v = q + 2 with (q + 1)/2 a prime power or
  the order of the core of a symmetric conference matrix (Seberry–Balonin's Corollary 1). Its remark on Miyamoto
  (J. Combin. Theory Ser. A 57 (1991), Theorem 5) concerns a different claim: a Hadamard matrix of order 4q "of
  symmetric block form" from a conference matrix of order q + 1 and a symmetric Hadamard matrix of order q − 1, which
  the census says the proof does not give.
- Ordinary Hadamard matrices of all these orders are classical (Miyamoto's theorem). Doubling gives symmetric matrices
  of orders 8q, 16q, ….

The sources searched are listed in [NOVELTY.md](../NOVELTY.md).

## Reproduction and verification

```
python3 search/sym_propus_miyamoto.py targets          # the 58 census values q = 1 (mod 8) and their inputs
python3 search/sym_propus_miyamoto.py certify 89       # writes symmetric/356/certificate.json
python3 verify/verify.py symmetric/356                 # rebuilds H and checks it
```

`verify/verify.py` (construction `symmetric-propus-miyamoto`) builds M from the description in the certificate (Paley
I or II over a prime field, Kronecker products with H₂, or another certificate of this repository, checked by its
SHA-256) and checks that it is a symmetric Hadamard matrix. It builds U and V over GF(q) (from exponent and logarithm
tables for prime powers), checks that the W_k are ±1, symmetric and pairwise amicable, with Σ W_k² = 4qI and W₃ = W₄,
and checks that H = Hᵀ, H Hᵀ = 4qI and the SHA-256 of the canonical CSV. `search/sym_propus_miyamoto.py` ran these
checks for all 48 orders.
