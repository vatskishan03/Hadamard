# Paley's matrix in circulant blocks, times T-matrices

Fourteen orders 4n of [FRONTIER.md](../FRONTIER.md) follow from two published results that the tables never combined.
For each we give an explicit Goethals–Seidel quadruple and the Hadamard matrix of order 4n. The same route gives 116
orders beyond Table 4 ([below](#beyond-table-4)).

| Order | n | q | s = (q + 1)/4 | t |
|------:|---:|---:|---:|---:|
| [3156](3156) | 789 | 1051 | 263 | 3 |
| [5260](5260) | 1315 | 1051 | 263 | 5 |
| [6036](6036) | 1509 | 2011 | 503 | 3 |
| [6356](6356) | 1589 | 907 | 227 | 7 |
| [6756](6756) | 1689 | 2251 | 563 | 3 |
| [7364](7364) | 1841 | 1051 | 263 | 7 |
| [7660](7660) | 1915 | 1531 | 383 | 5 |
| [8620](8620) | 2155 | 1723 | 431 | 5 |
| [10060](10060) | 2515 | 2011 | 503 | 5 |
| [10644](10644) | 2661 | 3547 | 887 | 3 |
| [10932](10932) | 2733 | 3643 | 911 | 3 |
| [11260](11260) | 2815 | 2251 | 563 | 5 |
| [11572](11572) | 2893 | 1051 | 263 | 11 |
| [11804](11804) | 2951 | 907 | 227 | 13 |

## Construction

**1. Paley's matrix as sixteen circulants.** Let q ≡ 3 (mod 8) be a prime and s = (q + 1)/4, which is odd. Paley's
skew-Hadamard matrix H = I + S of order q + 1 lives on the projective line PG(1, q), and PSL(2, q) acts on it by
signed permutations: H(gx, gy) = e(x) e(y) H(x, y). An element g of order s in a non-split torus has no fixed point,
so PG(1, q) splits into four g-orbits of length s. Choose signs E(x) along each orbit with E(gx) = c e(x) E(x); since
g^s = 1 the products of e over the four orbits are equal, and with c that product (c^s = c, s odd) the signs close up.
Then g is a plain automorphism of E H E, and listing each orbit along g makes E H E a 4 × 4 array of circulant s × s
blocks N_kj.

**2. Multiplication by T-matrices.** Let T₁, …, T₄ be T-matrices of order t: circulant, with entries 0, ±1, disjoint
supports covering ℤ_t and Σ T_k T_kᵀ = t I. Seberry and Yamada (J. Combin. Math. Combin. Comput. 7 (1990) 97–137,
Theorem 4; the 1992 survey, Theorems 3.10 and 8.1) put

    A_j = Σ_k T_k ⊗ N_kj        (j = 1, …, 4).

The A_j are ±1 matrices developed over ℤ_t × ℤ_s, and

    Σ_j A_j A_jᵀ = Σ_{k,l} T_k T_lᵀ ⊗ Σ_j N_kj N_ljᵀ = Σ_k T_k T_kᵀ ⊗ 4s I = 4st I,

because the block rows of E H E are orthogonal. No amicability of the blocks is needed. When gcd(s, t) = 1 the A_j
are circulants of length n = st, a Goethals–Seidel quadruple (Σ PAF = 0 off 0), and the Goethals–Seidel array gives
a Hadamard matrix of order 4n.

The T-matrices of orders 3, 5, 7, 11 and 13 come from base sequences BS(m + 1, m) by Turyn's construction; the
program finds the base sequences by a meet in the middle. Each certificate lists the T-matrices it uses.

**Reach.** We checked every entry of [FRONTIER.md](../FRONTIER.md) against every prime power q ≡ 3 (mod 8) with
(q + 1)/4 dividing n and every T-matrix order we know, allowing gcd(s, t) > 1 (the Goethals–Seidel array over
ℤ_t × ℤ_s). Exactly the fourteen orders above come out. The route also reaches seven further entries of Cati–Pasechnik
Table 4, which 2026 constructions have since settled: 1359 = 453 · 3 (q = 1811), 917 = 131 · 7, 1441 = 131 · 11,
1703 = 131 · 13, 2227 = 131 · 17, 2489 = 131 · 19 (q = 523) and 2335 = 467 · 5 (q = 1867). (An earlier version
listed six and missed 917; `paley_blocks.py quad 523 7` gives its quadruple.)

## Status

**New**, but obtainable classically. Cati–Pasechnik Table 4 lists each of the fourteen n as n(3): the smallest known
order 2ᵐ · n was 8n. We found no construction of these orders before ours. Both steps are published, but we found no
source that combines them:

- Step 1 appears in several forms. Ito's relative difference sets in the dicyclic group of order 8s, for every prime
  power 4s − 1, are equivalent to Hadamard matrices of order 4s whose sixteen blocks are ±A, ±B, ±C, ±D or their
  transposes, for s × s circulants A, B, C, D (B. Schmidt, Des. Codes Cryptogr. 17 (1999) 61–68, Theorems 2.2 and
  3.3). Yamada showed that Paley's matrix is equivalent to such a matrix of generalized quaternion type (M. Yamada,
  Discrete Math. 87 (1991) 187–196; J. Seberry, M. Yamada, J. Combin. Math. Combin. Comput. 13 (1993) 97–106,
  Theorem 1, where q + 1 = 4nN with N = 1 here). Balonin and Đoković computed the blocks for q = 139 (s = 35) and noted
  that they are not of Williamson type (arXiv:1508.00640, §9).
- Step 2 is Seberry–Yamada (1990), Theorem 4, for an OD(4s; 4s) made of sixteen circulant blocks.
- Seberry and Yamada (1993) multiply the generalized quaternion structure only under extra conditions (their Theorem 6
  (iv), which ordinary T-matrices do not satisfy) and with Williamson matrices, reaching 2ʳ m₁⋯mᵣ (q + 1). Schmidt
  uses the four-way split for q ≡ 3 (mod 8) only for a three-squares identity (his Remark 3.4), and the product theorems
  for Hadamard groups multiply only by Golay numbers 2ᵃ 10ᵇ 26ᶜ (Schmidt–Tan 2014, Result 1.2). Balonin and Đoković
  did not apply Cooper–Wallis because the blocks are not Williamson matrices; Theorem 4 does not need them to be.
- Williamson-type matrices of these orders s are not known to us, and Paley's quaternion-type blocks are not amicable
  ([NOTES.md](../NOTES.md) §21), so Cooper–Wallis does not apply directly.

The sources searched are listed in [NOVELTY.md](../NOVELTY.md); the audit of this route is summarised in
[NOTES.md](../NOTES.md) §22.

## Beyond Table 4

Above Cati–Pasechnik Table 4 the reference is Table 1.53 of the Handbook of Combinatorial Designs (2nd ed., 2007), which
covers odd n < 10000. After Đoković's updates (arXiv:0912.5091, arXiv:1008.2043), 691 of its entries with 3001 ≤ n ≤
9999 had no valid construction: 637 still have t ≥ 3, and 54 are recorded as known only through Seberry–Yamada's
Corollary 8.8(i), whose printed proof fails ([NOTES.md](../NOTES.md) §7) — 13 that Đoković settled in 2010, 33 of his
set Δ of 2009 and 8 that the Handbook lists with t = 2 ([below](#orders-recorded-only-through-corollary-88i)). The route
reaches 138 of the 691, with q ≡ 3 (mod 8) a prime or a prime power and T-matrices of an order we know. Among the 120
entries with t ≥ 3 or settled by Đoković in 2010 it reaches, we built 117, of orders 12052 to 39916. They are in the
[index](README.md#orders-beyond-the-tables), each with its README.

- 109 are new beyond the tables but obtainable classically: their T-matrices (orders 3 to 73) were known before 2026,
  those of order 73 from the near-normal sequences NN(36) of Đoković (2009). One of them, 32316 = 4 · 3 · 2693, is an
  order that Đoković claimed through Corollary 8.8(i), and is labelled First explicit construction.
- 6 need T-matrices of order 103, 151 or 163, known only since 2026, and are new beyond the tables: 16308, 21836,
  32012, 34556, 36844 and 39772.
- 16892 and 25132 are Known: SageMath's development version builds them by Cooper–Wallis since October 3, 2026, from
  the T-matrices of order 103 decoded from Schneider's matrix of order 2060.
- For seven of them gcd(s, t) > 1 (16300, 17532, 18036, 23724, 26604, 27036 and 37436), and the sequences live on ℤ_t
  × ℤ_s. For 25132, 36844 and 39772, q = 243 = 3⁵.
- Other routes reach some of them, and their READMEs say so. Cooper–Wallis with the T-matrices of order 131 of London
  and Kotsireas (2025) reaches 12052, 15196, 16244, 18340, 20436, 21484, 25676, 27772, 28820 and 38252 (for 18340 with
  Đoković's good matrices of order 35; we could not read that paper, which may list them), and Yang's multiplication of
  the same base sequences by 59 reaches 30916. With the Williamson-type matrices of order 47 of 2026 the T-matrices of
  order 131 also reach 24628. (18340, 28820, 30916 and 24628 were added by the closure audit of October 17,
  [NOTES.md](../NOTES.md) §28.) Cooper–Wallis with T-matrices posted in 2026 reaches 16300, 16308, 19932, 21836, 32012,
  33220, 34556, 36844 and 39772. Kiyasu's theorem reaches 15196 from our skew-Hadamard matrix of order 1900.
- 19028, 34196 and 39260 were already in this repository, by Cooper–Wallis.
- 28820 = 4 · 5 · 11 · 131 needs T-matrices of order 55, which SageMath does not store. We built them from base
  sequences BS(6, 5) by Yang's multiplication by 5 (Seberry–Yamada 1992, Theorem 5.23, with the Golay pair of length
  2; [`t-matrices/T55.txt`](../t-matrices/T55.txt)), so 28820 is among the 109. The same T-matrices give 38060
  ([below](#orders-recorded-only-through-corollary-88i)).

SageMath's development version (`hadamard_matrix(4n, existence=True)`, October 17, 2026) returns Unknown for all the
others. Each certificate is checked exactly when it is written (entries ±1, Σ PAF = 0 off 0, sequences rebuilt from
Paley's matrix by `verify/verify.py`); the dense check H Hᵀ = 4n I of `verify/verify.py` is running.

## Orders recorded only through Corollary 8.8(i)

Seberry and Yamada's Table A.1 (1992) obtains Williamson-type matrices of the primes 173, 293, 373, …, 1993 only from
their Corollary 8.8(i) (key wl, credited to Miyamoto's 1991 paper, which we could not read), and the printed proof of
that corollary does not give Williamson-type matrices ([NOTES.md](../NOTES.md) §7). Orders built by Cooper–Wallis from
those matrices are recorded as known by Seberry–Yamada's Table A.2 (and hence omitted by Cati–Pasechnik's Table 4), by
Đoković's set Δ (arXiv:0912.5091) or by Handbook Table 1.53 with t = 2, but have no other construction that we could
find. For s = 173 and 293, q = 4s − 1 = 691 and 1171 are primes ≡ 3 (mod 8), so the route reaches every multiple s · t
with T-matrices of order t; for s = 373 it reaches 9325 = 25 · 373 through q = 7459 (s = 1865, t = 5, over ℤ₅ × ℤ₁₈₆₅).
We built 26 of them; the T-matrices of orders 21 and 33 come from the Golay pairs of lengths 20 and 32 (Turyn), those of
order 43 from `search/tmat_mult.py 43 6`, and those of order 55 from Yang's multiplication (above). They are labelled
First explicit construction: a table records them as known, and the only construction we could trace is the corollary's.
Two later findings give the Williamson-type matrices themselves for s = 293 and 173 ([NOTES.md](../NOTES.md) §7): for
293, Seberry–Yamada's symmetric Williamson-type matrices of order 73 have a repeated block, and Lemma 26 in the
repeated-block form turns them into order 293 (checked); for 173, Verma's symmetric conference matrix of order 86
(arXiv:2610.08893, 6 October 2026, checked) gives them by Corollary 8.8(ii). So the orders with s = 293 were obtainable
classically, and those with s = 173 have a published route since October 6, 2026; the labels stay, as no one stated
these constructions.

| Order | n | s | q | t | Recorded by |
|------:|---:|---:|---:|---:|---|
| [2076](2076) | 519 = 3 · 173 | 173 | 691 | 3 | Seberry–Yamada Table A.2; omitted by Cati–Pasechnik Table 4 |
| [3460](3460) | 865 = 5 · 173 | 173 | 691 | 5 | Seberry–Yamada Table A.2; omitted by Cati–Pasechnik Table 4 |
| [4844](4844) | 1211 = 7 · 173 | 173 | 691 | 7 | Seberry–Yamada Table A.2; omitted by Cati–Pasechnik Table 4 |
| [5860](5860) | 1465 = 5 · 293 | 293 | 1171 | 5 | Seberry–Yamada Table A.2; omitted by Cati–Pasechnik Table 4 |
| [6228](6228) | 1557 = 3² · 173 | 173 | 691 | 9 | Seberry–Yamada Table A.2; omitted by Cati–Pasechnik Table 4 |
| [7612](7612) | 1903 = 11 · 173 | 173 | 691 | 11 | Seberry–Yamada Table A.2; omitted by Cati–Pasechnik Table 4 |
| [8204](8204) | 2051 = 7 · 293 | 293 | 1171 | 7 | Seberry–Yamada Table A.2; omitted by Cati–Pasechnik Table 4 |
| [8996](8996) | 2249 = 13 · 173 | 173 | 691 | 13 | Seberry–Yamada Table A.2; omitted by Cati–Pasechnik Table 4 |
| [12892](12892) | 3223 = 11 · 293 | 293 | 1171 | 11 | Handbook Table 1.53, t = 2 |
| [14532](14532) | 3633 = 3 · 7 · 173 | 173 | 691 | 21 | Handbook Table 1.53, t = 2 |
| [15236](15236) | 3809 = 13 · 293 | 293 | 1171 | 13 | Handbook Table 1.53, t = 2 |
| [15916](15916) | 3979 = 23 · 173 | 173 | 691 | 23 | Handbook Table 1.53, t = 2 |
| [19924](19924) | 4981 = 17 · 293 | 293 | 1171 | 17 | Đoković's set Δ (arXiv:0912.5091) |
| [20068](20068) | 5017 = 29 · 173 | 173 | 691 | 29 | Đoković's set Δ (arXiv:0912.5091) |
| [21452](21452) | 5363 = 31 · 173 | 173 | 691 | 31 | Đoković's set Δ (arXiv:0912.5091) |
| [22268](22268) | 5567 = 19 · 293 | 293 | 1171 | 19 | Đoković's set Δ (arXiv:0912.5091) |
| [22836](22836) | 5709 = 3 · 11 · 173 | 173 | 691 | 33 | Đoković's set Δ (arXiv:0912.5091) |
| [28372](28372) | 7093 = 41 · 173 | 173 | 691 | 41 | Đoković's set Δ (arXiv:0912.5091) |
| [29300](29300) | 7325 = 5² · 293 | 293 | 1171 | 25 | Đoković's set Δ (arXiv:0912.5091) |
| [29756](29756) | 7439 = 43 · 173 | 173 | 691 | 43 | Đoković's set Δ (arXiv:0912.5091) |
| [32524](32524) | 8131 = 47 · 173 | 173 | 691 | 47 | Đoković's set Δ (arXiv:0912.5091) |
| [33908](33908) | 8477 = 7² · 173 | 173 | 691 | 49 | Đoković's set Δ (arXiv:0912.5091) |
| [36332](36332) | 9083 = 31 · 293 | 293 | 1171 | 31 | Đoković's set Δ (arXiv:0912.5091) |
| [36676](36676) | 9169 = 53 · 173 | 173 | 691 | 53 | Đoković's set Δ (arXiv:0912.5091) |
| [37300](37300) | 9325 = 5² · 373 | 1865 | 7459 | 5 | Đoković's set Δ (arXiv:0912.5091) |
| [38060](38060) | 9515 = 5 · 11 · 173 | 173 | 691 | 55 | Đoković's set Δ (arXiv:0912.5091) |

The others (2611 in the range of Table 4; 3579, 3865, 3903, 4479, 4765, 5199, 5411, 5545, 5965, 5979, 6341, 6813, 8327,
8465, 8503, 8579, 8665, 9107 and 9965 above it) are not built here. 2427, 3065, 6743 and 8899 are built by Miyamoto's
matrix in bordered blocks, with a periodic Golay pair of length 404 or Williamson matrices of order 153
([MIYAMOTO_BLOCKS.md](MIYAMOTO_BLOCKS.md#orders-recorded-only-through-corollary-88i)); 3903 and 9107 would follow in the
same way from Williamson matrices of order 325, not implemented here; and 5259 is built by Corollary 8.8(ii) instead
([21036](21036)).

## Reproduction and verification

    python3 search/paley_blocks.py check 1051          # the sixteen circulant blocks for q = 1051
    python3 search/paley_blocks.py certify 1051 3      # writes hadamard/3156/certificate.json and sequences.txt
    python3 verify/verify.py 3156                      # rebuilds and checks H H^T = 3156 I

`verify/verify.py` (construction `paley-blocks`) builds Paley's matrix independently, from the symplectic form on
PG(1, q). It finds the blocks from the certificate's τ, power and c, rebuilds the four sequences, and compares them
with the certificate. It then checks Σ PAF = 0, assembles the Goethals–Seidel matrix, and checks H Hᵀ = 4n I and the
SHA-256 of its canonical CSV. The matrices are not stored.
