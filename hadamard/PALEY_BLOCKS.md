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
ℤ_t × ℤ_s). Exactly the fourteen orders above come out. The route also reaches six further entries of Cati–Pasechnik
Table 4, which 2026 constructions have since settled: 1359 = 453 · 3 (q = 1811), 1441 = 131 · 11, 1703 = 131 · 13,
2227 = 131 · 17, 2489 = 131 · 19 (q = 523) and 2335 = 467 · 5 (q = 1867).

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
covers odd n < 10000. After Đoković's updates (arXiv:0912.5091, arXiv:1008.2043), 671 of its entries with 3001 ≤ n ≤
9999 had no valid construction: 637 still have t ≥ 3, and 34 were settled only through Seberry–Yamada's Corollary
8.8(i), whose printed proof fails ([NOTES.md](../NOTES.md) §7): 13 by Đoković in 2010, 17 of his set Δ of 2009 and 4
that the Handbook lists with t = 2. One of the 17, 5259, is now settled through part (ii) of the corollary
([21036](21036)). The route reaches 120 of the 671, with q ≡ 3 (mod 8) a prime or a prime power and T-matrices of an
order we know. We built 116, of orders 12052 to 39916. They are in the [index](README.md#orders-beyond-the-tables), each
with its README.

- 108 are new beyond the tables but obtainable classically: their T-matrices (orders 3 to 73) were known before 2026,
  those of order 73 from the near-normal sequences NN(36) of Đoković (2009). One of them, 32316 = 4 · 3 · 2693, is an
  order that Đoković claimed through Corollary 8.8(i).
- 6 need T-matrices of order 103, 151 or 163, known only since 2026, and are new beyond the tables: 16308, 21836,
  32012, 34556, 36844 and 39772.
- 16892 and 25132 are Known: SageMath's development version builds them by Cooper–Wallis since October 3, 2026, from
  the T-matrices of order 103 decoded from Schneider's matrix of order 2060.
- For seven of them gcd(s, t) > 1 (16300, 17532, 18036, 23724, 26604, 27036 and 37436), and the sequences live on ℤ_t
  × ℤ_s. For 25132, 36844 and 39772, q = 243 = 3⁵.
- Other routes reach some of them, and their READMEs say so. Cooper–Wallis with the T-matrices of order 131 of London
  and Kotsireas (2025) reaches 12052, 15196, 16244, 20436, 21484, 25676, 27772 and 38252 (we could not read that
  paper, which may list them). Cooper–Wallis with T-matrices posted in 2026 reaches 16300, 16308, 19932, 21836,
  32012, 33220, 34556, 36844 and 39772. Kiyasu's theorem reaches 15196 from our skew-Hadamard matrix of order 1900.
- 19028, 34196 and 39260 were already in this repository, by Cooper–Wallis.
- Not built: 28820 = 4 · 131 · 55, for want of T-matrices of order 55 (`search/tmat_mult.py` found none over its
  multiplier subgroups); it also follows by Cooper–Wallis from T-matrices of order 131.

SageMath's development version (`hadamard_matrix(4n, existence=True)`, October 17, 2026) returns Unknown for all the
others. Each certificate is checked exactly when it is written (entries ±1, Σ PAF = 0 off 0, sequences rebuilt from
Paley's matrix by `verify/verify.py`); the dense check H Hᵀ = 4n I of `verify/verify.py` is running.

## Reproduction and verification

    python3 search/paley_blocks.py check 1051          # the sixteen circulant blocks for q = 1051
    python3 search/paley_blocks.py certify 1051 3      # writes hadamard/3156/certificate.json and sequences.txt
    python3 verify/verify.py 3156                      # rebuilds and checks H H^T = 3156 I

`verify/verify.py` (construction `paley-blocks`) builds Paley's matrix independently, from the symplectic form on
PG(1, q). It finds the blocks from the certificate's τ, power and c, rebuilds the four sequences, and compares them
with the certificate. It then checks Σ PAF = 0, assembles the Goethals–Seidel matrix, and checks H Hᵀ = 4n I and the
SHA-256 of its canonical CSV. The matrices are not stored.
