# T-matrices

T-matrices of order t are four circulant matrices with entries 0, ±1 whose supports partition ℤ_t and with Σ TᵢTᵢᵀ = t
I. By the Cooper–Wallis construction they give an orthogonal design OD(4t; t, t, t, t). Replacing its variables with
Williamson matrices of order w gives a Hadamard matrix of order 4tw.

Each file `T<t>.txt` lists the first rows of T₁, …, T₄, one per line, written over `0`, `+`, `-`. Order 71 is classical
(see the table). Orders 67, 101 and 131 come from `search/tmat_mult.py`: a multiplier subgroup, a random split of the
orbits into the supports of T₁ + T₂ and T₃ + T₄, and an exact join of the four complete sign lists (numpy for small
lists, `search/join4x.c` for large ones). The others were found by the search described in
[order 2060](../hadamard/2060): a small multiplier subgroup, parity conditions over GF(2), and meet-in-the-middle over
the signs.

| t | Earlier appearances known to us |
|----:|---|
| 67 | classical (T-sequences of length 67); this instance from `search/tmat_mult.py 67 29` (multiplier group of order 3) |
| 71 | classical: T-sequences from the Turyn-type sequences TT(24) of Best, Đoković, Kharaghani and Ramp (2013), as stored in SageMath's `t_sequences` module, through base sequences of lengths 47, 47, 24, 24 |
| 97 | none found; T-sequences of length 97 are open (Đoković, arXiv:1008.2043) |
| 101 | classical (T-sequences of length 1 + 10², Turyn 1974); this instance from `search/tmat_mult.py 101 95` (multiplier group of order 5) |
| 103 | decoded from the 2060 matrix in JD-Jones-ASES/Hadamard-T (2026-08-29) |
| 107 | classical: T-sequences from the Turyn-type sequences TT(36) of Kharaghani and Tayfeh-Rezaie (2005), through base sequences BS(71, 36), as SageMath's `t_sequences` module builds them (`T_sequences_smallcases(107)`); used for [20116](../hadamard/20116) |
| 109 | none found |
| 113 | Zuo–Xia–Xia (2013); Best–Đoković–Kharaghani–Ramp (2013) |
| 131 | Known: S. London and I. Kotsireas, New Turyn-type sequences, Cryptogr. Commun. 17 (2025) 1601–1610, construct TT(44), and Turyn's steps give base sequences BS(87, 44) and T-sequences of length 3 · 44 − 1 = 131. Not among the T-matrix orders listed by Seberry and Yamada (1992, §3: every t ≤ 210 except 73, 79, 83, 89, 97, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 183, 191, 193, 197 and 199). This instance from `search/tmat_mult.py 131 53` (multiplier group of order 5, row sums 5, 0, −9, 5) |
| 151 | posted by Schumer (2026-09-17) |
| 163 | JD-Jones-ASES/Hadamard-T (2026-08-29) |
| 181 | JD-Jones-ASES/Hadamard-T (2026-08-29) |
| 197 | none found |

For 97, 109, 151 and 163, which have the form 6m + 1, we could not check Lou, Xia and Xia, "On T-matrices of order
6m + 1" (Adv. Appl. Discrete Math. 21 (2019), paywalled). Treat them as possibly known. Against it: Seberry and Yamada's
book (2020, §1.16, p. 47), co-authored by Seberry, lists the known orders of T-matrices without 151 and 163, and Xia,
Zuo, Lou and Xia (Trans. Comb. 13 (2024)) mention neither order.

To check a file, the conditions are that every column of the 4 × t array has exactly one nonzero entry and that the
periodic autocorrelations of the four rows sum to t at shift 0 and to 0 elsewhere.

## Small orders used by the Paley-block orders

The orders beyond Table 4 in [hadamard/PALEY_BLOCKS.md](../hadamard/PALEY_BLOCKS.md) use T-matrices of small orders, all
classical: Seberry–Yamada (1992, §3) for every order but 73, and T-sequences of length 73 from the near-normal sequences
NN(36) of Đoković (arXiv:0912.5091, 2009). The files are the instances used there:

| t | Source of this instance |
|----:|---|
| 9, 11, 17, 21, 33 | T-sequences (1; 0), (0; (A + B)/2), (0; (A − B)/2), 0 from a Golay pair (A, B) of length t − 1 (Turyn 1974) |
| 15, 23, 25, 29, 35, 47, 53, 59 | T-sequences of length t as stored in SageMath's `t_sequences` module (`T_sequences_smallcases`) |
| 19, 31, 37, 41, 43, 49, 73 | `search/tmat_mult.py t g` with g = 7, 5, 10, 10, 6, 18, 8 (multiplier groups of orders 3, 3, 3, 5, 3, 3, 3); 37 is used by the Miyamoto-block order [39812](../hadamard/39812) |
| 55 | T-sequences of length 55 from Turyn's base sequences BS(6, 5) by Yang's multiplication by 5 (Seberry–Yamada 1992, Theorem 5.23, Golay pair of length 2): `search/yang5.py 5 --turyn`. SageMath stores none of length 55 |

Orders 3, 5, 7, 11 and 13 are also built by `search/paley_blocks.py` itself from base sequences; for t = 11 the order
11572 uses that built-in instance, the orders above 11804 the Golay instance of this directory.