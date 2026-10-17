# Hadamard matrix of order 21036

21036 = 4 · 5259, with 5259 = 3 · 1753.

## Construction

Cooper–Wallis with symmetric Williamson-type matrices of order 1753 that come from the symmetric Hadamard matrix of
order 876 of this repository ([symmetric/876](../../symmetric/876)).

**Williamson-type matrices of order 1753.** q = 1753 is a prime ≡ 1 (mod 8), and n = (q − 1)/2 = 876. Seberry and
Yamada's Corollary 8.8(ii) (1992; Corollary 29(ii) of JCMCC 7, 1990; after Miyamoto) turns a symmetric Hadamard
matrix M of order n into symmetric Williamson-type matrices of order q, through their Lemma 8.6 (Lemma 26 of JCMCC 7).
With g = 7 a primitive root mod q, s_i = g^(2i), t_j = g^(1 − 2j) (0 ≤ i, j < n) and χ the Legendre symbol:

    U₁ = −[χ(s_i − s_j)],   U₂ = [χ(s_i − t_j)],   U₃ = U₄ = 0,   V = (I, 0, M, M),
    S_k = U_k ⊗ [[1, 1], [1, 1]] + V_k ⊗ [[1, −1], [−1, 1]],
    W₁ = [[1, −e], [−eᵀ, S₁]],   W_k = [[1, e], [eᵀ, S_k]]  (k = 2, 3, 4).

The W_k are symmetric ±1 matrices of order 1753, pairwise amicable (W_k W_lᵀ = W_l W_kᵀ), with Σ W_k W_kᵀ = 4 · 1753 I.
`verify/verify.py` builds them and checks this exactly. Part (i) of the same corollary, with a non-symmetric input,
does not give Williamson-type matrices as printed ([NOTES.md](../../NOTES.md) §7); part (ii), used here, does.

**Cooper–Wallis.** T-matrices of order 3 (first rows `+00`, `0+0`, `00+`, `000`) give an orthogonal design
OD(12; 3, 3, 3, 3) through the Goethals–Seidel array. Substituting W₁, …, W₄ for its four variables gives the
Hadamard matrix of order 4 · 3 · 1753 = 21036 (Seberry–Yamada 1992, Theorem 3.8).

## Search

None beyond the symmetric matrix of order 876, which came from an exact propus chain search
([symmetric/876](../../symmetric/876)). `python3 search/miyamoto_wt.py cw 1753 symmetric/876/certificate.json T3.txt`,
with T3.txt holding the four rows above, writes the certificate.

## Status

**New, beyond the tables.** In Table 1.53 of the Handbook of Combinatorial Designs (2nd ed., 2007), n = 5259 has
t = 3: the smallest known Hadamard matrix with odd part 5259 had order 2³ · 5259. The order lies above the range of
Cati–Pasechnik Table 4, which stops at n = 2999.

- Đoković (arXiv:0912.5091) lists 5259 in his set Δ of entries that should have had t = 2. His proof for Δ takes
  Williamson-type matrices from Seberry–Yamada's Table A.1 (1992), and the only Williamson-type matrices of order
  1753 there are those of key wl, from Corollary 8.8(i) with Williamson-type matrices of order 438. The printed
  proof of Corollary 8.8(i) does not give Williamson-type matrices ([NOTES.md](../../NOTES.md) §7). Part (ii) needs a
  symmetric Hadamard matrix of order 876 = 4 · 219, and 219 was undecided in the Đoković–Seberry census of symmetric
  Hadamard matrices until [symmetric/876](../../symmetric/876).
- Paley's constructions do not apply: 4n − 1 = 21035 = 5 · 7 · 601 and 2n − 1 = 10517 = 13 · 809 are not prime
  powers.
- Miyamoto and Yamada–Kiyasu do not apply: 5259 is not a prime power, and neither is n − 2 = 5257 = 7 · 751.
- Turyn's Williamson-type matrices of order 1753 would need 2 · 1753 − 1 = 3505 = 5 · 701 to be a prime power.
  Cooper–Wallis with the other factorization would need T-matrices of order 1753, and Yang's products would need
  base sequences of total length 1753 or 5259; none is known to us. The Paley-block route does not reach 5259
  ([PALEY_BLOCKS.md](../PALEY_BLOCKS.md)).
- SageMath's development version (October 17, 2026) does not construct it: `hadamard_matrix(21036, existence=True)`
  returns Unknown.

We found no valid construction of order 21036 before this one. The sources searched are listed in
[NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, g, the input certificate and its SHA-256, the T-matrices and the SHA-256 of the matrix |

The matrix itself is not stored; `verify/verify.py` rebuilds it from the symmetric matrix of order 876.

## Verification

```
python3 verify/verify.py 21036
```

SHA-256 of `H21036.csv`: `66397e01284a7f4655ed19e1fe13e8940cf25318ab8bf5815f10a9ab33b1db79`
