# Skew-Hadamard matrix of order 7272

7272 = 8 · 909, with 909 = 3² · 101. 7272 = 3(q + 1) with q = 2423, a prime ≡ 7 (mod 8).

## Construction

A skew-Hadamard matrix is a Hadamard matrix H with H + Hᵀ = 2I. This one is the case q = 2423 of the construction in
[PALEY_SKEW.md](../PALEY_SKEW.md), which gives a skew-Hadamard matrix of order 3(q + 1) for every prime power q ≡ 3 (mod
4).

- **Paley's matrix in sixteen blocks.** Paley's skew-Hadamard matrix of order q + 1 = 2424 on PG(1, 2423), signed and
  listed along the four orbits, of length s = 606, of g = (x ↦ −1/(x + 3))² as in
  [PALEY_BLOCKS.md](../../hadamard/PALEY_BLOCKS.md) (sign c = 1), is a 4 × 4 array of negacyclic s × s blocks N_kj: s is
  even, and g^s acts as −I.
- **The pairing.** After the base point of orbit 3 is moved by 685 (`shifts`), the blocks of column 0 in rows 1 and 3
  satisfy N_30 = −N_10ᵀ. The pairing comes from an anti-automorphism of Paley's matrix that fixes two of the orbits and
  swaps the other two.
- **T-matrices of order 3.** T_0 = I, T_1 = Z, T_3 = Z² and T_2 = 0, with Z the cyclic shift (`T`, one row per block
  row). The matrices A_j = Σ_k T_k ⊗ N_kj are ±1 matrices on ℤ₃ × N₆₀₆ (circulant in ℤ₃, negacyclic of length 606), Σ_j
  A_j A_jᵀ = 7272 I, and A_0 = I ⊗ N_00 + Z ⊗ N_10 − Z² ⊗ N_10ᵀ is skew-type.
- **The matrix.** The Goethals–Seidel array with A_0 first and R: (i, a) ↦ (−i, 605 − a) is a skew-Hadamard matrix of
  order 7272. `sequences.txt` lists the four A_j as sequences on ℤ₃ × ℤ_s, at index i s + a, A_0 first.

## Search

None: the matrix is determined by q. `python3 search/paley_blocks.py skew 2423` writes the certificate, and `python3
search/skew_3q.py 2423` builds a matrix of the same family from the closed formula of [PALEY_SKEW.md](../PALEY_SKEW.md)
and checks it.

## Status

**New** as a skew-Hadamard matrix ([PALEY_SKEW.md](../PALEY_SKEW.md#status)).

- Cati–Pasechnik Table 5 gives 909(4): the smallest known skew order with odd part 909 was 16 · 909 = 14544. This matrix
  has order 8 · 909.
- The Đoković–Seberry census covers the orders 4v only; it lists 909 as undecided, and 3636 remains open.

Doubling (Seberry) gives skew-Hadamard matrices of all orders 2ᵏ · 7272. The sources searched are listed in
[NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, τ, power and c (the blocks), shifts and columns (the pairing), the T-matrices, the four sequences and the SHA-256 of the matrix |
| `sequences.txt` | the four sequences, one per line |

The matrix itself is not stored; `verify/verify.py` rebuilds it.

## Verification

```
python3 verify/verify.py skew/7272
```

SHA-256 of `H7272.csv`: `69867f9200a66e203108aaeea72d8d128710f35633f23a877bafa7500908dff2`
