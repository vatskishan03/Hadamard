# Skew-Hadamard matrix of order 7656

7656 = 8 · 957, with 957 = 3 · 11 · 29. 7656 = 3(q + 1) with q = 2551, a prime ≡ 7 (mod 8).

## Construction

A skew-Hadamard matrix is a Hadamard matrix H with H + Hᵀ = 2I. This one is the case q = 2551 of the construction in
[PALEY_SKEW.md](../PALEY_SKEW.md), which gives a skew-Hadamard matrix of order 3(q + 1) for every prime power q ≡ 3 (mod
4).

- **Paley's matrix in sixteen blocks.** Paley's skew-Hadamard matrix of order q + 1 = 2552 on PG(1, 2551), signed and
  listed along the four orbits, of length s = 638, of g = (x ↦ −1/(x + 4))² as in
  [PALEY_BLOCKS.md](../../hadamard/PALEY_BLOCKS.md) (sign c = 1), is a 4 × 4 array of negacyclic s × s blocks N_kj: s is
  even, and g^s acts as −I.
- **The pairing.** After the base point of orbit 3 is moved by 854 (`shifts`), the blocks of column 0 in rows 1 and 3
  satisfy N_30 = −N_10ᵀ. The pairing comes from an anti-automorphism of Paley's matrix that fixes two of the orbits and
  swaps the other two.
- **T-matrices of order 3.** T_0 = I, T_1 = Z, T_3 = Z² and T_2 = 0, with Z the cyclic shift (`T`, one row per block
  row). The matrices A_j = Σ_k T_k ⊗ N_kj are ±1 matrices on ℤ₃ × N₆₃₈ (circulant in ℤ₃, negacyclic of length 638), Σ_j
  A_j A_jᵀ = 7656 I, and A_0 = I ⊗ N_00 + Z ⊗ N_10 − Z² ⊗ N_10ᵀ is skew-type.
- **The matrix.** The Goethals–Seidel array with A_0 first and R: (i, a) ↦ (−i, 637 − a) is a skew-Hadamard matrix of
  order 7656. `sequences.txt` lists the four A_j as sequences on ℤ₃ × ℤ_s, at index i s + a, A_0 first.

## Search

None: the matrix is determined by q. `python3 search/paley_blocks.py skew 2551` writes the certificate, and `python3
search/skew_3q.py 2551` builds a matrix of the same family from the closed formula of [PALEY_SKEW.md](../PALEY_SKEW.md)
and checks it.

## Status

**New** as a skew-Hadamard matrix ([PALEY_SKEW.md](../PALEY_SKEW.md#status)).

- Cati–Pasechnik Table 5 gives 957(4): the smallest known skew order with odd part 957 was 16 · 957 = 15312. This matrix
  has order 8 · 957.
- The Đoković–Seberry census covers the orders 4v only; it lists 957 as undecided, and 3828 remains open.

Doubling (Seberry) gives skew-Hadamard matrices of all orders 2ᵏ · 7656. The sources searched are listed in
[NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, τ, power and c (the blocks), shifts and columns (the pairing), the T-matrices, the four sequences and the SHA-256 of the matrix |
| `sequences.txt` | the four sequences, one per line |

The matrix itself is not stored; `verify/verify.py` rebuilds it.

## Verification

```
python3 verify/verify.py skew/7656
```

SHA-256 of `H7656.csv`: `29bf759171beb4bbb141ac306ee33df8b503c1e05072c8106ccd20daa17dffd1`
