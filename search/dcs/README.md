# dcs: dihedral Golay sequences and the construction of Kamali and Kharaghani

Programs behind [search/NOTES.md](../NOTES.md) §104. Kamali and Kharaghani (Australas. J. Combin. 18 (1998) 139–145)
call a pair (a, b) of sequences over D = {±1, ±i, ±s, ±is} (s² = 1, is = −si) a dihedral Golay sequence (DCS) if the
polynomials f_a, f_b, f_a*, f_b* commute and f_a f_a* + f_b f_b* = 2n. A DCS of length g free of ±is, Golay lengths g₁
and g₂ and their Lemmas 1–3 and Theorem 7 give a Hadamard matrix of order 8(g₁g + g₂) (their Corollary 9).

Labels: PROVED = argument in NOTES §104; CHECKED = exact computation here.

## Files

| File | Contents |
|---|---|
| `dcs.c` | exhaustive search of the ε-palindromic ansatz: a real, b = r + i·ι + s·s with r, i ε-palindromic and s (−ε)-palindromic; then the DCS conditions reduce to N_a + N_u + N_i = 2nδ (u = r + s), solved by a meet in the middle on N_a |
| `dcscheck.py` | independent checker of a DCS in the 2 × 2 real model (1, i, s ↦ I, J, S; the involution is the transpose): commutation of the four polynomials, complementarity, no ±is entry |
| `kk_build.py` | Golay pairs of lengths 2^a 10^b 26^c, the D-circulants and the Craigen–Goethals–Seidel array of Kamali–Kharaghani's Lemma 2, Lemma 3 (DCS × Golay pair) |
| `kk_big.py` | the whole chain for a given DCS, g₁ and g₂: the Hadamard matrix of order 8(g₁g + g₂), its exact check and the SHA-256 of its CSV |

## Usage

```
gcc -O3 -march=native -o dcs dcs.c
./dcs 23 1                      # n = 23, eps = +1: no solution (CHECKED, about 20 s)
./dcs 25 -1 3                   # n = 25, eps = -1: 16 solutions, the first three printed
python3 dcscheck.py A B         # A, B comma-separated tokens 1, -1, i, -i, s, -s
python3 kk_big.py 1,1,1,-1,-1,1,-1 i,1,s,1,-s,1,i 1 1024     # H(8248) = 8 (7 + 2^10), H H^T = 8248 I
python3 kk_big.py 1,1,1,-1,-1,1,-1 i,1,s,1,-s,1,i 1 832      # H(6712) = 8 (7 + 2^5 26)
python3 kk_big.py 1,1,1,-1,-1,1,-1 i,1,s,1,-s,1,i 256 1      # H(14344) = 8 (2^8 7 + 1)
```

The DCS of length 7 above is one of the 24 that `./dcs 7 1` finds. The orders 6712 and 8248 were known (Kamali and
Kharaghani 1998; the Hadamard Matrix Explorer posts matrices); the computation checks the construction. SHA-256 of the
CSVs (rows of 1 and −1, comma-separated): 6712 `a09c3ddc1bc7b7a9325cc6ead799b82320df58fb9e529f4e85cdc28e42470c8b`,
8248 `a90312bda7140bbb0da5deba54dc4d1921e477c1ba503ac5f6e69b13bf74d5ff`.

## Results (CHECKED)

Solutions of the ansatz, up to the signs of u and i and with a up to negation: n = 3: 2; 5: 4 (ε = −1); 7: 24;
9: 16 (ε = −1); 11: 4; 13: 8 (ε = +1) and 4 (ε = −1); 15: 8; 17, 19, 21, 23: none; 25: 16 (ε = −1); 27: none. For odd n,
ε = −1 needs 2n to be a sum of two squares (the sum of the antipalindromic i vanishes). Length 25 is not a complex
Golay length.
