# cmblocks: Miyamoto's quaternary array in blocks, times T-matrices

Programs behind [hadamard/QUATERNARY_MIYAMOTO_BLOCKS.md](../../hadamard/QUATERNARY_MIYAMOTO_BLOCKS.md) and
[hadamard/20712](../../hadamard/20712) ([search/NOTES.md](../NOTES.md) §105). Only numpy is needed (sympy for
`scan_t4.py`).

| File | Contents |
|---|---|
| `paley_circ.py` | Paley's quaternary Hadamard matrix of order (p + 1)/2 (p ≡ 3 mod 8 a prime) in 2 × 2 block-circulant form, through the non-split torus of SL(2, p) |
| `conf_circ.py` | K = C + iI (C Paley's symmetric conference matrix of order p + 1, p ≡ 1 mod 4 a prime) in the same form; run alone, it checks the theorem for p = 5, 17, 29, 41 |
| `border_search.py` | all 4 × 4 unit border pairs (A, B) with A A* = B B* = 4I, A B* + 2 B R* = 0 and A − 2R = B (512 per sign variant of R) |
| `cmb.py` | the sixteen quaternary bordered blocks, T-matrices from base sequences, and checks of H H* = 4qI, commutation and R_B-transposition |
| `gsplug.py` | Turyn doubles of the blocks, the signed permutation R′, and the real Goethals–Seidel plug-in |
| `make_cert.py` | `python3 make_cert.py p t OUT.json [paley|conf]`: a certificate for verify/verify.py (construction `quaternary-miyamoto-blocks`, without the SHA-256) |
| `rebuild.py` | `python3 rebuild.py CERT.json`: rebuilds and checks a certificate independently of verify.py, including the intermediate claims (about 6 minutes and 1 GB for 20712); `--csv` streams the matrix |
| `test_small.py`, `test_gs.py`, `test_nega.py` | exact checks for q = 7, 11, 23 (both K sources, t = 3, 5) and for negacyclic K |
| `old_border_test.py` | the borders of COMPLEX_MIYAMOTO.md: 102 of the 120 block pairs do not commute |
| `cx_counter.py` | the linear (quaternary) Goethals–Seidel plug-in fails at q = 7, t = 3 (largest entry of G G* − 84I is 42) |
| `gs_search2.py` | exhaustive search of Goethals–Seidel-type arrays over ±X, ±Xᵀ, ±conj X, ±X* (times R): none for generic quaternary blocks, millions for real ones |
| `scan_t4.py` | the q ≡ 3 (mod 4) with 4q open and K of the right form, and the orders 8qt they give; Table 4 from `../data/cati_pasechnik_table4.txt` |

```
cd search/cmblocks
python3 conf_circ.py                          # theorem at q = 7, 19, 31, 43 (K = C + iI), t = 3
python3 make_cert.py 1723 3 cert.json         # the certificate of hadamard/20712 (q = 863, t = 3)
python3 rebuild.py ../../hadamard/20712/certificate.json
```
