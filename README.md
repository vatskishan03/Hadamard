# Hadamard

Explicit Hadamard matrices for orders that had no published construction, with
certificates that can be checked independently.

A Hadamard matrix of order n is an n × n matrix H with entries ±1 and HHᵀ = nI.
Such a matrix can exist only if n = 1, 2 or n ≡ 0 (mod 4), and the Hadamard
conjecture states that it exists for every multiple of 4. After the constructions
announced in August 2026 for the twelve unresolved orders below 2000 [1], the
smallest unresolved orders were 2060 and 2092.

## Results

| Order | n | Construction | Status |
|------:|---|---|---|
| [2092](2092) | 4 · 523 | Goethals–Seidel array, supplementary difference set in ℤ₅₂₃ invariant under a multiplier subgroup of order 9 | New |
| [2284](2284) | 4 · 571 | Goethals–Seidel, ℤ₅₇₁, multiplier subgroup of order 5; three blocks related by a multiplier of order 3 | New |
| [4492](4492) | 4 · 1123 | Goethals–Seidel, ℤ₁₁₂₃, multiplier subgroup of order 33 | New |
| [5932](5932) | 4 · 1483 | Goethals–Seidel, ℤ₁₄₈₃, multiplier subgroup of order 39 | New |
| [2060](2060) | 4 · 515 | Cooper–Wallis: T-matrices of order 103 and Williamson matrices of order 5 | Independent construction; the order was first settled publicly on 2026-08-23 [2] |
| [3532](3532) | 4 · 883 | Goethals–Seidel, ℤ₈₈₃, multiplier subgroup of order 21 | Order known [4]; new explicit matrix |

An order is marked new if it is listed as unresolved in the table of Cati and
Pasechnik [3], is not covered by the later updates of Đoković [4, 5] or by the
recursive constructions of Paley, Miyamoto, Yamada, Mathon and Spence, and no
construction was found in a literature search in October 2026. The table in [3]
is not complete (nine of its entries were settled in [4]), so these claims carry
that caveat.

Orders under investigation: [2148](2148), the smallest order we know to be
unresolved, and [4412](4412). [3676](3676) was settled by Schumer [8] while our
search for it was running.

## Skew-Hadamard results

A Hadamard matrix H is skew if H + Hᵀ = 2I. Skew-Hadamard matrices are conjectured to exist for every
order divisible by 4. Table 5 of Cati and Pasechnik [3] lists, for odd v ≤ 999, the smallest m for
which a skew-Hadamard matrix of order 2ᵐv is known. The orders below have m ≥ 3 there, so no skew
matrix of order 4v was known. Each matrix comes from a Goethals–Seidel difference family over ℤᵥ (for 1100, over the non-cyclic
group ℤ₅ × ℤ₅ × ℤ₁₁) with a skew first block. All four blocks are invariant under a multiplier subgroup H ≤ ℤᵥ* with −1 ∉ H.

| Order | v | H | Row sums | Table 5 entry |
|------:|---|---|---|---|
| [836](skew/836) | 209 (11 · 19) | ⟨26⟩, order 15 | 1, 9, 15, 23 | 209(4) |
| [1036](skew/1036) | 259 (7 · 37) | ⟨9⟩, order 9 | 1, 5, 7, 31 | 259(4) |
| [1100](skew/1100) | 275 (5² · 11), group ℤ₅ × ℤ₅ × ℤ₁₁ | Singer-3 on ℤ₅² × ⟨3⟩ on ℤ₁₁, order 15 | 1, 1, 3, 33 | 275(4) |
| [1148](skew/1148) | 287 (7 · 41) | ⟨16⟩, order 15 | 1, 3, 7, 33 | 287(4) |
| [1236](skew/1236) | 309 (3 · 103) | ⟨13⟩, order 17 | 1, 1, 3, 35 | 309(3) |
| [1364](skew/1364) | 341 (11 · 31) | ⟨14⟩, order 15 | 1, 9, 21, 29 | 341(4) |
| [1412](skew/1412) | 353 (prime) | ⟨22⟩, order 11 | 1, 21, 21, 23 | 353(4) |
| [1508](skew/1508) | 377 (13 · 29) | ⟨16⟩, order 21 | 1, 7, 27, 27 | 377(5) |
| [1612](skew/1612) | 403 (13 · 31) | ⟨100⟩, order 15 | 1, 5, 25, 31 | 403(5) |
| [1852](skew/1852) | 463 (prime) | ⟨33⟩, order 21 | 1, 1, 1, 43 | 463(7) |
| [1876](skew/1876) | 469 (7 · 67) | ⟨9⟩, order 33 | 1, 1, 5, 43 | 469(3) |
| [2404](skew/2404) | 601 (prime) | ⟨2⟩, order 25 | 1, 1, 1, 49 | 601(5) |
| [3556](skew/3556) | 889 (7 · 127) | ⟨2,25⟩, order 63 | 1, 5, 41, 43 | 889(4) |
| [3676](skew/3676) | 919 (prime) | ⟨58⟩, order 17; spin μ = 52 | 1, 35, 35, 35 | 919(3) |

Doubling (Seberry) then gives skew matrices of orders 2ᵏ · 4v for every k ≥ 0, which settles every
entry v(m) above at m = 2. The searches use `attack/skew_triage.py`, `attack/skew_rank.py`, `attack/skew_ktree.py` with
`attack/match4.c`, and, for the spin structure of 3676, `attack/spin_skew.py` with `attack/spin_ls.c`.
Each directory's README describes its search. `verify/verify.py` checks
H + Hᵀ = 2I for these certificates.

## Verification

```
python3 verify/verify.py
```

rebuilds every matrix from its `certificate.json`, checks the defining identity,
checks HHᵀ = nI exactly and compares the SHA-256 digest of the canonical CSV
encoding with the recorded value. It needs NumPy and takes about two minutes.

An independent checker in C verifies the shipped matrices directly:

```
cc -O2 -fopenmp verify/check_matrix.c -o check_matrix
gzip -dc 4492/H4492.csv.gz | ./check_matrix
```

## Layout

Each order has its own directory containing

- `certificate.json`: the construction data from which the matrix is rebuilt,
- `sequences.txt`: for Goethals–Seidel orders, the first rows of the four circulant blocks,
- `H<n>.csv.gz`: the matrix, one row per line, entries separated by commas,
- `README.md`: construction details and status.

Other top-level contents:

- [`skew`](skew): skew-Hadamard matrices, one directory per order, in the same layout.
- [`t-matrices`](t-matrices): the T-matrices of orders 97 to 197 found during this work.
- [`search`](search): the programs that found the Goethals–Seidel orders.
- [`METHOD.md`](METHOD.md): how the search works.
- [`FRONTIER.md`](FRONTIER.md): the unresolved orders up to 4 · 2999, with the cheapest structure for each and its predicted cost.

## References

1. L. Alpöge, P. Voinov, S. Reynolds-Haertle, Hadamard matrices of orders 668, 716,
   892, 1132, 1244, 1388, 1436, 1676, 1772, 1916, 1948 and 1964, announcement,
   August 12, 2026.
2. Hadamard matrix of order 2060, public gist by Schneider, August 23, 2026; decoded
   in JD-Jones-ASES/Hadamard-T.
3. M. Cati, D. V. Pasechnik, A database of constructions of Hadamard matrices,
   arXiv:2411.18897v2 (2025), Table 4.
4. D. Ž. Đoković, Small orders of Hadamard matrices and base sequences,
   Int. Math. Forum 6 (2011) 3061–3067, arXiv:1008.2043.
5. D. Ž. Đoković, Hadamard matrices of small order and Yang conjecture,
   J. Combin. Des. 18 (2010) 254–259, arXiv:0912.5091.
6. J.-M. Goethals, J. J. Seidel, Orthogonal matrices with zero diagonal,
   Canad. J. Math. 19 (1967) 1001–1010.
7. J. Cooper, J. Wallis, A construction for Hadamard arrays,
   Bull. Austral. Math. Soc. 7 (1972) 269–278.
8. Schumer, Hadamard matrices of orders 3676, 5436, 6644 and 7852, Zenodo,
   September 17, 2026, doi:10.5281/zenodo.22820997.
