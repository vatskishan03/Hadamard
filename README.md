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
unresolved, [3676](3676) and [4412](4412).

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
