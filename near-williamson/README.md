# Near-Williamson matrices

Near-Williamson matrices of order w are four circulant ±1 matrices A, B, C, D of order w with B, C, D symmetric and
A Aᵀ + B Bᵀ + C Cᵀ + D Dᵀ = 4w I (H. Kharaghani, A. Mohammadian, B. Tayfeh-Rezaie and V. Zaitsev, arXiv:2605.08661v2,
September 2026). With R the back-identity, (AR, B, C, D) are Williamson-type matrices: four ±1 matrices with
WₖWₗᵀ = WₗWₖᵀ and Σ WₖWₖᵀ = 4w I. Williamson-type matrices of order w and T-matrices of order t give a Hadamard
matrix of order 4tw by the Cooper–Wallis construction.

Each file `NW<w>.txt` lists the first rows of A, B, C, D, one per line, written over `+` and `-`. The files were
converted from the database of the authors (cs.uleth.ca/~zaitsevv/nw/NW<w>.txt, retrieved October 2026, where `1` stands
for `+`). The quadruples of orders 59, 65 and 67 are those printed in Tables 5, 6 and 7 of the paper; order 103 is from
the database only.

| w | Earlier appearances known to us |
|----:|---|
| 59 | arXiv:2605.08661v2 (September 2026, Table 5); v1 (May 2026) already contained near-Williamson matrices of order 59 and the quaternary Hadamard matrix of order 118; no Williamson matrices of order 59 exist |
| 65 | arXiv:2605.08661v2 (September 2026, Table 6) |
| 67 | arXiv:2605.08661v2 (September 2026, Table 7) |
| 103 | arXiv:2605.08661v2 database (cs.uleth.ca/~zaitsevv/nw/NW103.txt) |

To check a file, the conditions are that rows 2 to 4 are symmetric (b_k = b_(−k)) and that the periodic autocorrelations of
the four rows sum to 4w at shift 0 and to 0 elsewhere.
