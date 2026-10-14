# Near-Williamson matrices

Near-Williamson matrices of order w are four circulant ±1 matrices A, B, C, D of order w with B, C, D symmetric and
A Aᵀ + B Bᵀ + C Cᵀ + D Dᵀ = 4w I (H. Kharaghani, A. Mohammadian, B. Tayfeh-Rezaie and V. Zaitsev, arXiv:2605.08661v2,
September 2026). With R the back-identity, (AR, B, C, D) are Williamson-type matrices: four ±1 matrices with
WₖWₗᵀ = WₗWₖᵀ and Σ WₖWₖᵀ = 4w I. Williamson-type matrices of order w and T-matrices of order t give a Hadamard
matrix of order 4tw by the Cooper–Wallis construction.

Each file `NW<w>.txt` lists the first rows of A, B, C, D, one per line, written over `+` and `-`. The files were
converted from the database of the authors (cs.uleth.ca/~zaitsevv/nw/NW<w>.txt, retrieved October 2026, where `1` stands
for `+`). The quadruples of orders 59, 65 and 67 are those printed in Tables 5, 6 and 7 of the paper; order 103 is from
the database only. Order 111 is ours (see the table).

| w | Earlier appearances known to us |
|----:|---|
| 59 | arXiv:2605.08661v2 (September 2026, Table 5); v1 (May 2026) already contained near-Williamson matrices of order 59 and the quaternary Hadamard matrix of order 118; no Williamson matrices of order 59 exist |
| 65 | arXiv:2605.08661v2 (September 2026, Table 6) |
| 67 | arXiv:2605.08661v2 (September 2026, Table 7) |
| 103 | arXiv:2605.08661v2 database (cs.uleth.ca/~zaitsevv/nw/NW103.txt) |
| 111 | none found: **this repository** (October 2026). Not among the orders of arXiv:2605.08661v2: every odd w ≤ 69 plus 73, 79, 91, 93, 95, 97, 103, 109, 127, 133 and 181, including those it gives only electronically. Table A.1 of Seberry–Yamada (1992) records only Williamson-type matrices of order 222 for 111 (t = 3). Found by `attack/wt_twist.py search 111 10 38 --plain`, the order-3 twist of the paper's §4.4 with u = 10 and K′ = ⟨38⟩: A is constant on the orbits of ⟨10, 38⟩, B is symmetric and ⟨38⟩-invariant, C = B∘10 and D = B∘100. Rows 21, 1, 1, 1; the exact two-list join of this structure holds 36 quadruples with these row sums. |

To check a file, the conditions are that rows 2 to 4 are symmetric (b_k = b_(−k)) and that the periodic autocorrelations of
the four rows sum to 4w at shift 0 and to 0 elsewhere.
