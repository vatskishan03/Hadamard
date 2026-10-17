# Hadamard matrices

Each directory holds one order: `certificate.json`, the matrix `H<n>.csv.gz` where it is small enough to store, and a
README with the construction, the search and the evidence for its status. Status labels are defined in
[NOVELTY.md](../NOVELTY.md).

## Orders that were unresolved

Entries n of Cati–Pasechnik Table 4 (odd n ≤ 2999) for which no Hadamard matrix of order 4n was known.

| Order | n | Construction |
|------:|---|---|
| [2092](2092) | 523 | Goethals–Seidel over ℤ₅₂₃, multiplier subgroup of order 9 |
| [2148](2148) | 537 = 3 · 179 | bordered Goethals–Seidel array over ℤ₃ × ℤ₁₇₈ with a border of size 12, multiplier subgroup of order 11 |
| [2284](2284) | 571 | Goethals–Seidel over ℤ₅₇₁, multiplier subgroup of order 5, three blocks related by a multiplier of order 3 |
| [2292](2292) | 573 = 3 · 191 | Spence's bordered construction over ℤ₂ × ℤ₂ × ℤ₁₄₃, three blocks related by an automorphism of order 15 |
| [2636](2636) | 659 | bordered Goethals–Seidel array over ℤ₄ × ℤ₄ × ℤ₄₁ with a border of size 12, symmetry group of order 15 |
| [3732](3732) | 933 = 3 · 311 | Goethals–Seidel over ℤ₃ × GF(311), multiplier subgroup of order 31; the matrix is skew |
| [3964](3964) | 991 | bordered Goethals–Seidel array over ℤ₂ × ℤ₂ × ℤ₂₄₇ with a border of size 12, symmetry group of order 27 |
| [4492](4492) | 1123 | Goethals–Seidel over ℤ₁₁₂₃, multiplier subgroup of order 33 |
| [5020](5020) | 1255 = 5 · 251 | Goethals–Seidel over ℤ₅ × GF(251), multiplier subgroup of order 50 |
| [5932](5932) | 1483 | Goethals–Seidel over ℤ₁₄₈₃, multiplier subgroup of order 39 |
| [9172](9172) | 2293 | Miyamoto's theorem (q = 2293) from the matrix of order 2292 |

Doubling the matrices of orders 2636 and 3732 gives orders 5272 and 7464, which were open as well.

## Orders beyond the tables

These orders lie above the range of Cati–Pasechnik Table 4. They come from recursive theorems (Miyamoto,
Yamada–Kiyasu, Cooper–Wallis) applied to matrices of this repository or to objects published in 2026. Most of their
matrices are not stored (up to 7 GB as CSV); `python3 verify/verify.py <order>` rebuilds them.

| Order | n | Construction | Status |
|------:|---|---|---|
| [12596](12596) | 3149 = 47 · 67 | Cooper–Wallis: T-matrices of order 67, bordered Williamson-type matrices of order 47 | New, beyond the tables |
| [14708](14708) | 3677 | Miyamoto (q = 3677) from the skew-Hadamard matrix of order [3676](../skew/3676) | New, beyond the tables |
| [14932](14932) | 3733 | Miyamoto (q = 3733) from the matrix of order 3732 | First explicit construction |
| [19028](19028) | 4757 = 67 · 71 | Cooper–Wallis: T-sequences of length 71, Williamson-type matrices of order 67 from near-Williamson matrices | New, beyond the tables |
| [19364](19364) | 4841 = 47 · 103 | Cooper–Wallis: T-matrices of order 103, bordered Williamson-type matrices of order 47 | New, beyond the tables |
| [20084](20084) | 5021 | Miyamoto (q = 5021) from the matrix of order 5020 | First explicit construction |
| [21748](21748) | 5437 | Miyamoto (q = 5437) from a Cooper–Wallis matrix of order 5436 (T-matrices of order 151, Williamson matrices of order 9) | New, beyond the tables |
| [24124](24124) | 6031 | Yamada–Kiyasu (q = 6029) from a skew-Hadamard matrix of order 3016, the double of [skew 1508](../skew/1508) | New, beyond the tables |
| [24308](24308) | 6077 = 59 · 103 | Cooper–Wallis: T-matrices of order 103, near-Williamson matrices of order 59 | New, beyond the tables |
| [24412](24412) | 6103 | Yamada–Kiyasu (q = 6101) from [skew 3052](../skew/3052) | New, beyond the tables |
| [27604](27604) | 6901 = 67 · 103 | Cooper–Wallis: T-matrices of order 103, near-Williamson matrices of order 67 | New, beyond the tables |
| [29404](29404) | 7351 | Yamada–Kiyasu (q = 7349) from [skew 3676](../skew/3676) | First explicit construction |
| [31316](31316) | 7829 | Miyamoto (q = 7829) from a Cooper–Wallis matrix of order 7828 (T-matrices of order 103, Williamson matrices of order 19) | New, beyond the tables |
| [31412](31412) | 7853 | Miyamoto (q = 7853) from a Cooper–Wallis matrix of order 7852 (T-matrices of order 151, Williamson matrices of order 13) | New, beyond the tables |
| [34196](34196) | 8549 = 83 · 103 | Cooper–Wallis: T-matrices of order 103, bordered Williamson-type matrices of order 83 | New, beyond the tables |
| [35636](35636) | 8909 = 59 · 151 | Cooper–Wallis: T-matrices of order 151, near-Williamson matrices of order 59 | New, beyond the tables |
| [38468](38468) | 9617 = 59 · 163 | Cooper–Wallis: T-matrices of order 163, near-Williamson matrices of order 59 | New, beyond the tables |
| [39260](39260) | 9815 = 5 · 13 · 151 | Cooper–Wallis: T-matrices of order 151, near-Williamson matrices of order 65 | New, beyond the tables |
| [51580](51580) | 12895 = 5 · 2579 | Yamada–Kiyasu (q = 12893) from a skew-Hadamard matrix of order 6448, [skew 1612](../skew/1612) doubled twice | New, beyond the tables |

## Known orders, independent matrices

| Order | n | Construction | First construction |
|------:|---|---|---|
| [2060](2060) | 515 = 5 · 103 | Cooper–Wallis: T-matrices of order 103, Williamson matrices of order 5 | Schneider, 2026-08-23 |
| [3532](3532) | 883 | Goethals–Seidel over ℤ₈₈₃, multiplier subgroup of order 21 | Đoković, arXiv:1008.2043 (2010) |
