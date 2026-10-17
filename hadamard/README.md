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
| [3156](3156) | 789 = 3 · 263 | Paley's matrix of order 1052 in sixteen circulant blocks, T-matrices of order 3 |
| [3732](3732) | 933 = 3 · 311 | Goethals–Seidel over ℤ₃ × GF(311), multiplier subgroup of order 31; the matrix is skew |
| [3964](3964) | 991 | bordered Goethals–Seidel array over ℤ₂ × ℤ₂ × ℤ₂₄₇ with a border of size 12, symmetry group of order 27 |
| [4492](4492) | 1123 | Goethals–Seidel over ℤ₁₁₂₃, multiplier subgroup of order 33 |
| [5020](5020) | 1255 = 5 · 251 | Goethals–Seidel over ℤ₅ × GF(251), multiplier subgroup of order 50 |
| [5260](5260) | 1315 = 5 · 263 | Paley's matrix of order 1052 in sixteen circulant blocks, T-matrices of order 5 |
| [5932](5932) | 1483 | Goethals–Seidel over ℤ₁₄₈₃, multiplier subgroup of order 39 |
| [6036](6036) | 1509 = 3 · 503 | Paley's matrix of order 2012 in sixteen circulant blocks, T-matrices of order 3 |
| [6356](6356) | 1589 = 7 · 227 | Paley's matrix of order 908 in sixteen circulant blocks, T-matrices of order 7 |
| [6756](6756) | 1689 = 3 · 563 | Paley's matrix of order 2252 in sixteen circulant blocks, T-matrices of order 3 |
| [7364](7364) | 1841 = 7 · 263 | Paley's matrix of order 1052 in sixteen circulant blocks, T-matrices of order 7 |
| [7660](7660) | 1915 = 5 · 383 | Paley's matrix of order 1532 in sixteen circulant blocks, T-matrices of order 5 |
| [8620](8620) | 2155 = 5 · 431 | Paley's matrix of order 1724 in sixteen circulant blocks, T-matrices of order 5 |
| [9172](9172) | 2293 | Miyamoto's theorem (q = 2293) from the matrix of order 2292 |
| [10060](10060) | 2515 = 5 · 503 | Paley's matrix of order 2012 in sixteen circulant blocks, T-matrices of order 5 |
| [10644](10644) | 2661 = 3 · 887 | Paley's matrix of order 3548 in sixteen circulant blocks, T-matrices of order 3 |
| [10932](10932) | 2733 = 3 · 911 | Paley's matrix of order 3644 in sixteen circulant blocks, T-matrices of order 3 |
| [11260](11260) | 2815 = 5 · 563 | Paley's matrix of order 2252 in sixteen circulant blocks, T-matrices of order 5 |
| [11572](11572) | 2893 = 11 · 263 | Paley's matrix of order 1052 in sixteen circulant blocks, T-matrices of order 11 |
| [11804](11804) | 2951 = 13 · 227 | Paley's matrix of order 908 in sixteen circulant blocks, T-matrices of order 13 |

2284 is Concurrent: it also follows, by Turyn's doubling, from the quaternary Hadamard matrix of order 1142 on
Zaitsev's page ([2284](2284/README.md#status)). The others are New. The fourteen orders from Paley's matrix in
sixteen circulant blocks were obtainable from two published theorems that the tables did not combine
([PALEY_BLOCKS.md](PALEY_BLOCKS.md)).

Doubling the matrix of order 3732 gives order 7464, which was open as well, and doubling 9172 gives 18344, 36688 and
73376 (the table entry is 2293(6)).

## Orders beyond the tables

These orders lie above the range of Cati–Pasechnik Table 4. They come from recursive theorems (Miyamoto,
Yamada–Kiyasu, Cooper–Wallis) applied to matrices of this repository or to objects published in 2026, and,
for the Handbook entries from n = 3013 on, from Paley's matrix in sixteen circulant blocks times T-matrices
([PALEY_BLOCKS.md](PALEY_BLOCKS.md#beyond-table-4)). Most of their
matrices are not stored (up to 7 GB as CSV); `python3 verify/verify.py <order>` rebuilds them.

| Order | n | Construction | Status |
|------:|---|---|---|
| [12052](12052) | 3013 = 23 · 131 | Paley's matrix of order 524 in sixteen circulant blocks, T-matrices of order 23 | New, beyond the tables; obtainable classically ([PALEY_BLOCKS.md](PALEY_BLOCKS.md)) |
| [12068](12068) | 3017 = 7 · 431 | Paley's matrix of order 1724 in sixteen circulant blocks, T-matrices of order 7 | New, beyond the tables; obtainable classically ([PALEY_BLOCKS.md](PALEY_BLOCKS.md)) |
| [12596](12596) | 3149 = 47 · 67 | Cooper–Wallis: T-matrices of order 67, bordered Williamson-type matrices of order 47 | New, beyond the tables |
| [13076](13076) | 3269 = 7 · 467 | Paley's matrix of order 1868 in sixteen circulant blocks, T-matrices of order 7 | New, beyond the tables; obtainable classically ([PALEY_BLOCKS.md](PALEY_BLOCKS.md)) |
| [13092](13092) | 3273 = 3 · 1091 | Paley's matrix of order 4364 in sixteen circulant blocks, T-matrices of order 3 | New, beyond the tables; obtainable classically ([PALEY_BLOCKS.md](PALEY_BLOCKS.md)) |
| [13676](13676) | 3419 = 13 · 263 | Paley's matrix of order 1052 in sixteen circulant blocks, T-matrices of order 13 | New, beyond the tables; obtainable classically ([PALEY_BLOCKS.md](PALEY_BLOCKS.md)) |
| [13812](13812) | 3453 = 3 · 1151 | Paley's matrix of order 4604 in sixteen circulant blocks, T-matrices of order 3 | New, beyond the tables; obtainable classically ([PALEY_BLOCKS.md](PALEY_BLOCKS.md)) |
| [14708](14708) | 3677 | Miyamoto (q = 3677) from the skew-Hadamard matrix of order [3676](../skew/3676) | New, beyond the tables |
| [14860](14860) | 3715 = 5 · 743 | Paley's matrix of order 2972 in sixteen circulant blocks, T-matrices of order 5 | New, beyond the tables; obtainable classically ([PALEY_BLOCKS.md](PALEY_BLOCKS.md)) |
| [14932](14932) | 3733 | Miyamoto (q = 3733) from the matrix of order 3732 | First explicit construction |
| [15196](15196) | 3799 = 29 · 131 | Paley's matrix of order 524 in sixteen circulant blocks, T-matrices of order 29 | New, beyond the tables; obtainable classically ([PALEY_BLOCKS.md](PALEY_BLOCKS.md)) |
| [15324](15324) | 3831 = 3 · 1277 | Paley's matrix of order 5108 in sixteen circulant blocks, T-matrices of order 3 | New, beyond the tables; obtainable classically ([PALEY_BLOCKS.md](PALEY_BLOCKS.md)) |
| [15764](15764) | 3941 = 7 · 563 | Paley's matrix of order 2252 in sixteen circulant blocks, T-matrices of order 7 | New, beyond the tables; obtainable classically ([PALEY_BLOCKS.md](PALEY_BLOCKS.md)) |
| [15780](15780) | 3945 = 3 · 5 · 263 | Paley's matrix of order 1052 in sixteen circulant blocks, T-matrices of order 15 | New, beyond the tables; obtainable classically ([PALEY_BLOCKS.md](PALEY_BLOCKS.md)) |
| [16044](16044) | 4011 = 3 · 7 · 191 | Paley's matrix of order 5348 in sixteen circulant blocks, T-matrices of order 3 | New, beyond the tables; obtainable classically ([PALEY_BLOCKS.md](PALEY_BLOCKS.md)) |
| [16244](16244) | 4061 = 31 · 131 | Paley's matrix of order 524 in sixteen circulant blocks, T-matrices of order 31 | New, beyond the tables; obtainable classically ([PALEY_BLOCKS.md](PALEY_BLOCKS.md)) |
| [16300](16300) | 4075 = 5² · 163 | Paley's matrix of order 3260 in sixteen circulant blocks, T-matrices of order 5 (Goethals–Seidel array over ℤ₅ × ℤ₈₁₅) | New, beyond the tables; obtainable classically ([PALEY_BLOCKS.md](PALEY_BLOCKS.md)) |
| [16308](16308) | 4077 = 3³ · 151 | Paley's matrix of order 108 in sixteen circulant blocks, T-matrices of order 151 | New, beyond the tables |
| [16436](16436) | 4109 = 7 · 587 | Paley's matrix of order 2348 in sixteen circulant blocks, T-matrices of order 7 | New, beyond the tables; obtainable classically ([PALEY_BLOCKS.md](PALEY_BLOCKS.md)) |
| [16852](16852) | 4213 = 11 · 383 | Paley's matrix of order 1532 in sixteen circulant blocks, T-matrices of order 11 | New, beyond the tables; obtainable classically ([PALEY_BLOCKS.md](PALEY_BLOCKS.md)) |
| [17252](17252) | 4313 = 19 · 227 | Paley's matrix of order 908 in sixteen circulant blocks, T-matrices of order 19 | New, beyond the tables; obtainable classically ([PALEY_BLOCKS.md](PALEY_BLOCKS.md)) |
| [17532](17532) | 4383 = 3² · 487 | Paley's matrix of order 5844 in sixteen circulant blocks, T-matrices of order 3 (Goethals–Seidel array over ℤ₃ × ℤ₁₄₆₁) | New, beyond the tables; obtainable classically ([PALEY_BLOCKS.md](PALEY_BLOCKS.md)) |
| [17740](17740) | 4435 = 5 · 887 | Paley's matrix of order 3548 in sixteen circulant blocks, T-matrices of order 5 | New, beyond the tables; obtainable classically ([PALEY_BLOCKS.md](PALEY_BLOCKS.md)) |
| [18036](18036) | 4509 = 3³ · 167 | Paley's matrix of order 6012 in sixteen circulant blocks, T-matrices of order 3 (Goethals–Seidel array over ℤ₃ × ℤ₁₅₀₃) | New, beyond the tables; obtainable classically ([PALEY_BLOCKS.md](PALEY_BLOCKS.md)) |
| [18108](18108) | 4527 = 3² · 503 | Paley's matrix of order 2012 in sixteen circulant blocks, T-matrices of order 9 | New, beyond the tables; obtainable classically ([PALEY_BLOCKS.md](PALEY_BLOCKS.md)) |
| [18340](18340) | 4585 = 5 · 7 · 131 | Paley's matrix of order 524 in sixteen circulant blocks, T-matrices of order 35 | New, beyond the tables; obtainable classically ([PALEY_BLOCKS.md](PALEY_BLOCKS.md)) |
| [18964](18964) | 4741 = 11 · 431 | Paley's matrix of order 1724 in sixteen circulant blocks, T-matrices of order 11 | New, beyond the tables; obtainable classically ([PALEY_BLOCKS.md](PALEY_BLOCKS.md)) |
| [19028](19028) | 4757 = 67 · 71 | Cooper–Wallis: T-sequences of length 71, Williamson-type matrices of order 67 from near-Williamson matrices | New, beyond the tables; obtainable classically ([PALEY_BLOCKS.md](PALEY_BLOCKS.md)) |
| [19124](19124) | 4781 = 7 · 683 | Paley's matrix of order 2732 in sixteen circulant blocks, T-matrices of order 7 | New, beyond the tables; obtainable classically ([PALEY_BLOCKS.md](PALEY_BLOCKS.md)) |
| [19284](19284) | 4821 = 3 · 1607 | Paley's matrix of order 6428 in sixteen circulant blocks, T-matrices of order 3 | New, beyond the tables; obtainable classically ([PALEY_BLOCKS.md](PALEY_BLOCKS.md)) |
| [19364](19364) | 4841 = 47 · 103 | Cooper–Wallis: T-matrices of order 103, bordered Williamson-type matrices of order 47 | New, beyond the tables |
| [19932](19932) | 4983 = 3 · 11 · 151 | Paley's matrix of order 1812 in sixteen circulant blocks, T-matrices of order 11 | New, beyond the tables; obtainable classically ([PALEY_BLOCKS.md](PALEY_BLOCKS.md)) |
| [20084](20084) | 5021 | Miyamoto (q = 5021) from the matrix of order 5020 | First explicit construction |
| [20260](20260) | 5065 = 5 · 1013 | Paley's matrix of order 4052 in sixteen circulant blocks, T-matrices of order 5 | New, beyond the tables; obtainable classically ([PALEY_BLOCKS.md](PALEY_BLOCKS.md)) |
| [20436](20436) | 5109 = 3 · 13 · 131 | Paley's matrix of order 1572 in sixteen circulant blocks, T-matrices of order 13 | New, beyond the tables; obtainable classically ([PALEY_BLOCKS.md](PALEY_BLOCKS.md)) |
| [20884](20884) | 5221 = 23 · 227 | Paley's matrix of order 908 in sixteen circulant blocks, T-matrices of order 23 | New, beyond the tables; obtainable classically ([PALEY_BLOCKS.md](PALEY_BLOCKS.md)) |
| [21044](21044) | 5261 | Miyamoto (q = 5261) from the matrix of order [5260](5260) | New, beyond the tables; obtainable classically ([PALEY_BLOCKS.md](PALEY_BLOCKS.md)) |
| [21084](21084) | 5271 = 3 · 7 · 251 | Paley's matrix of order 7028 in sixteen circulant blocks, T-matrices of order 3 | New, beyond the tables; obtainable classically ([PALEY_BLOCKS.md](PALEY_BLOCKS.md)) |
| [21132](21132) | 5283 = 3² · 587 | Paley's matrix of order 2348 in sixteen circulant blocks, T-matrices of order 9 | New, beyond the tables; obtainable classically ([PALEY_BLOCKS.md](PALEY_BLOCKS.md)) |
| [21484](21484) | 5371 = 41 · 131 | Paley's matrix of order 524 in sixteen circulant blocks, T-matrices of order 41 | New, beyond the tables; obtainable classically ([PALEY_BLOCKS.md](PALEY_BLOCKS.md)) |
| [21732](21732) | 5433 = 3 · 1811 | Paley's matrix of order 7244 in sixteen circulant blocks, T-matrices of order 3 | New, beyond the tables; obtainable classically ([PALEY_BLOCKS.md](PALEY_BLOCKS.md)) |
| [21748](21748) | 5437 | Miyamoto (q = 5437) from a Cooper–Wallis matrix of order 5436 (T-matrices of order 151, Williamson matrices of order 9) | New, beyond the tables; obtainable classically (its input 5436 from [PALEY_BLOCKS.md](PALEY_BLOCKS.md)) |
| [21836](21836) | 5459 = 53 · 103 | Paley's matrix of order 212 in sixteen circulant blocks, T-matrices of order 103 | New, beyond the tables |
| [22132](22132) | 5533 = 11 · 503 | Paley's matrix of order 2012 in sixteen circulant blocks, T-matrices of order 11 | New, beyond the tables; obtainable classically ([PALEY_BLOCKS.md](PALEY_BLOCKS.md)) |
| [22412](22412) | 5603 = 13 · 431 | Paley's matrix of order 1724 in sixteen circulant blocks, T-matrices of order 13 | New, beyond the tables; obtainable classically ([PALEY_BLOCKS.md](PALEY_BLOCKS.md)) |
| [23020](23020) | 5755 = 5 · 1151 | Paley's matrix of order 4604 in sixteen circulant blocks, T-matrices of order 5 | New, beyond the tables; obtainable classically ([PALEY_BLOCKS.md](PALEY_BLOCKS.md)) |
| [23156](23156) | 5789 = 7 · 827 | Paley's matrix of order 3308 in sixteen circulant blocks, T-matrices of order 7 | New, beyond the tables; obtainable classically ([PALEY_BLOCKS.md](PALEY_BLOCKS.md)) |
| [23172](23172) | 5793 = 3 · 1931 | Paley's matrix of order 7724 in sixteen circulant blocks, T-matrices of order 3 | New, beyond the tables; obtainable classically ([PALEY_BLOCKS.md](PALEY_BLOCKS.md)) |
| [23724](23724) | 5931 = 3² · 659 | Paley's matrix of order 7908 in sixteen circulant blocks, T-matrices of order 3 (Goethals–Seidel array over ℤ₃ × ℤ₁₉₇₇) | New, beyond the tables; obtainable classically ([PALEY_BLOCKS.md](PALEY_BLOCKS.md)) |
| [24036](24036) | 6009 = 3 · 2003 | Paley's matrix of order 8012 in sixteen circulant blocks, T-matrices of order 3 | New, beyond the tables; obtainable classically ([PALEY_BLOCKS.md](PALEY_BLOCKS.md)) |
| [24124](24124) | 6031 | Yamada–Kiyasu (q = 6029) from a skew-Hadamard matrix of order 3016, the double of [skew 1508](../skew/1508) | New, beyond the tables |
| [24284](24284) | 6071 = 13 · 467 | Paley's matrix of order 1868 in sixteen circulant blocks, T-matrices of order 13 | New, beyond the tables; obtainable classically ([PALEY_BLOCKS.md](PALEY_BLOCKS.md)) |
| [24308](24308) | 6077 = 59 · 103 | Cooper–Wallis: T-matrices of order 103, near-Williamson matrices of order 59 | New, beyond the tables |
| [24412](24412) | 6103 | Yamada–Kiyasu (q = 6101) from [skew 3052](../skew/3052) | New, beyond the tables |
| [24628](24628) | 6157 = 47 · 131 | Paley's matrix of order 524 in sixteen circulant blocks, T-matrices of order 47 | New, beyond the tables; obtainable classically ([PALEY_BLOCKS.md](PALEY_BLOCKS.md)) |
| [24772](24772) | 6193 = 11 · 563 | Paley's matrix of order 2252 in sixteen circulant blocks, T-matrices of order 11 | New, beyond the tables; obtainable classically ([PALEY_BLOCKS.md](PALEY_BLOCKS.md)) |
| [24836](24836) | 6209 = 7 · 887 | Paley's matrix of order 3548 in sixteen circulant blocks, T-matrices of order 7 | New, beyond the tables; obtainable classically ([PALEY_BLOCKS.md](PALEY_BLOCKS.md)) |
| [25332](25332) | 6333 = 3 · 2111 | Paley's matrix of order 8444 in sixteen circulant blocks, T-matrices of order 3 | New, beyond the tables; obtainable classically ([PALEY_BLOCKS.md](PALEY_BLOCKS.md)) |
| [25508](25508) | 6377 = 7 · 911 | Paley's matrix of order 3644 in sixteen circulant blocks, T-matrices of order 7 | New, beyond the tables; obtainable classically ([PALEY_BLOCKS.md](PALEY_BLOCKS.md)) |
| [25676](25676) | 6419 = 7² · 131 | Paley's matrix of order 524 in sixteen circulant blocks, T-matrices of order 49 | New, beyond the tables; obtainable classically ([PALEY_BLOCKS.md](PALEY_BLOCKS.md)) |
| [25828](25828) | 6457 = 11 · 587 | Paley's matrix of order 2348 in sixteen circulant blocks, T-matrices of order 11 | New, beyond the tables; obtainable classically ([PALEY_BLOCKS.md](PALEY_BLOCKS.md)) |
| [25860](25860) | 6465 = 3 · 5 · 431 | Paley's matrix of order 5172 in sixteen circulant blocks, T-matrices of order 5 | New, beyond the tables; obtainable classically ([PALEY_BLOCKS.md](PALEY_BLOCKS.md)) |
| [26044](26044) | 6511 = 17 · 383 | Paley's matrix of order 1532 in sixteen circulant blocks, T-matrices of order 17 | New, beyond the tables; obtainable classically ([PALEY_BLOCKS.md](PALEY_BLOCKS.md)) |
| [26052](26052) | 6513 = 3 · 13 · 167 | Paley's matrix of order 2004 in sixteen circulant blocks, T-matrices of order 13 | New, beyond the tables; obtainable classically ([PALEY_BLOCKS.md](PALEY_BLOCKS.md)) |
| [26124](26124) | 6531 = 3 · 7 · 311 | Paley's matrix of order 8708 in sixteen circulant blocks, T-matrices of order 3 | New, beyond the tables; obtainable classically ([PALEY_BLOCKS.md](PALEY_BLOCKS.md)) |
| [26140](26140) | 6535 = 5 · 1307 | Paley's matrix of order 5228 in sixteen circulant blocks, T-matrices of order 5 | New, beyond the tables; obtainable classically ([PALEY_BLOCKS.md](PALEY_BLOCKS.md)) |
| [26156](26156) | 6539 = 13 · 503 | Paley's matrix of order 2012 in sixteen circulant blocks, T-matrices of order 13 | New, beyond the tables; obtainable classically ([PALEY_BLOCKS.md](PALEY_BLOCKS.md)) |
| [26300](26300) | 6575 = 5² · 263 | Paley's matrix of order 1052 in sixteen circulant blocks, T-matrices of order 25 | New, beyond the tables; obtainable classically ([PALEY_BLOCKS.md](PALEY_BLOCKS.md)) |
| [26332](26332) | 6583 = 29 · 227 | Paley's matrix of order 908 in sixteen circulant blocks, T-matrices of order 29 | New, beyond the tables; obtainable classically ([PALEY_BLOCKS.md](PALEY_BLOCKS.md)) |
| [26604](26604) | 6651 = 3² · 739 | Paley's matrix of order 8868 in sixteen circulant blocks, T-matrices of order 3 (Goethals–Seidel array over ℤ₃ × ℤ₂₂₁₇) | New, beyond the tables; obtainable classically ([PALEY_BLOCKS.md](PALEY_BLOCKS.md)) |
| [26740](26740) | 6685 = 5 · 7 · 191 | Paley's matrix of order 5348 in sixteen circulant blocks, T-matrices of order 5 | New, beyond the tables; obtainable classically ([PALEY_BLOCKS.md](PALEY_BLOCKS.md)) |
| [27036](27036) | 6759 = 3² · 751 | Paley's matrix of order 9012 in sixteen circulant blocks, T-matrices of order 3 (Goethals–Seidel array over ℤ₃ × ℤ₂₂₅₃) | New, beyond the tables; obtainable classically ([PALEY_BLOCKS.md](PALEY_BLOCKS.md)) |
| [27204](27204) | 6801 = 3 · 2267 | Paley's matrix of order 9068 in sixteen circulant blocks, T-matrices of order 3 | New, beyond the tables; obtainable classically ([PALEY_BLOCKS.md](PALEY_BLOCKS.md)) |
| [27524](27524) | 6881 = 7 · 983 | Paley's matrix of order 3932 in sixteen circulant blocks, T-matrices of order 7 | New, beyond the tables; obtainable classically ([PALEY_BLOCKS.md](PALEY_BLOCKS.md)) |
| [27604](27604) | 6901 = 67 · 103 | Cooper–Wallis: T-matrices of order 103, near-Williamson matrices of order 67 | New, beyond the tables |
| [27772](27772) | 6943 = 53 · 131 | Paley's matrix of order 524 in sixteen circulant blocks, T-matrices of order 53 | New, beyond the tables; obtainable classically ([PALEY_BLOCKS.md](PALEY_BLOCKS.md)) |
| [28148](28148) | 7037 = 31 · 227 | Paley's matrix of order 908 in sixteen circulant blocks, T-matrices of order 31 | New, beyond the tables; obtainable classically ([PALEY_BLOCKS.md](PALEY_BLOCKS.md)) |
| [28364](28364) | 7091 = 7 · 1013 | Paley's matrix of order 4052 in sixteen circulant blocks, T-matrices of order 7 | New, beyond the tables; obtainable classically ([PALEY_BLOCKS.md](PALEY_BLOCKS.md)) |
| [28932](28932) | 7233 = 3 · 2411 | Paley's matrix of order 9644 in sixteen circulant blocks, T-matrices of order 3 | New, beyond the tables; obtainable classically ([PALEY_BLOCKS.md](PALEY_BLOCKS.md)) |
| [29108](29108) | 7277 = 19 · 383 | Paley's matrix of order 1532 in sixteen circulant blocks, T-matrices of order 19 | New, beyond the tables; obtainable classically ([PALEY_BLOCKS.md](PALEY_BLOCKS.md)) |
| [29220](29220) | 7305 = 3 · 5 · 487 | Paley's matrix of order 9740 in sixteen circulant blocks, T-matrices of order 3 | New, beyond the tables; obtainable classically ([PALEY_BLOCKS.md](PALEY_BLOCKS.md)) |
| [29276](29276) | 7319 = 13 · 563 | Paley's matrix of order 2252 in sixteen circulant blocks, T-matrices of order 13 | New, beyond the tables; obtainable classically ([PALEY_BLOCKS.md](PALEY_BLOCKS.md)) |
| [29404](29404) | 7351 | Yamada–Kiyasu (q = 7349) from [skew 3676](../skew/3676) | First explicit construction |
| [29772](29772) | 7443 = 3² · 827 | Paley's matrix of order 3308 in sixteen circulant blocks, T-matrices of order 9 | New, beyond the tables; obtainable classically ([PALEY_BLOCKS.md](PALEY_BLOCKS.md)) |
| [30052](30052) | 7513 = 11 · 683 | Paley's matrix of order 2732 in sixteen circulant blocks, T-matrices of order 11 | New, beyond the tables; obtainable classically ([PALEY_BLOCKS.md](PALEY_BLOCKS.md)) |
| [30180](30180) | 7545 = 3 · 5 · 503 | Paley's matrix of order 2012 in sixteen circulant blocks, T-matrices of order 15 | New, beyond the tables; obtainable classically ([PALEY_BLOCKS.md](PALEY_BLOCKS.md)) |
| [30220](30220) | 7555 = 5 · 1511 | Paley's matrix of order 6044 in sixteen circulant blocks, T-matrices of order 5 | New, beyond the tables; obtainable classically ([PALEY_BLOCKS.md](PALEY_BLOCKS.md)) |
| [30460](30460) | 7615 = 5 · 1523 | Paley's matrix of order 6092 in sixteen circulant blocks, T-matrices of order 5 | New, beyond the tables; obtainable classically ([PALEY_BLOCKS.md](PALEY_BLOCKS.md)) |
| [30508](30508) | 7627 = 29 · 263 | Paley's matrix of order 1052 in sixteen circulant blocks, T-matrices of order 29 | New, beyond the tables; obtainable classically ([PALEY_BLOCKS.md](PALEY_BLOCKS.md)) |
| [30524](30524) | 7631 = 13 · 587 | Paley's matrix of order 2348 in sixteen circulant blocks, T-matrices of order 13 | New, beyond the tables; obtainable classically ([PALEY_BLOCKS.md](PALEY_BLOCKS.md)) |
| [30548](30548) | 7637 = 7 · 1091 | Paley's matrix of order 4364 in sixteen circulant blocks, T-matrices of order 7 | New, beyond the tables; obtainable classically ([PALEY_BLOCKS.md](PALEY_BLOCKS.md)) |
| [30916](30916) | 7729 = 59 · 131 | Paley's matrix of order 524 in sixteen circulant blocks, T-matrices of order 59 | New, beyond the tables; obtainable classically ([PALEY_BLOCKS.md](PALEY_BLOCKS.md)) |
| [31412](31412) | 7853 | Miyamoto (q = 7853) from a Cooper–Wallis matrix of order 7852 (T-matrices of order 151, Williamson matrices of order 13) | New, beyond the tables |
| [31932](31932) | 7983 = 3² · 887 | Paley's matrix of order 3548 in sixteen circulant blocks, T-matrices of order 9 | New, beyond the tables; obtainable classically ([PALEY_BLOCKS.md](PALEY_BLOCKS.md)) |
| [31956](31956) | 7989 = 3 · 2663 | Paley's matrix of order 10652 in sixteen circulant blocks, T-matrices of order 3 | New, beyond the tables; obtainable classically ([PALEY_BLOCKS.md](PALEY_BLOCKS.md)) |
| [32012](32012) | 8003 = 53 · 151 | Paley's matrix of order 212 in sixteen circulant blocks, T-matrices of order 151 | New, beyond the tables |
| [32228](32228) | 8057 = 7 · 1151 | Paley's matrix of order 4604 in sixteen circulant blocks, T-matrices of order 7 | New, beyond the tables; obtainable classically ([PALEY_BLOCKS.md](PALEY_BLOCKS.md)) |
| [32316](32316) | 8079 = 3 · 2693 | Paley's matrix of order 10772 in sixteen circulant blocks, T-matrices of order 3 | New, beyond the tables; obtainable classically ([PALEY_BLOCKS.md](PALEY_BLOCKS.md)) |
| [32692](32692) | 8173 = 11 · 743 | Paley's matrix of order 2972 in sixteen circulant blocks, T-matrices of order 11 | New, beyond the tables; obtainable classically ([PALEY_BLOCKS.md](PALEY_BLOCKS.md)) |
| [32756](32756) | 8189 = 19 · 431 | Paley's matrix of order 1724 in sixteen circulant blocks, T-matrices of order 19 | New, beyond the tables; obtainable classically ([PALEY_BLOCKS.md](PALEY_BLOCKS.md)) |
| [32796](32796) | 8199 = 3² · 911 | Paley's matrix of order 3644 in sixteen circulant blocks, T-matrices of order 9 | New, beyond the tables; obtainable classically ([PALEY_BLOCKS.md](PALEY_BLOCKS.md)) |
| [33068](33068) | 8267 = 7 · 1181 | Paley's matrix of order 4724 in sixteen circulant blocks, T-matrices of order 7 | New, beyond the tables; obtainable classically ([PALEY_BLOCKS.md](PALEY_BLOCKS.md)) |
| [33220](33220) | 8305 = 5 · 11 · 151 | Paley's matrix of order 3020 in sixteen circulant blocks, T-matrices of order 11 | New, beyond the tables; obtainable classically ([PALEY_BLOCKS.md](PALEY_BLOCKS.md)) |
| [34196](34196) | 8549 = 83 · 103 | Cooper–Wallis: T-matrices of order 103, bordered Williamson-type matrices of order 83 | New, beyond the tables |
| [34204](34204) | 8551 = 17 · 503 | Paley's matrix of order 2012 in sixteen circulant blocks, T-matrices of order 17 | New, beyond the tables; obtainable classically ([PALEY_BLOCKS.md](PALEY_BLOCKS.md)) |
| [34556](34556) | 8639 = 53 · 163 | Paley's matrix of order 212 in sixteen circulant blocks, T-matrices of order 163 | New, beyond the tables |
| [35236](35236) | 8809 = 23 · 383 | Paley's matrix of order 1532 in sixteen circulant blocks, T-matrices of order 23 | New, beyond the tables; obtainable classically ([PALEY_BLOCKS.md](PALEY_BLOCKS.md)) |
| [35484](35484) | 8871 = 3 · 2957 | Paley's matrix of order 11828 in sixteen circulant blocks, T-matrices of order 3 | New, beyond the tables; obtainable classically ([PALEY_BLOCKS.md](PALEY_BLOCKS.md)) |
| [35516](35516) | 8879 = 13 · 683 | Paley's matrix of order 2732 in sixteen circulant blocks, T-matrices of order 13 | New, beyond the tables; obtainable classically ([PALEY_BLOCKS.md](PALEY_BLOCKS.md)) |
| [35636](35636) | 8909 = 59 · 151 | Cooper–Wallis: T-matrices of order 151, near-Williamson matrices of order 59 | New, beyond the tables |
| [35756](35756) | 8939 = 7 · 1277 | Paley's matrix of order 5108 in sixteen circulant blocks, T-matrices of order 7 | New, beyond the tables; obtainable classically ([PALEY_BLOCKS.md](PALEY_BLOCKS.md)) |
| [36204](36204) | 9051 = 3 · 7 · 431 | Paley's matrix of order 5172 in sixteen circulant blocks, T-matrices of order 7 | New, beyond the tables; obtainable classically ([PALEY_BLOCKS.md](PALEY_BLOCKS.md)) |
| [36220](36220) | 9055 = 5 · 1811 | Paley's matrix of order 7244 in sixteen circulant blocks, T-matrices of order 5 | New, beyond the tables; obtainable classically ([PALEY_BLOCKS.md](PALEY_BLOCKS.md)) |
| [36388](36388) | 9097 = 11 · 827 | Paley's matrix of order 3308 in sixteen circulant blocks, T-matrices of order 11 | New, beyond the tables; obtainable classically ([PALEY_BLOCKS.md](PALEY_BLOCKS.md)) |
| [36596](36596) | 9149 = 7 · 1307 | Paley's matrix of order 5228 in sixteen circulant blocks, T-matrices of order 7 | New, beyond the tables; obtainable classically ([PALEY_BLOCKS.md](PALEY_BLOCKS.md)) |
| [36820](36820) | 9205 = 5 · 7 · 263 | Paley's matrix of order 1052 in sixteen circulant blocks, T-matrices of order 35 | New, beyond the tables; obtainable classically ([PALEY_BLOCKS.md](PALEY_BLOCKS.md)) |
| [36844](36844) | 9211 = 61 · 151 | Paley's matrix of order 244 in sixteen circulant blocks, T-matrices of order 151 | New, beyond the tables |
| [37204](37204) | 9301 = 71 · 131 | Paley's matrix of order 524 in sixteen circulant blocks, T-matrices of order 71 | New, beyond the tables; obtainable classically ([PALEY_BLOCKS.md](PALEY_BLOCKS.md)) |
| [37228](37228) | 9307 = 41 · 227 | Paley's matrix of order 908 in sixteen circulant blocks, T-matrices of order 41 | New, beyond the tables; obtainable classically ([PALEY_BLOCKS.md](PALEY_BLOCKS.md)) |
| [37436](37436) | 9359 = 7² · 191 | Paley's matrix of order 5348 in sixteen circulant blocks, T-matrices of order 7 (Goethals–Seidel array over ℤ₇ × ℤ₁₃₃₇) | New, beyond the tables; obtainable classically ([PALEY_BLOCKS.md](PALEY_BLOCKS.md)) |
| [37860](37860) | 9465 = 3 · 5 · 631 | Paley's matrix of order 12620 in sixteen circulant blocks, T-matrices of order 3 | New, beyond the tables; obtainable classically ([PALEY_BLOCKS.md](PALEY_BLOCKS.md)) |
| [38228](38228) | 9557 = 19 · 503 | Paley's matrix of order 2012 in sixteen circulant blocks, T-matrices of order 19 | New, beyond the tables; obtainable classically ([PALEY_BLOCKS.md](PALEY_BLOCKS.md)) |
| [38252](38252) | 9563 = 73 · 131 | Paley's matrix of order 524 in sixteen circulant blocks, T-matrices of order 73 | New, beyond the tables; obtainable classically ([PALEY_BLOCKS.md](PALEY_BLOCKS.md)) |
| [38292](38292) | 9573 = 3 · 3191 | Paley's matrix of order 12764 in sixteen circulant blocks, T-matrices of order 3 | New, beyond the tables; obtainable classically ([PALEY_BLOCKS.md](PALEY_BLOCKS.md)) |
| [38468](38468) | 9617 = 59 · 163 | Cooper–Wallis: T-matrices of order 163, near-Williamson matrices of order 59 | New, beyond the tables |
| [38636](38636) | 9659 = 13 · 743 | Paley's matrix of order 2972 in sixteen circulant blocks, T-matrices of order 13 | New, beyond the tables; obtainable classically ([PALEY_BLOCKS.md](PALEY_BLOCKS.md)) |
| [39012](39012) | 9753 = 3 · 3251 | Paley's matrix of order 13004 in sixteen circulant blocks, T-matrices of order 3 | New, beyond the tables; obtainable classically ([PALEY_BLOCKS.md](PALEY_BLOCKS.md)) |
| [39028](39028) | 9757 = 11 · 887 | Paley's matrix of order 3548 in sixteen circulant blocks, T-matrices of order 11 | New, beyond the tables; obtainable classically ([PALEY_BLOCKS.md](PALEY_BLOCKS.md)) |
| [39260](39260) | 9815 = 5 · 13 · 151 | Cooper–Wallis: T-matrices of order 151, near-Williamson matrices of order 65 | New, beyond the tables; obtainable classically ([PALEY_BLOCKS.md](PALEY_BLOCKS.md)) |
| [39540](39540) | 9885 = 3 · 5 · 659 | Paley's matrix of order 7908 in sixteen circulant blocks, T-matrices of order 5 | New, beyond the tables; obtainable classically ([PALEY_BLOCKS.md](PALEY_BLOCKS.md)) |
| [39772](39772) | 9943 = 61 · 163 | Paley's matrix of order 244 in sixteen circulant blocks, T-matrices of order 163 | New, beyond the tables |
| [39916](39916) | 9979 = 17 · 587 | Paley's matrix of order 2348 in sixteen circulant blocks, T-matrices of order 17 | New, beyond the tables; obtainable classically ([PALEY_BLOCKS.md](PALEY_BLOCKS.md)) |
| [40244](40244) | 10061 | Miyamoto (q = 10061) from the matrix of order [10060](10060) | New, beyond the tables; obtainable classically ([PALEY_BLOCKS.md](PALEY_BLOCKS.md)) |
| [45044](45044) | 11261 | Miyamoto (q = 11261) from the matrix of order [11260](11260) | New, beyond the tables; obtainable classically ([PALEY_BLOCKS.md](PALEY_BLOCKS.md)) |
| [51580](51580) | 12895 = 5 · 2579 | Yamada–Kiyasu (q = 12893) from a skew-Hadamard matrix of order 6448, [skew 1612](../skew/1612) doubled twice | New, beyond the tables |

## Known orders, independent matrices

| Order | n | Construction | First construction |
|------:|---|---|---|
| [2060](2060) | 515 = 5 · 103 | Cooper–Wallis: T-matrices of order 103, Williamson matrices of order 5 | Schneider, 2026-08-23 |
| [3532](3532) | 883 | Goethals–Seidel over ℤ₈₈₃, multiplier subgroup of order 21 | Đoković, arXiv:1008.2043 (2010) |
| [16892](16892) | 4223 = 41 · 103 | Paley's matrix of order 164 in sixteen circulant blocks, T-matrices of order 103 | SageMath development version (gh-42876, merged 2026-10-03) |
| [25132](25132) | 6283 = 61 · 103 | Paley's matrix of order 244 in sixteen circulant blocks, T-matrices of order 103 | SageMath development version (gh-42876, merged 2026-10-03) |
| [31316](31316) | 7829 | Miyamoto (q = 7829) from a Cooper–Wallis matrix of order 7828 (T-matrices of order 103, Williamson matrices of order 19) | SageMath development version (gh-42876, merged 2026-10-03) |
