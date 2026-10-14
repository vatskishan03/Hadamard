"""Census values reachable by Seberry's Theorem 7 (Ars Combin. 6 (1978)) with m >= 3: v = m k (4k - 1) with 4k an
amicable order (classical q + 1, q = 4k - 1 a prime power, or the new orders 52, 76, 100, 124) and m odd, among the
737 undecided skew values of the Đoković–Seberry census.  For each, the orthogonal design OD(4km; 1, m, 4km - m - 1)
that the theorem needs.  For m = 3, whether (t - 1)/2 (t = km) is a value of a^2 + ab + b^2 (necessary for the
Goethals–Seidel form, `attack/od13x.c`).

Usage: python3 attack/thm7_census.py CENSUS.txt   (pdftotext of the census PDF of 2025-05-28)
"""
import re, sys
from sympy import factorint, isprime
txt = open(sys.argv[1], encoding='utf-8', errors='replace').read()
hdr = r'CENSUS OF SYMMETRIC AND SKEW HADAMARD MATRICES OF ORDER 4v FOR ODD v < 2500\d*'
s = txt.index('Finally we list the 737 undecided cases of odd v < 2500:')
body = txt[s+len('Finally we list the 737 undecided cases of odd v < 2500:'):txt.index('4. The case v = 189')]
skew = set(int(x) for x in re.findall(r'\b\d+\b', re.sub(hdr, ' ', body)))
assert len(skew) == 737
def ppow(q):
    f = factorint(q); return len(f) == 1
amic = {}
for k in range(1, 200, 2):
    q = 4 * k - 1
    if ppow(q): amic[k] = 'q+1, q=%d' % q
for k in (13, 19, 25, 31): amic.setdefault(k, 'new (skew-propus %d)' % k)
def norm_ok(x):  # x = a^2 + ab + b^2 ?
    return any(a*a + a*b + b*b == x for a in range(0, 60) for b in range(-60, 60))
rows = []
for k in sorted(amic):
    base = k * (4 * k - 1)
    for m in range(3, 2500, 2):
        v = m * base
        if v >= 2500: break
        if v in skew:
            t = k * m
            extra = ''
            if m == 3: extra = 'singleton scheme norm (t-1)/2=%d: %s' % ((t - 1) // 2, norm_ok((t - 1) // 2))
            rows.append((v, k, m, 4 * k, amic[k], 4 * t, extra))
for r in sorted(rows):
    print('v=%4d  n=4k=%3d (%s)  m=%3d  needs OD(%d; 1, %d, %d)  %s' % (r[0], r[3], r[4], r[2], r[5], r[2], r[5] - r[2] - 1, r[6]))
print(len(rows), 'rows')
